/*
 * test_vote_tally: the vote tally over a read's hit trees equals the tally over every tree.
 *
 * tronko-assign counts a read's votes with vote_tally_hit_trees() over the distinct candidate
 * trees of the read (vote_hit_trees()), ascending, and clears them with vote_reset_hit_trees().
 * Before, it passed over every node of every tree in the reference; that pass is kept as
 * vote_tally_all_trees(). This test runs both on random forests and random reads and requires,
 * for every read:
 *   - the same number of votes (numMinNodes) and the same minNodes[0 .. numMinNodes);
 *   - the same vote count for every tree (zero outside the hit trees);
 *   - the same count of voting trees, maxRoot and max (strict ">" from 0, lowest index wins);
 *   - the same maxRoots list, in the same order;
 *   - voteRoot all zero again after vote_reset_hit_trees().
 * The reads' candidate lists (trees_search) come in BWA hit order: out of tree order, with
 * duplicate trees, with -1 and out-of-range ids among them, longer than the tree count, and
 * empty. Votes lie only in trees the list names, which is what placement guarantees; the
 * vote patterns include none, one node, sparse, dense, every node, and equal counts in several
 * trees (the tie rule). The buffers are reused from read to read, as in tronko-assign, so a
 * reset that missed a row would show up in the next read.
 *
 * It also counts, as evidence that the inputs exercise the rules, the reads whose result would
 * change if the hit list were taken in hit order instead of sorted, or without removing
 * duplicates. Deterministic (fixed seed), no input files, about a second.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vote_tally.h"

#define MAX_CAP 25 /* production builds use MAX_NUM_BWA_MATCHES 25, the repository 10 */

static uint64_t rng = 0x2545F4914F6CDD1DULL;
static uint64_t next_u64(void) /* splitmix64 */
{
	uint64_t z = (rng += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}
static int below(int n) { return (int)(next_u64() % (uint64_t)n); }

typedef struct {
	const char *name;
	int ntree;
	int *numspec;
	int **voteRoot;
	long nodes;    /* total nodes, sum of 2*numspec-1 */
	int maxnodes;  /* nodes of the largest tree */
} forest_t;

static void forest_init(forest_t *f, const char *name, int ntree, int min_spec, int max_spec, int big_spec)
{
	int t;
	f->name = name;
	f->ntree = ntree;
	f->numspec = malloc(ntree * sizeof(int));
	f->voteRoot = malloc(ntree * sizeof(int *));
	f->nodes = 0;
	f->maxnodes = 0;
	for (t = 0; t < ntree; t++) {
		int spec = min_spec + below(max_spec - min_spec + 1);
		if (t == ntree / 3 && big_spec > 0) spec = big_spec; /* one large tree, as in MiFish */
		f->numspec[t] = spec;
		f->voteRoot[t] = calloc(2 * spec - 1, sizeof(int));
		f->nodes += 2 * spec - 1;
		if (2 * spec - 1 > f->maxnodes) f->maxnodes = 2 * spec - 1;
	}
}

static void forest_free(forest_t *f)
{
	int t;
	for (t = 0; t < f->ntree; t++) free(f->voteRoot[t]);
	free(f->voteRoot);
	free(f->numspec);
}

static int all_zero(const forest_t *f)
{
	int t, k;
	for (t = 0; t < f->ntree; t++)
		for (k = 0; k < 2 * f->numspec[t] - 1; k++)
			if (f->voteRoot[t][k] != 0) return 0;
	return 1;
}

static int cmp_int(const void *a, const void *b)
{
	int x = *(const int *)a, y = *(const int *)b;
	return (x > y) - (x < y);
}

/* The hit list by another route: collect the valid ids of the scanned prefix, qsort, unique. */
static int hit_reference(const int *ts, int n, int ntree, int *out)
{
	int m = 0, i, u = 0;
	for (i = 0; i < n && i < ntree; i++)
		if (ts[i] >= 0 && ts[i] < ntree) out[m++] = ts[i];
	qsort(out, m, sizeof(int), cmp_int);
	for (i = 0; i < m; i++)
		if (u == 0 || out[u - 1] != out[i]) out[u++] = out[i];
	return u;
}

static long reads, multi_hit, multi_vote, unsorted, dups, invalid, longer_than_ntree, empty, ties;
static long order_matters, dedup_matters, failures;

static void fail(const forest_t *f, long read, const char *what)
{
	if (failures < 20) fprintf(stderr, "FAIL %s read %ld: %s\n", f->name, read, what);
	failures++;
}

/* One read: candidate list, votes, both tallies, comparison, reset. */
static void one_read(forest_t *f, int *ts, int *old_cv, int *old_mn, int *old_mr, int *new_mn)
{
	int ntree = f->ntree, i, k, t;
	int cap = (int[]){1, 2, 10, 25}[below(4)];
	int leaf_iter = below(cap + 1);
	int scan = leaf_iter < ntree ? leaf_iter : ntree;
	int style = below(8);
	int hit[MAX_CAP], ref[MAX_CAP], cv[MAX_CAP], new_mr[MAX_CAP], seen_bad = 0, seen_dup = 0;
	long id = reads++;

	/* trees_search: a prefix of scan entries in hit order, -1 after it (as tronko-assign
	   leaves it). */
	for (i = 0; i < ntree; i++) ts[i] = -1;
	for (i = 0; i < scan; i++) {
		int v = below(ntree);
		if (style == 0) v = i;                                           /* already ascending */
		else if (style <= 2 && i > 0 && below(3) == 0) v = ts[below(i)]; /* a duplicate */
		else if (style == 3 && below(4) == 0) v = (int[]){-1, -2, ntree, ntree + 7}[below(4)];
		ts[i] = v;
	}
	if (style == 4) /* a burst of neighbouring trees, the usual shape of a real hit list */
		for (i = 0, t = below(ntree); i < scan; i++) ts[i] = (t + below(5)) % ntree;
	int nhit = vote_hit_trees(ts, leaf_iter, ntree, hit);
	int nref = hit_reference(ts, leaf_iter, ntree, ref);
	if (nhit != nref || memcmp(hit, ref, nhit * sizeof(int)) != 0) fail(f, id, "hit list differs from sort+unique");
	int seen_unsorted = 0;
	for (i = 0; i < scan; i++) {
		if (ts[i] < 0 || ts[i] >= ntree) seen_bad = 1;
		for (k = 0; k < i; k++) if (ts[k] == ts[i] && ts[i] >= 0) seen_dup = 1;
		if (i > 0 && ts[i] < ts[i - 1]) seen_unsorted = 1;
	}
	invalid += seen_bad;
	dups += seen_dup;
	unsorted += seen_unsorted;
	longer_than_ntree += leaf_iter > ntree;
	empty += nhit == 0;
	multi_hit += nhit > 1;

	/* Votes, only in listed trees. */
	int pattern = below(6), same = 1 + below(4);
	for (k = 0; k < nhit; k++) {
		int tr = hit[k], nodes = 2 * f->numspec[tr] - 1, p = pattern == 5 ? below(5) : pattern;
		int *v = f->voteRoot[tr];
		switch (p) {
		case 0: break;                                              /* no votes */
		case 1: v[below(nodes)] = 1; break;                         /* one node */
		case 2: for (i = 0; i < nodes; i++) v[i] = below(10) == 0; break; /* sparse */
		case 3: for (i = 0; i < nodes; i++) v[i] = below(3) != 0; break;  /* dense */
		case 4: for (i = 0; i < nodes; i++) v[i] = 1; break;        /* every node */
		}
		if (pattern == 0 && nodes >= same) /* equal counts in every tree: the tie rule */
			for (i = 0; i < same; i++) v[i * (nodes / same)] = 1;
	}

	/* Old: every tree. minNodes starts at -1 as the old code left it. */
	for (i = 0; i < (int)f->nodes; i++) old_mn[i] = -1;
	vote_tally_t o = vote_tally_all_trees(f->voteRoot, f->numspec, ntree, old_cv, old_mn, old_mr);
	/* New: hit trees only, into buffers holding garbage. */
	for (i = 0; i < (int)f->nodes; i++) new_mn[i] = 0x5a5a5a5a;
	vote_tally_t n = vote_tally_hit_trees(f->voteRoot, f->numspec, hit, nhit, cv, new_mn, new_mr);

	if (o.numMinNodes != n.numMinNodes) fail(f, id, "numMinNodes differs");
	else if (memcmp(old_mn, new_mn, o.numMinNodes * sizeof(int)) != 0) fail(f, id, "minNodes differs");
	if (o.count != n.count) fail(f, id, "count differs");
	else if (memcmp(old_mr, new_mr, o.count * sizeof(int)) != 0) fail(f, id, "maxRoots differs");
	if (o.maxRoot != n.maxRoot) fail(f, id, "maxRoot differs");
	if (o.max != n.max) fail(f, id, "max differs");
	{
		long in_hit = 0, all = 0;
		for (k = 0; k < nhit; k++) {
			if (cv[k] != old_cv[hit[k]]) fail(f, id, "vote count of a hit tree differs");
			in_hit += cv[k];
		}
		for (t = 0; t < ntree; t++) all += old_cv[t];
		if (all != in_hit) fail(f, id, "votes outside the hit trees");
	}
	multi_vote += o.count > 1;
	if (o.count > 1) {
		int c0 = -1, tie = 0;
		for (t = 0; t < ntree; t++) if (old_cv[t] == o.max) { if (c0 >= 0) tie = 1; c0 = t; }
		ties += tie;
	}

	/* Would the result change without the sort, or without the duplicate removal? */
	{
		int raw[MAX_CAP], nraw = 0, cvr[MAX_CAP], mrr[MAX_CAP];
		int *mnr = malloc(((size_t)scan * f->maxnodes + 1) * sizeof(int)); /* a repeated tree counts twice */
		for (i = 0; i < scan; i++) if (ts[i] >= 0 && ts[i] < ntree) raw[nraw++] = ts[i];
		vote_tally_t r = vote_tally_hit_trees(f->voteRoot, f->numspec, raw, nraw, cvr, mnr, mrr);
		if (r.numMinNodes != o.numMinNodes || r.maxRoot != o.maxRoot || r.count != o.count ||
		    memcmp(mnr, old_mn, o.numMinNodes * sizeof(int)) != 0 || memcmp(mrr, old_mr, o.count * sizeof(int)) != 0) {
			int sorted_unique = 1;
			for (i = 1; i < nraw; i++) if (raw[i] <= raw[i - 1]) sorted_unique = 0;
			if (!sorted_unique) {
				int d = 0;
				for (i = 0; i < nraw; i++) for (k = 0; k < i; k++) if (raw[k] == raw[i]) d = 1;
				if (d) dedup_matters++; else order_matters++;
			}
		}
		free(mnr);
	}

	vote_reset_hit_trees(f->voteRoot, f->numspec, hit, nhit);
	if ((id & 15) == 0 || f->nodes < 20000) {
		if (!all_zero(f)) fail(f, id, "voteRoot not all zero after the reset");
	}
}

int main(void)
{
	static const struct { const char *name; int ntree, min_spec, max_spec, big_spec, reads; } shapes[] = {
		{"one tree", 1, 1, 400, 0, 20000},
		{"one-leaf trees", 12, 1, 1, 0, 20000},
		{"three trees", 3, 50, 200, 0, 30000},
		{"five trees", 5, 1, 60, 0, 30000},
		{"MiFish-like, 68 trees", 68, 2, 120, 3000, 20000},
		{"FWH-like, 4000 small trees", 4000, 1, 16, 0, 20000},
	};
	int s, r;
	for (s = 0; s < (int)(sizeof(shapes) / sizeof(shapes[0])); s++) {
		forest_t f;
		forest_init(&f, shapes[s].name, shapes[s].ntree, shapes[s].min_spec, shapes[s].max_spec, shapes[s].big_spec);
		/* The old tally needs room for every tree and node; the new one gets the same so an
		   overrun would not go unnoticed under AddressSanitizer either way. */
		int *ts = malloc(f.ntree * sizeof(int));
		int *old_cv = malloc(f.ntree * sizeof(int));
		int *old_mr = malloc(f.ntree * sizeof(int));
		int *old_mn = malloc(f.nodes * sizeof(int));
		int *new_mn = malloc(f.nodes * sizeof(int));
		for (r = 0; r < shapes[s].reads; r++) one_read(&f, ts, old_cv, old_mn, old_mr, new_mn);
		if (!all_zero(&f)) fail(&f, reads, "voteRoot not all zero at the end");
		free(ts);
		free(old_cv); free(old_mr); free(old_mn); free(new_mn);
		forest_free(&f);
	}
	printf("test_vote_tally: %ld reads on %d forests; %ld with more than one hit tree, %ld with votes in more than one tree, "
	       "%ld ties for the most votes, %ld lists out of tree order, %ld with duplicate trees, %ld with invalid ids, "
	       "%ld longer than the tree count, %ld empty; without the sort %ld results would change, without the "
	       "duplicate removal %ld; %ld differ\n",
	       reads, (int)(sizeof(shapes) / sizeof(shapes[0])), multi_hit, multi_vote, ties, unsorted, dups, invalid,
	       longer_than_ntree, empty, order_matters, dedup_matters, failures);
	return failures ? 1 : 0;
}
