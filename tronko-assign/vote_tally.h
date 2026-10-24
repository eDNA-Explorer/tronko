/*
 * vote_tally.h: one read's vote tally, over the trees the read was scored against.
 *
 * After placement, voteRoot[t][k] is 1 for every node k of tree t whose score lies inside the
 * Cinterval window of the read's best score, and 0 everywhere else. Placement writes votes only
 * into the trees of the read's accepted candidate slots, and each accepted tree is recorded in
 * trees_search[0 .. min(leaf_iter, ntree)). So a pass over those trees in ascending tree order
 * gives what a pass over every node of every tree in the reference gives:
 *   - the same vote count per tree (every other tree has none);
 *   - the same minNodes sequence: tree by tree in ascending tree index, node by node in
 *     ascending node index;
 *   - the same maxRoot under the strict ">" rule from 0: the lowest tree index with the most
 *     votes, or -1 when no tree has any;
 *   - the same count of voting trees, and the same ascending list of them (maxRoots).
 * Resetting those trees' rows leaves voteRoot all zero, as resetting every row did.
 *
 * For a large reference this matters: CO1_fwhF2_EPTDr2n has 33,656 trees and 7,086,092 nodes,
 * which the full pass read and rewrote for every read pair, while a read has at most
 * MAX_NUM_BWA_MATCHES candidate trees.
 */
#ifndef VOTE_TALLY_H
#define VOTE_TALLY_H

typedef struct {
	int numMinNodes; /* entries written to minNodes: the read's votes */
	int count;       /* trees with at least one vote */
	int max;         /* votes of maxRoot; 0 when no tree has votes */
	int maxRoot;     /* lowest tree index with the most votes; -1 when no tree has votes */
} vote_tally_t;

/* The distinct values of trees_search[0 .. min(n, ntree)) that are tree ids (0 <= t < ntree),
   ascending, into hit. Returns how many there are. trees_search is in BWA hit order, not tree
   order, so the sort is needed. hit needs room for min(n, ntree) entries. */
static inline int vote_hit_trees(const int *trees_search, int n, int ntree, int *hit)
{
	int nhit = 0, i, k, m;
	for (i = 0; i < n && i < ntree; i++) {
		int t = trees_search[i];
		if (t < 0 || t >= ntree) continue;
		for (k = nhit; k > 0 && hit[k-1] > t; k--) {}
		if (k > 0 && hit[k-1] == t) continue;
		for (m = nhit; m > k; m--) hit[m] = hit[m-1];
		hit[k] = t;
		nhit++;
	}
	return nhit;
}

/* The tally over the trees hit[0 .. nhit), which must be ascending and distinct and must
   include every tree that has a vote. countVotes[k] receives the votes of tree hit[k], minNodes
   the voted node indices, and maxRoots the trees with votes, ascending. countVotes and maxRoots
   need room for nhit entries, minNodes for all the votes. */
static inline vote_tally_t vote_tally_hit_trees(int *const *voteRoot, const int *numspec,
                                                const int *hit, int nhit, int *countVotes,
                                                int *minNodes, int *maxRoots)
{
	vote_tally_t r = {0, 0, 0, -1};
	int k, j;
	for (k = 0; k < nhit; k++) {
		const int *votes = voteRoot[hit[k]];
		int nodes = 2*numspec[hit[k]]-1;
		countVotes[k] = 0;
		for (j = 0; j < nodes; j++) {
			if (votes[j] == 1) {
				countVotes[k]++;
				minNodes[r.numMinNodes] = j;
				r.numMinNodes++;
			}
		}
	}
	for (k = 0; k < nhit; k++) {
		if (countVotes[k] > r.max) {
			r.max = countVotes[k];
			r.maxRoot = hit[k];
		}
		if (countVotes[k] > 0) {
			maxRoots[r.count] = hit[k];
			r.count++;
		}
	}
	return r;
}

/* Clears the vote rows of the trees hit[0 .. nhit). */
static inline void vote_reset_hit_trees(int **voteRoot, const int *numspec, const int *hit, int nhit)
{
	int k, j;
	for (k = 0; k < nhit; k++) {
		int *votes = voteRoot[hit[k]];
		int nodes = 2*numspec[hit[k]]-1;
		for (j = 0; j < nodes; j++) votes[j] = 0;
	}
}

/* The tally as tronko-assign computed it before the hit-tree change, a pass over every node of
   every tree, kept as the reference for tests/unit/test_vote_tally.c and the VOTE_SHADOW_CHECK
   build. countVotes and maxRoots need room for ntree entries, minNodes for all the votes. */
static inline vote_tally_t vote_tally_all_trees(int *const *voteRoot, const int *numspec, int ntree,
                                                int *countVotes, int *minNodes, int *maxRoots)
{
	vote_tally_t r = {0, 0, 0, -1};
	int i, j, count = 0, count2 = 0;
	for (i = 0; i < ntree; i++) {
		countVotes[i] = 0;
		for (j = 0; j < 2*numspec[i]-1; j++) {
			if (voteRoot[i][j] == 1) {
				countVotes[i]++;
				minNodes[count] = j;
				count++;
			}
		}
	}
	r.numMinNodes = count;
	for (i = 0; i < ntree; i++) {
		if (countVotes[i] > r.max) {
			r.max = countVotes[i];
			r.maxRoot = i;
		}
		if (countVotes[i] > 0) {
			r.count++;
		}
	}
	for (i = 0; i < ntree; i++) {
		if (countVotes[i] > 0) {
			maxRoots[count2] = i;
			count2++;
		}
	}
	return r;
}

#ifdef VOTE_SHADOW_CHECK
/* Debug build (make clean && make ARCH_FLAGS=-DVOTE_SHADOW_CHECK): after each read's tally,
   recompute the full pass and compare every value; after the reset, check that every vote row of
   every tree is zero. Any difference prints the read and aborts. At exit, one line on stderr:
     VOTE_SHADOW reads=<n> multi_hit=<n> multi_vote=<n> unsorted=<n> neg_in_prefix=<n>
                 dup_in_prefix=<n> cap_reached=<n> written_past_prefix=<n>
   multi_hit: reads with more than one hit tree; multi_vote: with votes in more than one tree;
   unsorted: trees_search out of tree order; neg_in_prefix: a -1 or other invalid id inside the
   scanned prefix (expected 0, reachable only through BWA's candidate-list over-read);
   dup_in_prefix: a tree listed twice (expected 0 unless the reference has fewer trees than
   MAX_NUM_BWA_MATCHES); cap_reached: leaf_iter == MAX_NUM_BWA_MATCHES; written_past_prefix:
   a tree id at or beyond min(leaf_iter, ntree) (expected 0). About as slow as the code before
   the change, plus one more pass over the vote rows per read. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static long vote_shadow_n[8];
__attribute__((destructor)) static void vote_shadow_report(void)
{
	fprintf(stderr, "VOTE_SHADOW reads=%ld multi_hit=%ld multi_vote=%ld unsorted=%ld neg_in_prefix=%ld "
	        "dup_in_prefix=%ld cap_reached=%ld written_past_prefix=%ld\n", vote_shadow_n[0], vote_shadow_n[1],
	        vote_shadow_n[2], vote_shadow_n[3], vote_shadow_n[4], vote_shadow_n[5], vote_shadow_n[6], vote_shadow_n[7]);
}
static void vote_shadow_count(int which) { __atomic_add_fetch(&vote_shadow_n[which], 1, __ATOMIC_RELAXED); }
static void vote_shadow_fail(int read, const char *what)
{
	fprintf(stderr, "VOTE_SHADOW MISMATCH read %d: %s\n", read, what);
	abort();
}
static void vote_shadow_check(int *const *voteRoot, const int *numspec, int ntree, const int *trees_search,
                              int leaf_iter, int cap, const int *hit, int nhit, const int *countVotes,
                              vote_tally_t t, const int *minNodes, const int *maxRoots, int read)
{
	static __thread int *cv, *mn, *mr;
	int i, k, scan = leaf_iter < ntree ? leaf_iter : ntree;
	if (!cv) {
		long nodes = 0;
		for (i = 0; i < ntree; i++) nodes += 2*numspec[i]-1;
		cv = malloc(ntree * sizeof(int));
		mr = malloc(ntree * sizeof(int));
		mn = malloc(nodes * sizeof(int));
		if (!cv || !mr || !mn) vote_shadow_fail(read, "out of memory");
	}
	vote_tally_t o = vote_tally_all_trees(voteRoot, numspec, ntree, cv, mn, mr);
	if (o.numMinNodes != t.numMinNodes) vote_shadow_fail(read, "numMinNodes");
	if (memcmp(mn, minNodes, o.numMinNodes * sizeof(int)) != 0) vote_shadow_fail(read, "minNodes");
	if (o.count != t.count) vote_shadow_fail(read, "count");
	if (memcmp(mr, maxRoots, o.count * sizeof(int)) != 0) vote_shadow_fail(read, "maxRoots");
	if (o.maxRoot != t.maxRoot || o.max != t.max) vote_shadow_fail(read, "maxRoot or max");
	for (k = 0; k < nhit; k++) {
		if (countVotes[k] != cv[hit[k]]) vote_shadow_fail(read, "votes of a hit tree");
		cv[hit[k]] = 0;
	}
	for (i = 0; i < ntree; i++)
		if (cv[i] != 0) vote_shadow_fail(read, "votes outside the hit trees");
	vote_shadow_count(0);
	if (nhit > 1) vote_shadow_count(1);
	if (o.count > 1) vote_shadow_count(2);
	for (i = 1; i < scan; i++) if (trees_search[i] < trees_search[i-1]) { vote_shadow_count(3); break; }
	for (i = 0; i < scan; i++) if (trees_search[i] < 0 || trees_search[i] >= ntree) { vote_shadow_count(4); break; }
	for (i = 1, k = 0; i < scan && !k; i++) { int m; for (m = 0; m < i; m++) if (trees_search[m] == trees_search[i] && trees_search[i] >= 0) k = 1; }
	if (k) vote_shadow_count(5);
	if (leaf_iter == cap) vote_shadow_count(6);
	for (i = scan; i < ntree; i++) if (trees_search[i] != -1) { vote_shadow_count(7); break; }
}
static void vote_shadow_check_clean(int *const *voteRoot, const int *numspec, int ntree, int read)
{
	int i, j;
	for (i = 0; i < ntree; i++)
		for (j = 0; j < 2*numspec[i]-1; j++)
			if (voteRoot[i][j] != 0) vote_shadow_fail(read, "voteRoot not all zero after the reset");
}
#endif

#endif
