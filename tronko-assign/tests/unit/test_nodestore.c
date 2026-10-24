/*
 * Unit test for the node store and the read-blocked scoring kernels (nodestore.c).
 *
 * The property the exactness argument rests on: for every job and every node, the new path
 * produces the same double, bit for bit, as production's loop (assignScores_Arr_paired and
 * getscore_Arr in assignment.c, linked here unchanged) accumulated into a zeroed nodeScores slot,
 * forward mate then reverse mate. Checked on random trees whose node counts straddle the tile
 * widths, with random posteriors that include sentinel columns, sentinels on A only, a C, G or T
 * value of 1, NaN and ties, and random rows that include -1 positions, '-', lowercase and other
 * characters, positions[0] == -1, one-mate and paired reads, candidates whose leaf string is empty,
 * blocks of 1, 7, 64 and 1000 reads and a block budget that closes blocks early, for every kernel
 * the CPU supports. Also: every leaf string from the store equals the legacy reader's, the store
 * equals the per-node arrays (ns_build's self-check), ns_reduce equals the reduction of
 * place_paired_with_nw, and trees the store path cannot reproduce are refused before anything
 * changes.
 *
 * Build and run: make test (in tronko-assign/), or
 *   gcc -O3 -std=gnu99 -I. -o tests/unit/test_nodestore tests/unit/test_nodestore.c nodestore.c \
 *       assignment.c getSequenceinRoot.c -lm -pthread && tests/unit/test_nodestore
 * Exit status 0 when every check passes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "global.h"
#include "assignment.h"
#include "getSequenceinRoot.h"
#include "nodestore.h"

/* the globals tronko-assign.c defines */
node **treeArr;
int *numspecArr, *numbaseArr, *rootArr;
char ****taxonomyArr;
type_of_PP Cinterval;

static long g_checks = 0, g_fail = 0, g_boundary_votes = 0;
#define CHECK(cond, ...) do { g_checks++; if (!(cond)) { g_fail++; if (g_fail <= 20) { fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); } } } while (0)

/* ---------------------------------------------------------------- random numbers */
static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void)
{
	uint64_t x = rng_state;
	x ^= x << 13; x ^= x >> 7; x ^= x << 17;
	return rng_state = x;
}
static int rndi(int n) { return (int)(rnd() % (uint64_t)n); }
static double rndu(void) { return (double)(rnd() >> 11) / 9007199254740992.0; }

/* bitwise equality, any two NaNs equal (their payloads may depend on operand order, which the
 * output never shows: a NaN score neither wins the maximum nor votes) */
static int same_double(double a, double b)
{
	if (isnan(a) && isnan(b)) return 1;
	return memcmp(&a, &b, sizeof(double)) == 0;
}

/* ---------------------------------------------------------------- trees */
typedef struct {
	int numspec, L;
	int empty_leaf;   /* a leaf whose every column is a sentinel, or -1 */
} tree_spec_t;

/* integer_values: posteriors of -1 to -4, so that sums tie exactly and nodes land exactly on the
 * vote window boundary maximum - Cinterval */
static int integer_values = 0;

static void make_tree(int t, int numspec, int L, int empty_leaf)
{
	int N = 2 * numspec - 1;
	numspecArr[t] = numspec;
	numbaseArr[t] = L;
	treeArr[t] = calloc((size_t)N, sizeof(node));
	for (int n = 0; n < N; n++) {
		treeArr[t][n].up[0] = treeArr[t][n].up[1] = -1;
		treeArr[t][n].down = -1;
		treeArr[t][n].posteriornc = calloc((size_t)L * 4 + 1, sizeof(type_of_PP));
	}
	/* topology: leaves are numspec-1 .. N-1, internal nodes 0 .. numspec-2, merged at random */
	int *pool = malloc(sizeof(int) * (size_t)N), np = 0;
	for (int n = numspec - 1; n < N; n++) pool[np++] = n;
	int *order = malloc(sizeof(int) * (size_t)(numspec > 1 ? numspec - 1 : 1));
	for (int i = 0; i < numspec - 1; i++) order[i] = i;
	for (int i = numspec - 2; i > 0; i--) { int j = rndi(i + 1); int x = order[i]; order[i] = order[j]; order[j] = x; }
	for (int i = 0; i < numspec - 1; i++) {
		int a = rndi(np); int ca = pool[a]; pool[a] = pool[--np];
		int b = rndi(np); int cb = pool[b]; pool[b] = pool[--np];
		int in = order[i];
		treeArr[t][in].up[0] = ca;
		treeArr[t][in].up[1] = cb;
		treeArr[t][ca].down = in;
		treeArr[t][cb].down = in;
		pool[np++] = in;
	}
	rootArr[t] = pool[0];
	free(pool);
	free(order);
	/* posteriors, as they are after store_PPs_Arr: log values, 1 for missing data */
	for (int n = 0; n < N; n++) {
		type_of_PP *pp = treeArr[t][n].posteriornc;
		int internal = n < numspec - 1;
		for (int p = 0; p < L; p++) {
			double *v = &pp[p * 4];
			int kind = rndi(1000);
			for (int b = 0; b < 4; b++) v[b] = integer_values ? -(double)(1 + rndi(4)) : log(0.0034 + 0.9866 * rndu());
			if (n == empty_leaf || kind < 350) {
				v[0] = v[1] = v[2] = v[3] = 1;              /* missing column */
			} else if (kind < 380) {
				v[0] = 1;                                   /* sentinel on A only */
			} else if (kind < 400) {
				v[1 + rndi(3)] = 1;                         /* a C, G or T of 1: added as 1.0 */
			} else if (kind < 460) {
				v[1 + rndi(3)] = v[0];                      /* ties for the leaf arg-max */
			} else if (kind < 470) {
				v[rndi(4)] = v[rndi(4)];
			} else if (kind < 475 && internal) {
				v[rndi(4)] = NAN;                           /* as in the CO1_fwhF2_EPTDr2n reference */
			}
		}
	}
}

static void free_tree(int t)
{
	int N = 2 * numspecArr[t] - 1;
	for (int n = 0; n < N; n++) free(treeArr[t][n].posteriornc);
	free(treeArr[t]);
}

/* ---------------------------------------------------------------- jobs and reads */
typedef struct {
	int alength;
	int *pos;
	char *q;
} job_t;

typedef struct {
	int leaf_iter;
	int coords[MAX_NUM_BWA_MATCHES][2];
	int has_job[MAX_NUM_BWA_MATCHES];   /* 0: empty leaf string, no job */
	job_t f[MAX_NUM_BWA_MATCHES], r[MAX_NUM_BWA_MATCHES];
	int paired;
	double *legacy[MAX_NUM_BWA_MATCHES];   /* nodeScores[k][tree][*] after both mates */
	double fwd_mm, rev_mm;
} read_t;

static const char QCHARS[] = "ACGTACGTACGTacgtacgt-NnR.";

static void make_job(job_t *j, int L)
{
	int kind = rndi(100);
	int cap = 2 * L + 8;
	j->pos = malloc(sizeof(int) * (size_t)(cap + 1));
	j->q = malloc((size_t)cap + 1);
	int a = 0;
	if (kind < 3) {                 /* alength 0: positions[0] stays -1 */
		j->pos[0] = -1;
		j->q[0] = '\0';
		j->alength = 0;
		return;
	}
	if (kind < 5) {                 /* positions[0] == -1 with rows after it: still the early exit */
		j->pos[a] = -1; j->q[a] = 'A'; a++;
	}
	int p = rndi(L > 4 ? L / 4 : 1);
	while (p < L && a < cap - 4) {
		int r = rndi(100);
		if (r < 12) { j->pos[a] = -1; j->q[a] = QCHARS[rndi(20)]; a++; continue; }   /* insertion */
		j->pos[a] = p;
		j->q[a] = QCHARS[rndi((int)sizeof(QCHARS) - 1)];
		a++;
		p += 1 + (rndi(10) == 0 ? rndi(3) : 0);
	}
	int tail = rndi(8) == 0 ? rndi(40) : 0;   /* the read runs past the leaf */
	for (int i = 0; i < tail && a < cap; i++) { j->pos[a] = -1; j->q[a] = "ACGT"[rndi(4)]; a++; }
	j->pos[a] = -1;
	j->q[a] = '\0';
	j->alength = a;
}

static void free_job(job_t *j) { free(j->pos); free(j->q); }

/* production's scoring of one candidate: zeroed slot, forward mate, then reverse mate */
static void legacy_scores(read_t *rd, int k, int ntrees)
{
	int t = rd->coords[k][0];
	int N = 2 * numspecArr[t] - 1;
	type_of_PP ***nodeScores = calloc(MAX_NUM_BWA_MATCHES, sizeof(type_of_PP **));
	nodeScores[k] = calloc((size_t)ntrees, sizeof(type_of_PP *));
	nodeScores[k][t] = calloc((size_t)N, sizeof(type_of_PP));
	if (rd->has_job[k]) {
		type_of_PP best = -9999999999999999; int strikes = 0;
		assignScores_Arr_paired(t, rootArr[t], rd->f[k].q, rd->f[k].pos, nodeScores, rd->f[k].alength, k, 0, NULL, NULL,
		                        0, &best, &strikes, 0, 0, 0, 0);
		if (rd->paired)
			assignScores_Arr_paired(t, rootArr[t], rd->r[k].q, rd->r[k].pos, nodeScores, rd->r[k].alength, k, 0, NULL, NULL,
			                        0, &best, &strikes, 0, 0, 0, 0);
	}
	rd->legacy[k] = nodeScores[k][t];
	free(nodeScores[k]);
	free(nodeScores);
}

/* the reduction of place_paired_with_nw (placement.c), on the legacy slots */
static void legacy_reduce(const read_t *rd, int **voteRoot, double *minimum)
{
	int i, j, k;
	type_of_PP maximum=-9999999999999999;
	for (i=0; i<rd->leaf_iter;i++){
		for (j=rd->coords[i][0]; j<rd->coords[i][0]+1; j++){
			for(k=0; k<2*numspecArr[j]-1; k++){
				if ( maximum < rd->legacy[i][k]){
					maximum=rd->legacy[i][k];
				}
			}
		}
	}
	for(i=0; i<rd->leaf_iter; i++){
		for(k=0; k<2*numspecArr[rd->coords[i][0]]-1; k++){
			if ( rd->legacy[i][k] >= (maximum-Cinterval) && rd->legacy[i][k] <= (maximum+Cinterval) ){
				voteRoot[rd->coords[i][0]][k]=1;
			}
		}
	}
	minimum[0] = maximum;
	minimum[1] = rd->fwd_mm;
	minimum[2] = rd->rev_mm;
}

/* ---------------------------------------------------------------- the blocked path, as the driver runs it */
static int **alloc_votes(int ntrees)
{
	int **v = malloc(sizeof(int *) * (size_t)ntrees);
	for (int t = 0; t < ntrees; t++) v[t] = calloc((size_t)(2 * numspecArr[t] - 1), sizeof(int));
	return v;
}
static void free_votes(int **v, int ntrees) { for (int t = 0; t < ntrees; t++) free(v[t]); free(v); }

static void compare_block(ns_block_t *B, read_t *reads, int ntrees, int **va, int **vb, const char *what)
{
	ns_block_score(B);
	int n = ns_block_nreads(B);
	for (int i = 0; i < n; i++) {
		ns_read_t *r = ns_block_read(B, i);
		read_t *rd = &reads[r->lineNumber];
		CHECK(r->leaf_iter == rd->leaf_iter, "%s: read %d leaf_iter", what, r->lineNumber);
		long bad = 0;
		for (int k = 0; k < rd->leaf_iter; k++) {
			const double *got = ns_read_scores(B, r, k);
			int N = 2 * numspecArr[rd->coords[k][0]] - 1;
			for (int x = 0; x < N; x++) if (!same_double(got[x], rd->legacy[k][x])) bad++;
		}
		CHECK(bad == 0, "%s: read %d: %ld node scores differ from the legacy loop", what, r->lineNumber, bad);
		/* the reduction */
		double ma[3], mb[3];
		for (int t = 0; t < ntrees; t++) {
			memset(va[t], 0, sizeof(int) * (size_t)(2 * numspecArr[t] - 1));
			memset(vb[t], 0, sizeof(int) * (size_t)(2 * numspecArr[t] - 1));
		}
		if (rd->leaf_iter > 0) {
			legacy_reduce(rd, va, ma);
			for (int k = 0; k < rd->leaf_iter; k++) {
				int N = 2 * numspecArr[rd->coords[k][0]] - 1;
				for (int x = 0; x < N; x++) if (rd->legacy[k][x] == ma[0] - Cinterval) g_boundary_votes++;
			}
			ns_reduce(B, r, vb, mb);
			CHECK(same_double(ma[0], mb[0]) && ma[1] == mb[1] && ma[2] == mb[2], "%s: read %d: reduction maximum", what, r->lineNumber);
			int vbad = 0;
			for (int t = 0; t < ntrees; t++)
				if (memcmp(va[t], vb[t], sizeof(int) * (size_t)(2 * numspecArr[t] - 1))) vbad++;
			CHECK(vbad == 0, "%s: read %d: votes differ", what, r->lineNumber);
		}
	}
	ns_block_reset(B);
}

static void run_blocked(read_t *reads, int nreads, int ntrees, int max_reads, size_t budget, const char *what)
{
	ns_block_t *B = ns_block_new(max_reads, budget);
	int **va = alloc_votes(ntrees), **vb = alloc_votes(ntrees);
	int *coord_rows[MAX_NUM_BWA_MATCHES];
	int coords[MAX_NUM_BWA_MATCHES][2];
	for (int k = 0; k < MAX_NUM_BWA_MATCHES; k++) coord_rows[k] = coords[k];
	for (int i = 0; i < nreads; i++) {
		read_t *rd = &reads[i];
		for (int k = 0; k < rd->leaf_iter; k++) { coords[k][0] = rd->coords[k][0]; coords[k][1] = rd->coords[k][1]; }
		if (!ns_block_fits(B, coord_rows, rd->leaf_iter)) compare_block(B, reads, ntrees, va, vb, what);
		ns_read_t *r = ns_block_add_read(B, i, i, coord_rows, rd->leaf_iter);
		/* place_paired_with_nw's order: every candidate's forward mate, then every reverse mate */
		for (int k = 0; k < rd->leaf_iter; k++)
			if (rd->has_job[k]) ns_sink_add_job(B, k, rd->f[k].pos, rd->f[k].q, rd->f[k].alength);
		if (rd->paired)
			for (int k = 0; k < rd->leaf_iter; k++)
				if (rd->has_job[k]) ns_sink_add_job(B, k, rd->r[k].pos, rd->r[k].q, rd->r[k].alength);
		r->fwd_mm = rd->fwd_mm;
		r->rev_mm = rd->rev_mm;
	}
	compare_block(B, reads, ntrees, va, vb, what);
	free_votes(va, ntrees);
	free_votes(vb, ntrees);
	ns_block_free(B);
}

/* ---------------------------------------------------------------- scenarios */

static int leaf_is_empty(int t, int leaf, char *seq, int *posr)
{
	getSequenceInNodeWithoutNs(t, leaf, seq, posr, 0, numbaseArr[t]);
	return seq[0] == '\0';
}

static void test_random(int seed_round)
{
	/* node counts around the tile widths 8, 32 and 64 */
	static const int specs[][2] = {
		{1, 40}, {2, 1}, {3, 50}, {4, 33}, {5, 70}, {9, 64}, {16, 90}, {17, 45}, {32, 120},
		{33, 80}, {48, 60}, {65, 100}, {100, 150}, {129, 70}, {300, 110},
	};
	const int ntrees = (int)(sizeof(specs) / sizeof(specs[0]));
	treeArr = calloc((size_t)ntrees, sizeof(node *));
	numspecArr = calloc((size_t)ntrees, sizeof(int));
	numbaseArr = calloc((size_t)ntrees, sizeof(int));
	rootArr = calloc((size_t)ntrees, sizeof(int));
	Cinterval = integer_values ? 2 : 5 + seed_round;
	for (int t = 0; t < ntrees; t++) {
		int numspec = specs[t][0], L = specs[t][1];
		int empty = (t % 3 == 1) ? numspec - 1 + rndi(numspec) : -1;   /* some trees get an all-missing leaf */
		make_tree(t, numspec, L, empty);
	}
	int maxL = 0;
	for (int t = 0; t < ntrees; t++) if (numbaseArr[t] > maxL) maxL = numbaseArr[t];
	char *seq = malloc((size_t)maxL + 1), *seq2 = malloc((size_t)maxL + 1);
	int *posr = malloc(sizeof(int) * ((size_t)maxL + 1)), *posr2 = malloc(sizeof(int) * ((size_t)maxL + 1));

	/* reads */
	const int nreads = 600;
	read_t *reads = calloc((size_t)nreads, sizeof(read_t));
	for (int i = 0; i < nreads; i++) {
		read_t *rd = &reads[i];
		rd->paired = rndi(4) != 0;
		int u = rndi(100);
		/* candidates are distinct trees: at most min(MAX_NUM_BWA_MATCHES, ntrees), whatever the cap */
		int maxc = MAX_NUM_BWA_MATCHES < ntrees ? MAX_NUM_BWA_MATCHES : ntrees;
		rd->leaf_iter = u < 4 ? 0 : u < 70 ? 1 : u < 90 ? 2 : 1 + rndi(maxc);
		int used[sizeof(specs) / sizeof(specs[0])] = {0};
		for (int k = 0; k < rd->leaf_iter; k++) {
			int t;
			do { t = rndi(100) < 50 ? ntrees - 1 - rndi(4) : rndi(ntrees); } while (used[t]);   /* distinct trees */
			used[t] = 1;
			int numspec = numspecArr[t];
			rd->coords[k][0] = t;
			rd->coords[k][1] = numspec - 1 + rndi(numspec);
			rd->has_job[k] = !leaf_is_empty(t, rd->coords[k][1], seq, posr);
			if (rd->has_job[k]) {
				make_job(&rd->f[k], numbaseArr[t]);
				if (rd->paired) make_job(&rd->r[k], numbaseArr[t]);
			}
			legacy_scores(rd, k, ntrees);
		}
		rd->fwd_mm = rndi(20);
		rd->rev_mm = rd->paired ? rndi(20) : 0;
	}

	/* legacy leaf strings for every leaf, before the store exists */
	int nleaves = 0;
	for (int t = 0; t < ntrees; t++) nleaves += numspecArr[t];
	char **leaf_seq = malloc(sizeof(char *) * (size_t)nleaves);
	int **leaf_pos = malloc(sizeof(int *) * (size_t)nleaves);
	for (int t = 0, li = 0; t < ntrees; t++)
		for (int n = numspecArr[t] - 1; n < 2 * numspecArr[t] - 1; n++, li++) {
			getSequenceInNodeWithoutNs(t, n, seq, posr, 0, numbaseArr[t]);
			leaf_seq[li] = strdup(seq);
			leaf_pos[li] = malloc(sizeof(int) * (strlen(seq) + 1));
			memcpy(leaf_pos[li], posr, sizeof(int) * strlen(seq));
		}

	/* the store, with its self-check; the per-node arrays are kept for the legacy reader */
	ns_build_stats_t st;
	int ok = ns_build(ntrees, 1 + seed_round % 3, 1, 1, &st);
	CHECK(ok == 1, "ns_build refused regular trees: %s", st.reason);
	if (!ok) exit(1);
	CHECK(st.checked && st.check_T_bad == 0 && st.check_S_bad == 0 && st.check_pad_bad == 0,
	      "store differs from the per-node arrays: T %ld S %ld padding %ld", st.check_T_bad, st.check_S_bad, st.check_pad_bad);
	CHECK(st.check_leaf_bad == 0 && st.check_leaves == nleaves, "self-check leaf strings: %ld of %ld differ", st.check_leaf_bad, st.check_leaves);
	CHECK(st.check_nan > 0, "the test posteriors should contain NaN");

	/* leaf strings from the store, whole length and random windows */
	long leaf_bad = 0, win_bad = 0;
	for (int t = 0, li = 0; t < ntrees; t++)
		for (int n = numspecArr[t] - 1; n < 2 * numspecArr[t] - 1; n++, li++) {
			ns_leaf_sequence(t, n, seq, posr, 0, numbaseArr[t]);
			size_t l = strlen(leaf_seq[li]);
			if (strcmp(seq, leaf_seq[li]) || memcmp(posr, leaf_pos[li], sizeof(int) * l)) leaf_bad++;
			int a = rndi(numbaseArr[t] + 1), b = a + rndi(numbaseArr[t] - a + 1);
			getSequenceInNodeWithoutNs(t, n, seq2, posr2, a, b);   /* legacy: arrays still present */
			ns_leaf_sequence(t, n, seq, posr, a, b);
			if (strcmp(seq, seq2) || memcmp(posr, posr2, sizeof(int) * strlen(seq))) win_bad++;
		}
	CHECK(leaf_bad == 0, "%ld of %d leaf strings from the store differ from the legacy reader", leaf_bad, nleaves);
	CHECK(win_bad == 0, "%ld leaf windows differ", win_bad);

	/* every kernel, every block shape */
	static const char *kernels[] = { "plain", "avx2", "avx512" };
	static const int blocks[] = { 1, 7, 64, 1000 };
	int nkernels = 0;
	for (size_t ki = 0; ki < sizeof(kernels) / sizeof(kernels[0]); ki++) {
		const char *kn = ns_select_kernel(kernels[ki]);
		if (!kn) { if (seed_round == 0) printf("  kernel %-6s not supported here, skipped\n", kernels[ki]); continue; }
		nkernels++;
		for (size_t bi = 0; bi < sizeof(blocks) / sizeof(blocks[0]); bi++) {
			char what[64];
			snprintf(what, sizeof what, "kernel %s block %d", kn, blocks[bi]);
			run_blocked(reads, nreads, ntrees, blocks[bi], (size_t)64 << 20, what);
		}
		char what[64];
		snprintf(what, sizeof what, "kernel %s budget 1 KiB", kn);
		run_blocked(reads, nreads, ntrees, 64, 1024, what);   /* blocks closed by the budget */
	}
	CHECK(nkernels >= 1, "no kernel available");

	/* release the per-node arrays as a real run does: the legacy reader now reads the store */
	ns_free_store(ntrees);
	ok = ns_build(ntrees, 2, 0, 0, &st);
	CHECK(ok == 1, "second ns_build failed");
	long null_bad = 0, redirect_bad = 0;
	for (int t = 0, li = 0; t < ntrees; t++)
		for (int n = 0; n < 2 * numspecArr[t] - 1; n++) {
			if (treeArr[t][n].posteriornc != NULL) null_bad++;
			if (n >= numspecArr[t] - 1) {
				getSequenceInNodeWithoutNs(t, n, seq, posr, 0, numbaseArr[t]);
				if (strcmp(seq, leaf_seq[li]) || memcmp(posr, leaf_pos[li], sizeof(int) * strlen(seq))) redirect_bad++;
				li++;
			}
		}
	CHECK(null_bad == 0, "%ld per-node arrays not released", null_bad);
	CHECK(redirect_bad == 0, "%ld leaf strings differ after the arrays were released", redirect_bad);
	ns_select_kernel(NULL);
	run_blocked(reads, nreads, ntrees, 64, (size_t)64 << 20, "default kernel after release");
	ns_free_store(ntrees);

	for (int i = 0; i < nleaves; i++) { free(leaf_seq[i]); free(leaf_pos[i]); }
	free(leaf_seq); free(leaf_pos);
	for (int i = 0; i < nreads; i++)
		for (int k = 0; k < reads[i].leaf_iter; k++) {
			if (reads[i].has_job[k]) { free_job(&reads[i].f[k]); if (reads[i].paired) free_job(&reads[i].r[k]); }
			free(reads[i].legacy[k]);
		}
	free(reads);
	for (int t = 0; t < ntrees; t++) free_tree(t);
	free(treeArr); free(numspecArr); free(numbaseArr); free(rootArr);
	free(seq); free(seq2); free(posr); free(posr2);
}

/* The reduction on crafted scores: nodes exactly on the lower edge of the vote window, one ulp
 * inside and outside it, NaN, and the maximum in the second candidate. Deterministic, so the
 * window's >= is covered whatever the random rounds happen to produce. */
static void test_reduce_edges(void)
{
	const int ntrees = 2;
	treeArr = calloc((size_t)ntrees, sizeof(node *));
	numspecArr = calloc((size_t)ntrees, sizeof(int));
	numbaseArr = calloc((size_t)ntrees, sizeof(int));
	rootArr = calloc((size_t)ntrees, sizeof(int));
	make_tree(0, 5, 10, -1);   /* 9 nodes */
	make_tree(1, 4, 10, -1);   /* 7 nodes */
	Cinterval = 2;
	ns_build_stats_t st;
	CHECK(ns_build(ntrees, 1, 1, 0, &st) == 1, "ns_build for the reduction test: %s", st.reason);
	static const double v0[9] = { -10.0, -12.0, -3.0, -4.0, -3.9999999999999996,
	                              NAN, -30.0, -7.0, -4.000000000000001 };
	static const double v1[7] = { -6.0, -4.0, -4.000000000000001, -2.0, -2.0000000000000004, NAN, -100.0 };
	read_t rd;
	memset(&rd, 0, sizeof rd);
	rd.leaf_iter = 2;
	rd.coords[0][0] = 0; rd.coords[0][1] = 4;
	rd.coords[1][0] = 1; rd.coords[1][1] = 3;
	rd.legacy[0] = (double *)v0;
	rd.legacy[1] = (double *)v1;
	rd.fwd_mm = 3; rd.rev_mm = 4;
	ns_block_t *B = ns_block_new(4, (size_t)1 << 20);
	int *rows[MAX_NUM_BWA_MATCHES], coords[MAX_NUM_BWA_MATCHES][2];
	for (int k = 0; k < MAX_NUM_BWA_MATCHES; k++) rows[k] = coords[k];
	coords[0][0] = 0; coords[0][1] = 4; coords[1][0] = 1; coords[1][1] = 3;
	ns_read_t *r = ns_block_add_read(B, 0, 0, rows, 2);
	r->fwd_mm = 3; r->rev_mm = 4;
	ns_block_score(B);   /* no jobs: both slots stay +0.0 */
	double *s0 = (double *)ns_read_scores(B, r, 0), *s1 = (double *)ns_read_scores(B, r, 1);
	for (int x = 0; x < 9; x++) CHECK(s0[x] == 0.0 && !signbit(s0[x]), "slot without jobs is not +0.0");
	for (int x = 0; x < 7; x++) CHECK(s1[x] == 0.0 && !signbit(s1[x]), "slot without jobs is not +0.0");
	memcpy(s0, v0, sizeof v0);
	memcpy(s1, v1, sizeof v1);
	int **va = alloc_votes(ntrees), **vb = alloc_votes(ntrees);
	double ma[3], mb[3];
	legacy_reduce(&rd, va, ma);
	ns_reduce(B, r, vb, mb);
	CHECK(same_double(ma[0], mb[0]) && ma[0] == -2.0 && ma[1] == mb[1] && ma[2] == mb[2], "reduction edges: maximum");
	int votes = 0, vbad = 0;
	for (int t = 0; t < ntrees; t++)
		for (int n = 0; n < 2 * numspecArr[t] - 1; n++) { votes += va[t][n]; vbad += va[t][n] != vb[t][n]; }
	CHECK(vbad == 0, "reduction edges: %d votes differ", vbad);
	CHECK(votes == 6, "reduction edges: %d votes, expected 6 (the window [-4, 0])", votes);
	if (va[1][1] == 1 && va[1][2] == 0) g_boundary_votes++;
	free_votes(va, ntrees);
	free_votes(vb, ntrees);
	ns_block_free(B);
	ns_free_store(ntrees);
	for (int t = 0; t < ntrees; t++) free_tree(t);
	free(treeArr); free(numspecArr); free(numbaseArr); free(rootArr);
}

/* Trees the store path cannot reproduce must be refused before anything changes. */
static void test_irregular(void)
{
	static const char *names[] = { "one-child node", "child index out of range", "node reached twice",
	                               "cycle", "unreachable node", "root out of range" };
	for (int kind = 0; kind < 6; kind++) {
		const int ntrees = 2;
		treeArr = calloc((size_t)ntrees, sizeof(node *));
		numspecArr = calloc((size_t)ntrees, sizeof(int));
		numbaseArr = calloc((size_t)ntrees, sizeof(int));
		rootArr = calloc((size_t)ntrees, sizeof(int));
		make_tree(0, 6, 20, -1);
		make_tree(1, 5, 20, -1);   /* nodes 0..3 internal, 4..8 leaves */
		node *T1 = treeArr[1];
		int r = rootArr[1];
		switch (kind) {
		case 0: T1[r].up[1] = -1; break;
		case 1: T1[r].up[0] = 9; break;
		case 2: { int c0 = T1[r].up[0], c1 = T1[r].up[1]; (void)c1; T1[r].up[1] = c0; } break;
		case 3: { int c = T1[r].up[0]; while (T1[c].up[0] != -1) c = T1[c].up[0]; T1[c].up[0] = r; T1[c].up[1] = r; } break;
		case 4: rootArr[1] = T1[r].up[0]; if (T1[rootArr[1]].up[0] == -1) rootArr[1] = T1[r].up[1]; break;
		case 5: rootArr[1] = 42; break;
		}
		type_of_PP *before = treeArr[0][0].posteriornc;
		ns_build_stats_t st;
		int ok = ns_build(ntrees, 2, 0, 0, &st);
		CHECK(ok == 0, "irregular tree (%s) was accepted", names[kind]);
		CHECK(ns_trees == NULL && treeArr[0][0].posteriornc == before, "irregular tree (%s): state changed", names[kind]);
		if (ok) ns_free_store(ntrees);
		rootArr[1] = r;
		for (int t = 0; t < ntrees; t++) free_tree(t);
		free(treeArr); free(numspecArr); free(numbaseArr); free(rootArr);
	}
}

int main(void)
{
	printf("test_nodestore: node store and blocked kernels against production's scoring loop\n");
	CHECK(ns_options_eligible(1, 0, 0, 0, 0, 0, 0) == 1, "production options not eligible");
	CHECK(ns_options_eligible(0, 0, 0, 0, 0, 0, 0) == 0, "WFA (no -w) must keep the legacy path");
	CHECK(ns_options_eligible(1, 1, 0, 0, 0, 0, 0) == 0, "-e must keep the legacy path");
	CHECK(ns_options_eligible(1, 0, 1, 0, 0, 0, 0) == 0, "print_all_nodes must keep the legacy path");
	CHECK(ns_options_eligible(1, 0, 0, 1, 0, 0, 0) == 0, "pruning must keep the legacy path");
	CHECK(ns_select_kernel("no-such-kernel") == NULL, "unknown kernel accepted");
	for (int round = 0; round < 4; round++) {
		rng_state = 0x9E3779B97F4A7C15ULL + (uint64_t)round * 0x632BE59BD9B4E019ULL;
		integer_values = round == 3;
		test_random(round);
	}
	test_irregular();
	test_reduce_edges();
	CHECK(g_boundary_votes > 0, "no node score fell on the vote window boundary");
	printf("test_nodestore: %ld checks, %ld failed: %s\n", g_checks, g_fail, g_fail ? "FAIL" : "PASS");
	return g_fail ? 1 : 0;
}
