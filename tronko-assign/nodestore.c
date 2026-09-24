/*
 * nodestore.c - node-innermost posterior store and read-blocked scoring
 *
 * What production computes (assignment.c, placement.c):
 * for a job, one mate of one read aligned to the BWA leaf of one candidate tree, giving rows
 * (positions[i], locQuery[i]) for i = 0 .. alength-1, every node n of the tree receives
 *
 *     s(n) = 9999999999                                 if positions[0] == -1
 *     s(n) = (...((+0.0 + x_0(n)) + x_1(n)) + ...) + x_{alength-1}(n)     otherwise
 *
 * with one IEEE binary64 addition per row, in row order, where x_i(n) is what getscore_Arr adds:
 *
 *     positions[i] == -1                           log(0.01)
 *     sigma(n,p) and locQuery[i] == '-'            0          (score = score + 0)
 *     sigma(n,p), any other character              log(0.01)
 *     A, C, G or T in either case                  V[n][p][b]  (posteriornc after store_PPs_Arr)
 *     '-'                                          log(0.25)
 *     any other character                          nothing is added
 *
 * sigma(n,p) is posteriornc[n][p*4+0] == 1, the missing-data sentinel, tested on A only. The
 * traversal (assignScores_Arr_paired) adds s(n) into nodeScores[match][tree][n], which is +0.0
 * when the read starts; the forward mate's jobs all run before the reverse mate's, so a slot ends
 * as (+0.0 + s_F(n)) + s_R(n). The reduction then reads only these values.
 *
 * What this module computes, and why it is the same double:
 *
 * - The store holds T[(p*4+b)*Npad + n] = sigma(n,p) ? log(0.01) : V[n][p][b], a copy with no
 *   arithmetic, and the bit S(p,n) = sigma(n,p) taken from the same double with the same == 1.
 *   A nucleotide row adds T, which is log(0.01) exactly where getscore_Arr adds log(0.01) for a
 *   nucleotide and V exactly where it adds V. A -1 row adds log(0.01). A '-' row adds 0 or
 *   log(0.25) and another character adds log(0.01) or nothing, both chosen by S and written as
 *   the literal operations of getscore_Arr. The constants are checked at start-up against what
 *   getscore_Arr itself returns (ns_check_constants); if they differ the path is not used.
 * - Each (job, node) accumulator starts at +0.0 and receives its rows in row order, one addition
 *   each. Tiles only partition the nodes; SIMD lanes are different nodes; the lane-wise additions
 *   of SSE2, AVX2 and AVX-512 round exactly as a scalar addsd under the default MXCSR, which
 *   nothing in the program changes. The kernels contain no multiplication, so contraction into
 *   fused multiply-add cannot occur, and no row is merged, reordered or pre-summed.
 * - A slot is zeroed when its read joins the block (+0.0, what nodeScores holds), then receives
 *   slot = slot + s for each of its jobs in the order they were produced: forward, then reverse.
 *   A candidate whose leaf string is empty has no job and keeps +0.0, as nodeScores does.
 * - Every node is scored exactly once only if the legacy traversal adds to every node exactly
 *   once. ns_build checks that for every tree before anything changes (ns_tree_regular) and
 *   otherwise leaves the whole run on the legacy path.
 * - The leaf string (getSequenceInNodeWithoutNs) takes sigma from S and the arg-max from T,
 *   which equals V wherever sigma is false; the scan starts from A with strict <, as the legacy
 *   one does. The self-check (TRONKO_NODESTORE_CHECK=1) compares every leaf of every tree.
 * - Which reads share a block, the block size and the thread count change no operand and no
 *   order: a slot's value depends only on its own jobs and the tree.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <malloc.h>
#include <sys/mman.h>
#include "nodestore.h"
#include "assignment.h"
#include "getSequenceinRoot.h"

#define NS_BIG 9999999999.0            /* getscore_Arr: score=9999999999 when positions[0] == -1 */
#define NS_ROUND_BYTES (1UL << 30)     /* per-node bytes transposed per round before release */
#define NS_TRIM_BYTES (2UL << 30)      /* released bytes between malloc_trim calls */
#define NS_ITEM_BYTES (2UL << 20)      /* target work per relayout item */
#define NS_SLAB_BYTES (256UL << 20)    /* small trees are carved from slabs of this size */
#define NS_OWN_MAP_BYTES (32UL << 20)  /* trees at least this large get their own mapping */

ns_tree_t *ns_trees = NULL;
int ns_max_npad = 0;
static double ns_log001, ns_log025;    /* getscore_Arr's log(0.01) and log(0.25) */

int ns_options_eligible(int use_nw, int use_leaf_portion, int print_all_nodes, int enable_pruning,
                        int early_termination, int print_alignments, int print_alignments_to_file)
{
	(void)hashmap_base;   /* hashmap_base.h defines this static in every file that includes global.h */
#if defined(OPTIMIZE_MEMORY) || defined(__FAST_MATH__)
	(void)use_nw; (void)use_leaf_portion; (void)print_all_nodes; (void)enable_pruning;
	(void)early_termination; (void)print_alignments; (void)print_alignments_to_file;
	return 0;   /* float posteriors, or a build that may reassociate sums: legacy only */
#else
	return use_nw == 1 && use_leaf_portion == 0 && print_all_nodes == 0 && enable_pruning == 0 &&
	       early_termination == 0 && print_alignments == 0 && print_alignments_to_file == 0;
#endif
}

/* ---------------------------------------------------------------- memory for the store */

typedef struct { void *addr; size_t len; } ns_map_t;
static ns_map_t *ns_maps = NULL;
static int ns_nmaps = 0, ns_cap_maps = 0;
static char *ns_slab = NULL;
static size_t ns_slab_left = 0;

static void *ns_map(size_t len, int huge)
{
	void *p = mmap(NULL, len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (p == MAP_FAILED) {
		fprintf(stderr, "tronko-assign: out of memory mapping %zu bytes for the node store\n", len);
		exit(1);
	}
#ifdef MADV_HUGEPAGE
	if (huge) (void)madvise(p, len, MADV_HUGEPAGE);   /* advice only; failure is harmless */
#else
	(void)huge;
#endif
	if (ns_nmaps == ns_cap_maps) {
		ns_cap_maps = ns_cap_maps ? 2 * ns_cap_maps : 64;
		ns_maps = realloc(ns_maps, (size_t)ns_cap_maps * sizeof(ns_map_t));
		if (!ns_maps) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
	}
	ns_maps[ns_nmaps].addr = p;
	ns_maps[ns_nmaps].len = len;
	ns_nmaps++;
	return p;
}

/* zero-filled, 64-byte aligned; bytes is a multiple of 64 */
static void *ns_alloc(size_t bytes)
{
	if (bytes >= NS_OWN_MAP_BYTES) return ns_map(bytes, 1);
	if (bytes > ns_slab_left) {
		ns_slab = ns_map(NS_SLAB_BYTES, 0);
		ns_slab_left = NS_SLAB_BYTES;
	}
	void *p = ns_slab;
	ns_slab += bytes;
	ns_slab_left -= bytes;
	return p;
}

void ns_free_store(int ntrees)
{
	(void)ntrees;
	for (int i = 0; i < ns_nmaps; i++) munmap(ns_maps[i].addr, ns_maps[i].len);
	free(ns_maps);
	ns_maps = NULL;
	ns_nmaps = ns_cap_maps = 0;
	ns_slab = NULL;
	ns_slab_left = 0;
	free(ns_trees);
	ns_trees = NULL;
}

/* ---------------------------------------------------------------- checks before any change */

/*
 * The legacy traversal (assignScores_Arr_paired, assignment.c) visits the root, and at every
 * visited node: a leaf (both children -1) gets one addition; an internal node (both children not
 * -1) gets one addition and both children are visited; a node with exactly one child -1 is scored
 * and discarded and its subtree is not visited. Replayed here with an explicit stack. The store
 * path scores every node, so it is taken only when every node receives exactly one addition.
 * Child indices outside the tree, cycles and shared subtrees make the tree irregular.
 */
static int ns_tree_regular(int t, int *stack, unsigned char *cnt)
{
	int N = 2 * numspecArr[t] - 1;
	int root = rootArr[t];
	if (N < 1 || root < 0 || root >= N) return 0;
	memset(cnt, 0, (size_t)N);
	long visits = 0;
	int sp = 0;
	stack[sp++] = root;
	while (sp > 0) {
		int node = stack[--sp];
		if (node < 0 || node >= N) return 0;
		if (++visits > N) return 0;
		int c0 = treeArr[t][node].up[0], c1 = treeArr[t][node].up[1];
		if (c0 == -1 && c1 == -1) {
			cnt[node]++;
		} else if (c0 != -1 && c1 != -1) {
			cnt[node]++;
			if (sp + 2 > N + 1) return 0;
			stack[sp++] = c1;
			stack[sp++] = c0;
		} else {
			return 0;
		}
	}
	for (int n = 0; n < N; n++)
		if (cnt[n] != 1) return 0;
	return 1;
}

/*
 * The kernel's constants must be the doubles getscore_Arr adds. Asked of getscore_Arr itself on a
 * real node: rows (p,'N'), (-1,'A') at a non-sentinel column p add nothing and then log(0.01);
 * row (p,'-') adds log(0.25). (A one-row job at position -1 would take the early exit.) If no
 * column of any tree is non-sentinel, log(0.25) can never be added and log(0.01) is checked with
 * (p,'-'), (-1,'A') at a sentinel column, which adds 0 and then log(0.01).
 */
static int ns_check_constants(int ntrees)
{
	int positions[4];
	char q[4];
	for (int t = 0; t < ntrees; t++) {
		int N = 2 * numspecArr[t] - 1, L = numbaseArr[t];
		for (int n = 0; n < N; n++) {
			const type_of_PP *pp = treeArr[t][n].posteriornc;
			for (int p = 0; p < L; p++) {
				if (pp[PP_IDX(p, 0)] == 1) continue;
				double v1, v2;
				positions[0] = p; positions[1] = -1; positions[2] = -1;
				q[0] = 'N'; q[1] = 'A'; q[2] = '\0';
				v1 = getscore_Arr(2, n, t, q, positions, 0, NULL, NULL);
				positions[0] = p; positions[1] = -1;
				q[0] = '-'; q[1] = '\0';
				v2 = getscore_Arr(1, n, t, q, positions, 0, NULL, NULL);
				return memcmp(&v1, &ns_log001, sizeof(double)) == 0 &&
				       memcmp(&v2, &ns_log025, sizeof(double)) == 0;
			}
		}
	}
	for (int t = 0; t < ntrees; t++) {
		if (numbaseArr[t] < 1) continue;
		double v1;
		positions[0] = 0; positions[1] = -1; positions[2] = -1;
		q[0] = '-'; q[1] = 'A'; q[2] = '\0';
		v1 = getscore_Arr(2, 0, t, q, positions, 0, NULL, NULL);
		return memcmp(&v1, &ns_log001, sizeof(double)) == 0;
	}
	return 1;   /* no column at all: no constant is ever added */
}

/* ---------------------------------------------------------------- relayout */

typedef struct { int t, p0, p1; } ns_item_t;
typedef struct {
	const ns_item_t *items;
	int nitems;
	int next;       /* shared counter, atomic */
} ns_work_t;

static size_t ns_old_bytes(int t)
{
	return (size_t)(2 * numspecArr[t] - 1) * (size_t)numbaseArr[t] * 4 * sizeof(type_of_PP);
}

/* Write positions [p0, p1) of tree t: T rows (p*4+b), padding lanes, and S rows p. Items of one
 * round own disjoint position ranges, so every byte of T and S has one writer. */
static void ns_transpose_item(const ns_item_t *it)
{
	ns_tree_t *s = &ns_trees[it->t];
	const int N = s->N, Npad = s->Npad, SW = s->SW;
	const size_t Np = (size_t)Npad;
	const double k01 = ns_log001;
	node *nodes = treeArr[it->t];
	for (int p = it->p0; p < it->p1; p++) {
		memset(s->S + (size_t)p * SW, 0, (size_t)SW);
		for (int b = 0; b < 4; b++)
			for (int n = N; n < Npad; n++) s->T[((size_t)p * 4 + b) * Np + n] = 0.0;
	}
	for (int n0 = 0; n0 < N; n0 += 64) {
		int n1 = n0 + 64 < N ? n0 + 64 : N;
		for (int n = n0; n < n1; n++) {
			const type_of_PP *src = nodes[n].posteriornc;
			for (int p = it->p0; p < it->p1; p++) {
				const type_of_PP *v = src + PP_IDX(p, 0);
				double *dst = s->T + (size_t)p * 4 * Np + n;
				if (v[0] == 1) {   /* the sentinel test of getscore_Arr, on the same double */
					s->S[(size_t)p * SW + (n >> 3)] |= (uint8_t)(1u << (n & 7));
					dst[0] = k01;
					dst[Np] = k01;
					dst[2 * Np] = k01;
					dst[3 * Np] = k01;
				} else {
					dst[0] = v[0];
					dst[Np] = v[1];
					dst[2 * Np] = v[2];
					dst[3 * Np] = v[3];
				}
			}
		}
	}
}

static void *ns_relayout_worker(void *arg)
{
	ns_work_t *w = arg;
	for (;;) {
		int i = __atomic_fetch_add(&w->next, 1, __ATOMIC_RELAXED);
		if (i >= w->nitems) break;
		ns_transpose_item(&w->items[i]);
	}
	return NULL;
}

static int ns_sbit(const ns_tree_t *s, int p, int n)
{
	return (s->S[(size_t)p * s->SW + ((size_t)n >> 3)] >> (n & 7)) & 1;
}

/* Self-check of one tree against its per-node arrays (TRONKO_NODESTORE_CHECK=1). */
static void ns_check_tree(int t, ns_build_stats_t *st, char *seq_a, int *pos_a, char *seq_b, int *pos_b)
{
	const ns_tree_t *s = &ns_trees[t];
	const size_t Np = (size_t)s->Npad;
	for (int n = 0; n < s->N; n++) {
		const type_of_PP *pp = treeArr[t][n].posteriornc;
		for (int p = 0; p < s->L; p++) {
			int sg = pp[PP_IDX(p, 0)] == 1;
			if (ns_sbit(s, p, n) != sg) st->check_S_bad++;
			for (int b = 0; b < 4; b++) {
				double want = sg ? ns_log001 : (double)pp[PP_IDX(p, b)];
				double got = s->T[((size_t)p * 4 + b) * Np + n];
				if (memcmp(&want, &got, sizeof(double)) != 0) st->check_T_bad++;
				if (isnan(pp[PP_IDX(p, b)])) st->check_nan++;
				st->check_values++;
			}
		}
	}
	for (int p = 0; p < s->L; p++) {
		for (int n = s->N; n < s->Npad; n++) {
			if (ns_sbit(s, p, n)) st->check_pad_bad++;
			for (int b = 0; b < 4; b++) {
				double z = s->T[((size_t)p * 4 + b) * Np + n];
				if (z != 0.0 || signbit(z)) st->check_pad_bad++;
			}
		}
	}
	/* every leaf: the legacy reader (per-node arrays still present) against the store */
	int numspec = numspecArr[t];
	for (int n = numspec - 1; n < 2 * numspec - 1; n++) {
		getSequenceInNodeWithoutNs(t, n, seq_a, pos_a, 0, s->L);
		ns_leaf_sequence(t, n, seq_b, pos_b, 0, s->L);
		size_t la = strlen(seq_a);
		if (la != strlen(seq_b) || memcmp(seq_a, seq_b, la) != 0 || memcmp(pos_a, pos_b, la * sizeof(int)) != 0)
			st->check_leaf_bad++;
		st->check_leaves++;
	}
}

int ns_build(int ntrees, int nthreads, int keep_old, int check, ns_build_stats_t *st)
{
	memset(st, 0, sizeof(*st));
	st->reason = "not attempted";
	ns_log001 = tronko_log001();
	ns_log025 = tronko_log025();
	if (ntrees < 1) { st->reason = "no trees"; return 0; }
	if (nthreads < 1) nthreads = 1;

	/* 1. every tree, before anything changes */
	int maxN = 0, maxL = 0;
	for (int t = 0; t < ntrees; t++) {
		if (numspecArr[t] < 1 || numbaseArr[t] < 0) { st->reason = "a tree has no species or a negative length"; return 0; }
		if (2 * numspecArr[t] - 1 > maxN) maxN = 2 * numspecArr[t] - 1;
		if (numbaseArr[t] > maxL) maxL = numbaseArr[t];
	}
	int *stack = malloc(((size_t)maxN + 2) * sizeof(int));
	unsigned char *cnt = malloc((size_t)maxN + 1);
	if (!stack || !cnt) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
	int regular = 1;
	for (int t = 0; t < ntrees && regular; t++) regular = ns_tree_regular(t, stack, cnt);
	free(stack);
	free(cnt);
	if (!regular) { st->reason = "a tree is not a proper binary tree reached once from its root"; return 0; }
	if (!ns_check_constants(ntrees)) { st->reason = "log(0.01) or log(0.25) differs from getscore_Arr"; return 0; }

	/* 2. shapes */
	ns_trees = calloc((size_t)ntrees, sizeof(ns_tree_t));
	if (!ns_trees) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
	ns_max_npad = 0;
	for (int t = 0; t < ntrees; t++) {
		ns_tree_t *s = &ns_trees[t];
		s->N = 2 * numspecArr[t] - 1;
		s->Npad = (s->N + 7) & ~7;
		s->L = numbaseArr[t];
		s->SW = s->Npad / 8;
		if (s->Npad > ns_max_npad) ns_max_npad = s->Npad;
	}

	char *seq_a = NULL, *seq_b = NULL;
	int *pos_a = NULL, *pos_b = NULL;
	if (check) {
		seq_a = malloc((size_t)maxL + 1); seq_b = malloc((size_t)maxL + 1);
		pos_a = malloc(((size_t)maxL + 1) * sizeof(int)); pos_b = malloc(((size_t)maxL + 1) * sizeof(int));
		if (!seq_a || !seq_b || !pos_a || !pos_b) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
	}

	/* 3. rounds: transpose a group of trees with all threads, then release their arrays */
	ns_item_t *items = NULL;
	int cap_items = 0;
	size_t since_trim = 0;
	pthread_t *th = malloc((size_t)nthreads * sizeof(pthread_t));
	if (!th) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
	for (int ta = 0, tb; ta < ntrees; ta = tb) {
		size_t bytes = 0;
		for (tb = ta; tb < ntrees && (tb == ta || bytes + ns_old_bytes(tb) <= NS_ROUND_BYTES); tb++)
			bytes += ns_old_bytes(tb);
		int nitems = 0;
		for (int t = ta; t < tb; t++) {
			ns_tree_t *s = &ns_trees[t];
			size_t tbytes = (size_t)4 * (size_t)s->L * (size_t)s->Npad * sizeof(double);
			size_t sbytes = ((size_t)s->L * (size_t)s->SW + 63) & ~(size_t)63;
			char *mem = ns_alloc(tbytes + sbytes > 0 ? tbytes + sbytes : 64);
			s->T = (double *)mem;
			s->S = (uint8_t *)(mem + tbytes);
			st->store_gib += (double)(tbytes + sbytes) / (1024.0 * 1024.0 * 1024.0);
			if (s->L == 0) continue;
			size_t per_pos = (size_t)s->N * 4 * sizeof(double);
			int chunk = (int)(NS_ITEM_BYTES / per_pos);
			if (chunk < 8) chunk = 8;
			for (int p0 = 0; p0 < s->L; p0 += chunk) {
				if (nitems == cap_items) {
					cap_items = cap_items ? 2 * cap_items : 1024;
					items = realloc(items, (size_t)cap_items * sizeof(ns_item_t));
					if (!items) { fprintf(stderr, "tronko-assign: out of memory\n"); exit(1); }
				}
				items[nitems].t = t;
				items[nitems].p0 = p0;
				items[nitems].p1 = p0 + chunk < s->L ? p0 + chunk : s->L;
				nitems++;
			}
		}
		ns_work_t w = { items, nitems, 0 };
		int nt = nthreads < nitems ? nthreads : nitems;
		int started = 0;
		for (int i = 1; i < nt; i++) {
			if (pthread_create(&th[i], NULL, ns_relayout_worker, &w) != 0) break;
			started = i;
		}
		ns_relayout_worker(&w);
		for (int i = 1; i <= started; i++) pthread_join(th[i], NULL);

		if (check) {
			for (int t = ta; t < tb; t++) ns_check_tree(t, st, seq_a, pos_a, seq_b, pos_b);
			st->checked = 1;
		}
		if (!keep_old) {
			for (int t = ta; t < tb; t++) {
				int N = 2 * numspecArr[t] - 1;
				for (int n = 0; n < N; n++) {
					free(treeArr[t][n].posteriornc);
					treeArr[t][n].posteriornc = NULL;
				}
				since_trim += ns_old_bytes(t);
				st->freed_gib += (double)ns_old_bytes(t) / (1024.0 * 1024.0 * 1024.0);
			}
			if (since_trim >= NS_TRIM_BYTES) {
				malloc_trim(0);
				st->trims++;
				since_trim = 0;
			}
		}
		st->rounds++;
	}
	if (!keep_old && since_trim > 0) {
		malloc_trim(0);
		st->trims++;
	}
	free(th);
	free(items);
	free(seq_a); free(seq_b); free(pos_a); free(pos_b);
	st->built = 1;
	st->reason = "";
	return 1;
}

/* getSequenceInNodeWithoutNs (getSequenceinRoot.c) on the store: a sentinel column emits nothing;
 * otherwise the arg-max over A, C, G, T starting from A with strict <, on the values of T, which
 * are the per-node values wherever the column is not a sentinel. */
void ns_leaf_sequence(int t, int node, char *seq, int *posr, int start_position, int end_position)
{
	const ns_tree_t *s = &ns_trees[t];
	const size_t Np = (size_t)s->Npad;
	int count = 0;
	for (int i = start_position; i < end_position; i++) {
		if (ns_sbit(s, i, node)) continue;
		const double *base = s->T + (size_t)i * 4 * Np + (size_t)node;
		double maximum = base[0];
		int index = 0;
		for (int j = 0; j < 4; j++) {
			if (maximum < base[(size_t)j * Np]) {
				maximum = base[(size_t)j * Np];
				index = j;
			}
		}
		seq[count] = "ACGT"[index];
		posr[count] = i;
		count++;
	}
	seq[count] = '\0';
}
