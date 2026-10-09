/*
 * nodestore.h - node-innermost posterior store and read-blocked scoring
 *
 * The production scoring loop (assignScores_Arr_paired and getscore_Arr in assignment.c) scores
 * one read against every node of a candidate tree, fetching each node's per-position quartet from
 * that node's own array. This module keeps, per tree, the same transformed posteriors with the
 * node index innermost, and scores a block of a thread's reads against tiles of nodes that stay in
 * cache. Every (read, node) score is the same sum of the same doubles in the same order, starting
 * from +0.0, so the output is byte-identical; the argument is written out at the top of
 * nodestore.c.
 *
 * The path is taken only for the production options (see ns_options_eligible) and only when every
 * tree passes the checks of ns_build; otherwise the run keeps the per-node arrays and the legacy
 * loop.
 */
#ifndef _NODESTORE_
#define _NODESTORE_

#include <stdint.h>
#include <stddef.h>
#include "global.h"

/* One tree in the node-innermost layout. */
typedef struct ns_tree {
	int N;        /* nodes, 2*numspecArr[t]-1 */
	int Npad;     /* N rounded up to a multiple of 8; lanes N..Npad-1 hold +0.0 */
	int L;        /* alignment columns, numbaseArr[t] */
	int SW;       /* bytes per sentinel row, Npad/8 */
	double *T;    /* T[(p*4+b)*Npad + n] = sigma(n,p) ? log(0.01) : posteriornc[n][p*4+b] */
	uint8_t *S;   /* S[p*SW + n/8] bit n%8 = sigma(n,p) = (posteriornc[n][p*4] == 1) */
} ns_tree_t;

extern ns_tree_t *ns_trees;   /* NULL unless ns_build succeeded */
extern int ns_max_npad;       /* largest Npad over all trees */

/* 1 when the command line is the production one this path reproduces. */
int ns_options_eligible(int use_nw, int use_leaf_portion, int print_all_nodes, int enable_pruning,
                        int early_termination, int print_alignments, int print_alignments_to_file);

typedef struct ns_build_stats {
	int built;              /* 1: the store exists and the per-node arrays are released */
	const char *reason;     /* why not, when built == 0 */
	const char *kernel;     /* kernel chosen for scoring */
	int rounds;             /* relayout rounds (groups of trees transposed together) */
	int trims;              /* malloc_trim calls */
	double store_gib;       /* bytes of T and S */
	double freed_gib;       /* bytes of per-node arrays released */
	int maps;               /* mmap regions holding the store (256 MiB slabs and large trees) */
	double t_transpose, t_free, t_trim;   /* seconds in each step of the relayout */
	/* self-check (TRONKO_NODESTORE_CHECK=1) */
	int checked;
	long check_values, check_T_bad, check_S_bad, check_pad_bad, check_nan;
	long check_leaves, check_leaf_bad;
} ns_build_stats_t;

/*
 * Build the store for trees 0..ntrees-1 from treeArr[t][n].posteriornc (after store_PPs_Arr), with
 * nthreads worker threads. Checks every tree first (a proper binary tree reached from rootArr[t],
 * every node added exactly once by the legacy traversal) and the scoring constants; if any check
 * fails nothing is changed and 0 is returned. On success the per-node arrays are freed and their
 * pointers set to NULL, unless keep_old is set (tests). check: compare the store with the per-node
 * arrays value by value and every leaf string with the legacy reader before releasing them.
 */
int ns_build(int ntrees, int nthreads, int keep_old, int check, ns_build_stats_t *st);
void ns_free_store(int ntrees);

/* getSequenceInNodeWithoutNs from the store: the sentinel from S, the arg-max from T. */
void ns_leaf_sequence(int t, int node, char *seq, int *posr, int start_position, int end_position);

/* ---------------------------------------------------------------- blocked scoring */

/* A (read, candidate) score vector: Npad doubles in the block arena, receiving its jobs in order. */
typedef struct ns_slot {
	size_t out;       /* offset in doubles into the arena, a multiple of 8 */
	int tree;
	int njob;         /* 0, 1 or 2: the forward mate's job, then the reverse mate's */
	int32_t code_off[2], ncode[2];
	int big[2];       /* positions[0] == -1: getscore_Arr returns 9999999999 for every node */
} ns_slot_t;

/* The per-read state that phase C needs, saved in phase A. */
typedef struct ns_read {
	int lineNumber, iter, leaf_iter;
	int coords[MAX_NUM_BWA_MATCHES][2];   /* results->leaf_coordinates[0..leaf_iter) */
	int slot[MAX_NUM_BWA_MATCHES];        /* slot index per candidate */
	type_of_PP fwd_mm, rev_mm;            /* forward and reverse mismatch counts */
} ns_read_t;

typedef struct ns_block ns_block_t;

typedef struct ns_counters {
	long blocks, reads, slots, jobs, rows, kernel_calls, budget_closes;
} ns_counters_t;

/* Kernel: every slot's jobs, tile by tile. Chosen once, before any placement thread starts. */
typedef void (*ns_kernel_fn)(const ns_tree_t *t, ns_slot_t *const *slots, int nslots,
                             const int32_t *codes, double *arena);
/* Select by name ("plain", "avx2", "avx512") or, for NULL or "", the best the CPU supports.
 * Returns the selected name, or NULL if the name is unknown or the CPU lacks it. */
const char *ns_select_kernel(const char *name);
extern ns_kernel_fn ns_kernel;

ns_block_t *ns_block_new(int max_reads, size_t budget_bytes);
void ns_block_free(ns_block_t *B);
void ns_block_reset(ns_block_t *B);
int ns_block_nreads(const ns_block_t *B);
ns_read_t *ns_block_read(ns_block_t *B, int i);
const ns_counters_t *ns_block_counters(const ns_block_t *B);
/* 1 if a read with these candidates may join: the block is empty, or it has room for one more
 * read and its slots fit the arena budget. */
int ns_block_fits(ns_block_t *B, int **leaf_coordinates, int leaf_iter);
/* Add a read: saves its candidates and gives each a zeroed slot. Its jobs arrive through
 * ns_sink_add_job until the next read is added. */
ns_read_t *ns_block_add_read(ns_block_t *B, int lineNumber, int iter, int **leaf_coordinates, int leaf_iter);
/* place_paired_with_nw's scoring call on this path: compile one mate's rows for candidate match. */
void ns_sink_add_job(ns_block_t *B, int match, const int *positions, const char *locQuery, int alength);
/* Phase B: score every slot of the block, tree by tree. */
void ns_block_score(ns_block_t *B);
/* The score vector of candidate i of read r. */
const double *ns_read_scores(const ns_block_t *B, const ns_read_t *r, int i);
/* Phase C: the reduction of place_paired_with_nw (placement.c), on the read's slots. */
void ns_reduce(const ns_block_t *B, const ns_read_t *r, int **voteRoot, type_of_PP *minimum_score);

#endif /* _NODESTORE_ */
