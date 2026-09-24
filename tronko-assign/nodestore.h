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

#endif /* _NODESTORE_ */
