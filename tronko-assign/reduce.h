/* The vote reduction every scoring path ends in: the node store (ns_reduce) and both legacy
 * placements (WFA and NW). Given each candidate tree's node scores, it finds the best score M
 * over the read's candidates and votes every node scoring within the window of it. One
 * implementation, so the paths cannot drift apart. */
#ifndef REDUCE_H
#define REDUCE_H

#include "global.h"

/* What the last reduction on this thread found. Filled only while g_reduce.collect is set
 * (--trace-read); finish_read, which runs right after the reduction on the same thread for
 * the same read, reads it. Leaves are nodes with no children; L2 is the best leaf of a
 * species other than L1's (-1 when there is none). Scores are raw, as the window sees them. */
typedef struct reduce_info {
	int ncand;
	type_of_PP maximum, window;
	int max_tree, max_node;
	int votes;
	int l1_tree, l1_node, l2_tree, l2_node;
	type_of_PP l1, l2;
} reduce_info_t;

typedef struct reduce_cfg {
	int collect;   /* fill reduce_last() */
} reduce_cfg_t;

/* Set once in main, before any placement thread starts; read-only after. */
extern reduce_cfg_t g_reduce;

/* Vote the read's candidates: sc[i] holds tree[i]'s 2*numspecArr[tree[i]]-1 node scores.
 * Candidates are visited in order and nodes in index order, and the first strictly greater
 * score wins the maximum, as every path did. Returns M. */
type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int **voteRoot);

/* This thread's last reduction (valid while g_reduce.collect). */
const reduce_info_t *reduce_last(void);

/* A leaf's species name (the deepest rank of its taxonomy row). */
const char *reduce_leaf_species(int tree, int node);

#endif
