/* The vote reduction every scoring path ends in: the node store (ns_reduce) and both legacy
 * placements (WFA and NW). Given each candidate tree's node scores, it finds the best score M
 * over the read's candidates and votes every node scoring within the window of it. One
 * implementation, so the paths cannot drift apart. */
#ifndef REDUCE_H
#define REDUCE_H

#include "global.h"

/* What the last reduction on this thread found. Filled only while g_reduce.collect is set
 * (--trace-read); finish_read, which runs right after the reduction on the same thread for
 * the same read, reads it. Leaves are nodes with no children; L1 is the best leaf (ties to
 * the lowest tree, then node), L2 the best leaf of a species other than L1's (-1 when there
 * is none). Scores are raw, as the window sees them. */
typedef struct reduce_info {
	int ncand;
	type_of_PP maximum, window;
	int max_tree, max_node;
	int votes;
	int l1_tree, l1_node, l2_tree, l2_node;
	type_of_PP l1, l2;
	int best_leaf;   /* 1 when --best-leaf-margin replaced the votes with L1 */
} reduce_info_t;

typedef struct reduce_cfg {
	int collect;     /* fill reduce_last() */
	int per_base;    /* --cinterval-per-base: window = -c x the read's A/C/G/T bases */
	int need_bases;  /* the callers count the read's bases (per_base, or a finite max_deficit) */
	int full_read;   /* --score-full-read: a read's leading overhang bases are scored rows too */
	/* --best-leaf-margin D (0 = off) and --best-leaf-max-deficit T (INFINITY = no bound): the
	 * votes become {L1} when L1 is voted, S(L1) - S(L2) >= D, and n*log_match - S(L1) <= T,
	 * n the read's A/C/G/T bases and log_match = log(1 - u), a base's score on a matching leaf. */
	type_of_PP margin, max_deficit, log_match;
} reduce_cfg_t;

/* getscore_Arr's score for a job with no aligned position (positions[0] == -1), on every
 * node: not evidence, so the new options never let it be the maximum. A paired slot adds the
 * other mate's score to it, so a slot holding it is any score above REDUCE_SENTINEL_FLOOR
 * (a real score is a sum of logs, at most 0). */
#define REDUCE_SENTINEL 9999999999.0
#define REDUCE_SENTINEL_FLOOR 1e9

/* Set once in main, before any placement thread starts; read-only after. */
extern reduce_cfg_t g_reduce;

/* Whether a scoring job of *alength* rows gets the sentinel: its row 0 is -1, which the row
 * builders never write, so a job with no aligned row. Under --score-full-read the rows start
 * with the read's leading overhang (-1 rows) and every read base has a row, so only a job with
 * no row at all. getscore_Arr and the node store both ask here. */
static inline int reduce_job_unaligned(const int *positions, int alength)
{
	return g_reduce.full_read ? alength == 0 : positions[0] == -1;
}

/* Vote the read's candidates: sc[i] holds tree[i]'s 2*numspecArr[tree[i]]-1 node scores.
 * Candidates are visited in order and nodes in index order, and the first strictly greater
 * score wins the maximum, as every path did. *read_bases* (A/C/G/T bases across both mates)
 * is read only when g_reduce.need_bases. Returns M. */
type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int read_bases, int **voteRoot);

/* A/C/G/T bases in a read (q2 may be NULL for single-end). */
int reduce_read_bases(const char *q1, const char *q2);

/* This thread's last reduction (valid while g_reduce.collect). */
const reduce_info_t *reduce_last(void);

/* A leaf's species name (the deepest rank of its taxonomy row). */
const char *reduce_leaf_species(int tree, int node);

#endif
