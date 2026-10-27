/* The vote reduction every scoring path ends in: the node store (ns_reduce) and both legacy
 * placements (WFA and NW). Given each candidate tree's node scores, it finds the best score M
 * over the read's candidates and votes every node scoring within the window of it. One
 * implementation, so the paths cannot drift apart. */
#ifndef REDUCE_H
#define REDUCE_H

#include "global.h"

/* Vote the read's candidates: sc[i] holds tree[i]'s 2*numspecArr[tree[i]]-1 node scores.
 * Candidates are visited in order and nodes in index order, and the first strictly greater
 * score wins the maximum, as every path did. Returns M. */
type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int **voteRoot);

#endif
