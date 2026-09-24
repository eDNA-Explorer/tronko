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

#endif
