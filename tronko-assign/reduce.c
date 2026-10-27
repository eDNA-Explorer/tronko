#include "reduce.h"

#include <ctype.h>
#include <math.h>
#include <string.h>

reduce_cfg_t g_reduce = {0};

static __thread reduce_info_t last;

const reduce_info_t *reduce_last(void)
{
	return &last;
}

static int is_leaf(int tree, int node)
{
	return treeArr[tree][node].up[0] == -1 && treeArr[tree][node].up[1] == -1;
}

const char *reduce_leaf_species(int tree, int node)
{
	return taxonomyArr[tree][treeArr[tree][node].taxIndex[0]][0];
}

int reduce_read_bases(const char *q1, const char *q2)
{
	int n = 0;
	const char *q;
	for (q = q1; q != NULL && *q; q++) {
		switch (toupper((unsigned char)*q)) {
		case 'A': case 'C': case 'G': case 'T': n++;
		}
	}
	for (q = q2; q != NULL && *q; q++) {
		switch (toupper((unsigned char)*q)) {
		case 'A': case 'C': case 'G': case 'T': n++;
		}
	}
	return n;
}

static int is_sentinel(type_of_PP s)
{
	return s >= REDUCE_SENTINEL;
}

/* Leaf (tree, node) scoring sc beats the current best (bt, bn, best): higher, or equal and
 * lower in (tree, node), so the choice does not depend on the order of the candidates. */
static int leaf_beats(type_of_PP sc, int tree, int node, type_of_PP best, int bt, int bn)
{
	if (sc != best) {
		return sc > best;
	}
	return bt < 0 || tree < bt || (tree == bt && node < bn);
}

/* L1, the best leaf, and L2, the best leaf naming another species: what the trace reports and
 * the best-leaf rule tests. The sentinel score is no leaf's evidence and never counts. */
static void best_leaves(int ncand, const type_of_PP *const *sc, const int *tree, reduce_info_t *info)
{
	int i, k;
	info->l1_tree = info->l1_node = info->l2_tree = info->l2_node = -1;
	info->l1 = info->l2 = -INFINITY;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if (is_leaf(tree[i], k) && !is_sentinel(sc[i][k])
			    && leaf_beats(sc[i][k], tree[i], k, info->l1, info->l1_tree, info->l1_node)){
				info->l1 = sc[i][k];
				info->l1_tree = tree[i];
				info->l1_node = k;
			}
		}
	}
	if (info->l1_tree < 0){
		return;
	}
	const char *species = reduce_leaf_species(info->l1_tree, info->l1_node);
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if (is_leaf(tree[i], k) && !is_sentinel(sc[i][k])
			    && leaf_beats(sc[i][k], tree[i], k, info->l2, info->l2_tree, info->l2_node)
			    && strcmp(reduce_leaf_species(tree[i], k), species) != 0){
				info->l2 = sc[i][k];
				info->l2_tree = tree[i];
				info->l2_node = k;
			}
		}
	}
}

/* --best-leaf-margin: L1 is voted, beats every other species' best leaf by the margin (a tie
 * never fires), and is within the deficit of a perfect match. */
static int best_leaf_fires(const reduce_info_t *info, type_of_PP maximum, type_of_PP window, int read_bases)
{
	if (info->l1_tree < 0 || info->l1 < maximum - window){
		return 0;
	}
	if (!(info->l1 - info->l2 >= g_reduce.margin)){
		return 0;
	}
	return !isfinite(g_reduce.max_deficit)
	       || (type_of_PP)read_bases * g_reduce.log_match - info->l1 <= g_reduce.max_deficit;
}

type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int read_bases, int **voteRoot)
{
	int i, k;
	type_of_PP maximum=-9999999999999999;
	int max_tree=-1, max_node=-1;
	int rule = g_reduce.margin > 0;
	/* Every new option takes the maximum over real scores only: the sentinel would win it. */
	int skip_sentinel = g_reduce.per_base || rule;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if ( maximum < sc[i][k] && !(skip_sentinel && is_sentinel(sc[i][k])) ){
				maximum=sc[i][k];
				max_tree=tree[i];
				max_node=k;
			}
		}
	}
	/* --cinterval-per-base: -c tolerated per read base, so one window per read (the same
	 * denominator for every candidate keeps the likelihood order). */
	type_of_PP window = g_reduce.per_base ? Cinterval * (type_of_PP)read_bases : Cinterval;
	reduce_info_t info;
	int fired = 0;
	if (rule || g_reduce.collect){
		best_leaves(ncand, sc, tree, &info);
		fired = rule && best_leaf_fires(&info, maximum, window, read_bases);
	}
	int votes=0;
	if (fired){
		voteRoot[info.l1_tree][info.l1_node]=1;
		votes=1;
	}else{
		for (i=0; i<ncand; i++){
			for (k=0; k<2*numspecArr[tree[i]]-1; k++){
				if ( sc[i][k] >= (maximum-window) && sc[i][k] <= (maximum+window) ){
					voteRoot[tree[i]][k]=1;
					votes++;
				}
			}
		}
	}
	if (g_reduce.collect){
		info.ncand = ncand;
		info.maximum = maximum;
		info.window = window;
		info.max_tree = max_tree;
		info.max_node = max_node;
		info.votes = votes;
		info.best_leaf = fired;
		last = info;
	}
	return maximum;
}
