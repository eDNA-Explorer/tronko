#include "reduce.h"

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

/* L1, the best leaf (ties to the earlier candidate, then the lower node), and L2, the best
 * leaf naming another species: what the trace reports and the best-leaf rule tests. */
static void best_leaves(int ncand, const type_of_PP *const *sc, const int *tree, reduce_info_t *info)
{
	int i, k;
	info->l1_tree = info->l1_node = info->l2_tree = info->l2_node = -1;
	info->l1 = info->l2 = -INFINITY;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if (is_leaf(tree[i], k) && sc[i][k] > info->l1){
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
			if (is_leaf(tree[i], k) && sc[i][k] > info->l2
			    && strcmp(reduce_leaf_species(tree[i], k), species) != 0){
				info->l2 = sc[i][k];
				info->l2_tree = tree[i];
				info->l2_node = k;
			}
		}
	}
}

type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int **voteRoot)
{
	int i, k;
	type_of_PP maximum=-9999999999999999;
	int max_tree=-1, max_node=-1;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if ( maximum < sc[i][k] ){
				maximum=sc[i][k];
				max_tree=tree[i];
				max_node=k;
			}
		}
	}
	int votes=0;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if ( sc[i][k] >= (maximum-Cinterval) && sc[i][k] <= (maximum+Cinterval) ){
				voteRoot[tree[i]][k]=1;
				votes++;
			}
		}
	}
	if (g_reduce.collect){
		last.ncand = ncand;
		last.maximum = maximum;
		last.window = Cinterval;
		last.max_tree = max_tree;
		last.max_node = max_node;
		last.votes = votes;
		best_leaves(ncand, sc, tree, &last);
	}
	return maximum;
}
