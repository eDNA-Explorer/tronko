#include "assignment.h"

/*type_of_PP **assignScores_Arr(int rootNum, int node, char *locQuery, int *positions, type_of_PP **scores, int alength){

	int child0 = treeArr[rootNum][node].up[0];
	int child1 = treeArr[rootNum][node].up[1];
	if (child0 == -1 && child1 == -1){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//printf("node %d, score %lf, root %d\n",node,tree[node].score[rootNum],nodesToCutMinVar[rootNum]);
		return scores;
	}else if (child0 != -1 && child1 != -1){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum);
		assignScores_Arr(rootNum, child0, locQuery, positions, scores, alength);
		assignScores_Arr(rootNum, child1, locQuery, positions, scores, alength);
	}*/
	/*if (child1 != -1 ){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum);
		assignScores_Arr(rootNum, child1, locQuery, positions, scores, alength);
	}*/
/*}*/
void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
    FILE* site_scores_file, char* readname,
    int early_termination, type_of_PP *best_score, int *strikes,
    type_of_PP strike_box, int max_strikes,
    int enable_pruning, type_of_PP pruning_threshold) {

	/*
	 * Iterative DFS traversal using an explicit stack.
	 * Replaces the recursive version to avoid ~2931 stack frames with 17 parameters each.
	 * The stack is allocated on the C stack (VLA) since max tree size is bounded.
	 */
	int stack[4096];  /* Max tree depth; binary tree with ~2931 nodes has depth <= ~12 */
	int sp = 0;
	stack[sp++] = node;

	type_of_PP *scores_arr = scores[search_number][rootNum];

	while (sp > 0) {
		int cur = stack[--sp];
		int child0 = treeArr[rootNum][cur].up[0];
		int child1 = treeArr[rootNum][cur].up[1];

		type_of_PP node_score = getscore_Arr(alength, cur, rootNum, locQuery, positions,
		                                      print_all_nodes, site_scores_file, readname);
		scores_arr[cur] += node_score;

		if (child0 == -1 && child1 == -1) {
			/* Leaf node — update best score tracking */
			if (early_termination && node_score > *best_score) {
				*best_score = node_score;
				*strikes = 0;
			}
		} else if (child0 != -1 && child1 != -1) {
			/* Internal node — pruning check, then push children */
			if (enable_pruning && *best_score > -9999999999999998) {
				if (node_score < *best_score - pruning_threshold) {
					continue;  /* Skip children entirely */
				}
			}
			/* Push child1 first so child0 is processed first (preserves DFS order) */
			stack[sp++] = child1;
			stack[sp++] = child0;
		}
	}
}
/*type_of_PP getscore(int alength, int node, int rootNum){
	type_of_PP score;
	int pos, i;
	score=0;
	if (when==2){
	  char* level;
	  level=malloc(sizeof(char)*20);
	  if (tree[node].taxIndex[1]==0){
	  level = "species";
	  }
	  if (tree[node].taxIndex[1]==1){
	  level = "genus";
	  }
	  if (tree[node].taxIndex[1]==2){
	  level="family";
	  }
	  if (tree[node].taxIndex[1]==3){
	  level="order";
	  }
	  if (tree[node].taxIndex[1]==4){
	  level="class";
	  }
	  if (tree[node].taxIndex[1]==5){
	  level="phylum";
	  }
	  printf("LCA Node %i, taxonomy (%s): %s, alength: %i %c\n",node,level,taxonomy[tree[node].taxIndex[0]][tree[node].taxIndex[1]],alength,locQuery[alength]);
	  }
	  if (alength < 40 ){
		  return 999999999999999;
	  }
	  int j=0;
	for (i=0; i<numbase; i++){//We assume that we have already ensured that the sequence only contains a, c, t, and g.  Missing data is not aligned
		pos=positions[j];//test if this is faster rather than referencing multiple times.  Compiler optimization may make the latter faster.
		if (i==pos){
			if (locQuery[j]=='-' && tree[node].taxIndex[1]==0 && node==rootNum){
				score=score+1;
			}else{
				//might be better simply to delete the parts of the query that does not align. This probably slows doww a lot.
				if (locQuery[j]=='a') score += PP[node][pos][0]; //test if a swtich statement is faster in this case. Might be compiler dependent.
				else if (locQuery[j]=='c') score += PP[node][pos][1];
				else if (locQuery[j]=='g') score += PP[node][pos][2];
				else if (locQuery[j]=='t') score += PP[node][pos][3];
				else score = score+1;
				//oldscore=score;
			}
		//}
			j++;
		}else{
			score = score+1;
		}
	}
	return score;
}*/
int checkPolyA(int rootNum, int node, int position){
	int i,j;
	int isPolyA = 0;
	for(i=0; i<position; i++){
		if ( treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1 ){
			isPolyA = 1;
		}else{
			isPolyA = 0;
			break;
		}
	}
	if ( isPolyA == 1 ){
		return 1;
	}
	for(i=numbaseArr[rootNum]-1; i>=position; i--){
		if ( treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1 ){
			isPolyA = 1;
		}else{
			isPolyA = 0;
			break;
		}
	}
	return isPolyA;
}
/*
 * Precomputed log constants and nucleotide lookup table.
 * Avoids recomputing log() on every iteration of the inner scoring loop.
 */
static const type_of_PP LOG_001 = -4.605170185988091;  /* log(0.01) */
static const type_of_PP LOG_025 = -1.3862943611198906; /* log(0.25) */

/*
 * Nucleotide character -> index lookup table.
 * A/a=0, C/c=1, G/g=2, T/t=3, '-'=4 (gap), everything else=5 (unknown).
 * Eliminates branch-heavy if/else chain in the inner loop.
 */
static const int NUC_LUT[256] = {
	[0 ... 255] = 5,  /* default: unknown */
	['A'] = 0, ['a'] = 0,
	['C'] = 1, ['c'] = 1,
	['G'] = 2, ['g'] = 2,
	['T'] = 3, ['t'] = 3,
	['-'] = 4
};

type_of_PP getscore_Arr(int alength, int node, int rootNum, char *locQuery, int *positions, int print_all_nodes, FILE* site_scores_file, char* readname){
	type_of_PP score;
	int i;
	score=0;
	if (positions[0]==-1){
		score=9999999999;
		return score;
	}

	/* Fast path: no debug output (the overwhelmingly common case) */
	if (print_all_nodes != 1) {
		const type_of_PP *pp = treeArr[rootNum][node].posteriornc;
		for (i=0; i<alength; i++){
			int pos = positions[i];
			if (pos==-1){
				score += LOG_001;
			}else{
				if (pp[PP_IDX(pos, 0)] == 1){
					/* Missing data position */
					if (locQuery[i] != '-'){
						score += LOG_001;
					}
					/* else: gap at missing data = +0, nothing to add */
				}else{
					int nuc = NUC_LUT[(unsigned char)locQuery[i]];
					if (nuc <= 3){
						score += pp[PP_IDX(pos, nuc)];
					}else if (nuc == 4){
						/* gap character */
						score += LOG_025;
					}
					/* nuc==5 (unknown): no contribution, same as original */
				}
			}
		}
		return score;
	}

	/* Slow path: with debug output (print_all_nodes == 1) */
	for (i=0; i<alength; i++){
		fprintf(site_scores_file,"%s\t%d\t%d\t%d\t",readname,rootNum,node,positions[i]);
		if (positions[i]==-1){
			score += LOG_001;
			fprintf(site_scores_file,"%lf\n",LOG_001);
		}else{
			if ( treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)]==1 && locQuery[i]=='-' ){
				score=score+0;
				fprintf(site_scores_file,"%lf\n",0.0);
			}else{
				if( treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)] == 1){
					score += LOG_001;
					fprintf(site_scores_file,"%lf\n",LOG_001);
				}else{
					int nuc = NUC_LUT[(unsigned char)locQuery[i]];
					if (nuc <= 3){
						type_of_PP pp_val = treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], nuc)];
						score += pp_val;
						fprintf(site_scores_file, PP_PRINT_FORMAT "\n", pp_val);
					}else if (nuc == 4){
						score += LOG_025;
						fprintf(site_scores_file,"%lf\n",LOG_025);
					}
				}
			}
		}
	}

	return score;
}
