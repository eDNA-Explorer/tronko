#include "reduce.h"

type_of_PP reduce_votes(int ncand, const type_of_PP *const *sc, const int *tree, int **voteRoot)
{
	int i, k;
	type_of_PP maximum=-9999999999999999;
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if ( maximum < sc[i][k] ){
				maximum=sc[i][k];
			}
		}
	}
	for (i=0; i<ncand; i++){
		for (k=0; k<2*numspecArr[tree[i]]-1; k++){
			if ( sc[i][k] >= (maximum-Cinterval) && sc[i][k] <= (maximum+Cinterval) ){
				voteRoot[tree[i]][k]=1;
			}
		}
	}
	return maximum;
}
