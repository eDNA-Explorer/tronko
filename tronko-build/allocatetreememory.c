#include "allocatetreememory.h"

void allocatetreememory_for_nucleotide_Arr(int numberOfTrees){
  int i, j, k;
  for(k=0; k<numberOfTrees; k++){
  	for (i=0; i<(numspecArr[k]*2-1); i++){
    	int nb = numbaseArr[k];
    	/* Allocate pointer arrays */
    	treeArr[k][i].likenc = malloc(nb*(sizeof(double *)));
    	treeArr[k][i].posteriornc = malloc(nb*(sizeof(double *)));
    	/* Allocate contiguous data blocks (one malloc per array instead of per-site) */
    	double *like_data = malloc(nb * 4 * sizeof(double));
    	double *post_data = malloc(nb * 4 * sizeof(double));
    	/* Set up row pointers into contiguous blocks */
    	for (j=0; j<nb; j++){
      		treeArr[k][i].likenc[j] = &like_data[j*4];
      		treeArr[k][i].posteriornc[j] = &post_data[j*4];
    	}
  	}
  }
  for (i=0;i<4;i++){
    PMATnc[0][i][4]=1.0;//missing data
    PMATnc[1][i][4]=1.0;//missing data
  }
}
void allocateTreeArrMemory(struct masterArr *m, int max_nodename){
	int i;
	for (i=0; i<m->numspec*2-1; i++){
		m->tree[i].like = malloc(STATESPACE*(sizeof(double)));
		m->tree[i].posterior = malloc(STATESPACE*(sizeof(double)));
		m->tree[i].name = malloc((max_nodename+1)*sizeof(char));
		m->tree[i].up[0] = -2;  // -2 is uninitialized / NULL value
		m->tree[i].up[1] = -2;
		m->tree[i].down = -2;
		m->tree[i].nd = -2;
		m->tree[i].depth = -2;
		m->tree[i].bl = -2.0;
		m->tree[i].likenc = NULL;
		m->tree[i].posteriornc = NULL;
		m->tree[i].s = -2;
		m->tree[i].numsites = -2;
		m->tree[i].spec = -2;
		m->tree[i].mrca = -2;
		memset(m->tree[i].name, '\0', max_nodename+1);
	}
}
void freeTreeMemory(int whichPartition){
	int i,j,k;
	for(i=0; i<2*numspecArr[whichPartition]-1; i++){
		/* Free contiguous data blocks (first row pointer = base of block) */
		if (treeArr[whichPartition][i].likenc && numbaseArr[whichPartition] > 0)
			free(treeArr[whichPartition][i].likenc[0]);
		if (treeArr[whichPartition][i].posteriornc && numbaseArr[whichPartition] > 0)
			free(treeArr[whichPartition][i].posteriornc[0]);
		free(treeArr[whichPartition][i].like);
		free(treeArr[whichPartition][i].posterior);
		free(treeArr[whichPartition][i].likenc);
		free(treeArr[whichPartition][i].posteriornc);
		free(treeArr[whichPartition][i].name);
	}
	free(treeArr[whichPartition]);
}
void *calloc_check(size_t nmemb, size_t size){
	void *ptr = calloc(nmemb, size);
	if(!ptr){
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	return ptr;
}
void allocateMemoryForTaxArr(int whichPartitions, int max_tax_name){
	int phylogenyLevels = 7;
	int max_tax_name_len = max_tax_name;
	int i,j,k;
	taxonomyArr = (char ****)calloc_check(whichPartitions,sizeof(char ***));
	for(k=0;k<whichPartitions;k++){
		taxonomyArr[k] = (char ***)calloc_check(numspecArr[k], sizeof(char **));
		for(i=0; i<numspecArr[k]; i++){
			taxonomyArr[k][i] = (char **)calloc_check(phylogenyLevels, sizeof(char *));
			for(j=0; j<phylogenyLevels; j++){
				taxonomyArr[k][i][j] = (char *)calloc_check(max_tax_name_len, sizeof(char));
			}
		}
	}
}
