#include "printtree.h"

/*old code for printing a tree*/
void printtreeArr(int whichRoot){
  int i,j;
  printf("\n");
  j=whichRoot;
  for (i=0; i<2*numspecArr[j]-1; i++){
    if (treeArr[j][i].up[0] != -1){
      printf("Node %i: up: (%i, %i) down: %i, name: %s, taxIndex[0]: %d, taxIndex[1]: %d, depth: %d, nd: %d, ",i,treeArr[j][i].up[0],treeArr[j][i].up[1],treeArr[j][i].down,treeArr[j][i].name,treeArr[j][i].taxIndex[0],treeArr[j][i].taxIndex[1],treeArr[j][i].depth,treeArr[j][i].nd);
	if ( treeArr[j][i].taxIndex[0] != -1 && treeArr[j][i].taxIndex[1] != -1 ){
		printf("taxonomy, %s",taxonomyArr[j][treeArr[j][i].taxIndex[0]][treeArr[j][i].taxIndex[1]]);
	}else{
		printf("taxonomy, Domain");
	}
    }else{
	printf("(Node %i): up: (%i, %i) down: %i, name: %s, taxIndex[0]: %d, taxIndex[1]: %d, depth: %d, nd: %d, ",i,treeArr[j][i].up[0],treeArr[j][i].up[1],treeArr[j][i].down,treeArr[j][i].name,treeArr[j][i].taxIndex[0],treeArr[j][i].taxIndex[1],treeArr[j][i].depth,treeArr[j][i].nd);
	if ( treeArr[j][i].taxIndex[0] != -1 && treeArr[j][i].taxIndex[1] != -1 ){
	       	printf("taxonomy, %s",taxonomyArr[j][treeArr[j][i].taxIndex[0]][treeArr[j][i].taxIndex[1]]);
	}else{
		printf("taxonomy, Domain");
	}
	}
    if (i != rootArr[j] ){
      printf(" (bl: %f)\n",treeArr[j][i].bl);
	}else {
		printf("\n");
	}
  }
}
void printTreeFile(int numberOfTrees, int max_nodename, int max_tax_name, int max_lineTaxonomy, Options opt){
	int i, j, k, l;
	char buf[BUFFER_SIZE];
	struct stat st = {0};
	if ( stat(opt.partitions_directory, &st) == -1){
		mkdir(opt.partitions_directory, 0700);
	}
	snprintf(buf,BUFFER_SIZE,"%s/reference_tree.txt",opt.partitions_directory);
	FILE *outputTree = fopen(buf,"w");
	if  ( outputTree == NULL ){ printf("Error opening reference tree file!\n"); exit(-1); }
	/* Opt 5: 256KB write buffer for bulk I/O */
	char *io_buf = malloc(256 * 1024);
	if (io_buf) setvbuf(outputTree, io_buf, _IOFBF, 256 * 1024);
	fprintf(outputTree,"%d\n",numberOfTrees);
	fprintf(outputTree,"%d\n",max_nodename);
	fprintf(outputTree,"%d\n",max_tax_name);
	fprintf(outputTree,"%d\n",max_lineTaxonomy);
	for (i=0; i<numberOfTrees;i++){
		fprintf(outputTree,"%d\t%d\t%d\n",numbaseArr[i],rootArr[i],numspecArr[i]);
	}
	for(i=0; i<numberOfTrees; i++){
		for(j=0; j<numspecArr[i]; j++){
			for(k=0; k<7; k++){
				if (k==6){
					fprintf(outputTree,"%s\n",taxonomyArr[i][j][k]);
				}else{
					fprintf(outputTree,"%s;",taxonomyArr[i][j][k]);
				}
			}
		}
	}
	/* Batch formatting: build lines in a local buffer to reduce fprintf call overhead.
	   Each posterior row is 4 doubles at ~24 chars each + tabs/newline ≈ 120 bytes max.
	   Use a 64KB buffer and flush when near-full. */
	{
		const int LBUF_SIZE = 65536;
		char *lbuf = malloc(LBUF_SIZE);
		if (!lbuf){ fprintf(stderr, "Out of memory for print buffer\n"); exit(1); }
		int lpos = 0;
		for(i=0; i<numberOfTrees;i++){
			int nb = numbaseArr[i];
			for(j=0;j<2*numspecArr[i]-1;j++){
				/* Node header line */
				if ( treeArr[i][j].up[0] != -1 && treeArr[i][j].up[1] != -1 ){
					lpos += snprintf(lbuf+lpos, LBUF_SIZE-lpos,
						"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t\n",
						i, j, treeArr[i][j].up[0], treeArr[i][j].up[1],
						treeArr[i][j].down, treeArr[i][j].depth,
						treeArr[i][j].taxIndex[0], treeArr[i][j].taxIndex[1]);
				}else{
					lpos += snprintf(lbuf+lpos, LBUF_SIZE-lpos,
						"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
						i, j, treeArr[i][j].up[0], treeArr[i][j].up[1],
						treeArr[i][j].down, treeArr[i][j].depth,
						treeArr[i][j].taxIndex[0], treeArr[i][j].taxIndex[1],
						treeArr[i][j].name);
				}
				if (lpos > LBUF_SIZE - 512){ fwrite(lbuf, 1, lpos, outputTree); lpos = 0; }
				/* Posterior probability rows */
				double **post = treeArr[i][j].posteriornc;
				for (k=0; k<nb; k++){
					lpos += snprintf(lbuf+lpos, LBUF_SIZE-lpos,
						"%.17g\t%.17g\t%.17g\t%.17g\n",
						post[k][0], post[k][1], post[k][2], post[k][3]);
					if (lpos > LBUF_SIZE - 512){ fwrite(lbuf, 1, lpos, outputTree); lpos = 0; }
				}
			}
		}
		if (lpos > 0) fwrite(lbuf, 1, lpos, outputTree);
		free(lbuf);
	}
	fclose(outputTree);
	free(io_buf);
}
void printTaxonomyArrToFile(int numberOfTrees){
	int i,j,k;
	FILE *outputTaxArr = fopen("/space/s2/lenore/partitions2/taxonomy_reference.txt","w");
	if ( outputTaxArr == NULL ){ printf("Error opening file!\n"); exit(1); }
	for(i=0; i<numberOfTrees;i++){
		fprintf(outputTaxArr,"TREE %d\n",i);
		for(j=0;j<numspecArr[i];j++){
			fprintf(outputTaxArr,"%d",j);
			for(k=0;k<7;k++){
				fprintf(outputTaxArr,"\t%s",taxonomyArr[i][j][k]);
			}
			fprintf(outputTaxArr,"\n");
		}
	}
	fclose(outputTaxArr);
}

void printPP(){
	int i,j;
	FILE *outputPP;
	outputPP=fopen("/space/s2/lenore/partitions2/initial_tree_PP.txt","w");
	if (outputPP==NULL){ printf("Error opening file!\n"); exit(1); }
	printf("numspec is %d\n",numspec);
	for(i=0; i<2*numspec-1; i++){
		printf("NODE\t%d\t%d\n",i,numbase);
		//fprintf(stderr,"NODE\t%d\n",i);
		for(j=0;j<numbase;j++){
			printf("%d\t%lf\t%lf\t%lf\t%lf\n",j,PP[i][j][0],PP[i][j][1],PP[i][j][2],PP[i][j][3]);
			fprintf(outputPP,"%d\t%lf\t%lf\t%lf\t%lf\n",j,PP[i][j][0],PP[i][j][1],PP[i][j][2],PP[i][j][3]);
		}
	}
	fclose(outputPP);
}
void printPP_Arr(int numberOfTrees){
	int i,j,k;
	FILE *outputPP;
	outputPP=fopen("/space/s2/lenore/partitions2/PP_Arr.txt","w");
	if (outputPP==NULL){ printf("Error opening file!\n"); exit(1); }
	for(i=0; i<numberOfTrees; i++){
		fprintf(outputPP,"TREE\t%d\n",i);
		for(j=0; j<2*numspecArr[i]-1; j++){
			fprintf(outputPP,"NODE\t%d\n",j);
			for(k=0; k<numbaseArr[i]; k++){
				fprintf(outputPP,"%d\t%lf\t%lf\t%lf\t%lf\n",k,PP_Arr[i][j][k][0],PP_Arr[i][j][k][1],PP_Arr[i][j][k][2],PP_Arr[i][j][k][3]);
			}
		}
	}
	fclose(outputPP);
}
