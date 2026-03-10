#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/wait.h>
#include "global.h"
#include "opt.h"
#include "options.h"
#include "hashmap.h"
#include "allocatetreememory.h"
#include "printtree.h"
#include "likelihood.h"
#include "getclade.h"
#include "readfasta.h"
#include "readreference.h"

HASHMAP(int, struct masterArr) mastermap;
struct node **treeArr;
int numspec, numbase, **seq, numundspec[MAXNUMBEROFINDINSPECIES+1];
char **nodeIDs;
char ***nodeIDsArr;
__thread int root,tip,comma=0; /*globals used to read in the tree*/
double Logfactorial[MAXNUMBEROFINDINSPECIES];
double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
double *statevector, UFC, *UFCnc, **templike_nc;
int COUNT2;
int COUNT;
double *localpi;
int localnode;
double parameters[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
type_of_PP ***PP;
type_of_PP ***PPcopy;
type_of_PP ****PP_Arr;
char ****taxonomyArr;
int *SPscoreArr;
int *numspecArr, *numbaseArr, ***seqArr, *rootArr;
double minVariance = 99999999999999;
int minVarNode = -1;
int returnNode=-1;
int *partitionSizes;
int *nodesToCut;
int *nodesToCutMinVar;
node **treeArr;

static pthread_mutex_t mastermap_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t spscore_mutex = PTHREAD_MUTEX_INITIALIZER;

/*This function calculates the number of  desdencents of each node in the tree stored in tree[node].nd*/
int get_number_descendantsArr(int node, struct masterArr *m){
	if (node==-1) return 0;
	if (m->tree[node].up[0]==-1) return (m->tree[node].nd=1);
	else return (m->tree[node].nd=(get_number_descendantsArr(m->tree[node].up[0],m)+get_number_descendantsArr(m->tree[node].up[1],m)));
}
void changePP_Arr (int node, int whichRoot){
	int child0 = treeArr[whichRoot][node].up[0];
	int child1 = treeArr[whichRoot][node].up[1];
	int i,j;
	if ( child0 < 0 || child1 < 0 ){
		/* Leaf node or degenerate node with missing child — mark gaps */
		for (i=0; i<numbaseArr[whichRoot]; i++){
			if ( seqArr[whichRoot][node-numspecArr[whichRoot]+1][i] == 4 ){
				for(j=0;j<4;j++){
					treeArr[whichRoot][node].posteriornc[i][j]=-1;
				}
			}
		}
		return;
	}else{
		changePP_Arr(child0,whichRoot);
		changePP_Arr(child1,whichRoot);
		return;
	}
}
void changePP_parents_Arr(int node, int whichRoot){
	int maxNodes = 2 * numspecArr[whichRoot] - 1;
	int iter = 0;
	while (iter < maxNodes) {
		int parent = treeArr[whichRoot][node].down;
		if ( parent < 0 ){
			return;
		}
		int child0 = treeArr[whichRoot][node].up[0];
		int child1 = treeArr[whichRoot][node].up[1];
		if ( child0 < 0 || child1 < 0 ){
			node = parent;
			iter++;
			continue;
		}
		int i,j;
		for( i=0; i<numbaseArr[whichRoot]; i++){
			for (j=0; j<4; j++){
				if ( treeArr[whichRoot][child0].posteriornc[i][j] == -1 && treeArr[whichRoot][child1].posteriornc[i][j] == -1 ){
					treeArr[whichRoot][node].posteriornc[i][j]=-1;
				}
			}
		}
		node = parent;
		iter++;
	}
}
void assignDepthArr(int node0, int node1, int depth, struct masterArr *m){
	if( node0 != -1 && node1 != -1){
		m->tree[node0].depth = depth;
		m->tree[node1].depth = depth;
		assignDepthArr(m->tree[node0].up[0], m->tree[node0].up[1],depth+1,m);
		assignDepthArr(m->tree[node1].up[0], m->tree[node1].up[1],depth+1,m);
	}
}
void findMaxTaxName(FILE* file, int* specs){
	int max_tax_name=0;
	char buffer[BUFFER_SIZE];
	char *s;
	char *lineTaxonomy;
	int max_line = 0;
	while( fgets(buffer,BUFFER_SIZE,file) != NULL ){
		s = strtok(buffer,"\t");
		lineTaxonomy = strtok(NULL,"\n");
		if (lineTaxonomy == NULL) continue;
		int taxsize = strlen(lineTaxonomy);
		if ( max_line < taxsize ){
			max_line = taxsize;
		}
		int j=6;
		while(j>-1){
			s = strtok_r(lineTaxonomy,";",&lineTaxonomy);
			if (s == NULL) break;
			int size = strlen(s);
			if ( max_tax_name < size ){
				max_tax_name = size;
			}
			j--;
		}
	}
	if ( specs[0] < max_tax_name ){
		specs[0] = max_tax_name;
	}
	if (specs[1] < max_line ){
		specs[1] = max_line;
	}
}
void assignTaxonomyToLeavesArr(char *tax,struct masterArr *m, int max_nodename, int max_tax_name){
	char buffer[BUFFER_SIZE];
	char *lineAccession, *lineTaxonomy, *taxLevelName;
	int i, j;
	FILE *taxonomy_file;

	/* Set internal nodes' taxIndex to -1 */
	for(i=0; i<m->numspec - 1; i++){
		m->tree[i].taxIndex[0] = -1;
		m->tree[i].taxIndex[1] = -1;
	}

	/* Opt 2: Build name->leaf-index hashmap for O(1) lookup */
	HASHMAP(char, int) name_map;
	hashmap_init(&name_map, hashmap_hash_string, strcmp);
	hashmap_set_key_alloc_funcs(&name_map, strdup, free);

	for(i = m->numspec - 1; i < 2*m->numspec - 1; i++){
		int *idx = malloc(sizeof(int));
		*idx = i;
		hashmap_put(&name_map, m->tree[i].name, idx);
	}

	/* Read taxonomy file ONCE (was re-opened N times before) */
	if ((taxonomy_file = fopen(tax,"r")) == (FILE *) NULL){
		fprintf(stderr, "*** taxonomy file could not be opened: %s\n", tax);
		exit(-1);
	}
	while(fgets(buffer, BUFFER_SIZE, taxonomy_file) != NULL){
		lineAccession = strtok(buffer, "\t");
		lineTaxonomy = strtok(NULL, "\n");
		if (lineAccession == NULL || lineTaxonomy == NULL) continue;

		int *leaf_ptr = hashmap_get(&name_map, lineAccession);
		if (leaf_ptr != NULL){
			i = *leaf_ptr;
			m->tree[i].taxIndex[0] = i - m->numspec + 1;
			m->tree[i].taxIndex[1] = 0;
			j = 6;
			while(j > -1){
				taxLevelName = strtok_r(lineTaxonomy, ";", &lineTaxonomy);
				if (taxLevelName == NULL) break;
				strncpy(m->taxonomy[i - m->numspec + 1][j], taxLevelName, max_tax_name);
				j--;
			}
		}
	}
	fclose(taxonomy_file);

	/* Cleanup hashmap */
	const char *hkey;
	int *hdata;
	hashmap_foreach(hkey, hdata, &name_map){
		free(hdata);
	}
	hashmap_cleanup(&name_map);
	/*if ( child0 == -1 && child1 == -1 ){
		if (( taxonomy_file = fopen(tax,"r")) == (FILE *) NULL) printf("*** taxonomy file could not be opened.\n");
		while( fgets(buffer,BUFFER_SIZE,taxonomy_file) != NULL){
			lineAccession = strtok(buffer,"\t");
			lineTaxonomy = strtok(NULL,"\n");
			assert(strlen(lineAccession) <= max_nodename);
			strncpy(&accessionID[0], lineAccession, max_nodename);
			if ( strcmp(accessionID,m->tree[node].name)==0 ){
				m->tree[node].taxIndex[0]=node-m->numspec+1;
				m->tree[node].taxIndex[1]=0;
				int j=6;
				while(j>-1){
					taxLevelName = strtok_r(lineTaxonomy,";",&lineTaxonomy);
					strncpy(m->taxonomy[node-m->numspec+1][j], taxLevelName, max_tax_name); 
					j--;
				}
			}
		}
		fclose(taxonomy_file);
		return;
	}else{
		m->tree[node].taxIndex[0] = -1;
		m->tree[node].taxIndex[1] = -1;
	}
	if(child0 != -1){
		assignTaxonomyToLeavesArr(child0,tax,m,max_nodename,max_tax_name);
	}
	if(child1 != -1){
		assignTaxonomyToLeavesArr(child1,tax,m,max_nodename,max_tax_name);
	}*/
}
/* Opt 3: Stack-allocated recursion — void + output parameter instead of malloc per call */
void getTaxonomyArr(int node, struct masterArr *m, int *out){
	int taxIndexA[2], taxIndexB[2];
	out[0] = -1;
	out[1] = -1;
	if ( m->tree[node].up[0] != -1 ){
		if ( m->tree[node].taxIndex[0] == -1 ){
			getTaxonomyArr(m->tree[node].up[0],m,taxIndexA);
			getTaxonomyArr(m->tree[node].up[1],m,taxIndexB);
			if ( taxIndexA[0]==-1 || taxIndexB[0]==-1 ){
				m->tree[node].taxIndex[0]=-1;
				m->tree[node].taxIndex[1]=-1;
			}else{
				int i=0;
				int phylogenyLevel=0;
				int maxABLevel = (taxIndexA[1] > taxIndexB[1]) ? taxIndexA[1] : taxIndexB[1];
				for(i=maxABLevel;i<7;i++){
					if ( strcmp(m->taxonomy[taxIndexA[0]][i],m->taxonomy[taxIndexB[0]][i])==0){
						phylogenyLevel = i;
						out[0] = taxIndexA[0];
						out[1] = phylogenyLevel;
						m->tree[node].taxIndex[0] = taxIndexA[0];
						m->tree[node].taxIndex[1] = phylogenyLevel;
						break;
					}
				}
			}
		}
	}else{
		out[0] = m->tree[node].taxIndex[0];
		out[1] = m->tree[node].taxIndex[1];
	}
}
int* getTaxonomyArr_UsePartitions(int node, int whichPartitions, char**** taxonomyArr_heap){
	int* taxIndexA = NULL;
	int* taxIndexB = NULL;
	int* taxonomyOfNode = malloc(sizeof(int)*2);
	taxonomyOfNode[0] = -1;
	taxonomyOfNode[1] = -1;
	if ( treeArr[whichPartitions][node].up[0] != -1 ){
		if ( treeArr[whichPartitions][node].taxIndex[0] == -1 ){
			taxIndexA = getTaxonomyArr_UsePartitions(treeArr[whichPartitions][node].up[0],whichPartitions,taxonomyArr_heap);
			taxIndexB = getTaxonomyArr_UsePartitions(treeArr[whichPartitions][node].up[1],whichPartitions,taxonomyArr_heap);
			if ( taxIndexA[0]==-1 || taxIndexB[0]==-1 ){
				treeArr[whichPartitions][node].taxIndex[0]=-1;
				treeArr[whichPartitions][node].taxIndex[1]=-1;
			}else{
				int i=0;
				int phylogenyLevel=0;
				int maxABLevel = (taxIndexA[1] > taxIndexB[1]) ? taxIndexA[1] : taxIndexB[1];
				//int maxABLevel = 0;
				for(i=maxABLevel;i<7;i++){
					if ( strcmp(taxonomyArr_heap[whichPartitions][taxIndexA[0]][i],taxonomyArr_heap[whichPartitions][taxIndexB[0]][i])==0){
						phylogenyLevel = i;
						taxonomyOfNode[0] = taxIndexA[0];
						taxonomyOfNode[1] = phylogenyLevel;
						treeArr[whichPartitions][node].taxIndex[0] = taxIndexA[0];
						treeArr[whichPartitions][node].taxIndex[1] = phylogenyLevel;
						break;
					}
				}
			}
		}
	}else{
		taxonomyOfNode[0] = treeArr[whichPartitions][node].taxIndex[0];
		taxonomyOfNode[1] = treeArr[whichPartitions][node].taxIndex[1];
	}
	if(taxIndexA != NULL) free(taxIndexA);
	if(taxIndexB != NULL) free(taxIndexB);
	return taxonomyOfNode;
}
void readSeqArr(gzFile partitionsFile, int maxname, struct masterArr *master){
	char buffer[FASTA_MAXLINE];
	int i, j, m, k=0, row=0;
	char c;
	int size;
	char nodename[maxname+1];
	char *s;
	int firstIter=1;
	int index=0;
	while( gzgets(partitionsFile,buffer,FASTA_MAXLINE) != NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( buffer[0] == '>'){
			if ( size > maxname ){ size = maxname; }
			for(i=1; buffer[i]!='\0' && (i-1) < maxname; i++){
				nodename[i-1]=buffer[i];
			}
			nodename[size-1]='\0';
			strcpy(master->names[index], nodename);
			index++;
		}else{
			if ( firstIter==1 ){ 
				int l=0;
				int length=0;
				for(l=0; buffer[l] != '\0'; l++){
					length++;
				}
				master->numbase=length;
				master->msa = (int**)malloc(master->numspec*sizeof(int*));
				for (i=0; i<master->numspec; i++){
					master->msa[i]=(int *)malloc(master->numbase*(sizeof(int)));
				}
				firstIter=0;
				i=0;
			}
			for( m=0; m<master->numbase; m++){
				c=tolower(buffer[m]);
				if ((c!=' ')&&(c!='\t')) {
					if (c=='a' || c=='A') master->msa[row][m] = 0;
					else if (c=='c' || c=='C') master->msa[row][m] = 1;
					else if (c=='g' || c=='G') master->msa[row][m] = 2;
					else if (c=='t' || c=='T') master->msa[row][m] = 3;
					else if (c=='n'||c=='-'||c=='~' || c=='N') master->msa[row][m] = 4;
					else{
						/* IUPAC ambiguity codes (Y,R,W,S,K,M,B,D,H,V,etc.) → treat as missing data */
						master->msa[row][m] = 4;
					}
				}
			}
			row++;
		}
	}
	if (master->numspec < 3) {printf("This is for more than two seq.s only!\n"); exit(-1);}
}
char* readNewickFile( FILE *file){
	fseek(file, 0, SEEK_END);        // Move to the end of the file
	long fileSize = ftell(file);     // Get the size of the file
	fseek(file, 0, SEEK_SET);        // Move back to the beginning
	char *buffer = malloc(fileSize + 1); // +1 for null terminator
	if (buffer == NULL) {
		fprintf(stderr, "Error: Memory allocation failed while reading Newick file.\n");
		exit(EXIT_FAILURE);
	}
	size_t bytesRead = fread(buffer, 1, fileSize, file);
	if (bytesRead != fileSize) {
		fprintf(stderr, "Error: Could not read the entire Newick file.\n");
		free(buffer);
		exit(EXIT_FAILURE);
	}
	buffer[bytesRead] = '\0';
	return buffer;
}
int *findLeavesOfMinVarArr(int node, int *leafNodeList, int size, struct masterArr *m, int local_minVarNode){
	if (node==-1){ return leafNodeList; }
	int child0 = m->tree[node].up[0];
	int child1 = m->tree[node].up[1];
	int i=0;
	int currentPos=-1;
	if ( node == local_minVarNode ){ return leafNodeList; }
	if (child0 ==-1 && child1 == -1){
		for(i=size-1;i>=0;i--){
			if (leafNodeList[i] == -1){
				currentPos = i;
			}
		}
		leafNodeList[currentPos]=node;
		return leafNodeList;
	}
	findLeavesOfMinVarArr(child0,leafNodeList,size,m,local_minVarNode);
	findLeavesOfMinVarArr(child1,leafNodeList,size,m,local_minVarNode);
	return leafNodeList;
}
char *removeDashArr(int *sequenceWithDash,int length, char *sequenceWithoutDash){
	int i,j;
	int countDash;
	countDash=0;
	for(i=0;i<length;i++){
		if( sequenceWithDash[i] < 4){
			if ( sequenceWithDash[i]==0 ){ sequenceWithoutDash[countDash]='A'; }
				else if (sequenceWithDash[i]==1 ){ sequenceWithoutDash[countDash]='C'; }
				else if (sequenceWithDash[i]==2 ){ sequenceWithoutDash[countDash]='G'; }
				else { sequenceWithoutDash[countDash]='T'; }
			countDash++;
		}
	}
	sequenceWithoutDash[countDash]='\0';
	return sequenceWithoutDash;
}
void printPartitionsToFileArr(int *partition1,int partition1size, int *partition2,int partition2size, int *partition3, int partition3size,int whichPartitions, int partitionCount, Options opt, struct masterArr *m){
	int i=0;
	int j=0;
	char buf_fasta[BUFFER_SIZE];
	char buf_tax[BUFFER_SIZE];
	int which = partitionCount;
	char *seqWithoutDash;
	struct stat st = {0};
	if ( stat(opt.partitions_directory, &st) == -1){
		mkdir(opt.partitions_directory, 0700);
	}
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);
	}else{
		snprintf(buf_fasta,BUFFER_SIZE,"%s/%spartition%d.fasta",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p1=fopen(buf_fasta,"w");
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	}else{
		snprintf(buf_tax,BUFFER_SIZE,"%s/%spartition%d_taxonomy.txt",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p1_tax=fopen(buf_tax,"w");
	seqWithoutDash = (char *)malloc((m->numbase+1)*sizeof(char));
	for(i=0; i<m->numbase+1; i++){
		seqWithoutDash[i]='\0';
	}
	for(i=0;i<partition1size;i++){
		fprintf(p1,">%s\n",m->names[partition1[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition1[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p1,"%s\n",seqWithoutDash);
		int taxnode=0;
		fprintf(p1_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition1[i]-m->numspec+1],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][0]);
		for (j=0; j<m->numbase; j++){
			seqWithoutDash[j]='\0';
		}
	}
	fclose(p1);
	fclose(p1_tax);
	which = partitionCount+1;
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);				
	}else{
		snprintf(buf_fasta,BUFFER_SIZE,"%s/%spartition%d.fasta",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p2=fopen(buf_fasta,"w");
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	}else{
		snprintf(buf_tax,BUFFER_SIZE,"%s/%spartition%d_taxonomy.txt",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p2_tax=fopen(buf_tax,"w");
	for(i=0;i<partition2size;i++){
		fprintf(p2,">%s\n",m->names[partition2[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition2[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p2,"%s\n",seqWithoutDash);
		fprintf(p2_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition2[i]-m->numspec+1],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][0]);
		for(j=0; j<m->numbase+1; j++){
			seqWithoutDash[j]='\0';
		}
	}
	fclose(p2);
	fclose(p2_tax);
	which=partitionCount+2; 
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);				
	}else{
		 snprintf(buf_fasta,BUFFER_SIZE,"%s/%spartition%d.fasta",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p3=fopen(buf_fasta,"w");
	if ( opt.prefix[0] == '\0' ){
		snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	}else{
		snprintf(buf_tax,BUFFER_SIZE,"%s/%spartition%d_taxonomy.txt",opt.partitions_directory,opt.prefix,which);
	}
	FILE *p3_tax=fopen(buf_tax,"w");
	for(i=0;i<partition3size;i++){
		fprintf(p3,">%s\n",m->names[partition3[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition3[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p3,"%s\n",seqWithoutDash);
		fprintf(p3_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition3[i]-m->numspec+1],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][0]);
		for(j=0; j<m->numbase; j++){
			seqWithoutDash[j]='\0';
		}
	}
	free(seqWithoutDash);
	fclose(p3);
	fclose(p3_tax);
}
double calculateSPArr(struct masterArr *m){
	int i,j,k;
	int numpairs=0;
	double SPscore=0;
	j=0;
	int *partition = (int*)malloc(sizeof(int)*m->numspec);
	for(i=m->numspec-1; i<2*m->numspec-1; i++){
		partition[j]=i;
		j++;
	}	
	for (i = 0; i<m->numbase; i++){
		for( j=0; j<m->numspec; j++){
			for( k=j+1; k<m->numspec; k++){
				if (m->msa[partition[j]-m->numspec+1][i]==m->msa[partition[k]-m->numspec+1][i]){
					SPscore=SPscore+3;
				}else if (m->msa[partition[j]-m->numspec+1][i] != 4 && m->msa[partition[k]-m->numspec+1][i] != 4 && m->msa[partition[j]-m->numspec+1][i]!=m->msa[partition[k]-m->numspec+1][i] ){
					SPscore=SPscore-2;
				}else{
					SPscore--;
				}
				numpairs++;
			}
		}
	}
	free(partition);
	//printf("raw SPscore: %lf\n",SPscore);
	//printf("numpairs: %d\n",numpairs);
	double exponent;
	exponent = (-5.0/2.0);
	//SPscore=SPscore/(m->numspec*exp(exponent));
	//SPscore=SPscore/(m->numspec*exp(-5/2));
	SPscore=SPscore/m->numspec;
	SPscore=SPscore/numpairs;
	printf("SPscore: %lf\n",SPscore);
	return SPscore;
}
void writeNewick(FILE* file, masterArr* m, int nodeIndex) {
	node* current = &m->tree[nodeIndex];
	if (current->up[0] == -1) {
		// Leaf node
		if (current->name != NULL) {
			fprintf(file, "%s", current->name);
		}else{
			fprintf(file, "UnnamedLeaf");
		}
	}else{
		// Internal node
		fprintf(file, "(");
		for (int i = 0; i < current->nd; i++) {
			if ( i>0 ){
				fprintf(file, ",");
			}
			writeNewick(file, m, current->up[i]);
		}
		fprintf(file, ")");

		// Optional: Output the internal node's name if it has one
		//if (current->name != NULL) {
		//	fprintf(file, "%s", current->name);
		//}
	}
		fprintf(file, ":%f", current->bl);
}
void exportTreeToNewick(struct masterArr* m, const char* filename){
	FILE* file = fopen(filename, "w");
	if ( file == NULL ){
		fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
		exit(EXIT_FAILURE);
	}
	writeNewick(file, m, m->root);
	fprintf(file, ";\n"); // Newick format ends with a semicolon
	fclose(file);
}
void findMinVarianceArr(int node, int size, struct masterArr *m, double *out_minVariance, int *out_minVarNode){
	if (node==-1){ return; }
	int child0 = m->tree[node].up[0];
	int child1 = m->tree[node].up[1];
	int parent = m->tree[node].down;
	int i=0;
	if ( parent == -1 ){
		findMinVarianceArr(child0,size,m,out_minVariance,out_minVarNode);
		findMinVarianceArr(child1,size,m,out_minVariance,out_minVarNode);
		return;
	}
	if (child0 == -1 ){ return; }
	if (child1 == -1 ){ return; }
	int num_children0 = m->tree[child0].nd;
	int num_children1 = m->tree[child1].nd;
	int num_ancestors = size-num_children1-num_children0;
	double mean = (double)(num_children0 + num_children1 + num_ancestors )/3;
	double variance = (double)((num_ancestors-mean)*(num_ancestors-mean) + (num_children0-mean)*(num_children0-mean) + (num_children1-mean)*(num_children1-mean))/3;
	if (*out_minVariance > variance){
		*out_minVariance = variance;
		*out_minVarNode = node;
	}
	findMinVarianceArr(child0,size,m,out_minVariance,out_minVarNode);
	findMinVarianceArr(child1,size,m,out_minVariance,out_minVarNode);
	return;
}
/* Unwrap multi-line FASTA sequences in-place (replaces gsed call) */
static void unwrap_fasta_inplace(const char *filepath){
	FILE *fp = fopen(filepath, "r");
	if (!fp){ fprintf(stderr, "unwrap_fasta_inplace: cannot open %s\n", filepath); return; }
	/* Read entire file into memory */
	fseek(fp, 0, SEEK_END);
	long fsize = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	char *data = (char *)malloc(fsize + 1);
	if (!data){ fclose(fp); return; }
	fsize = fread(data, 1, fsize, fp);
	data[fsize] = '\0';
	fclose(fp);
	/* Write back with continuation lines joined */
	fp = fopen(filepath, "w");
	if (!fp){ free(data); return; }
	int in_header = 0;
	long i;
	for (i = 0; i < fsize; i++){
		if (data[i] == '>'){
			in_header = 1;
			fputc(data[i], fp);
		} else if (data[i] == '\n'){
			if (in_header){
				/* End of header line: always keep this newline */
				fputc('\n', fp);
				in_header = 0;
			} else if (i + 1 >= fsize || data[i + 1] == '>'){
				/* End of sequence or EOF: keep newline */
				fputc('\n', fp);
			}
			/* Otherwise: skip newline to join continuation sequence lines */
		} else {
			fputc(data[i], fp);
		}
	}
	fclose(fp);
	free(data);
}
/* Convert FASTA to PHYLIP format (replaces fasta2phyml.pl) */
static void fasta_to_phylip(const char *fasta_path){
	FILE *fp = fopen(fasta_path, "r");
	if (!fp){ fprintf(stderr, "fasta_to_phylip: cannot open %s\n", fasta_path); return; }
	/* Build output filename: split on first '.' and append .phymlAln */
	char outpath[BUFFER_SIZE];
	strncpy(outpath, fasta_path, BUFFER_SIZE - 1);
	outpath[BUFFER_SIZE - 1] = '\0';
	char *dot = strchr(outpath, '.');
	if (dot) *dot = '\0';
	strncat(outpath, ".phymlAln", BUFFER_SIZE - strlen(outpath) - 1);
	/* Read all sequences */
	int capacity = 1024;
	int nseqs = 0;
	char **names = (char **)malloc(capacity * sizeof(char *));
	char **seqs = (char **)malloc(capacity * sizeof(char *));
	int maxNameLen = 0, maxSeqLen = 0;
	char line[65536];
	char *curName = NULL;
	char *curSeq = NULL;
	int curSeqLen = 0, curSeqCap = 0;
	while (fgets(line, sizeof(line), fp)){
		int len = strlen(line);
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
		if (line[0] == '>'){
			if (curName){
				if (nseqs >= capacity){ capacity *= 2; names = realloc(names, capacity * sizeof(char *)); seqs = realloc(seqs, capacity * sizeof(char *)); }
				names[nseqs] = curName;
				curSeq[curSeqLen] = '\0';
				seqs[nseqs] = curSeq;
				if (curSeqLen > maxSeqLen) maxSeqLen = curSeqLen;
				nseqs++;
			}
			int nameLen = len - 1;
			curName = (char *)malloc(nameLen + 1);
			memcpy(curName, line + 1, nameLen);
			curName[nameLen] = '\0';
			if (nameLen > maxNameLen) maxNameLen = nameLen;
			curSeqCap = 4096;
			curSeqLen = 0;
			curSeq = (char *)malloc(curSeqCap);
		} else {
			if (curSeq){
				while (curSeqLen + len + 1 > curSeqCap){ curSeqCap *= 2; curSeq = realloc(curSeq, curSeqCap); }
				memcpy(curSeq + curSeqLen, line, len);
				curSeqLen += len;
			}
		}
	}
	if (curName){
		if (nseqs >= capacity){ capacity *= 2; names = realloc(names, capacity * sizeof(char *)); seqs = realloc(seqs, capacity * sizeof(char *)); }
		names[nseqs] = curName;
		curSeq[curSeqLen] = '\0';
		seqs[nseqs] = curSeq;
		if (curSeqLen > maxSeqLen) maxSeqLen = curSeqLen;
		nseqs++;
	}
	fclose(fp);
	/* Write PHYLIP format */
	FILE *out = fopen(outpath, "w");
	if (!out){ fprintf(stderr, "fasta_to_phylip: cannot open %s for writing\n", outpath); goto cleanup; }
	fprintf(out, "  %d %d\n", nseqs, maxSeqLen);
	int i, j;
	for (i = 0; i < nseqs; i++){
		fputs(names[i], out);
		int padName = 4 + maxNameLen - (int)strlen(names[i]);
		for (j = 0; j < padName; j++) fputc(' ', out);
		fputs(seqs[i], out);
		int padSeq = maxSeqLen - (int)strlen(seqs[i]);
		for (j = 0; j < padSeq; j++) fputc('-', out);
		fputc('\n', out);
	}
	fclose(out);
cleanup:
	for (i = 0; i < nseqs; i++){ free(names[i]); free(seqs[i]); }
	free(names);
	free(seqs);
}
/* Copy a file using native C I/O (replaces system("cp ...")) */
static void copy_file(const char *src, const char *dst){
	FILE *in = fopen(src, "rb");
	if (!in){ fprintf(stderr, "copy_file: cannot open %s\n", src); return; }
	FILE *out = fopen(dst, "wb");
	if (!out){ fclose(in); fprintf(stderr, "copy_file: cannot open %s for writing\n", dst); return; }
	char buf[8192];
	size_t n;
	while ((n = fread(buf, 1, sizeof(buf), in)) > 0){
		fwrite(buf, 1, n, out);
	}
	fclose(in);
	fclose(out);
}
/* Get number of available CPU cores */
static int get_num_cores(void){
	long n = sysconf(_SC_NPROCESSORS_ONLN);
	return (n > 0) ? (int)n : 1;
}

/* Run the external tool pipeline for one partition:
   FAMSA -> unwrap -> [fasta2phyml] -> RAxML/VeryFastTree -> [nw_reroot]
   Returns 0 on success. */
static int run_partition_pipeline(int which, Options opt, int pipeline_threads){
	int status;
	char buf[BUFFER_SIZE];
	char buf2[BUFFER_SIZE];
	char buf3[BUFFER_SIZE];
	char famsa_threads_str[BUFFER_SIZE];

	/* Set thread limits for child processes */
	if (pipeline_threads > 0){
		char omp_str[32];
		snprintf(omp_str, sizeof(omp_str), "%d", pipeline_threads);
		setenv("OMP_NUM_THREADS", omp_str, 1);
	}

	/* Step 1: FAMSA alignment */
	if (opt.prefix[0] == '\0'){
		snprintf(buf2,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);
		snprintf(buf3,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
	}else{
		snprintf(buf2,BUFFER_SIZE,"%s/%spartition%d.fasta",opt.partitions_directory,opt.prefix,which);
		snprintf(buf3,BUFFER_SIZE,"%s/%spartition%d_MSA.fasta",opt.partitions_directory,opt.prefix,which);
	}
	int famsa_t = (opt.famsa_threads > 0) ? opt.famsa_threads : pipeline_threads;
	if (famsa_t <= 0) famsa_t = 1;
	snprintf(famsa_threads_str,BUFFER_SIZE,"%d",famsa_t);
	{
		pid_t pid = fork();
		if (pid == -1){
			fprintf(stderr, "can't fork for famsa, error occurred\n");
			return -1;
		}else if (pid == 0){
			int devnull = open("/dev/null", O_WRONLY);
			if (devnull >= 0){ dup2(devnull, STDERR_FILENO); close(devnull); }
			char *arguments[] = {"famsa","-t",famsa_threads_str,buf2,buf3,NULL};
			execvp("famsa",arguments);
			_exit(127);
		}else{
			waitpid(pid, &status, 0);
		}
	}

	/* Step 2: Unwrap FASTA */
	if (opt.prefix[0] == '\0'){
		snprintf(buf,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
	}else{
		snprintf(buf,BUFFER_SIZE,"%s/%spartition%d_MSA.fasta",opt.partitions_directory,opt.prefix,which);
	}
	unwrap_fasta_inplace(buf);

	/* Step 3: FASTA to PHYLIP conversion (RAxML only) */
	if (opt.fasttree == 0){
		fasta_to_phylip(buf);
	}

	/* Step 4: Tree inference */
	if (opt.fasttree == 0){
		/* RAxML path */
		if (opt.prefix[0] == '\0'){
			snprintf(buf,BUFFER_SIZE,"raxmlHPC-PTHREADS --silent -m GTRGAMMA -w %s/ -n partition%d -p 1234 -T %d -s %s/partition%d_MSA.phymlAln",
				opt.partitions_directory,which,pipeline_threads > 0 ? pipeline_threads : 8,opt.partitions_directory,which);
		}else{
			snprintf(buf,BUFFER_SIZE,"raxmlHPC-PTHREADS --silent -m GTRGAMMA -w %s/ -n %spartition%d -p 1234 -T %d -s %s/%spartition%d_MSA.phymlAln",
				opt.partitions_directory,opt.prefix,which,pipeline_threads > 0 ? pipeline_threads : 8,opt.partitions_directory,opt.prefix,which);
		}
		status = system(buf);
	}else{
		/* VeryFastTree path — fork/exec with stdout redirect instead of system() */
		char vft_input[BUFFER_SIZE], vft_output[BUFFER_SIZE], vft_threads[32];
		snprintf(vft_threads, sizeof(vft_threads), "%d",
			pipeline_threads > 0 ? pipeline_threads : get_num_cores());
		if (opt.prefix[0] == '\0'){
			snprintf(vft_input,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
			snprintf(vft_output,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d.reroot",opt.partitions_directory,which);
		}else{
			snprintf(vft_input,BUFFER_SIZE,"%s/%spartition%d_MSA.fasta",opt.partitions_directory,opt.prefix,which);
			snprintf(vft_output,BUFFER_SIZE,"%s/RAxML_bestTree.%spartition%d.reroot",opt.partitions_directory,opt.prefix,which);
		}
		pid_t vft_pid = fork();
		if (vft_pid == -1){
			fprintf(stderr, "can't fork for VeryFastTree\n");
			return -1;
		}else if (vft_pid == 0){
			int fd = open(vft_output, O_WRONLY | O_CREAT | O_TRUNC, 0666);
			if (fd >= 0){ dup2(fd, STDOUT_FILENO); close(fd); }
			int devnull = open("/dev/null", O_WRONLY);
			if (devnull >= 0){ dup2(devnull, STDERR_FILENO); close(devnull); }
			/* Try VeryFastTree first, fall back to FastTree */
			char *vft_args[] = {"VeryFastTree", "-gtr", "-gamma", "-nt", "-quiet", "-nosupport", "-threads", vft_threads, vft_input, NULL};
			execvp("VeryFastTree", vft_args);
			/* VeryFastTree not found — try FastTree (single-threaded, has support values) */
			char *ft_args[] = {"FastTree", "-gtr", "-gamma", "-nt", "-nosupport", vft_input, NULL};
			execvp("FastTree", ft_args);
			/* Neither found */
			fprintf(stderr, "Neither VeryFastTree nor FastTree found\n");
			_exit(127);
		}
		waitpid(vft_pid, &status, 0);
	}

	/* Step 5: Reroot tree (RAxML only) */
	if (opt.fasttree == 0){
		char buf4[BUFFER_SIZE];
		char buf5[BUFFER_SIZE];
		if (opt.prefix[0] == '\0'){
			snprintf(buf4,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d",opt.partitions_directory,which);
			snprintf(buf5,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d.reroot",opt.partitions_directory,which);
		}else{
			snprintf(buf4,BUFFER_SIZE,"%s/RAxML_bestTree.%spartition%d",opt.partitions_directory,opt.prefix,which);
			snprintf(buf5,BUFFER_SIZE,"%s/RAxML_bestTree.%spartition%d.reroot",opt.partitions_directory,opt.prefix,which);
		}
		pid_t pid = fork();
		if (pid == -1){
			fprintf(stderr, "can't fork for nw_reroot, error occurred\n");
			return -1;
		}else if (pid == 0){
			char *arguments[] = {"nw_reroot",buf4,NULL};
			int fd = open(buf5, O_WRONLY | O_CREAT, 0777);
			if (fd == -1){
				perror(buf5);
				_exit(-1);
			}
			fclose(stdout);
			dup2(fd, STDOUT_FILENO);
			close(fd);
			execvp(arguments[0],arguments);
			_exit(127);
		}else{
			waitpid(pid, &status, 0);
		}
	}

	return 0;
}

/* Remove partition files using unlink (replaces system("rm ...") calls) */
static void remove_partition_files(const char *partitions_directory, const char *tempfilename){
	char path[BUFFER_SIZE];
	snprintf(path, BUFFER_SIZE, "%s/%s_MSA.fasta", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/%s_taxonomy.txt", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_bestTree.%s", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_bestTree.%s.reroot", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/%s_MSA.phymlAln", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/%s_MSA.phymlAln.reduced", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_log.%s", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_result.%s", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_info.%s", partitions_directory, tempfilename);
	unlink(path);
	snprintf(path, BUFFER_SIZE, "%s/RAxML_parsimonyTree.%s", partitions_directory, tempfilename);
	unlink(path);
}
void createNewRoots(int rootCount, Options opt, int max_nodename, int max_lineTaxonomy, struct masterArr *m){
	int i,j,k,count;
	char buffer[BUFFER_SIZE];
	int *leaves;
	gzFile partition;
	double local_minVariance = 99999999999999;
	int local_minVarNode = -1;
	findMinVarianceArr(m->root,m->numspec,m,&local_minVariance,&local_minVarNode);
	int *partition1 = (int *)malloc(sizeof(int)*m->tree[m->tree[local_minVarNode].up[0]].nd);
	for(i=0; i<m->tree[m->tree[local_minVarNode].up[0]].nd; i++){
		partition1[i]=-1;
	}
	int *partition2 = (int *)malloc(sizeof(int)*m->tree[m->tree[local_minVarNode].up[1]].nd);
	for(i=0; i<m->tree[m->tree[local_minVarNode].up[1]].nd; i++){
		partition2[i]=-1;
	}
	int *partition3 = (int *)malloc(sizeof(int)*(m->tree[m->root].nd-m->tree[local_minVarNode].nd));
	for(i=0; i<(m->tree[m->root].nd-m->tree[local_minVarNode].nd); i++){
		partition3[i]=-1;
	}
	partition1 = findLeavesOfMinVarArr(m->tree[local_minVarNode].up[0],partition1,m->tree[m->tree[local_minVarNode].up[0]].nd,m,local_minVarNode);
	partition2 = findLeavesOfMinVarArr(m->tree[local_minVarNode].up[1],partition2,m->tree[m->tree[local_minVarNode].up[1]].nd,m,local_minVarNode);
	partition3 = findLeavesOfMinVarArr(m->root,partition3,(m->tree[m->root].nd-m->tree[local_minVarNode].nd),m,local_minVarNode);
	int *partitionSizes = (int *)malloc(sizeof(int)*3);
	partitionSizes[0]=m->tree[m->tree[local_minVarNode].up[0]].nd;
	partitionSizes[1]=m->tree[m->tree[local_minVarNode].up[1]].nd;
	partitionSizes[2]=(m->tree[m->root].nd-m->tree[local_minVarNode].nd);
	if ( partitionSizes[0] < 4 || partitionSizes[1] < 4 || partitionSizes[2] < 4){
		pthread_mutex_lock(&spscore_mutex);
		SPscoreArr[rootCount]=1;
		pthread_mutex_unlock(&spscore_mutex);
		return;
	}
	int partitionCount=-1;
	int freeSlots=0;
	pthread_mutex_lock(&spscore_mutex);
	for(i=MAX_NUMBEROFROOTS-1;i>=0;i--){
		if(SPscoreArr[i]==-1){
			partitionCount=i;
			freeSlots++;
		}
	}
	if (freeSlots < 3){
		pthread_mutex_unlock(&spscore_mutex);
		fprintf(stderr, "ERROR: partition slot limit reached (%d). Too many partitions for this dataset/threshold combination.\n", MAX_NUMBEROFROOTS);
		exit(EXIT_FAILURE);
	}
	SPscoreArr[partitionCount]=1;
	SPscoreArr[partitionCount+1]=1;
	SPscoreArr[partitionCount+2]=1;
	pthread_mutex_unlock(&spscore_mutex);
	partitionCount++;
	if (partitionCount > opt.restart){
		printPartitionsToFileArr(partition1,partitionSizes[0],partition2,partitionSizes[1],partition3,partitionSizes[2],rootCount,partitionCount,opt,m);
	}
	for(i=0; i<m->numspec; i++){
		for(j=0; j<7; j++){
			free(m->taxonomy[i][j]);
		}
		free(m->taxonomy[i]);
		free(m->msa[i]);
		free(m->names[i]);
	}
	for(i=0; i<2*m->numspec-1; i++){
		free(m->tree[i].like);
		free(m->tree[i].posterior);
		free(m->tree[i].name);
	}
	free(m->taxonomy);
	free(m->msa);
	free(m->names);
	free(m->tree);
	free(m->filename);
	pthread_mutex_lock(&mastermap_mutex);
	hashmap_remove(&mastermap,m->index);
	pthread_mutex_unlock(&mastermap_mutex);
	free(m);
	rootCount=partitionCount;
	int which,status;
	char buf[BUFFER_SIZE];
	/* Fork all 3 partition pipelines concurrently */
	int total_cores = get_num_cores();
	int pipeline_threads = total_cores / 3;
	if (pipeline_threads < 1) pipeline_threads = 1;
	pid_t pipeline_pids[3];
	int p;
	for(p = 0; p < 3; p++){
		which = rootCount + p;
		if (rootCount > opt.restart){
			pipeline_pids[p] = fork();
			if (pipeline_pids[p] == -1){
				fprintf(stderr, "can't fork pipeline for partition %d\n", which);
			}else if (pipeline_pids[p] == 0){
				/* Child process: run the entire pipeline */
				int rc = run_partition_pipeline(which, opt, pipeline_threads);
				_exit(rc);
			}
		}else{
			pipeline_pids[p] = -1;
		}
	}
	/* Wait for all 3 pipelines to complete */
	for(p = 0; p < 3; p++){
		if (pipeline_pids[p] > 0){
			waitpid(pipeline_pids[p], &status, 0);
		}
	}
	/* Sequential loading/parsing/scoring of the 3 partitions */
	for(which=rootCount;which<rootCount+3;which++){
		if ( opt.prefix[0] == '\0' ){
			snprintf(buf,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
		}else{
			snprintf(buf,BUFFER_SIZE,"%s/%spartition%d_MSA.fasta",opt.partitions_directory,opt.prefix,which);
		}
		printf("buffer: %s\n",buf);
		struct masterArr *t=malloc(sizeof(masterArr));
		t->tree = NULL;
		t->treeCapacity = 0;
		t->filename = (malloc)(300*sizeof(char));
		for(i=0; i<300; i++){
			t->filename[i] = '\0';
		}
		sprintf(t->index,"%d",which-1);
		if ((partition = gzopen(buf,"r")) == NULL ){
			fprintf(stderr, "Cannot open %s: %s\n", buf, strerror(errno));
			exit(EXIT_FAILURE);
		}
		t->numspec=setNumspecArr(partition);
		gzclose(partition);
		t->tree=(struct node*)malloc((2*t->numspec-1)*sizeof(struct node));
		t->names=(char**)malloc(sizeof(char*)*t->numspec);
		for(i=0; i<t->numspec; i++){
			t->names[i]=(char*)malloc(sizeof(char)*(max_nodename+1));
		}
		t->msa=(int**)malloc(t->numspec*sizeof(int*));
		t->taxonomy=(char***)malloc(t->numspec*sizeof(char**));
		if ((partition = gzopen(buf,"r")) == NULL ){
			fprintf(stderr, "Cannot open %s: %s\n", buf, strerror(errno));
			exit(EXIT_FAILURE);
		}
		readSeqArr(partition,max_nodename,t);
		gzclose(partition);
		if ( opt.prefix[0] == '\0' ){
			snprintf(buf,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d.reroot",opt.partitions_directory,which);
		}else{
			snprintf(buf,BUFFER_SIZE,"%s/RAxML_bestTree.%spartition%d.reroot",opt.partitions_directory,opt.prefix,which);
		}
		strcpy(t->filename,buf);
		if ( opt.fasttree == 0 ){
			allocateTreeArrMemory(t,max_nodename);
			comma=0;
			tip=0;
			t->root=getcladeArr_fast(buf,t,max_nodename)-1;
			t->tree[t->root].down = -1;
		}
		if ( opt.fasttree == 1 ){
			t->numspec = 0;
			t->numNodes = 0;
			FILE *fasttreefile = fopen(buf,"r");
			if (fasttreefile == NULL ){
				fprintf(stderr, "Error: Could not open Newick file %s.\n", buf);
				exit(EXIT_FAILURE);
			}
			char* newick = readNewickFile(fasttreefile);
			fclose(fasttreefile);
			srand(time(NULL));
			parseNewick(t, newick, max_nodename);
			free(newick);
			makeBinary(t,max_nodename);
			exportTreeToNewick(t,buf);
			/* Free parseNewick tree and re-allocate for getcladeArr_fast */
			int idx;
			for(idx=0; idx<t->numNodes; idx++){
				free(t->tree[idx].name);
				free(t->tree[idx].like);
				free(t->tree[idx].posterior);
			}
			free(t->tree);
			t->numspec = t->numspec;  /* parseNewick set this */
			t->tree = malloc((2*t->numspec-1) * sizeof(node));
			t->treeCapacity = 2*t->numspec-1;
			t->numNodes = 2*t->numspec-1;
			allocateTreeArrMemory(t,max_nodename);
			comma=0;
			tip=0;
			t->root=getcladeArr_fast(buf,t,max_nodename)-1;
			t->tree[t->root].down = -1;
		}
		get_number_descendantsArr(t->root,t);
		int child0 = t->tree[t->root].up[0];
		int child1 = t->tree[t->root].up[1];
		t->tree[t->root].depth = 0;
		assignDepthArr(child0,child1,1,t);
		if ( opt.prefix[0] == '\0' ){
			snprintf(buf,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);
		}else{
			snprintf(buf,BUFFER_SIZE,"%s/%spartition%d_taxonomy.txt",opt.partitions_directory,opt.prefix,which);
		}
		printf("buffer: %s\n",buf);
		for(j=0; j<t->numspec; j++){
			t->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
			for(k=0; k<7; k++){
				t->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
			}
		}
		assignTaxonomyToLeavesArr(buf,t,max_nodename,max_lineTaxonomy);
		{ int _tax_out[2]; getTaxonomyArr(t->root,t,_tax_out); }
		pthread_mutex_lock(&mastermap_mutex);
		hashmap_put(&mastermap, t->index, t);
		pthread_mutex_unlock(&mastermap_mutex);
		double initialSPscore=-1;
		if ( opt.use_spscore==1 && opt.use_min_leaves==0){
				initialSPscore = calculateSPArr(t);
				if (initialSPscore < opt.sp_score){
					pthread_mutex_lock(&spscore_mutex);
					SPscoreArr[which-1]=0;
					pthread_mutex_unlock(&spscore_mutex);
					createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
				}else{
					pthread_mutex_lock(&spscore_mutex);
					SPscoreArr[which-1]=1;
					pthread_mutex_unlock(&spscore_mutex);
				}
		}
		if ( opt.use_spscore==0 && opt.use_min_leaves==1  ){
			if ( t->numspec > opt.min_leaves ){
				pthread_mutex_lock(&spscore_mutex);
				SPscoreArr[which-1]=0;
				pthread_mutex_unlock(&spscore_mutex);
				createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
			}else{
				pthread_mutex_lock(&spscore_mutex);
				SPscoreArr[which-1]=1;
				pthread_mutex_unlock(&spscore_mutex);
			}
		}
		if ( opt.use_spscore==1 && opt.use_min_leaves==1 ){
			initialSPscore = calculateSPArr(t);
			if (initialSPscore < opt.sp_score && t->numspec > opt.min_leaves){
				createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
			}else{
				pthread_mutex_lock(&spscore_mutex);
				SPscoreArr[which-1]=1;
				pthread_mutex_unlock(&spscore_mutex);
			}
		}
	}
	free(partition1);
	free(partition2);
	free(partition3);
	free(partitionSizes);
}
void populate(int* nodes, char** taxnames, int node){
	int child0 = treeArr[0][node].up[0];
	int child1 = treeArr[0][node].up[1];
	int parent = treeArr[0][node].down;
	int i,j;
	int no_add=0;
	//if ( treeArr[0][node].taxIndex[1] == 2 ){ 
		for(i=0; i<100; i++){
			//if (strcmp(taxonomyArr[0][treeArr[0][node].taxIndex[0]][treeArr[0][node].taxIndex[1]],taxnames[i])==0){
			if (strcmp(taxonomyArr[0][treeArr[0][node].taxIndex[0]][2],taxnames[i])==0){
				if ( treeArr[0][node].depth > treeArr[0][nodes[i]].depth ){
					nodes[i]=node;
				}
				no_add=1;
			}
		}
		if ( no_add == 0 ){
			for(i=0; i<100; i++){
				if ( nodes[i]==-1 ){
					break;
				}
			}
			nodes[i] = node;
			//strcpy(taxnames[i],taxonomyArr[0][treeArr[0][node].taxIndex[0]][treeArr[0][node].taxIndex[1]]);
			strcpy(taxnames[i],taxonomyArr[0][treeArr[0][node].taxIndex[0]][2]);
		}
	//}else{
		if ( child0 != -1 && child1 != -1 ){
			populate(nodes,taxnames,child0);
			populate(nodes,taxnames,child1);
		}
	//}

}
/*void printFamilySeqs(){
	int* nodes = (int*)malloc(sizeof(int)*100);
	int i, j;
	for(i=0; i<100; i++){
		nodes[i]=-1;
	}
	char** taxnames = (char**)malloc(sizeof(char*)*100);
	for(i=0; i<100; i++){
		taxnames[i]=malloc(sizeof(char)*100);
		for(j=0; j<100; j++){
			taxnames[i][j]='\0';
		}
	}
	populate(nodes,taxnames,rootArr[0]);
	for(i=0; i<100; i++){
		if ( nodes[i] == -1 ){ break; }
		char* seq = (char*)malloc(sizeof(char)*1000);
		for(j=0; j<1000; j++){
			seq[j] = '\0';
		}
		printf(">%s\n",taxnames[i]);
		getSequenceInNode(0,nodes[i],seq);
		printf("%s\n",seq);
	}
}*/
int createNode( struct masterArr* m, int max_nodename){
	int newIndex = m->numNodes++;
	if (m->numNodes > m->treeCapacity){
		/* Double capacity to amortize realloc cost and reduce pointer invalidation */
		int newCap = m->treeCapacity * 2;
		if (newCap < m->numNodes) newCap = m->numNodes;
		if (newCap < 64) newCap = 64;
		node* newTree = realloc(m->tree, newCap * sizeof(node));
		if ( newTree == NULL ){
			fprintf(stderr, "Error: Memory allocation failed in createNode().\n");
			exit(EXIT_FAILURE);
		}
		m->tree = newTree;
		m->treeCapacity = newCap;
	}
	node* newNode = &m->tree[newIndex];
	int i;
	for(i=0; i<MAX_NODE_CHILDREN; i++){
		newNode->up[i] = -1;
	}
	newNode->down = -2;
	newNode->nd = 0;
	newNode->depth = 0;
	newNode->bl = 0.0;
	newNode->name = malloc((max_nodename+1)*sizeof(char));
	memset(newNode->name, '\0', max_nodename+1);  // name is 20 bytes
	newNode->like = malloc(STATESPACE*(sizeof(double)));
	newNode->posterior = malloc(STATESPACE*(sizeof(double)));
	newNode->likenc = NULL;
	newNode->posteriornc = NULL;
	return newIndex;
}
int addChild( struct masterArr* m, int parentIndex, int childIndex){
	node* parent = &m->tree[parentIndex];
	if (parent->nd >= MAX_NODE_CHILDREN - 1) {
		/* Overflow: create intermediate node to hold half the children,
		   keeping parent degree manageable for resolvePolytomy later */
		int intermediateIndex = createNode(m, 30);
		parent = &m->tree[parentIndex]; /* refresh after potential realloc */
		node* intermediate = &m->tree[intermediateIndex];
		int half = parent->nd / 2;
		int i;
		for (i = half; i < parent->nd; i++) {
			intermediate->up[intermediate->nd] = parent->up[i];
			m->tree[parent->up[i]].down = intermediateIndex;
			intermediate->nd++;
			parent->up[i] = -1;
		}
		parent->nd = half;
		intermediate->down = parentIndex;
		parent->up[parent->nd] = intermediateIndex;
		parent->nd++;
	}
	parent = &m->tree[parentIndex];
	node* child = &m->tree[childIndex];
	parent->up[parent->nd] = childIndex;
	parent->nd++;
	child->down = parentIndex;
}
void parseNewick( struct masterArr* m, const char* newick, int max_nodename ){
	int current = -1;
	int tipIndex = 0;
	int parent = -1;
	/* Pre-allocate tree capacity: count '(' and ',' to estimate nodes.
	   Each '(' is an internal node, each ',' + 1 gives leaf count within a clade.
	   Upper bound: 2 * (commas + 1) for a binary tree. Add 50% extra for polytomy resolution. */
	{
		int est_nodes = 0;
		const char *q = newick;
		while (*q){ if (*q == '(' || *q == ',') est_nodes++; q++; }
		est_nodes = est_nodes * 3; /* generous overallocation for resolvePolytomy */
		if (est_nodes < 64) est_nodes = 64;
		if (m->tree) free(m->tree);
		m->tree = malloc(est_nodes * sizeof(node));
		m->treeCapacity = est_nodes;
		m->numNodes = 0;
	}
	const char* p = newick;
	while( *p ){
		if ( *p == '(' ){
			int newIndex = createNode(m,max_nodename);
			if ( current != -1 ){
				addChild(m, current, newIndex);
			}else{
				m->root = newIndex;
			}
			current = newIndex;
			p++;
		}else if ( *p == ')' ){
			// Close current internal node and move up
			if ( current == -1 ){
				fprintf(stderr, "Error: Mismatched parentheses in Newick string.\n");
				exit(EXIT_FAILURE);
			}
			p++;
			// After ')', there may be a node label or branch length
			const char* start = p;
			// Parse node label if present
			if (*p && strchr("():,;", *p) == NULL && *p != ':') {
				// Read node label
				while (*p && strchr("():,;", *p) == NULL && *p != ':') {
					p++;
				}
				size_t len = p - start;
				if ( len > 0 ){
					//m->tree[current].name = malloc(len + 1);
					strncpy(m->tree[current].name, start, len);
					m->tree[current].name[len] = '\0';
				}
			}
			// Check for branch length
			if (*p == ':') {
				p++;
				m->tree[current].bl = strtod(p, (char**)&p);
			}
			// Check if we're at the root (no parent)
			if ( current == m->root ){
				// At root, no parent to move up to; simply continue
				current = -1;
			}else{
				// Otherwise, move up to the parent node
				current = m->tree[current].down;
			}
		}else if ( *p == ':' ){
			p++;
			m->tree[current].bl = strtod(p, (char**)&p);
		}else if ( *p == ','){
			// Move to sibling processing
			p++;
		} else if ( *p == ';' ){
			break;
		}else{
			//Parse node name
			const char* start = p;
			while (*p && strchr("():,;", *p) == NULL) {
				p++;
			}
			size_t len = p - start;
			if ( len > 0 ){
				int newIndex = createNode(m,max_nodename);
				//m->tree[newIndex].name = malloc( len + 1);
				strncpy(m->tree[newIndex].name, start, len);
				m->tree[newIndex].name[len] = '\0';
				m->tree[newIndex].nd = 0;
				m->tree[newIndex].depth = 1;
				m->numspec++;

				// Link to the current node
				addChild(m, current, newIndex);
				current = newIndex;
			}
			// Check for branch length
			if (*p == ':') {
				p++;
				m->tree[current].bl = strtod(p, (char**)&p);
			}
			// Move back up to parent
			current = m->tree[current].down;
		}
	}
	// Ensure we finish at the root
	if ( current != -1 ){
		fprintf(stderr, "Error: Unmatched opening parenthesis in Newick string.\n");
		exit(EXIT_FAILURE);
	}
}
void resolvePolytomy(struct masterArr* m, int nodeIndex, int max_nodename){
	while( m->tree[nodeIndex].nd > 2 ){
		node* current = &m->tree[nodeIndex];
		// Create a new internal node
		int newNodeIndex = createNode(m,max_nodename);
		current = &m->tree[nodeIndex];
		node* newNode = &m->tree[newNodeIndex];
		int nd = current->nd;
		//Randomly select two distinct children to combine
		int idx1 = rand() % nd;
		int idx2;
		do {
			idx2 = rand() % nd;
		}while(idx2 == idx1);

		//Get the child indices
		int childIndex1 = current->up[idx1];
		int childIndex2 = current->up[idx2];

		//Ensure idx1 and idx2 are in order to remove them correctly
		if ( idx1 > idx2 ){
			int tempIdx = idx1;
			idx1 = idx2;
			idx2 = tempIdx;

			int tempChild = childIndex1;
			childIndex1 = childIndex2;
			childIndex2 = tempChild;
		}

		// Remove the moved child from the current node's `up[]` array
		for (int i = idx1; i < current->nd - 1; i++) {
			current->up[i] = current->up[i + 1];
		}
		current->nd--;  // Decrement the number of descendants for the parent
	
		//Adjust idx2 since the array has shifted
		idx2--;
		for (int i = idx2; i < current->nd - 1; i++) {
			current->up[i] = current->up[i + 1];
		}
		current->nd--;

		// Clear the now unused slots
		current->up[current->nd] = -1;
		current->up[current->nd + 1] = -1;

		// Set the new internal node's children and parent relationships
		newNode->up[0] = childIndex1;
		newNode->up[1] = childIndex2;
		newNode->nd = 2;
		newNode->down = nodeIndex;

		// Update the moved children's parent relationships
		m->tree[childIndex1].down = newNodeIndex;
		m->tree[childIndex2].down = newNodeIndex;

		// Add the new internal node as a child of the current node
		current->up[current->nd] = newNodeIndex;
		current->nd++;  // Increment the number of descendants for the parent

		//printf("Resolving polytomy at node %d: combined children %d and %d into new internal node %d\n",nodeIndex, childIndex1, childIndex2, newNodeIndex);

		// Debugging: Print the current and newNode state
		//printf("Updated node %d: nd=%d, up={", nodeIndex, current->nd);
		//for (int i = 0; i < current->nd; i++) {
		//	printf("%d ", current->up[i]);
		//}
		//printf("}\n");
		//printf("New Node %d: nd=%d, up={%d, %d}, down=%d\n", newNodeIndex, newNode->nd, newNode->up[0], newNode->up[1], newNode->down);
	}
}
void makeBinary( struct masterArr* m, int max_nodename){
	/* Pre-allocate tree capacity for all nodes that resolvePolytomy will create.
	   Each polytomy of degree d creates d-2 new internal nodes. */
	int extraNodes = 0;
	int i;
	for(i=0; i<m->numNodes; i++){
		if (m->tree[i].nd > 2){
			extraNodes += m->tree[i].nd - 2;
		}
	}
	if (extraNodes > 0) {
		int needed = m->numNodes + extraNodes;
		if (needed > m->treeCapacity) {
			int newCap = needed + 64;
			node* newTree = realloc(m->tree, newCap * sizeof(node));
			if (newTree == NULL) {
				fprintf(stderr, "Error: Memory allocation failed in makeBinary pre-alloc.\n");
				exit(EXIT_FAILURE);
			}
			m->tree = newTree;
			m->treeCapacity = newCap;
		}
	}
	for(i=0; i<m->numNodes; i++){
		node* n = &m->tree[i];
		// Skip leaves
		if ( n->nd <= 2 ){
			continue;
		}
		// Resolve polytomy for nodes with more than 2 children
		resolvePolytomy(m,i,max_nodename);
	}
}
int main(int argc, char **argv){
	Options opt;
	opt.min_leaves=0;
	opt.use_partitions=0;
	opt.use_spscore=0;
	opt.use_min_leaves=0;
	opt.sp_score = 0.05;
	opt.fasttree=0;
	opt.missing_data=1;
	opt.restart = 0;
	opt.two_step = 0;
	opt.number_of_partitions = 0;
	opt.number_of_trees = 0;
	opt.famsa_threads = 0; /* 0 = auto-detect in pipeline */
	opt.remove_unused = 0;
	opt.export_subtrees = 0;
	int i, j, k, numberOfTrees;
	for(i=0; i<200; i++){
		opt.partitions_directory[i] = '\0';
		opt.readdir[i] = '\0';
	}
	for(i=0; i<2000; i++){
		opt.prefix[i] = '\0';
	}
	parse_options(argc, argv, &opt);
	if ( opt.partitions_directory[0] == '\0' ){
		fprintf(stderr, "You must specify an output directory with -d\n");
		exit(-1);
	}
	int max_nodename = 0;
	int max_tax_name = 0;
	int max_lineTaxonomy = 0;
	hashmap_init(&mastermap,hashmap_hash_string,strcmp);
	if (opt.number_of_trees==1 && opt.use_partitions==0){
		printf("Using a single tree... \n");
		numberOfTrees=1;
		struct masterArr *m = malloc(sizeof(masterArr));
		m->tree = malloc(sizeof(node *));
		m->treeCapacity = 0;
		int *specifications = (int*)malloc(3*sizeof(int));
		specifications[0]=0;
		specifications[1]=0;
		specifications[2]=0;
		gzFile infile;
		if (( infile = gzopen(opt.msa_file,"r")) == NULL) fprintf(stderr,"MSA file could not be opened.\n");
		setNumspec(infile,specifications);
		gzclose(infile);
		m->numspec = specifications[0];
		max_nodename = specifications[1];
		m->numbase = specifications[2];
		free(specifications);
		initlogfactorial();
		//nodeIDsArr = (char ***)malloc(sizeof(char**));
		//itoa(0,m->index,10);
		sprintf(m->index,"%d",0);
		m->tree = (struct node*)malloc((2*m->numspec-1)*sizeof(struct node));
		m->names = (char **)malloc(sizeof(char *)*m->numspec);
		for(i=0; i<m->numspec;i++) {
			m->names[i] = malloc((max_nodename+1)*sizeof(char));
		}
		//seqArr = (int ***)malloc(sizeof(int **));
		m->msa = (int **)malloc(m->numspec*sizeof(int *));
		for(i=0; i<m->numspec; i++){
			m->msa[i]=(int *)malloc(m->numbase*sizeof(int));
		}
		if (( infile = gzopen(opt.msa_file,"r")) == NULL) fprintf(stderr,"MSA file could not be opened.\n");
		readseq(infile,max_nodename,m);
		gzclose(infile);
		allocateTreeArrMemory(m,max_nodename);
		comma=0;
		tip=0;
		m->root=getcladeArr_fast(opt.tree_file,m,max_nodename)-1;
		m->tree[m->root].down = -1;
		get_number_descendantsArr(m->root,m);
		int child0 = m->tree[m->root].up[0];
		int child1 = m->tree[m->root].up[1];
		m->tree[m->root].depth=0;
		assignDepthArr(child0,child1,1,m);
		int* tax_specs = (int*)malloc(2*sizeof(int));
		tax_specs[0]=0;
		tax_specs[1]=0;
		FILE* taxfile;
		if (( taxfile = fopen(opt.taxonomy_file,"r")) == (FILE *) NULL) fprintf(stderr,"Taxonomy file could not be opened.\n");
		findMaxTaxName(taxfile,tax_specs);
		fclose(taxfile);
		max_tax_name = tax_specs[0];
		max_lineTaxonomy = tax_specs[1];
		free(tax_specs);
		m->taxonomy = (char ***)calloc_check(m->numspec, sizeof(char **));
		for(j=0; j<m->numspec; j++){
			m->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
			for(k=0; k<7; k++){
				m->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
			}
		}
		assignTaxonomyToLeavesArr(opt.taxonomy_file,m,max_nodename,max_tax_name);
		{ int _tax_out[2]; getTaxonomyArr(m->root,m,_tax_out); }
		hashmap_put(&mastermap,m->index,m);
		/*treeArr = malloc(sizeof(node*));
		treeArr[0] = m->tree;
		taxonomyArr = (char****)malloc(sizeof(char***));
		taxonomyArr[0] = m->taxonomy;
		seqArr = (int***)malloc(sizeof(int**));
		seqArr[0] = m->msa;
		numbaseArr = (int*)malloc(sizeof(int));
		numbaseArr[0] = m->numbase;
		numspecArr = (int*)malloc(sizeof(int));
		numspecArr[0] = m->numspec;
		rootArr = (int*)malloc(sizeof(int));
		rootArr[0] = m->root;*/
	}else{
		printf("Using %d trees... \n",opt.number_of_partitions);
		//treeArr = malloc(sizeof(node *)*opt.number_of_partitions);
		SPscoreArr = malloc(sizeof(int)*MAX_NUMBEROFROOTS);
		for(i=0;i<MAX_NUMBEROFROOTS;i++){
			SPscoreArr[i]=-1;
		}
		//numspecArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		//numbaseArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		//rootArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		partition_files* pf = malloc(sizeof(partition_files));
		pf->tree_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		pf->msa_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		pf->tax_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		for (i=0; i<opt.number_of_partitions; i++){
			pf->tree_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
			pf->msa_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
			pf->tax_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
		}
		readFilesInDir(opt.readdir,opt.number_of_partitions,pf);
		printf("Input files...\n");
		for (i=0; i<opt.number_of_partitions; i++){
			printf("%d : %s\n",i,pf->msa_files[i]);
			printf("%d : %s\n",i,pf->tree_files[i]);
			printf("%d : %s\n",i,pf->tax_files[i]);
		}
		int* tax_specs = (int*)malloc(2*sizeof(int));
		tax_specs[0]=0;
		tax_specs[1]=0;
		char buffer[BUFFER_SIZE];
		int* specifications = (int*)malloc(3*sizeof(int));
		specifications[0]=0;
		specifications[1]=0;
		specifications[2]=0;
		for(i=0; i<opt.number_of_partitions; i++){
			FILE *taxfiles;
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tax_files[i]);
			if (NULL==(taxfiles=fopen(buffer,"r"))){ puts ("Cannot open tax file!"); exit(-1); }
			findMaxTaxName(taxfiles,tax_specs);
			fclose(taxfiles);
			gzFile msafiles;
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->msa_files[i]);
			if ((msafiles=gzopen(buffer,"r"))==NULL){
				fprintf(stderr, "Cannot open %s: %s\n", buffer, strerror(errno));
				exit(EXIT_FAILURE);
			}
			setNumspec(msafiles,specifications);
			gzclose(msafiles);
		}
		max_tax_name = tax_specs[0];
		max_lineTaxonomy = tax_specs[1];
		free(tax_specs);
		int max_numspec = specifications[0];
		max_nodename = specifications[1];
		int max_numbase = specifications[2];
		free(specifications);
		gzFile partition;
		int status;
		int partition_count = opt.number_of_partitions;
		for(i=0; i<partition_count; i++){
			SPscoreArr[i]=0;
		}
		//taxonomyArr = (char ****)malloc(opt.number_of_partitions*sizeof(char ***));
		//nodeIDsArr = (char ***)malloc(opt.number_of_partitions*sizeof(char**));
		//seqArr = (int ***)malloc(opt.number_of_partitions*sizeof(int **));
		/*for(i=0; i<MAX_NUMBEROFROOTS; i++){
			nodeIDsArr[i] = (char**)malloc(sizeof(char*)*max_numspec);
			seqArr[i] = (int **)malloc(max_numspec*sizeof(int*));
			for(j=0; j<max_numspec; j++){
				nodeIDsArr[i][j]=malloc(max_nodename*sizeof(char));
				seqArr[i][j]=malloc(max_numbase*sizeof(int));
				for(k=0; k<max_numbase; k++){
					nodeIDsArr[i][j][k]='\0';
					seqArr[i][j][k]=0;
				}
			}
		}*/
		for(i=0; i<partition_count; i++){
			struct masterArr *m = malloc(sizeof(masterArr));
			m->tree = NULL;
			m->treeCapacity = 0;
			sprintf(m->index,"%d",i);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->msa_files[i]);
			if ((partition=gzopen(buffer,"r"))==NULL){
				fprintf(stderr, "Cannot open %s: %s\n", buffer, strerror(errno));
				exit(EXIT_FAILURE);	
			}
			m->numspec = setNumspecArr(partition);
			printf("m->numspec: %d\n",m->numspec);
			gzclose(partition);
			m->numNodes = 0;
			m->tree=(struct node*)malloc((2*m->numspec-1)*sizeof(struct node));
			m->msa=(int**)malloc(m->numspec*sizeof(int*));
			m->taxonomy=(char***)malloc(m->numspec*sizeof(char**));
			m->names=(char**)malloc(m->numspec*sizeof(char*));
			m->filename = (malloc)(300*sizeof(char));
			for(j=0; j<300; j++){
				m->filename[j] = '\0';
			}
			for(j=0;j<m->numspec;j++){
				m->names[j]=(char*)malloc(sizeof(char)*(max_nodename+1));
			}
			if ((partition=gzopen(buffer,"r"))==NULL){
				fprintf(stderr,"Cannot open %s: %s\n", buffer, strerror(errno));
				exit(EXIT_FAILURE);
			}
			readSeqArr(partition,max_nodename,m);
			gzclose(partition);
			allocateTreeArrMemory(m,max_nodename);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tree_files[i]);
			strcpy(m->filename,buffer);
			if ( opt.fasttree == 1 ){
				/* VeryFastTree can produce multifurcating trees.
				   Use parseNewick+makeBinary to resolve polytomies,
				   export the binary tree, then re-parse with getcladeArr_fast. */
				FILE *tmpTree = fopen(buffer, "r");
				if (!tmpTree){ fprintf(stderr, "*** tree file could not be opened: %s\n", buffer); exit(-1); }
				char* newick = readNewickFile(tmpTree);
				fclose(tmpTree);
				int saved_numspec = m->numspec;
				m->numspec = 0;
				m->numNodes = 0;
				srand(time(NULL));
				parseNewick(m, newick, max_nodename);
				free(newick);
				makeBinary(m, max_nodename);
				exportTreeToNewick(m, buffer);
				/* Free parseNewick tree and re-allocate for getcladeArr_fast */
				int idx;
				for(idx=0; idx<m->numNodes; idx++){
					free(m->tree[idx].name);
					free(m->tree[idx].like);
					free(m->tree[idx].posterior);
				}
				free(m->tree);
				m->numspec = saved_numspec;
				m->tree = malloc((2*m->numspec-1) * sizeof(node));
				m->treeCapacity = 2*m->numspec-1;
				m->numNodes = 2*m->numspec-1;
				allocateTreeArrMemory(m,max_nodename);
			}
			comma=0;
			tip=0;
			m->root=getcladeArr_fast(buffer,m,max_nodename)-1;
			m->tree[m->root].down = -1;
			get_number_descendantsArr(m->root,m);
			int child0 = m->tree[m->root].up[0];
			int child1 = m->tree[m->root].up[1];
			m->tree[m->root].depth=0;
			assignDepthArr(child0,child1,1,m);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tax_files[i]);
			m->taxonomy = (char ***)calloc_check(m->numspec, sizeof(char **));
			for(j=0; j<m->numspec; j++){
				m->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
				for(k=0; k<7; k++){
					m->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
				}
			}
			assignTaxonomyToLeavesArr(buffer,m,max_nodename,max_lineTaxonomy);
			{ int _tax_out[2]; getTaxonomyArr(m->root,m,_tax_out); }
			hashmap_put(&mastermap, m->index, m);
			if ( opt.use_partitions==1 && m->numspec > opt.min_leaves ){
				SPscoreArr[i]=0;
				createNewRoots(i,opt,max_nodename,max_lineTaxonomy,m);
			}
			free(pf->tax_files[i]);
			free(pf->msa_files[i]);
			free(pf->tree_files[i]);
		}
		free(pf->tax_files);
		free(pf->msa_files);
		free(pf->tree_files);
		free(pf);
	}
	free(SPscoreArr);
		if ( opt.two_step == 1 ){
			printf("Partitioning has finished and you've requested a two-step build ('-p') exiting...\n");
			int numberOfTrees=0;
			int key;
			struct masterArr* final;
			if ( opt.remove_unused == 1 ){
				struct dirent *de;
				DIR *dr = opendir(opt.partitions_directory);
				regex_t regex1;
				char msgbuf[100];
				int reti = regcomp(&regex1, "_MSA\\.fasta$",0);
				if (reti){
					fprintf(stderr, "Could not compile regex\n");
					exit(-1);
				}
				if (dr == NULL){
					printf("Could not open directory for reads\n");
					exit(-1);
				}
				while ((de = readdir(dr)) != NULL){
					reti = regexec(&regex1, de->d_name, 0, NULL, 0);
					if (!reti) {
						int index;
						char *tempfilename = (char*)malloc(100*sizeof(char));
						for(index=0; index<100; index++){
							tempfilename[index]='\0';
						}
						strcpy(tempfilename,de->d_name);
						for(index=strlen(tempfilename)-10; index<strlen(tempfilename); index++){
							tempfilename[index]='\0';
						}
						char* tempfilename2 = (char*)malloc(100*sizeof(char));
						for(index=0; index<100; index++){
							tempfilename2[index]='\0';
						}
						int index2=0;
						for(index=9; index<strlen(tempfilename); index++){
							tempfilename2[index2]=tempfilename[index];
							index2++;
						}
						char partition_buffer[BUFFER_SIZE];
						snprintf(partition_buffer,BUFFER_SIZE,"%s/RAxML_bestTree.%s.reroot",opt.partitions_directory,tempfilename);
						printf("%s\n",partition_buffer);
						int found = 0;
						hashmap_foreach(key,final,&mastermap){
							if ( strcmp(partition_buffer,final->filename) == 0 ){
								found = 1;
							}
						}
						if ( found == 0 ){
							printf("removing %s...\n",tempfilename);
							remove_partition_files(opt.partitions_directory, tempfilename);
						}
						free(tempfilename);
						free(tempfilename2);
					} else if ( reti != REG_NOMATCH){
						regerror(reti, &regex1, msgbuf, sizeof(msgbuf));
						fprintf(stderr, "Regex match failed: %s\n", msgbuf);
						exit(1);
					}
				}
				closedir(dr);
			}
			int index=0;
			char newick_buf[BUFFER_SIZE];
			snprintf(newick_buf,BUFFER_SIZE,"%s/tree_list.txt",opt.partitions_directory);
			FILE* list_newick_files = fopen(newick_buf,"w");
			if ( list_newick_files == NULL ){ printf("Error opening tree_list.txt file!\n"); exit(1); }
			char fp_buf2[BUFFER_SIZE];
			snprintf(fp_buf2,BUFFER_SIZE,"%s/final_partitions.txt",opt.partitions_directory);
			FILE* fp_file2 = fopen(fp_buf2,"w");
			if ( fp_file2 == NULL ){ printf("Error opening final_partitions.txt file!\n"); exit(1); }
			hashmap_foreach(key,final,&mastermap){
				fprintf(list_newick_files,"%s\n",final->filename);
				char directory[BUFFER_SIZE]; // Adjust size as needed
    				char *last_slash;
				char *filename;
				char *start_substring;
				char *end_substring;
				char substring[256];
				filename = strrchr(final->filename, '/');
				if (filename != NULL) {
					filename++; // Move past the last '/'
				} else {
					filename = final->filename; // No '/' found, use the whole path as filename
				}
				// Find the start of the substring
				start_substring = strstr(filename, "RAxML_bestTree.");
				if (start_substring != NULL) {
					start_substring += strlen("RAxML_bestTree."); // Move past "RAxML_bestTree."
					// Find the end of the substring
					end_substring = strstr(start_substring, ".reroot");
					if (end_substring != NULL) {
						size_t length = end_substring - start_substring;
						strncpy(substring, start_substring, length);
						substring[length] = '\0'; // Null-terminate the string
					} else {
						// ".reroot" not found, handle accordingly
						strcpy(substring, "");
					}
				} else {
					// "RAxML_bestTree." not found, handle accordingly
					strcpy(substring, "");
				}
				//printf("Filename: %s\n", filename);
				//printf("Substring: %s\n", substring);
    				// Find the last occurrence of '/'
    				last_slash = strrchr(final->filename, '/');
    				if (last_slash != NULL) {
        				// Calculate the length of the directory part
        				size_t directory_length = last_slash - final->filename;
        				// Copy the directory part to a new string
        				strncpy(directory, final->filename, directory_length);
        				directory[directory_length] = '\0'; // Null-terminate the string
					if ( strcmp(directory,opt.readdir)==0 ){
						char src_path[BUFFER_SIZE];
						char dst_path[BUFFER_SIZE];
						printf("copying %s/RAxML_bestTree.%s.reroot to %s...\n",opt.readdir,substring,opt.partitions_directory);
						snprintf(src_path,BUFFER_SIZE,"%s/RAxML_bestTree.%s.reroot",opt.readdir,substring);
						snprintf(dst_path,BUFFER_SIZE,"%s/RAxML_bestTree.%s.reroot",opt.partitions_directory,substring);
						copy_file(src_path, dst_path);
						printf("copying %s/%s_MSA.fasta to %s...\n",opt.readdir,substring,opt.partitions_directory);
						snprintf(src_path,BUFFER_SIZE,"%s/%s_MSA.fasta",opt.readdir,substring);
						snprintf(dst_path,BUFFER_SIZE,"%s/%s_MSA.fasta",opt.partitions_directory,substring);
						copy_file(src_path, dst_path);
						printf("copying %s/%s_taxonomy.txt to %s...\n",opt.readdir,substring,opt.partitions_directory);
						snprintf(src_path,BUFFER_SIZE,"%s/%s_taxonomy.txt",opt.readdir,substring);
						snprintf(dst_path,BUFFER_SIZE,"%s/%s_taxonomy.txt",opt.partitions_directory,substring);
						copy_file(src_path, dst_path);
					}
    				}
				/* Write partition number to final_partitions.txt */
				{
					char *pnum = strstr(substring, "partition");
					if (pnum) {
						pnum += strlen("partition");
						fprintf(fp_file2, "%s\n", pnum);
					}
				}
			}
			fclose(list_newick_files);
			fclose(fp_file2);
			exit(1);
		}
		printf("Using %d trees... \n",opt.number_of_partitions);
		//treeArr = malloc(sizeof(node *)*opt.number_of_partitions);
		SPscoreArr = malloc(sizeof(int)*MAX_NUMBEROFROOTS);
	numberOfTrees=0;
	int key;
	struct masterArr* final;
	hashmap_foreach(key,final,&mastermap){
		numberOfTrees++;
	}
	printf("Number of trees: %d\n",numberOfTrees);
	treeArr = malloc(numberOfTrees*sizeof(node*));
	taxonomyArr = (char****)malloc(numberOfTrees*sizeof(char***));
	seqArr = (int***)malloc(numberOfTrees*sizeof(int**));
	numbaseArr = (int*)malloc(numberOfTrees*sizeof(int));
	numspecArr = (int*)malloc(numberOfTrees*sizeof(int));
	rootArr = (int*)malloc(numberOfTrees*sizeof(int));
	int index=0;
	char newick_buf[BUFFER_SIZE];
	snprintf(newick_buf,BUFFER_SIZE,"%s/tree_list.txt",opt.partitions_directory);
	FILE* list_newick_files = fopen(newick_buf,"w");
	if ( list_newick_files == NULL ){ printf("Error opening list.txt file!\n"); exit(1); }
	char fp_buf[BUFFER_SIZE];
	snprintf(fp_buf,BUFFER_SIZE,"%s/final_partitions.txt",opt.partitions_directory);
	FILE* fp_file = fopen(fp_buf,"w");
	if ( fp_file == NULL ){ printf("Error opening final_partitions.txt file!\n"); exit(1); }
	hashmap_foreach(key,final,&mastermap){
		treeArr[index] = final->tree;
		taxonomyArr[index] =  final->taxonomy;
		seqArr[index] = final->msa;
		numbaseArr[index] = final->numbase;
		numspecArr[index] = final->numspec;
		rootArr[index] = final->root;
		//printtreeArr(index);
		index++;
		fprintf(list_newick_files,"%s\n",final->filename);
		/* Extract partition number from filename for final_partitions.txt */
		{
			char *fn = strrchr(final->filename, '/');
			fn = fn ? fn + 1 : final->filename;
			char *start = strstr(fn, "partition");
			if (start) {
				start += strlen("partition");
				char *end = strchr(start, '.');
				if (end) {
					fprintf(fp_file, "%.*s\n", (int)(end - start), start);
				}
			}
		}
	}
	fclose(list_newick_files);
	fclose(fp_file);

	/* Export final subtrees for ablation studies */
	if (opt.export_subtrees) {
		char export_dir[BUFFER_SIZE];
		snprintf(export_dir, BUFFER_SIZE, "%s/exported_subtrees", opt.partitions_directory);
		struct stat est = {0};
		if (stat(export_dir, &est) == -1) {
			mkdir(export_dir, 0755);
		}

		int export_idx = 0;
		struct masterArr* efinal;
		int ekey;
		hashmap_foreach(ekey, efinal, &mastermap) {
			/* Export Newick tree */
			char tree_path[BUFFER_SIZE];
			snprintf(tree_path, BUFFER_SIZE, "%s/RAxML_bestTree.%d.reroot", export_dir, export_idx);
			exportTreeToNewick(efinal, tree_path);

			/* Export unaligned FASTA */
			char fasta_path[BUFFER_SIZE];
			snprintf(fasta_path, BUFFER_SIZE, "%s/%d_MSA.fasta", export_dir, export_idx);
			FILE *efp = fopen(fasta_path, "w");
			if (efp == NULL) {
				fprintf(stderr, "Error: Cannot open %s for writing\n", fasta_path);
				exit(EXIT_FAILURE);
			}
			char *eseq = malloc((efinal->numbase + 1) * sizeof(char));
			for (int s = 0; s < efinal->numspec; s++) {
				fprintf(efp, ">%s\n", efinal->names[s]);
				removeDashArr(efinal->msa[s], efinal->numbase, eseq);
				fprintf(efp, "%s\n", eseq);
			}
			free(eseq);
			fclose(efp);

			/* Export taxonomy */
			char tax_path[BUFFER_SIZE];
			snprintf(tax_path, BUFFER_SIZE, "%s/%d_taxonomy.txt", export_dir, export_idx);
			FILE *etp = fopen(tax_path, "w");
			if (etp == NULL) {
				fprintf(stderr, "Error: Cannot open %s for writing\n", tax_path);
				exit(EXIT_FAILURE);
			}
			for (int s = 0; s < efinal->numspec; s++) {
				int ti = efinal->tree[s + efinal->numspec - 1].taxIndex[0];
				fprintf(etp, "%s\t%s;%s;%s;%s;%s;%s;%s\n", efinal->names[s],
					efinal->taxonomy[ti][6], efinal->taxonomy[ti][5],
					efinal->taxonomy[ti][4], efinal->taxonomy[ti][3],
					efinal->taxonomy[ti][2], efinal->taxonomy[ti][1],
					efinal->taxonomy[ti][0]);
			}
			fclose(etp);
			export_idx++;
		}
		printf("Exported %d subtrees to %s\n", export_idx, export_dir);
	}

	allocatetreememory_for_nucleotide_Arr(numberOfTrees);

	/* Open progress file for external monitoring */
	char progress_path[BUFFER_SIZE];
	snprintf(progress_path, BUFFER_SIZE, "%s/_progress.txt", opt.partitions_directory);
	FILE *progress_fp = fopen(progress_path, "w");
	if (progress_fp) {
		fprintf(progress_fp, "stage=posteriors\ntrees_total=%d\ntrees_done=0\n", numberOfTrees);
		fflush(progress_fp);
	}

	printf("Computing posteriors for %d trees...\n", numberOfTrees);
	fflush(stdout);
	int trees_done = 0;
	#pragma omp parallel for schedule(dynamic) private(i)
	for(i=0; i<numberOfTrees; i++){
		double local_params[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
		estimatenucparameters_Arr(local_params,i);
		getposterior_nc_Arr(local_params,i);
		int done;
		#pragma omp atomic capture
		done = ++trees_done;
		if (done == numberOfTrees || done % 10 == 0){
			printf("  Posteriors: %d/%d trees done\n", done, numberOfTrees);
			fflush(stdout);
			if (progress_fp) {
				fseek(progress_fp, 0, SEEK_SET);
				ftruncate(fileno(progress_fp), 0);
				fprintf(progress_fp, "stage=posteriors\ntrees_total=%d\ntrees_done=%d\n", numberOfTrees, done);
				fflush(progress_fp);
			}
		}
	}
	printf("Posterior computation complete.\n");
	if (progress_fp) {
		fseek(progress_fp, 0, SEEK_SET);
		ftruncate(fileno(progress_fp), 0);
		fprintf(progress_fp, "stage=done\ntrees_total=%d\ntrees_done=%d\n", numberOfTrees, numberOfTrees);
		fflush(progress_fp);
		fclose(progress_fp);
	}
	//FILE *for_monica = fopen("/space/s1/lenore/trout_copy/s2_copy/for_monica/OU061397_1_PP.txt","w");
	//if ( for_monica == NULL ){ printf("Error opening file!\n"); exit(1); }
	//for(i=0; i<numbaseArr[0]; i++){
	//	fprintf(for_monica,"%d",i);
	//	for(j=0; j<4; j++){
	//		fprintf(for_monica,"\t%.17g",treeArr[0][200].posteriornc[i][j]);
	//	}
	//	fprintf(for_monica,"\n");
	//}
	//fclose(for_monica);
	//exit(1);
	if (opt.missing_data==1){
		printf("Adjusting posteriors for missing data (%d trees)...\n", numberOfTrees);
		fflush(stdout);
		#pragma omp parallel for schedule(dynamic) private(i, j)
		for(i=0; i<numberOfTrees; i++){
			changePP_Arr(rootArr[i],i);
			for(j=numspecArr[i]-1;j<2*numspecArr[i]-1;j++){
				changePP_parents_Arr(j,i);
			}
		}
	}
	printf("Missing data adjustment complete.\n");
	fflush(stdout);
	printf("Writing reference_tree.txt...\n");
	fflush(stdout);
	printTreeFile(numberOfTrees,max_nodename,max_tax_name,max_lineTaxonomy,opt);
	printf("Done.\n");
	/*hashmap_foreach(key,final,&mastermap){
		for(i=0; i<final->numspec; i++){
			free(final->names[i]);
			for(j=0; j<7; j++){
				free(final->taxonomy[i][j]);
			}
			free(final->taxonomy[i]);
			free(final->msa[i]);
		}
		free(final->taxonomy);
		free(final->names);
		free(final->msa);
	}*/
	for(i=0; i<numberOfTrees; i++){
		freeTreeMemory(i);
		for(j=0; j<numspecArr[i]; j++){
			free(seqArr[i][j]);
			for(k=0; k<7; k++){
				free(taxonomyArr[i][j][k]);
			}
			free(taxonomyArr[i][j]);
		}
		free(taxonomyArr[i]);
		free(seqArr[i]);
	}
	free(treeArr);
	//hashmap_cleanup(&mastermap);
	free(taxonomyArr);
	free(seqArr);
	free(numbaseArr);
	free(numspecArr);
	free(rootArr);
}
