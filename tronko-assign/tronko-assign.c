#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "global.h"
#include "placement.h"
#include "needleman_wunsch.h"
#include "readreference.h"
#include "options.h"
#include "printAlignments.h"
#include "bwa_source_files_include.h"
#include "hashmap.h"
#include "allocateMemoryForResults.h"
#include "WFA2/wavefront_align.h"
#include "hashmap_base.h"
#include "logger.h"
#include "resource_monitor.h"
#include "crash_debug.h"
#include "symbol_resolver.h"
#include "tsv_memlog.h"
#ifdef ENABLE_PARQUET
#include "parquet_writer.h"

/*
 * Parse a taxonPath TSV string into an assignmentResult struct
 * Format: Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number
 * Returns 0 on success, -1 on parse error
 */
static int parse_tsv_to_result(const char *tsv_line, assignmentResult *result) {
    if (!tsv_line || !result) return -1;

    /* Make a copy since strtok modifies the string */
    char *line_copy = strdup(tsv_line);
    if (!line_copy) return -1;

    char *saveptr = NULL;
    char *token;
    int field = 0;

    token = strtok_r(line_copy, "\t", &saveptr);
    while (token && field < 7) {
        switch (field) {
            case 0: /* Readname */
                result->readname = strdup(token);
                break;
            case 1: /* Taxonomic_Path */
                result->taxonomic_path = strdup(token);
                break;
            case 2: /* Score */
                result->score = atof(token);
                break;
            case 3: /* Forward_Mismatch */
                result->forward_mismatch = atof(token);
                break;
            case 4: /* Reverse_Mismatch */
                result->reverse_mismatch = atof(token);
                break;
            case 5: /* Tree_Number */
                result->tree_number = atoi(token);
                break;
            case 6: /* Node_Number */
                result->node_number = atoi(token);
                break;
        }
        field++;
        token = strtok_r(NULL, "\t", &saveptr);
    }

    free(line_copy);

    /* Check if we got at least readname and taxonomic_path */
    if (field < 2) return -1;

    /* Handle unassigned case - may only have 2 fields */
    if (field == 2) {
        result->score = 0.0;
        result->forward_mismatch = 0.0;
        result->reverse_mismatch = 0.0;
        result->tree_number = -1;
        result->node_number = -1;
    }

    return 0;
}

static void free_assignment_result(assignmentResult *result) {
    if (result) {
        free(result->readname);
        free(result->taxonomic_path);
        result->readname = NULL;
        result->taxonomic_path = NULL;
    }
}
#endif

// Thread-local counter for dropped BWA matches (reset per read)
static __thread int dropped_matches_count = 0;

// Global diagnostic counters for measuring cap impact
static int g_overflow_read_count = 0;        // Reads that hit the cap
static int g_total_dropped_matches = 0;      // Total matches dropped across all reads
static int g_max_potential_matches = 0;      // Highest potential match count seen
static pthread_mutex_t g_overflow_stats_mutex = PTHREAD_MUTEX_INITIALIZER;

//int numspec, numbase, root/***seq, numundspec[MAXNUMBEROFINDINSPECIES+1]*/;
// Global TSV log file for use by readreference.c
FILE *g_tsv_log_file = NULL;
char ****taxonomyArr;
struct node **treeArr;
struct queryMatPaired *pairedQueryMat;
struct queryMatSingle *singleQueryMat;
type_of_PP Cinterval;
//struct hashmap_base map;
//struct map;
int *numspecArr, *numbaseArr, *rootArr;
struct timespec tstart={0,0}, tend={0.0};

char* strupr(char* s){
	unsigned char* t = (unsigned char*) s;
	while( *t){
		*t = toupper(*t);
		t++;
	}
	return s;
}
void store_PPs_Arr(int numberOfRoots, double c){
	int i, j, k, l;
	for(i=0; i<numberOfRoots; i++){
		for (j=0; j<2*numspecArr[i]-1; j++){
			for (k=0; k<numbaseArr[i]; k++){
				for (l=0; l<4; l++){
					if ( treeArr[i][j].posteriornc[PP_IDX(k, l)] == -1 ){   //Missing data
						treeArr[i][j].posteriornc[PP_IDX(k, l)]=1;
					}else{
							//treeArr[i][j].posteriornc[PP_IDX(k, l)] = log((c/3) + ((1-(((4*c)/3) * treeArr[i][j].posteriornc[PP_IDX(k, l)]))));
							double d = 1-c;
							double e = c/3;
							double f = d * treeArr[i][j].posteriornc[PP_IDX(k, l)];
							double g = e * (1-treeArr[i][j].posteriornc[PP_IDX(k, l)]);
							treeArr[i][j].posteriornc[PP_IDX(k, l)] = log( (f + g) );
					}
				}
			}
		}
	}
}
void assignDepthArr(int node0, int node1, int depth, int whichPartitions){
	if( node0 != -1 && node1 != -1){
		treeArr[whichPartitions][node0].depth = depth;
		treeArr[whichPartitions][node1].depth = depth;
		assignDepthArr(treeArr[whichPartitions][node0].up[0], treeArr[whichPartitions][node0].up[1],depth+1,whichPartitions);
		assignDepthArr(treeArr[whichPartitions][node1].up[0], treeArr[whichPartitions][node1].up[1],depth+1,whichPartitions);
	}
}
void printTreeInfo(int whichPartition, int node, FILE* file){
	int i,j;
	if (treeArr[whichPartition][node].up[0]==-1 && treeArr[whichPartition][node].up[1]==-1){
		fprintf(file,"%s\t%d\t%d\n",treeArr[whichPartition][node].name,whichPartition,node);
		return;
	}else{
		printTreeInfo(whichPartition,treeArr[whichPartition][node].up[0],file);
		printTreeInfo(whichPartition,treeArr[whichPartition][node].up[1],file);
	}
}
void run_bwa(int start, int end, bwaMatches* bwa_results, int concordant, int numberOfTrees, char *databasefile, int paired, int max_query_length, int max_readname_length, int max_acc_name){
	int i,j;
	int number_of_threads=1;
	if (paired != 0){
		main_mem(databasefile,end-start,number_of_threads, bwa_results, concordant, numberOfTrees, start, paired, start, end, max_query_length, max_readname_length, max_acc_name);
	}else{
		main_mem(databasefile,end-start,number_of_threads, bwa_results, concordant, numberOfTrees, start, paired, start, end, max_query_length, max_readname_length, max_acc_name);
	}
}
int getLCA_Arr(int node1, int node2, int whichRoot){
	if (node1 == node2){ return node1; }
	if (treeArr[whichRoot][node1].depth > treeArr[whichRoot][node2].depth){
		int tmp = node1;
		node1 = node2;
		node2 = tmp;
	}
	node2 = treeArr[whichRoot][node2].down;
	return getLCA_Arr(node1,node2,whichRoot);
}
int getKeysCount(int whichRoot, int node, int* minNodes, int matching_nodes, int* ancestors, int numMinNodes){
	int child0 = treeArr[whichRoot][node].up[0];
	int child1 = treeArr[whichRoot][node].up[1];
	if ( child0 != -1 && child1 != -1 ){
		matching_nodes += getKeysCount(whichRoot,child0,minNodes,matching_nodes,ancestors,numMinNodes) + getKeysCount(whichRoot,child1,minNodes,matching_nodes,ancestors,numMinNodes);
	}
	int i;
	for(i=0; i<numMinNodes; i++){
		if ( minNodes[i] == node ){
			matching_nodes++;
		}
	}
	if ( matching_nodes == numMinNodes ){
		for(i=0; i<2*numspecArr[whichRoot]-1; i++){
			if ( ancestors[i] == -1 ){
				break;
			}
		}
		ancestors[i] = node;
	}
	return matching_nodes;
}
int LCA_of_nodes(int whichRoot, int root_node, int* minNodes, int numMinNodes){
	int* ancestors = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
	int i;
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		ancestors[i] = -1;
	}
	int matching_nodes = 0;
	getKeysCount(whichRoot, root_node, minNodes, matching_nodes, ancestors, numMinNodes);
	int LCA = ancestors[0];
	free(ancestors);
	return LCA;
}
int getLCAofArray_Arr(int *minNodes,int whichRoot,int maxNumSpec, int number_of_total_nodes){
	int LCA = minNodes[0];
	int maxDepth = 1000000000;
	int i;
	for( i=1; i<number_of_total_nodes; i++){
		if ( minNodes[i] == -1 ){ return LCA; }
		if (treeArr[whichRoot][minNodes[i]].depth < maxDepth){
			LCA = getLCA_Arr(LCA, minNodes[i], whichRoot);
		}
	}
	return LCA;
}
int getLCAofArray_Arr_Multiple(int *voteroot,int whichRoot, int maxNumSpec, int number_of_total_nodes){
	/*int LCA = minNodes[0];
	int maxDepth = 1000000000;
	int i;
	for( i=1; i<number_of_total_nodes; i++){
		if ( minNodes[i] == 0 ){ return LCA; }
		if (treeArr[whichRoot][minNodes[i]].depth < maxDepth){
			LCA = getLCA_Arr(LCA, minNodes[i], whichRoot);
		}
	}
	return LCA;*/
	int i;
	int* minNodes = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		minNodes[i]=-1;
	}
	int count=0;
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		if ( voteroot[i] == 1 ){
			minNodes[count]=i;
			count++;
		}
	}
	int LCA = LCA_of_nodes(whichRoot,rootArr[whichRoot],minNodes,count);
	free(minNodes);
	return LCA;
}
void *runAssignmentOnChunk_WithBWA(void *ptr){
	struct mystruct *mstr = (mystruct *) ptr;
	resultsStruct *results=mstr->str;
	char **rootSeqs=mstr->rootSeqs;
	int numberOfTrees = mstr->ntree;
	char query_1[mstr->max_query_length];
	char query_2[mstr->max_query_length];
	int maxNumSpec = mstr->maxNumSpec;
	int iter = 0;
	int i,j,k,lineNumber;
	int end=mstr->end;
	int *minNodes=results->minNodes;
	char **LCAnames = results->LCAnames;
	int paired = mstr->paired;
	int print_alignments = results->print_alignments;
	int use_nw = mstr->use_nw;
	int print_alignments_to_file = mstr->print_alignments_to_file;
	int print_unassigned = mstr->print_unassigned;
	int use_leaf_portion = mstr->use_leaf_portion;
	int padding = mstr->padding;
	int max_query_length = mstr->max_query_length;
	int max_readname_length = mstr->max_readname_length;
	int max_acc_name = mstr->max_acc_name;
	int max_numbase = mstr->max_numbase;
	int number_of_total_nodes = mstr->number_of_total_nodes;
	int print_all_nodes = mstr->print_all_nodes;
	// Tier 1 optimization settings
	int early_termination = mstr->early_termination;
	type_of_PP strike_box = mstr->strike_box;
	int max_strikes = mstr->max_strikes;
	int enable_pruning = mstr->enable_pruning;
	type_of_PP pruning_factor = mstr->pruning_factor;
	/*affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch = 4,
		.gap_opening = 6,
		.gap_extension = 2,
	};*/
	/*mm_allocator_t* mm_allocator;
	affine_wavefronts_t* affine_wavefronts;
	char* pattern_alg;
	char* text_alg;*/
	/*if (use_nw==0){
		if (use_leaf_portion==1){
			mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
			pattern_alg=mm_allocator_calloc(mm_allocator,max_query_length+max_query_length+padding+50+1,char,true);
			text_alg=mm_allocator_calloc(mm_allocator,max_query_length+max_query_length+padding+50+1,char,true);
			affine_wavefronts = affine_wavefronts_new_complete(max_query_length+max_query_length+padding+50+1,max_query_length+max_query_length+padding+50+1,&affine_penalties,NULL,mm_allocator);
		}
	}*/
	char* resultsPath = (char*)malloc((max_readname_length+mstr->max_lineTaxonomy+120)*sizeof(char));
	bwaMatches* bwa_results = (bwaMatches *)malloc((end-(mstr->start))*sizeof(bwaMatches));
	for (i=0; i<end-mstr->start; i++){
		//bwa_results[i].readname = (char*)malloc((max_readname_length+1)*sizeof(char));
		bwa_results[i].concordant_matches_roots = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].concordant_matches_nodes = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].discordant_matches_roots = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].discordant_matches_nodes = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].use_portion = use_leaf_portion;
		if ( use_leaf_portion == 1 ){
			bwa_results[i].cigars_forward = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
			bwa_results[i].starts_forward = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
			if (paired == 1){
				bwa_results[i].cigars_reverse = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
				bwa_results[i].starts_reverse = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
			}
		}
		for(j=0; j<MAX_NUM_BWA_MATCHES; j++){
			bwa_results[i].discordant_matches_roots[j] = -1;
			bwa_results[i].discordant_matches_nodes[j] = -1;
			bwa_results[i].concordant_matches_roots[j] = -1;
			bwa_results[i].concordant_matches_nodes[j] = -1;
			if (use_leaf_portion == 1){
				bwa_results[i].starts_forward[j] = -1;
				if(paired==1){
					bwa_results[i].starts_reverse[j] = -1;
				}
			}
		}
		for(j=0; j<MAX_NUM_BWA_MATCHES; j++){
			if (use_leaf_portion == 1){
				bwa_results[i].cigars_forward[j] = (char *)malloc(MAX_CIGAR*sizeof(char));
				memset(bwa_results[i].cigars_forward[j],'\0',MAX_CIGAR);
				if (paired==1){
					bwa_results[i].cigars_reverse[j] = (char *)malloc(MAX_CIGAR*sizeof(char));
					memset(bwa_results[i].cigars_reverse[j],'\0',MAX_CIGAR);
				}
			}
		}
		bwa_results[i].n_matches = 0;
	}
	int trees_search[mstr->ntree];
	int *leaf_coord_arr;
	char *leaf_sequence;
	int *positionsInRoot;
	if (use_leaf_portion == 1){
		leaf_sequence = (char *)malloc((max_query_length+max_query_length+2*padding+1)*sizeof(char));
		positionsInRoot = (int *)malloc((max_query_length+max_query_length+2*padding+1)*sizeof(int));
	}else{
		leaf_sequence = (char *)malloc((max_query_length+mstr->max_numbase+1)*sizeof(char));
		positionsInRoot = (int *)malloc((max_query_length+mstr->max_numbase+1)*sizeof(int));
	}
	struct leafMap *leaf_map;	
	run_bwa(mstr->start, end, bwa_results, mstr->concordant, mstr->ntree, mstr->databasefile, paired, max_query_length, max_readname_length, max_acc_name);
	for ( lineNumber=mstr->start; lineNumber<end; lineNumber++){
		// Set crash context for current read processing
		char read_info[64];
		snprintf(read_info, sizeof(read_info), "read_%d", lineNumber);
		crash_set_current_read(read_info, 1, lineNumber);
		crash_set_processing_stage("BWA alignment and tree search");
		
		for(i=0; i<mstr->ntree; i++){
			trees_search[i]=-1;
		}
		j=0;
		int hashValue;
		int no_add=0;
		int leaf_iter=0;
		dropped_matches_count = 0;  // Reset dropped counter for this read
		if (bwa_results[iter].concordant_matches_roots[0]==-1 && mstr->concordant==1){
			for (i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].discordant_matches_roots[0] < 0 ){
					//for(j=0; j<mstr->ntree;j++){
					//	results->leaf_coordinates[j][0]=j;
					//	results->leaf_coordinates[j][1]=rootArr[j];
					//}
					//leaf_iter=mstr->ntree;
					//i=mstr->ntree;
					break;
				}else if ( bwa_results[iter].discordant_matches_roots[i]==-1){
					break;
				}else{
					if (leaf_iter < MAX_NUM_BWA_MATCHES) {
						results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
						results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
						if (use_leaf_portion==1){
							results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
							strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
							if ( paired==1){
								results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
								strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
							}
						}
					}
				}
				int index1=mstr->ntree-1;
				for(k=mstr->ntree-1; k>=0; k--){
					if (trees_search[k]==-1){
						index1=k;
					}
				}
				int found=0;
				for(k=0; k<index1; k++){
					if (leaf_iter < MAX_NUM_BWA_MATCHES && trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
						found=1;
					}
				}
				if (found==0){
					if (leaf_iter < MAX_NUM_BWA_MATCHES) {
						trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
						leaf_iter++;
					} else {
						dropped_matches_count++;
					}
				}
			}
		}else if (mstr->concordant==1){
			for(i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].concordant_matches_roots[i]==-1){
				//if (strlen(bwa_results[iter].concordant_leaf_matches[i])<=3){
					break;
				}else{
					if (leaf_iter < MAX_NUM_BWA_MATCHES) {
						results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].concordant_matches_roots[i];
						results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].concordant_matches_nodes[i];
						if (use_leaf_portion==1){
							results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
							strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
							if ( paired==1){
								results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
								strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
							}
						}
					}
					int index1=mstr->ntree-1;
					for(k=mstr->ntree-1; k>=0; k--){
						if (trees_search[k]==-1){
							index1=k;
						}
					}
					int found=0;
					for(k=0; k<index1; k++){
						if (leaf_iter < MAX_NUM_BWA_MATCHES && trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
							found=1;
						}
					}
					if (found==0){
						if (leaf_iter < MAX_NUM_BWA_MATCHES) {
							trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
							leaf_iter++;
						} else {
							dropped_matches_count++;
						}
					}
				}
			}
		}else{
			j=0;
			for(i=0; i<mstr->ntree; i++){
				//if(strlen(bwa_results[iter].discordant_leaf_matches[i])<=3){
				if(bwa_results[iter].discordant_matches_roots[i]==-1){
					break;
				}else{
					if (no_add==0){
						if (leaf_iter < MAX_NUM_BWA_MATCHES) {
							results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
							results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
							if (use_leaf_portion==1){
								results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
								strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
								if (paired==1){
									results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
									strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
								}
							}
						}
						int index1=mstr->ntree-1;
						for(k=mstr->ntree-1; k>=0; k--){
							if (trees_search[k]==-1){
								index1=k;
							}
						}
						int found=0;
						for(k=0; k<index1; k++){
							if (leaf_iter < MAX_NUM_BWA_MATCHES && trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
								found=1;
							}
						}
						if (found==0){
							if (leaf_iter < MAX_NUM_BWA_MATCHES) {
								trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
								leaf_iter++;
							} else {
								dropped_matches_count++;
							}
						}
						j++;
					}
					no_add=0;
				}
			}
			for(i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].concordant_matches_roots[i]==-1){
				//if (strlen(bwa_results[iter].concordant_leaf_matches[i])<=3){
					break;
				}else{
					if (no_add==0){
						if (leaf_iter < MAX_NUM_BWA_MATCHES) {
							results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
							results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
							if (use_leaf_portion == 1){
								results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
								strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
								if (paired==1){
									results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
									strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
								}
							}
						}
					int index1=mstr->ntree-1;
					for(k=mstr->ntree-1; k>=0; k--){
						if (trees_search[k]==-1){
							index1=k;
						}
					}
					int found=0;
					for(k=0; k<index1; k++){
						if (leaf_iter < MAX_NUM_BWA_MATCHES && trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
							found=1;
						}
					}
					if (found==0){
						if (leaf_iter < MAX_NUM_BWA_MATCHES) {
							trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
							leaf_iter++;
						} else {
							dropped_matches_count++;
						}
					}
						j++;
					}
					no_add=0;
				}
			}
		}

		// Update crash context and log if matches were dropped
		if (dropped_matches_count > 0) {
			int potential_matches = leaf_iter + dropped_matches_count;
			const char *read_name = paired ? pairedQueryMat->forward_name[lineNumber] : singleQueryMat->name[lineNumber];

			// Update global statistics (thread-safe)
			pthread_mutex_lock(&g_overflow_stats_mutex);
			g_overflow_read_count++;
			g_total_dropped_matches += dropped_matches_count;
			if (potential_matches > g_max_potential_matches) {
				g_max_potential_matches = potential_matches;
			}
			int local_overflow_count = g_overflow_read_count;
			pthread_mutex_unlock(&g_overflow_stats_mutex);

			crash_set_bwa_bounds_violation(leaf_iter, MAX_NUM_BWA_MATCHES, dropped_matches_count);

			// Log first 100 occurrences for debugging, then every 1000th
			if (local_overflow_count <= 100 || local_overflow_count % 1000 == 0) {
				LOG_WARN("Read %s: %d unique tree matches (capped at %d, dropped %d) [overflow #%d]",
				         read_name, potential_matches, MAX_NUM_BWA_MATCHES,
				         dropped_matches_count, local_overflow_count);
			}
		}

		results->minimum[0]=0;
		if (leaf_iter > 0 ){
			//strcmp(results->cigars_forward[0],"*")!=0 && strcmp(results->cigars_forward[0]," ")!=0){
		if (paired != 0){
			if ( use_nw==0 ){
				place_paired(pairedQueryMat->query1Mat[lineNumber],pairedQueryMat->query2Mat[lineNumber],rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,pairedQueryMat->forward_name[lineNumber],pairedQueryMat->reverse_name[lineNumber],print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes,early_termination,strike_box,max_strikes,enable_pruning,pruning_factor);
			}else{
				place_paired_with_nw(pairedQueryMat->query1Mat[lineNumber],pairedQueryMat->query2Mat[lineNumber],rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nw,results->aln,results->scoring,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,pairedQueryMat->forward_name[lineNumber],pairedQueryMat->reverse_name[lineNumber],print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes,early_termination,strike_box,max_strikes,enable_pruning,pruning_factor);
			}
		}else{
			if (use_nw==0){
				place_paired(singleQueryMat->queryMat[lineNumber],NULL,rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,singleQueryMat->name[lineNumber],NULL,print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes,early_termination,strike_box,max_strikes,enable_pruning,pruning_factor);
			}else{
				place_paired_with_nw(singleQueryMat->queryMat[lineNumber],NULL,rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nw,results->aln,results->scoring,results->nodeScores,results->voteRoot,leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,singleQueryMat->name[lineNumber],NULL,print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes,early_termination,strike_box,max_strikes,enable_pruning,pruning_factor);
			}
		}
		numberOfTrees = leaf_iter;
		memset(minNodes, -1, number_of_total_nodes * sizeof(int));
		for(i=0; i<leaf_iter; i++){
			memset(results->nodeScores[i][trees_search[i]], 0, (2*numspecArr[trees_search[i]]-1) * sizeof(type_of_PP));
			results->leaf_coordinates[i][0]=-1;
			results->leaf_coordinates[i][1]=-1;
		}
			if (use_leaf_portion==1){
				results->starts_forward[leaf_iter] = -1;
				for(j=0; j<MAX_CIGAR; j++){
					results->cigars_forward[leaf_iter][j]='\0';
				}
				if(paired==1){
					results->starts_reverse[leaf_iter]=-1;
					for(j=0; j<MAX_CIGAR; j++){
						results->cigars_reverse[leaf_iter][j]='\0';
					}
				}
			}
		}
		int countVotes[mstr->ntree];
		int count=0;
		for(i=0; i<mstr->ntree; i++){
			countVotes[i]=0;
			for(j=0;j<2*numspecArr[i]-1;j++){
				if (results->voteRoot[i][j]==1){
					countVotes[i]++;
					minNodes[count]=j;
					count++;
				}
			}
		}
		int numMinNodes=count;
		int max=0;
		int maxRoot=-1;
		count=0;
		for(i=0;i<mstr->ntree;i++){
			if (countVotes[i]>max){
				max=countVotes[i];
				maxRoot=i;
			}
			if (countVotes[i]>0){
				count++;
			}
		}
		int LCA, LCAs[count], maxRoots[count];
		int count2=0;
		for(i=0;i<mstr->ntree;i++){
			if ( countVotes[i]>0){
				maxRoots[count2]=i;
				count2++;
			}
		}
		int unassigned=0;
		int minLevel=0;
		int taxRoot,taxIndex0,taxIndex1,taxNode;
		if ( count == 1 ){
			// Set context for tree processing
			crash_set_current_tree(maxRoot);
			crash_set_processing_stage("LCA calculation and taxonomic assignment");
			clock_gettime(CLOCK_MONOTONIC, &tstart);
			//LCA=getLCAofArray_Arr(minNodes,maxRoot,maxNumSpec,number_of_total_nodes);
			LCA = LCA_of_nodes(maxRoot,rootArr[maxRoot],minNodes,numMinNodes);
		}else if (count != 0){
			for(i=0;i<count;i++){
				for(j=0; j<mstr->max_lineTaxonomy; j++){
					LCAnames[i][j]='\0';
				}
			}
			for(i=0;i<count;i++){
				LCAs[i]=getLCAofArray_Arr_Multiple(results->voteRoot[maxRoots[i]],maxRoots[i],maxNumSpec,number_of_total_nodes);
				if ( treeArr[maxRoots[i]][LCAs[i]].taxIndex[1]!=-1){ 
					strcpy(LCAnames[i],taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[1]]);
					if ( treeArr[maxRoots[i]][LCAs[i]].taxIndex[1] > minLevel ){
						minLevel=treeArr[maxRoots[i]][LCAs[i]].taxIndex[1];
					}
				}else{
					unassigned=1;
				}
			}
			LCA = LCAs[0];
			int correctTax=0;
			int stop=0;
			while(stop==0 && minLevel<=6){
			for(i=0; i<count;i++){
				if (treeArr[maxRoots[i]][LCAs[i]].taxIndex[0] != -1){
				for(j=i+1; j<count; j++){
					if ( treeArr[maxRoots[j]][LCAs[j]].taxIndex[0] != -1){
					if ( strcmp(taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][minLevel],taxonomyArr[maxRoots[j]][treeArr[maxRoots[j]][LCAs[j]].taxIndex[0]][minLevel])==0 && taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][minLevel] != "NA" ){
						correctTax++;
					}
					}
				}
				}
			}
			if (correctTax>=count-1){
				stop=1;
				taxRoot=maxRoots[0];
				taxNode=LCAs[0];
				taxIndex0=treeArr[maxRoots[0]][LCAs[0]].taxIndex[0];
				taxIndex1=minLevel;
			}else{
				minLevel++;
			}
			}
			if (minLevel>6){ unassigned=1;}
		}
		int resultsPathSize = max_readname_length+mstr->max_lineTaxonomy+120;
		resultsPath[0] = '\0';
		int print_un=1;
		if (count==1){
			const char *readname = (paired != 0) ?
				pairedQueryMat->forward_name[lineNumber] :
				singleQueryMat->name[lineNumber];
			int taxIndex1 = treeArr[maxRoot][LCA].taxIndex[1];

			/* Build taxonomy string into resultsPath using cursor */
			int pos = 0;
			pos += snprintf(resultsPath + pos, resultsPathSize - pos, "%s\t", readname);
			if ( taxIndex1 == -1 ){
				pos += snprintf(resultsPath + pos, resultsPathSize - pos, "unassigned Euk or Bac");
			}else{
				for(i=6;i>=taxIndex1;i--){
					pos += snprintf(resultsPath + pos, resultsPathSize - pos, "%s%s",
						taxonomyArr[maxRoot][treeArr[maxRoot][LCA].taxIndex[0]][i],
						(i==taxIndex1) ? "" : ";");
				}
			}
			pos += snprintf(resultsPath + pos, resultsPathSize - pos,
				"\t%lf\t%lf\t%lf\t%d\t%d",
				results->minimum[0], results->minimum[1], results->minimum[2],
				maxRoot, LCA);
			strcpy(results->taxonPath[iter], resultsPath);
		}else if (count==0){
			const char *readname = (paired != 0) ?
				pairedQueryMat->forward_name[lineNumber] :
				singleQueryMat->name[lineNumber];
			snprintf(resultsPath, resultsPathSize, "%s\tunassigned", readname);
			strcpy(results->taxonPath[iter],resultsPath);
		}else{
			const char *readname = (paired != 0) ?
				pairedQueryMat->forward_name[lineNumber] :
				singleQueryMat->name[lineNumber];
			int pos = 0;
			pos += snprintf(resultsPath + pos, resultsPathSize - pos, "%s\t", readname);
			if (unassigned==1 ){
				pos += snprintf(resultsPath + pos, resultsPathSize - pos, "unassigned Eukaryote or Bacteria");
			}else{
				int taxIndex1=minLevel;
				for(i=6;i>=taxIndex1;i--){
					pos += snprintf(resultsPath + pos, resultsPathSize - pos, "%s%s",
						taxonomyArr[taxRoot][treeArr[taxRoot][taxNode].taxIndex[0]][i],
						(i==taxIndex1) ? "" : ";");
				}
			}
			pos += snprintf(resultsPath + pos, resultsPathSize - pos,
				"\t%lf\t%lf\t%lf\t%d\t%d",
				results->minimum[0], results->minimum[1], results->minimum[2],
				maxRoot, LCA);
			strcpy(results->taxonPath[iter],resultsPath);
		}
		//free(bwa_results[iter].concordant_matches);
		//free(bwa_results[iter].discordant_matches);
		if (use_leaf_portion == 1){
			free(bwa_results[iter].starts_forward);
			if(paired==1){
				free(bwa_results[iter].starts_reverse);
			}
		}
		for(i=0; i<MAX_NUM_BWA_MATCHES;i++){
			//free(bwa_results[iter].discordant_leaf_matches);
			//free(bwa_results[iter].concordant_leaf_matches);
			if (use_leaf_portion == 1){
				free(bwa_results[iter].cigars_forward[i]);
				if(paired==1){
						free(bwa_results[iter].cigars_reverse[i]);
				}
			}
		}
		//free(bwa_results[iter].discordant_leaf_matches);
		//free(bwa_results[iter].concordant_leaf_matches);
		if (use_leaf_portion == 1){
			free(bwa_results[iter].cigars_forward);
			if (paired==1){
				free(bwa_results[iter].cigars_reverse);
			}
		}
		free(bwa_results[iter].concordant_matches_roots);
		free(bwa_results[iter].concordant_matches_nodes);
		free(bwa_results[iter].discordant_matches_roots);
		free(bwa_results[iter].discordant_matches_nodes);
		//free(bwa_results[iter].readname);
		crash_clear_bwa_context();  // Clear BWA context at end of read processing
		iter++;
		results->minimum[0] = -1;
		results->minimum[1] = -1;
		results->minimum[2] = -1;
		LCA = -1;
		if (leaf_iter > 0){
			for(i=0; i<mstr->ntree; i++){
				memset(results->voteRoot[i], 0, (2*numspecArr[i]-1) * sizeof(int));
			}
		}
	}
	/*if (use_nw == 0){
		affine_wavefronts_delete(affine_wavefronts);
		mm_allocator_free(mm_allocator,pattern_alg);
		mm_allocator_free(mm_allocator,text_alg);
		mm_allocator_delete(mm_allocator);
	}*/
	/*if (use_nw == 0){
		mm_allocator_delete(mm_allocator);
	}*/
	free(leaf_sequence);
	free(positionsInRoot);
	free(bwa_results);
	free(resultsPath);
	pthread_exit(NULL);
}
int main(int argc, char **argv){
	int i, j, k, numberOfTrees;
	Options opt;
	opt.use_nw=0;
	opt.reference_mode=1; // default use reference
	opt.reverse_single_read=0; // default no reverse single read
	opt.reverse_second_of_paired_read=0;
	opt.print_alignments = 0;
	opt.print_node_info[0] = '\0';
	opt.results_file[0] = '\0';
	opt.reference_file[0] = '\0';
	opt.print_trees_dir[0] = '\0';
	opt.fastq=0; //default is FASTA
	opt.unassigned=0; //don't print unassigned sequences
	opt.print_alignments_to_file=0; //don't print alignments to file
	opt.use_leaf_portion=0;
	opt.padding=0;
	opt.skip_build=0;
	opt.number_of_cores=1;
	opt.number_of_lines_to_read=50000;
	opt.score_constant = 0.01;
	opt.print_all_nodes=0;
	opt.verbose_level = -1;  // Disabled by default
	opt.log_file[0] = '\0';
	opt.enable_resource_monitoring = 0;
	opt.enable_timing = 0;
	opt.tsv_log_file[0] = '\0';
	// Tier 1 optimization defaults (disabled by default for baseline comparison)
	opt.early_termination = 0;
	opt.strike_box = 1.0;
	opt.max_strikes = 6;
	opt.enable_pruning = 0;
	opt.pruning_factor = 2.0;

	parse_options(argc, argv, &opt);
	
	// Initialize logging based on options
	if (opt.verbose_level >= 0) {
		log_level_t level = (log_level_t)(3 - opt.verbose_level);  // Convert to our enum (3=DEBUG, 2=INFO, etc.)
		const char* log_file = (opt.log_file[0] != '\0') ? opt.log_file : NULL;
		logger_init(level, log_file, 1);  // Always log to stderr
		
		if (opt.enable_resource_monitoring) {
			logger_enable_resource_monitoring(1);
			init_resource_monitoring();
		}
		
		if (opt.enable_timing) {
			logger_enable_timing(1);
		}
		
		// Initialize crash debugging system
		crash_clear_context();  // Initialize context tracking
		crash_config_t crash_config = {0};
		crash_config.enable_stack_trace = 1;
		crash_config.enable_register_dump = 1;
		crash_config.enable_memory_dump = 1;
		crash_config.enable_core_dump = 1;
		crash_config.enable_crash_log = 1;
		strcpy(crash_config.crash_log_directory, "/tmp");
		strcpy(crash_config.crash_log_prefix, "tronko_assign_crash");
		
		if (crash_debug_init(&crash_config) == 0) {
			LOG_DEBUG("Crash debugging system initialized");
		}
		
		// Initialize symbol resolver for enhanced crash reports
		symbol_resolver_config_t symbol_config = {0};
		symbol_config.enable_addr2line = 1;
		symbol_config.enable_source_lookup = 1;
		symbol_config.cache_symbols = 1;
		strcpy(symbol_config.debug_info_path, "/usr/lib/debug");
		strcpy(symbol_config.addr2line_path, "addr2line");
		
		if (symbol_resolver_init(&symbol_config) == 0) {
			LOG_DEBUG("Enhanced symbol resolution initialized");
		}
		
		// Log program startup (now that logging is initialized)
		LOG_MILESTONE_TIMED(MILESTONE_STARTUP);
		
		char milestone_info[256];
		snprintf(milestone_info, sizeof(milestone_info), 
			"Verbose logging enabled, level=%d, cores=%d, lines=%d", 
			opt.verbose_level, opt.number_of_cores, opt.number_of_lines_to_read);
		log_milestone_with_timing(MILESTONE_OPTIONS_PARSED, milestone_info);
		
		LOG_INFO("Full logging system initialized successfully");

#ifdef OPTIMIZE_MEMORY
		LOG_INFO("Memory optimization ENABLED: using float precision for posteriors");
#else
		LOG_INFO("Memory optimization DISABLED: using double precision for posteriors");
#endif

		// Set crash context for reference file loading
		crash_set_processing_stage("Loading reference database");
		crash_set_current_file(opt.reference_file);
	}

	// Initialize TSV memory log (independent of verbose logging)
	FILE *tsv_log = NULL;
	if (opt.tsv_log_file[0] != '\0') {
		tsv_log = fopen(opt.tsv_log_file, "w");
		if (tsv_log) {
			// Write header
			fprintf(tsv_log, "# tronko-assign memory log v1.0\n");
			fprintf(tsv_log, "wall_time\tphase\trss_mb\tvm_mb\tpeak_rss_mb\tcpu_user\tcpu_sys\textra_info\n");
			fflush(tsv_log);
			// Set global for readreference.c access
			g_tsv_log_file = tsv_log;
			// Initialize resource monitoring if not already done
			if (!opt.enable_resource_monitoring) {
				init_resource_monitoring();
			}
			TSV_LOG_SIMPLE(tsv_log, "STARTUP");
		} else {
			LOG_WARN("Could not open TSV log file: %s", opt.tsv_log_file);
		}
	}

	// Check if reference file is specified and exists
	if (opt.reference_file[0] == '\0') {
		printf("Error: Reference file not specified. Use -f to specify reference file. Exiting...\n");
		exit(-1);
	}
	
	struct stat st = {0};
	if ( stat(opt.reference_file, &st) == -1 ){
		printf("Cannot find reference_tree.txt file: %s. Exiting...\n", opt.reference_file);
		exit(-1);
	}
	if ( opt.fastq == 1 && opt.number_of_lines_to_read%4 != 0 ){
		printf("You chose FASTQ for your queries but the number of lines to read are not divisible by 4. Change -L to be divisible by 4. Exiting...\n");
	       exit(-1);	
	}
	if ( opt.fastq == 0 && opt.number_of_lines_to_read%2 != 0 ){
		printf("You chose FASTA for your queries but the number of lines to read are not divisible by 2. Change -L to be divisible by 2. Exiting...\n");
		exit(-1);
	}
	int* name_specs = (int*)malloc(3*sizeof(int));
	name_specs[0]=0;
	name_specs[1]=0;
	name_specs[2]=0;

	// Detect reference file format
	int ref_format = detect_reference_format(opt.reference_file);

	if (ref_format == FORMAT_BINARY) {
		// Load uncompressed binary format
		if (opt.verbose_level >= 0) {
			LOG_INFO("Loading binary format reference database: %s", opt.reference_file);
		}
		numberOfTrees = readReferenceBinary(opt.reference_file, name_specs);
		if (numberOfTrees < 0) {
			printf("Error: Failed to load binary reference file: %s. Exiting...\n", opt.reference_file);
			exit(-1);
		}
	} else if (ref_format == FORMAT_BINARY_GZIPPED) {
		// Load gzipped binary format
		if (opt.verbose_level >= 0) {
			LOG_INFO("Loading gzipped binary format reference database: %s", opt.reference_file);
		}
		numberOfTrees = readReferenceBinaryGzipped(opt.reference_file, name_specs);
		if (numberOfTrees < 0) {
			printf("Error: Failed to load gzipped binary reference file: %s. Exiting...\n", opt.reference_file);
			exit(-1);
		}
	} else if (ref_format == FORMAT_BINARY_ZSTD) {
		// Load zstd-compressed binary format
		if (opt.verbose_level >= 0) {
			LOG_INFO("Loading zstd-compressed binary format reference database: %s", opt.reference_file);
		}
		numberOfTrees = readReferenceBinaryZstd(opt.reference_file, name_specs);
		if (numberOfTrees < 0) {
			printf("Error: Failed to load zstd-compressed binary reference file: %s. Exiting...\n", opt.reference_file);
			exit(-1);
		}
	} else if (ref_format == FORMAT_TEXT) {
		// Load text format (existing code path)
		if (opt.verbose_level >= 0) {
			LOG_INFO("Loading text format reference database: %s", opt.reference_file);
		}
		gzFile referenceTree = gzopen(opt.reference_file, "r");
		if (referenceTree == Z_NULL) {
			printf("Error: Cannot open reference file: %s. Exiting...\n", opt.reference_file);
			exit(-1);
		}
		numberOfTrees = readReferenceTree(referenceTree, name_specs);
		gzclose(referenceTree);
	} else {
		printf("Error: Unknown reference file format: %s. Exiting...\n", opt.reference_file);
		exit(-1);
	}

	if (numberOfTrees <= 0) {
		printf("Error: No trees loaded from reference file: %s. Exiting...\n", opt.reference_file);
		exit(-1);
	}
	int max_nodename = name_specs[0];
	int max_taxname = name_specs[1];
	int max_lineTaxonomy = name_specs[2];
	free(name_specs);
	
	if (opt.verbose_level >= 0) {
		char ref_info[256];
		snprintf(ref_info, sizeof(ref_info),
			"Loaded %d trees, max_nodename=%d, max_taxname=%d",
			numberOfTrees, max_nodename, max_taxname);
		log_milestone_with_timing(MILESTONE_REFERENCE_LOADED, ref_info);
	}
	TSV_LOG(tsv_log, "REFERENCE_LOADED", "trees=%d", numberOfTrees);

	if ( opt.print_node_info[0] != '\0' ){
		printf("Printing Accession IDs, Tree numbers, and leaf numbers...\n");
		FILE* tree_info = fopen(opt.print_node_info,"w");
		if (tree_info == NULL ){ printf("Error opening node info file!\n"); exit(1); }
		for(i=0; i<numberOfTrees; i++){
			printTreeInfo(i,rootArr[i],tree_info);
		}
		fclose(tree_info);
		exit(1);
	}
	/*if (opt.print_trees_dir[0] = '\0' ){
		printf("Printing Newick trees used for assignment...\n");
		char newickbuffer[3000];
		for(i=0; i<numberOfTrees; i++){
			FILE *newick_out;
			snprintf(buffer,3000,"%s/%s.nwk",opt.print_trees_dir,i);
			if (NULL==(newick_out=fopen(buffer,"w"))){ puts("Cannot open newick file!\n"); exit(-1); }
			printNewick(newick_out,i,rootArr[i]);
			fclose(newick_out);
		}
		exit(1);
	}*/
	opt.print_leave_seqs=0;
	if (opt.print_leave_seqs == 1){
		printf("Printingleavesfile\n");
		FILE* leaves_file = fopen("leaves.fasta","w");
		if(leaves_file == NULL ){ printf("Error opening leaves file!\n"); exit(1); }
		for(i=0; i<numberOfTrees; i++){
			printLeaveSeqsToFile(leaves_file,rootArr[i],i,numbaseArr[i]);	
		}
		exit(1);
	}
	store_PPs_Arr(numberOfTrees,opt.score_constant);
	int maxNumSpec=0;
	int maxNumBase=0;
	int *specs = (int*)malloc(2*sizeof(int));
	specs[0]=0;
	specs[1]=0;
	int numspec_total=setNumbase_setNumspec(numberOfTrees,specs);
	maxNumSpec = specs[0];
	maxNumBase = specs[1];
	free(specs);
	Cinterval = opt.cinterval;
	//HASHMAP(char, leafMap) map;
	//hashmap_init(&map, hashmap_hash_string, strcmp);
	//for(i=0; i<numberOfTrees; i++){
	//	for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
	//		struct leafMap *l;
	//		l = malloc(sizeof(*l));
	//		l->name = treeArr[i][j].name;
	//		l->root = i;
	//		l->node = j;
	//		hashmap_put(&map,l->name,l);
	//	}
	//}
	//leaf_map = (leafMap *)malloc(numspec_total*sizeof(leafMap));
	//k=0;
	//for(i=0; i<numberOfTrees; i++){
	//	for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
	//		leaf_map[k].name = (char*)malloc(max_nodename*sizeof(char));
	//		strcpy(leaf_map[k].name,treeArr[i][j].name);
	//		leaf_map[k].root = i;
	//		leaf_map[k].node = j;
	//		k++;
	//	}
	//}
	int *read_specs = (int*)malloc(2*sizeof(int));
	read_specs[0] = 0;
	read_specs[1] = 0;
	int number_of_total_nodes = 0;
	for(i=0; i<numberOfTrees; i++){
		number_of_total_nodes += 2*numspecArr[i]-1;
	}
	int max_name_length = 0;
	int max_query_length = 0;
	int numberOfLinesToRead=opt.number_of_lines_to_read;
	mystruct mstr[opt.number_of_cores];//array of stuct that contains input and output for each thread
	if ( strcmp("single",opt.paired_or_single)==0){
		if (opt.skip_build==0){
			if (opt.verbose_level >= 0) {
				LOG_INFO("Building BWA index for: %s", opt.fasta_file);
			}
			bwa_index(2,opt.fasta_file);
			if (opt.verbose_level >= 0) {
				LOG_MILESTONE_TIMED(MILESTONE_BWA_INDEX_BUILT);
			}
			TSV_LOG_SIMPLE(tsv_log, "BWA_INDEX");
		} else {
			if (opt.verbose_level >= 0) {
				LOG_INFO("Skipping BWA index build");
			}
		}
		CompressedFile* reads_file = cf_open(opt.read1_file, "r");
		if ( reads_file == NULL ){
			printf("**reads file could not be opened.\n");
			exit(-1);
		}
		find_specs_for_reads_cf(read_specs, reads_file, opt.fastq);
		cf_close(reads_file);
		max_name_length = read_specs[0];
		max_query_length = read_specs[1];
		free(read_specs);
		
		if (opt.verbose_level >= 0) {
			char specs_info[256];
			snprintf(specs_info, sizeof(specs_info), 
				"Read specs detected: max_name=%d, max_query=%d", 
				max_name_length, max_query_length);
			log_milestone_with_timing(MILESTONE_READ_SPECS_DETECTED, specs_info);
		}
		singleQueryMat = malloc(sizeof(struct queryMatSingle));
		if ( opt.fastq == 0 ){
			singleQueryMat->queryMat = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/2));
			singleQueryMat->name = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/2));
			for (i=0; i<numberOfLinesToRead/2; i++){
				singleQueryMat->queryMat[i] = (char *)malloc(sizeof(char)*(max_query_length+1));
				singleQueryMat->name[i] = (char *)malloc(sizeof(char)*(max_name_length+1));
			}
		}else{
			singleQueryMat->queryMat = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/4));
			singleQueryMat->name = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/4));
			for (i=0; i<numberOfLinesToRead/4; i++){
				singleQueryMat->queryMat[i] = (char *)malloc(sizeof(char)*(max_query_length+1));
				singleQueryMat->name[i] = (char *)malloc(sizeof(char)*(max_name_length+1));
			}
		}
		FILE *results = NULL;
#ifdef ENABLE_PARQUET
		parquet_writer_t *parquet_writer = NULL;
		if (opt.parquet_enabled) {
			char parquet_filename[MAXFILENAME];
			snprintf(parquet_filename, MAXFILENAME, "%s.parquet", opt.parquet_prefix);
			char parquet_err[256];
			parquet_writer = parquet_writer_create(parquet_filename, parquet_err);
			if (!parquet_writer) {
				printf("Error creating Parquet writer: %s\n", parquet_err);
				exit(1);
			}
			printf("Writing results to Parquet: %s\n", parquet_filename);
		} else {
#endif
			results = fopen(opt.results_file,"w");
			if ( results == NULL ){ printf("Error opening output file!\n"); exit(1); }
			fprintf(results,"Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");
#ifdef ENABLE_PARQUET
		}
#endif
		int keepTrackOfReadLine=0;
		pthread_t threads[opt.number_of_cores];//array of our threads
		int divideFile, start, end;
		int returnLineNumber=0;
		CompressedFile* seqinfile = cf_open(opt.read1_file, "r");
		if (seqinfile == NULL){
			printf("*** fasta file could not be opened.\n");
			exit(-1);
		}
		int first_iter=1;
		for(i=0; i<opt.number_of_cores; i++){
			mstr[i].str = malloc(sizeof(struct resultsStruct));
			if ( opt.fastq == 0){
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,0,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}else{
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,0,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}
			mstr[i].concordant=0;
			mstr[i].maxNumSpec=maxNumSpec;
			mstr[i].databasefile=opt.fasta_file;
			mstr[i].numspec_total=numspec_total;
			mstr[i].use_nw = opt.use_nw;
			mstr[i].print_alignments_to_file = opt.print_alignments_to_file;
			mstr[i].print_unassigned = opt.unassigned;
			mstr[i].use_leaf_portion = opt.use_leaf_portion;
			mstr[i].padding = opt.padding;
			mstr[i].max_query_length = max_query_length;
			mstr[i].max_readname_length = max_name_length;
			mstr[i].max_acc_name = max_nodename;
			mstr[i].max_numbase = maxNumBase;
			mstr[i].max_lineTaxonomy = max_lineTaxonomy;
			mstr[i].number_of_total_nodes = number_of_total_nodes;
			mstr[i].print_all_nodes = opt.print_all_nodes;
			// Tier 1 optimization settings
			mstr[i].early_termination = opt.early_termination;
			mstr[i].strike_box = opt.strike_box;
			mstr[i].max_strikes = opt.max_strikes;
			mstr[i].enable_pruning = opt.enable_pruning;
			mstr[i].pruning_factor = opt.pruning_factor;
		}

		if (opt.verbose_level >= 0) {
			char mem_info[256];
			snprintf(mem_info, sizeof(mem_info),
				"Memory allocated: threads=%d, lines_per_batch=%d, total_nodes=%d",
				opt.number_of_cores, numberOfLinesToRead, number_of_total_nodes);
			log_milestone_with_timing(MILESTONE_MEMORY_ALLOCATED, mem_info);

			char thread_info[256];
			snprintf(thread_info, sizeof(thread_info),
				"Thread structures initialized for %d cores", opt.number_of_cores);
			log_milestone_with_timing(MILESTONE_THREADS_INITIALIZED, thread_info);
		}
		TSV_LOG(tsv_log, "THREADS_ALLOCATED", "threads=%d", opt.number_of_cores);

		int batch_count = 0;
		while (1){
			batch_count++;
			if (opt.verbose_level >= 0) {
				char start_info[128];
				snprintf(start_info, sizeof(start_info), "Starting batch %d", batch_count);
				log_milestone_with_timing(MILESTONE_BATCH_START, start_info);
			}
			TSV_LOG(tsv_log, "BATCH_START", "batch=%d", batch_count);

			if (opt.fastq==0){
				returnLineNumber=readInXNumberOfLines_cf(numberOfLinesToRead/2,seqinfile,0,opt,max_query_length,max_name_length);
			}else{
				returnLineNumber=readInXNumberOfLines_fastq_cf(numberOfLinesToRead/4,seqinfile,0,opt,max_query_length,max_name_length,first_iter);
			}
			if (returnLineNumber==0){
				break;
			}
			
			if (opt.verbose_level >= 0) {
				char batch_info[256];
				snprintf(batch_info, sizeof(batch_info),
					"Batch %d loaded: %d reads", batch_count, returnLineNumber);
				log_milestone_with_timing(MILESTONE_BATCH_LOADED, batch_info);
			}
			TSV_LOG(tsv_log, "BATCH_LOADED", "batch=%d,reads=%d", batch_count, returnLineNumber);
			divideFile = returnLineNumber/opt.number_of_cores;
			first_iter=0;
			j=0;
			for (i=0; i<opt.number_of_cores; i++){
				start=j;
				end=j+divideFile;
				if ( i==opt.number_of_cores-1){
					end=returnLineNumber;
				}
				mstr[i].start=start;
				mstr[i].end=end;
				mstr[i].paired=0;
				mstr[i].ntree = numberOfTrees;
				mstr[i].alignmentsdir = opt.print_alignments_dir;
				j=j+divideFile;
				mstr[i].str->taxonPath =(char**) malloc((end-start)*(sizeof(char *)));
				for(k=0; k<end-start; k++){
					mstr[i].str->taxonPath[k] = malloc((max_name_length+max_lineTaxonomy+120)*(sizeof(char)));
				}
			}
			
			if (opt.verbose_level >= 0) {
				LOG_DEBUG("Creating %d threads for batch processing", opt.number_of_cores);
			}
			
			for(i=0; i<opt.number_of_cores;i++){
				pthread_create(&threads[i], NULL, runAssignmentOnChunk_WithBWA, &mstr[i]);
			}
			for ( i=0; i<opt.number_of_cores;i++){
				pthread_join(threads[i], NULL);
			}
			
			if (opt.verbose_level >= 0) {
				LOG_MILESTONE_TIMED(MILESTONE_PLACEMENT_COMPLETE);
			}
			
#ifdef ENABLE_PARQUET
			if (opt.parquet_enabled) {
				/* Collect all results into a batch for Parquet */
				int total_results = 0;
				for (i = 0; i < opt.number_of_cores; i++) {
					total_results += (mstr[i].end - mstr[i].start);
				}
				if (total_results > 0) {
					assignmentResult *batch = calloc(total_results, sizeof(assignmentResult));
					int idx = 0;
					for (i = 0; i < opt.number_of_cores; i++) {
						for (j = 0; j < (mstr[i].end - mstr[i].start); j++) {
							parse_tsv_to_result(mstr[i].str->taxonPath[j], &batch[idx]);
							idx++;
						}
					}
					parquet_writer_write_batch(parquet_writer, batch, total_results);
					for (idx = 0; idx < total_results; idx++) {
						free_assignment_result(&batch[idx]);
					}
					free(batch);
				}
			} else {
#endif
				for ( i=0; i<opt.number_of_cores; i++){
					for ( j=0; j<(mstr[i].end-mstr[i].start); j++){
						fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
					}
				}
#ifdef ENABLE_PARQUET
			}
#endif

			if (opt.verbose_level >= 0) {
				char results_info[256];
				snprintf(results_info, sizeof(results_info),
					"Batch %d results written", batch_count);
				log_milestone_with_timing(MILESTONE_RESULTS_WRITTEN, results_info);

				char complete_info[256];
				snprintf(complete_info, sizeof(complete_info),
					"Batch %d completed: %d reads processed", batch_count, returnLineNumber);
				log_milestone_with_timing(MILESTONE_BATCH_COMPLETE, complete_info);
			}
			TSV_LOG(tsv_log, "BATCH_COMPLETE", "batch=%d", batch_count);
			for ( i=0; i<opt.number_of_cores; i++){
				for(j=0; j<mstr[i].end-mstr[i].start; j++){
					free(mstr[i].str->taxonPath[j]);
				}
				free(mstr[i].str->taxonPath);
			}
		}
		
		if (opt.verbose_level >= 0) {
			char cleanup_info[256];
			snprintf(cleanup_info, sizeof(cleanup_info),
				"Processing completed. Processed %d batches", batch_count);
			log_milestone_with_timing(MILESTONE_CLEANUP_START, cleanup_info);
		}

		// Log BWA overflow statistics summary
		if (g_overflow_read_count > 0) {
			LOG_WARN("=== BWA BOUNDS CAP SUMMARY ===");
			LOG_WARN("  Reads hitting cap: %d", g_overflow_read_count);
			LOG_WARN("  Total matches dropped: %d", g_total_dropped_matches);
			LOG_WARN("  Max potential matches seen: %d (cap is %d)",
			         g_max_potential_matches, MAX_NUM_BWA_MATCHES);
			LOG_WARN("  Consider increasing MAX_NUM_BWA_MATCHES if accuracy is affected");
		}

		// Close files
#ifdef ENABLE_PARQUET
		if (opt.parquet_enabled) {
			if (parquet_writer_close(parquet_writer) != 0) {
				printf("Warning: Error closing Parquet writer\n");
			}
		} else {
#endif
			fclose(results);
#ifdef ENABLE_PARQUET
		}
#endif
		cf_close(seqinfile);

		if (opt.verbose_level >= 0) {
			log_current_resource_usage("After closing files");
		}
		if ( opt.fastq == 0 ){
			for(i=0; i<numberOfLinesToRead/2; i++){
				free(singleQueryMat->queryMat[i]);
				free(singleQueryMat->name[i]);
			}
		}else{
			for(i=0; i<numberOfLinesToRead/4; i++){
				free(singleQueryMat->queryMat[i]);
				free(singleQueryMat->name[i]);
			}
		}
		free(singleQueryMat->queryMat);
		free(singleQueryMat->name);
		free(singleQueryMat);
		
		if (opt.verbose_level >= 0) {
			log_current_resource_usage("After freeing query matrices");
		}
		
		if (opt.verbose_level >= 0) {
			LOG_MILESTONE_TIMED(MILESTONE_CLEANUP_COMPLETE);
			LOG_MILESTONE_TIMED(MILESTONE_PROGRAM_END);
			
			// Cleanup crash debugging systems
			crash_debug_cleanup();
			symbol_resolver_cleanup();
			
			logger_cleanup();
			cleanup_resource_monitoring();
		}
	}else{
		CompressedFile* seqinfile_1 = cf_open(opt.read1_file, "r");
		CompressedFile* seqinfile_2 = cf_open(opt.read2_file, "r");
		if (seqinfile_1 == NULL){
			printf("*** fasta/fastq file could not be opened.\n");
			exit(-1);
		}
		if (seqinfile_2 == NULL){
			printf("*** fasta/fastq file could not be opened.\n");
			exit(-1);
		}
		find_specs_for_reads_cf(read_specs, seqinfile_1, opt.fastq);
		cf_close(seqinfile_1);
		find_specs_for_reads_cf(read_specs, seqinfile_2, opt.fastq);
		cf_close(seqinfile_2);
		max_name_length = read_specs[0];
		max_query_length = read_specs[1];
		free(read_specs);
		if (opt.skip_build==0){
			bwa_index(2,opt.fasta_file);
		}
		pairedQueryMat = malloc(sizeof(struct queryMatPaired));
		if (opt.fastq==0){
			pairedQueryMat->query1Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->query2Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->forward_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->reverse_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/2);
			for ( i=0; i<numberOfLinesToRead/2; i++){
				pairedQueryMat->query1Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->query2Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->forward_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				pairedQueryMat->reverse_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				for(j=0; j<max_name_length+1; j++){
					pairedQueryMat->forward_name[i][j]='\0';
					pairedQueryMat->reverse_name[i][j]='\0';
				}
				for(j=0; j<max_query_length+1; j++){
					pairedQueryMat->query1Mat[i][j] = '\0';
					pairedQueryMat->query2Mat[i][j] = '\0';
				}
			}
		}else{
			pairedQueryMat->query1Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->query2Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->forward_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->reverse_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/4);
			for ( i=0; i<numberOfLinesToRead/4; i++){
				pairedQueryMat->query1Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->query2Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->forward_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				 pairedQueryMat->reverse_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				for(j=0; j<max_name_length+1; j++){
					pairedQueryMat->forward_name[i][j]='\0';
					pairedQueryMat->reverse_name[i][j]='\0';
				}
				for(j=0; j<max_query_length+1; j++){
					pairedQueryMat->query1Mat[i][j]='\0';
					pairedQueryMat->query2Mat[i][j]='\0';
				}
			}
		}
		FILE *results = NULL;
#ifdef ENABLE_PARQUET
		parquet_writer_t *parquet_writer = NULL;
		if (opt.parquet_enabled) {
			char parquet_filename[MAXFILENAME];
			snprintf(parquet_filename, MAXFILENAME, "%s.parquet", opt.parquet_prefix);
			char parquet_err[256];
			parquet_writer = parquet_writer_create(parquet_filename, parquet_err);
			if (!parquet_writer) {
				printf("Error creating Parquet writer: %s\n", parquet_err);
				exit(1);
			}
			printf("Writing results to Parquet: %s\n", parquet_filename);
		} else {
#endif
			results = fopen(opt.results_file,"w");
			if ( results == NULL ){ printf("Error opening output file!\n"); exit(1); }
			fprintf(results,"Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");
#ifdef ENABLE_PARQUET
		}
#endif
		int keepTrackOfReadLine=0;
		pthread_t threads[opt.number_of_cores];//array of our thrads
		int divideFile, start, end;
		int returnLineNumber=0;//<- this is the number of records that we have read.
		int returnLineNumber2 =0;
		int concordant=1;
		seqinfile_1 = cf_open(opt.read1_file, "r");
		seqinfile_2 = cf_open(opt.read2_file, "r");
		if (seqinfile_1 == NULL){
			printf("*** fasta/fastq file could not be opened.\n");
			exit(-1);
		}
		if (seqinfile_2 == NULL){
			printf("*** fasta/fastq file could not be opened.\n");
			exit(-1);
		}
		int first_iter=1;
		for(i=0;i<opt.number_of_cores;i++){
			mstr[i].str = malloc(sizeof(struct resultsStruct));
			if ( opt.fastq == 0){
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,1,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}else{
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,1,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}
			mstr[i].concordant=concordant;
			mstr[i].maxNumSpec=maxNumSpec;
			mstr[i].databasefile=opt.fasta_file;
			mstr[i].numspec_total=numspec_total;
			mstr[i].use_nw = opt.use_nw;
			mstr[i].print_alignments_to_file = opt.print_alignments_to_file;
			mstr[i].print_unassigned = opt.unassigned;
			mstr[i].use_leaf_portion = opt.use_leaf_portion;
			mstr[i].padding = opt.padding;
			mstr[i].max_query_length = max_query_length;
			mstr[i].max_readname_length = max_name_length;
			mstr[i].max_acc_name = max_nodename;
			mstr[i].max_numbase = maxNumBase;
			mstr[i].max_lineTaxonomy = max_lineTaxonomy;
			mstr[i].number_of_total_nodes = number_of_total_nodes;
			mstr[i].print_all_nodes = opt.print_all_nodes;
			// Tier 1 optimization settings
			mstr[i].early_termination = opt.early_termination;
			mstr[i].strike_box = opt.strike_box;
			mstr[i].max_strikes = opt.max_strikes;
			mstr[i].enable_pruning = opt.enable_pruning;
			mstr[i].pruning_factor = opt.pruning_factor;
		}
		while (1){
			crash_set_processing_stage("Reading paired-end input files");
			crash_set_current_file(opt.read1_file);
			if (opt.fastq==0){
				returnLineNumber=readInXNumberOfLines_cf(numberOfLinesToRead/2,seqinfile_1,1,opt,max_query_length,max_name_length);
			}else{
				returnLineNumber=readInXNumberOfLines_fastq_cf(numberOfLinesToRead/4,seqinfile_1,1,opt,max_query_length,max_name_length,first_iter);
			}
			if (returnLineNumber==0)
				break;
			crash_set_current_file(opt.read2_file);
			if (opt.fastq==0){
				returnLineNumber2 = readInXNumberOfLines_cf(numberOfLinesToRead/2, seqinfile_2, 2, opt,max_query_length,max_name_length);
			}else{
				returnLineNumber2 = readInXNumberOfLines_fastq_cf(numberOfLinesToRead/4,seqinfile_2, 2, opt,max_query_length,max_name_length,first_iter);
			}
			returnLineNumber = returnLineNumber2;
			first_iter=0;
			divideFile= returnLineNumber/opt.number_of_cores;
			j=0;
			for ( i=0; i<opt.number_of_cores; i++){
				start=j;
				end=j+divideFile;
				if(i==opt.number_of_cores-1)
					end=returnLineNumber;
				mstr[i].start=start;
				mstr[i].end=end;
				mstr[i].paired = 1;
				mstr[i].ntree = numberOfTrees;
				mstr[i].alignmentsdir = opt.print_alignments_dir;
				j=j+divideFile;
				mstr[i].str->taxonPath =(char**) malloc((end-start)*(sizeof(char *)));
				for(k=0; k<end-start; k++){
					mstr[i].str->taxonPath[k] = malloc((max_name_length+max_lineTaxonomy+120)*(sizeof(char)));
				}
			}
			for (i=0; i<opt.number_of_cores;i++){
				pthread_create(&threads[i], NULL, runAssignmentOnChunk_WithBWA, &mstr[i]);
			}
			for ( i=0; i<opt.number_of_cores;i++){
				pthread_join(threads[i], NULL);
			}
#ifdef ENABLE_PARQUET
			if (opt.parquet_enabled) {
				/* Collect all results into a batch for Parquet */
				int total_results = 0;
				for (i = 0; i < opt.number_of_cores; i++) {
					total_results += (mstr[i].end - mstr[i].start);
				}
				if (total_results > 0) {
					assignmentResult *batch = calloc(total_results, sizeof(assignmentResult));
					int idx = 0;
					for (i = 0; i < opt.number_of_cores; i++) {
						for (j = 0; j < (mstr[i].end - mstr[i].start); j++) {
							parse_tsv_to_result(mstr[i].str->taxonPath[j], &batch[idx]);
							idx++;
						}
					}
					parquet_writer_write_batch(parquet_writer, batch, total_results);
					for (idx = 0; idx < total_results; idx++) {
						free_assignment_result(&batch[idx]);
					}
					free(batch);
				}
			} else {
#endif
				for ( i=0; i<opt.number_of_cores; i++){
					for ( j=0; j<(mstr[i].end-mstr[i].start); j++){
						fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
					}
				}
#ifdef ENABLE_PARQUET
			}
#endif
			for ( i=0; i<opt.number_of_cores; i++){
				for(j=0; j<mstr[i].end-mstr[i].start; j++){
					free(mstr[i].str->taxonPath[j]);
				}
				free(mstr[i].str->taxonPath);
			}
		}

		// Log BWA overflow statistics summary for paired-end
		if (g_overflow_read_count > 0) {
			LOG_WARN("=== BWA BOUNDS CAP SUMMARY ===");
			LOG_WARN("  Reads hitting cap: %d", g_overflow_read_count);
			LOG_WARN("  Total matches dropped: %d", g_total_dropped_matches);
			LOG_WARN("  Max potential matches seen: %d (cap is %d)",
			         g_max_potential_matches, MAX_NUM_BWA_MATCHES);
			LOG_WARN("  Consider increasing MAX_NUM_BWA_MATCHES if accuracy is affected");
		}

#ifdef ENABLE_PARQUET
		if (opt.parquet_enabled) {
			if (parquet_writer_close(parquet_writer) != 0) {
				printf("Warning: Error closing Parquet writer\n");
			}
		} else {
#endif
			fclose(results);
#ifdef ENABLE_PARQUET
		}
#endif
		cf_close(seqinfile_1);
		cf_close(seqinfile_2);
		if (opt.fastq==0){
			for(i=0; i<numberOfLinesToRead/2; i++){
				free(pairedQueryMat->query1Mat[i]);
				free(pairedQueryMat->query2Mat[i]);
				free(pairedQueryMat->forward_name[i]);
				free(pairedQueryMat->reverse_name[i]);
			}
		}else{
			for(i=0; i<numberOfLinesToRead/4; i++){
				free(pairedQueryMat->query1Mat[i]);
				free(pairedQueryMat->query2Mat[i]);
				free(pairedQueryMat->forward_name[i]);
				free(pairedQueryMat->reverse_name[i]);
			}
		}
		free(pairedQueryMat->query1Mat);
		free(pairedQueryMat->query2Mat);
		free(pairedQueryMat->forward_name);
		free(pairedQueryMat->reverse_name);
		free(pairedQueryMat);
	}
	// Free thread result structures
	if (opt.verbose_level >= 0) {
		log_current_resource_usage("Before freeing thread structures");
	}
	
	for(i=0; i<opt.number_of_cores; i++){
		if (opt.fastq == 0){
			if (strcmp("single",opt.paired_or_single)==0){
				freeMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, 0, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}else{
				freeMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, 1, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}
		}else{
			if (strcmp("single",opt.paired_or_single)==0){
				freeMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, 0, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}else{
				freeMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, 1, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}
		}
	}
	
	if (opt.verbose_level >= 0) {
		log_current_resource_usage("After freeing thread structures");
	}
	
	// Free tree data structures - this is where major memory should be freed
	if (opt.verbose_level >= 0) {
		char tree_info[256];
		snprintf(tree_info, sizeof(tree_info), "Before freeing tree data (%d trees)", numberOfTrees);
		log_current_resource_usage(tree_info);
	}
	
	for(i=0; i<numberOfTrees; i++){
		for(j=0; j<2*numspecArr[i]-1; j++){
			free(treeArr[i][j].posteriornc);  /* Single free per node (1D array) */
		}
		for(j=numspecArr[i]-1; j<(2*numspecArr[i]-1); j++){
			free(treeArr[i][j].name);
		}
		free(treeArr[i]);
	}
	free(treeArr);
	
	if (opt.verbose_level >= 0) {
		log_current_resource_usage("After freeing tree arrays");
	}
	
	// Free taxonomyArr - 4-dimensional array [trees][species][taxonomy_levels][taxonomy_names]
	if (taxonomyArr) {
		for(i=0; i<numberOfTrees; i++){
			if (taxonomyArr[i]) {
				for(j=0; j<numspecArr[i]; j++){
					if (taxonomyArr[i][j]) {
						for(k=0; k<7; k++){  // 7 phylogeny levels
							if (taxonomyArr[i][j][k]) {
								free(taxonomyArr[i][j][k]);
							}
						}
						free(taxonomyArr[i][j]);
					}
				}
				free(taxonomyArr[i]);
			}
		}
		free(taxonomyArr);
	}
	
	if (opt.verbose_level >= 0) {
		log_current_resource_usage("After freeing taxonomy arrays");
	}
	
	free(numbaseArr);
	free(rootArr);
	free(numspecArr);

	if (opt.verbose_level >= 0) {
		log_current_resource_usage("After freeing all tree data");
	}

	// Close TSV log file
	TSV_LOG_SIMPLE(tsv_log, "FINAL");
	if (tsv_log) {
		fclose(tsv_log);
		g_tsv_log_file = NULL;
	}
}
