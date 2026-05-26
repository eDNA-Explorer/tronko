#ifndef _ALLOCFORPTHREADS_
#define _ALLOCFORPTHREADS_

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "needleman_wunsch.h"
#include "alignment.h"
#include "alignment_scoring.h"

void allocateMemForResults( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int print_alignments, int maxNumSpec, int paired, int use_nw, int max_lineTaxonomy, int max_name_length, int max_query_length, int max_numbase, int use_portion, int padding_size, int number_of_total_nodes, int max_leaf_matches, int normalize_scores);
void freeMemForResults ( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int paired, int use_nw, int use_portion, int maxNumSpec, int number_of_total_nodes, int max_leaf_matches);
#endif /* _ALLOCFORPTHREADS_ */
