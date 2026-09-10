#include "allocateMemoryForResults.h"
#include "candidate_workspace.h"

#include <stdint.h>
#include <string.h>

int allocateMemForResults( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int print_alignments, int maxNumSpec, int paired, int use_nw, int max_lineTaxonomy, int max_name_length, int max_query_length, int max_numbase, int use_portion, int padding_size, int number_of_total_nodes){
	size_t scratch_length;
	(void)sizeOfChunk;
	(void)num_threads;
	(void)numberOfTrees;
	(void)maxNumSpec;
	(void)paired;
	(void)max_lineTaxonomy;
	(void)max_name_length;
	(void)number_of_total_nodes;
	if (results == NULL || max_query_length < 0 || max_numbase < 0 ||
		padding_size < 0) return -1;
	memset(results, 0, sizeof(*results));
	if (use_portion==1){
		scratch_length = (size_t)max_query_length * 2 +
			(size_t)padding_size * 2 + 1;
	}else{
		scratch_length = (size_t)max_query_length + (size_t)max_numbase + 1;
	}
	if (scratch_length > SIZE_MAX / sizeof(int)) return -1;
	results->positions = malloc(scratch_length * sizeof(*results->positions));
	results->locQuery = malloc(scratch_length * sizeof(*results->locQuery));
	results->workspace = malloc(sizeof(*results->workspace));
	if (results->positions == NULL || results->locQuery == NULL ||
		results->workspace == NULL) goto allocation_failed;
	candidate_workspace_init(results->workspace);
	candidate_workspace_reset(results->workspace);
	if ( use_nw == 1 ){
		results->nw=needleman_wunsch_new();
		results->aln=alignment_create(max_query_length+max_numbase+1);
		results->scoring=malloc(sizeof(scoring_t));
		if (results->nw == NULL || results->aln == NULL ||
			results->scoring == NULL) goto allocation_failed;
	}
	results->minimum=(type_of_PP *)malloc(3*sizeof(type_of_PP));
	if (results->minimum == NULL) goto allocation_failed;
	if ( print_alignments == 1){
		results->print_alignments=1;
	}else{
		results->print_alignments=0;
	}
	if (use_nw==1){
		int match=2;
		int mismatch=-1;
		int gap_open=-3;
		int gap_extend=-1;
		int num_of_mismatches=0;
		int num_of_indels = 0;
		bool no_start_gap_penalty=true;
		bool no_end_gap_penalty=true;
		bool no_gaps_in_a = false;
		bool no_gaps_in_b = false;
		bool no_mismatches = false;
		bool case_sensitive=false;
		scoring_init(results->scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	}
	return 0;

allocation_failed:
	free(results->positions);
	free(results->locQuery);
	if (results->workspace != NULL) {
		candidate_workspace_destroy(results->workspace);
		free(results->workspace);
	}
	if (results->nw != NULL) needleman_wunsch_free(results->nw);
	if (results->aln != NULL) alignment_free(results->aln);
	free(results->scoring);
	free(results->minimum);
	memset(results, 0, sizeof(*results));
	return -1;
}
void freeMemForResults ( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int paired, int use_nw, int use_portion, int maxNumSpec, int number_of_total_nodes){
	(void)sizeOfChunk;
	(void)num_threads;
	(void)numberOfTrees;
	(void)paired;
	(void)use_portion;
	(void)maxNumSpec;
	(void)number_of_total_nodes;
	if (results == NULL) return;
	free(results->positions);
	free(results->locQuery);
	if (results->workspace != NULL)
		candidate_workspace_destroy(results->workspace);
	free(results->workspace);
	free(results->minimum);
	if (use_nw==1){
		if (results->nw != NULL) needleman_wunsch_free(results->nw);
		if (results->aln != NULL) alignment_free(results->aln);
		free(results->scoring);
	}
	free(results);
}
