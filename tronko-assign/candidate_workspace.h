#ifndef TRONKO_CANDIDATE_WORKSPACE_H
#define TRONKO_CANDIDATE_WORKSPACE_H

#include <stdbool.h>
#include <stddef.h>

#include "global.h"

typedef struct CandidatePlacement {
	int tree_id;
	int leaf_node;
	size_t score_offset;
	size_t node_count;
	int forward_start;
	int reverse_start;
	bool concordant;
	char forward_cigar[MAX_CIGAR];
	char reverse_cigar[MAX_CIGAR];
	int winning_count;
	int lca_node;
	bool has_winner;
} CandidatePlacement;

typedef struct CandidateWorkspace {
	CandidatePlacement candidates[MAX_NUM_BWA_MATCHES];
	size_t candidate_count;
	type_of_PP *score_arena;
	size_t score_capacity;
	size_t score_length;
} CandidateWorkspace;

typedef int (*CandidateLcaFunction)(int node1, int node2, int tree_id);

void candidate_workspace_init(CandidateWorkspace *workspace);
void candidate_workspace_reset(CandidateWorkspace *workspace);
void candidate_workspace_destroy(CandidateWorkspace *workspace);

int candidate_workspace_add(CandidateWorkspace *workspace, int tree_id,
	int leaf_node, int forward_start, const char *forward_cigar,
	int reverse_start, const char *reverse_cigar);

int candidate_workspace_select(CandidateWorkspace *workspace,
	const bwaMatches *matches, int concordant, int paired, int use_portion,
	size_t *dropped_count);

int candidate_workspace_prepare_scores(CandidateWorkspace *workspace,
	const int *species_counts, size_t number_of_trees);

type_of_PP *candidate_workspace_scores(CandidateWorkspace *workspace,
	size_t candidate_index);
const type_of_PP *candidate_workspace_scores_const(
	const CandidateWorkspace *workspace, size_t candidate_index);
int candidate_workspace_summarize_winners(CandidateWorkspace *workspace,
	type_of_PP maximum, type_of_PP interval, CandidateLcaFunction get_lca);
size_t candidate_workspace_collect_positive(const CandidateWorkspace *workspace,
	size_t indices[MAX_NUM_BWA_MATCHES]);

#ifdef TRONKO_TESTING
void candidate_workspace_fail_next_growth(void);
#endif

#endif
