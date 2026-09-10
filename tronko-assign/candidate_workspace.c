#include "candidate_workspace.h"

#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef TRONKO_TESTING
static int fail_next_growth;

void candidate_workspace_fail_next_growth(void)
{
	fail_next_growth = 1;
}
#endif

void candidate_workspace_init(CandidateWorkspace *workspace)
{
	if (workspace != NULL) memset(workspace, 0, sizeof(*workspace));
}

void candidate_workspace_reset(CandidateWorkspace *workspace)
{
	size_t i;
	if (workspace == NULL) return;
	workspace->candidate_count = 0;
	workspace->score_length = 0;
	for (i = 0; i < MAX_NUM_BWA_MATCHES; ++i) {
		CandidatePlacement *candidate = &workspace->candidates[i];
		candidate->tree_id = -1;
		candidate->leaf_node = -1;
		candidate->score_offset = 0;
		candidate->node_count = 0;
		candidate->forward_start = -1;
		candidate->reverse_start = -1;
		candidate->concordant = false;
		candidate->forward_cigar[0] = '\0';
		candidate->reverse_cigar[0] = '\0';
		candidate->winning_count = 0;
		candidate->lca_node = -1;
		candidate->has_winner = false;
	}
}

void candidate_workspace_destroy(CandidateWorkspace *workspace)
{
	if (workspace == NULL) return;
	free(workspace->score_arena);
	memset(workspace, 0, sizeof(*workspace));
}

static void copy_cigar(char destination[MAX_CIGAR], const char *source)
{
	if (source == NULL) {
		destination[0] = '\0';
		return;
	}
	strncpy(destination, source, MAX_CIGAR - 1);
	destination[MAX_CIGAR - 1] = '\0';
}

int candidate_workspace_add(CandidateWorkspace *workspace, int tree_id,
	int leaf_node, int forward_start, const char *forward_cigar,
	int reverse_start, const char *reverse_cigar)
{
	size_t i;
	CandidatePlacement *candidate;
	if (workspace == NULL || tree_id < 0 || leaf_node < 0) return -EINVAL;
	for (i = 0; i < workspace->candidate_count; ++i) {
		if (workspace->candidates[i].tree_id == tree_id) return 1;
	}
	if (workspace->candidate_count >= MAX_NUM_BWA_MATCHES) return 2;
	candidate = &workspace->candidates[workspace->candidate_count++];
	candidate->tree_id = tree_id;
	candidate->leaf_node = leaf_node;
	candidate->forward_start = forward_start;
	candidate->reverse_start = reverse_start;
	copy_cigar(candidate->forward_cigar, forward_cigar);
	copy_cigar(candidate->reverse_cigar, reverse_cigar);
	return 0;
}

static int add_match_list(CandidateWorkspace *workspace, const int *roots,
	const int *nodes, const int *forward_starts, char *const *forward_cigars,
	const int *reverse_starts, char *const *reverse_cigars, int paired,
	int use_portion, bool concordant, size_t *dropped_count)
{
	size_t i;
	if (roots == NULL || nodes == NULL) return -EINVAL;
	for (i = 0; i < MAX_NUM_BWA_MATCHES && roots[i] >= 0; ++i) {
		int status = candidate_workspace_add(workspace, roots[i], nodes[i],
			use_portion && forward_starts != NULL ? forward_starts[i] : -1,
			use_portion && forward_cigars != NULL ? forward_cigars[i] : NULL,
			use_portion && paired && reverse_starts != NULL ? reverse_starts[i] : -1,
			use_portion && paired && reverse_cigars != NULL ? reverse_cigars[i] : NULL);
		if (status < 0) return status;
		if (status == 0)
			workspace->candidates[workspace->candidate_count - 1].concordant =
				concordant;
		if (status == 2 && dropped_count != NULL) (*dropped_count)++;
	}
	return 0;
}

int candidate_workspace_select(CandidateWorkspace *workspace,
	const bwaMatches *matches, int concordant, int paired, int use_portion,
	size_t *dropped_count)
{
	int status;
	if (workspace == NULL || matches == NULL) return -EINVAL;
	candidate_workspace_reset(workspace);
	if (dropped_count != NULL) *dropped_count = matches->dropped_matches;

	if (concordant && matches->concordant_matches_roots[0] >= 0) {
		return add_match_list(workspace, matches->concordant_matches_roots,
			matches->concordant_matches_nodes, matches->starts_forward,
			matches->cigars_forward, matches->starts_reverse,
			matches->cigars_reverse, paired, use_portion, true, dropped_count);
	}
	status = add_match_list(workspace, matches->discordant_matches_roots,
		matches->discordant_matches_nodes, matches->starts_forward,
		matches->cigars_forward, matches->starts_reverse,
		matches->cigars_reverse, paired, use_portion, false, dropped_count);
	if (status != 0 || concordant) return status;
	return add_match_list(workspace, matches->concordant_matches_roots,
		matches->concordant_matches_nodes, matches->starts_forward,
		matches->cigars_forward, matches->starts_reverse,
		matches->cigars_reverse, paired, use_portion, true, dropped_count);
}

int candidate_workspace_prepare_scores(CandidateWorkspace *workspace,
	const int *species_counts, size_t number_of_trees)
{
	size_t i;
	size_t required = 0;
	type_of_PP *grown;
	if (workspace == NULL || species_counts == NULL) return -EINVAL;
	for (i = 0; i < workspace->candidate_count; ++i) {
		CandidatePlacement *candidate = &workspace->candidates[i];
		size_t species_count;
		size_t node_count;
		if (candidate->tree_id < 0 || (size_t)candidate->tree_id >= number_of_trees)
			return -ERANGE;
		if (species_counts[candidate->tree_id] <= 0) return -ERANGE;
		species_count = (size_t)species_counts[candidate->tree_id];
		if (species_count > SIZE_MAX / 2 + 1) return -EOVERFLOW;
		node_count = species_count * 2 - 1;
		if (required > SIZE_MAX - node_count) return -EOVERFLOW;
		candidate->score_offset = required;
		candidate->node_count = node_count;
		required += node_count;
	}
	if (required > SIZE_MAX / sizeof(type_of_PP)) return -EOVERFLOW;
	if (required > workspace->score_capacity) {
	#ifdef TRONKO_TESTING
		if (fail_next_growth) {
			fail_next_growth = 0;
			return -ENOMEM;
		}
	#endif
		grown = realloc(workspace->score_arena, required * sizeof(type_of_PP));
		if (grown == NULL && required != 0) return -ENOMEM;
		workspace->score_arena = grown;
		workspace->score_capacity = required;
	}
	workspace->score_length = required;
	if (required != 0)
		memset(workspace->score_arena, 0, required * sizeof(type_of_PP));
	return 0;
}

type_of_PP *candidate_workspace_scores(CandidateWorkspace *workspace,
	size_t candidate_index)
{
	if (workspace == NULL || candidate_index >= workspace->candidate_count)
		return NULL;
	return workspace->score_arena + workspace->candidates[candidate_index].score_offset;
}

const type_of_PP *candidate_workspace_scores_const(
	const CandidateWorkspace *workspace, size_t candidate_index)
{
	if (workspace == NULL || candidate_index >= workspace->candidate_count)
		return NULL;
	return workspace->score_arena + workspace->candidates[candidate_index].score_offset;
}

int candidate_workspace_summarize_winners(CandidateWorkspace *workspace,
	type_of_PP maximum, type_of_PP interval, CandidateLcaFunction get_lca)
{
	size_t candidate_index;
	if (workspace == NULL || interval < 0) return -EINVAL;
	for (candidate_index = 0; candidate_index < workspace->candidate_count;
		++candidate_index) {
		CandidatePlacement *candidate = &workspace->candidates[candidate_index];
		const type_of_PP *scores = candidate_workspace_scores_const(
			workspace, candidate_index);
		size_t node_id;
		if (scores == NULL && candidate->node_count != 0) return -EINVAL;
		candidate->winning_count = 0;
		candidate->lca_node = -1;
		candidate->has_winner = false;
		for (node_id = 0; node_id < candidate->node_count; ++node_id) {
			if (scores[node_id] >= maximum - interval &&
				scores[node_id] <= maximum + interval) {
				candidate->winning_count++;
				if (!candidate->has_winner) {
					candidate->lca_node = (int)node_id;
					candidate->has_winner = true;
				} else {
					if (get_lca == NULL) return -EINVAL;
					candidate->lca_node = get_lca(candidate->lca_node,
						(int)node_id, candidate->tree_id);
				}
			}
		}
	}
	return 0;
}

size_t candidate_workspace_collect_positive(const CandidateWorkspace *workspace,
	size_t indices[MAX_NUM_BWA_MATCHES])
{
	size_t count = 0;
	size_t i;
	if (workspace == NULL || indices == NULL) return 0;
	for (i = 0; i < workspace->candidate_count; ++i) {
		size_t insert_at;
		if (!workspace->candidates[i].has_winner) continue;
		insert_at = count;
		while (insert_at > 0 &&
			workspace->candidates[indices[insert_at - 1]].tree_id >
				workspace->candidates[i].tree_id) {
			indices[insert_at] = indices[insert_at - 1];
			insert_at--;
		}
		indices[insert_at] = i;
		count++;
	}
	return count;
}
