#include "unity.h"

#include <errno.h>
#include <stddef.h>
#include <string.h>

#include "candidate_workspace.h"

static CandidateWorkspace workspace;

void setUp(void)
{
	candidate_workspace_init(&workspace);
	candidate_workspace_reset(&workspace);
}

void tearDown(void)
{
	candidate_workspace_destroy(&workspace);
}

static int lowest_node_lca(int node1, int node2, int tree_id)
{
	(void)tree_id;
	return node1 < node2 ? node1 : node2;
}

static void fill_match_arrays(int roots[MAX_NUM_BWA_MATCHES],
	int nodes[MAX_NUM_BWA_MATCHES])
{
	size_t i;
	for (i = 0; i < MAX_NUM_BWA_MATCHES; ++i) {
		roots[i] = -1;
		nodes[i] = -1;
	}
}

static void test_reset_starts_empty(void)
{
	TEST_ASSERT_EQUAL_UINT64(0, workspace.candidate_count);
	TEST_ASSERT_EQUAL_UINT64(0, workspace.score_length);
}

static void test_add_one_candidate(void)
{
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 3, 17,
		12, "8M", 20, "7M1S"));
	TEST_ASSERT_EQUAL_UINT64(1, workspace.candidate_count);
	TEST_ASSERT_EQUAL_INT(3, workspace.candidates[0].tree_id);
	TEST_ASSERT_EQUAL_INT(17, workspace.candidates[0].leaf_node);
	TEST_ASSERT_EQUAL_STRING("8M", workspace.candidates[0].forward_cigar);
}

static void test_candidate_cap_is_bounded(void)
{
	int i;
	for (i = 0; i < MAX_NUM_BWA_MATCHES; ++i) {
		TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, i, i,
			-1, NULL, -1, NULL));
	}
	TEST_ASSERT_EQUAL_UINT64(MAX_NUM_BWA_MATCHES, workspace.candidate_count);
	TEST_ASSERT_EQUAL_INT(2, candidate_workspace_add(&workspace,
		MAX_NUM_BWA_MATCHES, 99, -1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_UINT64(MAX_NUM_BWA_MATCHES, workspace.candidate_count);
}

static void test_duplicate_tree_keeps_first_hit(void)
{
	int concordant_roots[MAX_NUM_BWA_MATCHES];
	int concordant_nodes[MAX_NUM_BWA_MATCHES];
	int discordant_roots[MAX_NUM_BWA_MATCHES];
	int discordant_nodes[MAX_NUM_BWA_MATCHES];
	bwaMatches matches;
	size_t dropped = 0;

	memset(&matches, 0, sizeof(matches));
	fill_match_arrays(concordant_roots, concordant_nodes);
	fill_match_arrays(discordant_roots, discordant_nodes);
	concordant_roots[0] = 2;
	concordant_nodes[0] = 20;
	concordant_roots[1] = 2;
	concordant_nodes[1] = 21;
	concordant_roots[2] = 1;
	concordant_nodes[2] = 10;
	matches.concordant_matches_roots = concordant_roots;
	matches.concordant_matches_nodes = concordant_nodes;
	matches.discordant_matches_roots = discordant_roots;
	matches.discordant_matches_nodes = discordant_nodes;

	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_select(&workspace, &matches,
		1, 0, 0, &dropped));
	TEST_ASSERT_EQUAL_UINT64(2, workspace.candidate_count);
	TEST_ASSERT_EQUAL_INT(2, workspace.candidates[0].tree_id);
	TEST_ASSERT_EQUAL_INT(20, workspace.candidates[0].leaf_node);
	TEST_ASSERT_EQUAL_INT(1, workspace.candidates[1].tree_id);
	TEST_ASSERT_EQUAL_UINT64(0, dropped);
}

static void test_score_arena_uses_candidate_node_sum_and_reuses_capacity(void)
{
	int species_counts[] = {2, 3, 4};
	type_of_PP *first_arena;
	size_t first_capacity;

	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 0, 2,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 2, 6,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_prepare_scores(&workspace,
		species_counts, 3));
	TEST_ASSERT_EQUAL_UINT64(10, workspace.score_length);
	TEST_ASSERT_EQUAL_UINT64(0, workspace.candidates[0].score_offset);
	TEST_ASSERT_EQUAL_UINT64(3, workspace.candidates[0].node_count);
	TEST_ASSERT_EQUAL_UINT64(3, workspace.candidates[1].score_offset);
	TEST_ASSERT_EQUAL_UINT64(7, workspace.candidates[1].node_count);
	first_arena = workspace.score_arena;
	first_capacity = workspace.score_capacity;
	workspace.score_arena[0] = 42;

	candidate_workspace_reset(&workspace);
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 1, 4,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_prepare_scores(&workspace,
		species_counts, 3));
	TEST_ASSERT_EQUAL_PTR(first_arena, workspace.score_arena);
	TEST_ASSERT_EQUAL_UINT64(first_capacity, workspace.score_capacity);
	TEST_ASSERT_EQUAL_UINT64(5, workspace.score_length);
	TEST_ASSERT_TRUE(workspace.score_arena[0] == (type_of_PP)0);
}

static void test_winner_summary_includes_interval_endpoints_and_sorts_trees(void)
{
	int species_counts[] = {2, 2, 2};
	size_t indices[MAX_NUM_BWA_MATCHES];
	type_of_PP *scores;

	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 2, 2,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 0, 2,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 1, 2,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_prepare_scores(&workspace,
		species_counts, 3));

	scores = candidate_workspace_scores(&workspace, 0);
	scores[0] = 9.0;
	scores[1] = 8.99;
	scores[2] = 11.0;
	scores = candidate_workspace_scores(&workspace, 1);
	scores[0] = 10.0;
	scores = candidate_workspace_scores(&workspace, 2);
	scores[0] = 8.0;

	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_summarize_winners(&workspace,
		10.0, 1.0, lowest_node_lca));
	TEST_ASSERT_EQUAL_INT(2, workspace.candidates[0].winning_count);
	TEST_ASSERT_EQUAL_INT(0, workspace.candidates[0].lca_node);
	TEST_ASSERT_EQUAL_INT(1, workspace.candidates[1].winning_count);
	TEST_ASSERT_FALSE(workspace.candidates[2].has_winner);
	TEST_ASSERT_EQUAL_UINT64(2,
		candidate_workspace_collect_positive(&workspace, indices));
	TEST_ASSERT_EQUAL_INT(0, workspace.candidates[indices[0]].tree_id);
	TEST_ASSERT_EQUAL_INT(2, workspace.candidates[indices[1]].tree_id);
}

static void test_failed_growth_preserves_previous_arena(void)
{
	int species_counts[] = {2, 6};
	type_of_PP *original_arena;
	size_t original_capacity;

	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 0, 2,
		-1, NULL, -1, NULL));
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_prepare_scores(&workspace,
		species_counts, 2));
	original_arena = workspace.score_arena;
	original_capacity = workspace.score_capacity;

	candidate_workspace_reset(&workspace);
	TEST_ASSERT_EQUAL_INT(0, candidate_workspace_add(&workspace, 1, 10,
		-1, NULL, -1, NULL));
	candidate_workspace_fail_next_growth();
	TEST_ASSERT_EQUAL_INT(-ENOMEM, candidate_workspace_prepare_scores(&workspace,
		species_counts, 2));
	TEST_ASSERT_EQUAL_PTR(original_arena, workspace.score_arena);
	TEST_ASSERT_EQUAL_UINT64(original_capacity, workspace.score_capacity);
	TEST_ASSERT_EQUAL_UINT64(0, workspace.score_length);
}

int main(void)
{
	UNITY_BEGIN();
	RUN_TEST(test_reset_starts_empty);
	RUN_TEST(test_add_one_candidate);
	RUN_TEST(test_candidate_cap_is_bounded);
	RUN_TEST(test_duplicate_tree_keeps_first_hit);
	RUN_TEST(test_score_arena_uses_candidate_node_sum_and_reuses_capacity);
	RUN_TEST(test_winner_summary_includes_interval_endpoints_and_sorts_trees);
	RUN_TEST(test_failed_growth_preserves_previous_arena);
	return UNITY_END();
}
