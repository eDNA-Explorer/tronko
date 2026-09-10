#include "bwa_context.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "hashmap.h"

struct tronko_bwa_context {
	bwaidx_t *index;
	HASHMAP(char, leafMap) leaf_map;
	leafMap *coordinates;
	size_t coordinate_count;
	size_t unique_count;
	size_t duplicate_count;
	int map_initialized;
};

static void destroy_partial_context(tronko_bwa_context *context)
{
	if (context == NULL) return;
	if (context->map_initialized) hashmap_cleanup(&context->leaf_map);
	free(context->coordinates);
	bwa_idx_destroy(context->index);
	free(context);
}

int tronko_bwa_context_create(const char *database_file,
	int number_of_trees, tronko_bwa_context **out_context)
{
	tronko_bwa_context *context;
	size_t leaf_count = 0;
	size_t coordinate_index = 0;
	int tree_id;

	if (database_file == NULL || out_context == NULL || number_of_trees <= 0)
		return -EINVAL;
	*out_context = NULL;

	for (tree_id = 0; tree_id < number_of_trees; ++tree_id) {
		if (numspecArr[tree_id] < 0 ||
			leaf_count > SIZE_MAX - (size_t)numspecArr[tree_id])
			return -EOVERFLOW;
		leaf_count += (size_t)numspecArr[tree_id];
	}
	if (leaf_count > SIZE_MAX / sizeof(leafMap)) return -EOVERFLOW;

	context = calloc(1, sizeof(*context));
	if (context == NULL) return -ENOMEM;

	/* BWA diagnostics are process configuration, never worker-owned state. */
	bwa_verbose = 1;
	context->index = bwa_idx_load_from_shm(database_file);
	if (context->index == NULL)
		context->index = bwa_idx_load(database_file, BWA_IDX_ALL);
	if (context->index == NULL) {
		destroy_partial_context(context);
		return -ENOENT;
	}

	hashmap_init(&context->leaf_map, hashmap_hash_string, strcmp);
	context->map_initialized = 1;
	if (hashmap_reserve(&context->leaf_map, leaf_count) != 0) {
		destroy_partial_context(context);
		return -ENOMEM;
	}

	context->coordinates = calloc(leaf_count, sizeof(*context->coordinates));
	if (context->coordinates == NULL && leaf_count != 0) {
		destroy_partial_context(context);
		return -ENOMEM;
	}
	context->coordinate_count = leaf_count;

	for (tree_id = 0; tree_id < number_of_trees; ++tree_id) {
		int node_id;
		for (node_id = numspecArr[tree_id] - 1;
			node_id < 2 * numspecArr[tree_id] - 1; ++node_id) {
			leafMap *coordinate = &context->coordinates[coordinate_index++];
			int status;
			coordinate->name = treeArr[tree_id][node_id].name;
			coordinate->root = tree_id;
			coordinate->node = node_id;
			status = hashmap_put(&context->leaf_map, coordinate->name, coordinate);
			if (status == 0) {
				context->unique_count++;
			} else if (status == -EEXIST) {
				/* Preserve the baseline's first-seen coordinate. */
				context->duplicate_count++;
			} else {
				destroy_partial_context(context);
				return status;
			}
		}
	}

	*out_context = context;
	return 0;
}

void tronko_bwa_context_destroy(tronko_bwa_context *context)
{
	destroy_partial_context(context);
}

const bwaidx_t *tronko_bwa_context_index(const tronko_bwa_context *context)
{
	return context == NULL ? NULL : context->index;
}

const leafMap *tronko_bwa_context_lookup_leaf(
	const tronko_bwa_context *context, const char *accession)
{
	if (context == NULL || accession == NULL) return NULL;
	return hashmap_base_get(&context->leaf_map.map_base, accession);
}

size_t tronko_bwa_context_leaf_count(const tronko_bwa_context *context)
{
	return context == NULL ? 0 : context->unique_count;
}

size_t tronko_bwa_context_duplicate_count(const tronko_bwa_context *context)
{
	return context == NULL ? 0 : context->duplicate_count;
}
