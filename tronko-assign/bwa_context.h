#ifndef TRONKO_BWA_CONTEXT_H
#define TRONKO_BWA_CONTEXT_H

#include <stddef.h>

#include "global.h"
#include "bwa_source_files/bwa.h"

typedef struct tronko_bwa_context tronko_bwa_context;

typedef struct tronko_bwa_chunk_options {
	int number_of_seqs;
	int concordant;
	int startline;
	int paired;
	int start;
	int end;
	int max_query_length;
	int max_readname_length;
	int max_acc_name;
} tronko_bwa_chunk_options;

int tronko_bwa_context_create(const char *database_file,
	int number_of_trees, tronko_bwa_context **out_context);
void tronko_bwa_context_destroy(tronko_bwa_context *context);

int tronko_bwa_align_chunk(const tronko_bwa_context *context,
	bwaMatches *results, const tronko_bwa_chunk_options *options);

/* Read-only accessors used by the embedded BWA adapter. */
const bwaidx_t *tronko_bwa_context_index(const tronko_bwa_context *context);
const leafMap *tronko_bwa_context_lookup_leaf(
	const tronko_bwa_context *context, const char *accession);
size_t tronko_bwa_context_leaf_count(const tronko_bwa_context *context);
size_t tronko_bwa_context_duplicate_count(const tronko_bwa_context *context);

#endif
