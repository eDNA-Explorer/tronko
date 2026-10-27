#ifndef TRONKO_SAM_PARSE_H
#define TRONKO_SAM_PARSE_H
#include "global.h"

/* What the SAM parse needs to know about the batch (sam_parse.c). */
typedef struct {
	bwaMatches *results;     /* one entry per read (single) or pair (paired) of the batch */
	int concordant;
	int ntree;               /* number of trees: the leaf map covers trees 0..ntree-1 */
	int startline;           /* index of the batch's first read in singleQueryMat/pairedQueryMat */
	int paired;
	int max_readname_length;
	int max_acc_name;
} sam_parse_ctx_t;

void sam_parse_batch(const sam_parse_ctx_t *ctx, int n_seqs, char **sam);

/* The tree and node of the leaf named name (a reference name in the SAM), from the leaf map the
 * first batch built; 0 if there is no such leaf. Read-only, so placement threads may call it. */
int sam_leaf_lookup(const char *name, int *tree, int *node);

#endif
