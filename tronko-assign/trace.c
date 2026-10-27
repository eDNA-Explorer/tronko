#include "trace.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static FILE *trace_fp = NULL;
static char *trace_pattern = NULL;

int trace_init(const char *pattern, const char *path)
{
	if (pattern == NULL || pattern[0] == '\0') {
		return 0;
	}
	if (path == NULL || path[0] == '\0' || strcmp(path, "-") == 0) {
		trace_fp = stderr;
	} else if ((trace_fp = fopen(path, "w")) == NULL) {
		return -1;
	}
	trace_pattern = strdup(pattern);
	fputs("#trace\tread\tcandidates\ttrees\tmaximum\tmax_tree:node\twindow\tvotes"
	      "\tbest_leaf\tbest_leaf_species\tbest_leaf_score"
	      "\tother_species_leaf\tother_species\tother_species_score\tbest_leaf_rule\trow\n", trace_fp);
	return 0;
}

int trace_wanted(const char *name)
{
	if (trace_pattern == NULL) {
		return 0;
	}
	return strcmp(trace_pattern, "*") == 0 || strcmp(trace_pattern, name) == 0;
}

static void put_leaf(FILE *f, int tree, int node, type_of_PP score)
{
	if (tree < 0) {
		fputs("\t-\t-\t-", f);
		return;
	}
	fprintf(f, "\t%d:%d\t%s\t%.17g", tree, node, reduce_leaf_species(tree, node), score);
}

void trace_write(const char *name, int ncand, const int *trees, const reduce_info_t *info, const char *row)
{
	char *line = NULL;
	size_t len = 0;
	FILE *f = open_memstream(&line, &len);
	if (f == NULL) {
		return;
	}
	int i;
	fprintf(f, "trace\t%s\t%d\t", name, ncand);
	for (i = 0; i < ncand; i++) {
		fprintf(f, i ? ",%d" : "%d", trees[i]);
	}
	if (ncand == 0) {
		/* no candidates: no trees and nothing reduced (11 columns) */
		fputs("-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-", f);
	} else {
		fprintf(f, "\t%.17g\t%d:%d\t%.17g\t%d", info->maximum, info->max_tree, info->max_node,
		        info->window, info->votes);
		put_leaf(f, info->l1_tree, info->l1_node, info->l1);
		put_leaf(f, info->l2_tree, info->l2_node, info->l2);
		fprintf(f, "\t%d", info->best_leaf);
	}
	fprintf(f, "\t%s\n", row);
	fclose(f);
	/* One record per call, whole: placement threads finish reads concurrently. */
	flockfile(trace_fp);
	fwrite(line, 1, len, trace_fp);
	funlockfile(trace_fp);
	free(line);
}

void trace_close(void)
{
	if (trace_fp != NULL && trace_fp != stderr) {
		fclose(trace_fp);
	}
	trace_fp = NULL;
	free(trace_pattern);
	trace_pattern = NULL;
}
