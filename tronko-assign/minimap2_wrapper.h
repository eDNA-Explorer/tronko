#ifndef MINIMAP2_WRAPPER_H
#define MINIMAP2_WRAPPER_H

#include "global.h"

/*
 * Map reads using minimap2 and populate bwaMatches struct.
 * Same contract as run_bwa() — caller uses bwaMatches identically regardless
 * of which aligner produced them.
 *
 * The minimap2 index is built once on first call and cached for subsequent calls.
 *
 * Parameters match run_bwa() plus mm2_kmer/mm2_window for minimap2 tuning.
 */
void run_minimap2(int start, int end, bwaMatches *bwa_results, int concordant,
                  int numberOfTrees, char *databaseFile, int paired,
                  int max_query_length, int max_readname_length,
                  int max_acc_name, int max_bwa_matches,
                  int mm2_kmer, int mm2_window);

#endif /* MINIMAP2_WRAPPER_H */
