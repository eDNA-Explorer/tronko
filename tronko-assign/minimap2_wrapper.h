#ifndef MINIMAP2_WRAPPER_H
#define MINIMAP2_WRAPPER_H

#include "global.h"

void run_minimap2(int start, int end, bwaMatches *bwa_results, int concordant,
                  int numberOfTrees, char *databaseFile, int paired,
                  int max_query_length, int max_readname_length,
                  int max_acc_name, int max_leaf_matches,
                  int mm2_kmer, int mm2_window);

#endif /* MINIMAP2_WRAPPER_H */
