#include "bwa_source_files/bntseq.h"
#include "bwa_source_files/bwa.h"
#include "bwa_source_files/bwamem.h"
#include "bwa_source_files/bwt.h"

// Forward declaration of main_mem from fastmap.c
int main_mem(char* databaseFile, int number_of_seqs, int number_of_threads,
             bwaMatches* bwa_results, int concordant, int numberOfTrees,
             int startline, int paired, int start, int end,
             int max_query_length, int max_readname_length, int max_acc_name);
