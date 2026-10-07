#ifndef TRONKO_BWAMEM3_H
#define TRONKO_BWAMEM3_H
#include "global.h"

/* The candidate search of tronko-assign: BWA-MEM3 v0.14.0, run unchanged as a separate program,
 * once per batch (bwamem3.c). */

/* The release tronko-assign accepts (`bwa-mem3 version`, first line). */
#define BWAMEM3_VERSION "0.14.0"

/* Prepares the aligner for this process: finds the bwa-mem3 binary (bin_option, else
 * $TRONKO_BWA_MEM3, else bwa-mem3 beside the tronko-assign executable, else bwa-mem3 on PATH) and
 * checks its version; builds the BWA-MEM3 index of fasta beside it unless skip_build (-6), in
 * which case the index must exist; makes the process's temporary directory; and, when use_shm,
 * stages the index in shared memory. Exits with a message on any failure. */
void bwamem3_init(const char *bin_option, const char *fasta, int threads, int skip_build, int use_shm);

/* Aligns the batch held in singleQueryMat (paired == 0) or pairedQueryMat (paired == 1), entries
 * 0..n-1, and fills results[0..n-1] through the SAM parse. Reorients Tronko's copy of each read as
 * the fork did (rule H8). Exits with a message on any failure. */
void bwamem3_run_batch(int n, int paired, bwaMatches *results, int concordant, int ntree,
                       int max_readname_length, int max_acc_name);

/* Drops the shared-memory index if this process staged it and removes the temporary directory.
 * Registered with atexit by bwamem3_init; safe to call more than once. */
void bwamem3_finish(void);

#endif
