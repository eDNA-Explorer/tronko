#ifndef _READ_REF_
#define _READ_REF_

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <dirent.h>
#include <regex.h>
#include "global.h"
#include "allocatetreememory.h"
int readInXNumberOfLines_fastq(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length,int first_iter);
int readInXNumberOfLines(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length);
void shiftUp(int iter, int jump, int numberOfLinesToRead);
int readReferenceTree( gzFile referenceTree,int* name_specs);
int setNumbase_setNumspec(int numberOfPartitions,int* specs);
void find_specs_for_reads(int* specs, gzFile file, int format);

// Format detection for reference files
int detect_reference_format(const char *filename);

// Binary format reader (mirrors readReferenceTree but for .trkb files)
int readReferenceBinary(const char *filename, int *name_specs);

// Gzipped binary format reader
int readReferenceBinaryGzipped(const char *filename, int *name_specs);

#endif /* _READ_REF_ */
