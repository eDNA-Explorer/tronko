/*
 * minimap2_shim.c
 *
 * Provides minimap2-specific global variables and utility functions that
 * would normally come from minimap2's misc.c, but can't be compiled from
 * misc.c directly because it also defines cputime() and realtime() which
 * conflict with BWA's utils.c.
 *
 * This file is compiled as part of the main tronko-assign binary.
 * The minimap2 library (libminimap2.a) resolves these symbols at link time.
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "minimap2_src/minimap.h"
#include "minimap2_src/mmpriv.h"
#include "minimap2_src/ksort.h"

/* Global variables needed by minimap2 */
int mm_verbose = 1;
int mm_dbg_flag = 0;
double mm_realtime0;

/* peakrss - not defined in BWA */
long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

/* Error-checked I/O wrappers used by minimap2 internally */
void mm_err_puts(const char *str)
{
	int ret;
	ret = puts(str);
	if (ret == EOF) {
		perror("[ERROR] failed to write the results");
		exit(EXIT_FAILURE);
	}
}

void mm_err_fwrite(const void *p, size_t size, size_t nitems, FILE *fp)
{
	size_t ret;
	ret = fwrite(p, size, nitems, fp);
	if (ret != nitems) {
		perror("[ERROR] failed to write data");
		exit(EXIT_FAILURE);
	}
}

void mm_err_fread(void *p, size_t size, size_t nitems, FILE *fp)
{
	size_t ret;
	ret = fread(p, size, nitems, fp);
	if (ret != nitems) {
		perror("[ERROR] failed to read data");
		exit(EXIT_FAILURE);
	}
}

/* Radix sort instantiations used by minimap2 internally (from misc.c) */
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8)

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT_GENERIC(uint64_t)
