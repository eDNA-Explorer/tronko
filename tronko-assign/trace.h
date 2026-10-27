/* --trace-read: one diagnostic record per traced read, written as the read is finished.
 * Read-only: it never changes what is scored or assigned, and it is off (one branch per read)
 * unless asked for. */
#ifndef TRACE_H
#define TRACE_H

#include "reduce.h"

/* Turn tracing on for reads named *pattern* ("*" = every read), writing to *path* ("" or "-"
 * = stderr). Returns 0, or -1 if the file cannot be opened. Called once in main. */
int trace_init(const char *pattern, const char *path);

/* 1 when tracing is on and *name* is traced. */
int trace_wanted(const char *name);

/* Write the record of one read: its candidates, what the reduction found (when it had
 * candidates) and its output row. Each record is one line, written whole. */
void trace_write(const char *name, int ncand, const int *trees, const reduce_info_t *info, const char *row);

void trace_close(void);

#endif
