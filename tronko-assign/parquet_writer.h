/*
 * parquet_writer.h
 * Parquet output support for tronko-assign using Carquet library
 */
#ifndef _PARQUET_WRITER_H_
#define _PARQUET_WRITER_H_

#ifdef ENABLE_PARQUET

#include <stdint.h>

/* Forward declaration - actual struct defined in parquet_writer.c */
typedef struct parquet_writer parquet_writer_t;

/* Forward declaration of assignmentResult from global.h */
struct assignmentResult;

/*
 * Create a new Parquet writer for assignment results
 *
 * @param filename  Output filename (e.g., "results_0.parquet")
 * @param err       Error message buffer (caller-allocated, min 256 bytes)
 * @return          Writer handle, or NULL on error
 */
parquet_writer_t* parquet_writer_create(const char *filename, char *err);

/*
 * Write a batch of assignment results
 *
 * @param writer    Writer handle from parquet_writer_create
 * @param results   Array of assignment results
 * @param count     Number of results in array
 * @return          0 on success, -1 on error
 */
int parquet_writer_write_batch(parquet_writer_t *writer,
                                const struct assignmentResult *results,
                                int count);

/*
 * Close writer and finalize Parquet file
 * Writes footer metadata and releases resources
 *
 * @param writer    Writer handle (freed after call)
 * @return          0 on success, -1 on error
 */
int parquet_writer_close(parquet_writer_t *writer);

#endif /* ENABLE_PARQUET */
#endif /* _PARQUET_WRITER_H_ */
