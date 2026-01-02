#ifndef FORMAT_BINARY_H
#define FORMAT_BINARY_H

#include "format_common.h"

tronko_db_t *load_binary(const char *filename, int verbose);
tronko_db_t *load_binary_zstd(const char *filename, int verbose);
int write_binary(tronko_db_t *db, const char *filename, int verbose);
int write_binary_zstd(tronko_db_t *db, const char *filename, int compression_level, int verbose);

#endif
