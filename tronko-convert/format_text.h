#ifndef FORMAT_TEXT_H
#define FORMAT_TEXT_H

#include "format_common.h"

tronko_db_t *load_text(const char *filename, int verbose);
int write_text(tronko_db_t *db, const char *filename, int verbose);

#endif
