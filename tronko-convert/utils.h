#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

// Error-checked I/O (modeled after BWA's utils.c:124-174)
FILE *err_xopen(const char *fn, const char *mode);
gzFile err_xzopen(const char *fn, const char *mode);
size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t err_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
int err_fseek(FILE *stream, long offset, int whence);
long err_ftell(FILE *stream);
int err_fclose(FILE *stream);
int err_gzclose(gzFile fp);

// Little-endian write helpers
void write_u16(FILE *fp, uint16_t v);
void write_u32(FILE *fp, uint32_t v);
void write_u64(FILE *fp, uint64_t v);
void write_i32(FILE *fp, int32_t v);
void write_float(FILE *fp, float v);

// Little-endian read helpers
uint16_t read_u16(FILE *fp);
uint32_t read_u32(FILE *fp);
uint64_t read_u64(FILE *fp);
int32_t read_i32(FILE *fp);
float read_float(FILE *fp);

#endif
