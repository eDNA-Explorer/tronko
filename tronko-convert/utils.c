#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>

static void fatal(const char *func, const char *msg) {
    fprintf(stderr, "[FATAL] %s: %s\n", func, msg);
    exit(1);
}

FILE *err_xopen(const char *fn, const char *mode) {
    FILE *fp = fopen(fn, mode);
    if (!fp) {
        fprintf(stderr, "[FATAL] Cannot open '%s': %s\n", fn, strerror(errno));
        exit(1);
    }
    return fp;
}

gzFile err_xzopen(const char *fn, const char *mode) {
    gzFile fp = gzopen(fn, mode);
    if (!fp) {
        fprintf(stderr, "[FATAL] Cannot open '%s': %s\n", fn,
                errno ? strerror(errno) : "Out of memory");
        exit(1);
    }
    return fp;
}

size_t err_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t ret = fwrite(ptr, size, nmemb, stream);
    if (ret != nmemb) fatal("fwrite", strerror(errno));
    return ret;
}

size_t err_fread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t ret = fread(ptr, size, nmemb, stream);
    if (ret != nmemb) {
        fatal("fread", ferror(stream) ? strerror(errno) : "Unexpected end of file");
    }
    return ret;
}

int err_fseek(FILE *stream, long offset, int whence) {
    int ret = fseek(stream, offset, whence);
    if (ret != 0) fatal("fseek", strerror(errno));
    return ret;
}

long err_ftell(FILE *stream) {
    long ret = ftell(stream);
    if (ret == -1) fatal("ftell", strerror(errno));
    return ret;
}

int err_fclose(FILE *stream) {
    int ret = fclose(stream);
    if (ret != 0) fatal("fclose", strerror(errno));
    return ret;
}

int err_gzclose(gzFile fp) {
    int ret = gzclose(fp);
    if (ret != Z_OK) fatal("gzclose", "gzip close failed");
    return ret;
}

// Little-endian write helpers
void write_u16(FILE *fp, uint16_t v) {
    uint8_t b[2] = { v & 0xFF, (v >> 8) & 0xFF };
    err_fwrite(b, 1, 2, fp);
}

void write_u32(FILE *fp, uint32_t v) {
    uint8_t b[4] = { v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF, (v >> 24) & 0xFF };
    err_fwrite(b, 1, 4, fp);
}

void write_u64(FILE *fp, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (v >> (i * 8)) & 0xFF;
    err_fwrite(b, 1, 8, fp);
}

void write_i32(FILE *fp, int32_t v) {
    write_u32(fp, (uint32_t)v);
}

void write_float(FILE *fp, float v) {
    err_fwrite(&v, sizeof(float), 1, fp);
}

// Little-endian read helpers
uint16_t read_u16(FILE *fp) {
    uint8_t b[2];
    err_fread(b, 1, 2, fp);
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

uint32_t read_u32(FILE *fp) {
    uint8_t b[4];
    err_fread(b, 1, 4, fp);
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

uint64_t read_u64(FILE *fp) {
    uint8_t b[8];
    err_fread(b, 1, 8, fp);
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) v |= ((uint64_t)b[i] << (i * 8));
    return v;
}

int32_t read_i32(FILE *fp) {
    return (int32_t)read_u32(fp);
}

float read_float(FILE *fp) {
    float v;
    err_fread(&v, sizeof(float), 1, fp);
    return v;
}
