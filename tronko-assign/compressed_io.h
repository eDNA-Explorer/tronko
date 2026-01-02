#ifndef COMPRESSED_IO_H
#define COMPRESSED_IO_H

/**
 * Unified compressed file I/O wrapper for tronko-assign
 *
 * Provides transparent decompression for:
 * - Plain text files
 * - Gzip compressed files (.gz)
 * - Zstandard compressed files (.zst)
 *
 * Format is auto-detected from magic bytes, not file extension.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <zstd.h>

// Compression format types
#define COMPRESS_FORMAT_PLAIN 0
#define COMPRESS_FORMAT_GZIP  1
#define COMPRESS_FORMAT_ZSTD  2

// Unified compressed file handle
typedef struct {
    int format;
    union {
        gzFile gz;          // For gzip and plain (gzopen handles both)
        struct {
            FILE* file;
            ZSTD_DCtx* dctx;
            uint8_t* in_buf;
            uint8_t* out_buf;
            size_t in_buf_size;
            size_t out_buf_size;
            size_t out_pos;
            size_t out_end;
            int eof;
        } zstd;
    } handle;
} CompressedFile;

/**
 * Detect compression format from magic bytes
 *
 * @param filename Path to the file to detect
 * @return COMPRESS_FORMAT_PLAIN, COMPRESS_FORMAT_GZIP, or COMPRESS_FORMAT_ZSTD
 */
static inline int detect_compression_format(const char* filename) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) return COMPRESS_FORMAT_PLAIN;

    uint8_t magic[4];
    size_t bytes_read = fread(magic, 1, 4, fp);
    fclose(fp);

    if (bytes_read < 2) return COMPRESS_FORMAT_PLAIN;

    // Gzip: 0x1f 0x8b
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return COMPRESS_FORMAT_GZIP;
    }

    // Zstd: 0x28 0xb5 0x2f 0xfd
    if (bytes_read >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 &&
        magic[2] == 0x2f && magic[3] == 0xfd) {
        return COMPRESS_FORMAT_ZSTD;
    }

    return COMPRESS_FORMAT_PLAIN;
}

/**
 * Open compressed file for reading
 * Auto-detects format based on magic bytes.
 *
 * @param filename Path to the file to open
 * @param mode Open mode (currently only "r" is supported)
 * @return CompressedFile handle or NULL on failure
 */
static inline CompressedFile* cf_open(const char* filename, const char* mode) {
    (void)mode;  // Currently only read mode is supported

    CompressedFile* cf = (CompressedFile*)calloc(1, sizeof(CompressedFile));
    if (!cf) return NULL;

    cf->format = detect_compression_format(filename);

    if (cf->format == COMPRESS_FORMAT_ZSTD) {
        cf->handle.zstd.file = fopen(filename, "rb");
        if (!cf->handle.zstd.file) {
            free(cf);
            return NULL;
        }
        cf->handle.zstd.dctx = ZSTD_createDCtx();
        if (!cf->handle.zstd.dctx) {
            fclose(cf->handle.zstd.file);
            free(cf);
            return NULL;
        }
        cf->handle.zstd.in_buf_size = ZSTD_DStreamInSize();
        cf->handle.zstd.out_buf_size = ZSTD_DStreamOutSize();
        cf->handle.zstd.in_buf = (uint8_t*)malloc(cf->handle.zstd.in_buf_size);
        cf->handle.zstd.out_buf = (uint8_t*)malloc(cf->handle.zstd.out_buf_size);
        if (!cf->handle.zstd.in_buf || !cf->handle.zstd.out_buf) {
            if (cf->handle.zstd.in_buf) free(cf->handle.zstd.in_buf);
            if (cf->handle.zstd.out_buf) free(cf->handle.zstd.out_buf);
            ZSTD_freeDCtx(cf->handle.zstd.dctx);
            fclose(cf->handle.zstd.file);
            free(cf);
            return NULL;
        }
        cf->handle.zstd.out_pos = 0;
        cf->handle.zstd.out_end = 0;
        cf->handle.zstd.eof = 0;
    } else {
        // gzopen handles both gzip and plain text transparently
        cf->handle.gz = gzopen(filename, "r");
        if (!cf->handle.gz) {
            free(cf);
            return NULL;
        }
    }

    return cf;
}

/**
 * Refill zstd output buffer from compressed input
 * Internal helper function.
 *
 * @param cf CompressedFile handle (must be ZSTD format)
 * @return Number of bytes decompressed, 0 on EOF, -1 on error
 */
static inline int zstd_refill(CompressedFile* cf) {
    if (cf->handle.zstd.eof) return 0;

    size_t bytes_read = fread(cf->handle.zstd.in_buf, 1,
                        cf->handle.zstd.in_buf_size, cf->handle.zstd.file);
    if (bytes_read == 0) {
        cf->handle.zstd.eof = 1;
        return 0;
    }

    ZSTD_inBuffer in = { cf->handle.zstd.in_buf, bytes_read, 0 };
    ZSTD_outBuffer out = { cf->handle.zstd.out_buf, cf->handle.zstd.out_buf_size, 0 };

    while (in.pos < in.size) {
        size_t ret = ZSTD_decompressStream(cf->handle.zstd.dctx, &out, &in);
        if (ZSTD_isError(ret)) {
            fprintf(stderr, "ZSTD decompression error: %s\n", ZSTD_getErrorName(ret));
            return -1;
        }
        if (ret == 0) break; // Frame complete
    }

    cf->handle.zstd.out_pos = 0;
    cf->handle.zstd.out_end = out.pos;
    return (int)out.pos;
}

/**
 * Read line from compressed file (gzgets equivalent)
 * Reads until newline, EOF, or buffer is full.
 *
 * @param buf Buffer to store the line
 * @param size Maximum number of bytes to read (including null terminator)
 * @param cf CompressedFile handle
 * @return buf on success, NULL on EOF or error
 */
static inline char* cf_gets(char* buf, int size, CompressedFile* cf) {
    if (cf->format != COMPRESS_FORMAT_ZSTD) {
        return gzgets(cf->handle.gz, buf, size);
    }

    // zstd streaming read
    int i = 0;
    while (i < size - 1) {
        if (cf->handle.zstd.out_pos >= cf->handle.zstd.out_end) {
            int refill_result = zstd_refill(cf);
            if (refill_result <= 0) break;
        }
        char c = (char)cf->handle.zstd.out_buf[cf->handle.zstd.out_pos++];
        buf[i++] = c;
        if (c == '\n') break;
    }

    if (i == 0) return NULL;
    buf[i] = '\0';
    return buf;
}

/**
 * Close compressed file and free resources
 *
 * @param cf CompressedFile handle
 * @return 0 on success, -1 on error
 */
static inline int cf_close(CompressedFile* cf) {
    if (!cf) return -1;

    if (cf->format == COMPRESS_FORMAT_ZSTD) {
        ZSTD_freeDCtx(cf->handle.zstd.dctx);
        free(cf->handle.zstd.in_buf);
        free(cf->handle.zstd.out_buf);
        fclose(cf->handle.zstd.file);
    } else {
        gzclose(cf->handle.gz);
    }

    free(cf);
    return 0;
}

#endif // COMPRESSED_IO_H
