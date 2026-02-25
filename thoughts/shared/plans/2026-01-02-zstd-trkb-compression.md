# ZSTD Compression for TRKB Binary Format Implementation Plan

## Overview

Add zstd compression support to the TRKB binary format:
- **tronko-convert**: Output zstd-compressed `.trkb` files by default
- **tronko-assign**: Read zstd-compressed `.trkb` files with auto-detection

This provides 3.9x smaller files than gzipped text and 2.5x faster decompression than gzip.

## Current State Analysis

### tronko-convert
- Outputs uncompressed binary via `write_binary()` in `format_binary.c:16`
- Links `-lz -lm` only (no zstd) - `Makefile:4`
- Uses `FILE*` with `fwrite()` for sequential output

### tronko-assign
- Complete zstd streaming infrastructure in `compressed_io.h:1-218`
- Already links `-lzstd` in `Makefile:12`
- Has `readReferenceBinary()` (line 701) and `readReferenceBinaryGzipped()` (line 991)
- Format detection in `detect_reference_format()` at line 89

### Key Discoveries:
- `compressed_io.h` provides `cf_open()`, `cf_gets()`, `cf_close()` for text-based streaming
- Binary reading needs bulk `fread()`/`gzread()` equivalent - not provided by current `cf_*` functions
- Gzip reader uses `gz_read_*` helper functions for little-endian binary values

## Desired End State

After implementation:
1. `tronko-convert -i ref.txt -o ref.trkb` produces zstd-compressed output by default
2. `tronko-convert -i ref.txt -o ref.trkb -u` produces uncompressed output
3. `tronko-assign -f ref.trkb ...` auto-detects and reads zstd-compressed files
4. All existing functionality preserved (text, gzip support)

### Verification:
```bash
# Round-trip test
tronko-convert -i ref.txt -o ref.trkb
tronko-convert -i ref.trkb -o roundtrip.txt -t
diff ref.txt roundtrip.txt

# Assignment test
tronko-assign -f ref.trkb -s -g query.fasta -o results.txt
```

## What We're NOT Doing

- Separate `.trkb.zst` extension (use magic detection)
- Streaming writes in tronko-build (future enhancement)
- Memory-mapped loading (future enhancement)

## Design Note: Compression Level

The `-c` compression level flag only applies to **tronko-convert** (the writer). Zstd decompression doesn't require knowing the compression level - the format is self-describing. This means:
- Users of tronko-assign don't need to specify anything about compression
- tronko-assign will correctly decompress files regardless of what level was used to create them

## Implementation Approach

Use streaming compression/decompression for memory efficiency with large databases. The zstd streaming API (`ZSTD_compressStream2`/`ZSTD_decompressStream`) allows processing data in chunks without holding the entire file in memory.

---

## Phase 1: tronko-convert Writer Updates

### Overview
Add zstd compression as the default output format for tronko-convert.

### Changes Required:

#### 1. Makefile
**File**: `tronko-convert/Makefile`
**Changes**: Add `-lzstd` to linker flags

```makefile
LDFLAGS = -lz -lm -lzstd
```

#### 2. format_common.h Constants
**File**: `tronko-convert/format_common.h`
**Changes**: Add zstd-related constants

```c
// Add after line 19 (FORMAT_BINARY_GZIPPED definition)
#define FORMAT_BINARY_ZSTD     4  // Zstd-compressed binary format

// Add after TAXONOMY_LEVELS
#define ZSTD_COMPRESSION_LEVEL 19  // Max compression (write-once, read-many)
```

#### 3. format_binary.h Header
**File**: `tronko-convert/format_binary.h`
**Changes**: Add function prototype for zstd writer

```c
// Add after write_binary() declaration
int write_binary_zstd(tronko_db_t *db, const char *filename, int compression_level, int verbose);
```

#### 4. format_binary.c Implementation
**File**: `tronko-convert/format_binary.c`
**Changes**: Add `write_binary_zstd()` function using streaming compression

```c
#include <zstd.h>

// Streaming zstd writer context
typedef struct {
    FILE *fp;
    ZSTD_CCtx *cctx;
    uint8_t *out_buf;
    size_t out_buf_size;
} zstd_writer_t;

static zstd_writer_t *zstd_writer_open(const char *filename, int level) {
    zstd_writer_t *w = calloc(1, sizeof(zstd_writer_t));
    if (!w) return NULL;

    w->fp = fopen(filename, "wb");
    if (!w->fp) { free(w); return NULL; }

    w->cctx = ZSTD_createCCtx();
    if (!w->cctx) { fclose(w->fp); free(w); return NULL; }

    ZSTD_CCtx_setParameter(w->cctx, ZSTD_c_compressionLevel, level);
    ZSTD_CCtx_setParameter(w->cctx, ZSTD_c_checksumFlag, 1);

    w->out_buf_size = ZSTD_CStreamOutSize();
    w->out_buf = malloc(w->out_buf_size);
    if (!w->out_buf) {
        ZSTD_freeCCtx(w->cctx);
        fclose(w->fp);
        free(w);
        return NULL;
    }

    return w;
}

static int zstd_writer_write(zstd_writer_t *w, const void *data, size_t size) {
    ZSTD_inBuffer in = { data, size, 0 };

    while (in.pos < in.size) {
        ZSTD_outBuffer out = { w->out_buf, w->out_buf_size, 0 };
        size_t ret = ZSTD_compressStream2(w->cctx, &out, &in, ZSTD_e_continue);
        if (ZSTD_isError(ret)) {
            fprintf(stderr, "ZSTD compression error: %s\n", ZSTD_getErrorName(ret));
            return -1;
        }
        if (out.pos > 0) {
            if (fwrite(w->out_buf, 1, out.pos, w->fp) != out.pos) return -1;
        }
    }
    return 0;
}

static int zstd_writer_close(zstd_writer_t *w) {
    if (!w) return -1;

    // Flush remaining data
    ZSTD_inBuffer in = { NULL, 0, 0 };
    size_t ret;
    do {
        ZSTD_outBuffer out = { w->out_buf, w->out_buf_size, 0 };
        ret = ZSTD_compressStream2(w->cctx, &out, &in, ZSTD_e_end);
        if (ZSTD_isError(ret)) {
            fprintf(stderr, "ZSTD flush error: %s\n", ZSTD_getErrorName(ret));
            break;
        }
        if (out.pos > 0) {
            fwrite(w->out_buf, 1, out.pos, w->fp);
        }
    } while (ret > 0);

    ZSTD_freeCCtx(w->cctx);
    free(w->out_buf);
    fclose(w->fp);
    free(w);
    return 0;
}

// Helper to write little-endian values to zstd stream
static int zstd_write_u16(zstd_writer_t *w, uint16_t v) {
    uint8_t b[2] = { v & 0xFF, (v >> 8) & 0xFF };
    return zstd_writer_write(w, b, 2);
}

static int zstd_write_u32(zstd_writer_t *w, uint32_t v) {
    uint8_t b[4] = { v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF, (v >> 24) & 0xFF };
    return zstd_writer_write(w, b, 4);
}

static int zstd_write_i32(zstd_writer_t *w, int32_t v) {
    return zstd_write_u32(w, (uint32_t)v);
}

static int zstd_write_u64(zstd_writer_t *w, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (v >> (i * 8)) & 0xFF;
    return zstd_writer_write(w, b, 8);
}

int write_binary_zstd(tronko_db_t *db, const char *filename, int compression_level, int verbose) {
    zstd_writer_t *w = zstd_writer_open(filename, compression_level);
    if (!w) {
        fprintf(stderr, "Error: Cannot open '%s' for zstd writing\n", filename);
        return -1;
    }

    if (verbose) {
        fprintf(stderr, "  Writing zstd-compressed binary format v%d.%d (level %d)\n",
                TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR, compression_level);
    }

    // Calculate section offsets (same as write_binary)
    uint64_t header_end = 64 + 16;  // FILE_HEADER_SIZE + GLOBAL_META_SIZE
    uint64_t tree_meta_size = db->num_trees * 12;  // TREE_META_SIZE
    uint64_t taxonomy_offset = header_end + tree_meta_size;

    // Calculate taxonomy section size
    uint64_t taxonomy_size = db->num_trees * 8;
    for (int t = 0; t < db->num_trees; t++) {
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                taxonomy_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
    }

    uint64_t node_offset = taxonomy_offset + taxonomy_size;

    // Calculate node section size
    uint64_t node_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        node_size += 4;
        node_size += db->trees[t].num_nodes * 32;  // NODE_RECORD_SIZE
        for (int n = db->trees[t].numspec - 1; n < db->trees[t].num_nodes; n++) {
            if (db->trees[t].nodes[n].name) {
                node_size += 2 + strlen(db->trees[t].nodes[n].name) + 1;
            }
        }
    }

    uint64_t posterior_offset = node_offset + node_size;

    // Calculate posterior section size
    uint64_t posterior_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        posterior_size += (uint64_t)db->trees[t].num_nodes * db->trees[t].numbase * 4 * sizeof(float);
    }

    uint64_t total_size = posterior_offset + posterior_size + 8;

    if (verbose) {
        fprintf(stderr, "  Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset);
    }

    // === Write File Header (64 bytes) ===
    uint8_t header_bytes[8] = {
        TRONKO_MAGIC_0, TRONKO_MAGIC_1, TRONKO_MAGIC_2, TRONKO_MAGIC_3,
        TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR,
        0x01,  // Little-endian
        0x01   // Float precision
    };
    zstd_writer_write(w, header_bytes, 8);

    // Header CRC (of bytes 0-7)
    uint32_t header_crc = tronko_crc32(0, header_bytes, 8);
    zstd_write_u32(w, header_crc);
    zstd_write_u32(w, 0);  // Reserved

    zstd_write_u64(w, taxonomy_offset);
    zstd_write_u64(w, node_offset);
    zstd_write_u64(w, posterior_offset);
    zstd_write_u64(w, total_size);

    // Padding to 64 bytes
    uint8_t reserved[16] = {0};
    zstd_writer_write(w, reserved, 16);

    // === Write Global Metadata (16 bytes) ===
    zstd_write_i32(w, db->num_trees);
    zstd_write_i32(w, db->max_nodename);
    zstd_write_i32(w, db->max_tax_name);
    zstd_write_i32(w, db->max_line_taxonomy);

    // === Write Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < db->num_trees; i++) {
        zstd_write_i32(w, db->trees[i].numbase);
        zstd_write_i32(w, db->trees[i].root);
        zstd_write_i32(w, db->trees[i].numspec);
    }

    // === Write Taxonomy Section ===
    for (int t = 0; t < db->num_trees; t++) {
        uint32_t tree_tax_size = 0;
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                tree_tax_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
        zstd_write_u32(w, tree_tax_size);
        zstd_write_u32(w, 0);  // Reserved

        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                const char *str = db->trees[t].taxonomy[s][l];
                uint16_t len = strlen(str) + 1;
                zstd_write_u16(w, len);
                zstd_writer_write(w, str, len);
            }
        }
    }

    // === Write Node Section ===
    for (int t = 0; t < db->num_trees; t++) {
        zstd_write_u32(w, db->trees[t].num_nodes);

        uint32_t name_offset = 0;

        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            zstd_write_i32(w, node->up[0]);
            zstd_write_i32(w, node->up[1]);
            zstd_write_i32(w, node->down);
            zstd_write_i32(w, node->depth);
            zstd_write_i32(w, node->taxIndex[0]);
            zstd_write_i32(w, node->taxIndex[1]);

            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                zstd_write_u32(w, name_offset + 1);
                name_offset += 2 + strlen(node->name) + 1;
            } else {
                zstd_write_u32(w, 0);
            }
            zstd_write_u32(w, 0);  // Reserved
        }

        // Write name table
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                uint16_t len = strlen(node->name) + 1;
                zstd_write_u16(w, len);
                zstd_writer_write(w, node->name, len);
            }
        }
    }

    // === Write Posterior Section ===
    if (verbose) {
        fprintf(stderr, "  Writing posteriors...\n");
    }

    for (int t = 0; t < db->num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            zstd_writer_write(w, db->trees[t].nodes[n].posteriors,
                             db->trees[t].numbase * 4 * sizeof(float));
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: wrote %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // === Write Footer ===
    // Note: For zstd, we skip the CRC calculation since zstd has its own checksumming
    zstd_write_u32(w, 0);  // data_crc placeholder
    zstd_write_u32(w, TRONKO_FOOTER_MAGIC);

    zstd_writer_close(w);

    if (verbose) {
        fprintf(stderr, "  Zstd compression complete\n");
    }

    return 0;
}
```

#### 5. tronko-convert.c CLI Updates
**File**: `tronko-convert/tronko-convert.c`
**Changes**: Make zstd default, add `-u` flag for uncompressed

```c
// Update print_usage()
static void print_usage(const char *prog) {
    fprintf(stderr, "tronko-convert: Convert tronko reference database formats\n\n");
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "  -i <file>   Input file (required)\n");
    fprintf(stderr, "  -o <file>   Output file (required)\n");
    fprintf(stderr, "  -t          Output as text (default: zstd-compressed binary)\n");
    fprintf(stderr, "  -u          Output as uncompressed binary\n");
    fprintf(stderr, "  -c <level>  Zstd compression level 1-19 (default: 19, max compression)\n");
    fprintf(stderr, "              Lower = faster compression, larger file\n");
    fprintf(stderr, "              Higher = slower compression, smaller file\n");
    fprintf(stderr, "  -v          Verbose output\n");
    fprintf(stderr, "  -h          Show this help\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb\n", prog);
    fprintf(stderr, "  %s -i reference_tree.trkb -o reference_tree.txt -t\n", prog);
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb -u  # uncompressed\n", prog);
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb -c 3  # fast compression\n", prog);
}

// Update main()
int main(int argc, char *argv[]) {
    char *input_file = NULL;
    char *output_file = NULL;
    int output_text = 0;
    int output_uncompressed = 0;
    int compression_level = ZSTD_COMPRESSION_LEVEL;  // Default: 19
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "i:o:tuc:vh")) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'o': output_file = optarg; break;
            case 't': output_text = 1; break;
            case 'u': output_uncompressed = 1; break;
            case 'c':
                compression_level = atoi(optarg);
                if (compression_level < 1 || compression_level > 19) {
                    fprintf(stderr, "Error: Compression level must be 1-19\n");
                    return 1;
                }
                break;
            case 'v': verbose = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    // ... (input validation unchanged)

    // Update output format reporting
    if (verbose) {
        fprintf(stderr, "Input format: %s\n",
                input_format == FORMAT_BINARY ? "binary (.trkb)" :
                input_format == FORMAT_BINARY_ZSTD ? "zstd binary (.trkb)" : "text");
        fprintf(stderr, "Output format: %s\n",
                output_text ? "text" :
                output_uncompressed ? "uncompressed binary (.trkb)" :
                "zstd binary (.trkb)");
    }

    // ... (loading unchanged)

    // Update write logic
    int result;
    if (output_text) {
        result = write_text(db, output_file, verbose);
    } else if (output_uncompressed) {
        result = write_binary(db, output_file, verbose);  // Existing function
    } else {
        result = write_binary_zstd(db, output_file, compression_level, verbose);  // NEW default
    }

    // ... (rest unchanged)
}
```

### Success Criteria:

#### Automated Verification:
- [x] `make clean && make` succeeds in tronko-convert directory
- [x] `./tronko-convert -i ../tronko-build/example_datasets/single_tree/reference_tree.txt -o /tmp/test.trkb -v` creates zstd-compressed file
- [x] `file /tmp/test.trkb` reports "Zstandard compressed data"
- [x] `./tronko-convert -i ../tronko-build/example_datasets/single_tree/reference_tree.txt -o /tmp/test_uncomp.trkb -u -v` creates uncompressed file
- [x] `./tronko-convert -i ../tronko-build/example_datasets/single_tree/reference_tree.txt -o /tmp/test_fast.trkb -c 3 -v` creates zstd file with level 3
- [x] File size comparison: level 19 < level 3 < uncompressed
- [x] Invalid compression level (`-c 25`) produces error

#### Manual Verification:
- [x] `./tronko-convert -h` shows updated usage with `-u` and `-c` flags

---

## Phase 2: tronko-assign Reader Updates

### Overview
Add zstd decompression support for reading `.trkb` files in tronko-assign.

### Changes Required:

#### 1. global.h Format Constant
**File**: `tronko-assign/global.h`
**Changes**: Add FORMAT_BINARY_ZSTD constant

```c
// Update after line 54 (FORMAT_BINARY_GZIPPED)
#define FORMAT_BINARY_ZSTD     4  // Zstd-compressed binary format (.trkb)
```

#### 2. readreference.h Header
**File**: `tronko-assign/readreference.h`
**Changes**: Add function prototype

```c
// Add prototype for zstd reader
int readReferenceBinaryZstd(const char *filename, int *name_specs);
```

#### 3. readreference.c Format Detection
**File**: `tronko-assign/readreference.c`
**Changes**: Update `detect_reference_format()` to detect zstd-compressed TRKB

```c
// In detect_reference_format(), after line 130 (gzip detection):

    // Check for zstd magic: 0x28 0xb5 0x2f 0xfd
    if (bytes_read >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 &&
        magic[2] == 0x2f && magic[3] == 0xfd) {
        // Zstd compressed - decompress first bytes to check for TRKB magic
        FILE *zfp = fopen(filename, "rb");
        if (zfp) {
            // Read compressed data
            size_t comp_size = 256;  // Small buffer for header detection
            uint8_t *comp_buf = malloc(comp_size);
            size_t decomp_size = 64;
            uint8_t *decomp_buf = malloc(decomp_size);

            if (comp_buf && decomp_buf) {
                size_t n_read = fread(comp_buf, 1, comp_size, zfp);
                size_t result = ZSTD_decompress(decomp_buf, decomp_size, comp_buf, n_read);

                if (!ZSTD_isError(result) && result >= 4) {
                    if (decomp_buf[0] == TRONKO_MAGIC_0 && decomp_buf[1] == TRONKO_MAGIC_1 &&
                        decomp_buf[2] == TRONKO_MAGIC_2 && decomp_buf[3] == TRONKO_MAGIC_3) {
                        LOG_DEBUG("Detected zstd-compressed binary format: %s", filename);
                        free(comp_buf);
                        free(decomp_buf);
                        fclose(zfp);
                        return FORMAT_BINARY_ZSTD;
                    }
                }
            }

            if (comp_buf) free(comp_buf);
            if (decomp_buf) free(decomp_buf);
            fclose(zfp);
        }
        LOG_WARN("Zstd file does not contain TRKB data: %s", filename);
        return FORMAT_UNKNOWN;
    }
```

#### 4. readreference.c Zstd Reader Implementation
**File**: `tronko-assign/readreference.c`
**Changes**: Add `readReferenceBinaryZstd()` function (mirroring `readReferenceBinaryGzipped()`)

```c
// Zstd streaming read helpers
typedef struct {
    FILE *file;
    ZSTD_DCtx *dctx;
    uint8_t *in_buf;
    uint8_t *out_buf;
    size_t in_buf_size;
    size_t out_buf_size;
    size_t in_pos;
    size_t in_end;
    size_t out_pos;
    size_t out_end;
    int eof;
} zstd_reader_t;

static zstd_reader_t *zstd_reader_open(const char *filename) {
    zstd_reader_t *r = calloc(1, sizeof(zstd_reader_t));
    if (!r) return NULL;

    r->file = fopen(filename, "rb");
    if (!r->file) { free(r); return NULL; }

    r->dctx = ZSTD_createDCtx();
    if (!r->dctx) { fclose(r->file); free(r); return NULL; }

    r->in_buf_size = ZSTD_DStreamInSize();
    r->out_buf_size = ZSTD_DStreamOutSize();
    r->in_buf = malloc(r->in_buf_size);
    r->out_buf = malloc(r->out_buf_size);

    if (!r->in_buf || !r->out_buf) {
        if (r->in_buf) free(r->in_buf);
        if (r->out_buf) free(r->out_buf);
        ZSTD_freeDCtx(r->dctx);
        fclose(r->file);
        free(r);
        return NULL;
    }

    return r;
}

static int zstd_reader_refill(zstd_reader_t *r) {
    if (r->eof) return 0;

    // Move remaining input data to start
    if (r->in_pos > 0 && r->in_pos < r->in_end) {
        memmove(r->in_buf, r->in_buf + r->in_pos, r->in_end - r->in_pos);
        r->in_end -= r->in_pos;
        r->in_pos = 0;
    } else {
        r->in_pos = 0;
        r->in_end = 0;
    }

    // Read more compressed data
    size_t space = r->in_buf_size - r->in_end;
    if (space > 0) {
        size_t n = fread(r->in_buf + r->in_end, 1, space, r->file);
        r->in_end += n;
        if (n == 0 && r->in_end == 0) {
            r->eof = 1;
            return 0;
        }
    }

    ZSTD_inBuffer in = { r->in_buf, r->in_end, r->in_pos };
    ZSTD_outBuffer out = { r->out_buf, r->out_buf_size, 0 };

    size_t ret = ZSTD_decompressStream(r->dctx, &out, &in);
    if (ZSTD_isError(ret)) {
        LOG_ERROR("ZSTD decompression error: %s", ZSTD_getErrorName(ret));
        return -1;
    }

    r->in_pos = in.pos;
    r->out_pos = 0;
    r->out_end = out.pos;

    return (int)out.pos;
}

static size_t zstd_reader_read(zstd_reader_t *r, void *buf, size_t size) {
    size_t total = 0;
    uint8_t *dst = buf;

    while (total < size) {
        if (r->out_pos >= r->out_end) {
            int ret = zstd_reader_refill(r);
            if (ret <= 0) break;
        }

        size_t avail = r->out_end - r->out_pos;
        size_t need = size - total;
        size_t copy = (avail < need) ? avail : need;

        memcpy(dst + total, r->out_buf + r->out_pos, copy);
        r->out_pos += copy;
        total += copy;
    }

    return total;
}

static void zstd_reader_close(zstd_reader_t *r) {
    if (!r) return;
    ZSTD_freeDCtx(r->dctx);
    free(r->in_buf);
    free(r->out_buf);
    fclose(r->file);
    free(r);
}

// Little-endian read helpers for zstd reader
static uint16_t zstd_read_u16(zstd_reader_t *r) {
    uint8_t b[2];
    if (zstd_reader_read(r, b, 2) != 2) {
        LOG_ERROR("Unexpected end of file reading uint16 (zstd)");
        return 0;
    }
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t zstd_read_u32(zstd_reader_t *r) {
    uint8_t b[4];
    if (zstd_reader_read(r, b, 4) != 4) {
        LOG_ERROR("Unexpected end of file reading uint32 (zstd)");
        return 0;
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t zstd_read_i32(zstd_reader_t *r) {
    return (int32_t)zstd_read_u32(r);
}

static uint64_t zstd_read_u64(zstd_reader_t *r) {
    uint8_t b[8];
    if (zstd_reader_read(r, b, 8) != 8) {
        LOG_ERROR("Unexpected end of file reading uint64 (zstd)");
        return 0;
    }
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}

/**
 * Read reference database from zstd-compressed binary format (.trkb)
 * Populates the same global structures as readReferenceTree()
 *
 * @param filename Path to zstd-compressed .trkb file
 * @param name_specs Output array: [max_nodename, max_tax_name, max_line_taxonomy]
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinaryZstd(const char *filename, int *name_specs) {
    zstd_reader_t *r = zstd_reader_open(filename);
    if (!r) {
        LOG_ERROR("Cannot open zstd-compressed binary reference file: %s", filename);
        return -1;
    }

    // === Validate Header ===
    uint8_t magic[4];
    if (zstd_reader_read(r, magic, 4) != 4) {
        LOG_ERROR("Failed to read magic number (zstd)");
        zstd_reader_close(r);
        return -1;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        LOG_ERROR("Invalid binary format magic number (zstd)");
        zstd_reader_close(r);
        return -1;
    }

    // Read version
    uint8_t version_bytes[2];
    zstd_reader_read(r, version_bytes, 2);
    uint8_t version_major = version_bytes[0];
    uint8_t version_minor = version_bytes[1];
    LOG_DEBUG("Binary format version: %d.%d (zstd)", version_major, version_minor);

    if (version_major > 1) {
        LOG_ERROR("Unsupported binary format version: %d.%d", version_major, version_minor);
        zstd_reader_close(r);
        return -1;
    }

    // Read flags
    uint8_t flags[2];
    zstd_reader_read(r, flags, 2);
    uint8_t endianness = flags[0];
    uint8_t precision = flags[1];

    if (endianness != 0x01) {
        LOG_ERROR("Only little-endian binary format is supported");
        zstd_reader_close(r);
        return -1;
    }

    if (precision != 0x01) {
        LOG_ERROR("Only float precision binary format is supported");
        zstd_reader_close(r);
        return -1;
    }

    // Skip header CRC and reserved
    zstd_read_u32(r);  // header_crc
    zstd_read_u32(r);  // reserved

    // Read section offsets (for validation/logging only - we read sequentially)
    uint64_t taxonomy_offset = zstd_read_u64(r);
    uint64_t node_offset = zstd_read_u64(r);
    uint64_t posterior_offset = zstd_read_u64(r);
    uint64_t total_size = zstd_read_u64(r);

    LOG_DEBUG("Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu",
              (unsigned long)taxonomy_offset, (unsigned long)node_offset,
              (unsigned long)posterior_offset, (unsigned long)total_size);

    // Skip reserved bytes (16 bytes to complete 64-byte header)
    uint8_t skip[16];
    zstd_reader_read(r, skip, 16);

    // === Read Global Metadata (16 bytes) ===
    int32_t numberOfTrees = zstd_read_i32(r);
    int32_t max_nodename = zstd_read_i32(r);
    int32_t max_tax_name = zstd_read_i32(r);
    int32_t max_lineTaxonomy = zstd_read_i32(r);

    name_specs[0] = max_nodename;
    name_specs[1] = max_tax_name;
    name_specs[2] = max_lineTaxonomy;

    LOG_DEBUG("Database: %d trees, max_nodename=%d, max_tax=%d (zstd)",
              numberOfTrees, max_nodename, max_tax_name);

    // === Allocate Global Arrays ===
    numbaseArr = (int*)malloc(numberOfTrees * sizeof(int));
    rootArr = (int*)malloc(numberOfTrees * sizeof(int));
    numspecArr = (int*)malloc(numberOfTrees * sizeof(int));

    if (!numbaseArr || !rootArr || !numspecArr) {
        LOG_ERROR("Failed to allocate tree metadata arrays");
        zstd_reader_close(r);
        return -1;
    }

    // === Read Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < numberOfTrees; i++) {
        numbaseArr[i] = zstd_read_i32(r);
        rootArr[i] = zstd_read_i32(r);
        numspecArr[i] = zstd_read_i32(r);

        LOG_DEBUG("Tree %d: numbase=%d, root=%d, numspec=%d",
                  i, numbaseArr[i], rootArr[i], numspecArr[i]);
    }

    // === Read Taxonomy Section ===
    allocateMemoryForTaxArr(numberOfTrees, max_tax_name);

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t tree_tax_size = zstd_read_u32(r);
        zstd_read_u32(r);  // reserved
        (void)tree_tax_size;

        for (int s = 0; s < numspecArr[t]; s++) {
            for (int l = 0; l < 7; l++) {
                uint16_t len = zstd_read_u16(r);
                if (len > 0 && len <= (uint16_t)(max_tax_name + 1)) {
                    zstd_reader_read(r, taxonomyArr[t][s][l], len);
                } else if (len > 0) {
                    LOG_WARN("Taxonomy string length %d exceeds max %d", len, max_tax_name);
                    // Skip oversized string
                    char *skip_buf = malloc(len);
                    if (skip_buf) {
                        zstd_reader_read(r, skip_buf, len);
                        free(skip_buf);
                    }
                }
            }
        }

        LOG_DEBUG("Tree %d: read %d taxonomy entries (zstd)", t, numspecArr[t]);
    }

    // === Read Node Section ===
    treeArr = malloc(numberOfTrees * sizeof(struct node *));
    if (!treeArr) {
        LOG_ERROR("Failed to allocate tree array");
        zstd_reader_close(r);
        return -1;
    }

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t num_nodes = zstd_read_u32(r);
        int expected_nodes = 2 * numspecArr[t] - 1;

        if ((int)num_nodes != expected_nodes) {
            LOG_ERROR("Node count mismatch for tree %d: got %d, expected %d",
                      t, num_nodes, expected_nodes);
            zstd_reader_close(r);
            return -1;
        }

        allocateTreeArrMemory(t, max_nodename);

        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) {
            LOG_ERROR("Failed to allocate name offset array");
            zstd_reader_close(r);
            return -1;
        }

        // Read node records
        for (int n = 0; n < (int)num_nodes; n++) {
            treeArr[t][n].up[0] = zstd_read_i32(r);
            treeArr[t][n].up[1] = zstd_read_i32(r);
            treeArr[t][n].down = zstd_read_i32(r);
            treeArr[t][n].depth = zstd_read_i32(r);
            treeArr[t][n].taxIndex[0] = zstd_read_i32(r);
            treeArr[t][n].taxIndex[1] = zstd_read_i32(r);
            name_offsets[n] = zstd_read_u32(r);
            zstd_read_u32(r);  // reserved
        }

        // For zstd, we read names sequentially (can't seek)
        // Names are stored in order after node records
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                uint16_t len = zstd_read_u16(r);
                if (len > 0 && len <= (uint16_t)(max_nodename + 1)) {
                    zstd_reader_read(r, treeArr[t][n].name, len);
                }
            }
        }

        free(name_offsets);

        LOG_DEBUG("Tree %d: read %d node structures (zstd)", t, num_nodes);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t, (int)num_nodes, numbaseArr[t]);
            fflush(g_tsv_log_file);
        }
    }

    // === Read Posterior Section ===
    LOG_INFO("Loading posteriors from zstd-compressed binary format...");

    for (int t = 0; t < numberOfTrees; t++) {
        int num_nodes = 2 * numspecArr[t] - 1;
        int numbase = numbaseArr[t];

        for (int n = 0; n < num_nodes; n++) {
            size_t count = numbase * 4;

#ifdef OPTIMIZE_MEMORY
            if (zstd_reader_read(r, treeArr[t][n].posteriornc, count * sizeof(float)) != count * sizeof(float)) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (zstd)", t, n);
                zstd_reader_close(r);
                return -1;
            }
#else
            float *temp = malloc(count * sizeof(float));
            if (!temp) {
                LOG_ERROR("Failed to allocate temp buffer for posterior conversion");
                zstd_reader_close(r);
                return -1;
            }
            if (zstd_reader_read(r, temp, count * sizeof(float)) != count * sizeof(float)) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (zstd)", t, n);
                free(temp);
                zstd_reader_close(r);
                return -1;
            }
            for (size_t i = 0; i < count; i++) {
                treeArr[t][n].posteriornc[i] = (double)temp[i];
            }
            free(temp);
#endif
        }

        LOG_DEBUG("Tree %d: read posteriors for %d nodes (zstd)", t, num_nodes);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    t);
            fflush(g_tsv_log_file);
        }
    }

    zstd_reader_close(r);

    LOG_INFO("Loaded %d trees from zstd-compressed binary format", numberOfTrees);

    return numberOfTrees;
}
```

#### 5. tronko-assign.c Main Handler
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Add case for FORMAT_BINARY_ZSTD in format switch

Find the section where format is handled (search for `FORMAT_BINARY_GZIPPED`) and add:

```c
case FORMAT_BINARY_ZSTD:
    LOG_INFO("Loading zstd-compressed binary reference database...");
    numberOfPartitions = readReferenceBinaryZstd(opt.reference_file, name_specs);
    break;
```

### Success Criteria:

#### Automated Verification:
- [x] `make clean && make` succeeds in tronko-assign directory
- [x] Create test file: `cd tronko-convert && ./tronko-convert -i ../tronko-build/example_datasets/single_tree/reference_tree.txt -o /tmp/test.trkb`
- [x] `./tronko-assign -f /tmp/test.trkb -s -g ../example_datasets/16S_Bacteria/query_single.fasta -o /tmp/results.txt -w` completes successfully
- [x] Results match assignment with uncompressed file

#### Manual Verification:
- [x] Verbose output shows "zstd-compressed binary format" detection

---

## Phase 3: Round-Trip Testing

### Overview
Comprehensive testing of the full read/write cycle.

### Test Scripts:

#### 1. Basic Round-Trip
```bash
#!/bin/bash
set -e

cd tronko-convert

# Text -> Zstd binary
./tronko-convert -i ../tronko-build/example_datasets/single_tree/reference_tree.txt \
                 -o /tmp/test.trkb -v

# Verify it's zstd
file /tmp/test.trkb | grep -q "Zstandard"

# Zstd binary -> Text
./tronko-convert -i /tmp/test.trkb -o /tmp/roundtrip.txt -t -v

# Compare
diff ../tronko-build/example_datasets/single_tree/reference_tree.txt /tmp/roundtrip.txt

echo "Round-trip test PASSED"
```

#### 2. Assignment Comparison
```bash
#!/bin/bash
set -e

# Create both formats
cd tronko-convert
./tronko-convert -i ../example_datasets/16S_Bacteria/reference_tree.txt \
                 -o /tmp/ref.trkb
./tronko-convert -i ../example_datasets/16S_Bacteria/reference_tree.txt \
                 -o /tmp/ref_uncomp.trkb -u

cd ../tronko-assign

# Run assignment with zstd
./tronko-assign -f /tmp/ref.trkb -s \
    -g ../example_datasets/16S_Bacteria/query.fasta \
    -o /tmp/results_zstd.txt -w

# Run assignment with uncompressed
./tronko-assign -f /tmp/ref_uncomp.trkb -s \
    -g ../example_datasets/16S_Bacteria/query.fasta \
    -o /tmp/results_uncomp.txt -w

# Compare results
diff /tmp/results_zstd.txt /tmp/results_uncomp.txt

echo "Assignment comparison PASSED"
```

#### 3. Size Comparison
```bash
#!/bin/bash
cd tronko-convert

# Create all formats
./tronko-convert -i ../example_datasets/16S_Bacteria/reference_tree.txt \
                 -o /tmp/ref.trkb
./tronko-convert -i ../example_datasets/16S_Bacteria/reference_tree.txt \
                 -o /tmp/ref_uncomp.trkb -u

# Compress text with gzip for comparison
gzip -k -f ../example_datasets/16S_Bacteria/reference_tree.txt -c > /tmp/ref.txt.gz

echo "=== File Size Comparison ==="
ls -lh ../example_datasets/16S_Bacteria/reference_tree.txt /tmp/ref.txt.gz \
       /tmp/ref_uncomp.trkb /tmp/ref.trkb
```

### Success Criteria:

#### Automated Verification:
- [x] Round-trip test produces identical text files
- [x] Assignment results identical between zstd and uncompressed formats
- [x] Zstd file is smaller than gzipped text

#### Manual Verification:
- [x] File sizes match expected compression ratios (~3.9x vs gzipped text)
- [x] Loading time is acceptable (faster than gzip)

---

## Testing Strategy

### Unit Tests:
- Zstd writer produces valid zstd stream
- Zstd reader reads back what writer produces
- Format detection correctly identifies zstd-compressed TRKB
- Error handling for corrupted zstd data

### Integration Tests:
- Round-trip: text -> zstd binary -> text (exact match)
- Assignment results match between formats
- Large file handling (16S_Bacteria dataset)

### Manual Testing Steps:
1. Build both modules with `make clean && make`
2. Run round-trip test on example dataset
3. Run assignment test comparing zstd vs uncompressed results
4. Verify file sizes match expected compression ratios
5. Test error cases (invalid file, truncated file)

## Performance Considerations

- Streaming compression/decompression keeps memory usage bounded
- Zstd level 19 compression is slow but produces smallest files
- Decompression is very fast (~2.5x faster than gzip)
- Consider adding `-c <level>` flag for users who want faster compression (future)

## Migration Notes

- Existing `.trkb` files (uncompressed) will continue to work
- New default produces zstd-compressed files with same `.trkb` extension
- Use `-u` flag to produce uncompressed files for debugging/mmap future
- No version bump needed - format is auto-detected by magic bytes

## References

- Original research: `thoughts/shared/research/2026-01-02-zstd-trkb-compression-support.md`
- TRKB format spec: `tronko-convert/FORMAT_SPECIFICATION.md`
- Existing gzip reader: `tronko-assign/readreference.c:991-1313`
- Zstd streaming API: `tronko-assign/compressed_io.h`
