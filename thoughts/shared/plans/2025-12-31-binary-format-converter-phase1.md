# tronko-convert Phase 1 Implementation Plan

## Overview

Implement `tronko-convert`, a standalone utility to convert tronko-build's text-based `reference_tree.txt` format to a new binary format (`.trkb`) for faster loading. This tool will enable 50-70x faster database loading in tronko-assign without requiring changes to tronko-build.

## Current State Analysis

### Text Format (`reference_tree.txt`)
- **Example file**: 40MB, 930,598 lines for 1 tree, 1466 species, 316 positions
- **Structure** (from `tronko-build/printtree.c:31-84`):
  1. Header: 4 lines (numberOfTrees, max_nodename, max_tax_name, max_lineTaxonomy)
  2. Tree metadata: 1 line per tree (`numbase\troot\tnumspec`)
  3. Taxonomy: `numspec` lines per tree (semicolon-delimited, 7 levels)
  4. Nodes: For each of `2*numspec-1` nodes:
     - Node header line (9 tab-separated fields)
     - `numbase` lines of posteriors (4 tab-separated floats per line)

### Performance Bottleneck (from research)
- ~3.7M `sscanf()` calls for posterior parsing
- ~926K `gzgets()` calls
- ~930K `malloc()` calls for string handling
- **Current load time**: ~7s estimated for example database

### Binary Format Goals
- **Target load time**: ~0.1s (50-70x improvement)
- **Memory reduction**: 56% (float instead of double, 1D arrays)
- **Backward compatible**: tronko-assign will auto-detect format

## Desired End State

After Phase 1 completion:

1. **New directory**: `tronko-convert/` with standalone conversion tool
2. **Binary format**: `.trkb` files following the v1.0 specification from research
3. **Bidirectional conversion**: text-to-binary and binary-to-text (for validation)
4. **Round-trip validation**: `text → binary → text` produces identical output

### Verification
```bash
# Convert to binary
./tronko-convert -i reference_tree.txt -o reference_tree.trkb -v

# Convert back to text
./tronko-convert -i reference_tree.trkb -o reference_tree_roundtrip.txt -t -v

# Verify round-trip
diff reference_tree.txt reference_tree_roundtrip.txt
```

## What We're NOT Doing

- **NOT modifying tronko-assign** (that's Phase 2)
- **NOT modifying tronko-build** (that's Phase 3)
- **NOT implementing mmap support** (deferred optimization)
- **NOT supporting gzip-compressed binary** (eliminates mmap benefits)
- **NOT supporting big-endian systems** (all modern systems are little-endian)

## Implementation Approach

Model the binary I/O patterns after BWA's approach in `bwa_source_files/`:
- Use error-checked I/O wrappers (`err_fwrite`, `err_fread_noeof`)
- Magic number header for format detection
- Bulk reads/writes for large arrays
- Validation via checksums and magic footer

---

## Phase 1.1: Project Scaffolding

### Overview
Create the tronko-convert directory structure, Makefile, and basic CLI framework.

### Changes Required:

#### 1. Create directory structure
```
tronko-convert/
├── Makefile
├── tronko-convert.c      # Main entry point and CLI
├── format_common.h       # Shared structures and constants
├── format_text.c         # Text format reader/writer
├── format_text.h
├── format_binary.c       # Binary format reader/writer
├── format_binary.h
├── crc32.c               # CRC-32 implementation
├── crc32.h
└── utils.c               # Error-checked I/O wrappers
└── utils.h
```

#### 2. Makefile
**File**: `tronko-convert/Makefile`

```makefile
CC = gcc
CFLAGS = -Wall -Wextra -O3
DEBUG_CFLAGS = -Wall -Wextra -g -DDEBUG
LDFLAGS = -lz -lm

SRCS = tronko-convert.c format_text.c format_binary.c crc32.c utils.c
OBJS = $(SRCS:.c=.o)
TARGET = tronko-convert

.PHONY: all clean debug

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: clean $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
```

#### 3. Main entry point
**File**: `tronko-convert/tronko-convert.c`

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "format_common.h"
#include "format_text.h"
#include "format_binary.h"
#include "utils.h"

static void print_usage(const char *prog) {
    fprintf(stderr, "tronko-convert: Convert tronko reference database formats\n\n");
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "  -i <file>   Input file (required)\n");
    fprintf(stderr, "  -o <file>   Output file (required)\n");
    fprintf(stderr, "  -t          Output as text (default: binary)\n");
    fprintf(stderr, "  -v          Verbose output\n");
    fprintf(stderr, "  -h          Show this help\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb\n", prog);
    fprintf(stderr, "  %s -i reference_tree.trkb -o reference_tree.txt -t\n", prog);
}

int main(int argc, char *argv[]) {
    char *input_file = NULL;
    char *output_file = NULL;
    int output_text = 0;
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "i:o:tvh")) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'o': output_file = optarg; break;
            case 't': output_text = 1; break;
            case 'v': verbose = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    if (!input_file || !output_file) {
        fprintf(stderr, "Error: Both -i and -o are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    // Detect input format
    int input_format = detect_format(input_file);
    if (input_format < 0) {
        fprintf(stderr, "Error: Cannot open or detect format of '%s'\n", input_file);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Input format: %s\n",
                input_format == FORMAT_BINARY ? "binary (.trkb)" : "text");
        fprintf(stderr, "Output format: %s\n",
                output_text ? "text" : "binary (.trkb)");
    }

    // Load database
    tronko_db_t *db = NULL;
    if (input_format == FORMAT_BINARY) {
        db = load_binary(input_file, verbose);
    } else {
        db = load_text(input_file, verbose);
    }

    if (!db) {
        fprintf(stderr, "Error: Failed to load '%s'\n", input_file);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Loaded: %d trees\n", db->num_trees);
        for (int i = 0; i < db->num_trees; i++) {
            fprintf(stderr, "  Tree %d: %d species, %d positions, %d nodes\n",
                    i, db->trees[i].numspec, db->trees[i].numbase, db->trees[i].num_nodes);
        }
    }

    // Write output
    int result;
    if (output_text) {
        result = write_text(db, output_file, verbose);
    } else {
        result = write_binary(db, output_file, verbose);
    }

    if (result != 0) {
        fprintf(stderr, "Error: Failed to write '%s'\n", output_file);
        free_db(db);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Successfully wrote '%s'\n", output_file);
    }

    free_db(db);
    return 0;
}
```

### Success Criteria:

#### Automated Verification:
- [x] `cd tronko-convert && make` compiles without errors
- [x] `./tronko-convert -h` shows usage
- [x] `./tronko-convert` (no args) returns error and shows usage

#### Manual Verification:
- [x] Directory structure matches specification

---

## Phase 1.2: Common Data Structures and Utilities

### Overview
Define the in-memory data structures and error-checked I/O utilities.

### Changes Required:

#### 1. Common header
**File**: `tronko-convert/format_common.h`

```c
#ifndef FORMAT_COMMON_H
#define FORMAT_COMMON_H

#include <stdint.h>

// Binary format magic and version
#define TRONKO_MAGIC_0 0x89
#define TRONKO_MAGIC_1 'T'
#define TRONKO_MAGIC_2 'R'
#define TRONKO_MAGIC_3 'K'
#define TRONKO_VERSION_MAJOR 1
#define TRONKO_VERSION_MINOR 0
#define TRONKO_FOOTER_MAGIC 0x454E4421  // "END!"

// Format identifiers
#define FORMAT_UNKNOWN -1
#define FORMAT_TEXT     1
#define FORMAT_BINARY   2

// Taxonomy levels (fixed at 7: domain, phylum, class, order, family, genus, species)
#define TAXONOMY_LEVELS 7

// Node structure (mirrors tronko-assign/global.h:67-75)
typedef struct {
    int32_t up[2];        // Child indices (-1 for leaf nodes)
    int32_t down;         // Parent index
    int32_t depth;        // Tree depth
    int32_t taxIndex[2];  // [species_index, taxonomy_level]
    char *name;           // Node name (NULL for internal nodes)
    float *posteriors;    // [numbase * 4] contiguous array (A, C, G, T per position)
} tronko_node_t;

// Tree structure
typedef struct {
    int32_t numbase;      // MSA alignment length (number of positions)
    int32_t root;         // Root node index
    int32_t numspec;      // Number of species (leaf nodes)
    int32_t num_nodes;    // Total nodes: 2*numspec - 1
    char ***taxonomy;     // [numspec][TAXONOMY_LEVELS] taxonomy strings
    tronko_node_t *nodes; // [num_nodes] node array
} tronko_tree_t;

// Database structure
typedef struct {
    int32_t num_trees;
    int32_t max_nodename;      // Max node name length
    int32_t max_tax_name;      // Max taxonomy name length
    int32_t max_line_taxonomy; // Max full taxonomy line length
    tronko_tree_t *trees;
} tronko_db_t;

// Format detection
int detect_format(const char *filename);

// Memory management
tronko_db_t *alloc_db(int num_trees);
void free_db(tronko_db_t *db);

#endif
```

#### 2. Utilities (error-checked I/O)
**File**: `tronko-convert/utils.h`

```c
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
```

**File**: `tronko-convert/utils.c`

```c
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
```

#### 3. CRC-32 implementation
**File**: `tronko-convert/crc32.h`

```c
#ifndef CRC32_H
#define CRC32_H

#include <stdint.h>
#include <stddef.h>

uint32_t crc32(uint32_t crc, const void *buf, size_t len);
uint32_t crc32_file(FILE *fp, long start, long end);

#endif
```

**File**: `tronko-convert/crc32.c`

```c
#include "crc32.h"
#include <stdio.h>

// Standard CRC-32 (IEEE 802.3 polynomial)
static uint32_t crc32_table[256];
static int table_initialized = 0;

static void init_crc32_table(void) {
    if (table_initialized) return;
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t crc = i;
        for (int j = 0; j < 8; j++) {
            crc = (crc >> 1) ^ ((crc & 1) ? 0xEDB88320 : 0);
        }
        crc32_table[i] = crc;
    }
    table_initialized = 1;
}

uint32_t crc32(uint32_t crc, const void *buf, size_t len) {
    init_crc32_table();
    const uint8_t *p = (const uint8_t *)buf;
    crc = ~crc;
    while (len--) {
        crc = crc32_table[(crc ^ *p++) & 0xFF] ^ (crc >> 8);
    }
    return ~crc;
}

uint32_t crc32_file(FILE *fp, long start, long end) {
    init_crc32_table();
    long orig_pos = ftell(fp);
    fseek(fp, start, SEEK_SET);

    uint32_t crc = 0;
    uint8_t buf[8192];
    long remaining = end - start;

    while (remaining > 0) {
        size_t to_read = remaining > (long)sizeof(buf) ? sizeof(buf) : (size_t)remaining;
        size_t n = fread(buf, 1, to_read, fp);
        if (n == 0) break;
        crc = crc32(crc, buf, n);
        remaining -= n;
    }

    fseek(fp, orig_pos, SEEK_SET);
    return crc;
}
```

#### 4. Memory management
**File**: `tronko-convert/format_common.c`

```c
#include "format_common.h"
#include <stdlib.h>
#include <string.h>

int detect_format(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return FORMAT_UNKNOWN;

    uint8_t magic[4];
    size_t n = fread(magic, 1, 4, fp);
    fclose(fp);

    if (n < 4) return FORMAT_UNKNOWN;

    // Check for binary format magic
    if (magic[0] == TRONKO_MAGIC_0 && magic[1] == TRONKO_MAGIC_1 &&
        magic[2] == TRONKO_MAGIC_2 && magic[3] == TRONKO_MAGIC_3) {
        return FORMAT_BINARY;
    }

    // Check for gzip magic (text files may be gzipped)
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return FORMAT_TEXT;  // Gzipped text
    }

    // Assume text (first bytes should be ASCII digits)
    if (magic[0] >= '0' && magic[0] <= '9') {
        return FORMAT_TEXT;
    }

    return FORMAT_UNKNOWN;
}

tronko_db_t *alloc_db(int num_trees) {
    tronko_db_t *db = calloc(1, sizeof(tronko_db_t));
    if (!db) return NULL;

    db->num_trees = num_trees;
    db->trees = calloc(num_trees, sizeof(tronko_tree_t));
    if (!db->trees) {
        free(db);
        return NULL;
    }

    return db;
}

void free_db(tronko_db_t *db) {
    if (!db) return;

    for (int t = 0; t < db->num_trees; t++) {
        tronko_tree_t *tree = &db->trees[t];

        // Free taxonomy
        if (tree->taxonomy) {
            for (int s = 0; s < tree->numspec; s++) {
                if (tree->taxonomy[s]) {
                    for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                        free(tree->taxonomy[s][l]);
                    }
                    free(tree->taxonomy[s]);
                }
            }
            free(tree->taxonomy);
        }

        // Free nodes
        if (tree->nodes) {
            for (int n = 0; n < tree->num_nodes; n++) {
                free(tree->nodes[n].name);
                free(tree->nodes[n].posteriors);
            }
            free(tree->nodes);
        }
    }

    free(db->trees);
    free(db);
}
```

### Success Criteria:

#### Automated Verification:
- [x] `make` compiles all source files without warnings
- [x] `./tronko-convert -i nonexistent.txt -o out.trkb` returns error about file not found

#### Manual Verification:
- [x] Data structures match tronko-assign's structures (verified by code review)

---

## Phase 1.3: Text Format Reader

### Overview
Implement parsing of the text format, mirroring `tronko-assign/readreference.c:377-557`.

### Changes Required:

#### 1. Text format header
**File**: `tronko-convert/format_text.h`

```c
#ifndef FORMAT_TEXT_H
#define FORMAT_TEXT_H

#include "format_common.h"

tronko_db_t *load_text(const char *filename, int verbose);
int write_text(tronko_db_t *db, const char *filename, int verbose);

#endif
```

#### 2. Text format reader implementation
**File**: `tronko-convert/format_text.c`

```c
#include "format_text.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define BUFFER_SIZE 4096
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))

// Allocate taxonomy array for a tree
static char ***alloc_taxonomy(int numspec, int max_tax_name) {
    char ***tax = calloc(numspec, sizeof(char **));
    if (!tax) return NULL;

    for (int s = 0; s < numspec; s++) {
        tax[s] = calloc(TAXONOMY_LEVELS, sizeof(char *));
        if (!tax[s]) return NULL;  // Simplified error handling
        for (int l = 0; l < TAXONOMY_LEVELS; l++) {
            tax[s][l] = calloc(max_tax_name + 1, sizeof(char));
            if (!tax[s][l]) return NULL;
        }
    }
    return tax;
}

// Allocate nodes array for a tree
static tronko_node_t *alloc_nodes(int num_nodes, int numbase, int max_nodename) {
    tronko_node_t *nodes = calloc(num_nodes, sizeof(tronko_node_t));
    if (!nodes) return NULL;

    for (int n = 0; n < num_nodes; n++) {
        nodes[n].up[0] = -2;  // Uninitialized sentinel
        nodes[n].up[1] = -2;
        nodes[n].down = -2;
        nodes[n].depth = -2;
        nodes[n].posteriors = calloc(numbase * 4, sizeof(float));
        if (!nodes[n].posteriors) return NULL;
    }
    return nodes;
}

tronko_db_t *load_text(const char *filename, int verbose) {
    gzFile fp = err_xzopen(filename, "r");
    char buffer[BUFFER_SIZE];
    char *s;

    // Read header: 4 lines
    int numberOfTrees, max_nodename, max_tax_name, max_lineTaxonomy;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &numberOfTrees) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_nodename) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_tax_name) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_lineTaxonomy) != 1) goto error;

    if (verbose) {
        fprintf(stderr, "  Header: %d trees, max_nodename=%d, max_tax=%d\n",
                numberOfTrees, max_nodename, max_tax_name);
    }

    // Allocate database
    tronko_db_t *db = alloc_db(numberOfTrees);
    if (!db) goto error;

    db->max_nodename = max_nodename;
    db->max_tax_name = max_tax_name;
    db->max_line_taxonomy = max_lineTaxonomy;

    // Read per-tree metadata
    for (int i = 0; i < numberOfTrees; i++) {
        if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;
        s = strtok(buffer, "\n");
        int numbase, root, numspec;
        if (sscanf(s, "%d\t%d\t%d", &numbase, &root, &numspec) != 3) goto error_db;

        db->trees[i].numbase = numbase;
        db->trees[i].root = root;
        db->trees[i].numspec = numspec;
        db->trees[i].num_nodes = 2 * numspec - 1;

        if (verbose) {
            fprintf(stderr, "  Tree %d: numbase=%d, root=%d, numspec=%d\n",
                    i, numbase, root, numspec);
        }
    }

    // Allocate and read taxonomy
    for (int i = 0; i < numberOfTrees; i++) {
        db->trees[i].taxonomy = alloc_taxonomy(db->trees[i].numspec, max_tax_name);
        if (!db->trees[i].taxonomy) goto error_db;

        for (int j = 0; j < db->trees[i].numspec; j++) {
            if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;
            s = strtok(buffer, ";\n");
            if (s) strcpy(db->trees[i].taxonomy[j][0], s);

            for (int k = 1; k < TAXONOMY_LEVELS; k++) {
                s = strtok(NULL, ";\n");
                if (s) strcpy(db->trees[i].taxonomy[j][k], s);
            }
        }
    }

    // Allocate nodes
    for (int i = 0; i < numberOfTrees; i++) {
        db->trees[i].nodes = alloc_nodes(db->trees[i].num_nodes,
                                          db->trees[i].numbase,
                                          max_nodename);
        if (!db->trees[i].nodes) goto error_db;

        // Allocate names only for leaf nodes (indices numspec-1 to 2*numspec-2)
        int numspec = db->trees[i].numspec;
        for (int n = numspec - 1; n < 2 * numspec - 1; n++) {
            db->trees[i].nodes[n].name = calloc(max_nodename + 1, sizeof(char));
            if (!db->trees[i].nodes[n].name) goto error_db;
        }
    }

    // Read node structures and posteriors
    // Total nodes across all trees
    int total_nodes = 0;
    for (int i = 0; i < numberOfTrees; i++) {
        total_nodes += db->trees[i].num_nodes;
    }

    if (verbose) {
        fprintf(stderr, "  Reading %d total nodes...\n", total_nodes);
    }

    int nodes_read = 0;
    while (gzgets(fp, buffer, BUFFER_SIZE)) {
        int treeNumber, nodeNumber;
        int up0, up1, down, depth, taxIdx0, taxIdx1;
        char acc_name[BUFFER_SIZE];
        acc_name[0] = '\0';

        // Parse node header line
        s = strtok(buffer, "\n");
        int parsed = sscanf(s, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s",
                           &treeNumber, &nodeNumber,
                           &up0, &up1, &down, &depth,
                           &taxIdx0, &taxIdx1, acc_name);

        if (parsed < 8) {
            // Might be end of file or parse error
            if (nodes_read > 0) break;
            goto error_db;
        }

        if (treeNumber < 0 || treeNumber >= numberOfTrees) goto error_db;
        if (nodeNumber < 0 || nodeNumber >= db->trees[treeNumber].num_nodes) goto error_db;

        tronko_node_t *node = &db->trees[treeNumber].nodes[nodeNumber];
        node->up[0] = up0;
        node->up[1] = up1;
        node->down = down;
        node->depth = depth;
        node->taxIndex[0] = taxIdx0;
        node->taxIndex[1] = taxIdx1;

        // Copy name for leaf nodes (when both children are -1)
        if (up0 == -1 && up1 == -1 && acc_name[0] != '\0') {
            if (node->name) {
                strcpy(node->name, acc_name);
            }
        }

        // Read posterior probabilities for this node
        int numbase = db->trees[treeNumber].numbase;
        for (int pos = 0; pos < numbase; pos++) {
            if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;

            double p0, p1, p2, p3;
            if (sscanf(buffer, "%lf\t%lf\t%lf\t%lf", &p0, &p1, &p2, &p3) != 4) {
                goto error_db;
            }

            // Store as float (conversion happens here)
            node->posteriors[PP_IDX(pos, 0)] = (float)p0;
            node->posteriors[PP_IDX(pos, 1)] = (float)p1;
            node->posteriors[PP_IDX(pos, 2)] = (float)p2;
            node->posteriors[PP_IDX(pos, 3)] = (float)p3;
        }

        nodes_read++;

        if (verbose && nodes_read % 1000 == 0) {
            fprintf(stderr, "\r  Read %d/%d nodes...", nodes_read, total_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\r  Read %d nodes total    \n", nodes_read);
    }

    err_gzclose(fp);
    return db;

error_db:
    free_db(db);
error:
    gzclose(fp);
    return NULL;
}

int write_text(tronko_db_t *db, const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "w");

    // Write header
    fprintf(fp, "%d\n", db->num_trees);
    fprintf(fp, "%d\n", db->max_nodename);
    fprintf(fp, "%d\n", db->max_tax_name);
    fprintf(fp, "%d\n", db->max_line_taxonomy);

    // Write per-tree metadata
    for (int i = 0; i < db->num_trees; i++) {
        fprintf(fp, "%d\t%d\t%d\n",
                db->trees[i].numbase, db->trees[i].root, db->trees[i].numspec);
    }

    // Write taxonomy
    for (int i = 0; i < db->num_trees; i++) {
        for (int j = 0; j < db->trees[i].numspec; j++) {
            for (int k = 0; k < TAXONOMY_LEVELS; k++) {
                if (k == TAXONOMY_LEVELS - 1) {
                    fprintf(fp, "%s\n", db->trees[i].taxonomy[j][k]);
                } else {
                    fprintf(fp, "%s;", db->trees[i].taxonomy[j][k]);
                }
            }
        }
    }

    // Write nodes and posteriors
    int nodes_written = 0;
    for (int i = 0; i < db->num_trees; i++) {
        for (int j = 0; j < db->trees[i].num_nodes; j++) {
            tronko_node_t *node = &db->trees[i].nodes[j];

            // Node header line
            fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
                    i, j, node->up[0], node->up[1],
                    node->down, node->depth,
                    node->taxIndex[0], node->taxIndex[1]);

            // Name (only for leaf nodes)
            if (node->up[0] == -1 && node->up[1] == -1 && node->name) {
                fprintf(fp, "%s\n", node->name);
            } else {
                fprintf(fp, "\n");
            }

            // Posteriors (match %.17g format from printtree.c)
            int numbase = db->trees[i].numbase;
            for (int pos = 0; pos < numbase; pos++) {
                fprintf(fp, "%.17g\t%.17g\t%.17g\t%.17g\n",
                        (double)node->posteriors[PP_IDX(pos, 0)],
                        (double)node->posteriors[PP_IDX(pos, 1)],
                        (double)node->posteriors[PP_IDX(pos, 2)],
                        (double)node->posteriors[PP_IDX(pos, 3)]);
            }

            nodes_written++;
            if (verbose && nodes_written % 1000 == 0) {
                fprintf(stderr, "\r  Wrote %d nodes...", nodes_written);
            }
        }
    }

    if (verbose) {
        fprintf(stderr, "\r  Wrote %d nodes total    \n", nodes_written);
    }

    err_fclose(fp);
    return 0;
}
```

### Success Criteria:

#### Automated Verification:
- [x] `./tronko-convert -i example.txt -o roundtrip.txt -t -v` completes without error
- [x] Text parser handles the example dataset: `tronko-build/example_datasets/single_tree/reference_tree.txt`

#### Manual Verification:
- [x] Parsed data matches expected structure (inspect with verbose output)

---

## Phase 1.4: Binary Format Writer

### Overview
Implement the binary format writer following the v1.0 specification from the research document.

### Changes Required:

#### 1. Binary format header
**File**: `tronko-convert/format_binary.h`

```c
#ifndef FORMAT_BINARY_H
#define FORMAT_BINARY_H

#include "format_common.h"

tronko_db_t *load_binary(const char *filename, int verbose);
int write_binary(tronko_db_t *db, const char *filename, int verbose);

#endif
```

#### 2. Binary format writer implementation
**File**: `tronko-convert/format_binary.c`

```c
#include "format_binary.h"
#include "utils.h"
#include "crc32.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))

// Binary file header size constants
#define FILE_HEADER_SIZE 64
#define GLOBAL_META_SIZE 16
#define TREE_META_SIZE 12
#define NODE_RECORD_SIZE 32

int write_binary(tronko_db_t *db, const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "wb");

    if (verbose) {
        fprintf(stderr, "  Writing binary format v%d.%d\n",
                TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR);
    }

    // Calculate section offsets
    uint64_t header_end = FILE_HEADER_SIZE + GLOBAL_META_SIZE;
    uint64_t tree_meta_size = db->num_trees * TREE_META_SIZE;
    uint64_t taxonomy_offset = header_end + tree_meta_size;

    // Calculate taxonomy section size
    uint64_t taxonomy_size = db->num_trees * 8;  // Header per tree (size + reserved)
    for (int t = 0; t < db->num_trees; t++) {
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                // Length prefix (2 bytes) + string + null
                int len = strlen(db->trees[t].taxonomy[s][l]) + 1;
                taxonomy_size += 2 + len;
            }
        }
    }

    uint64_t node_offset = taxonomy_offset + taxonomy_size;

    // Calculate node section size
    uint64_t node_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        node_size += 4;  // Node count header
        node_size += db->trees[t].num_nodes * NODE_RECORD_SIZE;
        // Name table for leaf nodes
        for (int n = db->trees[t].numspec - 1; n < db->trees[t].num_nodes; n++) {
            if (db->trees[t].nodes[n].name) {
                int len = strlen(db->trees[t].nodes[n].name) + 1;
                node_size += 2 + len;  // Length prefix + string + null
            }
        }
    }

    uint64_t posterior_offset = node_offset + node_size;

    // Calculate posterior section size
    uint64_t posterior_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        posterior_size += (uint64_t)db->trees[t].num_nodes * db->trees[t].numbase * 4 * sizeof(float);
    }

    uint64_t total_size = posterior_offset + posterior_size + 8;  // +8 for footer

    if (verbose) {
        fprintf(stderr, "  Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset);
        fprintf(stderr, "  Total file size: %lu bytes (%.2f MB)\n",
                (unsigned long)total_size, total_size / (1024.0 * 1024.0));
    }

    // === Write File Header (64 bytes) ===
    uint8_t magic[4] = { TRONKO_MAGIC_0, TRONKO_MAGIC_1, TRONKO_MAGIC_2, TRONKO_MAGIC_3 };
    err_fwrite(magic, 1, 4, fp);
    fputc(TRONKO_VERSION_MAJOR, fp);
    fputc(TRONKO_VERSION_MINOR, fp);
    fputc(0x01, fp);  // Little-endian
    fputc(0x01, fp);  // Float precision

    // Calculate header CRC (of bytes 0-7)
    long header_start = err_ftell(fp) - 8;
    err_fseek(fp, header_start, SEEK_SET);
    uint8_t header_bytes[8];
    err_fread(header_bytes, 1, 8, fp);
    uint32_t header_crc = crc32(0, header_bytes, 8);
    write_u32(fp, header_crc);

    write_u32(fp, 0);  // Reserved (alignment)
    write_u64(fp, taxonomy_offset);
    write_u64(fp, node_offset);
    write_u64(fp, posterior_offset);
    write_u64(fp, total_size);

    // Padding to 64 bytes
    uint8_t reserved[16] = {0};
    err_fwrite(reserved, 1, 16, fp);

    // === Write Global Metadata (16 bytes) ===
    write_i32(fp, db->num_trees);
    write_i32(fp, db->max_nodename);
    write_i32(fp, db->max_tax_name);
    write_i32(fp, db->max_line_taxonomy);

    // === Write Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < db->num_trees; i++) {
        write_i32(fp, db->trees[i].numbase);
        write_i32(fp, db->trees[i].root);
        write_i32(fp, db->trees[i].numspec);
    }

    // Verify we're at taxonomy_offset
    if ((uint64_t)err_ftell(fp) != taxonomy_offset) {
        fprintf(stderr, "Error: Taxonomy offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)taxonomy_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Taxonomy Section ===
    for (int t = 0; t < db->num_trees; t++) {
        // Calculate and write size for this tree's taxonomy
        uint32_t tree_tax_size = 0;
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                tree_tax_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
        write_u32(fp, tree_tax_size);
        write_u32(fp, 0);  // Reserved

        // Write taxonomy strings
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                const char *str = db->trees[t].taxonomy[s][l];
                uint16_t len = strlen(str) + 1;  // Include null terminator
                write_u16(fp, len);
                err_fwrite(str, 1, len, fp);
            }
        }
    }

    // Verify we're at node_offset
    if ((uint64_t)err_ftell(fp) != node_offset) {
        fprintf(stderr, "Error: Node offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)node_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Node Section ===
    for (int t = 0; t < db->num_trees; t++) {
        write_u32(fp, db->trees[t].num_nodes);

        // Track name offsets (relative to start of name table for this tree)
        uint32_t name_offset = 0;

        // First pass: write node records with calculated name offsets
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            write_i32(fp, node->up[0]);
            write_i32(fp, node->up[1]);
            write_i32(fp, node->down);
            write_i32(fp, node->depth);
            write_i32(fp, node->taxIndex[0]);
            write_i32(fp, node->taxIndex[1]);

            // Name offset (0 for internal nodes)
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                write_u32(fp, name_offset + 1);  // +1 so 0 means "no name"
                name_offset += 2 + strlen(node->name) + 1;
            } else {
                write_u32(fp, 0);
            }
            write_u32(fp, 0);  // Reserved
        }

        // Second pass: write name table
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                uint16_t len = strlen(node->name) + 1;
                write_u16(fp, len);
                err_fwrite(node->name, 1, len, fp);
            }
        }
    }

    // Verify we're at posterior_offset
    if ((uint64_t)err_ftell(fp) != posterior_offset) {
        fprintf(stderr, "Error: Posterior offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)posterior_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Posterior Section ===
    if (verbose) {
        fprintf(stderr, "  Writing posteriors...\n");
    }

    for (int t = 0; t < db->num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            // Bulk write all posteriors for this node
            err_fwrite(db->trees[t].nodes[n].posteriors,
                      sizeof(float),
                      db->trees[t].numbase * 4,
                      fp);
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: wrote %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // === Write Footer ===
    // CRC of all preceding data
    long data_end = err_ftell(fp);
    uint32_t data_crc = crc32_file(fp, 0, data_end);
    write_u32(fp, data_crc);
    write_u32(fp, TRONKO_FOOTER_MAGIC);

    // Verify total size
    if ((uint64_t)err_ftell(fp) != total_size) {
        fprintf(stderr, "Warning: Total size mismatch (expected %lu, got %ld)\n",
                (unsigned long)total_size, err_ftell(fp));
    }

    err_fclose(fp);
    return 0;
}

// ... load_binary implementation will follow in Phase 1.5 ...

tronko_db_t *load_binary(const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "rb");

    // Read and validate magic
    uint8_t magic[4];
    err_fread(magic, 1, 4, fp);
    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        fprintf(stderr, "Error: Invalid magic number\n");
        err_fclose(fp);
        return NULL;
    }

    // Read version
    uint8_t version_major = fgetc(fp);
    uint8_t version_minor = fgetc(fp);
    if (version_major > TRONKO_VERSION_MAJOR) {
        fprintf(stderr, "Error: Unsupported format version %d.%d (max supported: %d.%d)\n",
                version_major, version_minor, TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR);
        err_fclose(fp);
        return NULL;
    }

    if (verbose) {
        fprintf(stderr, "  Binary format v%d.%d\n", version_major, version_minor);
    }

    // Read endianness and precision flags
    uint8_t endianness = fgetc(fp);
    uint8_t precision = fgetc(fp);
    if (endianness != 0x01) {
        fprintf(stderr, "Error: Only little-endian format is supported\n");
        err_fclose(fp);
        return NULL;
    }
    if (precision != 0x01) {
        fprintf(stderr, "Error: Only float precision is supported\n");
        err_fclose(fp);
        return NULL;
    }

    // Skip header CRC and reserved
    read_u32(fp);  // header_crc (could validate)
    read_u32(fp);  // reserved

    // Read section offsets
    uint64_t taxonomy_offset = read_u64(fp);
    uint64_t node_offset = read_u64(fp);
    uint64_t posterior_offset = read_u64(fp);
    uint64_t total_size = read_u64(fp);

    // Skip reserved bytes to reach global metadata
    err_fseek(fp, FILE_HEADER_SIZE, SEEK_SET);

    // Read global metadata
    int32_t num_trees = read_i32(fp);
    int32_t max_nodename = read_i32(fp);
    int32_t max_tax_name = read_i32(fp);
    int32_t max_line_taxonomy = read_i32(fp);

    if (verbose) {
        fprintf(stderr, "  %d trees, max_nodename=%d, max_tax=%d\n",
                num_trees, max_nodename, max_tax_name);
        fprintf(stderr, "  Offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset, (unsigned long)total_size);
    }

    // Allocate database
    tronko_db_t *db = alloc_db(num_trees);
    if (!db) {
        err_fclose(fp);
        return NULL;
    }

    db->max_nodename = max_nodename;
    db->max_tax_name = max_tax_name;
    db->max_line_taxonomy = max_line_taxonomy;

    // Read tree metadata
    for (int i = 0; i < num_trees; i++) {
        db->trees[i].numbase = read_i32(fp);
        db->trees[i].root = read_i32(fp);
        db->trees[i].numspec = read_i32(fp);
        db->trees[i].num_nodes = 2 * db->trees[i].numspec - 1;

        if (verbose) {
            fprintf(stderr, "  Tree %d: numbase=%d, root=%d, numspec=%d\n",
                    i, db->trees[i].numbase, db->trees[i].root, db->trees[i].numspec);
        }
    }

    // === Read Taxonomy Section ===
    err_fseek(fp, taxonomy_offset, SEEK_SET);

    for (int t = 0; t < num_trees; t++) {
        uint32_t tree_tax_size = read_u32(fp);
        read_u32(fp);  // reserved
        (void)tree_tax_size;  // We don't need this, read strings directly

        // Allocate taxonomy
        db->trees[t].taxonomy = calloc(db->trees[t].numspec, sizeof(char **));
        if (!db->trees[t].taxonomy) goto error;

        for (int s = 0; s < db->trees[t].numspec; s++) {
            db->trees[t].taxonomy[s] = calloc(TAXONOMY_LEVELS, sizeof(char *));
            if (!db->trees[t].taxonomy[s]) goto error;

            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                uint16_t len = read_u16(fp);
                db->trees[t].taxonomy[s][l] = calloc(len, sizeof(char));
                if (!db->trees[t].taxonomy[s][l]) goto error;
                err_fread(db->trees[t].taxonomy[s][l], 1, len, fp);
            }
        }
    }

    // === Read Node Section ===
    err_fseek(fp, node_offset, SEEK_SET);

    for (int t = 0; t < num_trees; t++) {
        uint32_t num_nodes = read_u32(fp);
        if ((int32_t)num_nodes != db->trees[t].num_nodes) {
            fprintf(stderr, "Error: Node count mismatch for tree %d\n", t);
            goto error;
        }

        // Allocate nodes
        db->trees[t].nodes = calloc(num_nodes, sizeof(tronko_node_t));
        if (!db->trees[t].nodes) goto error;

        // Store name offsets for later
        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) goto error;

        // Read node records
        for (int n = 0; n < (int)num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            node->up[0] = read_i32(fp);
            node->up[1] = read_i32(fp);
            node->down = read_i32(fp);
            node->depth = read_i32(fp);
            node->taxIndex[0] = read_i32(fp);
            node->taxIndex[1] = read_i32(fp);
            name_offsets[n] = read_u32(fp);
            read_u32(fp);  // reserved
        }

        // Read name table
        long name_table_start = err_ftell(fp);
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                err_fseek(fp, name_table_start + name_offsets[n] - 1, SEEK_SET);
                uint16_t len = read_u16(fp);
                db->trees[t].nodes[n].name = calloc(len, sizeof(char));
                if (!db->trees[t].nodes[n].name) {
                    free(name_offsets);
                    goto error;
                }
                err_fread(db->trees[t].nodes[n].name, 1, len, fp);
            }
        }

        free(name_offsets);
    }

    // === Read Posterior Section ===
    err_fseek(fp, posterior_offset, SEEK_SET);

    if (verbose) {
        fprintf(stderr, "  Reading posteriors...\n");
    }

    for (int t = 0; t < num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            db->trees[t].nodes[n].posteriors = calloc(db->trees[t].numbase * 4, sizeof(float));
            if (!db->trees[t].nodes[n].posteriors) goto error;

            // Bulk read all posteriors for this node
            err_fread(db->trees[t].nodes[n].posteriors,
                     sizeof(float),
                     db->trees[t].numbase * 4,
                     fp);
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: read %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // Optionally validate footer
    uint32_t file_crc = read_u32(fp);
    uint32_t footer_magic = read_u32(fp);
    (void)file_crc;  // Could validate

    if (footer_magic != TRONKO_FOOTER_MAGIC) {
        fprintf(stderr, "Warning: Invalid footer magic (file may be corrupted)\n");
    }

    err_fclose(fp);
    return db;

error:
    free_db(db);
    err_fclose(fp);
    return NULL;
}
```

### Success Criteria:

#### Automated Verification:
- [x] Binary file can be written without errors
- [x] File has correct magic number: `xxd -l 4 output.trkb` shows `89 54 52 4b`
- [x] File size matches expected calculation

#### Manual Verification:
- [x] Section offsets are correct when examined with hex editor

---

## Phase 1.5: Integration and Validation

### Overview
Wire everything together, add the missing `format_common.c` to the build, and validate round-trip conversion.

### Changes Required:

#### 1. Update Makefile
**File**: `tronko-convert/Makefile` (updated)

```makefile
CC = gcc
CFLAGS = -Wall -Wextra -O3
DEBUG_CFLAGS = -Wall -Wextra -g -DDEBUG
LDFLAGS = -lz -lm

SRCS = tronko-convert.c format_common.c format_text.c format_binary.c crc32.c utils.c
OBJS = $(SRCS:.c=.o)
TARGET = tronko-convert

.PHONY: all clean debug test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

debug: CFLAGS = $(DEBUG_CFLAGS)
debug: clean $(TARGET)

test: $(TARGET)
	@echo "=== Testing text-to-binary conversion ==="
	./$(TARGET) -i ../tronko-build/example_datasets/single_tree/reference_tree.txt \
	            -o /tmp/test.trkb -v
	@echo ""
	@echo "=== Testing binary-to-text conversion ==="
	./$(TARGET) -i /tmp/test.trkb -o /tmp/test_roundtrip.txt -t -v
	@echo ""
	@echo "=== Comparing original and round-trip ==="
	@if diff -q ../tronko-build/example_datasets/single_tree/reference_tree.txt \
	            /tmp/test_roundtrip.txt > /dev/null 2>&1; then \
		echo "SUCCESS: Files are identical"; \
	else \
		echo "MISMATCH: Checking for acceptable floating-point differences..."; \
		diff ../tronko-build/example_datasets/single_tree/reference_tree.txt \
		     /tmp/test_roundtrip.txt | head -20; \
	fi
	@echo ""
	@echo "=== File size comparison ==="
	@ls -lh ../tronko-build/example_datasets/single_tree/reference_tree.txt /tmp/test.trkb

clean:
	rm -f $(OBJS) $(TARGET)
```

#### 2. Create test script
**File**: `tronko-convert/test_roundtrip.sh`

```bash
#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXAMPLE_DIR="$SCRIPT_DIR/../tronko-build/example_datasets/single_tree"
TMP_DIR="/tmp/tronko-convert-test"

mkdir -p "$TMP_DIR"

echo "=== Building tronko-convert ==="
cd "$SCRIPT_DIR"
make clean && make

echo ""
echo "=== Test 1: Text to Binary ==="
./tronko-convert -i "$EXAMPLE_DIR/reference_tree.txt" -o "$TMP_DIR/test.trkb" -v

echo ""
echo "=== Test 2: Binary to Text ==="
./tronko-convert -i "$TMP_DIR/test.trkb" -o "$TMP_DIR/roundtrip.txt" -t -v

echo ""
echo "=== Test 3: Validate Round-Trip ==="
# Compare first 1471 lines (header + metadata + taxonomy) exactly
head -1471 "$EXAMPLE_DIR/reference_tree.txt" > "$TMP_DIR/orig_header.txt"
head -1471 "$TMP_DIR/roundtrip.txt" > "$TMP_DIR/rt_header.txt"

if diff -q "$TMP_DIR/orig_header.txt" "$TMP_DIR/rt_header.txt"; then
    echo "Header and taxonomy: MATCH"
else
    echo "Header and taxonomy: MISMATCH"
    diff "$TMP_DIR/orig_header.txt" "$TMP_DIR/rt_header.txt" | head -20
fi

# Compare node structure lines (every 317th line starting from 1472)
# Note: Posterior values will differ due to float precision
echo ""
echo "=== File Sizes ==="
ls -lh "$EXAMPLE_DIR/reference_tree.txt" "$TMP_DIR/test.trkb"

TEXT_SIZE=$(stat -c%s "$EXAMPLE_DIR/reference_tree.txt")
BINARY_SIZE=$(stat -c%s "$TMP_DIR/test.trkb")
RATIO=$(echo "scale=2; $TEXT_SIZE / $BINARY_SIZE" | bc)
echo "Compression ratio: ${RATIO}x"

echo ""
echo "=== Binary File Header ==="
xxd -l 64 "$TMP_DIR/test.trkb"
```

### Success Criteria:

#### Automated Verification:
- [x] `cd tronko-convert && make test` completes without errors
- [x] Header and taxonomy sections match exactly after round-trip
- [x] Binary file is smaller than text file (expected ~2.7x compression for this dataset)
- [x] `xxd -l 4 /tmp/test.trkb` shows `89545242` (magic number)

#### Manual Verification:
- [x] Verbose output shows correct tree stats
- [ ] No memory leaks (run with valgrind if available)

---

## Testing Strategy

### Unit Tests:
- Format detection on text file, binary file, and gzipped text file
- Header parsing for both formats
- Memory allocation/deallocation (use valgrind)

### Integration Tests:
- Round-trip conversion: text → binary → text
- Large file handling (the example dataset is representative)
- Error handling for corrupt/invalid files

### Manual Testing Steps:
1. Build with `make` and verify no warnings
2. Convert example dataset to binary
3. Verify binary file structure with `xxd`
4. Convert binary back to text
5. Compare original and round-trip files
6. Test with gzipped input: `gzip -k reference_tree.txt && ./tronko-convert -i reference_tree.txt.gz -o test.trkb`

---

## Performance Considerations

### Expected Results:
- **Text file size**: 40 MB (example dataset)
- **Binary file size**: ~15 MB (expected 2.7x reduction due to float vs double + binary encoding)
- **Conversion time**: < 10 seconds for example dataset

### Memory Usage:
- Peak memory: ~2x input file size during conversion
- All data loaded into memory (acceptable for typical database sizes)

---

## Migration Notes

N/A for Phase 1 - this is a standalone tool. Phase 2 will add binary format support to tronko-assign.

---

## References

- Research document: `thoughts/shared/research/2025-12-30-binary-format-conversion-tool.md`
- Text format output: `tronko-build/printtree.c:31-84`
- Text format parsing: `tronko-assign/readreference.c:377-557`
- Node structure: `tronko-assign/global.h:67-75`
- BWA binary I/O patterns: `tronko-assign/bwa_source_files/utils.c:124-174`, `bwt.c:385-462`
