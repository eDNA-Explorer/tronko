---
date: 2025-12-30T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: eDNA-Explorer/tronko
topic: "Binary format conversion tool for tronko reference database"
tags: [research, codebase, tronko-assign, tronko-build, binary-format, optimization, performance]
status: complete
last_updated: 2025-12-30
last_updated_by: Claude
---

# Research: Binary Format Conversion Tool for Tronko Reference Database

**Date**: 2025-12-30
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: eDNA-Explorer/tronko

## Research Question

How could we introduce a conversion tool to convert tronko-build's output into a binary format for faster loading? This follows up on optimization opportunities identified in the reference database loading research.

## Summary

The current text-based `reference_tree.txt` format requires ~3.7 million `sscanf` calls for a typical database (1 tree, 1466 species, 316 positions), causing significant loading overhead. A binary format could achieve **10-50x faster loading** by:

1. Eliminating text parsing overhead
2. Enabling direct memory reads (`fread`) or memory-mapping (`mmap`)
3. Supporting bulk allocation patterns

This document specifies a portable binary format (`.trkb`) and details two implementation approaches: a standalone conversion tool (`tronko-convert`) and native binary output in tronko-build.

## Binary Format Specification (v1.0)

### Design Principles

1. **Magic number identification** - Non-ASCII start byte + readable signature
2. **Little-endian throughout** - Modern standard, Intel/AMD compatible
3. **Explicit serialization** - No struct slurping for portability
4. **Section-based layout** - Offset pointers enable mmap compatibility
5. **Version-aware** - Forward compatibility for future enhancements
6. **Float precision** - Use `float` (4 bytes) instead of `double` (8 bytes) for posteriors

### File Structure Overview

```
+---------------------------+
| File Header (64 bytes)    |
+---------------------------+
| Global Metadata (16 bytes)|
+---------------------------+
| Tree Metadata Section     |
+---------------------------+
| Taxonomy String Table     |
+---------------------------+
| Node Structure Section    |
+---------------------------+
| Posterior Data Section    |
+---------------------------+
| [Optional] Checksum       |
+---------------------------+
```

### Section 1: File Header (64 bytes, fixed)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint8[4] | Magic: `0x89 'T' 'R' 'K'` |
| 4 | 1 | uint8 | Version major (1) |
| 5 | 1 | uint8 | Version minor (0) |
| 6 | 1 | uint8 | Endianness: 0x01=little, 0x02=big |
| 7 | 1 | uint8 | Precision flag: 0x01=float, 0x02=double |
| 8 | 4 | uint32 | Header CRC-32 (of bytes 0-7) |
| 12 | 4 | uint32 | Reserved (alignment) |
| 16 | 8 | uint64 | Taxonomy section offset |
| 24 | 8 | uint64 | Node section offset |
| 32 | 8 | uint64 | Posterior section offset |
| 40 | 8 | uint64 | Total file size |
| 48 | 16 | - | Reserved for future use |

### Section 2: Global Metadata (16 bytes, fixed)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | Number of trees |
| 4 | 4 | int32 | Max node name length |
| 8 | 4 | int32 | Max taxonomy name length |
| 12 | 4 | int32 | Max taxonomy line length |

### Section 3: Tree Metadata (12 bytes per tree)

For each tree `i` in `[0, numberOfTrees)`:

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | `numbaseArr[i]` - MSA alignment length |
| 4 | 4 | int32 | `rootArr[i]` - Root node index |
| 8 | 4 | int32 | `numspecArr[i]` - Number of species |

### Section 4: Taxonomy String Table

**Header** (8 bytes per tree):
| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint32 | Taxonomy data size for this tree |
| 4 | 4 | uint32 | Reserved |

**Data** (per tree, per species):
```
For each species in tree:
  For each of 7 taxonomy levels:
    uint16 length (including null terminator)
    char[length] taxonomy string (null-terminated)
```

This length-prefixed format allows:
- Fast seeking to any taxonomy entry
- No fixed buffer size assumptions
- Efficient memory allocation

### Section 5: Node Structure Section

**Per-tree header** (4 bytes):
| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint32 | Number of nodes (should be `2*numspec-1`) |

**Per-node record** (32 bytes, fixed):
| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | `up[0]` - First child index |
| 4 | 4 | int32 | `up[1]` - Second child index |
| 8 | 4 | int32 | `down` - Parent index |
| 12 | 4 | int32 | `depth` - Tree depth |
| 16 | 4 | int32 | `taxIndex[0]` - Taxonomy row |
| 20 | 4 | int32 | `taxIndex[1]` - Taxonomy level |
| 24 | 4 | uint32 | Name offset (0 if internal node) |
| 28 | 4 | uint32 | Reserved |

**Node name table** (after all node records):
```
For each leaf node:
  uint16 length
  char[length] name (null-terminated)
```

### Section 6: Posterior Data Section

**Per-tree posterior block**:
```
For each tree:
  For each node (2*numspec-1):
    For each position (numbase):
      float[4] posteriors (A, C, G, T)
```

Total size per tree: `(2*numspec-1) * numbase * 4 * sizeof(float)` bytes

For the example database (1466 species, 316 positions):
- Nodes: 2931
- Posteriors: 2931 × 316 × 4 × 4 = **14.84 MB** (vs 28.26 MB with double)

### Section 7: Optional Checksum (8 bytes)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint32 | CRC-32 of all preceding data |
| 4 | 4 | uint32 | Magic footer: `0x454E4421` ("END!") |

## Implementation Approaches

### Approach A: Standalone Conversion Tool (Recommended for Phase 1)

Create `tronko-convert` as a separate utility:

```bash
# Convert text to binary
tronko-convert -i reference_tree.txt -o reference_tree.trkb

# Convert binary to text (for debugging)
tronko-convert -i reference_tree.trkb -o reference_tree.txt -t
```

**Advantages:**
- No changes to tronko-build required
- Can be added incrementally
- Easier to test and validate
- Users can choose format

**File structure:**

```
tronko-convert/
├── Makefile
├── tronko-convert.c      # Main entry point
├── format_text.c         # Text format reader/writer
├── format_binary.c       # Binary format reader/writer
├── format_common.h       # Shared structures
└── crc32.c               # CRC-32 implementation
```

### Approach B: Native Binary Output in tronko-build

Add `-b` flag to tronko-build:

```bash
tronko-build -l -m msa.fasta -t tree.nwk -x taxonomy.txt -d output/ -b
# Outputs: output/reference_tree.trkb
```

**Implementation:**
1. Add binary output option to `tronko-build/options.c`
2. Create `tronko-build/printbinary.c` alongside `printtree.c`
3. Call `printBinaryFile()` when `-b` flag is set

### Approach C: Hybrid (tronko-assign auto-detection)

Modify tronko-assign to auto-detect format:

```c
// In tronko-assign.c
int format = detect_format(opt.reference_file);
if (format == FORMAT_BINARY) {
    numberOfTrees = readReferenceBinary(filename, name_specs);
} else {
    numberOfTrees = readReferenceTree(referenceTree, name_specs);
}
```

## Implementation Plan

### Phase 1: Standalone Converter (2-3 days)

**Step 1.1: Create tronko-convert scaffolding**

```c
// tronko-convert/tronko-convert.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include "format_common.h"

void print_usage(void) {
    fprintf(stderr, "Usage: tronko-convert [options]\n");
    fprintf(stderr, "  -i <file>   Input file (required)\n");
    fprintf(stderr, "  -o <file>   Output file (required)\n");
    fprintf(stderr, "  -t          Output as text (default: binary)\n");
    fprintf(stderr, "  -v          Verbose output\n");
    fprintf(stderr, "  -h          Show this help\n");
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
            case 'h': print_usage(); return 0;
            default: print_usage(); return 1;
        }
    }

    if (!input_file || !output_file) {
        print_usage();
        return 1;
    }

    // Detect input format
    int input_format = detect_format(input_file);

    // Load data
    tronko_db_t *db = NULL;
    if (input_format == FORMAT_BINARY) {
        db = load_binary(input_file);
    } else {
        db = load_text(input_file);
    }

    if (!db) {
        fprintf(stderr, "Error: Failed to load %s\n", input_file);
        return 1;
    }

    // Write output
    int result;
    if (output_text) {
        result = write_text(db, output_file);
    } else {
        result = write_binary(db, output_file);
    }

    free_db(db);
    return result;
}
```

**Step 1.2: Define common data structures**

```c
// tronko-convert/format_common.h
#ifndef FORMAT_COMMON_H
#define FORMAT_COMMON_H

#include <stdint.h>

#define TRONKO_MAGIC_0 0x89
#define TRONKO_MAGIC_1 'T'
#define TRONKO_MAGIC_2 'R'
#define TRONKO_MAGIC_3 'K'
#define TRONKO_VERSION_MAJOR 1
#define TRONKO_VERSION_MINOR 0

#define FORMAT_TEXT   1
#define FORMAT_BINARY 2

typedef struct {
    int32_t up[2];
    int32_t down;
    int32_t depth;
    int32_t taxIndex[2];
    char *name;           // NULL for internal nodes
    float *posteriors;    // [numbase * 4] contiguous array
} tronko_node_t;

typedef struct {
    int32_t numbase;
    int32_t root;
    int32_t numspec;
    int32_t num_nodes;    // 2*numspec - 1
    char ***taxonomy;     // [numspec][7]
    tronko_node_t *nodes; // [num_nodes]
} tronko_tree_t;

typedef struct {
    int32_t num_trees;
    int32_t max_nodename;
    int32_t max_tax_name;
    int32_t max_line_taxonomy;
    tronko_tree_t *trees;
} tronko_db_t;

// Format detection
int detect_format(const char *filename);

// Text format I/O
tronko_db_t *load_text(const char *filename);
int write_text(tronko_db_t *db, const char *filename);

// Binary format I/O
tronko_db_t *load_binary(const char *filename);
int write_binary(tronko_db_t *db, const char *filename);

// Memory management
void free_db(tronko_db_t *db);

#endif
```

**Step 1.3: Implement binary writer**

```c
// tronko-convert/format_binary.c (excerpt)
#include "format_common.h"
#include <stdio.h>
#include <string.h>

// Little-endian write helpers
static void write_u32(FILE *fp, uint32_t v) {
    uint8_t b[4] = { v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF, (v >> 24) & 0xFF };
    fwrite(b, 1, 4, fp);
}

static void write_u64(FILE *fp, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (v >> (i * 8)) & 0xFF;
    fwrite(b, 1, 8, fp);
}

static void write_float(FILE *fp, float v) {
    fwrite(&v, sizeof(float), 1, fp);  // Native float format is portable enough
}

int write_binary(tronko_db_t *db, const char *filename) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) return -1;

    // Calculate section offsets
    uint64_t header_size = 64 + 16;  // File header + global metadata
    uint64_t tree_meta_size = db->num_trees * 12;
    uint64_t taxonomy_offset = header_size + tree_meta_size;

    // ... calculate other offsets ...

    // Write file header (64 bytes)
    uint8_t magic[4] = { TRONKO_MAGIC_0, TRONKO_MAGIC_1, TRONKO_MAGIC_2, TRONKO_MAGIC_3 };
    fwrite(magic, 1, 4, fp);
    fputc(TRONKO_VERSION_MAJOR, fp);
    fputc(TRONKO_VERSION_MINOR, fp);
    fputc(0x01, fp);  // Little-endian
    fputc(0x01, fp);  // Float precision
    write_u32(fp, 0); // CRC placeholder
    write_u32(fp, 0); // Reserved
    write_u64(fp, taxonomy_offset);
    write_u64(fp, node_offset);
    write_u64(fp, posterior_offset);
    write_u64(fp, total_size);

    // Padding to 64 bytes
    uint8_t reserved[16] = {0};
    fwrite(reserved, 1, 16, fp);

    // Write global metadata (16 bytes)
    write_u32(fp, db->num_trees);
    write_u32(fp, db->max_nodename);
    write_u32(fp, db->max_tax_name);
    write_u32(fp, db->max_line_taxonomy);

    // Write tree metadata
    for (int i = 0; i < db->num_trees; i++) {
        write_u32(fp, db->trees[i].numbase);
        write_u32(fp, db->trees[i].root);
        write_u32(fp, db->trees[i].numspec);
    }

    // Write taxonomy section
    // ... length-prefixed strings ...

    // Write node structures
    // ... fixed 32-byte records ...

    // Write posterior data (bulk)
    for (int t = 0; t < db->num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            fwrite(db->trees[t].nodes[n].posteriors,
                   sizeof(float),
                   db->trees[t].numbase * 4,
                   fp);
        }
    }

    fclose(fp);
    return 0;
}
```

**Step 1.4: Implement binary reader**

```c
tronko_db_t *load_binary(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return NULL;

    // Validate magic
    uint8_t magic[4];
    fread(magic, 1, 4, fp);
    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        fclose(fp);
        return NULL;
    }

    // Read version
    uint8_t version_major = fgetc(fp);
    uint8_t version_minor = fgetc(fp);
    if (version_major > TRONKO_VERSION_MAJOR) {
        fprintf(stderr, "Error: Unsupported format version %d.%d\n",
                version_major, version_minor);
        fclose(fp);
        return NULL;
    }

    // Read rest of header and allocate structures...
    // ... (similar pattern to BWA's bwt_restore_bwt)

    return db;
}
```

### Phase 2: Update tronko-assign to support binary format

**Step 2.1: Add format detection**

```c
// tronko-assign/readreference.c (new function)
int detect_reference_format(const char *filename) {
    // Try to open uncompressed first
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        // Try gzipped
        gzFile gz = gzopen(filename, "rb");
        if (!gz) return -1;

        char buf[4];
        gzread(gz, buf, 4);
        gzclose(gz);

        // Check if it's binary (after gunzip)
        if (buf[0] == 0x89 && buf[1] == 'T' && buf[2] == 'R' && buf[3] == 'K') {
            return FORMAT_BINARY;  // Gzipped binary
        }
        return FORMAT_TEXT;  // Gzipped text
    }

    uint8_t magic[4];
    fread(magic, 1, 4, fp);
    fclose(fp);

    if (magic[0] == 0x89 && magic[1] == 'T' && magic[2] == 'R' && magic[3] == 'K') {
        return FORMAT_BINARY;
    }
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return FORMAT_TEXT;  // Gzipped text
    }
    return FORMAT_TEXT;  // Plain text
}
```

**Step 2.2: Add binary reading function**

```c
// tronko-assign/readreference_binary.c (new file)
int readReferenceBinary(const char *filename, int *name_specs) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return -1;

    // Skip header validation (already done in detect_format)
    fseek(fp, 64, SEEK_SET);  // Skip file header

    // Read global metadata
    int32_t buf[4];
    fread(buf, sizeof(int32_t), 4, fp);
    int numberOfTrees = buf[0];
    name_specs[0] = buf[1];  // max_nodename
    name_specs[1] = buf[2];  // max_tax_name
    name_specs[2] = buf[3];  // max_line_taxonomy

    // Allocate global arrays
    numbaseArr = malloc(numberOfTrees * sizeof(int));
    rootArr = malloc(numberOfTrees * sizeof(int));
    numspecArr = malloc(numberOfTrees * sizeof(int));

    // Read tree metadata
    for (int i = 0; i < numberOfTrees; i++) {
        fread(buf, sizeof(int32_t), 3, fp);
        numbaseArr[i] = buf[0];
        rootArr[i] = buf[1];
        numspecArr[i] = buf[2];
    }

    // ... continue reading sections using fread ...
    // Key optimization: read posterior data in bulk

    fclose(fp);
    return numberOfTrees;
}
```

**Step 2.3: Update main() to use auto-detection**

```c
// tronko-assign/tronko-assign.c (modified)
int format = detect_reference_format(opt.reference_file);
if (format == FORMAT_BINARY) {
    numberOfTrees = readReferenceBinary(opt.reference_file, name_specs);
} else {
    gzFile referenceTree = gzopen(opt.reference_file, "r");
    numberOfTrees = readReferenceTree(referenceTree, name_specs);
    gzclose(referenceTree);
}
```

### Phase 3: Native binary output in tronko-build (Optional)

Add `-b` flag to output binary directly, eliminating the conversion step for new databases.

## Performance Projections

### Current Text Format Loading

For example database (1 tree, 1466 species, 316 positions):

| Operation | Count | Estimated Time |
|-----------|-------|----------------|
| gzgets() calls | ~926,000 | ~2.0s |
| strtok() calls | ~3,800,000 | ~1.5s |
| sscanf() calls | ~3,700,000 | ~3.0s |
| malloc() calls | ~930,000 | ~0.5s |
| **Total** | | **~7.0s** |

### Binary Format Loading (Projected)

| Operation | Count | Estimated Time |
|-----------|-------|----------------|
| fread() for posteriors | 1 bulk read | ~0.05s |
| fread() for nodes | 1 bulk read | ~0.01s |
| fread() for taxonomy | ~10,000 | ~0.02s |
| malloc() calls | ~10 | ~0.001s |
| **Total** | | **~0.1s** |

**Projected speedup: 50-70x**

### Memory Comparison

| Component | Text (double) | Binary (float) | Savings |
|-----------|--------------|----------------|---------|
| Posteriors | 28.26 MB | 14.13 MB | 50% |
| Pointers | 7.07 MB | 0 MB (1D array) | 100% |
| Total | ~39 MB | ~17 MB | 56% |

## Code References

- `tronko-build/printtree.c:31-84` - Current text output implementation
- `tronko-assign/readreference.c:377-557` - Current text parsing implementation
- `tronko-assign/allocatetreememory.c:17-40` - Node memory allocation
- `tronko-assign/bwa_source_files/bwt.c:385-462` - BWA binary I/O patterns (model)
- `tronko-assign/bwa_source_files/utils.c:124-140` - Error-checked I/O wrappers
- `tronko-assign/global.h:60-68` - Node structure definition

## Architecture Insights

1. **BWA Pattern**: The BWA code provides an excellent model for binary file I/O in this codebase, using error-checked wrappers and bulk reads.

2. **Backward Compatibility**: By auto-detecting format, existing text databases continue to work while new binary databases provide performance benefits.

3. **Conversion Tool First**: Starting with a standalone converter reduces risk - users can validate the binary format produces identical results before modifying tronko-assign.

4. **Float Precision**: The existing `USE_FLOAT_PP` compile flag in tronko-assign shows float precision was already considered adequate for assignment.

## Related Research

- `thoughts/shared/research/2025-12-29-reference-database-loading.md` - Original optimization analysis

## Open Questions

1. **Compression**: Should binary format support gzip compression? (Reduces disk space but eliminates mmap benefits)

2. **Memory mapping**: Full mmap support requires page-aligned sections. Worth the complexity?

3. **Endianness**: Little-endian is standard, but should we support big-endian systems with runtime byte-swapping?

4. **Validation**: Add checksums per-section or just file-level?

## Recommendations

1. **Start with tronko-convert** - Low risk, validates format before modifying core tools
2. **Use float precision** - Already validated as adequate, cuts memory in half
3. **Little-endian only** - All modern x86/ARM systems, simplifies implementation
4. **Bulk posterior reads** - Single fread() for the entire posterior section per tree
5. **Skip mmap initially** - Add later if needed; fread() is simpler and nearly as fast
