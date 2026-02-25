# Gzipped Binary Format Support Implementation Plan

## Overview

Add support for loading gzipped binary reference trees (`.trkb.gz`) in tronko-assign. This enables 2.15x smaller file sizes compared to gzipped text while maintaining the performance benefits of the binary format (no string parsing, bulk reads).

## Current State Analysis

### Format Detection (`readreference.c:52-89`)
The `detect_reference_format()` function checks magic bytes to identify file format:
- Binary magic: `0x89 'T' 'R' 'K'` → returns `FORMAT_BINARY`
- Gzip magic: `0x1f 0x8b` → returns `FORMAT_TEXT` (incorrect for `.trkb.gz`)
- ASCII digit → returns `FORMAT_TEXT`

**Problem**: All gzipped files are assumed to be text format.

### Binary Reader (`readreference.c:649-929`)
Uses standard `FILE*` with `fopen()`, `fread()`, `fseek()`, `fclose()`. Cannot read gzipped content.

### Text Reader (`readreference.c:459-639`)
Uses zlib's `gzFile` with `gzopen()`, `gzgets()`, `gzclose()`. Already handles gzipped text.

### Dispatch Logic (`tronko-assign.c:900-927`)
```c
int ref_format = detect_reference_format(opt.reference_file);
if (ref_format == FORMAT_BINARY) {
    numberOfTrees = readReferenceBinary(opt.reference_file, name_specs);
} else if (ref_format == FORMAT_TEXT) {
    gzFile referenceTree = gzopen(opt.reference_file, "r");
    numberOfTrees = readReferenceTree(referenceTree, name_specs);
}
```

### Key Discoveries
- Binary format uses forward-only seeks to section offsets (compatible with `gzseek()`)
- Name table reading has local backward seeks that need restructuring
- All helper functions (`read_u16`, `read_u32`, etc.) use `FILE*`

## Desired End State

After implementation:
1. `tronko-assign -f reference_tree.trkb.gz ...` works correctly
2. Format detection correctly identifies gzipped binary files
3. Performance remains good (no excessive decompression overhead)

### Verification
```bash
# Create test gzipped binary file
gzip -k example_datasets/single_tree/reference_tree.trkb

# Test assignment with gzipped binary
tronko-assign -r -f example_datasets/single_tree/reference_tree.trkb.gz \
  -a example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_results_gz.txt -w

# Compare output with uncompressed binary
diff /tmp/test_results.txt /tmp/test_results_gz.txt
```

## What We're NOT Doing

- Adding zstd support (future enhancement)
- Memory-mapped decompression (optimization for later)
- Changing the binary format itself
- Modifying tronko-convert (already outputs uncompressed binary; users can gzip externally)

## Implementation Approach

Use zlib's `gzFile` for the binary reader, mirroring the text reader's approach. This is simpler than memory decompression and leverages existing zlib integration.

---

## Phase 1: Add FORMAT_BINARY_GZIPPED Constant

### Overview
Add a new format constant to distinguish gzipped binary from gzipped text.

### Changes Required

#### 1. Update global.h
**File**: `tronko-assign/global.h`
**Changes**: Add new format constant

```c
// Reference file format constants (must match tronko-convert/format_common.h)
#define FORMAT_UNKNOWN        -1
#define FORMAT_TEXT            1
#define FORMAT_BINARY          2
#define FORMAT_BINARY_GZIPPED  3  // NEW: Gzipped binary format
```

#### 2. Update format_common.h (for consistency)
**File**: `tronko-convert/format_common.h`
**Changes**: Add matching constant

```c
// Format identifiers
#define FORMAT_UNKNOWN        -1
#define FORMAT_TEXT            1
#define FORMAT_BINARY          2
#define FORMAT_BINARY_GZIPPED  3  // NEW: Gzipped binary format
```

### Success Criteria

#### Automated Verification:
- [x] Code compiles without warnings: `cd tronko-assign && make clean && make`
- [x] No conflicts with existing constants

#### Manual Verification:
- [x] Constants match between tronko-assign and tronko-convert

---

## Phase 2: Update Format Detection

### Overview
Modify `detect_reference_format()` to peek inside gzipped files and detect their actual content type.

### Changes Required

#### 1. Update detect_reference_format()
**File**: `tronko-assign/readreference.c`
**Changes**: Replace gzip detection logic (lines 75-79)

```c
int detect_reference_format(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        LOG_ERROR("Cannot open file for format detection: %s", filename);
        return FORMAT_UNKNOWN;
    }

    uint8_t magic[4];
    size_t n = fread(magic, 1, 4, fp);
    fclose(fp);

    if (n < 4) {
        LOG_ERROR("File too small for format detection: %s", filename);
        return FORMAT_UNKNOWN;
    }

    // Check for binary format magic: 0x89 'T' 'R' 'K'
    if (magic[0] == TRONKO_MAGIC_0 && magic[1] == TRONKO_MAGIC_1 &&
        magic[2] == TRONKO_MAGIC_2 && magic[3] == TRONKO_MAGIC_3) {
        LOG_DEBUG("Detected binary format: %s", filename);
        return FORMAT_BINARY;
    }

    // Check for gzip magic: 0x1f 0x8b
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        // Gzipped file - decompress first 4 bytes to check content type
        gzFile gz = gzopen(filename, "rb");
        if (gz) {
            uint8_t inner_magic[4];
            int bytes_read = gzread(gz, inner_magic, 4);
            gzclose(gz);

            if (bytes_read == 4) {
                if (inner_magic[0] == TRONKO_MAGIC_0 && inner_magic[1] == TRONKO_MAGIC_1 &&
                    inner_magic[2] == TRONKO_MAGIC_2 && inner_magic[3] == TRONKO_MAGIC_3) {
                    LOG_DEBUG("Detected gzipped binary format: %s", filename);
                    return FORMAT_BINARY_GZIPPED;
                }
            }
        }
        LOG_DEBUG("Detected gzipped text format: %s", filename);
        return FORMAT_TEXT;  // Gzipped text
    }

    // Assume text if first byte is ASCII digit (numberOfTrees line)
    if (magic[0] >= '0' && magic[0] <= '9') {
        LOG_DEBUG("Detected plain text format: %s", filename);
        return FORMAT_TEXT;
    }

    LOG_WARN("Unknown file format: %s", filename);
    return FORMAT_UNKNOWN;
}
```

### Success Criteria

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make`
- [x] Format detection test:
  ```bash
  # Should return different values for each format
  gzip -c example_datasets/single_tree/reference_tree.trkb > /tmp/test.trkb.gz
  # Verify detection works (will need a test harness or debug output)
  ```

#### Manual Verification:
- [x] `.trkb` files return `FORMAT_BINARY`
- [x] `.trkb.gz` files return `FORMAT_BINARY_GZIPPED`
- [x] `.txt.gz` files return `FORMAT_TEXT`
- [x] `.txt` files return `FORMAT_TEXT`

---

## Phase 3: Create Gzip-Aware Binary Reader Helpers

### Overview
Create gzip-compatible versions of the read helper functions using zlib's `gzFile`.

### Changes Required

#### 1. Add gzip read helpers
**File**: `tronko-assign/readreference.c`
**Changes**: Add after existing `read_*` helpers (around line 47)

```c
// Gzip-compatible little-endian read helpers for binary format
static uint16_t gz_read_u16(gzFile gz) {
    uint8_t b[2];
    if (gzread(gz, b, 2) != 2) {
        LOG_ERROR("Unexpected end of file reading uint16 (gzip)");
        return 0;
    }
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t gz_read_u32(gzFile gz) {
    uint8_t b[4];
    if (gzread(gz, b, 4) != 4) {
        LOG_ERROR("Unexpected end of file reading uint32 (gzip)");
        return 0;
    }
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t gz_read_i32(gzFile gz) {
    return (int32_t)gz_read_u32(gz);
}

static uint64_t gz_read_u64(gzFile gz) {
    uint8_t b[8];
    if (gzread(gz, b, 8) != 8) {
        LOG_ERROR("Unexpected end of file reading uint64 (gzip)");
        return 0;
    }
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}
```

### Success Criteria

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make`
- [x] No compiler warnings about unused functions (will be used in Phase 4)

---

## Phase 4: Implement readReferenceBinaryGzipped()

### Overview
Create a gzip-aware version of `readReferenceBinary()` that uses `gzFile` instead of `FILE*`.

### Changes Required

#### 1. Add new function
**File**: `tronko-assign/readreference.c`
**Changes**: Add after `readReferenceBinary()` (around line 930)

```c
/**
 * Read reference database from gzipped binary format (.trkb.gz)
 * Populates the same global structures as readReferenceTree()
 *
 * @param filename Path to .trkb.gz file
 * @param name_specs Output array: [max_nodename, max_tax_name, max_line_taxonomy]
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinaryGzipped(const char *filename, int *name_specs) {
    gzFile gz = gzopen(filename, "rb");
    if (!gz) {
        LOG_ERROR("Cannot open gzipped binary reference file: %s", filename);
        return -1;
    }

    // === Validate Header ===
    uint8_t magic[4];
    if (gzread(gz, magic, 4) != 4) {
        LOG_ERROR("Failed to read magic number (gzip)");
        gzclose(gz);
        return -1;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        LOG_ERROR("Invalid binary format magic number (gzip)");
        gzclose(gz);
        return -1;
    }

    // Read version
    uint8_t version_major, version_minor;
    gzread(gz, &version_major, 1);
    gzread(gz, &version_minor, 1);
    LOG_DEBUG("Binary format version: %d.%d (gzip)", version_major, version_minor);

    if (version_major > 1) {
        LOG_ERROR("Unsupported binary format version: %d.%d", version_major, version_minor);
        gzclose(gz);
        return -1;
    }

    // Read flags
    uint8_t endianness, precision;
    gzread(gz, &endianness, 1);
    gzread(gz, &precision, 1);

    if (endianness != 0x01) {
        LOG_ERROR("Only little-endian binary format is supported");
        gzclose(gz);
        return -1;
    }

    if (precision != 0x01) {
        LOG_ERROR("Only float precision binary format is supported");
        gzclose(gz);
        return -1;
    }

    // Skip header CRC and reserved
    gz_read_u32(gz);  // header_crc
    gz_read_u32(gz);  // reserved

    // Read section offsets (we still need these for validation, but we read sequentially)
    uint64_t taxonomy_offset = gz_read_u64(gz);
    uint64_t node_offset = gz_read_u64(gz);
    uint64_t posterior_offset = gz_read_u64(gz);
    uint64_t total_size = gz_read_u64(gz);

    LOG_DEBUG("Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu",
              (unsigned long)taxonomy_offset, (unsigned long)node_offset,
              (unsigned long)posterior_offset, (unsigned long)total_size);

    // Skip to global metadata (seek to byte 64)
    gzseek(gz, BINARY_FILE_HEADER_SIZE, SEEK_SET);

    // === Read Global Metadata (16 bytes) ===
    int32_t numberOfTrees = gz_read_i32(gz);
    int32_t max_nodename = gz_read_i32(gz);
    int32_t max_tax_name = gz_read_i32(gz);
    int32_t max_lineTaxonomy = gz_read_i32(gz);

    name_specs[0] = max_nodename;
    name_specs[1] = max_tax_name;
    name_specs[2] = max_lineTaxonomy;

    LOG_DEBUG("Database: %d trees, max_nodename=%d, max_tax=%d (gzip)",
              numberOfTrees, max_nodename, max_tax_name);

    // === Allocate Global Arrays ===
    numbaseArr = (int*)malloc(numberOfTrees * sizeof(int));
    rootArr = (int*)malloc(numberOfTrees * sizeof(int));
    numspecArr = (int*)malloc(numberOfTrees * sizeof(int));

    if (!numbaseArr || !rootArr || !numspecArr) {
        LOG_ERROR("Failed to allocate tree metadata arrays");
        gzclose(gz);
        return -1;
    }

    // === Read Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < numberOfTrees; i++) {
        numbaseArr[i] = gz_read_i32(gz);
        rootArr[i] = gz_read_i32(gz);
        numspecArr[i] = gz_read_i32(gz);

        LOG_DEBUG("Tree %d: numbase=%d, root=%d, numspec=%d",
                  i, numbaseArr[i], rootArr[i], numspecArr[i]);
    }

    // === Read Taxonomy Section ===
    gzseek(gz, taxonomy_offset, SEEK_SET);

    // Allocate taxonomy array
    allocateMemoryForTaxArr(numberOfTrees, max_tax_name);

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t tree_tax_size = gz_read_u32(gz);
        gz_read_u32(gz);  // reserved
        (void)tree_tax_size;

        for (int s = 0; s < numspecArr[t]; s++) {
            for (int l = 0; l < 7; l++) {  // 7 taxonomy levels
                uint16_t len = gz_read_u16(gz);
                if (len > 0 && len <= (uint16_t)(max_tax_name + 1)) {
                    if (gzread(gz, taxonomyArr[t][s][l], len) != (int)len) {
                        LOG_ERROR("Failed to read taxonomy string (gzip)");
                        gzclose(gz);
                        return -1;
                    }
                } else if (len > 0) {
                    LOG_WARN("Taxonomy string length %d exceeds max %d", len, max_tax_name);
                    // Skip oversized string by reading and discarding
                    char *skip_buf = malloc(len);
                    if (skip_buf) {
                        gzread(gz, skip_buf, len);
                        free(skip_buf);
                    }
                }
            }
        }

        LOG_DEBUG("Tree %d: read %d taxonomy entries (gzip)", t, numspecArr[t]);
    }

    // === Read Node Section ===
    gzseek(gz, node_offset, SEEK_SET);

    // Allocate tree array
    treeArr = malloc(numberOfTrees * sizeof(struct node *));
    if (!treeArr) {
        LOG_ERROR("Failed to allocate tree array");
        gzclose(gz);
        return -1;
    }

    for (int t = 0; t < numberOfTrees; t++) {
        uint32_t num_nodes = gz_read_u32(gz);
        int expected_nodes = 2 * numspecArr[t] - 1;

        if ((int)num_nodes != expected_nodes) {
            LOG_ERROR("Node count mismatch for tree %d: got %d, expected %d",
                      t, num_nodes, expected_nodes);
            gzclose(gz);
            return -1;
        }

        // Allocate nodes for this tree
        allocateTreeArrMemory(t, max_nodename);

        // Read node records into temporary array first (for name offset handling)
        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) {
            LOG_ERROR("Failed to allocate name offset array");
            gzclose(gz);
            return -1;
        }

        // Read all node records (32 bytes each)
        for (int n = 0; n < (int)num_nodes; n++) {
            treeArr[t][n].up[0] = gz_read_i32(gz);
            treeArr[t][n].up[1] = gz_read_i32(gz);
            treeArr[t][n].down = gz_read_i32(gz);
            treeArr[t][n].depth = gz_read_i32(gz);
            treeArr[t][n].taxIndex[0] = gz_read_i32(gz);
            treeArr[t][n].taxIndex[1] = gz_read_i32(gz);
            name_offsets[n] = gz_read_u32(gz);
            gz_read_u32(gz);  // reserved
        }

        // Read name table - needs special handling for gzip
        // Store current position, then read names sequentially
        z_off_t name_table_start = gztell(gz);

        // For gzip, we can't efficiently seek backwards, so we read names in offset order
        // First, find all non-zero offsets and their indices
        typedef struct { uint32_t offset; int node_idx; } offset_entry_t;
        int name_count = 0;
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) name_count++;
        }

        offset_entry_t *sorted_offsets = malloc(name_count * sizeof(offset_entry_t));
        if (!sorted_offsets && name_count > 0) {
            LOG_ERROR("Failed to allocate offset sort array");
            free(name_offsets);
            gzclose(gz);
            return -1;
        }

        int idx = 0;
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                sorted_offsets[idx].offset = name_offsets[n];
                sorted_offsets[idx].node_idx = n;
                idx++;
            }
        }

        // Sort by offset (simple bubble sort - name_count is small)
        for (int i = 0; i < name_count - 1; i++) {
            for (int j = 0; j < name_count - i - 1; j++) {
                if (sorted_offsets[j].offset > sorted_offsets[j + 1].offset) {
                    offset_entry_t tmp = sorted_offsets[j];
                    sorted_offsets[j] = sorted_offsets[j + 1];
                    sorted_offsets[j + 1] = tmp;
                }
            }
        }

        // Read names in order
        for (int i = 0; i < name_count; i++) {
            gzseek(gz, name_table_start + sorted_offsets[i].offset - 1, SEEK_SET);
            uint16_t len = gz_read_u16(gz);
            int n = sorted_offsets[i].node_idx;
            if (len > 0 && len <= (uint16_t)(max_nodename + 1)) {
                if (gzread(gz, treeArr[t][n].name, len) != (int)len) {
                    LOG_ERROR("Failed to read node name (gzip)");
                    free(sorted_offsets);
                    free(name_offsets);
                    gzclose(gz);
                    return -1;
                }
            }
        }

        free(sorted_offsets);
        free(name_offsets);

        LOG_DEBUG("Tree %d: read %d node structures (gzip)", t, num_nodes);

        // TSV logging
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
    gzseek(gz, posterior_offset, SEEK_SET);

    LOG_INFO("Loading posteriors from gzipped binary format...");

    for (int t = 0; t < numberOfTrees; t++) {
        int num_nodes = 2 * numspecArr[t] - 1;
        int numbase = numbaseArr[t];

        for (int n = 0; n < num_nodes; n++) {
            size_t count = numbase * 4;

#ifdef OPTIMIZE_MEMORY
            // Direct read - posteriors are float
            if (gzread(gz, treeArr[t][n].posteriornc, count * sizeof(float)) != (int)(count * sizeof(float))) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (gzip)", t, n);
                gzclose(gz);
                return -1;
            }
#else
            // Need to convert float to double
            float *temp = malloc(count * sizeof(float));
            if (!temp) {
                LOG_ERROR("Failed to allocate temp buffer for posterior conversion");
                gzclose(gz);
                return -1;
            }
            if (gzread(gz, temp, count * sizeof(float)) != (int)(count * sizeof(float))) {
                LOG_ERROR("Failed to read posteriors for tree %d node %d (gzip)", t, n);
                free(temp);
                gzclose(gz);
                return -1;
            }
            for (size_t i = 0; i < count; i++) {
                treeArr[t][n].posteriornc[i] = (double)temp[i];
            }
            free(temp);
#endif
        }

        LOG_DEBUG("Tree %d: read posteriors for %d nodes (gzip)", t, num_nodes);

        // TSV logging
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

    gzclose(gz);

    LOG_INFO("Loaded %d trees from gzipped binary format", numberOfTrees);

    return numberOfTrees;
}
```

#### 2. Add function declaration
**File**: `tronko-assign/readreference.h`
**Changes**: Add after `readReferenceBinary` declaration (line 26)

```c
// Binary format reader (mirrors readReferenceTree but for .trkb files)
int readReferenceBinary(const char *filename, int *name_specs);

// Gzipped binary format reader
int readReferenceBinaryGzipped(const char *filename, int *name_specs);
```

### Success Criteria

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] No memory leaks detected with valgrind (basic test)

#### Manual Verification:
- [x] Function signature matches header declaration

---

## Phase 5: Update Dispatch Logic

### Overview
Update the main dispatch logic in tronko-assign.c to handle the new format type.

### Changes Required

#### 1. Update format dispatch
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Update lines 900-927 to handle `FORMAT_BINARY_GZIPPED`

```c
// Detect reference file format
int ref_format = detect_reference_format(opt.reference_file);

if (ref_format == FORMAT_BINARY) {
    // Load uncompressed binary format
    if (opt.verbose_level >= 0) {
        LOG_INFO("Loading binary format reference database: %s", opt.reference_file);
    }
    numberOfTrees = readReferenceBinary(opt.reference_file, name_specs);
    if (numberOfTrees < 0) {
        printf("Error: Failed to load binary reference file: %s. Exiting...\n", opt.reference_file);
        exit(-1);
    }
} else if (ref_format == FORMAT_BINARY_GZIPPED) {
    // Load gzipped binary format
    if (opt.verbose_level >= 0) {
        LOG_INFO("Loading gzipped binary format reference database: %s", opt.reference_file);
    }
    numberOfTrees = readReferenceBinaryGzipped(opt.reference_file, name_specs);
    if (numberOfTrees < 0) {
        printf("Error: Failed to load gzipped binary reference file: %s. Exiting...\n", opt.reference_file);
        exit(-1);
    }
} else if (ref_format == FORMAT_TEXT) {
    // Load text format (existing code path)
    if (opt.verbose_level >= 0) {
        LOG_INFO("Loading text format reference database: %s", opt.reference_file);
    }
    gzFile referenceTree = gzopen(opt.reference_file, "r");
    if (referenceTree == Z_NULL) {
        printf("Error: Cannot open reference file: %s. Exiting...\n", opt.reference_file);
        exit(-1);
    }
    numberOfTrees = readReferenceTree(referenceTree, name_specs);
    gzclose(referenceTree);
} else {
    printf("Error: Unknown reference file format: %s. Exiting...\n", opt.reference_file);
    exit(-1);
}
```

### Success Criteria

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Assignment with gzipped binary produces correct output

#### Manual Verification:
- [x] Log messages correctly identify format type

---

## Phase 6: Testing

### Overview
Create and run tests to verify the implementation works correctly.

### Test Cases

#### 1. Create test gzipped binary file
```bash
cd /home/jimjeffers/Work/tronko

# First ensure we have a binary file
ls -la tronko-build/example_datasets/single_tree/reference_tree.trkb

# Create gzipped version
gzip -c tronko-build/example_datasets/single_tree/reference_tree.trkb > /tmp/reference_tree.trkb.gz
```

#### 2. Test format detection
```bash
# Run with verbose logging to see format detection
./tronko-assign/tronko-assign -v 0 -r \
  -f /tmp/reference_tree.trkb.gz \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_results_gzbin.txt -w 2>&1 | head -20
```

#### 3. Compare results
```bash
# Run with uncompressed binary
./tronko-assign/tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.trkb \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_results_bin.txt -w

# Compare outputs
diff /tmp/test_results_bin.txt /tmp/test_results_gzbin.txt
echo "Exit code: $?"  # Should be 0 (no differences)
```

#### 4. Performance comparison
```bash
# Time uncompressed binary
time ./tronko-assign/tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.trkb \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_bin.txt -w

# Time gzipped binary
time ./tronko-assign/tronko-assign -r \
  -f /tmp/reference_tree.trkb.gz \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_gzbin.txt -w
```

### Success Criteria

#### Automated Verification:
- [x] All test commands complete without errors
- [x] `diff` shows no differences between gzipped and uncompressed binary results
- [ ] Valgrind shows no memory leaks: `valgrind --leak-check=full ./tronko-assign/tronko-assign ...`

#### Manual Verification:
- [x] Performance overhead is acceptable (< 2x slower than uncompressed)
- [x] Memory usage is reasonable during decompression

---

## Performance Considerations

- **Decompression overhead**: `gzread()` adds CPU overhead for decompression. For large files, this is offset by reduced I/O.
- **Seek behavior**: `gzseek()` is efficient for forward seeks but slow for backwards. The name table handling sorts offsets to ensure forward-only access.
- **Memory usage**: Gzip decompression uses a small internal buffer (~32KB). No additional large allocations are needed.
- **Buffer size**: zlib uses default buffer sizes that work well for streaming reads.

---

## Migration Notes

No migration needed. Users can:
1. Continue using existing `.trkb` or `.txt.gz` files
2. Optionally gzip their `.trkb` files: `gzip reference_tree.trkb`
3. The new format is detected automatically

---

## References

- Original research: `thoughts/shared/research/2025-12-31-gzipped-binary-format-support.md`
- Format detection: `tronko-assign/readreference.c:52`
- Binary reader: `tronko-assign/readreference.c:649`
- Text reader: `tronko-assign/readreference.c:459`
- Format constants: `tronko-assign/global.h:51-53`
- Dispatch logic: `tronko-assign/tronko-assign.c:900-927`
