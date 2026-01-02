# TRKB Binary Format Specification

**Version**: 1.0
**Status**: Implemented
**Last Updated**: 2026-01-02

This document provides a comprehensive specification of the `.trkb` binary format used by tronko-convert, including design rationale and performance characteristics.

## Table of Contents

1. [Overview](#overview)
2. [Performance Benefits](#performance-benefits)
3. [File Structure](#file-structure)
4. [Section Specifications](#section-specifications)
5. [Design Decisions](#design-decisions)
6. [Implementation Notes](#implementation-notes)

---

## Overview

The `.trkb` (Tronko Binary) format is a custom binary serialization format for tronko reference databases. It replaces the text-based `reference_tree.txt` format produced by tronko-build, providing significant improvements in loading performance and storage efficiency.

### Key Characteristics

| Property | Value |
|----------|-------|
| Magic Number | `0x89 'T' 'R' 'K'` (4 bytes) |
| Byte Order | Little-endian (all multi-byte values) |
| Posterior Precision | 32-bit float (IEEE 754) |
| Checksum | CRC-32 (IEEE 802.3 polynomial) |
| Compression | None (designed for mmap compatibility) |

### Format Detection

Files are automatically detected by magic bytes:

| First Bytes | Format |
|-------------|--------|
| `0x89 0x54 0x52 0x4B` | Binary (.trkb) |
| `0x1f 0x8b` | Gzipped text |
| ASCII digits (`0x30`-`0x39`) | Plain text |

---

## Performance Benefits

The binary format achieves substantial performance improvements over text format:

### 1. File Size Reduction (~2.7x smaller)

**Root Cause**: Posterior probabilities dominate the database size. Each alignment position stores 4 probability values (A, C, G, T) for every node in the tree.

| Format | Single Position Storage | For 4 Posteriors |
|--------|------------------------|------------------|
| Text | `0.12345678901234567` (~18 chars) | ~72 bytes + delimiters |
| Binary | 4-byte float | 16 bytes |

**Additional savings**:
- No tab/newline delimiters between values
- 32-bit float vs 64-bit double (text uses `%.17g` format)
- Length-prefixed strings vs delimiter scanning

**Example** (1466 species, 316 positions):
- Text: ~40 MB
- Binary: ~15 MB
- Compression ratio: 2.7x

### 2. Loading Time Reduction (50-70x faster)

**Root Cause**: Text parsing requires millions of function calls for string-to-number conversion.

#### Text Format Loading Path
```
For each of ~926,000 lines:
  gzgets()     → read line into buffer
  strtok()     → tokenize (called ~4x per line)
  sscanf()     → parse 4 doubles from ASCII
  (float)cast  → convert double to float
```

For the example database:
- ~926,000 `gzgets()` calls
- ~3,700,000 `strtok()` calls
- ~3,700,000 `sscanf()` calls (the bottleneck)
- ~930,000 `malloc()` calls

#### Binary Format Loading Path
```
For each node:
  fread()  → bulk read entire posterior array (single syscall)
```

For the example database:
- ~2,931 `fread()` calls for posteriors (one per node)
- ~10 large allocations

**Why this matters**: `sscanf()` with floating-point format specifiers is computationally expensive. It must:
1. Skip whitespace
2. Parse sign
3. Parse integer part digit-by-digit
4. Parse decimal point
5. Parse fractional part digit-by-digit
6. Parse exponent (if present)
7. Assemble into IEEE 754 representation

Binary format bypasses all of this - the bytes are already in IEEE 754 format.

### 3. Memory Efficiency

#### During Loading

| Aspect | Text | Binary |
|--------|------|--------|
| Line buffer | 4 KB | Not needed |
| Temporary doubles | 4 x 8 bytes per line | Not needed |
| zlib buffers | ~32 KB | Not needed |
| String tokenization | Dynamic | Not needed |

#### Runtime (In-Memory)

| Component | Text (double) | Binary (float) | Savings |
|-----------|--------------|----------------|---------|
| Posteriors/node | numbase × 4 × 8 bytes | numbase × 4 × 4 bytes | 50% |

For example database (316 positions, 2931 nodes):
- Text: 2931 × 316 × 4 × 8 = ~29.6 MB
- Binary: 2931 × 316 × 4 × 4 = ~14.8 MB

### Performance Summary

| Metric | Text Format | Binary Format | Improvement |
|--------|-------------|---------------|-------------|
| File size | ~40 MB | ~15 MB | 2.7x smaller |
| Load time | ~7 seconds | ~0.1 seconds | 50-70x faster |
| Peak memory | ~50 MB | ~20 MB | 2.5x less |
| I/O syscalls | ~1M | ~3K | 300x fewer |

---

## File Structure

```
+---------------------------+
| File Header (64 bytes)    |  Fixed size, contains offsets
+---------------------------+
| Global Metadata (16 bytes)|  Database-wide parameters
+---------------------------+
| Tree Metadata Section     |  12 bytes per tree
+---------------------------+
| Taxonomy String Table     |  Variable size, length-prefixed
+---------------------------+
| Node Structure Section    |  32 bytes per node + name table
+---------------------------+
| Posterior Data Section    |  Bulk float arrays (largest section)
+---------------------------+
| Footer (8 bytes)          |  CRC-32 + magic
+---------------------------+
```

---

## Section Specifications

### File Header (64 bytes)

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | uint8[4] | magic | `0x89 'T' 'R' 'K'` |
| 4 | 1 | uint8 | version_major | Format version (1) |
| 5 | 1 | uint8 | version_minor | Format version (0) |
| 6 | 1 | uint8 | endianness | `0x01` = little-endian |
| 7 | 1 | uint8 | precision | `0x01` = float32 |
| 8 | 4 | uint32 | header_crc | CRC-32 of bytes 0-7 |
| 12 | 4 | uint32 | reserved | Alignment padding |
| 16 | 8 | uint64 | taxonomy_offset | Byte offset to taxonomy section |
| 24 | 8 | uint64 | node_offset | Byte offset to node section |
| 32 | 8 | uint64 | posterior_offset | Byte offset to posterior section |
| 40 | 8 | uint64 | total_size | Total file size in bytes |
| 48 | 16 | - | reserved | Future use |

**Design note**: The magic number starts with `0x89` (non-ASCII) to distinguish binary from text files, followed by readable "TRK" for easy identification in hex editors.

### Global Metadata (16 bytes)

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | int32 | num_trees | Number of trees in database |
| 4 | 4 | int32 | max_nodename | Maximum node name length |
| 8 | 4 | int32 | max_tax_name | Maximum taxonomy name length |
| 12 | 4 | int32 | max_line_taxonomy | Maximum taxonomy line length |

### Tree Metadata (12 bytes per tree)

For each tree `i` in `[0, num_trees)`:

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | int32 | numbase | MSA alignment length |
| 4 | 4 | int32 | root | Root node index |
| 8 | 4 | int32 | numspec | Number of species (leaf nodes) |

**Derived value**: `num_nodes = 2 * numspec - 1` (binary tree property)

### Taxonomy Section

Per tree:

| Size | Type | Description |
|------|------|-------------|
| 4 | uint32 | Taxonomy data size for this tree |
| 4 | uint32 | Reserved |

Per species (numspec entries), per level (7 levels):

| Size | Type | Description |
|------|------|-------------|
| 2 | uint16 | String length (including null terminator) |
| N | char[N] | Null-terminated taxonomy string |

**Taxonomy levels**: domain, phylum, class, order, family, genus, species

### Node Section

Per tree:

| Size | Type | Description |
|------|------|-------------|
| 4 | uint32 | Number of nodes (should equal `2 * numspec - 1`) |

Per node (32 bytes fixed):

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | int32 | up[0] | First child index (-1 for leaf) |
| 4 | 4 | int32 | up[1] | Second child index (-1 for leaf) |
| 8 | 4 | int32 | down | Parent index |
| 12 | 4 | int32 | depth | Tree depth |
| 16 | 4 | int32 | taxIndex[0] | Species/taxonomy row index |
| 20 | 4 | int32 | taxIndex[1] | Taxonomy level |
| 24 | 4 | uint32 | name_offset | Offset into name table (0 = no name) |
| 28 | 4 | uint32 | reserved | Future use |

**Name table** (follows node records):

For each leaf node with a name:

| Size | Type | Description |
|------|------|-------------|
| 2 | uint16 | Name length (including null) |
| N | char[N] | Null-terminated name |

### Posterior Section

This section contains the bulk of the file data. Per tree, per node, per position:

```
float[4] posteriors;  // A, C, G, T probabilities
```

**Memory layout**: `posteriors[node][position][nucleotide]`
- Total floats per tree: `num_nodes × numbase × 4`
- Total bytes per tree: `num_nodes × numbase × 16`

**Why floats instead of doubles**:
- Phylogenetic placement algorithms don't require double precision
- tronko-assign already has `USE_FLOAT_PP` compile option
- Halves the posterior storage size

### Footer (8 bytes)

| Offset | Size | Type | Field | Description |
|--------|------|------|-------|-------------|
| 0 | 4 | uint32 | data_crc | CRC-32 of all preceding bytes |
| 4 | 4 | uint32 | footer_magic | `0x454E4421` ("END!") |

---

## Design Decisions

### Why No Compression?

The format intentionally omits compression to enable:

1. **Memory mapping (mmap)**: Future optimization can map the posterior section directly into memory without decompression
2. **Random access**: Section offsets allow seeking to any part without scanning
3. **Simplicity**: No compression library dependencies beyond what's already used

For storage efficiency, files can be externally compressed (gzip, zstd) and decompressed on load.

### Why Little-Endian Only?

- All modern x86/x64 and ARM processors are little-endian
- Eliminates runtime byte-swapping overhead
- Simplifies implementation

Big-endian systems (rare in bioinformatics) can use the text format.

### Why Fixed-Size Node Records?

32-byte fixed records enable:
- Direct indexing: `offset = base + node_index × 32`
- Predictable memory layout
- Future mmap compatibility

### Why Separate Name Table?

Internal nodes don't have names. Storing names inline would waste space:
- Internal nodes: would need placeholder or variable-size records
- Leaf nodes only: ~50% of nodes have names

The offset-based approach stores names compactly after all node records.

### Why CRC-32?

- Fast to compute
- Standard algorithm (IEEE 802.3 polynomial)
- Sufficient for detecting corruption
- No external library required

---

## Implementation Notes

### Reading Binary Files

```c
// Efficient bulk read for posteriors
for (int n = 0; n < num_nodes; n++) {
    fread(node->posteriors, sizeof(float), numbase * 4, fp);
}
```

### Writing Binary Files

```c
// Bulk write posteriors
fwrite(node->posteriors, sizeof(float), numbase * 4, fp);
```

### Calculating File Size

```c
uint64_t estimate_size(tronko_db_t *db) {
    uint64_t size = 64 + 16;  // Header + global meta
    size += db->num_trees * 12;  // Tree metadata

    for (int t = 0; t < db->num_trees; t++) {
        // Taxonomy: ~avg_name_len * 7 * numspec
        // Nodes: 32 * num_nodes + leaf_names
        // Posteriors: num_nodes * numbase * 16
        size += db->trees[t].num_nodes * db->trees[t].numbase * 16;
    }

    size += 8;  // Footer
    return size;
}
```

### Validation

```c
// Magic number check
if (magic[0] != 0x89 || magic[1] != 'T' ||
    magic[2] != 'R' || magic[3] != 'K') {
    return FORMAT_UNKNOWN;
}

// Version compatibility
if (version_major > SUPPORTED_VERSION_MAJOR) {
    return ERROR_UNSUPPORTED_VERSION;
}

// Footer validation
if (footer_magic != 0x454E4421) {
    return ERROR_CORRUPTED;
}
```

---

## References

### Source Files

| File | Description |
|------|-------------|
| `format_binary.c` | Binary format reader/writer |
| `format_text.c` | Text format reader/writer |
| `format_common.h` | Data structures and constants |
| `crc32.c` | CRC-32 implementation |
| `utils.c` | Error-checked I/O helpers |

### Related Documentation

| Document | Description |
|----------|-------------|
| `README.md` | User documentation and quick reference |
| `thoughts/shared/research/2025-12-30-binary-format-conversion-tool.md` | Original research and analysis |
| `thoughts/shared/plans/2025-12-31-binary-format-converter-phase1.md` | Implementation plan |

### External References

| Reference | Relevance |
|-----------|-----------|
| IEEE 754 | Float representation standard |
| IEEE 802.3 | CRC-32 polynomial |
| tronko-assign `global.h:60-68` | Original node structure |
| tronko-build `printtree.c` | Text format output |
| tronko-assign `readreference.c` | Text format parsing |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2025-12-31 | Initial specification |

---

## Future Considerations

### Potential Enhancements (not implemented)

1. **Memory-mapped loading**: Map posterior section directly without copying
2. **Optional compression**: Per-section compression flags in header
3. **Streaming writes**: Write posteriors during tronko-build computation
4. **Index section**: Enable partial tree loading for very large databases

These would require version 1.1+ and new header flags.
