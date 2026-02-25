---
date: 2025-12-31T22:55:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Gzipped Binary Format Support for Reference Databases"
tags: [research, codebase, tronko-assign, binary-format, compression, performance]
status: complete
last_updated: 2025-12-31
last_updated_by: Claude
---

# Research: Gzipped Binary Format Support for Reference Databases

**Date**: 2025-12-31T22:55:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Our existing benchmarks tested tronko-assign on uncompressed reference trees vs the new binary format (.trkb). However, production often uses gzipped text format (.txt.gz). Questions:
1. Will binary format still be more performant than gzipped text?
2. Can we support loading gzipped binary format (.trkb.gz)?
3. How much smaller will gzipped binary be compared to gzipped text?

## Summary

**Gzipped binary is 2.15x smaller than gzipped text** (352 MB vs 758 MB for the 16S_Bacteria dataset). Binary format will still outperform gzipped text because it requires less decompression and no string parsing. Supporting gzipped binary is feasible with straightforward zlib integration changes.

## Detailed Findings

### File Size Comparison

Using the 16S_Bacteria reference database:

| Format | Size | Relative to gzipped text |
|--------|------|--------------------------|
| Gzipped text (`.txt.gz`) | 758 MB | 1.00x |
| Uncompressed binary (`.trkb`) | 5.8 GB | 7.65x |
| **Gzipped binary (`.trkb.gz`)** | **352 MB** | **0.46x** |

The binary format compresses better because:
- Raw floats/integers have better entropy for compression than ASCII text
- No redundant whitespace or delimiters
- Structured data layout allows more efficient compression patterns

### Performance Analysis

**Will binary still be faster when both formats are gzipped?**

Yes, for these reasons:

1. **Less data to decompress**: 352 MB vs 758 MB = 53% less decompression work
2. **No string parsing**: Binary reads integers/floats directly via `fread()` vs text requiring `sscanf()` for every line
3. **Bulk reads**: Binary reads posteriors in large contiguous chunks vs line-by-line text parsing with string tokenization
4. **CPU efficiency**: Decompression + direct memory copy vs decompression + tokenization + parsing + type conversion

### Current Implementation Status

**Text format reader** (`readReferenceTree()` in `readreference.c:459`):
- Uses `gzFile` from zlib
- Transparently handles both plain and gzipped text
- Already production-ready for `.txt.gz` files

**Binary format reader** (`readReferenceBinary()` in `readreference.c:649`):
- Uses standard `FILE*` with `fopen()`
- Uses `fread()` for data, `fseek()` for section navigation
- **No gzip support currently**

**Format detection** (`detect_reference_format()` in `readreference.c:52`):
- Checks magic bytes at file start
- Binary magic: `0x89 'T' 'R' 'K'`
- Gzip magic: `0x1f 0x8b`
- **Problem**: Currently assumes all gzipped files are text format (line 76-78)

```c
// readreference.c:76-78 - Current logic
if (magic[0] == 0x1f && magic[1] == 0x8b) {
    LOG_DEBUG("Detected gzipped text format: %s", filename);
    return FORMAT_TEXT;  // Assumes gzip = text, wrong for .trkb.gz!
}
```

### Implementation Plan for Gzipped Binary Support

#### 1. Update Format Detection

Need to decompress first 4 bytes of gzipped files to check the actual content:

```c
if (magic[0] == 0x1f && magic[1] == 0x8b) {
    // Gzipped file - decompress first 4 bytes to check content
    gzFile gz = gzopen(filename, "rb");
    if (gz) {
        uint8_t inner_magic[4];
        if (gzread(gz, inner_magic, 4) == 4) {
            gzclose(gz);
            if (inner_magic[0] == TRONKO_MAGIC_0 && inner_magic[1] == TRONKO_MAGIC_1 &&
                inner_magic[2] == TRONKO_MAGIC_2 && inner_magic[3] == TRONKO_MAGIC_3) {
                return FORMAT_BINARY_GZIPPED;  // New format type
            }
            return FORMAT_TEXT;  // Gzipped text
        }
        gzclose(gz);
    }
}
```

#### 2. Modify Binary Reader for Gzip Support

Convert `readReferenceBinary()` to use zlib:

| Current | New |
|---------|-----|
| `FILE *fp = fopen(...)` | `gzFile gz = gzopen(...)` |
| `fread(buf, 1, n, fp)` | `gzread(gz, buf, n)` |
| `fseek(fp, offset, SEEK_SET)` | `gzseek(gz, offset, SEEK_SET)` |
| `fclose(fp)` | `gzclose(gz)` |

**Seek considerations**:
- `gzseek()` is efficient for forward seeks
- All current seeks in `readReferenceBinary()` are forward seeks to section offsets
- The name table reading has some local seeking that may need restructuring

#### 3. Alternative: Memory-Mapped Decompression

For maximum performance, could decompress entire file to memory first:
```c
// Decompress to memory buffer, then use buffer like mmap
void *buf = decompress_entire_file(filename, &size);
// Use memcpy instead of fread, pointer arithmetic instead of fseek
```

This avoids repeated `gzread()` overhead but requires more memory.

## Code References

- `tronko-assign/readreference.c:52` - `detect_reference_format()` function
- `tronko-assign/readreference.c:76-78` - Gzip detection logic (needs update)
- `tronko-assign/readreference.c:459` - `readReferenceTree()` text format reader
- `tronko-assign/readreference.c:649` - `readReferenceBinary()` binary format reader
- `tronko-assign/tronko-assign.c:900-927` - Format detection and dispatch logic
- `tronko-assign/readreference.h:23-26` - Function declarations

## Architecture Insights

The binary format uses a sectioned layout with offset table in the header:
1. **Header** (64 bytes): Magic, version, flags, section offsets
2. **Global metadata**: Tree count, max name lengths
3. **Tree metadata**: Per-tree numbase, root, numspec
4. **Taxonomy section**: Length-prefixed strings
5. **Node section**: Fixed-size node records + name table
6. **Posterior section**: Bulk float arrays

This design enables efficient seeking - the reader jumps directly to each section using the offset table. For gzip support, forward-only seeking works well since sections are read in order.

## Open Questions

1. **Memory vs streaming tradeoff**: For very large databases, is it better to stream decompress or load entire compressed file to memory first?

2. **Compression level**: Should we recommend a specific gzip compression level for `.trkb.gz` files? Higher compression = smaller file but slower decompression.

3. **zstd alternative**: zstd offers better compression ratios and faster decompression than gzip. Worth considering as a future enhancement.

## Recommendation

**Implement gzipped binary support**. Benefits clearly outweigh implementation effort:
- **Storage**: 2.15x smaller than gzipped text (significant for cloud storage costs)
- **Performance**: Still faster than gzipped text (less to decompress + no parsing)
- **Implementation**: Straightforward zlib conversion, ~50-100 lines of code changes
