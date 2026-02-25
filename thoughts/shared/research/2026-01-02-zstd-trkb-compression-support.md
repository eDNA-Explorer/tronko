---
date: 2026-01-02T10:30:00-08:00
researcher: Claude
git_commit: 8fa99b63db2052b40a62e61b9726a6bb6960c389
branch: experimental
repository: tronko
topic: "ZSTD compression support for TRKB binary format"
tags: [research, codebase, tronko-convert, tronko-assign, zstd, compression, binary-format, performance]
status: complete
last_updated: 2026-01-02
last_updated_by: Claude
---

# Research: ZSTD Compression Support for TRKB Binary Format

**Date**: 2026-01-02T10:30:00-08:00
**Researcher**: Claude
**Git Commit**: 8fa99b63db2052b40a62e61b9726a6bb6960c389
**Branch**: experimental
**Repository**: tronko

## Research Question

What would be involved in:
1. Updating tronko-convert to output zstd-compressed `.trkb` files by default
2. Updating tronko-assign to stream decompress zstd-compressed `.trkb` files

## Summary

The infrastructure for zstd support **already exists** in tronko-assign via `compressed_io.h`. Adding zstd compression to the TRKB format requires:

- **tronko-convert**: ~100 lines - add zstd compression wrapper and make it the default output
- **tronko-assign**: ~300 lines - add `readReferenceBinaryZstd()` function mirroring the existing gzip version
- **Total**: ~400 lines of new code

The changes are straightforward because the zstd streaming patterns are already implemented for query file handling.

## Motivation: Benchmark Results

Testing with the 16S_Bacteria reference database showed compelling results:

| Format | Size | vs Text.gz | Decompress Time |
|--------|------|------------|-----------------|
| Text + gzip | 758 MB | 1.0x | 2.53s |
| Binary + gzip | 352 MB | 2.2x smaller | 2.53s |
| Binary + zstd -3 | 420 MB | 1.8x smaller | 2.09s |
| **Binary + zstd -19** | **195 MB** | **3.9x smaller** | **1.01s** |

Key findings:
- Zstd -19 produces files **74% smaller** than gzipped text
- Decompression is **2.5x faster** than gzip despite smaller size
- Binary format compresses better because posterior floats have repetitive byte patterns

## Detailed Findings

### Current TRKB Format Support in tronko-assign

The codebase already supports multiple reference database formats:

| Format | Constant | Reader Function | File |
|--------|----------|-----------------|------|
| Text (.txt) | `FORMAT_TEXT` | `readReferenceTree()` | `readreference.c:511` |
| Text gzipped (.txt.gz) | `FORMAT_TEXT` | `readReferenceTree()` | (gzFile handles transparently) |
| Binary (.trkb) | `FORMAT_BINARY` | `readReferenceBinary()` | `readreference.c:701` |
| Binary gzipped (.trkb.gz) | `FORMAT_BINARY_GZIPPED` | `readReferenceBinaryGzipped()` | `readreference.c:991` |
| **Binary zstd (.trkb.zst)** | - | - | **Not implemented** |

### Existing ZSTD Infrastructure

`compressed_io.h` provides complete zstd streaming decompression:

```c
// CompressedFile abstraction (compressed_io.h:28-44)
typedef struct {
    int format;
    union {
        gzFile gz;
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
```

Key functions:
- `cf_open()` - Opens file with auto-detection (line 84)
- `cf_gets()` - Line-by-line reading with streaming decompression (line 174)
- `zstd_refill()` - Refills decompression buffer (line 138)
- `cf_close()` - Cleanup (line 202)

The zstd magic detection is already implemented:
```c
// Zstd: 0x28 0xb5 0x2f 0xfd (compressed_io.h:67-70)
if (bytes_read >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 &&
    magic[2] == 0x2f && magic[3] == 0xfd) {
    return COMPRESS_FORMAT_ZSTD;
}
```

### tronko-convert Current State

The converter currently:
- Links only `-lz -lm` (no zstd) - `Makefile:4`
- Outputs uncompressed binary via `write_binary()` - `format_binary.c:16`
- Uses `FILE*` with `fwrite()` for output

## Code References

### tronko-assign
- `tronko-assign/compressed_io.h:1-218` - Complete zstd streaming infrastructure
- `tronko-assign/readreference.c:89-141` - Format detection function
- `tronko-assign/readreference.c:701-981` - `readReferenceBinary()` (uncompressed)
- `tronko-assign/readreference.c:991-1313` - `readReferenceBinaryGzipped()` (template for zstd version)
- `tronko-assign/global.h:53-60` - Format constants and TRKB magic
- `tronko-assign/Makefile:12` - Already links `-lzstd`

### tronko-convert
- `tronko-convert/format_binary.c:16-242` - Binary writer
- `tronko-convert/Makefile:4` - Currently missing `-lzstd`
- `tronko-convert/tronko-convert.c:23-102` - Main entry point

## Implementation Requirements

### Phase 1: tronko-convert (Writer)

**Changes:**

| File | Change | Lines |
|------|--------|-------|
| `Makefile` | Add `-lzstd` to `LDFLAGS` | 1 |
| `format_common.h` | Add `ZSTD_COMPRESSION_LEVEL 19` constant | 2 |
| `format_binary.c` | Add `write_binary_zstd()` function | ~60 |
| `tronko-convert.c` | Make zstd default, add `-u` for uncompressed | ~20 |

**Writer approach options:**

1. **Buffer-then-compress** (simpler):
   ```c
   // Write to memory buffer, then ZSTD_compress() entire buffer
   void* buffer = malloc(estimated_size);
   // ... write to buffer ...
   size_t compressed_size = ZSTD_compress(compressed, compressed_bound, buffer, size, 19);
   fwrite(compressed, 1, compressed_size, fp);
   ```

2. **Streaming compression** (memory-efficient for large DBs):
   ```c
   ZSTD_CCtx* cctx = ZSTD_createCCtx();
   ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, 19);
   // Use ZSTD_compressStream2() for each section
   ```

Recommendation: Use streaming for memory efficiency with large databases.

### Phase 2: tronko-assign (Reader)

**Changes:**

| File | Change | Lines |
|------|--------|-------|
| `global.h` | Add `FORMAT_BINARY_ZSTD 4` | 1 |
| `readreference.h` | Add `readReferenceBinaryZstd()` prototype | 1 |
| `readreference.c` | Update `detect_reference_format()` for zstd | ~15 |
| `readreference.c` | Add `readReferenceBinaryZstd()` function | ~300 |
| `tronko-assign.c` | Add case for `FORMAT_BINARY_ZSTD` | 5 |

**Format detection update:**
```c
// In detect_reference_format():
// Check for zstd magic: 0x28 0xb5 0x2f 0xfd
if (magic[0] == 0x28 && magic[1] == 0xb5 &&
    magic[2] == 0x2f && magic[3] == 0xfd) {
    // Need to decompress first 4 bytes to check for TRKB magic
    // Use ZSTD_decompress() on small buffer
    return FORMAT_BINARY_ZSTD;
}
```

**Reader approach:**

Two options for `readReferenceBinaryZstd()`:

1. **Extend CompressedFile** - Add `cf_read()` for bulk binary reads:
   ```c
   size_t cf_read(void* buf, size_t size, size_t nmemb, CompressedFile* cf);
   ```
   This would allow reusing the existing zstd infrastructure.

2. **Direct ZSTD_DCtx** - Mirror the gzip version but with zstd:
   ```c
   // Similar to readReferenceBinaryGzipped() but using:
   ZSTD_DCtx* dctx = ZSTD_createDCtx();
   ZSTD_inBuffer in = { compressed_buf, bytes_read, 0 };
   ZSTD_outBuffer out = { output_buf, output_size, 0 };
   ZSTD_decompressStream(dctx, &out, &in);
   ```

Recommendation: Option 2 is more straightforward - copy the gzip version and adapt.

## Architecture Insights

### Why Binary Compresses Better Than Text

The posterior section dominates file size (~99%). Binary compresses better because:

1. **Repetitive float patterns**: Many probabilities are exactly 0.0 or 1.0
2. **Related nodes have similar values**: Tree structure creates redundancy
3. **IEEE 754 byte patterns**: Binary floats have predictable byte sequences

Text format suffers from:
- Variable-length decimal representations
- Random digit patterns
- Redundant ASCII delimiters

### Design Decision: Compression Level

Zstd -19 (max) is recommended because:
- Reference databases are write-once, read-many
- Compression time is negligible vs. build time
- Decompression is still fast (faster than gzip)
- Storage/transfer savings are substantial

Can add `-c <level>` flag for users who want faster compression.

## Testing Strategy

1. **Round-trip validation**:
   ```bash
   tronko-convert -i ref.txt -o ref.trkb.zst
   tronko-convert -i ref.trkb.zst -o roundtrip.txt -t
   diff ref.txt roundtrip.txt
   ```

2. **Format detection**: Test all four format combinations
3. **Large file test**: Use 16S_Bacteria database
4. **Memory profiling**: Ensure streaming doesn't exceed memory limits

## Open Questions

1. **Should we support multiple compression levels?**
   - Recommendation: Default to 19, add `-c` flag for customization

2. **Should uncompressed `.trkb` still be supported?**
   - Recommendation: Yes, with `-u` flag for debugging/mmap future

3. **File extension convention?**
   - Option A: `.trkb` always (detect by magic)
   - Option B: `.trkb.zst` for compressed
   - Recommendation: Option A (simpler, follows gzip pattern for text)

## Related Research

- `thoughts/shared/research/2025-12-30-binary-format-conversion-tool.md` - Original TRKB format design
- `thoughts/shared/plans/2025-12-31-binary-format-converter-phase1.md` - tronko-convert implementation

## References

- `tronko-convert/FORMAT_SPECIFICATION.md` - Binary format specification
- [Zstd documentation](https://facebook.github.io/zstd/zstd_manual.html)
- [Zstd streaming API](https://github.com/facebook/zstd/blob/dev/lib/zstd.h)
