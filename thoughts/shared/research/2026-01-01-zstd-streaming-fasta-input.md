---
date: 2026-01-01T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Streaming zstd-compressed FASTA files in tronko-assign"
tags: [research, codebase, zstd, compression, fasta, streaming, performance]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Streaming zstd-compressed FASTA files in tronko-assign

**Date**: 2026-01-01T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

What would be involved in streaming zstd-compressed FASTA inputs to tronko-assign instead of pre-decompressing files to disk?

## Summary

**Good news**: tronko-assign already has a clean abstraction layer using zlib's `gzFile` type and `gzgets()` for transparent gzip decompression. Adding zstd support is feasible with **three main approaches** ranging from minimal effort to full custom implementation:

| Approach | Effort | Pros | Cons |
|----------|--------|------|------|
| **1. zstd's zlib wrapper** | Low (~1 day) | Drop-in replacement, minimal code changes | Adds both libraries as deps, ~10% overhead |
| **2. Auto-detect wrapper** | Medium (~2-3 days) | Best UX, handles all formats | More complex file opening logic |
| **3. Custom zstd streaming** | High (~1 week) | Most efficient, full control | Requires implementing line buffering |

**Recommendation**: Start with **Approach 1** (zstd zlib wrapper) for quick wins, then migrate to **Approach 2** for better UX.

## Detailed Findings

### Current Architecture: How FASTA Reading Works

tronko-assign uses a **two-pass streaming pattern**:

1. **First pass**: Scan entire file to determine max sequence/name lengths
2. **Allocate**: Pre-allocate memory based on those dimensions
3. **Second pass**: Read batches of sequences into pre-allocated buffers

All file I/O uses zlib's `gzFile` type for transparent gzip decompression:

```c
// Current pattern in tronko-assign.c:1123
gzFile reads_file = gzopen(opt.read1_file, "r");
find_specs_for_reads(read_specs, reads_file, opt.fastq);
gzclose(reads_file);

// Main reading uses gzgets() - readreference.c:351
while(gzgets(query_reads, buffer, buffer_size) != NULL) {
    // Process line...
}
```

### Entry Points Requiring Changes

**Primary entry points (query FASTA/FASTQ files):**

| Location | Mode | Function |
|----------|------|----------|
| `tronko-assign.c:1123` | Single-end | First open for spec scanning |
| `tronko-assign.c:1163` | Single-end | Main processing open |
| `tronko-assign.c:1356-1357` | Paired-end | First open for spec scanning |
| `tronko-assign.c:1423-1424` | Paired-end | Main processing open |

**Secondary entry points (reference FASTA for BWA index - less common to compress):**

| Location | Component |
|----------|-----------|
| `bwa_source_files/bwtindex.c:271, 305` | BWA index building |
| `bwa_source_files/utils.c:84` | xzopen wrapper |

### Approach 1: zstd's zlib Wrapper (Recommended Starting Point)

zstd includes a **drop-in zlib wrapper** that provides `gzopen/gzread/gzgets` equivalents with automatic format detection.

**Implementation:**

1. Add zstd as a dependency in Makefile
2. Replace header include:

```c
// Change from:
#include <zlib.h>

// To:
#include <zstd_zlibwrapper.h>
```

3. Link against both zlib and zstd:

```makefile
LDFLAGS = -lz -lzstd
```

**Makefile changes:**

```makefile
# Add to LDFLAGS
LDFLAGS += -lzstd

# Or if using bundled zstd with zlib wrapper:
ZSTD_DIR = zstd
CFLAGS += -I$(ZSTD_DIR)/lib -I$(ZSTD_DIR)/zlibWrapper
LDFLAGS += -L$(ZSTD_DIR)/lib -lzstd
```

**Pros:**
- Minimal code changes (just header swap)
- Auto-detects zstd vs gzip during decompression
- Keeps gzip compatibility

**Cons:**
- Adds zstd library dependency
- Slight overhead from format detection
- Some zlib functions unsupported: `deflateCopy`, `inflateSync`, `inflateBack`

### Approach 2: Auto-Detect Wrapper (Best UX)

Create a unified file handle that auto-detects format based on magic bytes:

```c
// zstd_file.h - Unified compression wrapper
typedef struct {
    int format;           // 0=plain, 1=gzip, 2=zstd
    union {
        FILE* plain;
        gzFile gz;
        struct {
            FILE* file;
            ZSTD_DCtx* dctx;
            char* in_buf;
            char* out_buf;
            size_t in_size, in_pos;
            size_t out_size, out_pos;
        } zstd;
    } handle;
    char* line_buf;       // For line buffering
    size_t line_buf_size;
} CompressedFile;

CompressedFile* cf_open(const char* path);
char* cf_gets(char* buf, int size, CompressedFile* cf);
void cf_close(CompressedFile* cf);
```

**Format detection using magic bytes:**

| Format | Magic Bytes | Detection |
|--------|-------------|-----------|
| Gzip | `0x1f 0x8b` | First 2 bytes |
| Zstd | `0x28 0xb5 0x2f 0xfd` | First 4 bytes |
| Plain | ASCII text | Default fallback |

**Integration:**

Replace all `gzopen`/`gzgets`/`gzclose` calls with `cf_open`/`cf_gets`/`cf_close`.

### Approach 3: Native zstd Streaming (Maximum Performance)

Implement direct zstd streaming with line buffering:

```c
#include <zstd.h>

typedef struct {
    FILE* file;
    ZSTD_DCtx* dctx;
    void* in_buf;
    void* out_buf;
    size_t in_buf_size;
    size_t out_buf_size;
    size_t out_pos;
    size_t out_end;
    int eof;
} ZstdFile;

ZstdFile* zstd_fopen(const char* path) {
    ZstdFile* zf = calloc(1, sizeof(ZstdFile));
    zf->file = fopen(path, "rb");
    zf->dctx = ZSTD_createDCtx();
    zf->in_buf_size = ZSTD_DStreamInSize();
    zf->out_buf_size = ZSTD_DStreamOutSize();
    zf->in_buf = malloc(zf->in_buf_size);
    zf->out_buf = malloc(zf->out_buf_size);
    return zf;
}

// Refill output buffer from compressed input
static int zstd_refill(ZstdFile* zf) {
    if (zf->eof) return 0;

    size_t read = fread(zf->in_buf, 1, zf->in_buf_size, zf->file);
    if (read == 0) { zf->eof = 1; return 0; }

    ZSTD_inBuffer in = { zf->in_buf, read, 0 };
    ZSTD_outBuffer out = { zf->out_buf, zf->out_buf_size, 0 };

    while (in.pos < in.size) {
        size_t ret = ZSTD_decompressStream(zf->dctx, &out, &in);
        if (ZSTD_isError(ret)) return -1;
    }

    zf->out_pos = 0;
    zf->out_end = out.pos;
    return out.pos;
}

// Read line (gzgets equivalent)
char* zstd_gets(char* buf, int size, ZstdFile* zf) {
    int i = 0;
    while (i < size - 1) {
        if (zf->out_pos >= zf->out_end) {
            if (zstd_refill(zf) <= 0) break;
        }
        char c = ((char*)zf->out_buf)[zf->out_pos++];
        buf[i++] = c;
        if (c == '\n') break;
    }
    if (i == 0) return NULL;
    buf[i] = '\0';
    return buf;
}

void zstd_fclose(ZstdFile* zf) {
    ZSTD_freeDCtx(zf->dctx);
    free(zf->in_buf);
    free(zf->out_buf);
    fclose(zf->file);
    free(zf);
}
```

### Memory Requirements

**zstd streaming decompression allocates internally:**
- History buffer: `Window_Size + 2 * Block_Maximum_Size`
- Default max window: 128 MB (can be limited via `ZSTD_d_maxWindowLog`)

**Recommended buffer sizes:**
```c
size_t in_size = ZSTD_DStreamInSize();   // ~128 KB
size_t out_size = ZSTD_DStreamOutSize(); // ~128 KB
```

**Limiting memory usage:**
```c
ZSTD_DCtx* dctx = ZSTD_createDCtx();
ZSTD_DCtx_setParameter(dctx, ZSTD_d_maxWindowLog, 23);  // Max 8 MB
```

### Existing Patterns to Follow

The codebase already has compression handling patterns in:

| Location | Pattern |
|----------|---------|
| `readreference.c:89-141` | Magic byte format detection |
| `readreference.c:48-83` | `gz_read_*` binary helpers |
| `readreference.c:991-1313` | Full gzipped binary reader |
| `bwa_source_files/utils.c:75-88` | Error-checked `xzopen` wrapper |

**The format detection pattern is particularly relevant:**

```c
// Existing pattern in readreference.c:89
int detect_reference_format(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    uint8_t magic[4];
    fread(magic, 1, 4, fp);
    fclose(fp);

    // Gzip: 0x1f 0x8b
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return FORMAT_GZIP;
    }
    // Add: Zstd: 0x28 0xb5 0x2f 0xfd
    if (magic[0] == 0x28 && magic[1] == 0xb5 &&
        magic[2] == 0x2f && magic[3] == 0xfd) {
        return FORMAT_ZSTD;
    }
    return FORMAT_PLAIN;
}
```

## Code References

| File | Lines | Description |
|------|-------|-------------|
| `tronko-assign/tronko-assign.c` | 1123, 1163 | Single-end FASTA opens |
| `tronko-assign/tronko-assign.c` | 1356-1357, 1423-1424 | Paired-end FASTA opens |
| `tronko-assign/readreference.c` | 314-500 | `readInXNumberOfLines()` batch reader |
| `tronko-assign/readreference.c` | 1335-1371 | `find_specs_for_reads()` first-pass scanner |
| `tronko-assign/readreference.c` | 89-141 | `detect_reference_format()` magic byte detection |
| `tronko-assign/global.h` | 45 | `FASTA_MAXLINE` constant (40000) |
| `tronko-assign/Makefile` | - | Build configuration |

## Implementation Roadmap

### Phase 1: zlib Wrapper (Quick Win)

1. Add zstd as git submodule or system dependency
2. Update Makefile to link zstd
3. Replace `#include <zlib.h>` with wrapper header
4. Test with `.zst` compressed FASTA files

### Phase 2: Auto-Detection (Better UX)

1. Create `compressed_file.h/c` wrapper module
2. Implement magic byte detection
3. Replace `gzopen`/`gzgets`/`gzclose` calls
4. Add format detection to command line help

### Phase 3: Optimization (If Needed)

1. Profile decompression overhead
2. Consider parallel decompression if bottleneck
3. Pre-buffer decompression in separate thread

## Open Questions

1. **Reference FASTA compression**: Should BWA index building also support zstd? (Less common use case)
2. **Compression level**: What zstd compression level provides best size/speed tradeoff for FASTA?
3. **Seekable format**: Would random access be valuable? (zstd has a seekable format extension)
4. **CI/CD**: Add zstd-compressed test files to the test suite?

## Related Research

- `2025-12-29-tronko-assign-streaming-architecture.md` - Overall streaming architecture
- `2025-12-31-gzipped-binary-format-support.md` - Similar work for reference database compression
