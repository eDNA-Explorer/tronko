# zstd-Compressed FASTA Input Support Implementation Plan

## Overview

Add support for reading zstd-compressed FASTA/FASTQ query files in tronko-assign using zstd's zlib wrapper for minimal code changes.

## Current State Analysis

tronko-assign uses zlib's `gzFile` type and `gzgets()` for transparent gzip decompression of input files. The pattern is:

```c
gzFile reads_file = gzopen(opt.read1_file, "r");
while(gzgets(reads_file, buffer, buffer_size) != NULL) { ... }
gzclose(reads_file);
```

### Key Files Using zlib for Query Reads

| File | Lines | Purpose |
|------|-------|---------|
| `tronko-assign.c` | 3 | `#include <zlib.h>` |
| `tronko-assign.c` | 1123, 1163 | Single-end FASTA opens |
| `tronko-assign.c` | 1356-1357, 1423-1424 | Paired-end FASTA opens |
| `readreference.h` | 6 | `#include <zlib.h>` |
| `readreference.c` | 115, 178, 351, 1342 | Spec scanning and batch reading |

### BWA Files (Not Modified)

BWA source files also use zlib but for reference FASTA indexing, not query reads. These remain unchanged:
- `bwa_source_files/bntseq.c`
- `bwa_source_files/bwtindex.c`
- `bwa_source_files/utils.c`
- etc.

## Desired End State

Users can pass `.zst` compressed FASTA/FASTQ files directly:

```bash
tronko-assign -f reference.txt -s -g reads.fasta.zst -o results.txt
```

The decompression happens automatically with ~256 KB of buffer memory, streaming through the file without loading it entirely into memory.

### Verification

```bash
# Compress test file
zstd example.fasta -o example.fasta.zst

# Run with compressed input
./tronko-assign -f ref.txt -s -g example.fasta.zst -o /tmp/out.txt

# Verify output matches uncompressed run
diff /tmp/out.txt /tmp/out_uncompressed.txt
```

## What We're NOT Doing

- NOT modifying BWA source files (they handle reference FASTAs, not query reads)
- NOT supporting zstd compression for reference database files (`.trkb`)
- NOT bundling zstd source - using system library
- NOT implementing native zstd streaming (using wrapper for simplicity)

## Implementation Approach

Use zstd's [zlibWrapper](https://github.com/facebook/zstd/tree/dev/zlibWrapper) which provides drop-in replacements for `gzopen`/`gzgets`/`gzclose` that auto-detect and decompress both gzip and zstd formats.

**Key insight**: For decompression, the wrapper auto-detects the format. No code changes needed beyond swapping headers and linking.

---

## Phase 1: Add zstd Dependency

### Overview
Install zstd development headers and update the Makefile to link against libzstd.

### Changes Required

#### 1. Document Dependency
**File**: `README.md`
**Changes**: Add zstd to build requirements

```markdown
## Build Requirements

- gcc
- zlib development headers (`zlib1g-dev` on Debian/Ubuntu)
- zstd development headers (`libzstd-dev` on Debian/Ubuntu)
```

#### 2. Update Makefile
**File**: `tronko-assign/Makefile`
**Changes**: Add `-lzstd` to LIBS

```makefile
# Change from:
LIBS = -lm -pthread -lz -lrt -std=gnu99

# To:
LIBS = -lm -pthread -lz -lzstd -lrt -std=gnu99
```

### Success Criteria

#### Automated Verification:
- [ ] `pkg-config --exists libzstd && echo "zstd found"` returns success
- [ ] `make clean && make` compiles without errors
- [ ] Existing gzip and uncompressed FASTA tests still work

#### Manual Verification:
- [ ] Verify libzstd-dev is documented in build instructions

---

## Phase 2: Create Compression Wrapper Header

### Overview
Create a unified header that wraps zlib functions and adds zstd support via format auto-detection.

### Changes Required

#### 1. Create Wrapper Header
**File**: `tronko-assign/compressed_io.h` (new file)

```c
#ifndef COMPRESSED_IO_H
#define COMPRESSED_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// Detect format from magic bytes
static inline int detect_compression_format(const char* filename) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) return COMPRESS_FORMAT_PLAIN;

    uint8_t magic[4];
    size_t read = fread(magic, 1, 4, fp);
    fclose(fp);

    if (read < 2) return COMPRESS_FORMAT_PLAIN;

    // Gzip: 0x1f 0x8b
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return COMPRESS_FORMAT_GZIP;
    }

    // Zstd: 0x28 0xb5 0x2f 0xfd
    if (read >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 &&
        magic[2] == 0x2f && magic[3] == 0xfd) {
        return COMPRESS_FORMAT_ZSTD;
    }

    return COMPRESS_FORMAT_PLAIN;
}

// Open compressed file for reading
static inline CompressedFile* cf_open(const char* filename, const char* mode) {
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
        cf->handle.zstd.in_buf_size = ZSTD_DStreamInSize();
        cf->handle.zstd.out_buf_size = ZSTD_DStreamOutSize();
        cf->handle.zstd.in_buf = (uint8_t*)malloc(cf->handle.zstd.in_buf_size);
        cf->handle.zstd.out_buf = (uint8_t*)malloc(cf->handle.zstd.out_buf_size);
        cf->handle.zstd.out_pos = 0;
        cf->handle.zstd.out_end = 0;
        cf->handle.zstd.eof = 0;
    } else {
        // gzopen handles both gzip and plain text
        cf->handle.gz = gzopen(filename, mode);
        if (!cf->handle.gz) {
            free(cf);
            return NULL;
        }
    }

    return cf;
}

// Refill zstd output buffer
static inline int zstd_refill(CompressedFile* cf) {
    if (cf->handle.zstd.eof) return 0;

    size_t read = fread(cf->handle.zstd.in_buf, 1,
                        cf->handle.zstd.in_buf_size, cf->handle.zstd.file);
    if (read == 0) {
        cf->handle.zstd.eof = 1;
        return 0;
    }

    ZSTD_inBuffer in = { cf->handle.zstd.in_buf, read, 0 };
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
    return out.pos;
}

// Read line from compressed file (gzgets equivalent)
static inline char* cf_gets(char* buf, int size, CompressedFile* cf) {
    if (cf->format != COMPRESS_FORMAT_ZSTD) {
        return gzgets(cf->handle.gz, buf, size);
    }

    // zstd streaming read
    int i = 0;
    while (i < size - 1) {
        if (cf->handle.zstd.out_pos >= cf->handle.zstd.out_end) {
            if (zstd_refill(cf) <= 0) break;
        }
        char c = (char)cf->handle.zstd.out_buf[cf->handle.zstd.out_pos++];
        buf[i++] = c;
        if (c == '\n') break;
    }

    if (i == 0) return NULL;
    buf[i] = '\0';
    return buf;
}

// Close compressed file
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
```

### Success Criteria

#### Automated Verification:
- [ ] Header compiles without errors: `gcc -c -x c -fsyntax-only compressed_io.h`
- [ ] No memory leaks with valgrind test

#### Manual Verification:
- [ ] Code review for edge cases (NULL handling, error paths)

---

## Phase 3: Integrate Wrapper into tronko-assign

### Overview
Replace `gzFile`/`gzopen`/`gzgets`/`gzclose` calls with the new `CompressedFile` wrapper for query FASTA/FASTQ files.

### Changes Required

#### 1. Update tronko-assign.c
**File**: `tronko-assign/tronko-assign.c`

**Add include** (after line 3):
```c
#include "compressed_io.h"
```

**Single-end spec scanning** (lines ~1123-1128):
```c
// Change from:
gzFile reads_file = gzopen(opt.read1_file,"r");
find_specs_for_reads(read_specs, reads_file, opt.fastq);
gzclose(reads_file);

// To:
CompressedFile* reads_file = cf_open(opt.read1_file, "r");
find_specs_for_reads_cf(read_specs, reads_file, opt.fastq);
cf_close(reads_file);
```

**Single-end main processing** (lines ~1163, 1320):
```c
// Change from:
gzFile seqinfile = gzopen(opt.read1_file,"r");
// ... processing ...
gzclose(seqinfile);

// To:
CompressedFile* seqinfile = cf_open(opt.read1_file, "r");
// ... processing (update readInXNumberOfLines calls) ...
cf_close(seqinfile);
```

**Paired-end spec scanning** (lines ~1356-1367):
```c
// Change from:
gzFile seqinfile_1 = gzopen(opt.read1_file,"r");
gzFile seqinfile_2 = gzopen(opt.read2_file,"r");
// ...
gzclose(seqinfile_1);
gzclose(seqinfile_2);

// To:
CompressedFile* seqinfile_1 = cf_open(opt.read1_file, "r");
CompressedFile* seqinfile_2 = cf_open(opt.read2_file, "r");
// ...
cf_close(seqinfile_1);
cf_close(seqinfile_2);
```

**Paired-end main processing** (lines ~1423-1424, 1528-1529):
Same pattern as above.

#### 2. Update readreference.h
**File**: `tronko-assign/readreference.h`

**Add include and update function signatures**:
```c
#include "compressed_io.h"

// Add new function signatures that accept CompressedFile*
void find_specs_for_reads_cf(struct readSpecs *specs, CompressedFile *file, int fastq);
int readInXNumberOfLines_cf(CompressedFile *seqinfile, ...);
```

#### 3. Update readreference.c
**File**: `tronko-assign/readreference.c`

**Add wrapper functions** that call `cf_gets` instead of `gzgets`:

```c
// New function: find_specs_for_reads_cf
void find_specs_for_reads_cf(struct readSpecs *specs, CompressedFile *file, int fastq) {
    char buffer[FASTA_MAXLINE];
    specs->max_name_length = 0;
    specs->max_seq_length = 0;
    specs->count = 0;

    while(cf_gets(buffer, FASTA_MAXLINE, file) != NULL) {
        // Same logic as original, using cf_gets
        // ...
    }
}

// New function: readInXNumberOfLines_cf
int readInXNumberOfLines_cf(CompressedFile *seqinfile, ...) {
    // Same logic as original, using cf_gets instead of gzgets
    // ...
}
```

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make` compiles without errors
- [ ] Test with uncompressed FASTA: `./tronko-assign -f ref.txt -s -g test.fasta -o /tmp/out.txt`
- [ ] Test with gzip FASTA: `./tronko-assign -f ref.txt -s -g test.fasta.gz -o /tmp/out.txt`
- [ ] Test with zstd FASTA: `./tronko-assign -f ref.txt -s -g test.fasta.zst -o /tmp/out.txt`
- [ ] Output matches for all three formats

#### Manual Verification:
- [ ] Test with paired-end reads (both files compressed)
- [ ] Test with mixed compression (read1.fasta.zst, read2.fasta.gz)
- [ ] Verify memory usage stays constant during large file processing

---

## Phase 4: Update Build System and Documentation

### Overview
Ensure CI/CD builds work and documentation is updated.

### Changes Required

#### 1. Update GitHub Actions
**File**: `.github/workflows/build.yml`

Add zstd installation step:
```yaml
- name: Install dependencies
  run: |
    sudo apt-get update
    sudo apt-get install -y zlib1g-dev libzstd-dev
```

#### 2. Update README
**File**: `README.md`

Document supported input formats:
```markdown
## Supported Input Formats

tronko-assign accepts FASTA/FASTQ files in the following formats:
- Plain text (`.fasta`, `.fastq`)
- Gzip compressed (`.fasta.gz`, `.fastq.gz`)
- Zstandard compressed (`.fasta.zst`, `.fastq.zst`)

Compression format is auto-detected from file magic bytes, not file extension.
```

### Success Criteria

#### Automated Verification:
- [ ] GitHub Actions build passes
- [ ] `make` works on fresh Ubuntu with only documented dependencies

#### Manual Verification:
- [ ] README accurately describes the feature
- [ ] Build instructions are complete

---

## Testing Strategy

### Unit Tests

Create `test_compressed_io.c`:
- Test format detection with gzip, zstd, and plain files
- Test reading lines from each format
- Test error handling for corrupted files
- Test memory cleanup on close

### Integration Tests

1. **Single-end zstd**: Compress existing test FASTA, run assignment, compare output
2. **Paired-end zstd**: Both files compressed
3. **Mixed formats**: One gzip, one zstd
4. **Large file**: 1GB+ compressed file to verify streaming (no OOM)

### Manual Testing Steps

1. Create test files:
   ```bash
   zstd example.fasta -o example.fasta.zst
   gzip -k example.fasta
   ```

2. Run with each format and compare outputs:
   ```bash
   ./tronko-assign -f ref.txt -s -g example.fasta -o plain.txt
   ./tronko-assign -f ref.txt -s -g example.fasta.gz -o gzip.txt
   ./tronko-assign -f ref.txt -s -g example.fasta.zst -o zstd.txt
   diff plain.txt gzip.txt && diff plain.txt zstd.txt && echo "All match!"
   ```

3. Monitor memory usage:
   ```bash
   /usr/bin/time -v ./tronko-assign -f ref.txt -s -g large.fasta.zst -o out.txt
   # Check "Maximum resident set size" stays reasonable
   ```

## Performance Considerations

- Decompression adds minimal overhead (~5-10% CPU time)
- Memory usage: ~256 KB for zstd buffers + existing allocations
- I/O is typically the bottleneck, not decompression
- zstd decompression is faster than gzip in most cases

## References

- Research document: `thoughts/shared/research/2026-01-01-zstd-streaming-fasta-input.md`
- zstd zlibWrapper: https://github.com/facebook/zstd/tree/dev/zlibWrapper
- zstd streaming API: https://facebook.github.io/zstd/zstd_manual.html
