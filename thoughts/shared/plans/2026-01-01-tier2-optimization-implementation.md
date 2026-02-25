# Tier 2 Optimization Implementation Plan

## Overview

This plan implements the Tier 2 optimizations identified in the research prioritization matrix. All optimizations will be **feature-flagged** via command-line options to enable benchmarking and A/B testing of their effects on performance and accuracy.

## Current State Analysis

### Key Findings

1. **Binary database loading** (`readreference.c:701-981`): Uses sequential `fread()` calls to load posterior probabilities. Each node's posteriors are read individually with `fread(treeArr[t][n].posteriornc, sizeof(float), count, fp)`.

2. **BWA index loading** (`bwa_source_files/bwt.c:443-462`): The `bwt_restore_bwt()` function uses `fread_fix()` which reads in 16MB chunks but still copies data to malloc'd memory.

3. **Maximum score finding** (`placement.c:889-904`): Triple-nested loop iterating over matches, trees, and nodes to find maximum score. Pure sequential scan.

4. **Cinterval range comparison** (`placement.c:927-937`): Another triple-nested loop checking if scores fall within `[maximum-Cinterval, maximum+Cinterval]` range.

5. **Two-phase screening**: Not currently implemented. All nodes in matched trees are scored.

### Existing Infrastructure

- **Feature flags already exist** for Tier 1 optimizations: `--early-termination`, `--enable-pruning`, etc. in `options.c:42-49`
- **Options struct** in `global.h:190-238` can be extended with new Tier 2 flags
- **Binary format support** already exists in `readreference.c:701-981` via `readReferenceBinary()`

## Desired End State

After implementation:
1. New command-line flags: `--enable-mmap-db`, `--enable-mmap-bwa`, `--enable-simd`, `--enable-two-phase`
2. Memory-mapped loading for `.trkb` files when `--enable-mmap-db` is set
3. Memory-mapped loading for BWA index when `--enable-mmap-bwa` is set
4. SIMD-accelerated max-finding and range comparison when `--enable-simd` is set
5. Two-phase candidate screening when `--enable-two-phase` is set

### Verification

- All existing test cases pass with default settings
- Each feature can be enabled independently and in combination
- Benchmarks show measurable performance improvement
- No accuracy degradation (100% identical output files)

## What We're NOT Doing

- GPU acceleration (Tier 4)
- Pre-computed bounds in database (requires tronko-build changes)
- A*/IDA* search algorithm changes (Tier 3)
- Changes to the binary format specification
- Multi-threading changes beyond OpenMP

---

## Phase 1: Feature Flag Infrastructure

### Overview
Establish command-line options and data structures for all Tier 2 optimizations before implementing any.

### Changes Required

#### 1. Options Struct Extension
**File**: `tronko-assign/global.h`
**Lines**: ~232-238 (after existing Tier 1 flags)

Add new fields:
```c
// Tier 2 optimization toggles
int enable_mmap_db;         // Enable mmap for binary database (default: 0)
int enable_mmap_bwa;        // Enable mmap for BWA index (default: 0)
int enable_simd;            // Enable SIMD vectorization (default: 0)
                            // NOTE: Requires compile-time opt-in: make ENABLE_SIMD=1
int enable_two_phase;       // Enable two-phase screening (default: 0)
double two_phase_threshold; // Multiplier of Cinterval for candidate selection (default: 3.0)
                            // Higher = more conservative (keeps more candidates, better accuracy)
                            // Lower = more aggressive (fewer candidates, better performance)
```

#### 2. Command-Line Options
**File**: `tronko-assign/options.c`
**Lines**: ~49 (after existing long_options)

Add new options:
```c
{"enable-mmap-db",no_argument,0,0},
{"disable-mmap-db",no_argument,0,0},
{"enable-mmap-bwa",no_argument,0,0},
{"disable-mmap-bwa",no_argument,0,0},
{"enable-simd",no_argument,0,0},
{"disable-simd",no_argument,0,0},
{"enable-two-phase",no_argument,0,0},
{"disable-two-phase",no_argument,0,0},
{"two-phase-threshold",required_argument,0,0},
```

#### 3. Option Parsing
**File**: `tronko-assign/options.c`
**Lines**: ~152-158 (inside case 0 switch)

Add parsing logic for new options.

#### 4. Help Text Update
**File**: `tronko-assign/options.c`
**Lines**: ~92 (before closing of usage string)

Add documentation for new flags.

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make` compiles without errors
- [ ] `./tronko-assign --help` shows new options
- [ ] `./tronko-assign --enable-mmap-db --help` doesn't error

#### Manual Verification:
- [ ] Options appear in help text with descriptions
- [ ] Invalid option combinations produce helpful error messages

---

## Phase 2: Memory-Mapped Binary Database Loading

### Overview
Implement mmap-based loading for `.trkb` binary format files, providing zero-copy access to posterior probability data.

### Changes Required

#### 1. New mmap Loading Function
**File**: `tronko-assign/readreference.c`

Create new function `readReferenceBinaryMmap()` that:
1. Opens file and validates magic/version header
2. Uses `mmap()` to map entire file into memory
3. Sets up pointer arithmetic to access sections directly
4. For posteriors: points `treeArr[t][n].posteriornc` directly into mmap'd region

```c
/**
 * Read reference database from binary format using mmap
 * @param filename Path to .trkb file
 * @param name_specs Output array
 * @param mmap_handle Output: pointer to mmap'd region (caller must munmap)
 * @param mmap_size Output: size of mmap'd region
 * @return Number of trees loaded, or -1 on error
 */
int readReferenceBinaryMmap(const char *filename, int *name_specs,
                            void **mmap_handle, size_t *mmap_size);
```

#### 2. Global mmap Handle Storage
**File**: `tronko-assign/global.h`

Add globals to track mmap'd regions for cleanup:
```c
extern void *g_db_mmap_handle;
extern size_t g_db_mmap_size;
```

#### 3. Cleanup Function
**File**: `tronko-assign/readreference.c`

Add `cleanupMmapDb()` function to `munmap()` on program exit.

#### 4. Integration with Main
**File**: `tronko-assign/tronko-assign.c`

Modify database loading to check `opt.enable_mmap_db` flag and call appropriate function.

### Key Implementation Details

```c
// In readReferenceBinaryMmap:
int fd = open(filename, O_RDONLY);
struct stat st;
fstat(fd, &st);
void *mapped = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);

// Advise kernel about access pattern
madvise(mapped, st.st_size, MADV_SEQUENTIAL);

// For posterior access (example):
// Instead of allocating and reading, point directly:
size_t posterior_base = posterior_offset + tree_posterior_offset + node_offset;
treeArr[t][n].posteriornc = (float*)(((char*)mapped) + posterior_base);
```

### Caveats

- Cannot use mmap with gzipped `.trkb.gz` files (fall back to regular loading)
- Need to track which nodes have mmap'd posteriors vs allocated posteriors
- Add flag to node struct or use separate tracking array

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make` compiles without errors
- [ ] Test with small dataset: identical output with and without `--enable-mmap-db`
- [ ] No memory leaks detected with valgrind

#### Manual Verification:
- [ ] Memory usage is significantly lower (check with `top` or `/proc/self/status`)
- [ ] Startup time is faster for large databases
- [ ] Multiple processes can share the same mmap'd file

---

## Phase 3: Memory-Mapped BWA Index Loading

### Overview
Implement mmap-based loading for BWA index files (`.bwt` and `.sa`), enabling lazy loading and shared memory across processes.

### Design Decision: Wrapper Module Approach

**Rationale**: The BWA code in `bwa_source_files/` is third-party code from bwa-mem. To minimize risk and maintain the ability to easily update the upstream BWA code, we will:

1. **NOT modify** the original `bwt.c` or `bwt.h` files
2. **Create a new wrapper module** `bwt_mmap_wrapper.c` that provides mmap alternatives
3. **Use conditional compilation** to select between original and mmap implementations
4. **Allow easy rollback** by simply not using the wrapper

### Changes Required

#### 1. New Wrapper Module
**File**: `tronko-assign/bwt_mmap_wrapper.h` (new file)

```c
#ifndef _BWT_MMAP_WRAPPER_H_
#define _BWT_MMAP_WRAPPER_H_

#include "bwa_source_files/bwt.h"

// Extended bwt_t wrapper that tracks mmap state externally
// This avoids modifying the original bwt_t structure
typedef struct {
    bwt_t *bwt;           // Original BWA structure
    int is_mmaped;        // Flag: 1 if using mmap
    void *bwt_mmap;       // mmap handle for BWT data
    size_t bwt_mmap_size;
    void *sa_mmap;        // mmap handle for SA data
    size_t sa_mmap_size;
} bwt_mmap_t;

/**
 * Restore BWT from file using mmap (wrapper around original)
 * Returns a bwt_mmap_t that wraps the standard bwt_t
 */
bwt_mmap_t *bwt_restore_bwt_mmap_wrapper(const char *fn);

/**
 * Restore SA from file using mmap
 */
void bwt_restore_sa_mmap_wrapper(const char *fn, bwt_mmap_t *bwt_wrap);

/**
 * Destroy bwt_mmap_t, properly cleaning up mmap'd regions
 */
void bwt_mmap_destroy(bwt_mmap_t *bwt_wrap);

/**
 * Get the underlying bwt_t pointer for use with existing BWA functions
 */
static inline bwt_t *bwt_mmap_get_bwt(bwt_mmap_t *bwt_wrap) {
    return bwt_wrap ? bwt_wrap->bwt : NULL;
}

#endif
```

#### 2. Wrapper Implementation
**File**: `tronko-assign/bwt_mmap_wrapper.c` (new file)

```c
#include "bwt_mmap_wrapper.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

bwt_mmap_t *bwt_restore_bwt_mmap_wrapper(const char *fn) {
    bwt_mmap_t *wrap = calloc(1, sizeof(bwt_mmap_t));
    if (!wrap) return NULL;

    // Open file and get size
    int fd = open(fn, O_RDONLY);
    if (fd < 0) {
        free(wrap);
        return NULL;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        close(fd);
        free(wrap);
        return NULL;
    }

    // mmap the file
    void *mapped = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);  // Can close fd after mmap

    if (mapped == MAP_FAILED) {
        free(wrap);
        return NULL;
    }

    // Advise kernel about access pattern
    madvise(mapped, st.st_size, MADV_SEQUENTIAL);

    // Allocate bwt_t structure
    wrap->bwt = calloc(1, sizeof(bwt_t));
    if (!wrap->bwt) {
        munmap(mapped, st.st_size);
        free(wrap);
        return NULL;
    }

    // Parse header from mmap'd region (same format as original)
    // BWT file format: primary (8 bytes), L2[1-4] (32 bytes), then BWT data
    bwtint_t *header = (bwtint_t*)mapped;
    wrap->bwt->primary = header[0];
    memcpy(wrap->bwt->L2 + 1, header + 1, 4 * sizeof(bwtint_t));

    // Point bwt array directly into mmap'd region
    wrap->bwt->bwt_size = (st.st_size - sizeof(bwtint_t) * 5) >> 2;
    wrap->bwt->bwt = (uint32_t*)((char*)mapped + sizeof(bwtint_t) * 5);
    wrap->bwt->seq_len = wrap->bwt->L2[4];

    // Generate count table (this is computed, not stored)
    bwt_gen_cnt_table(wrap->bwt);

    // Store mmap info for cleanup
    wrap->is_mmaped = 1;
    wrap->bwt_mmap = mapped;
    wrap->bwt_mmap_size = st.st_size;

    return wrap;
}

void bwt_mmap_destroy(bwt_mmap_t *bwt_wrap) {
    if (!bwt_wrap) return;

    if (bwt_wrap->is_mmaped) {
        // Don't free bwt->bwt - it points into mmap'd region
        if (bwt_wrap->bwt_mmap) {
            munmap(bwt_wrap->bwt_mmap, bwt_wrap->bwt_mmap_size);
        }
        if (bwt_wrap->sa_mmap) {
            munmap(bwt_wrap->sa_mmap, bwt_wrap->sa_mmap_size);
        }
        // Free SA if it was allocated separately
        if (bwt_wrap->bwt && bwt_wrap->bwt->sa && !bwt_wrap->sa_mmap) {
            free(bwt_wrap->bwt->sa);
        }
        if (bwt_wrap->bwt) {
            free(bwt_wrap->bwt);
        }
    } else {
        // Use original destroy for non-mmap case
        bwt_destroy(bwt_wrap->bwt);
    }
    free(bwt_wrap);
}
```

#### 3. Integration in Main Code
**File**: `tronko-assign/tronko-assign.c` (or BWA calling code)

Use conditional logic to select implementation:
```c
#include "bwt_mmap_wrapper.h"

// Global to track which mode we're using
static bwt_mmap_t *g_bwt_mmap = NULL;
static bwt_t *g_bwt = NULL;

void load_bwa_index(const char *prefix, int use_mmap) {
    char fn[1024];
    snprintf(fn, sizeof(fn), "%s.bwt", prefix);

    if (use_mmap) {
        g_bwt_mmap = bwt_restore_bwt_mmap_wrapper(fn);
        if (!g_bwt_mmap) {
            LOG_WARN("mmap BWA loading failed, falling back to standard loading");
            g_bwt = bwt_restore_bwt(fn);
        }
    } else {
        g_bwt = bwt_restore_bwt(fn);
    }
}

bwt_t *get_bwt(void) {
    if (g_bwt_mmap) return bwt_mmap_get_bwt(g_bwt_mmap);
    return g_bwt;
}

void cleanup_bwa_index(void) {
    if (g_bwt_mmap) {
        bwt_mmap_destroy(g_bwt_mmap);
        g_bwt_mmap = NULL;
    }
    if (g_bwt) {
        bwt_destroy(g_bwt);
        g_bwt = NULL;
    }
}
```

#### 4. Makefile Updates
**File**: `tronko-assign/Makefile`

Add wrapper source conditionally:
```makefile
# BWA mmap wrapper (use: make ENABLE_BWA_MMAP=1)
ifdef ENABLE_BWA_MMAP
    BWA_MMAP_SOURCES = bwt_mmap_wrapper.c
    BWA_MMAP_FLAGS = -DENABLE_BWA_MMAP
else
    BWA_MMAP_SOURCES =
    BWA_MMAP_FLAGS =
endif

SOURCES = ... $(BWA_MMAP_SOURCES)
```

### Advantages of This Approach

1. **Zero modifications to upstream BWA code** - can update `bwa_source_files/` without conflicts
2. **Easy rollback** - just don't compile with `ENABLE_BWA_MMAP=1`
3. **Runtime fallback** - if mmap fails, automatically falls back to original implementation
4. **Clear separation** - mmap logic is isolated in its own module

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make` compiles without errors
- [ ] Test assignment produces identical results with/without `--enable-mmap-bwa`
- [ ] No memory leaks

#### Manual Verification:
- [ ] BWA alignment still works correctly
- [ ] Memory footprint reduced during BWA index load phase
- [ ] Index loading is faster (especially for repeat runs)

---

## Phase 4: SIMD Vectorized Max-Finding

### Overview
Implement AVX2/AVX-512 vectorized maximum finding for the score comparison loops in placement.c.

### Changes Required

#### 1. SIMD Detection Header
**File**: `tronko-assign/simd_utils.h` (new file)

```c
#ifndef _SIMD_UTILS_H_
#define _SIMD_UTILS_H_

#include <stdbool.h>

// Runtime SIMD capability detection
bool simd_has_avx2(void);
bool simd_has_avx512(void);

// SIMD-accelerated max finding
// Returns index of maximum element and stores max value
int simd_find_max_float(const float *arr, int n, float *max_value);

// SIMD-accelerated range comparison
// Sets output[i] = 1 if low <= arr[i] <= high, else 0
void simd_range_compare_float(const float *arr, int n,
                               float low, float high, int *output);

#endif
```

#### 2. SIMD Implementation
**File**: `tronko-assign/simd_utils.c` (new file)

Implement using AVX2 intrinsics:
```c
#include <immintrin.h>
#include "simd_utils.h"

int simd_find_max_float(const float *arr, int n, float *max_value) {
    if (n == 0) {
        *max_value = -INFINITY;
        return -1;
    }

    int max_idx = 0;
    float max_val = arr[0];

#ifdef __AVX2__
    if (simd_has_avx2() && n >= 8) {
        __m256 max_vec = _mm256_set1_ps(-INFINITY);
        // Process 8 elements at a time
        int i;
        for (i = 0; i + 7 < n; i += 8) {
            __m256 vals = _mm256_loadu_ps(&arr[i]);
            max_vec = _mm256_max_ps(max_vec, vals);
        }
        // Horizontal max reduction
        float temp[8];
        _mm256_storeu_ps(temp, max_vec);
        for (int j = 0; j < 8; j++) {
            if (temp[j] > max_val) {
                max_val = temp[j];
            }
        }
        // Find actual index (need second pass for index)
        for (int j = 0; j < n; j++) {
            if (arr[j] == max_val) {
                max_idx = j;
                break;
            }
        }
        // Handle remainder
        for (; i < n; i++) {
            if (arr[i] > max_val) {
                max_val = arr[i];
                max_idx = i;
            }
        }
    } else
#endif
    {
        // Scalar fallback
        for (int i = 1; i < n; i++) {
            if (arr[i] > max_val) {
                max_val = arr[i];
                max_idx = i;
            }
        }
    }

    *max_value = max_val;
    return max_idx;
}
```

#### 3. Integrate into Placement
**File**: `tronko-assign/placement.c`
**Lines**: ~889-904 (max finding loop)

Add conditional SIMD path:
```c
// After collecting all scores into a flat array:
if (enable_simd && simd_has_avx2()) {
    // Flatten nodeScores into contiguous array for SIMD
    // Call simd_find_max_float()
} else {
    // Existing triple-nested loop
}
```

#### 4. Range Comparison SIMD
**File**: `tronko-assign/placement.c`
**Lines**: ~927-937

```c
if (enable_simd && simd_has_avx2()) {
    simd_range_compare_float(scores_flat, total_nodes,
                              maximum - Cinterval, maximum + Cinterval,
                              voteRoot_flat);
} else {
    // Existing loop
}
```

#### 5. Makefile Updates
**File**: `tronko-assign/Makefile`

Add SIMD source and conditional AVX flags:
```makefile
# SIMD support (use: make ENABLE_SIMD=1)
ifdef ENABLE_SIMD
    SIMD_FLAGS = -mavx2 -mfma
    SIMD_SOURCES = simd_utils.c
else
    SIMD_FLAGS =
    SIMD_SOURCES =
endif

SOURCES = ... $(SIMD_SOURCES)
```

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make ENABLE_SIMD=1` compiles without errors
- [ ] Test produces identical results with/without `--enable-simd`
- [ ] SIMD path is taken when flag is enabled (verify via debug logging)

#### Manual Verification:
- [ ] Benchmark shows speedup on max-finding operations
- [ ] Works on machines without AVX2 (graceful fallback)

---

## Phase 5: Two-Phase Screening

### Overview
Implement a two-phase scoring approach where Phase 1 computes approximate scores at strategic "checkpoint" nodes to identify candidate trees, and Phase 2 only fully scores nodes in the top candidates.

### Design Decision: Accuracy-First Approach

**Rationale**: Per user requirements, we prioritize accuracy over performance. The two-phase screening should:

1. Use a **conservative default threshold** (3.0 * Cinterval) that keeps most candidates
2. Emit **warnings** when enabled, reminding users this is an approximation
3. Include **accuracy tracking** to log when results might have differed
4. Provide **easy tuning** via `--two-phase-threshold` for users who want more aggressive screening

### Changes Required

#### 1. Checkpoint Node Identification
**File**: `tronko-assign/assignment.h`

Add structure to track checkpoint nodes:
```c
typedef struct {
    int tree_id;
    int checkpoint_node_ids[16];  // Up to 16 checkpoints per tree
    int num_checkpoints;
    type_of_PP checkpoint_scores[16];
} TreeCheckpoints;
```

#### 2. Phase 1: Quick Screening
**File**: `tronko-assign/assignment.c`

Add function:
```c
/**
 * Compute approximate scores using checkpoint nodes only
 * @return Number of candidate trees (trees within threshold)
 */
int computeCheckpointScores(int rootNum, char *locQuery, int *positions,
                            int alength, TreeCheckpoints *checkpoints);
```

Checkpoints should be:
- Root node
- Nodes at depth 1, 2, 3 (if they exist)
- A sample of leaf nodes

#### 3. Phase 2: Full Scoring (existing)
Only call `assignScores_Arr_paired()` for trees that pass Phase 1 screening.

#### 4. Integration in Placement
**File**: `tronko-assign/placement.c`

Modify the tree scoring section:
```c
if (enable_two_phase) {
    // Emit warning about approximation
    static int two_phase_warned = 0;
    if (!two_phase_warned) {
        LOG_WARN("Two-phase screening enabled (threshold=%.1f*Cinterval). "
                 "This is an approximation that may affect accuracy.", two_phase_threshold);
        two_phase_warned = 1;
    }

    // Phase 1: Quick screening
    int *candidate_trees = malloc(numberOfTotalRoots * sizeof(int));
    type_of_PP *approx_scores = malloc(numberOfTotalRoots * sizeof(type_of_PP));
    int num_candidates = 0;
    type_of_PP best_approx = -9999999999999999;

    for (int t = 0; t < numberOfTotalRoots; t++) {
        approx_scores[t] = computeCheckpointScore(t, ...);
        if (approx_scores[t] > best_approx) {
            best_approx = approx_scores[t];
        }
    }

    // Select candidates (conservative: keep anything within threshold of best)
    type_of_PP cutoff = best_approx - (two_phase_threshold * Cinterval);
    for (int t = 0; t < numberOfTotalRoots; t++) {
        if (approx_scores[t] >= cutoff) {
            candidate_trees[num_candidates++] = t;
        }
    }

    LOG_DEBUG("Two-phase screening: %d/%d trees passed Phase 1 (%.1f%%)",
              num_candidates, numberOfTotalRoots,
              100.0 * num_candidates / numberOfTotalRoots);

    // Phase 2: Full scoring only on candidates
    for (int c = 0; c < num_candidates; c++) {
        int t = candidate_trees[c];
        assignScores_Arr_paired(t, ...);
    }

    free(candidate_trees);
    free(approx_scores);
} else {
    // Existing full scoring for all trees
}
```

#### 5. Accuracy Tracking and Validation
**File**: `tronko-assign/placement.c`

Add debug-mode validation that compares two-phase result to full result:
```c
#ifdef DEBUG_TWO_PHASE
    // In debug builds, also run full scoring and compare
    // This allows validation that two-phase isn't missing the best match
    int full_best_tree = run_full_scoring(...);
    if (full_best_tree != two_phase_best_tree) {
        LOG_WARN("Two-phase screening missed best tree! "
                 "Two-phase selected tree %d, full scoring selected tree %d",
                 two_phase_best_tree, full_best_tree);
    }
#endif
```

Add logging to track:
- Number of trees screened in Phase 1
- Number of candidates passed to Phase 2
- Percentage of trees eliminated
- (Debug mode) Whether the final assignment would have been different without screening

### Success Criteria

#### Automated Verification:
- [ ] `make clean && make` compiles without errors
- [ ] Test with known dataset produces same results (may need threshold tuning)
- [ ] Logging shows Phase 1/Phase 2 statistics

#### Manual Verification:
- [ ] Performance improvement measurable (should reduce scoring calls significantly)
- [ ] Accuracy validation: compare output with/without two-phase on test datasets
- [ ] Threshold parameter behaves as expected

---

## Testing Strategy

### Unit Tests

1. **SIMD functions**: Test `simd_find_max_float()` and `simd_range_compare_float()` with known inputs
2. **mmap functions**: Test loading small binary files with both methods, compare results
3. **Two-phase screening**: Test checkpoint selection logic

### Integration Tests

1. **Full pipeline tests** with each optimization enabled individually
2. **Combination tests** with multiple optimizations enabled
3. **Regression tests** comparing output files byte-for-byte

### Benchmark Tests

Create benchmark script:
```bash
#!/bin/bash
# benchmark_tier2.sh

DB="reference_tree.trkb"
READS="test_reads.fasta"

echo "=== Tier 2 Optimization Benchmarks ==="
echo ""

echo "Baseline (no optimizations):"
time ./tronko-assign -r -f $DB -a ref.fasta -s -g $READS -o /tmp/baseline.txt

echo ""
echo "With mmap-db:"
time ./tronko-assign -r -f $DB -a ref.fasta -s -g $READS -o /tmp/mmap.txt --enable-mmap-db

echo ""
echo "With SIMD (requires: make ENABLE_SIMD=1):"
time ./tronko-assign -r -f $DB -a ref.fasta -s -g $READS -o /tmp/simd.txt --enable-simd

echo ""
echo "With two-phase (conservative threshold=3.0):"
time ./tronko-assign -r -f $DB -a ref.fasta -s -g $READS -o /tmp/twophase.txt \
    --enable-two-phase --two-phase-threshold 3.0

echo ""
echo "With all Tier 2 (excluding two-phase for accuracy):"
time ./tronko-assign -r -f $DB -a ref.fasta -s -g $READS -o /tmp/all_safe.txt \
    --enable-mmap-db --enable-mmap-bwa --enable-simd

echo ""
echo "=== Verification ==="

# These MUST be identical (no accuracy impact)
echo -n "mmap-db: "
diff /tmp/baseline.txt /tmp/mmap.txt && echo "IDENTICAL" || echo "DIFFERS (BUG!)"

echo -n "SIMD: "
diff /tmp/baseline.txt /tmp/simd.txt && echo "IDENTICAL" || echo "DIFFERS (BUG!)"

echo -n "All safe optimizations: "
diff /tmp/baseline.txt /tmp/all_safe.txt && echo "IDENTICAL" || echo "DIFFERS (BUG!)"

# Two-phase accuracy check
echo ""
echo "=== Two-Phase Accuracy Check ==="
TOTAL=$(wc -l < /tmp/baseline.txt)
MATCHING=$(comm -12 <(sort /tmp/baseline.txt) <(sort /tmp/twophase.txt) | wc -l)
ACCURACY=$(echo "scale=4; $MATCHING / $TOTAL * 100" | bc)
echo "Two-phase accuracy: $MATCHING / $TOTAL ($ACCURACY%)"

if [ "$MATCHING" -eq "$TOTAL" ]; then
    echo "Two-phase: IDENTICAL (100% accuracy)"
else
    DIFF_COUNT=$((TOTAL - MATCHING))
    echo "Two-phase: $DIFF_COUNT reads differ from baseline"
    echo "Consider increasing --two-phase-threshold if accuracy is critical"
fi
```

### Manual Testing Steps

1. Run with `16S_Bacteria` test dataset
2. Compare memory usage with `htop` or `/proc/self/status`
3. Run multiple times to measure variance
4. Test on different architectures (Intel vs AMD, x86 vs ARM)

---

## Performance Considerations

### Memory
- mmap reduces peak memory by ~50% for posterior data (stays on disk until accessed)
- SIMD requires temporary buffers for flattened arrays
- Two-phase uses minimal additional memory (just candidate list)

### CPU
- SIMD requires AVX2 support (most modern CPUs)
- Fallback paths ensure correctness on older hardware
- Two-phase reduces total computation significantly

### I/O
- mmap trades memory for page faults on first access
- Sequential access pattern works well with kernel prefetching
- Multiple processes can share mmap'd files

---

## Migration Notes

- No database format changes required
- All optimizations are backward compatible (off by default)
- Existing scripts continue to work unchanged
- Users can gradually enable optimizations for their workflows

---

## References

- Original research: `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md`
- SIMD details: `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`
- Memory access patterns: `thoughts/shared/research/2026-01-01-memory-access-pattern-optimization.md`
- Algorithm research: `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`

---

## Implementation Order

1. **Phase 1**: Feature flag infrastructure (1 day)
2. **Phase 2**: mmap binary database (3 days)
3. **Phase 3**: mmap BWA index (2 days)
4. **Phase 4**: SIMD vectorization (3 days)
5. **Phase 5**: Two-phase screening (4 days)

Each phase should be implemented, tested, and merged independently to allow incremental validation.
