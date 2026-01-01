---
date: 2026-01-01T00:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "SIMD Vectorization Opportunities in tronko-assign"
tags: [research, codebase, simd, avx2, performance, optimization, vectorization]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: SIMD Vectorization Opportunities in tronko-assign

**Date**: 2026-01-01T00:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Can SIMD vectorization (AVX2/AVX-512) be applied to speed up the posterior score comparisons across nodes in tronko-assign? Specifically:
- Process 4-8 floats simultaneously instead of one at a time
- Identify hot loops and data structures suitable for vectorization
- Assess current autovectorization status and gaps

## Summary

**Yes, there are significant vectorization opportunities**, but they require different strategies:

1. **Low-hanging fruit**: The maximum-finding loop in `placement.c:873-888` iterates over contiguous float/double arrays and is immediately vectorizable with SIMD max operations.

2. **Medium effort**: The Cinterval range comparison in `placement.c:911-921` can be vectorized with SIMD comparison and mask operations.

3. **High effort**: The score accumulation in `getscore_Arr` (`assignment.c:154-207`) requires restructuring due to indirect indexing and character branching, but has the highest potential impact.

4. **Quick win**: Adding `-march=native` to the Makefile would enable autovectorization and allow the compiler to use AVX2/AVX-512 where applicable.

5. **Existing SIMD**: BWA's `ksw.c` already uses explicit SSE2 intrinsics; WFA2 has OpenMP pragmas but `-fopenmp` is not passed to the compiler.

## Detailed Findings

### 1. Current Compiler Configuration

**Location**: `tronko-assign/Makefile`

The current build uses minimal optimization flags:

```makefile
OPTIMIZATION = -O3
LIBS = -lm -pthread -lz -lrt -std=gnu99
```

**Missing flags for SIMD**:
- No `-march=native` or specific architecture targeting
- No `-mavx2` or `-mavx512f`
- No `-fopenmp` (WFA2 has OpenMP pragmas that are currently non-functional)
- No `-ftree-vectorize` hints (implicitly enabled at -O3 but limited without architecture flags)

**Recommendation**: Add `-march=native -mtune=native` to enable autovectorization for the host CPU.

### 2. Hot Loop #1: Maximum Score Finding (High Vectorization Potential)

**Location**: `tronko-assign/placement.c:873-888`

```c
type_of_PP maximum = -9999999999999999;
int minRoot = 0, minNode = 0, match_number = 0;

for (i = 0; i < number_of_matches; i++) {
    for (j = leaf_coordinates[i][0]; j < leaf_coordinates[i][0]+1; j++) {
        for (k = 0; k < 2*numspecArr[j]-1; k++) {
            if (maximum < nodeScores[i][j][k]) {
                maximum = nodeScores[i][j][k];
                match_number = i;
                minRoot = j;
                minNode = k;
            }
        }
    }
}
```

**Analysis**:
- Inner loop iterates over contiguous `nodeScores[i][j][k]` array (allocated at `allocateMemoryForResults.c:12-21`)
- `2*numspecArr[j]-1` nodes per tree (typically thousands)
- Simple max-finding pattern

**SIMD Strategy (AVX2 - 8 floats or 4 doubles)**:
```c
// Pseudocode for AVX2 vectorized max-finding
__m256 max_vec = _mm256_set1_ps(-FLT_MAX);
__m256i idx_vec = _mm256_set_epi32(7,6,5,4,3,2,1,0);
__m256i max_idx = idx_vec;

for (k = 0; k < num_nodes - 7; k += 8) {
    __m256 scores = _mm256_loadu_ps(&nodeScores[i][j][k]);
    __m256 cmp = _mm256_cmp_ps(scores, max_vec, _CMP_GT_OQ);
    max_vec = _mm256_blendv_ps(max_vec, scores, cmp);
    max_idx = _mm256_blendv_epi8(max_idx, idx_vec, (__m256i)cmp);
    idx_vec = _mm256_add_epi32(idx_vec, _mm256_set1_epi32(8));
}
// Horizontal reduction to find max and index
```

**Estimated Speedup**: 4-8x for the inner loop (depends on float vs double)

### 3. Hot Loop #2: Cinterval Range Comparison (Medium Vectorization Potential)

**Location**: `tronko-assign/placement.c:911-921`

```c
for (i = 0; i < number_of_matches; i++) {
    for (k = 0; k < 2*numspecArr[leaf_coordinates[i][0]]-1; k++) {
        if (nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) &&
            nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval)) {
            voteRoot[leaf_coordinates[i][0]][k] = 1;
            index++;
        }
    }
}
```

**Analysis**:
- Range comparison on contiguous array
- Binary output to `voteRoot` array

**SIMD Strategy (AVX2)**:
```c
// Pseudocode
__m256 lower = _mm256_set1_ps(maximum - Cinterval);
__m256 upper = _mm256_set1_ps(maximum + Cinterval);

for (k = 0; k < num_nodes - 7; k += 8) {
    __m256 scores = _mm256_loadu_ps(&nodeScores[i][tree][k]);
    __m256 in_range = _mm256_and_ps(
        _mm256_cmp_ps(scores, lower, _CMP_GE_OQ),
        _mm256_cmp_ps(scores, upper, _CMP_LE_OQ)
    );
    // Convert mask to integers and store to voteRoot
    int mask = _mm256_movemask_ps(in_range);
    // Process mask bits...
}
```

**Estimated Speedup**: 4-8x for the comparison loop

### 4. Hot Loop #3: Score Accumulation (Challenging but High Impact)

**Location**: `tronko-assign/assignment.c:154-207`

```c
for (i = 0; i < alength; i++) {
    if (positions[i] == -1) {
        score = score + log(0.01);
    } else {
        if (locQuery[i] == 'a' || locQuery[i] == 'A') {
            score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
        } else if (locQuery[i] == 'c' || locQuery[i] == 'C') {
            score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
        } else if (locQuery[i] == 'g' || locQuery[i] == 'G') {
            score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 2)];
        } else if (locQuery[i] == 't' || locQuery[i] == 'T') {
            score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 3)];
        }
        // ...
    }
}
```

**Vectorization Challenges**:
1. **Indirect indexing**: `positions[i]` creates non-contiguous memory access (gather pattern)
2. **Character branching**: `locQuery[i]` requires mapping A/C/G/T to indices 0/1/2/3
3. **Conditional logic**: Multiple code paths for missing data, gaps, etc.

**SIMD Strategy (Requires AVX2 Gather)**:

```c
// Step 1: Pre-compute nucleotide index array (0-3 for ACGT)
int8_t nuc_idx[alength];  // Precompute from locQuery

// Step 2: Use gather to load non-contiguous posterior values
for (i = 0; i < alength - 7; i += 8) {
    // Load 8 position indices
    __m256i pos_vec = _mm256_loadu_si256((__m256i*)&positions[i]);

    // Load 8 nucleotide indices
    // (need to compute: pos*4 + nuc_idx)
    __m256i full_idx = _mm256_add_epi32(
        _mm256_slli_epi32(pos_vec, 2),  // pos * 4
        _mm256_cvtepi8_epi32(_mm_loadl_epi64((__m128i*)&nuc_idx[i]))
    );

    // Gather 8 posterior values
    __m256 posteriors = _mm256_i32gather_ps(posteriornc, full_idx, sizeof(float));

    // Horizontal sum (or accumulate to vector and sum at end)
    score_vec = _mm256_add_ps(score_vec, posteriors);
}
```

**Prerequisites**:
- Precompute `nuc_idx[]` lookup table from `locQuery`
- Handle the `positions[i] == -1` case with masked operations
- Use `OPTIMIZE_MEMORY` build for float (8 floats/vector vs 4 doubles)

**Estimated Speedup**: 2-4x (limited by gather latency)

### 5. Data Structure Considerations

**Posterior Probability Storage** (`global.h:29`):
```c
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))
```

The 1D layout `[pos0_A, pos0_C, pos0_G, pos0_T, pos1_A, ...]` is actually good for SIMD gather operations because all 4 nucleotide values for a position are contiguous.

**type_of_PP Configuration** (`global.h:14-22`):
```c
#ifdef OPTIMIZE_MEMORY
    #define type_of_PP float   // 8 values per AVX2 register
#else
    #define type_of_PP double  // 4 values per AVX2 register
#endif
```

**Recommendation**: Build with `OPTIMIZE_MEMORY=1` for better SIMD utilization (8 floats vs 4 doubles per vector).

### 6. Existing SIMD Usage in Codebase

**BWA ksw.c** (SSE2 intrinsics):
- Location: `tronko-assign/bwa_source_files/ksw.c:149-212`
- Already uses `_mm_*` intrinsics for Smith-Waterman alignment
- Processes 16 x 8-bit scores simultaneously

**WFA2** (OpenMP pragmas):
- Location: `tronko-assign/WFA2/commons.h:256-260`
- Uses `PRAGMA_LOOP_VECTORIZE` for compiler hints
- Currently non-functional because `-fopenmp` not passed

### 7. Implementation Priority and Effort Matrix

| Hot Loop | File:Lines | Impact | Effort | Priority |
|----------|------------|--------|--------|----------|
| Maximum finding | placement.c:873-888 | Medium | Low | 1 |
| Cinterval comparison | placement.c:911-921 | Medium | Low | 2 |
| Score accumulation | assignment.c:154-207 | High | High | 3 |
| Compiler flags | Makefile | Medium | Trivial | 0 |

## Code References

- `tronko-assign/assignment.c:143-210` - `getscore_Arr` function (score accumulation)
- `tronko-assign/assignment.c:24-65` - `assignScores_Arr_paired` (recursive tree traversal)
- `tronko-assign/placement.c:862-891` - Maximum score finding loop
- `tronko-assign/placement.c:911-921` - Cinterval range comparison
- `tronko-assign/global.h:14-22` - `type_of_PP` definition (float vs double)
- `tronko-assign/global.h:29` - `PP_IDX` macro
- `tronko-assign/Makefile:13` - Optimization flags
- `tronko-assign/bwa_source_files/ksw.c:29-33` - SSE2 includes
- `tronko-assign/WFA2/commons.h:253-261` - Vectorization pragma macros

## Architecture Insights

1. **Memory Layout is SIMD-Friendly**: The `nodeScores[match][tree][node]` 3D array has contiguous inner dimension, enabling direct vectorized loads for the node loop.

2. **Posterior Storage is Gather-Ready**: The `PP_IDX(pos, nuc)` macro creates a layout suitable for AVX2 gather operations.

3. **Build Configuration Affects Vector Width**: `OPTIMIZE_MEMORY=1` uses float (8/register) vs double (4/register).

4. **Recursive Tree Traversal Limits Parallelism**: The `assignScores_Arr_paired` function traverses the tree recursively, which limits vectorization at that level. Vectorization should focus on the inner scoring loop.

## Recommended Next Steps

### Step 0: Enable Autovectorization (Immediate Win)
```makefile
# In tronko-assign/Makefile
OPTIMIZATION = -O3 -march=native -mtune=native -ffast-math
```

This enables:
- AVX2/AVX-512 instruction generation
- Compiler autovectorization
- Fast math for floating-point optimizations

### Step 1: Vectorize Maximum Finding
Create a SIMD version of the max-finding loop using AVX2 intrinsics. This is a well-understood pattern with predictable speedup.

### Step 2: Vectorize Cinterval Comparison
Apply similar SIMD pattern for range comparison and mask generation.

### Step 3: Restructure Score Accumulation
This requires:
- Pre-computing nucleotide index lookup table
- Implementing AVX2 gather-based scoring
- Handling edge cases (missing data, gaps)

### Step 4: Enable OpenMP in WFA2
```makefile
LIBS = -lm -pthread -lz -lrt -std=gnu99 -fopenmp
```

This would activate the existing OpenMP pragmas in WFA2.

## Open Questions

1. **Profiling Data**: What percentage of runtime is spent in each hot loop? Profiling with `perf` would help prioritize.

2. **CPU Target**: What minimum CPU should be supported? (SSE4.1, AVX2, AVX-512)

3. **Memory Bandwidth**: For large trees, are the loops compute-bound or memory-bound?

4. **Accuracy Impact**: Would `-ffast-math` affect result reproducibility?

5. **ARM Support**: BWA has SSE2-to-NEON translation. Would new SIMD code need NEON equivalents?
