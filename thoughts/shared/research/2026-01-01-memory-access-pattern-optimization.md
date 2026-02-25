---
date: 2026-01-01T00:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Memory Access Pattern Optimization Strategies for tronko-assign"
tags: [research, codebase, tronko-assign, performance, memory-access, cache-optimization, simd, mmap]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Memory Access Pattern Optimization Strategies for tronko-assign

**Date**: 2026-01-01
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Can we improve tronko-assign's algorithmic performance through memory access pattern optimizations, specifically:
1. Cache-oblivious algorithms for better memory locality
2. Structure-of-Arrays (SoA) vs Array-of-Structures (AoS) for better SIMD utilization
3. Memory-mapped I/O to let OS handle paging for huge databases

## Summary

This research identifies **significant optimization opportunities** in tronko-assign's memory access patterns:

| Optimization | Current State | Improvement Potential | Implementation Effort |
|--------------|---------------|----------------------|----------------------|
| **SoA for `posteriornc`** | Already 1D layout | N/A - done | Completed |
| **mmap for Binary DB** | fread/fseek | HIGH - eliminate copies | Medium |
| **mmap for BWA Index** | bulk fread | HIGH - 3 large files | Medium |
| **SIMD in scoring** | Scalar loops | MEDIUM - 4-8x potential | Medium-High |
| **Cache-oblivious tree** | Recursive DFS | LOW - not bottleneck | High |

**Key finding**: The `posteriornc` data structure (90% of memory) has already been optimized to 1D contiguous layout. The next highest-impact opportunities are **memory-mapped I/O** for database loading and **SIMD vectorization** of the scoring loop.

---

## Detailed Findings

### 1. Structure-of-Arrays Analysis

#### Current Data Layout Status

The codebase has **already been optimized** to use a 1D array layout for the critical `posteriornc` data:

**Node Structure** (`tronko-assign/global.h:85-93`):
```c
typedef struct node{
    int up[2];           // 8 bytes - child indices
    int down;            // 4 bytes - parent index
    int nd;              // 4 bytes - node descriptor
    int depth;           // 4 bytes - tree depth
    type_of_PP *posteriornc;  // 8 bytes - 1D array pointer
    char *name;          // 8 bytes - leaf name
    int taxIndex[2];     // 8 bytes - taxonomy indices
}node;
```

**1D Layout with Index Macro** (`global.h:29`):
```c
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))
// Access: posteriornc[PP_IDX(position, nucleotide)]
// Layout: [pos0_A, pos0_C, pos0_G, pos0_T, pos1_A, pos1_C, ...]
```

This is effectively **AoS at the node level, but SoA-like within each node's probability data** - the nucleotide values for each position are contiguous, enabling potential SIMD vectorization.

#### SoA Transformation Opportunity: Node Metadata

A **full SoA transformation** for the node structure could separate tree-navigation fields from scoring fields:

```c
// Theoretical SoA layout for tree navigation
struct TreeStructure {
    int *up0;           // [num_nodes]
    int *up1;           // [num_nodes]
    int *down;          // [num_nodes]
    int *depth;         // [num_nodes]
};

// Separate SoA for scoring (already contiguous per-node)
struct PosteriorData {
    type_of_PP *data;   // [num_nodes * numbase * 4]
};
```

**Assessment**: LOW PRIORITY - Tree traversal is not the bottleneck; scoring loop iteration dominates. The current layout is acceptable.

### 2. Cache-Oblivious Algorithm Opportunities

#### Current Hot Loops

**Pattern 1: Posterior Probability Transformation** (`tronko-assign.c:50-69`)
```c
for(i=0; i<numberOfRoots; i++){
    for (j=0; j<2*numspecArr[i]-1; j++){
        for (k=0; k<numbaseArr[i]; k++){
            for (l=0; l<4; l++){
                treeArr[i][j].posteriornc[PP_IDX(k, l)] = log( (f + g) );
            }
        }
    }
}
```
- **Access pattern**: Sequential within each node's `posteriornc` (good locality)
- **Cache behavior**: Already optimal - stride-1 access through contiguous memory
- **Optimization**: None needed for cache; SIMD opportunity exists

**Pattern 2: Scoring Loop** (`assignment.c:143-210`)
```c
for (i=0; i<alength; i++){
    if (locQuery[i]=='A'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
    } else if (locQuery[i]=='C'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
    }
    // ... G, T
}
```
- **Access pattern**: Indirect indexing via `positions[i]` (random access within node)
- **Cache behavior**: Moderate - positions may be non-sequential
- **Optimization opportunity**: Prefetching based on positions array

**Pattern 3: Recursive Tree Traversal** (`assignment.c:24-52`)
```c
void assignScores_Arr_paired(int rootNum, int node, ...) {
    int child0 = treeArr[rootNum][node].up[0];
    int child1 = treeArr[rootNum][node].up[1];
    // Recursive calls for each child
    assignScores_Arr_paired(rootNum, child0, ...);
    assignScores_Arr_paired(rootNum, child1, ...);
}
```
- **Access pattern**: Depth-first traversal (tree order, not memory order)
- **Cache behavior**: Poor - jumps between arbitrary memory locations

#### Cache-Oblivious Considerations

**Van Emde Boas Tree Layout**: For cache-oblivious tree traversal, nodes could be reordered in memory to match traversal patterns:

```
Traditional array layout:  [root, level1_a, level1_b, level2_a, level2_b, ...]
Van Emde Boas layout:      [root, subtree1_recursive, subtree2_recursive, ...]
```

**Assessment**: LOW PRIORITY for tronko-assign because:
1. The scoring loop (Pattern 2) dominates runtime, not tree traversal
2. Tree structure is fixed after loading; reordering would complicate indexing
3. Each query touches only nodes on the path from matched leaf to root (small subset)

### 3. Memory-Mapped I/O Analysis

#### Current I/O Patterns

| File Type | Function | Current Method | mmap Potential |
|-----------|----------|----------------|----------------|
| Text Reference DB | `readReferenceTree()` | gzgets + sscanf | None (gzip) |
| **Binary Reference DB** | `readReferenceBinary()` | fread + fseek | **HIGH** |
| Gzipped Binary DB | `readReferenceBinaryGzipped()` | gzread + gzseek | None (gzip) |
| **BWA BWT Index** | `bwt_restore_bwt()` | bulk fread | **HIGH** |
| **BWA SA Index** | `bwt_restore_sa()` | bulk fread | **HIGH** |
| **BWA PAC File** | `bwa_idx_load_from_disk()` | bulk fread | **HIGH** |
| FASTA/FASTQ Reads | `readInXNumberOfLines()` | gzgets batch | Low (streaming) |

#### High-Impact mmap Candidates

**1. Binary Reference Database** (`readreference.c:701-981`)

The uncompressed binary format (`.trkb`) is an excellent mmap candidate:
- **File structure**: Header (64 bytes) + Taxonomy section + Node section + Posterior section
- **Posterior section**: Contiguous `float` array per tree
- **Current loading**: Multiple `fseek` + `fread` calls
- **mmap benefit**: Zero-copy access, OS-managed paging, reduced syscalls

```c
// Current approach
fseek(fp, posterior_offset, SEEK_SET);
fread(treeArr[t][n].posteriornc, sizeof(float), count, fp);

// With mmap
void *mapped = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
treeArr[t][n].posteriornc = (float*)(mapped + posterior_offset + node_offset);
```

**Implementation considerations**:
- Posteriors can point directly into mapped region (no copy needed)
- Log transformation (`store_PPs_Arr`) still needed - could use copy-on-write or transform in-place
- Deallocation simplified: single `munmap` vs thousands of `free` calls

**2. BWA Index Files** (`bwa_source_files/bwt.c:443-462`)

BWA loads three large index files entirely into memory:
- `.bwt` file: BWT array (potentially gigabytes for large genomes)
- `.sa` file: Suffix array
- `.pac` file: Packed reference sequence

Current loading:
```c
// bwt.c:456 - bulk read of entire BWT
fread_fix(fp, bwt->bwt_size<<2, bwt->bwt);
```

**mmap benefits**:
- Lazy loading: Only pages actually accessed are loaded
- Reduced startup time for large indexes
- Shared memory: Multiple tronko-assign processes can share the same mapped pages
- Better memory accounting: Memory attributed to file, not process heap

### 4. SIMD Vectorization Opportunities

#### Existing SIMD in Codebase

**BWA ksw.c** - SSE2 Smith-Waterman (`bwa_source_files/ksw.c:29-374`):
```c
#ifdef __aarch64__
#include "sse2neon.h"
#else
#include <emmintrin.h>
#endif

// SIMD-accelerated SW cells
__m128i e, h, t, f = zero, max = zero;
h = _mm_adds_epu8(h, _mm_load_si128(S + j));
h = _mm_subs_epu8(h, shift);
```
- Processes 16 `uint8_t` or 8 `int16_t` values in parallel
- Has ARM NEON fallback via `sse2neon.h`

**WFA2** - Compiler-assisted vectorization (`WFA2/commons.h:255-261`):
```c
#if defined(__clang__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("GCC ivdep")
#endif
```

**WFA2 Packed Matching** (`WFA2/wavefront_extend.c:172-196`):
```c
// Compare 8 characters at once using 64-bit XOR
uint64_t cmp = *pattern_blocks ^ *text_blocks;
while (__builtin_expect(cmp==0,0)) {
    offset += 8;
    ++pattern_blocks; ++text_blocks;
    cmp = *pattern_blocks ^ *text_blocks;
}
const int equal_chars = __builtin_ctzl(cmp) / 8;
```

#### Unvectorized Scoring Loop

The main scoring function (`assignment.c:143-210`) is currently scalar:

```c
for (i=0; i<alength; i++){
    if (locQuery[i]=='A'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
    } else if (locQuery[i]=='C'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
    } else if (locQuery[i]=='G'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 2)];
    } else if (locQuery[i]=='T'){
        score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 3)];
    }
}
```

**SIMD Vectorization Strategy**:

1. **Pre-encode query as nucleotide indices**:
```c
// Convert 'A','C','G','T' -> 0,1,2,3
int8_t query_indices[alength];
for (i=0; i<alength; i++) {
    query_indices[i] = nucleotide_to_index[locQuery[i]];
}
```

2. **Gather posterior values using SIMD**:
```c
// With AVX2/AVX-512 gather instructions
__m256 posteriors = _mm256_i32gather_ps(
    node_posteriors,           // base pointer
    positions_with_offsets,    // indices: positions[i]*4 + query_indices[i]
    sizeof(float)
);
```

3. **Horizontal sum for score**:
```c
// Reduce 8 floats to single sum
score = _mm256_reduce_add_ps(posteriors);
```

**Estimated speedup**: 4-8x for scoring loop (depending on alignment length and gather efficiency)

**Challenges**:
- Gather instructions have moderate latency (~8-12 cycles on Skylake)
- Indirect indexing via `positions[]` limits pure SIMD benefit
- Branch elimination provides some benefit even without full SIMD

### 5. Needleman-Wunsch Alignment Matrix

The N-W alignment (`alignment.c:89-167`) is already structured for potential vectorization:

```c
for(seq_j = 0; seq_j < len_j; seq_j++) {
    for(seq_i = 0; seq_i < len_i; seq_i++) {
        match_scores[index] = MAX4(
            match_scores[index_upleft] + substitution_penalty,
            gap_a_scores[index_upleft] + substitution_penalty,
            gap_b_scores[index_upleft] + substitution_penalty,
            min);
        // ... similar for gap_a, gap_b
    }
}
```

**Anti-diagonal SIMD**: Standard DP can be vectorized by processing anti-diagonals:
- Each anti-diagonal's cells are independent
- Can compute 4-16 cells in parallel

**Alternative**: The WFA2 library already provides optimized alignment; consider using it for all alignments rather than N-W.

---

## Recommendations

### Priority 1: Memory-Mapped Binary Database (HIGH IMPACT, MEDIUM EFFORT)

**Implementation plan**:
1. Create `readReferenceBinaryMmap()` function in `readreference.c`
2. Use `mmap()` with `MAP_PRIVATE` for copy-on-write if log transform needed in-place
3. Modify `store_PPs_Arr()` to work with mapped memory
4. Add cleanup to use `munmap()` instead of per-node `free()`

**Expected benefits**:
- Eliminate memory copy during load (currently ~14-28 MB for example database)
- Faster startup (10-50x for large databases)
- Better multi-process memory sharing
- Simplified memory management

### Priority 2: Memory-Mapped BWA Index (HIGH IMPACT, MEDIUM EFFORT)

**Implementation plan**:
1. Modify `bwt_restore_bwt()` to use mmap
2. Modify `bwt_restore_sa()` to use mmap
3. Modify `bwa_idx_load_from_disk()` to use mmap for PAC

**Expected benefits**:
- Lazy loading of large index files
- Shared memory across processes
- Reduced peak memory for partial index access

### Priority 3: SIMD Scoring Loop (MEDIUM IMPACT, MEDIUM-HIGH EFFORT)

**Implementation plan**:
1. Pre-encode query sequences to nucleotide indices
2. Implement AVX2/SSE2 gather-based scoring
3. Fall back to scalar for non-aligned portions
4. Add compile-time detection for SIMD support

**Expected benefits**:
- 4-8x speedup in scoring (the dominant operation)
- Better utilization of modern CPU capabilities

### Priority 4: Compiler Hints (LOW EFFORT, LOW-MEDIUM IMPACT)

**Quick wins**:
```c
// Add to hot loops
PRAGMA_LOOP_VECTORIZE
for (i=0; i<alength; i++) { ... }

// Add restrict qualifiers
type_of_PP getscore_Arr(..., type_of_PP * restrict posteriornc, ...)

// Add prefetch hints
__builtin_prefetch(&posteriornc[PP_IDX(positions[i+4], 0)], 0, 3);
```

### Not Recommended

1. **Cache-oblivious tree layout**: High complexity, low benefit for tronko's access patterns
2. **Full SoA transformation**: Already have 1D layout where it matters most
3. **Partial database loading**: Algorithm requires access to all nodes (LCA computation)

---

## Code References

- `tronko-assign/global.h:29` - `PP_IDX` macro for 1D access
- `tronko-assign/global.h:85-93` - Node structure definition
- `tronko-assign/allocatetreememory.c:12-28` - Tree allocation (already uses `calloc` for 1D)
- `tronko-assign/readreference.c:701-981` - Binary format loading (mmap candidate)
- `tronko-assign/assignment.c:143-210` - Scoring loop (SIMD candidate)
- `tronko-assign/bwa_source_files/bwt.c:443-462` - BWA index loading (mmap candidate)
- `tronko-assign/bwa_source_files/ksw.c:29-374` - Existing SSE2 SIMD example
- `tronko-assign/WFA2/commons.h:255-261` - Vectorization pragma macros

## Historical Context

Prior research has already addressed some memory optimizations:
- `thoughts/shared/plans/2025-12-30-memory-optimizations.md` - float vs double optimization (implemented)
- `thoughts/shared/plans/2025-12-30-bulk-allocation-1d-layout.md` - 1D array layout (implemented)
- `thoughts/shared/research/2025-12-29-reference-database-loading.md` - mmap considered but not yet implemented

## Related Research

- `thoughts/shared/research/2025-12-29-rust-port-feasibility.md` - SIMD opportunities via block-aligner
- `thoughts/shared/research/2025-12-31-benchmark-binary-vs-text-format.md` - Binary format performance

## Open Questions

1. **mmap vs copy-on-write**: Should log transformation happen in-place (requiring COW pages) or copy to separate buffer?
2. **SIMD portability**: Target SSE2 (universal) or AVX2 (better performance, narrower compatibility)?
3. **Prefetch distance**: What's the optimal prefetch distance for the scoring loop given typical alignment lengths?
4. **BWA mmap sharing**: How to handle multiple tronko-assign processes sharing BWA index mappings?
