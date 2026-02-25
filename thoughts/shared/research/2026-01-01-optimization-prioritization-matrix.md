---
date: 2026-01-01T14:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Optimization Prioritization Matrix for tronko-assign"
tags: [research, codebase, prioritization, performance, optimization, planning]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Optimization Prioritization Matrix for tronko-assign

**Date**: 2026-01-01T14:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Based on today's comprehensive research into performance optimization strategies for tronko-assign, what is the optimal prioritization of efforts based on implementation effort versus expected performance gains?

## Summary

This document consolidates findings from four research efforts conducted on 2026-01-01 and presents a prioritized matrix of optimization opportunities. The analysis reveals that **algorithmic improvements (early termination, pruning) offer the best effort-to-reward ratio**, followed by memory-mapped I/O, SIMD vectorization, and finally GPU acceleration for maximum performance.

**Key finding**: A phased approach starting with low-effort algorithmic changes can achieve 3-15x speedup within a week, with cumulative gains potentially reaching 100-500x with full GPU implementation.

## Source Research Documents

This prioritization synthesizes findings from:

1. `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`
2. `thoughts/shared/research/2026-01-01-memory-access-pattern-optimization.md`
3. `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`
4. `thoughts/shared/research/2026-01-01-gpu-acceleration-tronko-assign.md`

---

## Quick Reference: Effort vs. Impact Matrix

```
                        | LOW EFFORT | MEDIUM EFFORT | HIGH EFFORT | VERY HIGH EFFORT
------------------------|------------|---------------|-------------|------------------
HIGHEST IMPACT (10-50x) |            |               | GPU Score   | Full GPU Pipeline
                        |            |               | Kernel      |
------------------------|------------|---------------|-------------|------------------
HIGH IMPACT (5-20x)     | Early Term | mmap Binary   | Pre-computed| WFA-GPU
                        | + Pruning  | DB            | Bounds      | Integration
                        |            | mmap BWA Idx  |             |
------------------------|------------|---------------|-------------|------------------
MEDIUM IMPACT (2-5x)    | Compiler   | Two-phase     | SIMD Score  |
                        | Flags      | Screening     | Accumulation|
                        | SIMD Max   | A*/IDA* Search|             |
                        | SIMD Range |               |             |
------------------------|------------|---------------|-------------|------------------
LOW IMPACT (<2x)        | OpenMP     |               |             |
                        | Enable     |               |             |
```

---

## Tier 1: Immediate Wins (Do First)

**Source**: `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`
**Source**: `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`

| # | Optimization | File(s) | Effort | Expected Gain | Risk |
|---|--------------|---------|--------|---------------|------|
| **1** | **Compiler flags** (`-march=native -mtune=native`) | Makefile | 5 min | 10-30% | None |
| **2** | **Early termination** (strike counter) | assignment.c:24-65 | 1-2 days | **2-5x** | Low* |
| **3** | **Subtree pruning** (skip bad branches) | assignment.c:24-65 | 1-2 days | **2-10x** | Low* |
| **4** | **Enable OpenMP** for WFA2 | Makefile | 5 min | 10-20% | None |

*With conservative thresholds (2x Cinterval) to ensure no false negatives

### Tier 1 Details

#### 1. Compiler Flags (Trivial)

From SIMD research: Current Makefile lacks architecture-specific flags.

```makefile
# Current
OPTIMIZATION = -O3

# Recommended
OPTIMIZATION = -O3 -march=native -mtune=native -ffast-math
```

This enables autovectorization and allows the compiler to use AVX2/AVX-512 where applicable.

#### 2-3. Early Termination + Subtree Pruning

From algorithm research: The "baseball heuristic" from pplacer can terminate after finding nodes that are clearly suboptimal.

**Implementation approach**:
- Modify `assignScores_Arr_paired()` to accept `best_score` and `strikes` parameters
- Add command-line options: `--strike-box`, `--max-strikes`
- If `node_score > best_score + 2*Cinterval`, skip entire subtree

**Expected impact**: pplacer reports 2-10x speedup with default parameters.

#### 4. Enable OpenMP

From SIMD research: WFA2 has OpenMP pragmas but `-fopenmp` is not passed to the compiler.

```makefile
LIBS = -lm -pthread -lz -lrt -std=gnu99 -fopenmp
```

### Tier 1 Combined Impact

**Potentially 3-15x speedup with minimal code changes**

---

## Tier 2: Medium-Term Gains (Next Priority)

**Source**: `thoughts/shared/research/2026-01-01-memory-access-pattern-optimization.md`
**Source**: `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`
**Source**: `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`

| # | Optimization | File(s) | Effort | Expected Gain | Risk |
|---|--------------|---------|--------|---------------|------|
| **5** | **mmap Binary DB** | readreference.c:701-981 | 3-5 days | **2-5x startup**, shared mem | Low |
| **6** | **mmap BWA Index** | bwa_source_files/bwt.c:443-462 | 3-5 days | **2-5x startup**, lazy load | Low |
| **7** | **SIMD max finding** | placement.c:873-888 | 2-3 days | 4-8x (inner loop) | None |
| **8** | **SIMD range comparison** | placement.c:911-921 | 2-3 days | 4-8x (inner loop) | None |
| **9** | **Two-phase screening** | assignment.c, placement.c | 3-5 days | **3-10x** | Low** |

**~99.9% accuracy vs 100% - needs validation

### Tier 2 Details

#### 5-6. Memory-Mapped I/O

From memory access research: The uncompressed binary format (`.trkb`) and BWA index files are excellent mmap candidates.

**Benefits**:
- Zero-copy access, OS-managed paging
- Reduced syscalls
- Shared memory across processes
- Lazy loading of large index files

**Implementation**:
```c
// Current approach
fseek(fp, posterior_offset, SEEK_SET);
fread(treeArr[t][n].posteriornc, sizeof(float), count, fp);

// With mmap
void *mapped = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
treeArr[t][n].posteriornc = (float*)(mapped + posterior_offset + node_offset);
```

#### 7-8. SIMD Vectorization (Hot Loops)

From SIMD research: Maximum-finding loop and Cinterval range comparison are immediately vectorizable.

**Maximum finding** (`placement.c:873-888`):
- Inner loop iterates over contiguous `nodeScores[i][j][k]` array
- Simple max-finding pattern maps perfectly to AVX2 `_mm256_max_ps`

**Cinterval comparison** (`placement.c:911-921`):
- Range comparison on contiguous array
- Binary output to `voteRoot` array
- Maps to SIMD comparison and mask operations

**Estimated speedup**: 4-8x for these inner loops

#### 9. Two-Phase Screening

From algorithm research (EPA-ng approach):

**Phase 1**: Compute approximate scores at strategic "checkpoint" nodes
**Phase 2**: Only score nodes in candidate trees (top 10%)

**Expected impact**: EPA-ng reports reducing candidate branches from thousands to <10.

### Tier 2 Combined Impact

**Additional 2-5x on top of Tier 1**

---

## Tier 3: Substantial Investment

**Source**: `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`
**Source**: `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md`

| # | Optimization | File(s) | Effort | Expected Gain | Risk |
|---|--------------|---------|--------|---------------|------|
| **10** | **SIMD score accumulation** (AVX2 gather) | assignment.c:143-210 | 1-2 weeks | 4-8x (scoring) | Medium |
| **11** | **Pre-computed bounds** (DB build change) | tronko-build, global.h | 1 week | **5-20x** | Low |
| **12** | **A*/IDA* search** | assignment.c | 1 week | 2-10x | Medium |

### Tier 3 Details

#### 10. SIMD Score Accumulation

From SIMD research: The score accumulation in `getscore_Arr` requires restructuring due to indirect indexing and character branching.

**Challenges**:
- Indirect indexing: `positions[i]` creates non-contiguous memory access (gather pattern)
- Character branching: `locQuery[i]` requires mapping A/C/G/T to indices 0/1/2/3

**Strategy**:
1. Pre-compute nucleotide index lookup table from `locQuery`
2. Implement AVX2 gather-based scoring
3. Handle edge cases (missing data, gaps)

**Estimated speedup**: 2-4x (limited by gather latency)

#### 11. Pre-computed Bounds

From algorithm research: Store `max_descendant_score` in node structure during database build.

**Implementation**:
- Add auxiliary fields to node structure in `global.h`
- Compute bounds during `tronko-build`
- Use bounds for aggressive branch-and-bound pruning

**Research finding**: Analytical upper bounds can eliminate 92-98% of tree space.

#### 12. A*/IDA* Search

From algorithm research: Replace recursive DFS with iterative best-first search.

**Benefits**:
- Process nodes in score order rather than tree order
- Early termination when f_score exceeds best + Cinterval
- IDA* variant uses O(d) memory vs O(n) for standard A*

---

## Tier 4: Major Architecture Changes

**Source**: `thoughts/shared/research/2026-01-01-gpu-acceleration-tronko-assign.md`

| # | Optimization | Effort | Expected Gain | Risk | Notes |
|---|--------------|--------|---------------|------|-------|
| **13** | **GPU score kernel** (CUDA) | 2-3 weeks | **10-50x** (scoring) | Medium | NVIDIA only |
| **14** | **GPU parallel reduction** (Thrust) | 1 week | 100x (max finding) | Low | Builds on #13 |
| **15** | **WFA-GPU integration** | 2-3 weeks | 5-20x (alignment) | Medium | Replaces WFA2 |
| **16** | **Full GPU pipeline** | 1-2 months | **5-20x overall** | High | Multi-GPU support |

### Tier 4 Details

#### 13-14. GPU Score Kernel + Parallel Reduction

From GPU research: Score calculation and maximum finding are ideal GPU workloads.

**Score kernel**:
- Each position lookup is independent
- Maps perfectly to GPU parallel reduction
- Launch config: Grid `(num_reads, ceil(num_nodes/256))`, Block `(256)` threads

**Parallel reduction**:
- Using Thrust for maximum finding across 17K+ trees
- `thrust::max_element` provides optimized parallel reduction

#### 15. WFA-GPU Integration

From GPU research: WFA-GPU is a CUDA implementation of Wavefront Alignment.

- Allocates central diagonals of wavefronts in shared memory
- Batch-based processing
- Could replace WFA2 in `placement.c`

#### 16. Full GPU Pipeline

From GPU research: Complete pipeline with async streams for transfer/compute overlap.

**Considerations**:
- Tree data (~80GB for large databases) may exceed GPU memory
- Multi-GPU support needed for very large databases
- Consider BarraCUDA for BWA step

### Tier 4 Combined Impact

**5-20x overall speedup (on top of earlier tiers)**

---

## Recommended Implementation Order

```
Week 1:  [Tier 1: #1-4] --------------------------> ~3-15x speedup
         |
Week 2-3:[Tier 2: #5-6, #7-8] -------------------> Additional 2-3x
         |
Week 4:  [Tier 2: #9] ---------------------------> Additional 2-3x
         |
         +-- CHECKPOINT: Benchmark & validate accuracy
         |
Week 5-6:[Tier 3: #10-11] -----------------------> Additional 2-5x
         |
         +-- DECISION POINT: Is GPU worth it?
         |
If yes:  [Tier 4: #13-16] -----------------------> 5-20x more
```

---

## Expected Cumulative Speedup

| After | Optimizations | Cumulative Speedup | Accuracy |
|-------|---------------|--------------------|----------|
| Tier 1 | Compiler + Algorithm pruning | **3-15x** | 100% |
| Tier 2 | mmap + SIMD + Two-phase | **10-50x** | ~99.9% |
| Tier 3 | Full SIMD + Pre-computed | **20-100x** | 100% |
| Tier 4 | GPU acceleration | **100-500x** | 100% |

---

## Key Trade-offs

| Factor | Conservative Path | Aggressive Path |
|--------|-------------------|-----------------|
| **Accuracy** | Tiers 1-2 with 100% accuracy | Tier 2 with two-phase (~99.9%) |
| **Portability** | Avoid GPU, use SIMD only | Full GPU (NVIDIA lock-in) |
| **Maintenance** | Small changes, low risk | Architecture overhaul |
| **Time to value** | 1-2 weeks for 10x | 2+ months for 100x+ |

---

## Code References

### Tier 1 Targets
- `tronko-assign/Makefile:13` - Optimization flags
- `tronko-assign/assignment.c:24-65` - Tree traversal (early termination)
- `tronko-assign/WFA2/commons.h:255-261` - OpenMP pragmas

### Tier 2 Targets
- `tronko-assign/readreference.c:701-981` - Binary format loading (mmap)
- `tronko-assign/bwa_source_files/bwt.c:443-462` - BWA index loading (mmap)
- `tronko-assign/placement.c:873-888` - Maximum score finding (SIMD)
- `tronko-assign/placement.c:911-921` - Cinterval range comparison (SIMD)

### Tier 3 Targets
- `tronko-assign/assignment.c:143-210` - Score accumulation (SIMD gather)
- `tronko-assign/global.h:85-93` - Node structure (pre-computed bounds)

### Tier 4 Targets
- New `kernels.cu` file for CUDA kernels
- `tronko-assign/WFA2/` - Replace with WFA-GPU

---

## Related Research

- `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md` - Detailed SIMD analysis
- `thoughts/shared/research/2026-01-01-memory-access-pattern-optimization.md` - Memory optimization details
- `thoughts/shared/research/2026-01-01-gpu-acceleration-tronko-assign.md` - GPU implementation details
- `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md` - Algorithm improvements

---

## Open Questions

1. **Benchmark baseline**: What is the current performance on standardized test datasets?
2. **Accuracy validation**: How to validate that pruning maintains 100% accuracy?
3. **GPU hardware**: What GPU hardware is available for deployment?
4. **Cross-platform**: Is ARM/NEON support needed for portability?
5. **Database compatibility**: Do optimizations require database rebuild?

---

## Recommendation

Start with **Tier 1** (compiler flags + algorithm pruning) for immediate 3-15x gains with minimal risk. This can be done in a week and establishes a benchmark baseline for evaluating further optimizations.
