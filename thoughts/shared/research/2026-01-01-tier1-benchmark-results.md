# Tier 1 Optimizations Benchmark Results

**Date**: 2026-01-01
**Dataset**: single_tree (164 reads, 1 tree, ~2931 nodes)
**Platform**: Linux 6.17.9-arch1-1

## Executive Summary

| Optimization | Status | Speedup | Accuracy |
|-------------|--------|---------|----------|
| NATIVE_ARCH | Working | ~0% | 100% |
| FAST_MATH | Working | ~2% | 100% |
| Subtree Pruning | Working | ~0% | 100% |
| Early Termination | **BROKEN** | ~15% | 0% |

**Key Finding**: Early termination is broken and produces incorrect results. It should not be used until the implementation is fixed.

## Compile-Time Optimizations

All compile-time optimizations produce **byte-identical output** to baseline.

| Configuration | Avg Time (s) | Min | Max | Speedup |
|--------------|--------------|-----|-----|---------|
| baseline | 0.524 | 0.515 | 0.532 | 1.00x |
| NATIVE_ARCH | 0.526 | 0.522 | 0.530 | 1.00x |
| FAST_MATH | 0.513 | 0.507 | 0.518 | 1.02x |
| NATIVE_ARCH + FAST_MATH | 0.521 | 0.507 | 0.538 | 1.01x |

**Observations**:
- Compile-time flags show minimal benefit on this small dataset
- Differences are within margin of error (~2%)
- FAST_MATH provides slight improvement
- All outputs are identical to baseline (verified with diff)

## Runtime Optimizations

### Subtree Pruning (Working)

| Configuration | Avg Time (s) | Accuracy |
|--------------|--------------|----------|
| No runtime opts | 0.515 | 100% |
| --enable-pruning | 0.514 | 100% |
| --pruning-factor 1.5 | 0.529 | 100% |
| --pruning-factor 2.0 | 0.510 | 100% |
| --pruning-factor 2.5 | 0.511 | 100% |
| --pruning-factor 3.0 | 0.514 | 100% |

**Observations**:
- Pruning maintains 100% accuracy as expected
- Minimal speedup on small tree (tree is not deep enough for pruning to help significantly)
- All pruning factor values produce identical, correct results

### Early Termination (BROKEN - DO NOT USE)

| Configuration | Avg Time (s) | Accuracy |
|--------------|--------------|----------|
| --early-termination | 0.445 | 0% |
| --early-termination --strike-box 0.5 | 0.441 | 0% |
| --early-termination --strike-box 50.0 --max-strikes 100 | 0.443 | 0% |

**Critical Issue**: Early termination produces completely incorrect results regardless of parameter settings. All reads are assigned to the order level (Charadriiformes, node 16) instead of species level.

**Root Cause Analysis**:
The early termination implementation terminates traversal too aggressively. Even with very conservative settings (strike-box=50.0, max-strikes=100), it stops at high-level nodes instead of exploring to leaf nodes.

**Example Output Comparison**:
```
# Baseline (correct):
GU572157.1_0  Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge  -29.976211

# Early termination (incorrect):
GU572157.1_0  Eukaryota;Chordata;Aves;Charadriiformes  0.000000
```

## Recommendations

### For Production Use

```bash
# Recommended build (safe, slight improvement)
make NATIVE_ARCH=1 FAST_MATH=1

# DO NOT use early termination flags until fixed:
# --early-termination (BROKEN)
# --strike-box (BROKEN)
# --max-strikes (BROKEN)

# Pruning is safe but shows minimal benefit on small trees:
# --enable-pruning (safe, 100% accuracy)
```

### For Maximum Accuracy (Current Recommendation)
```bash
make NATIVE_ARCH=1 FAST_MATH=1
./tronko-assign ...  # No runtime optimization flags
```

### For Large Trees (When Pruning May Help)
```bash
make NATIVE_ARCH=1 FAST_MATH=1
./tronko-assign --enable-pruning --pruning-factor 2.0 ...
```

## Next Steps

1. **Fix Early Termination Bug**: The implementation needs debugging. Likely issues:
   - Strike threshold calculation relative to Cinterval
   - Best score tracking across the DFS traversal
   - Possible issue with how scores are compared (sign/magnitude)

2. **Test on Larger Datasets**: The single_tree dataset is too small to see significant optimization benefits. Test on 16S_Bacteria dataset.

3. **Profile Hotspots**: Use perf/valgrind to identify actual bottlenecks before further optimization.

## Large Dataset Testing (16S_Bacteria)

**Status**: BLOCKED - Critical bug prevents testing

Attempted to benchmark with 16S_Bacteria dataset (17,868 trees, ~1.27M nodes) but encountered a segmentation fault. See `thoughts/shared/bugs/2026-01-01-16S-bacteria-segfault.md` for details.

| Dataset | Trees | Nodes | Memory | Status |
|---------|-------|-------|--------|--------|
| single_tree | 1 | 2,931 | ~40 MB | Working |
| 16S_Bacteria | 17,868 | 1,265,572 | ~13 GB | **CRASH** |

The crash occurs during read processing after successfully loading the reference database. This is a separate issue from the early termination bug and needs investigation before large-scale benchmarks can proceed.

## Raw Data

Full benchmark results saved to: `tronko-assign/benchmark_results.tsv`

## Environment

- GCC version: (system default)
- CPU: (run `lscpu` for details)
- RAM: (run `free -h` for details)
- Iterations per test: 3
