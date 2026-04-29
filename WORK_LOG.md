# Tronko Fork — Work Log

All work done on the `optimize-tronko-build` branch, organized by area. This branch represents ~3 months of performance optimization, bug fixing, accuracy research, and infrastructure work on the tronko phylogenetic assignment system.

---

## Table of Contents

1. [Performance: tronko-build (6.3x speedup)](#1-performance-tronko-build)
2. [Performance: tronko-assign (28-35% speedup)](#2-performance-tronko-assign)
3. [Binary Database Format (.trkb) with zstd Compression](#3-binary-database-format)
4. [End-to-End Build Pipeline (build-tronko-db.sh)](#4-build-pipeline)
5. [Ablation System (ablate-tronko-db.sh)](#5-ablation-system)
6. [AncestralClust Optimizations](#6-ancestralclust-optimizations)
7. [Bug Fixes](#7-bug-fixes)
8. [Correctness Verification](#8-correctness-verification)
9. [Accuracy Improvements (tronko-assign)](#9-accuracy-improvements)
10. [Research & Analysis](#10-research-and-analysis)
11. [Tools & Infrastructure](#11-tools-and-infrastructure)
12. [Performance Summary](#12-performance-summary)

---

## 1. Performance: tronko-build

**Result: 6.3x speedup** on single-tree builds (23.1s -> 3.7s on Charadriiformes). All optimizations produce byte-identical output, verified against golden references.

### Optimizations

| # | Optimization | What Changed | Files |
|---|---|---|---|
| 1 | **O(n) in-memory Newick parser** | Replaced O(n^2) FILE-based `getcladeArr()` with single-pass in-memory `getcladeArr_fast()` using pointer arithmetic | `getclade.c`, `getclade.h` |
| 2 | **Contiguous memory allocation** | Single malloc per node for all `likenc`/`posteriornc` data instead of per-site-per-nucleotide allocations (millions of mallocs eliminated) | `allocatetreememory.c` |
| 3 | **Batch I/O** | 64KB `snprintf` buffer + bulk `fwrite` instead of per-field `fprintf` for `reference_tree.txt` output | `printtree.c` |
| 4 | **Hashmap taxonomy lookup** | O(1) hashmap lookup instead of O(N) re-opening taxonomy file for every leaf | `tronko-build.c` |
| 5 | **Stack-allocated recursion** | `getTaxonomyArr()` uses output param instead of malloc-per-recursive-call | `tronko-build.c` |
| 6 | **Doubling realloc strategy** | `createNode()` doubles tree capacity instead of realloc-per-node (O(n) -> O(log n) reallocs) | `tronko-build.c`, `global.h` |
| 7 | **FastTree fork/exec** | Direct `fork/execvp` instead of `system()` shell invocation | `tronko-build.c` |
| 8 | **Parallel partition pipelines** | FAMSA -> FastTree -> nw_reroot runs 3 partitions concurrently via `fork()` | `tronko-build.c` |
| 9 | **OpenMP parallel posteriors** | `#pragma omp parallel for` on outer tree loop for multi-partition databases. Each tree's ML optimization runs on a separate thread with `threadprivate` globals | `tronko-build.c`, `global.h`, `opt.c` |

### Benchmarks

| Test Case | Baseline | Optimized | Speedup |
|---|---|---|---|
| Single tree (Charadriiformes, 1,413 taxa, 317bp) | 23.1s | 3.7s | **6.3x** |
| Docker (Charadriiformes, gcc-12) | 10.33s | 4.50s | **2.3x** |
| Docker (CO1 Metazoa, 2,998 seqs, 362bp) | 35.79s | 11.60s | **3.1x** |

---

## 2. Performance: tronko-assign

**Result: 28-35% speedup** with byte-identical output verified via golden output testing across 6 configurations.

### Optimizations

| Optimization | Description |
|---|---|
| **Precomputed log constants** | Constants `log(0.25)`, `log(0.01)` computed once at startup instead of per-call |
| **256-entry nucleotide LUT** | Branchless nucleotide indexing in `getscore_Arr` hot loop |
| **Fast/slow path split** | Debug checks removed from hot loop; separate path when tracing is active |
| **WFA aligner reuse** | Single WFA aligner instance reused across all matches per thread (was create/destroy per match) |
| **Iterative tree DFS** | Explicit stack replaces recursive tree traversal (prevents stack overflow on deep trees) |
| **Output formatting** | `snprintf` cursor pattern replaces `asprintf`/`strcat`/`free` chains |
| **Bulk array clearing** | `memset` replaces per-element zeroing loops |

### Benchmarks (10K reads, single tree, 1 thread)

| Aligner | Before | After | Speedup |
|---|---|---|---|
| Needleman-Wunsch | 17.78s | 12.78s | **28.1%** |
| WFA | 14.54s | 9.48s | **34.8%** |

---

## 3. Binary Database Format

**Result: 30x compression** vs uncompressed binary, **4x smaller** than gzipped text.

Added `tronko-convert` tool and native binary format support (`.trkb`) with zstd streaming compression.

| Format | Size (16S_Bacteria, 17,868 trees) | Load Time |
|---|---|---|
| Text (gzipped) | 758 MB | — |
| Binary (uncompressed) | 5.8 GB | 3.14s |
| Binary (zstd compressed) | 192 MB | 3.77s |

Accuracy validation on 100K reads: 99.72% identical assignments (0.28% difference is normal multi-threaded floating-point variance).

### Files

- `tronko-convert/` — Format conversion tool with zstd support
- `tronko-convert/FORMAT_SPECIFICATION.md` — Binary format spec
- `tronko-assign/compressed_io.h` — Zstd streaming decompression

---

## 4. Build Pipeline

`build-tronko-db.sh` — Complete end-to-end pipeline for building CRUXv2-style tronko databases from raw FASTA + taxonomy.

### Pipeline Steps

```
Input FASTA + Taxonomy
    |
    v
1. AncestralClust clustering (sequence similarity-based)
    |
    v
2. Per-cluster FAMSA alignment
    |
    v
3. Per-cluster FastTree tree inference
    |
    v
4. Tree rerooting (nw_reroot)
    |
    v
5. tronko-build (partitioning + posterior computation)
    |
    v
6. marker.fasta concatenation (gap-free sequences)
    |
    v
7. BWA indexing
    |
    v
8. tronko-convert to .trkb binary format
```

### Features

- **Checkpoint caching** (`.cache/` directory) — resume on failure without recomputing
- Colon-to-underscore renaming in FASTA headers (FastTree compatibility)
- Configurable clustering parameters (cutoff, bin size)
- `final_partitions.txt` manifest for identifying leaf-level partitions
- Automatic deduplication handling for FastTree name parsing

---

## 5. Ablation System

`ablate-tronko-db.sh` — Remove sequences from a tronko database without a full rebuild. Skips alignment and tree building; only recomputes posteriors.

| Step | Full Rebuild (101K seqs) | Ablation |
|---|---|---|
| FAMSA + FastTree | ~30 min | Skipped |
| SP-score partitioning | ~15 min | Skipped |
| Posterior computation | ~5 min | ~5 min |
| BWA indexing | ~2 min | ~2 min |
| **Total** | **~50 min** | **~7 min** |

### How It Works

1. Master build exports subtrees via `-E` flag (one-time)
2. Ablation prunes specified accessions from each subtree's tree (via `nw_prune`), filters MSA and taxonomy
3. Posteriors recomputed from scratch on pruned trees
4. New `reference_tree.txt` + BWA index generated

Full documentation in `ABLATION.md` with Python API examples and batch ablation patterns.

---

## 6. AncestralClust Optimizations

Optimized fork of [lpipes/AncestralClust](https://github.com/lpipes/AncestralClust) for the clustering step:

| Optimization | Description |
|---|---|
| **qsort replaces bubble sort** | O(n log n) instead of O(n^2) for distance sorting |
| **Nucleotide lookup table** | Branchless base conversion |
| **Flat contiguous distance matrix** | Single allocation instead of pointer-per-row |
| **Balanced thread partitioning** | Even work distribution across threads |
| **O(N) tree diameter** | Replaces O(N^2) all-pairs computation |
| **WFA aligner reuse** | Single aligner instance across all pairwise alignments |
| **Direct CIGAR-to-JC distance** | Skips string materialization of alignment |
| **ML estimation malloc elimination** | Stack allocation in hot loop |

---

## 7. Bug Fixes

### Critical Fixes

| Bug | Impact | Fix | Commit |
|---|---|---|---|
| **tronko-assign segfault on ablated databases** | Crash when many sequences pruned | 3 fixes: NaN posterior guard, unifurcation suppression, uninitialized filename | `1ad46ba` |
| **marker.fasta MSA gap bug** | BWA can't align against gapped references; silent misassignment | Build pipeline now concatenates raw (gap-free) FASTAs instead of MSA-aligned ones | `2ad1af6` |
| **Infinite loop in changePP_parents_Arr** | Hang on degenerate trees | Converted recursive to iterative with `parent < 0` guard | `2ad1af6` |
| **ZSTD streaming buffer overflow** | All reads silently dropped from compressed FASTA | Preserve unconsumed input across refill calls | `71f6ec3` |
| **Buffer overflow in path buffers** | Silent memory corruption with long output paths | `char[200]` -> `BUFFER_SIZE` for all path buffers | `4aa29fc` |

### Robustness Improvements

| Fix | Description |
|---|---|
| **Likelihood underflow guard** | Flat 0.25 priors when max < 1e-300 (prevents NaN posteriors) |
| **Eigen convergence fallback** | Identity matrix fallback instead of `exit(-1)` on degenerate trees |
| **IUPAC ambiguity codes** | Treated as missing data instead of crashing |
| **Null taxonomy guards** | Graceful skip on malformed taxonomy lines |
| **`getopt` return type** | `int` not `char` (correct per POSIX, fixes aarch64) |
| **Bounds check in readSeqArr** | Guard against node name buffer overflow |
| **macOS build fix** | `_XOPEN_SOURCE` + `_DARWIN_C_SOURCE` for `ucontext.h`/`Dl_info` |

---

## 8. Correctness Verification

Rigorous verification that optimizations produce identical output to upstream.

### tronko-build

| Experiment | Result |
|---|---|
| Same compiler (both gcc-15), same input | **Byte-identical** output (SHA-256 match) |
| Multi-threaded (OMP 1 vs 4 threads) | **Byte-identical** output |
| Different compilers (Clang vs gcc-15) | Max relative diff 1.45e-06 (compiler FP difference only) |

### tronko-assign

- Golden output tests across **6 configurations** (NW/WFA x single-end/paired-end x single/multi-thread)
- All optimizations produce byte-identical output

### Scripts

- `tests/verify_no_logic_change.sh` — Automated branch comparison
- `tests/test_tronko_build_golden.sh` — Golden output regression tests
- `tests/CORRECTNESS_VERIFICATION.md` — Detailed analysis with SHA-256 hashes

---

## 9. Accuracy Improvements

### Implemented

| Feature | Description | Flag |
|---|---|---|
| **minimap2 aligner** | Alternative to BWA with better sensitivity at high divergence (>5%) | `--aligner minimap2` |
| **minimap2 tuning** | Configurable k-mer and window size for amplicon optimization | `--minimap2-kmer`, `--minimap2-window` |
| **Best-leaf override** | Skip LCA voting when a single leaf clearly dominates | `--best-leaf-threshold`, `--best-leaf-max-votes` |
| **Trace diagnostics** | Per-read scoring visibility for debugging assignments | `--trace-read READNAME` |
| **Max leaf matches** | Runtime-configurable leaf match cap (was compile-time) | `--max-leaf-matches` (alias: `--max-bwa-matches`) |
| **NUMCAT override** | Compile-time gamma rate category control | `make NUMCAT=4` |

### Removed (dead parameters)

Config sweep showed `<0.001 F1 range`, effectively no impact:
- `--consensus-threshold`
- `--soft-voting`
- `--vote-temperature`

### Planned / Researched

#### Leaf-Centered Voting & Alignment Quality Filtering

**Problem**: Internal tree nodes outscore leaves because their smoothed posteriors don't pay the full -5.7 log-penalty per mismatch site. This causes over-conservative assignments (genus when species is correct).

**Proposed solution**: Two-signal adaptive voting using `S_leaf_max` (best leaf score) and `delta_12` (gap to second-best leaf):

- **Case 1** (good leaf match, clear winner): Leaf-centered voting with narrow window -> species-level
- **Case 2** (good leaf match, multiple candidates): Leaf-centered, standard width -> genus-level
- **Case 3** (poor leaf match): Standard global-max voting -> lets internal nodes set the right level

New CLI parameters: `--leaf-quality-threshold`, `--narrow-factor`

Full analysis with worked examples in `thoughts/shared/plans/2026-03-11-leaf-centered-voting-and-alignment-filtering.md`.

#### "NA" Consensus Bug Fix

Two bugs in `tronko-assign.c:753-780` that partially cancel each other:

1. **Pointer comparison against "NA"** — compares `char*` addresses, never matches heap-allocated strings. "NA"="NA" incorrectly counts as taxonomic agreement.
2. **correctTax accumulates across levels** — not reset per taxonomic level, inflating consensus counts.

Combined effect: over-inflated consensus, accepting more specific levels than warranted.

Full analysis in `thoughts/shared/plans/2026-03-11-minimap2-integration-and-na-consensus-fix.md`.

#### Alignment Quality Filtering

- **MAPQ filtering**: BWA computes MAPQ but tronko-assign discards it (`%*d` in SAM parsing). Adding `--min-mapq` would filter ambiguous matches before scoring.
- **Post-alignment mismatch threshold**: Filter matches with >15% mismatch rate after WFA alignment.

---

## 10. Research & Analysis

### Major Research Documents

| Document | Topic |
|---|---|
| `thoughts/shared/research/2026-03-11-tronko-assign-amplicon-assignment-improvements.md` | Comprehensive pipeline analysis with prioritized improvements |
| `thoughts/shared/research/2026-03-11-ablation-setup-deep-dive.md` | Ablation system correctness and accuracy |
| `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md` | Phased optimization plan |
| `thoughts/shared/research/2026-01-03-non-determinism-mitigation-strategies.md` | ~3% multi-threading variance analysis and mitigation |
| `thoughts/shared/research/2026-01-01-algorithm-optimization-branch-bound-early-termination.md` | Branch-and-bound early termination design |
| `thoughts/shared/research/2026-01-01-bwa-multithreading-feasibility.md` | BWA threading analysis |
| `thoughts/shared/research/2026-01-01-gpu-acceleration-tronko-assign.md` | GPU feasibility study |
| `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md` | SIMD vectorization analysis |
| `thoughts/shared/research/2025-12-29-rust-port-feasibility.md` | Rust rewrite feasibility |

### Explored but Not Pursued

| Area | Reason |
|---|---|
| **Float precision (OPTIMIZE_MEMORY)** | 42% memory reduction but ~5% of reads lose ~2.5 taxonomic levels. Recommended: double for production |
| **Rust port** | Feasible but massive effort for C codebase with embedded BWA/WFA |
| **GPU acceleration** | Posterior scoring is memory-bound, not compute-bound. Limited benefit. |
| **SIMD vectorization** | Scoring loop has data-dependent branching that limits vectorization |
| **Full Bayesian assignment model** | Existing analysis shows simple thresholds work nearly as well |

---

## 11. Tools & Infrastructure

| Tool/Feature | Description |
|---|---|
| `tronko-convert` | Binary format converter (text <-> binary, with zstd compression) |
| `patch_taxonomy.py` | Fix taxonomy mismatches without full database rebuild |
| `sweep_tronko_build.py` | Parameter sweep across SP-score thresholds and min-leaf configs |
| `test_determinism.sh` | Statistical variance analysis for multi-threaded reproducibility |
| `confidence_analysis.py` / `confidence_analysis2.py` | Tested 16+ candidate confidence signals |
| Golden output test framework | 6 configurations, 10K benchmark reads |
| Dockerfile | Reproducible database builds |
| CodeQL security analysis | Automated C/C++ security scanning |
| Reproducibility documentation | 3 modes: single-threaded (deterministic), multi-threaded (~3% variance), consensus voting |

---

## 12. Performance Summary

| Component | Metric | Result |
|---|---|---|
| tronko-build (single tree) | Speedup | **6.3x** |
| tronko-build (2,998 seqs) | Speedup | **3.1x** |
| tronko-assign (NW) | Speedup | **28%** |
| tronko-assign (WFA) | Speedup | **35%** |
| Database size | Compression | **30x** vs uncompressed, **4x** vs gzipped text |
| Memory (float precision) | Reduction | **42%** (with accuracy tradeoff) |
| Ablation turnaround | Time savings | **7 min** vs **50 min** full rebuild |

All performance optimizations verified byte-identical to upstream output.

---

## Chronological Commit History

| Date | Commit | Summary |
|---|---|---|
| Dec 2025 | `1b47897` | Initial tronko-build optimizations + ancestralclust + build pipeline |
| Dec 2025 | `640f33e` | tronko-assign 28-35% speedup with golden output validation |
| Jan 2026 | `a3edf44` | AncestralClust optimization, FastTree switch, checkpoint caching |
| Jan 2026 | `00cc57e` | Binary reference tree format + memory optimizations |
| Jan 2026 | `d47ed01` | tronko-convert native zstd streaming |
| Jan 2026 | `71f6ec3` | Fix ZSTD streaming decompression buffer overflow |
| Jan 2026 | `e89bd07` | Reproducibility documentation + variance testing |
| Feb 2026 | `978cf2e` | Performance benchmarks + build-tronko-db.sh large dataset fixes |
| Mar 2026 | `2ad1af6` | Fix marker.fasta MSA gap bug + final_partitions.txt manifest |
| Mar 2026 | `1ad46ba` | Fix segfault on heavily ablated databases |
| Mar 2026 | `9e7fe1e` | Correctness verification docs + ablation research |
| Mar 2026 | `c0a6538` | Accuracy tuning parameters (trace, max-matches, consensus, soft voting) |
| Mar 2026 | `4aa29fc` | minimap2 aligner, best-leaf override, remove dead voting params |
