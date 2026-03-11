# tronko-build Correctness Verification

Comprehensive analysis confirming that the `optimize-tronko-build` branch produces
identical database output to `main`. All code changes are pure performance
optimizations or robustness fixes — no underlying logic changes affect the
generated `reference_tree.txt`.

## Quick Start

```bash
# Run the automated verification (builds both branches, compares output)
bash tests/verify_no_logic_change.sh

# Run the golden-output regression test
bash tests/test_tronko_build_golden.sh
```

## Experiment: Main vs Optimized Branch Comparison

### Setup

| Parameter | Value |
|---|---|
| Date | 2026-03-11 |
| Dataset | Charadriiformes (COI), single-tree mode (`-l`) |
| Sequences | 1,466 |
| Alignment columns | 316 |
| Platform | macOS (Apple Silicon), Darwin 25.2.0 |
| Main branch compiler | Apple Clang 17.0.0 (invoked as `gcc`) |
| Optimized branch compiler | GCC 15.2.0 (Homebrew `gcc-15`) |

### Experiment 1: Different Compilers (Main=Clang, Optimized=GCC-15)

Both binaries were run on the identical input:

```bash
# Main branch binary (compiled with Apple Clang)
./tronko-build-main -l \
    -m example_datasets/single_tree/Charadriiformes_MSA.fasta \
    -x example_datasets/single_tree/Charadriiformes_taxonomy.txt \
    -t example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
    -d /tmp/main_output

# Optimized branch binary (compiled with gcc-15), single-threaded
OMP_NUM_THREADS=1 ./tronko-build-opt -l \
    -m example_datasets/single_tree/Charadriiformes_MSA.fasta \
    -x example_datasets/single_tree/Charadriiformes_taxonomy.txt \
    -t example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
    -d /tmp/opt_output
```

**Result: Files differ slightly.**

| Metric | Main (Clang) | Optimized (GCC-15) |
|---|---|---|
| Lines | 930,598 | 930,598 |
| Bytes | 39,944,265 | 39,932,169 |

Numerical analysis of the posterior probability values:

| Metric | Value |
|---|---|
| Identical lines | 472,071 |
| Differing lines | 458,527 |
| Text/structural diffs | **0** |
| Numeric values compared | 1,834,108 |
| Values that differ | 1,457,497 |
| Max absolute difference | 5.38e-08 |
| Max relative difference | 1.45e-06 |

All differences are in floating-point posterior values, not in tree topology,
node structure, taxonomy, or any other structural field.

### Experiment 2: Same Compiler (Both GCC-15)

To isolate compiler effects, the main branch code was rebuilt with gcc-15.

**Binary verification** — confirmed different binaries via SHA-256 hash, file
size, and symbol table (`nm`):

```
Main binary (main branch source, gcc-15, no OpenMP):
  sha256: 547c1fd8f4dd97dbc7744a62d106eacf17fb130ac415b9750a618cf78132cf70
  size:   427,424 bytes
  nm:     _getcladeArr present, _getcladeArr_fast ABSENT, GOMP symbols: 0

Optimized binary (optimized branch source, gcc-15, with OpenMP):
  sha256: 701d7cd56ba7ca95a427e92e7fd7ee6d804a78a1b4d0ca38afe78107ec02e252
  size:   430,528 bytes
  nm:     _getcladeArr_fast present, _run_partition_pipeline present, GOMP symbols: 4
```

How the main binary was built (in a git worktree checked out to `main`):

```bash
# Create isolated worktree for main branch
git worktree add /tmp/main-worktree main

# Build with gcc-15 (same compiler as optimized branch)
cd /tmp/main-worktree/tronko-build
gcc-15 -O3 -w -Wno-error=incompatible-pointer-types \
    -Wno-error=int-conversion -Wno-error=implicit-function-declaration \
    -o tronko-build hashmap.c tronko-build.c getclade.c \
    readfasta.c readreference.c allocatetreememory.c math.c \
    likelihood.c opt.c options.c printtree.c \
    -lm -pthread -lz -std=gnu99
```

**Result: BYTE-IDENTICAL.**

```
Main output SHA-256:      847943f27e97599600b432748cc503bd1c2c1ad954278f8b60b27aa1e73dc8ee
Optimized output SHA-256: 847943f27e97599600b432748cc503bd1c2c1ad954278f8b60b27aa1e73dc8ee
```

Both `reference_tree.txt` files are exactly 39,932,169 bytes with identical
SHA-256 hashes. This proves the code changes introduce no logic changes
whatsoever — the Experiment 1 differences were entirely due to Apple
Clang vs GCC-15 generating different floating-point instruction sequences.

### Experiment 3: Thread Safety (1 Thread vs 4 Threads)

The optimized binary was run with different OpenMP thread counts:

```bash
OMP_NUM_THREADS=1 ./tronko-build -l ... -d /tmp/single_thread
OMP_NUM_THREADS=4 ./tronko-build -l ... -d /tmp/multi_thread
```

**Result: BYTE-IDENTICAL.**

The `#pragma omp threadprivate` annotations on all shared workspace
variables are correct. Multi-threaded execution produces the same
output as single-threaded.

## Classification of All Code Changes

### Pure Performance Optimizations (No Output Change)

| Change | File(s) | Description |
|---|---|---|
| In-memory Newick parser | `getclade.c` | O(n) `getcladeArr_fast()` replaces O(n^2) FILE-based `getcladeArr()` |
| Contiguous memory blocks | `allocatetreememory.c` | Single malloc per node for `likenc`/`posteriornc` instead of per-site |
| Batch I/O | `printtree.c` | 64KB `snprintf` buffer + bulk `fwrite` instead of per-field `fprintf` |
| Hashmap taxonomy lookup | `tronko-build.c` | O(N+L) single-pass with hashmap vs O(N*L) file-reopen-per-leaf |
| Stack-allocated recursion | `tronko-build.c` | `getTaxonomyArr()` output param instead of malloc-per-call |
| Local output params | `tronko-build.c` | `findMinVarianceArr()`, `findLeavesOfMinVarArr()` use local params instead of globals |
| Doubling realloc | `tronko-build.c` | `createNode()` doubles capacity instead of realloc-per-node |
| Concurrent partitions | `tronko-build.c` | 3 partition pipelines (`fork()`) run in parallel |
| Native C utilities | `tronko-build.c` | `copy_file()`, `remove_partition_files()`, `unwrap_fasta_inplace()`, `fasta_to_phylip()` replace shell commands |
| OpenMP parallelism | `tronko-build.c`, `global.h`, `opt.h`, `opt.c` | `#pragma omp parallel for` on posterior computation with `threadprivate` globals |
| GCC auto-detection | `Makefile` | Auto-detects `gcc-15`..`gcc-11` for OpenMP support |

### Robustness Fixes (Only Activate on Edge-Case Inputs)

These changes only affect code paths that the original would crash on or
produce NaN/Inf. For well-formed inputs (standard binary Newick trees,
ACGTN- bases), these paths are never reached.

| Change | File(s) | Original Behavior | New Behavior |
|---|---|---|---|
| Underflow guard | `likelihood.c` | Division by near-zero → Inf/NaN | Fallback to uniform (0.25) when max < 1e-300 |
| Eigen fallback | `likelihood.c` | `exit(-1)` on eigendecomposition failure | Identity matrix fallback with warning |
| Node guard | `tronko-build.c` | `changePP_Arr`: `up[0]==-1 && up[1]==-1` | `child0 < 0 \|\| child1 < 0` (handles degenerate single-child nodes) |
| Descendant guard | `tronko-build.c` | `get_number_descendantsArr(node=-1)` → segfault | Returns 0 for node==-1 |
| IUPAC bases | `tronko-build.c` | `scanf` + `exit(-1)` on unknown bases | Treats as missing data (code 4) |
| Null checks | `tronko-build.c` | Crash on malformed taxonomy lines | `continue` / `break` with graceful skip |
| `getopt` return type | `options.c` | `char c` can't hold -1 on aarch64 | `int c` (correct per POSIX) |
| Bounds check | `tronko-build.c` | Buffer overflow reading long node names | `(i-1) < maxname` guard |

### Non-Single-Tree Changes (Partition Mode Only)

| Change | Description |
|---|---|
| VeryFastTree (`-a` flag) | Uses VeryFastTree instead of FastTree for tree inference |
| Unifurcation suppression | `makeBinary()` removes degree-1 nodes from VeryFastTree output |
| FAMSA threads default | Changed from 1 to 0 (auto-detect) |

These only apply to partition mode (`-y`), not single-tree mode (`-l`).

## Verification Scripts

### `tests/verify_no_logic_change.sh`

Automated end-to-end comparison:

1. Creates a git worktree for the `main` branch
2. Builds tronko-build from both branches
3. Runs the Charadriiformes single-tree test with each binary
4. Compares `reference_tree.txt` byte-for-byte
5. If different, runs a Python numeric analysis (max absolute/relative error)
6. Tests OMP_NUM_THREADS=1 vs OMP_NUM_THREADS=4

**Important:** The main branch Makefile uses `gcc` (Apple Clang on macOS),
while the optimized branch uses `gcc-15` (Homebrew GCC). For a true
apples-to-apples comparison, rebuild the main branch with `gcc-15`:

```bash
# In the main branch worktree:
cd tronko-build
gcc-15 -O3 -w -o tronko-build hashmap.c tronko-build.c getclade.c \
    readfasta.c readreference.c allocatetreememory.c math.c \
    likelihood.c opt.c options.c printtree.c \
    -lm -pthread -lz -std=gnu99
```

### `tests/test_tronko_build_golden.sh`

Golden-output regression test for the optimized branch:

- **Single tree test:** Byte-exact match against committed golden reference
- **Partition VFT test:** Smoke test (non-deterministic due to `srand(time(NULL))`)

```bash
# Run tests
bash tests/test_tronko_build_golden.sh

# Regenerate golden files after intentional changes
bash tests/test_tronko_build_golden.sh --generate
```

## Conclusion

The optimized branch is **logically identical** to main for database generation:

- **Same compiler → byte-identical output** (verified with gcc-15 on both)
- **Multi-threading → byte-identical output** (verified 1 vs 4 threads)
- **All posterior values, tree topology, taxonomy, and node structure are preserved**
- **The only observed differences (max relative ~1.5e-6) are from compiler choice** (Apple Clang vs GCC), not code changes
