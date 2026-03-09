# Tronko Performance Optimizations

## tronko-build

### Summary

The tronko-build pipeline was optimized for a **6.3x speedup** on single-tree builds (23.1s → 3.7s on the Charadriiformes example dataset). All optimizations produce byte-for-byte identical output, verified against golden reference files.

### Optimizations

#### 1. O(n) In-Memory Tree Parser (`getclade.c`)

**Before:** `getcladeArr()` used FILE-based I/O with `fgetc()`. For every internal node, `getnodenumbArr()` scanned ahead through the Newick file using `fgetc`/`fsetpos` to count commas and determine node numbering — making the parser O(n^2) in the number of nodes.

**After:** `getcladeArr_fast()` reads the entire Newick file into memory once, then parses it in a single O(n) pass using pointer arithmetic. Node numbering is computed by capturing the global comma count at the correct point (after the separator between children), matching the original `getnodenumbArr` formula without lookahead.

**Files:** `getclade.c`, `getclade.h`, `tronko-build.c` (4 call sites updated)

#### 2. Contiguous Memory Allocation (`allocatetreememory.c`)

**Before:** Each node's `likenc` and `posteriornc` arrays were allocated as separate small mallocs — two 32-byte allocations per site per node, totaling millions of allocations for large trees.

**After:** Each node gets a single contiguous block for all its `likenc` data and another for `posteriornc`. Row pointers index into the block, preserving the existing `[site][nucleotide]` access pattern while dramatically reducing malloc overhead and improving cache locality.

**Files:** `allocatetreememory.c`

#### 3. Batch I/O for Reference Tree Output (`printtree.c`)

**Before:** `printTreeFile()` called `fprintf("%.17g")` individually for every posterior probability value — 4 calls per site per node, millions of calls total for large trees.

**After:** Lines are formatted into a 64KB buffer using `snprintf`, then flushed to disk with `fwrite` when the buffer is near-full. Combined with the existing 256KB `setvbuf`, this reduces per-call overhead by orders of magnitude.

**Files:** `printtree.c`

#### 4. VeryFastTree fork/exec (`tronko-build.c`)

**Before:** VeryFastTree was invoked via `system()`, which spawns a shell to parse the command string and handle stdout redirection.

**After:** Direct `fork/execvp` with manual `dup2` for stdout redirection. Eliminates shell startup overhead and is cleaner for process management in the parallel pipeline.

**Files:** `tronko-build.c` (`run_partition_pipeline()`)

#### 5. VeryFastTree realloc Bug Fix (`tronko-build.c`, `global.h`)

**Before:** `createNode()` called `realloc(m->tree, numNodes * sizeof(node))` on every single node creation during `parseNewick()` and `resolvePolytomy()`. Each realloc could move the entire tree array, invalidating any existing pointers into it and corrupting leaf names.

**After:** Added a `treeCapacity` field to `masterArr`. `parseNewick()` pre-allocates capacity estimated from the Newick string (3x the number of commas + parentheses). `createNode()` only reallocs when capacity is exceeded, using a doubling strategy. This reduces realloc frequency from O(n) to O(log n) and eliminates pointer invalidation during tree construction.

**Files:** `global.h`, `tronko-build.c`

#### 6. Parallel Partition Pipelines (`tronko-build.c`)

The 3-partition external tool pipeline (FAMSA → VeryFastTree/RAxML → nw_reroot) now runs all 3 partitions in parallel using `fork()`, with thread counts auto-divided across available cores. Previously these ran sequentially.

**Files:** `tronko-build.c` (`run_partition_pipeline()`, `createNewRoots()`)

### Benchmark

| Test Case | Baseline | Optimized | Speedup |
|-----------|----------|-----------|---------|
| Single tree (Charadriiformes, 1413 taxa, 317bp) | 23.1s | 3.7s | 6.3x |

Measured on Apple Silicon (M-series), `gcc-15 -O3 -fopenmp`.

### Testing

Golden output regression tests are in `tests/test_tronko_build_golden.sh`:

```bash
# Run tests
bash tests/test_tronko_build_golden.sh

# Regenerate golden files (after intentional output changes)
bash tests/test_tronko_build_golden.sh --generate
```

- **Single tree test:** Exact byte-for-byte match against golden reference
- **Partition VFT test:** Smoke test (checks output exists with >100 lines) because `srand(time(NULL))` in polytomy resolution makes output non-deterministic

## tronko-assign

See commit `640f33e` for tronko-assign optimizations (28-35% speedup with golden output validation).
