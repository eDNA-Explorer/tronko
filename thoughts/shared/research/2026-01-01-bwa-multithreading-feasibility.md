---
date: 2026-01-01T22:30:00-08:00
researcher: Claude
git_commit: 8fa99b63db2052b40a62e61b9726a6bb6960c389
branch: experimental
repository: tronko
topic: "Multi-threaded BWA Alignment Feasibility Assessment"
tags: [research, codebase, bwa, threading, performance, optimization]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Multi-threaded BWA Alignment Feasibility Assessment

**Date**: 2026-01-01T22:30:00-08:00
**Researcher**: Claude
**Git Commit**: 8fa99b63db2052b40a62e61b9726a6bb6960c389
**Branch**: experimental
**Repository**: tronko

## Research Question

How feasible is enabling multi-threaded BWA alignment in tronko-assign? Currently the thread count is hardcoded to 1 at `fastmap.c:697`. What would need to change, and what performance gains can we expect?

## Summary

**Feasibility: MEDIUM** - The threading infrastructure already exists and is functional, but enabling it requires fixing thread-safety issues in the tronko-specific code that was added on top of BWA.

**Expected Impact: LOW (10-20%)** - According to the optimization prioritization matrix, this is categorized as a low-impact optimization. Other optimizations (early termination, mmap, SIMD) offer significantly better returns.

**Key Blockers**:
1. Broken mutex pattern in `bwamem.c:877-911` (creates local mutex per call)
2. Global query matrix modifications during alignment
3. Result collection in Step 2 writes to shared `bwaMatches` array

**Recommendation**: This is a valid optimization but should be deprioritized in favor of Tier 1-2 optimizations that offer 2-10x gains. If pursued, it can be done in ~2-3 days with careful attention to the thread-safety issues.

---

## Detailed Findings

### The Hardcoded Constraint

**Location**: `tronko-assign/bwa_source_files/fastmap.c:696-697`

```c
//	if (opt->n_threads < 1) opt->n_threads = 1;
	opt->n_threads=1;
```

The original conditional check is commented out and replaced with a hardcoded `n_threads=1`.

Additionally, at line 842, the pipeline is hardcoded to single-threaded:
```c
//kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);
kt_pipeline(1, process, &aux, 3);
```

### Existing Threading Infrastructure

The good news: **BWA's threading infrastructure is complete and functional**. The implementation includes:

#### 1. kt_pipeline - Pipeline Parallelism (`kthread.c:63-151`)

A 3-stage producer-consumer pipeline:
- **Step 0**: Read sequences from query matrices → `ktp_data_t`
- **Step 1**: BWA alignment via `mem_process_seqs()`
- **Step 2**: Parse SAM output, populate `bwaMatches` results

Uses `pthread_mutex_t` and `pthread_cond_t` for ordered stage execution.

#### 2. kt_for - Parallel For Loop (`kthread.c:6-61`)

Work-stealing parallel-for using `__sync_fetch_and_add` atomics. Used inside `mem_process_seqs()` to parallelize alignment across sequences.

#### 3. Thread-Local Buffers (`bwamem.c:1302-1309`)

BWA properly allocates per-thread `smem_aux_t` buffers:
```c
w.aux = malloc(opt->n_threads * sizeof(smem_aux_t));
for (i = 0; i < opt->n_threads; ++i)
    w.aux[i] = smem_aux_init();
kt_for(opt->n_threads, worker1, &w, ...);
```

#### 4. Worker Functions (`bwamem.c:1251-1286`)

`worker1()` and `worker2()` access disjoint sequence indices, avoiding data races on the sequence arrays themselves.

---

### Thread Safety Issues (Why It's Disabled)

#### Issue 1: Broken Mutex Pattern (`bwamem.c:877-911`)

The code creates a **local mutex for each operation**, providing no cross-thread protection:

```c
// This is repeated 3 times in the function
if ( p->is_rev == 1 && paired == 0){
    pthread_mutex_t lock;                    // LOCAL mutex!
    pthread_mutex_init(&lock, NULL);
    pthread_mutex_lock(&lock);
    strcpy(seq_copy,singleQueryMat->queryMat[startline+seq_index]);
    reverse_complement(seq_copy);
    strcpy(singleQueryMat->queryMat[startline+seq_index],seq_copy);
    pthread_mutex_unlock(&lock);
}
```

**Problem**: Each thread creates its own mutex instance. Different threads could simultaneously write to overlapping indices of `singleQueryMat->queryMat` or `pairedQueryMat->query1Mat/query2Mat`.

**Fix Required**: Either:
- Use a global mutex for query matrix access, OR
- Make query matrices read-only (copy-on-write for reverse complement), OR
- Ensure each thread only accesses its assigned indices (currently true if `startline+seq_index` ranges don't overlap)

#### Issue 2: Global Query Matrix Modifications

The following global arrays are **modified** during alignment:

| Global | Location | Modification |
|--------|----------|--------------|
| `singleQueryMat->queryMat` | `bwamem.c:880` | Reverse complement in-place |
| `pairedQueryMat->query1Mat` | `bwamem.c:893` | Reverse complement in-place |
| `pairedQueryMat->query2Mat` | `bwamem.c:907` | Reverse complement in-place |

These same arrays are **read** in Step 0 (`bwa.c:75-93`) and Step 2 (`fastmap.c:187-536`).

#### Issue 3: Step 2 Result Collection (`fastmap.c:105-553`)

Step 2 writes to `aux->results[j]` based on parsed SAM output. With multiple pipeline workers:
- Multiple workers could be in Step 2 simultaneously
- The hashmap is created locally per Step 2 invocation (safe)
- But `aux->results` is shared across all workers

The index `j` is computed from sequence matching logic and could overlap if workers process interleaved batches.

#### Issue 4: global_bns Static Variable (`bwamem.c:45`)

```c
static const bntseq_t *global_bns; // for debugging only
```

Set at line 1297 before `kt_for`. With single-threaded execution this is safe. With multiple `kt_for` invocations from different threads, this would race.

---

### What Would Need to Change

#### Minimal Fix (Estimated: 2-3 days)

1. **Replace local mutexes with global mutex** in `bwamem.c:877-911`:
   ```c
   // At file scope
   static pthread_mutex_t g_query_matrix_mutex = PTHREAD_MUTEX_INITIALIZER;

   // In function
   pthread_mutex_lock(&g_query_matrix_mutex);
   // ... modify query matrix ...
   pthread_mutex_unlock(&g_query_matrix_mutex);
   ```

2. **Verify index isolation** - Confirm that each thread's `startline+seq_index` range is disjoint (should be true by design).

3. **Add thread parameter to run_bwa()** (`tronko-assign.c:173`):
   ```c
   void run_bwa(int start, int end, bwaMatches* bwa_results, int concordant,
                int numberOfTrees, char *databasefile, int paired,
                int max_query_length, int max_readname_length, int max_acc_name,
                int num_threads) {  // New parameter
       main_mem(databasefile, end-start, num_threads, ...);
   }
   ```

4. **Remove hardcoding** in `fastmap.c:697`:
   ```c
   if (opt->n_threads < 1) opt->n_threads = 1;  // Restore original check
   ```

5. **Enable pipeline threading** in `fastmap.c:842`:
   ```c
   kt_pipeline(no_mt_io ? 1 : 2, process, &aux, 3);  // Restore original
   ```

#### Cleaner Fix (Estimated: 3-5 days)

Instead of adding mutexes, make query matrices read-only:

1. **Pre-compute reverse complements** before BWA:
   - Create `singleQueryMat->queryMatRC` for all reverse complement sequences
   - Pass flag to indicate which version to use

2. **Remove in-place modifications** in `bwamem.c:877-911`

3. Benefits:
   - No mutexes needed
   - Better cache behavior (no false sharing)
   - Cleaner code

---

### Expected Performance Impact

From the optimization prioritization matrix (`thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md`):

| Optimization | Effort | Expected Gain |
|--------------|--------|---------------|
| **BWA Threading (OpenMP Enable)** | Low | **10-20%** |
| Early Termination + Pruning | Low-Medium | **2-10x** |
| mmap Binary DB | Medium | **2-5x startup** |
| SIMD Max Finding | Medium | **4-8x (inner loop)** |
| Two-Phase Screening | Medium | **3-10x** |

**Key insight**: BWA threading is the **lowest impact** optimization in the matrix. The actual alignment step is already relatively fast compared to:
- Database loading
- Score calculation across 17K+ trees
- Phylogenetic placement

The bottleneck is not in BWA alignment itself, but in the downstream score computation.

---

### Architecture Insight: Why Limited Gains?

The current tronko-assign threading model is:

```
Main Thread
    ├── Thread 1: runAssignmentOnChunk_WithBWA(chunk 1)
    │       └── run_bwa (single-threaded BWA)
    │       └── place_paired / phylogenetic placement
    │
    ├── Thread 2: runAssignmentOnChunk_WithBWA(chunk 2)
    │       └── run_bwa (single-threaded BWA)
    │       └── place_paired / phylogenetic placement
    │
    └── Thread N: runAssignmentOnChunk_WithBWA(chunk N)
```

Each **chunk thread** already runs its own BWA alignment. Adding threading **inside** BWA (`kt_for` parallelization) would only help if:
1. You're running with fewer chunk threads than cores, OR
2. The alignment of individual chunks is the bottleneck

In practice, with `--number-of-cores N`, you already have N concurrent BWA operations. Adding internal BWA threading would mostly cause thread contention.

---

## Recommendations

### If You Want Quick Wins

Skip BWA threading. Instead:

1. **Enable OpenMP for WFA2** - 5 min change, 10-20% gain on alignment step:
   ```makefile
   LIBS = -lm -pthread -lz -lrt -std=gnu99 -fopenmp
   ```

2. **Add compiler flags** - 5 min change, 10-30% overall:
   ```makefile
   OPTIMIZATION = -O3 -march=native -mtune=native -ffast-math
   ```

3. **Implement early termination** - 1-2 days, 2-5x gain

### If You Want to Enable BWA Threading Anyway

1. Use the "Minimal Fix" approach above
2. Add a `--bwa-threads N` command-line option (default: 1 for backward compatibility)
3. Test thoroughly with paired and single-end reads
4. Benchmark to confirm actual gains in your workload

---

## Code References

### Current Threading Hardcoding
- `tronko-assign/bwa_source_files/fastmap.c:697` - `opt->n_threads=1`
- `tronko-assign/bwa_source_files/fastmap.c:842` - `kt_pipeline(1, ...)`

### Broken Mutex Pattern
- `tronko-assign/bwa_source_files/bwamem.c:877-884` - Single-end reverse complement
- `tronko-assign/bwa_source_files/bwamem.c:890-897` - Paired-end read1 reverse complement
- `tronko-assign/bwa_source_files/bwamem.c:904-911` - Paired-end read2 reverse complement

### Threading Infrastructure
- `tronko-assign/bwa_source_files/kthread.c:49-61` - `kt_for()` implementation
- `tronko-assign/bwa_source_files/kthread.c:123-151` - `kt_pipeline()` implementation
- `tronko-assign/bwa_source_files/bwamem.c:1288-1320` - `mem_process_seqs()` with threading

### Global State
- `tronko-assign/global.h:112-113` - Query matrix declarations
- `tronko-assign/bwa_source_files/bwamem.c:45` - `global_bns` static variable

---

## Related Research

- `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md` - Full optimization priority analysis
- `thoughts/shared/plans/2026-01-01-tier2-optimization-implementation.md` - Tier 2 implementation plan (includes mmap BWA wrapper)
- `thoughts/shared/research/2025-12-29-tronko-assign-streaming-architecture.md` - I/O and data flow overview

---

## Open Questions

1. **Contention testing**: How much do chunk threads contend when BWA internal threading is enabled?
2. **Optimal thread ratio**: What's the best ratio of chunk threads to BWA threads?
3. **Memory pressure**: Does internal BWA threading increase memory usage significantly?
4. **Accuracy verification**: Does enabling threading change any output (it shouldn't, but needs testing)?

---

## Conclusion

Enabling multi-threaded BWA alignment is **feasible but low priority**. The infrastructure exists, but tronko-specific modifications introduced thread-safety issues that need fixing. Expected gains are modest (10-20%) compared to other optimizations that offer 2-10x improvements.

**Bottom line**: Fix if you have spare time, but focus on Tier 1 optimizations first.
