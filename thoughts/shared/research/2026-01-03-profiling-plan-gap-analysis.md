---
date: 2026-01-03T12:00:00-08:00
researcher: Claude
git_commit: d47ed018e85e993edf7124a67692377f2987b3b8
branch: experimental
repository: tronko
topic: "Gap Analysis: CPU Profiling Plan for tronko-assign Bottleneck Identification"
tags: [research, performance, profiling, bottleneck, optimization]
status: complete
last_updated: 2026-01-03
last_updated_by: Claude
---

# Research: Gap Analysis for CPU Profiling Plan

**Date**: 2026-01-03T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: d47ed018e85e993edf7124a67692377f2987b3b8
**Branch**: experimental
**Repository**: tronko

## Research Question

What gaps exist in the current CPU profiling plan (`thoughts/shared/plans/2026-01-03-cpu-profiling-bottleneck-analysis.md`) for identifying bottlenecks in tronko-assign, and what additional profiling strategies should be considered?

## Summary

The current profiling plan is solid for high-level CPU profiling with `perf`, but lacks detail on several critical areas discovered through codebase analysis:

1. **Character Branching in Scoring** - The innermost loop uses chained if-else for nucleotide comparison, causing branch misprediction
2. **Memory Allocation Churn** - Specific allocation sites in LCA computation and batch processing cause per-read/per-batch malloc/free
3. **Recursive Call Overhead** - Tree traversal uses recursion with ~12 parameters per call
4. **Data Layout Inefficiency** - Posterior probability arrays use strided access (4-element stride) that may cause cache issues
5. **WFA2/BWA Internal Behavior** - Both libraries have internal SIMD and threading that need separate profiling
6. **Reference Loading I/O** - Binary format uses many small reads (2-8 bytes) vs bulk reads

---

## Detailed Findings

### 1. Hot Path Analysis - What the Plan Gets Right

The plan correctly identifies the key functions:
- `getscore_Arr` / `getscore_Arr_ncbi` - score accumulation
- `assignScores_Arr_paired` / `assignScores_Arr_single` - tree traversal
- `wavefront_align` (WFA2) - detailed alignment
- `bwa_*` functions - BWA alignment

**Verified Call Hierarchy:**
```
main() [tronko-assign.c:902]
  └→ BATCH LOOP [lines 1318-1430]
       └→ pthread_create(runAssignmentOnChunk_WithBWA) [line 1369]
            └→ run_bwa() → main_mem() → mem_process_seqs() [BWA alignment]
            └→ FOR EACH READ [line 369]:
                 └→ place_paired() [placement.c:60]
                      └→ FOR EACH BWA MATCH:
                           └→ wavefront_align() [placement.c:116]
                           └→ assignScores_Arr_paired() [assignment.c:24]
                                └→ RECURSIVE FOR EACH NODE:
                                     └→ getscore_Arr() [assignment.c:172] ← INNERMOST LOOP
```

### 2. Critical Gap: Character Branching in getscore_Arr

**Location:** `assignment.c:206-225`

The scoring loop uses **chained if-else** for nucleotide comparison:

```c
if (locQuery[i]=='a' || locQuery[i]=='A'){
    score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
}else if (locQuery[i]=='c' || locQuery[i]=='C'){
    score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
}else if (locQuery[i]=='g' || locQuery[i]=='G'){
    score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 2)];
}else if (locQuery[i]=='t' || locQuery[i]=='T'){
    score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 3)];
}
```

**Why This Matters:**
- Each iteration has 4 conditional branches
- Branch predictor cannot learn pattern (nucleotides are ~25% each)
- Expected branch misprediction rate: ~75% per iteration
- With 150-300bp queries, this is 450-900+ mispredictions per node

**Profiling Addition:**
```bash
# Measure branch mispredictions specifically
perf stat -e branch-misses,branches \
    ./tronko-assign [args...]

# Annotate to see if mispredictions cluster in getscore_Arr
perf annotate getscore_Arr --symbol-filter=getscore_Arr
```

**Optimization Opportunity:**
Replace with lookup table:
```c
// Pre-computed: 'A'->0, 'a'->0, 'C'->1, 'c'->1, etc.
static const int8_t nuc_to_idx[256] = { ... };

// In loop:
int idx = nuc_to_idx[(unsigned char)locQuery[i]];
if (idx >= 0) {
    score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], idx)];
}
```

---

### 3. Critical Gap: Memory Allocation Churn Sites

The plan mentions "time in malloc/free" but doesn't identify specific allocation sites. Research found:

#### Per-LCA-Call Allocations (HIGH FREQUENCY)

**Location:** `tronko-assign.c:215,223,250,262`

```c
// LCA_of_nodes() - called for EVERY read
int* ancestors = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
// ... computation ...
free(ancestors);

// getLCAofArray_Arr_Multiple() - called for EVERY read
int* minNodes = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
// ... computation ...
free(minNodes);
```

**Frequency:** 2 malloc/free pairs per read processed
**Size:** Variable, proportional to tree size (thousands of ints)

#### Per-Batch Allocations (MEDIUM FREQUENCY)

**Location:** `tronko-assign.c:1358-1361`

```c
// Allocated EVERY batch
mstr[i].str->taxonPath = (char**) malloc((end-start)*(sizeof(char *)));
for(k=0; k<end-start; k++){
    mstr[i].str->taxonPath[k] = malloc((max_name_length+max_lineTaxonomy+120)*(sizeof(char)));
}
// Freed at lines 1425-1428
```

**Frequency:** batch_size * num_threads allocations per batch

**Profiling Addition:**
```bash
# Memory allocation tracing
perf record -e syscalls:sys_enter_brk,syscalls:sys_enter_mmap \
    ./tronko-assign [args...]

# Or use LD_PRELOAD with malloc tracing:
LD_PRELOAD=/usr/lib/libtcmalloc.so HEAPPROFILE=/tmp/heap \
    ./tronko-assign [args...]
```

---

### 4. Critical Gap: Recursive Call Overhead

**Location:** `assignment.c:24-93`

The `assignScores_Arr_paired` function is recursive with **12+ parameters**:

```c
void assignScores_Arr_paired(
    int rootNum,                    // 4 bytes
    int node,                       // 4 bytes
    char *locQuery,                 // 8 bytes (pointer)
    int *positions,                 // 8 bytes (pointer)
    type_of_PP ***scores,           // 8 bytes (pointer)
    int alength,                    // 4 bytes
    int search_number,              // 4 bytes
    int print_all_nodes,            // 4 bytes
    FILE* site_scores_file,         // 8 bytes (pointer)
    char* readname,                 // 8 bytes (pointer)
    int early_termination,          // 4 bytes
    type_of_PP *best_score,         // 8 bytes (pointer)
    int *strikes,                   // 8 bytes (pointer)
    type_of_PP strike_box,          // 8 bytes (double)
    int max_strikes,                // 4 bytes
    int enable_pruning,             // 4 bytes
    type_of_PP pruning_threshold    // 8 bytes (double)
)
```

**Stack Frame Size:** ~96 bytes per recursive call
**Tree Depth:** Typically 10-20 levels for balanced trees
**Per-Node Overhead:** ~960-1920 bytes of stack activity per tree traversal

**Profiling Addition:**
```bash
# Profile call/ret instruction overhead
perf stat -e instructions,cycles,call,ret \
    ./tronko-assign [args...]

# Check if stack operations dominate
perf annotate assignScores_Arr_paired
```

---

### 5. Critical Gap: Data Layout and Cache Effects

**Location:** `global.h:29`, `assignment.c:206-225`

The `PP_IDX` macro:
```c
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))
```

Creates this memory layout:
```
[pos0_A, pos0_C, pos0_G, pos0_T, pos1_A, pos1_C, pos1_G, pos1_T, ...]
```

**Access Pattern in getscore_Arr:**
- Iterates through positions sequentially
- But only accesses ONE nucleotide per position (based on query character)
- Results in stride-1 access within 4-element groups

**Potential Issue:**
- Cache line (64 bytes) holds 8 doubles or 16 floats
- Each access uses only 1 of 4 values per position
- 75% of loaded cache data is unused per iteration

**Alternative Layout (not currently implemented):**
```
// Structure-of-Arrays instead of Array-of-Structures
posteriornc_A[numbase]
posteriornc_C[numbase]
posteriornc_G[numbase]
posteriornc_T[numbase]
```

**Profiling Addition:**
```bash
# Cache efficiency metrics
perf stat -e L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
    ./tronko-assign [args...]

# Memory bandwidth
perf stat -e mem_load_retired.l3_miss,mem_load_retired.fb_hit \
    ./tronko-assign [args...]
```

---

### 6. Gap: BWA Internal Threading and SIMD

**Location:** `bwa_source_files/kthread.c`, `bwa_source_files/ksw.c`

BWA has its own threading model with work stealing:
```c
// kthread.c:31,40 - Lock-free work stealing
k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
```

And explicit SSE2 SIMD:
```c
// ksw.c:161-176 - SSE2 Smith-Waterman
_mm_max_epu8(H, _mm_adds_epu8(E, vgape));
_mm_max_epu8(H, _mm_subs_epu8(_mm_adds_epu8(Hd, vgapo), vgape1));
```

**The Current Plan Treats BWA as a Black Box**

**Profiling Addition:**
```bash
# Profile BWA functions specifically
perf record -g --call-graph dwarf ./tronko-assign [args...]
perf report --symbol-filter='bwa\|mem_\|ksw_'

# Check if BWA SIMD is actually running
perf stat -e fp_arith_inst_retired.128b_packed_single,fp_arith_inst_retired.scalar_single \
    ./tronko-assign [args...]
```

---

### 7. Gap: WFA2 Vectorization Status

**Location:** `WFA2/commons.h:253-261`, `WFA2/wavefront_compute_affine.c:62`

WFA2 has vectorization hints:
```c
#define PRAGMA_LOOP_VECTORIZE _Pragma("GCC ivdep")

// Used in compute loops:
PRAGMA_LOOP_VECTORIZE
for (k=lo;k<=hi;++k) { ... }
```

**But:** The Makefile doesn't pass `-fopenmp` by default (lines 31-36 only enable it with `ENABLE_OPENMP=1`)

**Profiling Addition:**
```bash
# Check vectorization effectiveness
perf stat -e fp_arith_inst_retired.256b_packed_double,fp_arith_inst_retired.scalar_double \
    ./tronko-assign [args...]

# Compare with OpenMP enabled
make clean && make ENABLE_OPENMP=1
perf stat [same counters] ./tronko-assign [args...]
```

---

### 8. Gap: Reference Database I/O Pattern

**Location:** `readreference.c:12-46` (binary format readers)

Binary format uses many small reads:
```c
// read_u16() - 2-byte reads
fread(b, 1, 2, fp);

// read_u32() - 4-byte reads
fread(b, 1, 4, fp);

// read_u64() - 8-byte reads
fread(b, 1, 8, fp);
```

**Bulk reads only for posteriors:**
```c
// readreference.c:968-978
fread(treeArr[t][n].posteriornc, sizeof(float), count, fp);
```

**Profiling Addition:**
```bash
# I/O tracing
strace -c ./tronko-assign [args...] 2>&1 | grep -E 'read|write'

# Or with perf
perf stat -e block:block_rq_issue,block:block_rq_complete \
    ./tronko-assign [args...]
```

---

### 9. Gap: Synchronization Contention

**Location:** `tronko-assign.c:111`, `logger.c:52`, `crash_debug.c:45`

Three mutex types exist:
```c
static pthread_mutex_t g_overflow_stats_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t log_mutex;
pthread_mutex_t context_mutex;
```

**Profiling Addition:**
```bash
# Lock contention analysis
perf lock record ./tronko-assign -C 8 [args...]
perf lock report

# Or with mutrace
LD_PRELOAD=/usr/lib/libmutrace.so ./tronko-assign -C 8 [args...]
```

---

## Revised Profiling Checklist

### Phase 1: Basic CPU Characteristics (as in original plan)
- [x] IPC (instructions per cycle)
- [x] Cache miss rates
- [x] Branch misprediction rates

### Phase 2: Function Hotspots (as in original plan)
- [x] `perf record -g`
- [x] `perf report`
- [x] Flame graph

### Phase 3: NEW - Branch Misprediction Analysis
- [ ] Annotate `getscore_Arr` for branch instructions
- [ ] Measure branch miss rate specifically in scoring loop
- [ ] Test lookup table optimization

### Phase 4: NEW - Memory Allocation Profiling
- [ ] Profile malloc/free frequency
- [ ] Identify if LCA temp arrays dominate allocation time
- [ ] Test pre-allocated buffers

### Phase 5: NEW - Cache Efficiency Analysis
- [ ] L1/L2/L3 miss rates
- [ ] Memory bandwidth utilization
- [ ] Evaluate data layout alternatives

### Phase 6: NEW - Library-Specific Profiling
- [ ] BWA function breakdown
- [ ] WFA2 vectorization verification
- [ ] Compare with ENABLE_OPENMP=1

### Phase 7: NEW - Multi-threaded Analysis (if using -C > 1)
- [ ] Lock contention with `perf lock`
- [ ] Thread idle time
- [ ] Work distribution balance

---

## Code References

- `assignment.c:172-239` - getscore_Arr (innermost scoring loop)
- `assignment.c:24-93` - assignScores_Arr_paired (recursive tree traversal)
- `tronko-assign.c:215,223,250,262` - LCA temp allocation sites
- `tronko-assign.c:1358-1361` - Per-batch result allocation
- `global.h:29` - PP_IDX macro (data layout)
- `bwa_source_files/ksw.c:149-212` - BWA SSE2 SIMD
- `WFA2/commons.h:253-261` - WFA2 vectorization pragmas
- `readreference.c:12-46` - Binary format small read functions

## Architecture Insights

1. **Execution is dominated by tree traversal depth** - Each query visits every node in matching trees
2. **Recursion with many parameters** - Stack pressure from 12+ parameter recursive calls
3. **Branch-heavy innermost loop** - Character comparisons prevent vectorization
4. **Embedded libraries have hidden complexity** - BWA and WFA2 each have internal optimizations

## Historical Context

From `EXPERIMENTS_LOG.md`:
- Tier 2 optimizations (mmap, SIMD max-finding, two-phase screening) showed **no measurable impact**
- This strongly suggests the bottleneck is NOT in areas previously targeted
- Branch misprediction and memory allocation are more likely culprits

From `thoughts/shared/research/2026-01-01-tier1-benchmark-results.md`:
- Early termination is **broken** (0% accuracy)
- Subtree pruning shows ~0% speedup on small datasets
- FAST_MATH provides only ~2% improvement

## Related Research

- `thoughts/shared/plans/2026-01-03-cpu-profiling-bottleneck-analysis.md` - Original profiling plan
- `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md` - Optimization tiers
- `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md` - SIMD opportunities

## Open Questions

1. What is the actual branch misprediction rate in `getscore_Arr`?
2. Does the recursive call overhead dominate over the scoring computation itself?
3. Would converting tree traversal to iterative (explicit stack) help?
4. Is the 4-element stride in PP_IDX causing cache line waste?
5. How much time is spent in BWA vs WFA2 vs tronko-assign scoring?
