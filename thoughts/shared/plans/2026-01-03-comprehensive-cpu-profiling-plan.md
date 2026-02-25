# Comprehensive CPU Profiling Plan for tronko-assign Bottleneck Identification

## Overview

This plan extends the original CPU profiling plan (`2026-01-03-cpu-profiling-bottleneck-analysis.md`) based on gap analysis research. It addresses seven specific bottleneck areas identified through codebase analysis that the original plan missed:

1. Branch misprediction in nucleotide scoring
2. Memory allocation churn in LCA computation
3. Recursive call overhead in tree traversal
4. Cache inefficiency from data layout
5. BWA/WFA2 library internal behavior
6. Reference database I/O patterns
7. Multi-threaded lock contention

## Current State Analysis

### What We Know
- Tier 2 optimizations (mmap, SIMD max-finding, two-phase screening) had **no measurable impact**
- This indicates bottlenecks are NOT in the areas we previously targeted
- The innermost loop (`getscore_Arr`) uses chained if-else for nucleotide comparison
- Per-read memory allocations occur in LCA computation
- Tree traversal is recursive with 12+ parameters per call

### Key Code Locations
- `assignment.c:172-239` - getscore_Arr (innermost scoring loop)
- `assignment.c:24-93` - assignScores_Arr_paired (recursive tree traversal)
- `tronko-assign.c:215,223,250,262` - LCA temp allocation sites
- `tronko-assign.c:1358-1361` - Per-batch result allocation
- `global.h:29` - PP_IDX macro (data layout)

## Desired End State

After completing this profiling plan, we will have:
1. Quantitative data on where CPU time is spent (function-level and line-level)
2. Branch misprediction rates for the nucleotide scoring loop
3. Memory allocation frequency and overhead measurements
4. Cache efficiency metrics (L1/L2/L3 miss rates)
5. BWA vs WFA2 vs tronko-assign time breakdown
6. Lock contention data for multi-threaded runs
7. Clear prioritization for which optimizations to pursue next

### Verification
- All profiling commands documented with expected output formats
- Results saved to `/tmp/profiling/` directory with timestamped files
- Summary report generated comparing actual bottlenecks to prior assumptions

## What We're NOT Doing

- **Not implementing optimizations** - This plan is purely for measurement
- **Not modifying tronko-assign code** - We're profiling existing behavior
- **Not profiling tronko-build** - Focus is on assignment bottlenecks
- **Not doing GPU profiling** - CPU-only analysis for now
- **Not micro-benchmarking individual functions** - Focus on realistic workloads

## Implementation Approach

Use a layered profiling approach:
1. Start with high-level CPU statistics to identify broad categories
2. Drill down with function-level profiling
3. Investigate specific bottlenecks (branches, cache, allocations)
4. Profile library behavior (BWA, WFA2) separately
5. Analyze multi-threaded behavior

All phases use `perf` (Linux performance counters) with supplemental tools for specific analyses.

---

## Phase 1: Environment Setup and Baseline

### Overview
Prepare the profiling environment, build debug-enabled binaries, and establish baseline measurements.

### Changes Required

#### 1. Create Profiling Directory Structure
```bash
mkdir -p /tmp/profiling/{perf,flamegraph,memory,reports}
```

#### 2. Build Debug Binary with Symbols
```bash
cd tronko-assign
make clean
make debug
# This adds -g flag for debug symbols while keeping -O3 optimization
```

#### 3. Install Required Tools
```bash
# On Arch Linux
sudo pacman -S perf linux-tools flamegraph

# Verify FlameGraph scripts are available
which stackcollapse-perf.pl flamegraph.pl
```

#### 4. Verify perf Access
```bash
# Check if perf can access hardware counters
perf stat ls

# If permission denied, run as root or adjust kernel settings:
sudo sysctl -w kernel.perf_event_paranoid=-1
```

### Test Command
```bash
# Baseline run to verify setup
cd /home/jimjeffers/Work/tronko/tronko-assign
./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/baseline_results.txt
```

### Success Criteria

#### Automated Verification:
- [ ] Debug binary builds successfully: `make debug` exits 0
- [ ] Output file created: `ls /tmp/profiling/baseline_results.txt`
- [ ] perf can record: `perf stat -e cycles ./tronko-assign --help 2>&1 | grep cycles`

#### Manual Verification:
- [ ] Debug symbols present: `nm tronko-assign | grep getscore_Arr` shows symbols
- [ ] perf annotate works: `perf annotate` shows source code (not just assembly)

---

## Phase 2: High-Level CPU Characteristics

### Overview
Capture overall CPU behavior metrics to determine if we're CPU-bound, memory-bound, or branch-bound.

### Commands

#### 2.1 Basic CPU Counters
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase2_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase2_basic_counters.txt
```

#### 2.2 Extended Metrics (if available)
```bash
perf stat -e cycles,instructions,L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,branch-loads,branch-load-misses,\
context-switches,cpu-migrations \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase2b_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase2_extended_counters.txt
```

### Analysis Criteria

| Metric | Healthy | Concerning | Action |
|--------|---------|------------|--------|
| IPC (instructions/cycle) | > 1.5 | < 1.0 | Investigate stalls |
| Cache miss rate | < 3% | > 5% | Phase 5 (cache analysis) |
| Branch miss rate | < 2% | > 5% | Phase 3 (branch analysis) |
| Context switches | Low | High | Phase 7 (threading) |

### Success Criteria

#### Automated Verification:
- [ ] Counter file created: `test -s /tmp/profiling/perf/phase2_basic_counters.txt`
- [ ] Contains IPC data: `grep -q "instructions" /tmp/profiling/perf/phase2_basic_counters.txt`

#### Manual Verification:
- [ ] Calculate IPC from output (instructions / cycles)
- [ ] Calculate cache miss rate ((cache-misses / cache-references) * 100)
- [ ] Calculate branch miss rate ((branch-misses / branches) * 100)
- [ ] Document these values in `/tmp/profiling/reports/phase2_summary.md`

---

## Phase 3: Function-Level Hotspot Analysis

### Overview
Identify which functions consume the most CPU time using sampling-based profiling.

### Commands

#### 3.1 Record with Call Graph (DWARF)
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

perf record -g --call-graph dwarf -o /tmp/profiling/perf/phase3_callgraph.data \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase3_results.txt
```

#### 3.2 Generate Report
```bash
perf report -i /tmp/profiling/perf/phase3_callgraph.data --sort=dso,symbol \
    > /tmp/profiling/perf/phase3_hotspots.txt 2>&1

# Top 20 functions only
perf report -i /tmp/profiling/perf/phase3_callgraph.data --stdio --sort=symbol \
    | head -50 > /tmp/profiling/perf/phase3_top_functions.txt
```

#### 3.3 Generate Flame Graph
```bash
perf script -i /tmp/profiling/perf/phase3_callgraph.data | \
    stackcollapse-perf.pl | \
    flamegraph.pl > /tmp/profiling/flamegraph/phase3_flamegraph.svg
```

### Expected Hotspot Candidates

Based on the research, we expect to see:
1. `getscore_Arr` / `getscore_Arr_ncbi` - scoring loop
2. `assignScores_Arr_paired` - tree traversal
3. `wavefront_align` / WFA2 functions - alignment
4. `bwa_*` / `mem_*` / `ksw_*` - BWA alignment
5. `malloc` / `free` - memory allocation

**Key Finding:** If the actual hotspot differs from this list, our optimization assumptions were wrong.

### Success Criteria

#### Automated Verification:
- [ ] perf.data file created: `test -s /tmp/profiling/perf/phase3_callgraph.data`
- [ ] Hotspot report created: `test -s /tmp/profiling/perf/phase3_hotspots.txt`
- [ ] Flame graph created: `test -s /tmp/profiling/flamegraph/phase3_flamegraph.svg`

#### Manual Verification:
- [ ] Open flame graph in browser: `xdg-open /tmp/profiling/flamegraph/phase3_flamegraph.svg`
- [ ] Identify top 5 functions by CPU time percentage
- [ ] Document findings in `/tmp/profiling/reports/phase3_summary.md`

---

## Phase 4: Branch Misprediction Analysis (Critical Gap)

### Overview
The research identified that `getscore_Arr` uses chained if-else for nucleotide comparison, potentially causing ~75% branch misprediction per iteration.

### Commands

#### 4.1 Branch-Specific Counters
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

perf stat -e branches,branch-misses,branch-loads,branch-load-misses \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase4_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase4_branch_counters.txt
```

#### 4.2 Record Branch Events
```bash
perf record -e branch-misses -g --call-graph dwarf \
    -o /tmp/profiling/perf/phase4_branch_misses.data \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase4b_results.txt
```

#### 4.3 Annotate getscore_Arr
```bash
perf annotate -i /tmp/profiling/perf/phase4_branch_misses.data \
    --symbol=getscore_Arr \
    > /tmp/profiling/perf/phase4_getscore_annotated.txt 2>&1
```

### Analysis Focus

Look for branch instructions clustering around lines 206-225 in assignment.c:
```c
if (locQuery[i]=='a' || locQuery[i]=='A'){  // Branch 1
    ...
}else if (locQuery[i]=='c' || locQuery[i]=='C'){  // Branch 2
    ...
}else if (locQuery[i]=='g' || locQuery[i]=='G'){  // Branch 3
    ...
}else if (locQuery[i]=='t' || locQuery[i]=='T'){  // Branch 4
    ...
}
```

**Hypothesis:** With ~25% probability for each nucleotide, branch predictor cannot learn pattern, causing ~75% misprediction rate in this loop.

### Success Criteria

#### Automated Verification:
- [ ] Branch counter file exists: `test -s /tmp/profiling/perf/phase4_branch_counters.txt`
- [ ] Annotation file exists: `test -s /tmp/profiling/perf/phase4_getscore_annotated.txt`

#### Manual Verification:
- [ ] Calculate overall branch misprediction rate
- [ ] Check if `getscore_Arr` appears in top branch misprediction sources
- [ ] Identify if branch mispredictions cluster at nucleotide comparison lines
- [ ] Document findings in `/tmp/profiling/reports/phase4_summary.md`
- [ ] **Key Question Answered:** Is branch misprediction rate > 5% in getscore_Arr?

---

## Phase 5: Memory Allocation Profiling (Critical Gap)

### Overview
Research identified per-read malloc/free cycles in LCA computation at `tronko-assign.c:215,223,250,262`.

### Commands

#### 5.1 Record malloc/free Activity
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# Method 1: perf probe on malloc/free (if available)
perf record -e 'probe_libc:*' -g \
    -o /tmp/profiling/perf/phase5_malloc.data \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase5_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase5_malloc_output.txt
```

#### 5.2 Alternative: ltrace for Allocation Counts
```bash
# Count malloc/free calls
ltrace -e malloc+free -c \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase5b_results.txt \
    2>&1 | tee /tmp/profiling/memory/phase5_allocation_counts.txt
```

#### 5.3 Strace for System Calls
```bash
# Count memory-related system calls
strace -c -e brk,mmap,munmap \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase5c_results.txt \
    2>&1 | tee /tmp/profiling/memory/phase5_syscall_counts.txt
```

### Allocation Sites to Monitor

| Location | Function | Frequency | Purpose |
|----------|----------|-----------|---------|
| `tronko-assign.c:215` | `LCA_of_nodes` | Per-read | ancestors array |
| `tronko-assign.c:250` | `getLCAofArray_Arr_Multiple` | Per-read | minNodes array |
| `tronko-assign.c:1358-1361` | Batch loop | Per-batch | taxonPath strings |

### Success Criteria

#### Automated Verification:
- [ ] Allocation count file exists: `test -s /tmp/profiling/memory/phase5_allocation_counts.txt`
- [ ] Syscall count file exists: `test -s /tmp/profiling/memory/phase5_syscall_counts.txt`

#### Manual Verification:
- [ ] Count total malloc/free calls
- [ ] Estimate if allocation frequency matches expected (2 per read for LCA)
- [ ] Calculate percentage of CPU time in allocation functions (from Phase 3)
- [ ] Document findings in `/tmp/profiling/reports/phase5_summary.md`
- [ ] **Key Question Answered:** Is > 5% of CPU time spent in malloc/free?

---

## Phase 6: Cache Efficiency Analysis (Critical Gap)

### Overview
The PP_IDX macro creates 4-element interleaved layout. Scoring loop accesses only 1 of 4 values, potentially wasting 75% of cache line loads.

### Commands

#### 6.1 L1/L2/L3 Cache Metrics
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

perf stat -e L1-dcache-loads,L1-dcache-load-misses,\
L1-icache-load-misses,\
LLC-loads,LLC-load-misses \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase6_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase6_cache_counters.txt
```

#### 6.2 Memory Bandwidth Metrics
```bash
# Note: These events may not be available on all CPUs
perf stat -e mem_load_retired.l3_miss,mem_load_retired.fb_hit,\
mem_inst_retired.all_loads,mem_inst_retired.all_stores \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase6b_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase6_memory_bandwidth.txt
```

#### 6.3 Record Cache Misses by Location
```bash
perf record -e L1-dcache-load-misses -g --call-graph dwarf \
    -o /tmp/profiling/perf/phase6_cache_misses.data \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase6c_results.txt

perf report -i /tmp/profiling/perf/phase6_cache_misses.data --stdio \
    > /tmp/profiling/perf/phase6_cache_miss_functions.txt
```

### Current Data Layout
```
PP_IDX(pos, nuc) = pos * 4 + nuc

Memory: [pos0_A, pos0_C, pos0_G, pos0_T, pos1_A, pos1_C, pos1_G, pos1_T, ...]
```

Access pattern in `getscore_Arr`: accesses only ONE of the 4 nucleotide values per position based on query character.

### Success Criteria

#### Automated Verification:
- [ ] Cache counter file exists: `test -s /tmp/profiling/perf/phase6_cache_counters.txt`
- [ ] Cache miss function report exists: `test -s /tmp/profiling/perf/phase6_cache_miss_functions.txt`

#### Manual Verification:
- [ ] Calculate L1 cache miss rate
- [ ] Calculate LLC (L3) cache miss rate
- [ ] Check if `getscore_Arr` or related functions dominate cache misses
- [ ] Document findings in `/tmp/profiling/reports/phase6_summary.md`
- [ ] **Key Question Answered:** Is cache miss rate > 5% and attributable to scoring loop?

---

## Phase 7: Library-Specific Profiling (BWA and WFA2)

### Overview
BWA and WFA2 are treated as black boxes in current profiling. They have internal SIMD and threading that may dominate runtime.

### Commands

#### 7.1 BWA Function Breakdown
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# Record and filter for BWA functions
perf record -g --call-graph dwarf \
    -o /tmp/profiling/perf/phase7_libraries.data \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase7_results.txt

# Filter for BWA symbols
perf report -i /tmp/profiling/perf/phase7_libraries.data --stdio \
    | grep -E 'bwa|mem_|ksw_' \
    > /tmp/profiling/perf/phase7_bwa_functions.txt

# Filter for WFA2 symbols
perf report -i /tmp/profiling/perf/phase7_libraries.data --stdio \
    | grep -E 'wavefront|wfa|WFA' \
    > /tmp/profiling/perf/phase7_wfa2_functions.txt
```

#### 7.2 Check SIMD Utilization
```bash
# Check if vectorized instructions are being used
perf stat -e fp_arith_inst_retired.128b_packed_single,\
fp_arith_inst_retired.256b_packed_single,\
fp_arith_inst_retired.scalar_single \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase7b_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase7_simd_counters.txt
```

#### 7.3 Compare with OpenMP-Enabled WFA2
```bash
# Build with OpenMP
cd /home/jimjeffers/Work/tronko/tronko-assign
make clean
make debug ENABLE_OPENMP=1

# Profile with OpenMP
perf stat -e cycles,instructions \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase7c_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase7_openmp_enabled.txt

# Restore non-OpenMP build
make clean
make debug
```

### Success Criteria

#### Automated Verification:
- [ ] BWA function file exists: `test -s /tmp/profiling/perf/phase7_bwa_functions.txt`
- [ ] WFA2 function file exists: `test -s /tmp/profiling/perf/phase7_wfa2_functions.txt`

#### Manual Verification:
- [ ] Calculate BWA time as percentage of total
- [ ] Calculate WFA2 time as percentage of total
- [ ] Calculate tronko-assign scoring time as percentage of total
- [ ] Document the three-way breakdown in `/tmp/profiling/reports/phase7_summary.md`
- [ ] **Key Question Answered:** What percentage of time is in BWA vs WFA2 vs tronko scoring?

---

## Phase 8: Multi-threaded Analysis

### Overview
Analyze lock contention and thread balance when using multiple cores (-C option).

### Commands

#### 8.1 Run with Multiple Threads
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# 4-core run with thread tracing
perf stat -e context-switches,cpu-migrations,\
sched:sched_switch,sched:sched_wakeup \
    ./tronko-assign -r -C 4 \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase8_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase8_threading_counters.txt
```

#### 8.2 Lock Contention Analysis
```bash
# Record lock events
perf lock record -o /tmp/profiling/perf/phase8_lock.data \
    ./tronko-assign -r -C 4 \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase8b_results.txt

perf lock report -i /tmp/profiling/perf/phase8_lock.data \
    > /tmp/profiling/perf/phase8_lock_report.txt 2>&1
```

#### 8.3 Per-Thread Flame Graph
```bash
# Record with thread info
perf record -g --call-graph dwarf -T \
    -o /tmp/profiling/perf/phase8_threaded.data \
    ./tronko-assign -r -C 4 \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase8c_results.txt

perf script -i /tmp/profiling/perf/phase8_threaded.data | \
    stackcollapse-perf.pl --pid | \
    flamegraph.pl > /tmp/profiling/flamegraph/phase8_threaded_flamegraph.svg
```

### Known Mutex Locations

| Mutex | Location | Purpose |
|-------|----------|---------|
| `g_overflow_stats_mutex` | `tronko-assign.c:111` | Overflow stats |
| `log_mutex` | `logger.c:52` | Logging |
| `context_mutex` | `crash_debug.c:45` | Crash context |

### Success Criteria

#### Automated Verification:
- [ ] Threading counter file exists: `test -s /tmp/profiling/perf/phase8_threading_counters.txt`
- [ ] Lock report exists: `test -s /tmp/profiling/perf/phase8_lock_report.txt`
- [ ] Threaded flame graph exists: `test -s /tmp/profiling/flamegraph/phase8_threaded_flamegraph.svg`

#### Manual Verification:
- [ ] Identify if any mutex has significant wait time
- [ ] Check if threads have balanced work (via flame graph width)
- [ ] Document findings in `/tmp/profiling/reports/phase8_summary.md`
- [ ] **Key Question Answered:** Is lock contention > 1% of runtime?

---

## Phase 9: I/O Pattern Analysis

### Overview
Check if reference loading or output writing causes bottlenecks.

### Commands

#### 9.1 I/O Tracing
```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# System call tracing for I/O
strace -c -e read,write,pread64,pwrite64 \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase9_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase9_io_syscalls.txt
```

#### 9.2 Block I/O Events
```bash
perf stat -e block:block_rq_issue,block:block_rq_complete \
    ./tronko-assign -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/profiling/phase9b_results.txt \
    2>&1 | tee /tmp/profiling/perf/phase9_block_io.txt
```

### Success Criteria

#### Automated Verification:
- [ ] I/O syscall file exists: `test -s /tmp/profiling/perf/phase9_io_syscalls.txt`
- [ ] Block I/O file exists: `test -s /tmp/profiling/perf/phase9_block_io.txt`

#### Manual Verification:
- [ ] Check if I/O time is significant vs compute time
- [ ] Identify if many small reads (reference loading) cause overhead
- [ ] Document findings in `/tmp/profiling/reports/phase9_summary.md`

---

## Phase 10: Generate Final Report

### Overview
Synthesize all findings into actionable recommendations.

### Create Summary Report

After completing all phases, create `/tmp/profiling/reports/FINAL_PROFILING_REPORT.md`:

```markdown
# tronko-assign Profiling Results Summary

## Date: [DATE]
## Test Dataset: single_tree example

## Key Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| IPC | [X.XX] | [CPU/Memory/Branch bound] |
| Cache Miss Rate | [X.X%] | [Low/Medium/High] |
| Branch Miss Rate | [X.X%] | [Low/Medium/High] |
| Time in malloc/free | [X.X%] | [Low/Medium/High] |

## CPU Time Breakdown

| Component | Percentage |
|-----------|------------|
| BWA Alignment | [XX%] |
| WFA2 Alignment | [XX%] |
| Tree Scoring (getscore_Arr) | [XX%] |
| Tree Traversal (assignScores_*) | [XX%] |
| Memory Allocation | [XX%] |
| Other | [XX%] |

## Bottleneck Ranking

1. [Primary bottleneck] - [XX%] of time
2. [Secondary bottleneck] - [XX%] of time
3. [Tertiary bottleneck] - [XX%] of time

## Recommended Optimizations

Based on profiling data, prioritize:

1. **[Optimization]** - Targets [XX%] of runtime
2. **[Optimization]** - Targets [XX%] of runtime
3. **[Optimization]** - Targets [XX%] of runtime

## Tier 3 Optimization Applicability

| Optimization | Target | Actual Bottleneck? | Proceed? |
|--------------|--------|-------------------|----------|
| #10 SIMD score accumulation | getscore_Arr | [Yes/No] | [Yes/No] |
| #11 Pre-computed bounds | assignScores | [Yes/No] | [Yes/No] |
| #12 A*/IDA* search | Tree traversal | [Yes/No] | [Yes/No] |

## Files Generated

- `/tmp/profiling/flamegraph/*.svg` - Flame graphs
- `/tmp/profiling/perf/*.txt` - Raw profiling data
- `/tmp/profiling/reports/*.md` - Phase summaries
```

### Success Criteria

#### Automated Verification:
- [ ] All phase summary files exist in `/tmp/profiling/reports/`
- [ ] Final report exists: `test -s /tmp/profiling/reports/FINAL_PROFILING_REPORT.md`

#### Manual Verification:
- [ ] All open questions from research document answered
- [ ] Clear optimization priority established
- [ ] Decision made on whether to proceed with Tier 3 optimizations

---

## Testing Strategy

### Datasets to Profile

1. **Small (single_tree)**: Baseline for quick iteration
2. **Medium (16S_Bacteria subset)**: Representative production workload
3. **Large (full production)**: Scale-dependent bottlenecks

### Commands for Larger Dataset (if available)

```bash
# Replace [LARGE_REFERENCE] and [LARGE_QUERY] with actual paths
perf record -g --call-graph dwarf \
    -o /tmp/profiling/perf/large_dataset.data \
    ./tronko-assign -r -C 4 \
    -f [LARGE_REFERENCE] \
    -a [ALIGNMENT_FILE] \
    -p -1 [PAIRED_F] -2 [PAIRED_R] \
    -o /tmp/profiling/large_results.txt
```

---

## Performance Considerations

- **perf overhead**: ~5-10% slowdown during profiling is normal
- **DWARF call graphs**: Higher overhead but accurate stack traces
- **Lock tracing**: May have significant overhead; run separately
- **Flame graphs**: Generate after profiling completes, not during

---

## References

- Original profiling plan: `thoughts/shared/plans/2026-01-03-cpu-profiling-bottleneck-analysis.md`
- Gap analysis research: `thoughts/shared/research/2026-01-03-profiling-plan-gap-analysis.md`
- Optimization tiers: `thoughts/shared/research/2026-01-01-optimization-prioritization-matrix.md`
- SIMD opportunities: `thoughts/shared/research/2026-01-01-simd-vectorization-tronko-assign.md`
- Experiments log: `EXPERIMENTS_LOG.md`
