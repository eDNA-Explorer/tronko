# CPU Profiling Plan: Finding the Real Bottleneck

## Context

Tier 2 optimizations (mmap, SIMD max-finding, two-phase screening) had no measurable performance impact. This indicates our assumptions about where time is spent were wrong. Before investing in Tier 3 optimizations (#10 SIMD score accumulation, #11 pre-computed bounds, #12 A*/IDA* search), we need actual profiling data.

## Current Tracking vs. What We Need

### What Memory/Performance Logging Already Tracks

From `thoughts/shared/plans/2025-12-29-memory-performance-logging.md`:

| Metric | Granularity | Purpose |
|--------|-------------|---------|
| RSS/VM memory | Per-milestone, per-tree | Memory growth tracking |
| Wall time | Per-milestone | Phase duration |
| CPU time (user+sys) | Per-milestone | Total CPU consumed |
| Peak RSS | Cumulative | Memory high-water mark |

**Milestones tracked:**
- STARTUP → OPTIONS_PARSED → REFERENCE_LOADED → BWA_INDEX_BUILT
- THREADS_INITIALIZED → BATCH_START → BWA_ALIGNMENT_COMPLETE
- DETAILED_ALIGNMENT_COMPLETE → PLACEMENT_COMPLETE → LCA_COMPLETE
- RESULTS_WRITTEN → BATCH_COMPLETE → CLEANUP → PROGRAM_END

### What Current Logging Does NOT Tell Us

| Gap | Why It Matters |
|-----|----------------|
| **Which functions consume CPU** | Milestones show phase duration, not internal breakdown |
| **Hot loops within functions** | Can't see if time is in `getscore_Arr` vs `WFA2_align` vs BWA |
| **Cache miss rates** | Memory access patterns may be the real bottleneck |
| **Branch misprediction** | Character branching in scoring could cause pipeline stalls |
| **Lock contention** | Multi-threaded sections may have hidden synchronization costs |
| **I/O wait vs compute** | Can't distinguish blocked I/O from CPU-bound work |

---

## Profiling Approach

### Phase 1: CPU Time Distribution (perf stat)

**Goal:** Understand high-level CPU characteristics

```bash
cd tronko-assign

# Basic CPU counters
perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
    ./tronko-assign -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/test_results.txt
```

**What to look for:**
- IPC (instructions per cycle) < 1.0 suggests memory-bound or branch-bound code
- Cache miss rate > 5% suggests memory access pattern issues
- Branch miss rate > 2% suggests unpredictable branches (character comparisons?)

### Phase 2: Function-Level Hotspots (perf record + report)

**Goal:** Identify which functions consume the most CPU time

```bash
# Record with call graph
perf record -g --call-graph dwarf \
    ./tronko-assign -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/test_results.txt

# View hotspots
perf report --sort=dso,symbol
```

**Expected candidates (ranked by Tier 2 assumptions):**
1. `getscore_Arr` / `getscore_Arr_ncbi` - score accumulation
2. `assignScores_Arr_paired` / `assignScores_Arr_single` - tree traversal
3. `wavefront_align` (WFA2) - detailed alignment
4. `bwa_*` functions - BWA alignment
5. `findMax` / max-finding loops in placement.c

**If the actual hotspot is different, our optimization targets were wrong.**

### Phase 3: Flame Graph Visualization

**Goal:** Visualize call stack to see where time accumulates

```bash
# Generate flame graph (requires FlameGraph tools)
perf script | stackcollapse-perf.pl | flamegraph.pl > tronko_flame.svg
```

**What to look for:**
- Wide bars = functions consuming significant time
- Deep stacks = many levels of function calls (overhead?)
- Unexpected wide bars = functions we didn't consider optimizing

### Phase 4: Line-Level Annotation (perf annotate)

**Goal:** Find hot lines within identified functions

```bash
# After perf record, annotate specific function
perf annotate getscore_Arr
perf annotate assignScores_Arr_paired
```

**What to look for:**
- Specific loops consuming disproportionate time
- Memory access instructions with high cycle counts (cache misses)
- Branch instructions with high counts (mispredictions)

### Phase 5: Memory Access Analysis (optional, if CPU profile suggests memory-bound)

```bash
# L1/L2/L3 cache behavior
perf stat -e L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
    ./tronko-assign [args...]

# Memory bandwidth
perf stat -e mem_load_retired.l3_miss,mem_load_retired.fb_hit \
    ./tronko-assign [args...]
```

---

## Hypotheses to Test

Based on Tier 2's failure, here are revised hypotheses:

| Hypothesis | How to Confirm | Implication |
|------------|----------------|-------------|
| **BWA alignment dominates** | >50% time in `bwa_*` functions | Tier 3 won't help; need BWA alternatives |
| **WFA2 alignment dominates** | >50% time in `wavefront_*` | Consider WFA-GPU or simpler alignment |
| **Memory allocation churn** | Significant time in `malloc`/`free` | Pre-allocate buffers, use arenas |
| **String parsing overhead** | Time in `sscanf`/`strtok`/string ops | Optimize input parsing |
| **Lock contention** | Time in `pthread_mutex_*` | Reduce synchronization |
| **I/O blocking** | Low CPU utilization during phases | Async I/O, prefetching |

---

## Integration with Existing Logging

### Use TSV logs to correlate with profiling

```bash
# Run with both TSV logging and perf
perf record -g --call-graph dwarf \
    ./tronko-assign -V2 -R --tsv-log /tmp/profile_run.tsv \
    -r -f reference.txt -a reads.fasta -s -g query.fasta -o /tmp/out.txt

# TSV shows wall time per phase
# perf shows CPU breakdown within phases
```

### Add phase markers to perf (optional enhancement)

Could add `perf_event_open` markers at milestone boundaries to correlate perf data with phases. This is optional - manual correlation usually sufficient.

---

## Decision Matrix After Profiling

| If Hotspot Is... | Then Do... | Skip... |
|------------------|------------|---------|
| `getscore_Arr` | Tier 3 #10 (SIMD score) | - |
| `assignScores_Arr_*` | Tier 3 #12 (A*/IDA*) | - |
| BWA functions | Explore BarraCUDA/GPU-BWA | All Tier 3 |
| WFA2 functions | WFA-GPU or simpler aligner | Tier 3 #10, #12 |
| `malloc`/`free` | Memory pooling | All Tier 3 |
| `pthread_*` | Lock-free structures | - |
| I/O functions | Async I/O, larger buffers | All Tier 3 |

---

## Test Datasets

Profile with representative workloads:

1. **Small (baseline):** `single_tree` example dataset
2. **Medium:** 16S_Bacteria dataset (if available)
3. **Large:** Production-scale reference database

Performance characteristics may differ by scale - profile all three.

---

## Success Criteria

Profiling is complete when we can answer:

1. [ ] Which function(s) consume >50% of CPU time?
2. [ ] Is the code CPU-bound, memory-bound, or I/O-bound?
3. [ ] Are cache miss rates or branch mispredictions significant?
4. [ ] Does the hotspot change between small/medium/large datasets?
5. [ ] Which specific Tier 3 optimization (if any) targets the actual bottleneck?

---

## Next Steps

1. Run Phase 1 (`perf stat`) to get baseline CPU characteristics
2. Run Phase 2 (`perf record`) to identify function hotspots
3. Generate flame graph for visualization
4. Compare hotspots against Tier 3 optimization targets
5. Update optimization priority based on actual data
