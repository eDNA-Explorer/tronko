---
date: 2025-12-29T12:00:00-08:00
researcher: Claude
git_commit: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
branch: main
repository: tronko
topic: "Memory Performance Logging Strategy for tronko-assign"
tags: [research, codebase, tronko-assign, performance, memory, logging]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
---

# Research: Memory Performance Logging Strategy for tronko-assign

**Date**: 2025-12-29
**Researcher**: Claude
**Git Commit**: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
**Branch**: main
**Repository**: tronko

## Research Question
Based on the data flow analysis, where should we add memory performance logging to track usage during a run? How can we output logs for cross-run comparison?

## Summary

Memory logging should be added at **7 key instrumentation points** that correspond to major memory state changes. The recommended approach uses Linux's `/proc/self/status` for accurate RSS tracking, outputs to a separate TSV log file (controlled by a new `-M` flag), and includes timestamps for timeline analysis.

## Recommended Instrumentation Points

### Phase Diagram with Logging Points

```
┌────────────────────────────────────────────────────────────────────────────┐
│                    MEMORY LOGGING INSTRUMENTATION POINTS                    │
└────────────────────────────────────────────────────────────────────────────┘

[1] STARTUP ─────────────────────────────────────────────────────────────────
    │ Log: baseline memory before any allocations
    │ Location: tronko-assign.c:744 (start of main)
    ▼
[2] AFTER REFERENCE LOAD ────────────────────────────────────────────────────
    │ Log: memory after loading all trees + posteriors (BIGGEST JUMP)
    │ Location: tronko-assign.c:788 (after readReferenceTree returns)
    │
    │ Key metrics to capture:
    │   - Total posteriors loaded
    │   - Number of trees, total nodes
    │   - Taxonomy array size
    ▼
[3] AFTER BWA INDEX ─────────────────────────────────────────────────────────
    │ Log: memory after BWA FM-index construction
    │ Location: tronko-assign.c:872 (single) / 1079 (paired)
    ▼
[4] AFTER PER-THREAD ALLOCATION ─────────────────────────────────────────────
    │ Log: memory after resultsStruct allocation for all threads
    │ Location: tronko-assign.c:934 (after allocateMemForResults loop)
    │
    │ Key metrics:
    │   - Number of threads
    │   - nodeScores size per thread
    ▼
[5] BATCH LOOP START (per batch) ────────────────────────────────────────────
    │ Log: memory at start of each batch iteration
    │ Location: tronko-assign.c:937 (inside while loop, after readInXLines)
    │
    │ Key metrics:
    │   - Batch number
    │   - Reads in batch
    ▼
[6] PEAK BATCH (per batch) ──────────────────────────────────────────────────
    │ Log: memory after thread execution completes
    │ Location: tronko-assign.c:970 (after pthread_join loop)
    │
    │ Captures peak per-batch memory including:
    │   - taxonPath arrays
    │   - Thread-local bwa_results
    ▼
[7] FINAL CLEANUP ───────────────────────────────────────────────────────────
    │ Log: memory after all cleanup
    │ Location: tronko-assign.c:1205 (end of main)
```

## Implementation Design

### 1. Memory Measurement Function

Add to a new file `memlog.c` or directly in `tronko-assign.c`:

```c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    long rss_kb;        // Resident Set Size in KB
    long vm_size_kb;    // Virtual Memory size in KB
    long vm_peak_kb;    // Peak Virtual Memory
    long rss_peak_kb;   // Peak RSS (VmHWM)
    double timestamp;   // Seconds since program start
} MemorySnapshot;

static struct timespec program_start;
static int memlog_initialized = 0;

void memlog_init() {
    clock_gettime(CLOCK_MONOTONIC, &program_start);
    memlog_initialized = 1;
}

MemorySnapshot get_memory_snapshot() {
    MemorySnapshot snap = {0};

    // Get elapsed time
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    snap.timestamp = (now.tv_sec - program_start.tv_sec) +
                     (now.tv_nsec - program_start.tv_nsec) * 1e-9;

    // Read /proc/self/status for memory info
    FILE* status = fopen("/proc/self/status", "r");
    if (status) {
        char line[256];
        while (fgets(line, sizeof(line), status)) {
            if (strncmp(line, "VmRSS:", 6) == 0)
                sscanf(line + 6, "%ld", &snap.rss_kb);
            else if (strncmp(line, "VmSize:", 7) == 0)
                sscanf(line + 7, "%ld", &snap.vm_size_kb);
            else if (strncmp(line, "VmPeak:", 7) == 0)
                sscanf(line + 7, "%ld", &snap.vm_peak_kb);
            else if (strncmp(line, "VmHWM:", 6) == 0)
                sscanf(line + 6, "%ld", &snap.rss_peak_kb);
        }
        fclose(status);
    }
    return snap;
}
```

### 2. Logging Output Format

**TSV format for easy parsing and comparison:**

```
# tronko-assign memory log
# version: 1.0
# reference_file: /path/to/reference.txt
# num_trees: 100
# num_cores: 8
# batch_size: 50000
# timestamp	phase	rss_mb	vm_mb	peak_rss_mb	extra_info
0.000	STARTUP	12.5	45.2	12.5
0.234	REFERENCE_LOADED	9876.3	10234.5	9876.3	trees=100,nodes=198000
0.456	BWA_INDEX	10234.5	10890.2	10234.5
0.512	THREADS_ALLOCATED	10456.7	11234.5	10456.7	threads=8
0.634	BATCH_START	10456.7	11234.5	10456.7	batch=1,reads=25000
1.234	BATCH_COMPLETE	10567.8	11345.6	10567.8	batch=1
1.256	BATCH_START	10456.7	11234.5	10567.8	batch=2,reads=25000
...
45.678	FINAL	234.5	456.7	10567.8	total_reads=1000000
```

### 3. Logging Macro

```c
#define MEMLOG(file, phase, extra_fmt, ...) \
    do { \
        if (file) { \
            MemorySnapshot _snap = get_memory_snapshot(); \
            fprintf(file, "%.3f\t%s\t%.1f\t%.1f\t%.1f\t" extra_fmt "\n", \
                _snap.timestamp, phase, \
                _snap.rss_kb / 1024.0, \
                _snap.vm_size_kb / 1024.0, \
                _snap.rss_peak_kb / 1024.0, \
                ##__VA_ARGS__); \
            fflush(file); \
        } \
    } while(0)
```

### 4. Command-Line Flag Addition

Add to `options.c` and `options.h`:

```c
// In options.h, add to struct options:
char *memlog_file;      // Memory log output file (NULL = disabled)

// In options.c, add case:
case 'M':
    opt->memlog_file = optarg;
    break;
```

## Specific Code Insertion Points

### Point 1: Startup
**File**: `tronko-assign.c:744`
```c
int main(int argc, char* argv[]) {
    // ... existing declarations ...
    FILE* memlog = NULL;

    memlog_init();  // Initialize timing baseline

    parseOptions(argc, argv, &opt);

    if (opt.memlog_file) {
        memlog = fopen(opt.memlog_file, "w");
        // Write header
        fprintf(memlog, "# tronko-assign memory log v1.0\n");
        fprintf(memlog, "timestamp\tphase\trss_mb\tvm_mb\tpeak_rss_mb\textra_info\n");
        MEMLOG(memlog, "STARTUP", "");
    }
```

### Point 2: After Reference Load
**File**: `tronko-assign.c:788` (after readReferenceTree call)
```c
readReferenceTree(/*...*/);

MEMLOG(memlog, "REFERENCE_LOADED", "trees=%d,total_nodes=%d",
       numberOfTrees, number_of_total_nodes);
```

### Point 3: After BWA Index
**File**: `tronko-assign.c:872` (single-end) and `1079` (paired-end)
```c
bwa_idx = bwa_idx_load_new(/*...*/);

MEMLOG(memlog, "BWA_INDEX", "");
```

### Point 4: After Per-Thread Allocation
**File**: `tronko-assign.c:934` (after the allocateMemForResults loop)
```c
for (i=0; i<opt.number_of_cores; i++) {
    mstr[i].str = malloc(sizeof(struct resultsStruct));
    allocateMemForResults(/*...*/);
    // ... rest of setup ...
}

MEMLOG(memlog, "THREADS_ALLOCATED", "threads=%d", opt.number_of_cores);
```

### Point 5 & 6: Batch Loop
**File**: `tronko-assign.c:937` and `970`
```c
int batch_num = 0;
while (1) {
    returnLineNumber = readInXNumberOfLines(/*...*/);
    if (returnLineNumber == 0) break;
    batch_num++;

    MEMLOG(memlog, "BATCH_START", "batch=%d,reads=%d",
           batch_num, returnLineNumber);

    // ... thread setup and execution ...

    for (i=0; i<opt.number_of_cores; i++) {
        pthread_join(threads[i], NULL);
    }

    MEMLOG(memlog, "BATCH_COMPLETE", "batch=%d", batch_num);

    // ... result collection and cleanup ...
}
```

### Point 7: Final Cleanup
**File**: `tronko-assign.c:1205` (end of main)
```c
// After all cleanup...
MEMLOG(memlog, "FINAL", "");

if (memlog) fclose(memlog);
return 0;
```

## Comparison Analysis Scripts

### Basic Comparison Script (bash)

```bash
#!/bin/bash
# compare_memlogs.sh - Compare memory usage across runs

echo "Run,Peak_RSS_MB,Time_to_Peak_s,Final_RSS_MB,Total_Time_s"
for log in "$@"; do
    run_name=$(basename "$log" .memlog)
    peak_rss=$(awk -F'\t' 'NR>2 {print $5}' "$log" | sort -n | tail -1)
    peak_time=$(awk -F'\t' -v peak="$peak_rss" 'NR>2 && $5==peak {print $1; exit}' "$log")
    final_rss=$(awk -F'\t' '/FINAL/ {print $3}' "$log")
    total_time=$(awk -F'\t' '/FINAL/ {print $1}' "$log")
    echo "$run_name,$peak_rss,$peak_time,$final_rss,$total_time"
done
```

### Python Analysis Example

```python
import pandas as pd
import matplotlib.pyplot as plt

def load_memlog(path):
    return pd.read_csv(path, sep='\t', comment='#',
                       names=['timestamp', 'phase', 'rss_mb', 'vm_mb', 'peak_rss_mb', 'extra'])

def compare_runs(log_files):
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    for path in log_files:
        df = load_memlog(path)
        label = path.split('/')[-1].replace('.memlog', '')

        # RSS over time
        axes[0].plot(df['timestamp'], df['rss_mb'], label=label)

        # Phase comparison
        phases = df[df['phase'].str.contains('LOADED|INDEX|ALLOCATED')]
        axes[1].bar(phases['phase'], phases['rss_mb'], alpha=0.7, label=label)

    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('RSS (MB)')
    axes[0].legend()
    axes[0].set_title('Memory Usage Over Time')

    axes[1].set_ylabel('RSS (MB)')
    axes[1].legend()
    axes[1].set_title('Memory by Phase')

    plt.tight_layout()
    plt.savefig('memory_comparison.png')
```

## Expected Memory Profile

Based on the data flow analysis, a typical run should show:

```
Phase                    Expected RSS Change     Cumulative %
─────────────────────────────────────────────────────────────
STARTUP                  ~10-50 MB               <1%
REFERENCE_LOADED         +5-15 GB (dominant!)    80-95%
BWA_INDEX               +100-500 MB              +2-5%
THREADS_ALLOCATED       +50-500 MB               +1-5%
BATCH_START             +50-200 MB               +1-2%
BATCH_COMPLETE          Same or slightly lower
FINAL                   ~50-200 MB (cleanup)
```

## Additional Metrics to Consider

### Per-Tree Breakdown (Optional)
For debugging which trees consume most memory:

```c
// In readreference.c after each tree is loaded
MEMLOG(memlog, "TREE_LOADED", "tree=%d,nodes=%d,bases=%d,est_mb=%.1f",
       i, 2*numspecArr[i]-1, numbaseArr[i],
       (2*numspecArr[i]-1) * numbaseArr[i] * 32.0 / (1024*1024));
```

### Thread-Level Tracking (Advanced)
For per-thread memory in multi-threaded section:

```c
// Requires thread-local storage or pthread_getspecific
// More complex to implement but useful for thread scaling analysis
```

## Code References

- `tronko-assign/tronko-assign.c:744` - main() entry point
- `tronko-assign/tronko-assign.c:788` - After readReferenceTree
- `tronko-assign/tronko-assign.c:872` - BWA index (single-end)
- `tronko-assign/tronko-assign.c:1079` - BWA index (paired-end)
- `tronko-assign/tronko-assign.c:912-934` - Thread allocation loop
- `tronko-assign/tronko-assign.c:935-981` - Batch processing loop
- `tronko-assign/readreference.c:311-433` - Reference tree loading
- `tronko-assign/allocatetreememory.c:17-40` - Posterior allocation
- `tronko-assign/allocateMemoryForResults.c:3-89` - Per-thread results
- `tronko-assign/options.c:70-217` - Command-line parsing

## Architecture Insights

1. **Existing Infrastructure**: The codebase already uses `clock_gettime(CLOCK_MONOTONIC)` for timing (mostly commented out in placement.c), so the pattern is familiar

2. **BWA Verbose Pattern**: BWA uses a `bwa_verbose` global with levels 1-4. A similar `memlog_enabled` or passing the FILE* would fit the existing style

3. **No Memory Tracking Currently**: There's no existing memory measurement code - this would be new infrastructure

4. **Linux-Specific**: The `/proc/self/status` approach is Linux-specific. For portability, could add `getrusage(RUSAGE_SELF, &r)` fallback (gives `ru_maxrss`)

## Open Questions

1. Should per-tree memory be logged during reference loading (many log lines but useful for debugging)?
2. Should the memory log include CPU time alongside wall time?
3. Is JSON output preferable to TSV for some use cases?
4. Should there be a "verbose memory" mode that logs every allocation?
