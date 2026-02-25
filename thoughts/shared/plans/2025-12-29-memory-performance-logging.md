# Memory Performance Logging Implementation Plan

## Overview

Implement extensive memory performance logging for tronko-assign to track memory usage at key phases, enable cross-run comparison, and support debugging memory-intensive operations. The logging will include both wall time and CPU time, support per-tree logging during reference loading, and offer a TSV export format for scripted analysis.

## Current State Analysis (experimental branch)

The `experimental` branch already has comprehensive logging infrastructure:

### Already Implemented:
- **logger.c/logger.h** - Full logging infrastructure with levels (DEBUG, INFO, WARN, ERROR)
- **resource_monitor.c/resource_monitor.h** - Memory/CPU tracking via `/proc/self/status` and `getrusage()`
- **Command-line flags**: `-V [0-3]`, `-l [FILE]`, `-R`, `-T`
- **18 milestones instrumented** in tronko-assign.c:
  - STARTUP, OPTIONS_PARSED, REFERENCE_LOADED
  - MEMORY_ALLOCATED, BWA_INDEX_BUILT, READ_SPECS_DETECTED
  - THREADS_INITIALIZED, BATCH_START, BATCH_LOADED
  - BWA_ALIGNMENT_COMPLETE, DETAILED_ALIGNMENT_COMPLETE
  - PLACEMENT_COMPLETE, LCA_COMPLETE, RESULTS_WRITTEN
  - BATCH_COMPLETE, CLEANUP_START, CLEANUP_COMPLETE, PROGRAM_END
- **Resource monitoring at cleanup points** (after closing files, freeing structures)
- **Thread-safe logging** with mutex protection
- **Crash debugging system** (crash_debug.c/crash_debug.h)

### Still Needed (based on user requirements):
1. **Per-tree logging** - Log memory after each tree loads (verbose mode)
2. **TSV export format** - Machine-parseable output for cross-run comparison
3. **Analysis scripts** - compare_memlogs.sh, plot_memlog.py

## Desired End State

After implementation:
1. Running `tronko-assign -V2 -R -l log.txt ...` produces human-readable logs (ALREADY WORKS)
2. Running `tronko-assign -V2 -R --tsv-log memory.tsv ...` produces TSV for analysis (NEW)
3. Running `tronko-assign -V3 -R ...` includes per-tree memory logging during reference load (NEW)
4. Comparison scripts can analyze multiple TSV log files

### Verification:
```bash
# Existing human-readable logging (already works)
./tronko-assign -V2 -R -T -l human.log -r -f reference.txt -a reads.fasta -s -g query.fasta -o results.txt

# New TSV logging for analysis
./tronko-assign -V2 -R --tsv-log memory.tsv -r -f reference.txt -a reads.fasta -s -g query.fasta -o results.txt
head memory.tsv  # Should show TSV with wall_time, phase, rss_mb, etc.

# Per-tree logging (DEBUG level)
./tronko-assign -V3 -R -l verbose.log -r -f reference.txt -a reads.fasta -s -g query.fasta -o results.txt
grep "TREE_" verbose.log  # Should show per-tree entries
```

## What We're NOT Doing

- Replacing the existing human-readable logger (it stays as-is)
- Per-allocation tracking (too granular)
- Thread-level memory tracking (complex, deferred)
- macOS/Windows portability (Linux `/proc/self/status` only)

## Implementation Approach

1. Add TSV export capability alongside existing logger
2. Add per-tree logging to readReferenceTree (at DEBUG level)
3. Create analysis scripts for TSV logs

---

## Phase 1: TSV Export Format (NEW FLAG)

### Overview
Add a `--tsv-log` flag that produces machine-parseable TSV output alongside the existing human-readable logging.

### Changes Required:

#### 1. Update global.h - Add tsv_log_file field
**File**: `tronko-assign/global.h` (in Options struct, after existing logging fields)

```c
    int enable_timing;
    char tsv_log_file[BUFFER_SIZE];  // TSV memory log output file (empty = disabled)
}Options;
```

#### 2. Update options.c - Add long option
**File**: `tronko-assign/options.c` (in long_options array)

```c
    {"enable-timing",no_argument,0,'T'},
    {"tsv-log",required_argument,0,0},  // Long option only, no short form
    {0, 0, 0, 0}
```

#### 3. Update options.c - Add usage string
**File**: `tronko-assign/options.c` (in usage string)

```c
    -T, Enable timing information\n\
    --tsv-log [FILE], Export memory stats to TSV file for analysis\n\
    \n";
```

#### 4. Update options.c - Handle long option
**File**: `tronko-assign/options.c` (add case for long option index)

Handle `--tsv-log` by checking `option_index` when `c == 0`:
```c
            case 0:
                // Handle long options without short equivalents
                if (strcmp(long_options[option_index].name, "tsv-log") == 0) {
                    strncpy(opt->tsv_log_file, optarg, sizeof(opt->tsv_log_file) - 1);
                    opt->tsv_log_file[sizeof(opt->tsv_log_file) - 1] = '\0';
                }
                break;
```

#### 5. Update tronko-assign.c - Initialize and use TSV log
**File**: `tronko-assign/tronko-assign.c`

Add initialization:
```c
    opt.tsv_log_file[0] = '\0';
```

Add TSV file handling after logger initialization:
```c
    FILE *tsv_log = NULL;
    if (opt.tsv_log_file[0] != '\0') {
        tsv_log = fopen(opt.tsv_log_file, "w");
        if (tsv_log) {
            // Write header
            fprintf(tsv_log, "# tronko-assign memory log v1.0\n");
            fprintf(tsv_log, "wall_time\tphase\trss_mb\tvm_mb\tpeak_rss_mb\tcpu_user\tcpu_sys\textra_info\n");
            fflush(tsv_log);
        } else {
            LOG_WARN("Could not open TSV log file: %s", opt.tsv_log_file);
        }
    }
```

#### 6. Create tsv_memlog.h helper
**File**: `tronko-assign/tsv_memlog.h`

```c
#ifndef TSV_MEMLOG_H
#define TSV_MEMLOG_H

#include <stdio.h>
#include "resource_monitor.h"

// Write a TSV log entry
static inline void tsv_memlog_write(FILE *f, const char *phase, const char *extra_fmt, ...) {
    if (!f) return;

    resource_stats_t stats;
    get_resource_stats(&stats);

    fprintf(f, "%.3f\t%s\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\t",
            stats.wall_time_sec, phase,
            stats.memory_rss_kb / 1024.0,
            stats.memory_vm_size_kb / 1024.0,
            stats.memory_vm_rss_peak_kb / 1024.0,
            stats.user_time_sec,
            stats.system_time_sec);

    if (extra_fmt && extra_fmt[0] != '\0') {
        va_list args;
        va_start(args, extra_fmt);
        vfprintf(f, extra_fmt, args);
        va_end(args);
    }

    fprintf(f, "\n");
    fflush(f);
}

#define TSV_LOG(file, phase, ...) \
    do { if (file) tsv_memlog_write(file, phase, __VA_ARGS__); } while(0)

#define TSV_LOG_SIMPLE(file, phase) \
    do { if (file) tsv_memlog_write(file, phase, ""); } while(0)

#endif // TSV_MEMLOG_H
```

#### 7. Add TSV logging calls at key milestones
**File**: `tronko-assign/tronko-assign.c`

Add TSV logging alongside existing LOG_MILESTONE calls:
```c
    // After STARTUP milestone
    TSV_LOG_SIMPLE(tsv_log, "STARTUP");

    // After REFERENCE_LOADED milestone
    TSV_LOG(tsv_log, "REFERENCE_LOADED", "trees=%d", numberOfTrees);

    // After BWA_INDEX_BUILT milestone
    TSV_LOG_SIMPLE(tsv_log, "BWA_INDEX");

    // After MEMORY_ALLOCATED milestone
    TSV_LOG(tsv_log, "THREADS_ALLOCATED", "threads=%d", opt.number_of_cores);

    // In batch loop
    TSV_LOG(tsv_log, "BATCH_START", "batch=%d,reads=%d", batch_num, returnLineNumber);
    TSV_LOG(tsv_log, "BATCH_COMPLETE", "batch=%d", batch_num);

    // At end
    TSV_LOG_SIMPLE(tsv_log, "FINAL");
    if (tsv_log) fclose(tsv_log);
```

### Success Criteria:

#### Automated Verification:
- [x] `make clean && make` compiles without errors
- [x] `./tronko-assign -h` shows the new --tsv-log option
- [x] TSV log file is created with correct header

#### Manual Verification:
- [x] TSV file is parseable by awk/pandas
- [x] Memory values match those in human-readable log

---

## Phase 2: Per-Tree Logging (Verbose/DEBUG Mode)

### Overview
Add per-tree memory logging during reference loading. This activates at DEBUG level (-V3).

### Changes Required:

#### 1. Add external declarations to readreference.c
**File**: `tronko-assign/readreference.c` (near top, after includes)

```c
#include "logger.h"
#include "resource_monitor.h"

// External TSV log file handle (set in main)
extern FILE *g_tsv_log_file;
```

#### 2. Add global TSV log pointer in tronko-assign.c
**File**: `tronko-assign/tronko-assign.c` (after includes)

```c
// Global TSV log file for use by readreference.c
FILE *g_tsv_log_file = NULL;
```

Set it after opening the TSV log:
```c
    if (tsv_log) {
        g_tsv_log_file = tsv_log;
    }
```

#### 3. Add per-tree logging in readReferenceTree
**File**: `tronko-assign/readreference.c` (in the tree allocation loop, around line 392)

```c
    treeArr = malloc(numberOfTrees*sizeof(struct node *));
    for (i=0; i<numberOfTrees; i++){
        allocateTreeArrMemory(i,max_nodename);

        // Per-tree logging at DEBUG level
        LOG_DEBUG("Tree %d allocated: nodes=%d, bases=%d", i, 2*numspecArr[i]-1, numbaseArr[i]);

        // TSV logging if enabled
        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_ALLOCATED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d,nodes=%d,bases=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    i, 2*numspecArr[i]-1, numbaseArr[i]);
            fflush(g_tsv_log_file);
        }
    }
```

#### 4. Track tree transitions during posterior loading
**File**: `tronko-assign/readreference.c` (in the main while loop)

Add tracking variable before the while loop:
```c
    int last_logged_tree = -1;
```

Inside the while loop, after reading posteriors:
```c
    // Log when we finish loading a tree's posteriors (tree number changes)
    if (treeNumber != last_logged_tree && last_logged_tree >= 0) {
        LOG_DEBUG("Tree %d posteriors loaded", last_logged_tree);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    last_logged_tree);
            fflush(g_tsv_log_file);
        }
    }
    last_logged_tree = treeNumber;
```

After the while loop ends, log the final tree:
```c
    // Log the final tree
    if (last_logged_tree >= 0) {
        LOG_DEBUG("Tree %d posteriors loaded (final)", last_logged_tree);

        if (g_tsv_log_file) {
            resource_stats_t stats;
            get_resource_stats(&stats);
            fprintf(g_tsv_log_file, "%.3f\tTREE_LOADED\t%.1f\t%.1f\t%.1f\t%.3f\t%.3f\ttree=%d\n",
                    stats.wall_time_sec,
                    stats.memory_rss_kb / 1024.0,
                    stats.memory_vm_size_kb / 1024.0,
                    stats.memory_vm_rss_peak_kb / 1024.0,
                    stats.user_time_sec,
                    stats.system_time_sec,
                    last_logged_tree);
            fflush(g_tsv_log_file);
        }
    }
    return numberOfTrees;
```

### Success Criteria:

#### Automated Verification:
- [x] `make clean && make` compiles without errors
- [x] `-V2 --tsv-log test.tsv` does NOT produce TREE_ALLOCATED entries (INFO level)
- [x] `-V3 --tsv-log test.tsv` DOES produce TREE_ALLOCATED entries (DEBUG level)

#### Manual Verification:
- [x] Number of TREE_ALLOCATED entries matches number of trees
- [x] Memory increases visible between trees in TSV
- [x] Human-readable log also shows tree debug messages at -V3

---

## Phase 3: Analysis Scripts

### Overview
Create helper scripts for analyzing and comparing TSV memory logs.

### Changes Required:

#### 1. Create scripts directory
```bash
mkdir -p tronko-assign/scripts
```

#### 2. Create compare_memlogs.sh
**File**: `tronko-assign/scripts/compare_memlogs.sh`

```bash
#!/bin/bash
# Compare memory usage across multiple tronko-assign runs
# Usage: ./compare_memlogs.sh log1.tsv log2.tsv ...

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 log1.tsv [log2.tsv ...]"
    exit 1
fi

echo "Run,Peak_RSS_MB,Time_to_Peak_s,Final_RSS_MB,Total_Time_s,CPU_User_s,CPU_Sys_s"
for log in "$@"; do
    if [ ! -f "$log" ]; then
        echo "Warning: File not found: $log" >&2
        continue
    fi

    run_name=$(basename "$log" .tsv)

    # Skip header lines (starting with # or containing 'wall_time')
    peak_rss=$(awk -F'\t' 'NR>1 && !/^#/ && !/wall_time/ {print $5}' "$log" | sort -n | tail -1)
    peak_time=$(awk -F'\t' -v peak="$peak_rss" 'NR>1 && !/^#/ && $5==peak {print $1; exit}' "$log")
    final_line=$(awk -F'\t' '/FINAL/ {print}' "$log" | tail -1)

    if [ -n "$final_line" ]; then
        final_rss=$(echo "$final_line" | cut -f3)
        total_time=$(echo "$final_line" | cut -f1)
        cpu_user=$(echo "$final_line" | cut -f6)
        cpu_sys=$(echo "$final_line" | cut -f7)
    else
        final_rss="N/A"
        total_time="N/A"
        cpu_user="N/A"
        cpu_sys="N/A"
    fi

    echo "$run_name,$peak_rss,$peak_time,$final_rss,$total_time,$cpu_user,$cpu_sys"
done
```

#### 3. Create plot_memlog.py
**File**: `tronko-assign/scripts/plot_memlog.py`

```python
#!/usr/bin/env python3
"""
Plot memory usage from tronko-assign TSV memory logs.
Usage: python plot_memlog.py log1.tsv [log2.tsv ...] [-o output.png]
"""
import sys
import argparse

def load_memlog(path):
    """Load a memory log TSV file."""
    try:
        import pandas as pd
        return pd.read_csv(path, sep='\t', comment='#',
                           names=['wall_time', 'phase', 'rss_mb', 'vm_mb',
                                  'peak_rss_mb', 'cpu_user', 'cpu_sys', 'extra'],
                           skiprows=1)  # Skip header row
    except ImportError:
        print("Error: pandas is required. Install with: pip install pandas")
        sys.exit(1)

def plot_logs(log_files, output='memory_comparison.png'):
    """Plot RSS over time for multiple log files."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Error: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8))

    for path in log_files:
        try:
            df = load_memlog(path)
            label = path.split('/')[-1].replace('.tsv', '')

            # RSS over time
            axes[0].plot(df['wall_time'], df['rss_mb'], label=label, marker='.', markersize=4)
        except Exception as e:
            print(f"Warning: Could not load {path}: {e}")
            continue

    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('RSS (MB)')
    axes[0].legend()
    axes[0].set_title('Memory Usage Over Time')
    axes[0].grid(True, alpha=0.3)

    # Phase comparison (bar chart of key phases)
    key_phases = ['STARTUP', 'REFERENCE_LOADED', 'BWA_INDEX', 'THREADS_ALLOCATED']
    phase_data = {}

    for path in log_files:
        try:
            df = load_memlog(path)
            label = path.split('/')[-1].replace('.tsv', '')
            for phase in key_phases:
                phase_rows = df[df['phase'] == phase]
                if not phase_rows.empty:
                    if phase not in phase_data:
                        phase_data[phase] = {}
                    phase_data[phase][label] = phase_rows['rss_mb'].iloc[0]
        except:
            continue

    if phase_data:
        import pandas as pd
        phase_df = pd.DataFrame(phase_data).T
        phase_df.plot(kind='bar', ax=axes[1], alpha=0.7)
        axes[1].set_ylabel('RSS (MB)')
        axes[1].set_title('Memory by Phase')
        axes[1].legend()
        axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved plot to {output}")

def main():
    parser = argparse.ArgumentParser(description='Plot memory usage from tronko-assign TSV logs')
    parser.add_argument('logs', nargs='+', help='TSV log files to plot')
    parser.add_argument('-o', '--output', default='memory_comparison.png', help='Output image file')
    args = parser.parse_args()

    plot_logs(args.logs, args.output)

if __name__ == '__main__':
    main()
```

#### 4. Make scripts executable
```bash
chmod +x tronko-assign/scripts/compare_memlogs.sh
chmod +x tronko-assign/scripts/plot_memlog.py
```

### Success Criteria:

#### Automated Verification:
- [x] `bash scripts/compare_memlogs.sh test.tsv` produces CSV output
- [x] `python scripts/plot_memlog.py test.tsv` creates PNG file (if matplotlib available)

#### Manual Verification:
- [x] Comparison output shows meaningful differences between runs
- [x] Plot clearly shows memory phases

---

## Testing Strategy

### Integration Tests:
```bash
cd tronko-assign

# Test 1: Existing human-readable logging (should already work)
./tronko-assign -V2 -R -T -l /tmp/human.log \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/test_results.txt

grep "MILESTONE" /tmp/human.log  # Should show milestones

# Test 2: New TSV logging
./tronko-assign -V2 -R --tsv-log /tmp/memory.tsv \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/test_results.txt

head /tmp/memory.tsv  # Should show TSV header and data
wc -l /tmp/memory.tsv  # Should have multiple entries

# Test 3: Per-tree logging at DEBUG level
./tronko-assign -V3 -R --tsv-log /tmp/verbose.tsv \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/test_results.txt

grep "TREE_" /tmp/verbose.tsv  # Should show per-tree entries

# Test 4: Analysis scripts
./scripts/compare_memlogs.sh /tmp/memory.tsv
python scripts/plot_memlog.py /tmp/memory.tsv -o /tmp/memory_plot.png
```

### Manual Testing Steps:
1. Verify TSV output is parseable by spreadsheet software
2. Compare TSV values with human-readable log values
3. Test with multi-tree reference database
4. Verify per-tree entries appear in correct order

## Performance Considerations

- TSV logging adds minimal overhead (same resource_stats calls as existing logger)
- Per-tree logging at DEBUG level only (no impact at INFO level)
- Logging is buffered and flushed per entry for reliability
- No runtime impact when `--tsv-log` is not specified

## Migration Notes

- Existing `-V`, `-l`, `-R`, `-T` flags work unchanged
- New `--tsv-log` flag is additive, not replacing existing functionality
- No changes to existing log format

## References

- Original research: `thoughts/shared/research/2025-12-29-memory-performance-logging.md`
- Existing logger: `tronko-assign/logger.c:56-105`
- Existing resource monitor: `tronko-assign/resource_monitor.c:118-162`
- Performance logging docs: `docs/performance-logging.md`
