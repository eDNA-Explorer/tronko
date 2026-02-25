---
date: 2025-12-30T16:45:00-08:00
researcher: Claude
git_commit: experimental
branch: experimental
repository: tronko
topic: "Memory Logging Usage Guide for tronko-assign"
tags: [documentation, tronko-assign, performance, memory, logging, tsv]
status: complete
last_updated: 2025-12-30
last_updated_by: Claude
---

# Memory Logging Usage Guide for tronko-assign

**Date**: 2025-12-30
**Branch**: experimental

## Overview

tronko-assign includes a TSV memory logging feature that tracks memory usage at key phases during execution. This enables cross-run comparison, performance debugging, and identifying memory-intensive operations.

## Quick Start

### Basic TSV Logging

```bash
./tronko-assign -V2 -R --tsv-log memory.tsv \
    -r -f reference_tree.txt \
    -a reference.fasta \
    -s -g query.fasta \
    -o results.txt
```

This produces a TSV file with memory statistics at each milestone.

### Verbose Per-Tree Logging

For detailed per-tree memory tracking during reference loading:

```bash
./tronko-assign -V3 -R --tsv-log verbose.tsv \
    -r -f reference_tree.txt \
    -a reference.fasta \
    -s -g query.fasta \
    -o results.txt
```

At DEBUG level (-V3), the log includes TREE_ALLOCATED and TREE_LOADED entries for each tree.

## Command-Line Options

| Option | Description |
|--------|-------------|
| `--tsv-log [FILE]` | Output memory statistics to TSV file |
| `-V [LEVEL]` | Verbose logging level: 0=ERROR, 1=WARN, 2=INFO, 3=DEBUG |
| `-R` | Enable resource monitoring (memory/CPU usage) |
| `-T` | Enable timing information |
| `-l [FILE]` | Human-readable log file path |

## TSV Log Format

### Header

```
# tronko-assign memory log v1.0
wall_time	phase	rss_mb	vm_mb	peak_rss_mb	cpu_user	cpu_sys	extra_info
```

### Columns

| Column | Description |
|--------|-------------|
| `wall_time` | Seconds since program start |
| `phase` | Current processing phase |
| `rss_mb` | Resident Set Size in MB (physical memory) |
| `vm_mb` | Virtual Memory size in MB |
| `peak_rss_mb` | Peak RSS observed so far |
| `cpu_user` | User CPU time in seconds |
| `cpu_sys` | System CPU time in seconds |
| `extra_info` | Phase-specific details (e.g., trees=5, batch=3) |

### Phases Logged

| Phase | Description | Extra Info |
|-------|-------------|------------|
| `STARTUP` | Program initialization | - |
| `TREE_ALLOCATED` | After each tree's memory allocated (DEBUG only) | tree=N, nodes=N, bases=N |
| `TREE_LOADED` | After each tree's posteriors loaded (DEBUG only) | tree=N |
| `REFERENCE_LOADED` | All reference data loaded | trees=N |
| `BWA_INDEX` | BWA FM-index constructed | - |
| `THREADS_ALLOCATED` | Thread structures initialized | threads=N |
| `BATCH_START` | Beginning of batch processing | batch=N |
| `BATCH_LOADED` | Reads loaded into memory | batch=N, reads=N |
| `BATCH_COMPLETE` | Batch processing finished | batch=N |
| `FINAL` | Program end | - |

## Example Output

### INFO Level (-V2)

```
# tronko-assign memory log v1.0
wall_time	phase	rss_mb	vm_mb	peak_rss_mb	cpu_user	cpu_sys	extra_info
0.002	STARTUP	2.4	4.3	2.4	0.001	0.000
0.298	REFERENCE_LOADED	52.8	54.6	52.8	0.273	0.023	trees=1
0.359	BWA_INDEX	53.8	55.3	57.2	0.326	0.025
0.360	THREADS_ALLOCATED	55.6	57.0	57.2	0.326	0.025	threads=1
0.360	BATCH_START	55.6	57.0	57.2	0.326	0.025	batch=1
0.360	BATCH_LOADED	55.6	57.0	57.2	0.326	0.025	batch=1,reads=100
0.450	BATCH_COMPLETE	55.6	57.0	57.2	0.400	0.030	batch=1
0.500	FINAL	55.5	56.8	57.2	0.450	0.035
```

### DEBUG Level (-V3)

Includes per-tree entries:

```
0.017	TREE_ALLOCATED	52.8	54.6	52.8	0.002	0.014	tree=0,nodes=2931,bases=316
0.302	TREE_LOADED	52.7	54.6	52.7	0.289	0.012	tree=0
```

## Analysis Scripts

Two analysis scripts are provided in `tronko-assign/scripts/`:

### compare_memlogs.sh

Compare memory usage across multiple runs:

```bash
./scripts/compare_memlogs.sh run1.tsv run2.tsv run3.tsv
```

Output (CSV format):
```
Run,Peak_RSS_MB,Time_to_Peak_s,Final_RSS_MB,Total_Time_s,CPU_User_s,CPU_Sys_s
run1,9876.5,0.456,234.5,45.678,42.123,3.456
run2,9901.2,0.512,245.6,47.890,44.567,3.678
run3,9856.7,0.423,228.9,43.567,40.789,3.234
```

### plot_memlog.py

Generate visual comparison plots:

```bash
python scripts/plot_memlog.py run1.tsv run2.tsv -o comparison.png
```

Requirements:
- pandas
- matplotlib

Install with: `pip install pandas matplotlib`

## Use Cases

### 1. Performance Baseline

Establish a memory baseline for your reference database:

```bash
./tronko-assign -V2 -R --tsv-log baseline.tsv \
    -r -f my_reference.txt -a my_ref.fasta \
    -s -g sample.fasta -o results.txt
```

### 2. Compare Different Configurations

Test impact of thread count on memory:

```bash
# 4 threads
./tronko-assign -V2 -R -C 4 --tsv-log 4threads.tsv ...

# 8 threads
./tronko-assign -V2 -R -C 8 --tsv-log 8threads.tsv ...

# Compare
./scripts/compare_memlogs.sh 4threads.tsv 8threads.tsv
```

### 3. Debug Large Reference Databases

Identify which trees consume most memory:

```bash
./tronko-assign -V3 -R --tsv-log trees.tsv -r -f large_reference.txt ...

# Find trees with highest memory delta
awk -F'\t' '/TREE_ALLOCATED/ {print $8, $3}' trees.tsv | sort -t= -k2 -n
```

### 4. Monitor Batch Processing

Track memory stability across batches:

```bash
./tronko-assign -V2 -R --tsv-log batches.tsv -L 10000 ...

# Extract batch memory
grep "BATCH_COMPLETE" batches.tsv | cut -f1,3,8
```

## Interpreting Results

### Expected Memory Profile

```
Phase                   Typical % of Peak    Notes
────────────────────────────────────────────────────
STARTUP                 <1%                  Minimal baseline
REFERENCE_LOADED        80-95%               Largest allocation
BWA_INDEX              +2-5%                 FM-index structure
THREADS_ALLOCATED      +1-5%                 Per-thread arrays
BATCH_*                +1-2%                 Temporary read data
FINAL                  ~1-5%                 After cleanup
```

### Key Metrics

1. **Peak RSS** - Maximum physical memory used during run
2. **Time to Peak** - When peak memory is reached (usually after REFERENCE_LOADED)
3. **Final RSS** - Memory after cleanup (should be much lower than peak)
4. **Memory Delta** - Difference between phases indicates allocation sizes

### Red Flags

- **Peak at BATCH_COMPLETE instead of REFERENCE_LOADED**: Memory leak in processing
- **Final RSS >> STARTUP**: Incomplete cleanup
- **TREE_LOADED memory inconsistent**: Possible posterior loading issues
- **Memory growth across batches**: Leak in batch processing

## Combining with Human-Readable Logs

For comprehensive debugging, use both TSV and human-readable logs:

```bash
./tronko-assign -V3 -R -T \
    -l human.log \
    --tsv-log memory.tsv \
    -r -f reference.txt -a ref.fasta \
    -s -g query.fasta -o results.txt
```

This gives you:
- `memory.tsv` - Machine-parseable memory data
- `human.log` - Detailed human-readable log with milestones and context
- `stderr` - Real-time progress output

## Technical Details

### How Memory is Measured

The logging uses Linux's `/proc/self/status` for memory metrics and `getrusage()` for CPU time:

- **VmRSS**: Resident Set Size (physical RAM)
- **VmSize**: Total virtual memory
- **VmHWM**: High-water mark (peak RSS)
- **ru_utime/ru_stime**: User/system CPU time

### Performance Impact

- TSV logging adds minimal overhead (<1% runtime)
- Per-tree logging at DEBUG level has slightly higher overhead
- Logging is buffered and flushed per entry for reliability
- No impact when `--tsv-log` is not specified

### Portability

The memory measurement is Linux-specific (uses `/proc/self/status`). On other systems:
- macOS: Would need `task_info()` or `getrusage()` fallback
- Windows: Would need `GetProcessMemoryInfo()`

Currently only Linux is supported.

## Code References

- TSV logging macros: `tronko-assign/tsv_memlog.h`
- Resource monitoring: `tronko-assign/resource_monitor.c`
- Main instrumentation: `tronko-assign/tronko-assign.c`
- Per-tree logging: `tronko-assign/readreference.c`
- Analysis scripts: `tronko-assign/scripts/`

## Related Documentation

- Research: `thoughts/shared/research/2025-12-29-memory-performance-logging.md`
- Implementation Plan: `thoughts/shared/plans/2025-12-29-memory-performance-logging.md`
