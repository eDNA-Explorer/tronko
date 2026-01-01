# Bug Report: Segmentation Fault with 16S_Bacteria Dataset

**Date**: 2026-01-01
**Severity**: Critical
**Status**: Open

## Summary

tronko-assign crashes with SIGSEGV when processing reads against the 16S_Bacteria reference database. The crash occurs during batch processing after successfully loading the reference database.

## Environment

- Platform: Linux 6.17.9-arch1-1
- Kernel: 6.17.9-arch1-1
- tronko-assign: experimental branch (commit 77ade9e)
- GCC: System default
- Dataset: 16S_Bacteria (17,868 trees, ~1.27M total nodes)

## Complete Reproduction Steps

### 1. Prepare Test Data

The test reads were extracted from compressed production data files:

```bash
# Source files (in /home/jimjeffers/Work/tronko/)
# clk06lckc0001jn0f9hw7rb86-16S_Bacteria-paired_F.fasta.zst (251,950,316 bytes)
# clk06lckc0001jn0f9hw7rb86-16S_Bacteria-paired_R.fasta.zst (230,948,368 bytes)

# Decompress full files
cd /home/jimjeffers/Work/tronko
zstd -d -k clk06lckc0001jn0f9hw7rb86-16S_Bacteria-paired_F.fasta.zst -o /tmp/16S_reads_F.fasta
zstd -d -k clk06lckc0001jn0f9hw7rb86-16S_Bacteria-paired_R.fasta.zst -o /tmp/16S_reads_R.fasta

# Full files contain 9,496,090 reads each
# Create 1000-read subset (4000 lines = 1000 reads, 4 lines per FASTA entry)
head -4000 /tmp/16S_reads_F.fasta > /tmp/16S_subset_1k_F.fasta
head -4000 /tmp/16S_reads_R.fasta > /tmp/16S_subset_1k_R.fasta
```

### 2. Reference Database Files Used

Located in `/home/jimjeffers/Work/tronko/example_datasets/16S_Bacteria/`:

| File | Size | Description |
|------|------|-------------|
| `reference_tree.trkb.gz` | 369,096,550 bytes | Gzipped binary format |
| `reference_tree.txt.gz` | 793,999,337 bytes | Gzipped text format |
| `16S_Bacteria.fasta` | 180,496,222 bytes | Reference sequences |
| `16S_Bacteria.fasta.bwt` | 172,148,652 bytes | BWA index (pre-built) |
| `16S_Bacteria.fasta.sa` | 86,074,328 bytes | BWA index (pre-built) |
| `16S_Bacteria.fasta.pac` | 43,037,140 bytes | BWA index (pre-built) |
| `16S_Bacteria.fasta.ann` | 22,695,885 bytes | BWA index (pre-built) |
| `16S_Bacteria.fasta.amb` | 608,745 bytes | BWA index (pre-built) |

### 3. Build tronko-assign

```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# Clean all build artifacts
make clean
# Output: rm -f tronko-assign

# Build baseline (no optimization flags)
make
# Output: gcc -O3 -fno-stack-protector -o tronko-assign [sources...] -lm -pthread -lz -lrt -std=gnu99

# Copy to baseline binary
cp tronko-assign tronko-assign-baseline
```

### 4. Run Failing Command

```bash
cd /home/jimjeffers/Work/tronko/tronko-assign

# Attempt 1: Binary format with skip-bwa-build
./tronko-assign-baseline -r \
    -f ../example_datasets/16S_Bacteria/reference_tree.trkb.gz \
    -a ../example_datasets/16S_Bacteria/16S_Bacteria.fasta \
    -s -g /tmp/16S_subset_1k_F.fasta \
    -o /tmp/baseline_16s.txt -w -6

# Result: Segmentation fault (core dumped)

# Attempt 2: Text format with skip-bwa-build
./tronko-assign-baseline -r \
    -f ../example_datasets/16S_Bacteria/reference_tree.txt.gz \
    -a ../example_datasets/16S_Bacteria/16S_Bacteria.fasta \
    -s -g /tmp/16S_subset_1k_F.fasta \
    -o /tmp/baseline_16s.txt -w -6

# Result: Segmentation fault (core dumped)

# Attempt 3: With verbose logging
./tronko-assign-baseline -r \
    -f ../example_datasets/16S_Bacteria/reference_tree.trkb.gz \
    -a ../example_datasets/16S_Bacteria/16S_Bacteria.fasta \
    -s -g /tmp/16S_subset_1k_F.fasta \
    -o /tmp/baseline_16s.txt -w -6 -V 3

# Result: Segmentation fault with detailed logging (see below)
```

## Verbose Output Before Crash

```
[2026-01-01 13:35:49.972] [INFO] Skipping BWA index build
[2026-01-01 13:35:49.972] [INFO] MILESTONE: READ_SPECS_DETECTED - Read specs detected: max_name=26, max_query=60 [11.441s total, 4.804s since last]
[2026-01-01 13:35:50.094] [INFO] MILESTONE: MEMORY_ALLOCATED - Memory allocated: threads=1, lines_per_batch=50000, total_nodes=1265572 [11.563s total, 0.122s since last]
[2026-01-01 13:35:50.094] [INFO] MILESTONE: THREADS_INITIALIZED - Thread structures initialized for 1 cores [11.563s total, 0.000s since last]
[2026-01-01 13:35:50.094] [INFO] MILESTONE: BATCH_START - Starting batch 1 [11.563s total, 0.000s since last]
[2026-01-01 13:35:50.094] [INFO] MILESTONE: BATCH_LOADED - Batch 1 loaded: 3291 reads [11.563s total, 0.000s since last]
[2026-01-01 13:35:51.910] [ERROR] === CRASH DETECTED ===
[2026-01-01 13:35:51.910] [ERROR] Process: 1154924, Signal: SIGSEGV
[2026-01-01 13:35:51.910] [ERROR] Message: CRASH: SIGSEGV (Segmentation fault) at address (nil) in process 1154924
```

## Crash Report

Full crash report from `/tmp/tronko_assign_crash_1165173_1767303641.crash`:

```
=== TRONKO CRASH REPORT ===
Crash Time: Thu Jan  1 13:40:41 2026
Process ID: 1165173
Program: /home/jimjeffers/Work/tronko/tronko-assign/tronko-assign-baseline
Working Directory: /home/jimjeffers/Work/tronko/tronko-assign
Command Line: unknown

CRASH: SIGSEGV (Segmentation fault) at address (nil) in process 1165173

Application Context:
  Processing Stage: Loading reference database
  Current File: ../example_datasets/16S_Bacteria/reference_tree.txt.gz (line 4000)

Stack Trace:
#0  ./tronko-assign-baseline(+0x680db) [0x55a34dc5c0db]
#1  /usr/lib/libc.so.6(+0x3e4d0) [0x7f5a7c03e4d0]
#2  ./tronko-assign-baseline(+0x30382) [0x55a34dc24382]
#3  ./tronko-assign-baseline(+0x36ac2) [0x55a34dc2aac2]
#4  /usr/lib/libc.so.6(+0x9698b) [0x7f5a7c09698b]
#5  /usr/lib/libc.so.6(+0x11a9cc) [0x7f5a7c11a9cc]

Register state:
RAX: 0x00007f5a7a15d480  RBX: 0x00007f5a7b9c2010
RCX: 0x0000000000000000  RDX: 0x0000000000000000
RSI: 0x00007f5a620ae970  RDI: 0x000055a6652ca330
RBP: 0x00007f5a620aee70  RSP: 0x00007f5a620ae950
RIP: 0x000055a34dc24382

Memory state:
RSS: 13410768 KB (~13 GB)
VM Size: 13530092 KB (~13 GB)
CPU Usage: 0.00%

=== END CRASH REPORT ===
```

## Working Comparison (single_tree dataset)

The same binary works correctly with the smaller single_tree dataset:

```bash
# This works
./tronko-assign-baseline -r \
    -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/baseline_single.txt -w -6

# Output: 165 lines (164 reads + header), completes in ~0.5 seconds
```

## Dataset Comparison

| Dataset | Trees | Total Nodes | Memory Usage | Result |
|---------|-------|-------------|--------------|--------|
| single_tree | 1 | 2,931 | ~40 MB | Works |
| 16S_Bacteria | 17,868 | 1,265,572 | ~13 GB | **CRASH** |

## Analysis

### Key Observations

1. **Crash timing**: Occurs after successful reference loading, during first read batch processing
2. **NULL pointer**: Address `(nil)` indicates dereferencing a NULL pointer
3. **Memory**: 13 GB allocated successfully before crash
4. **Batch loaded**: 3,291 reads loaded (note: subset file has 1,000 reads but logging shows 3,291)
5. **Both formats crash**: Issue occurs with both `.trkb.gz` and `.txt.gz` formats

### Possible Causes

1. **Array bounds overflow**: With 17,868 trees, index calculations may overflow int limits
2. **Memory allocation failure**: Some allocation may fail silently, returning NULL
3. **Uninitialized pointer**: Data structure not properly initialized for large tree counts
4. **Tree/node indexing bug**: `leaf_coordinates` or `nodeScores` arrays may be incorrectly sized

### Relevant Code Locations

Based on stack trace offset `+0x30382`, likely in:
- `placement.c` - `runAssignmentsSingle()` or read processing loop
- `assignment.c` - `assignScores_Arr_paired()`

## Files Summary

### Input Files
- `/tmp/16S_subset_1k_F.fasta` - 1,000 forward reads (created from zst archive)
- `example_datasets/16S_Bacteria/reference_tree.trkb.gz` - Binary reference database
- `example_datasets/16S_Bacteria/16S_Bacteria.fasta` - Reference sequences

### Build Artifacts
- `tronko-assign` - Compiled binary (removed by `make clean`, rebuilt by `make`)
- `tronko-assign-baseline` - Copy of compiled binary for testing

### Output Files
- `/tmp/baseline_16s.txt` - Empty (crash before output)
- `/tmp/tronko_assign_crash_*.crash` - Crash reports

## Workaround

None currently available. The 16S_Bacteria dataset cannot be processed.
