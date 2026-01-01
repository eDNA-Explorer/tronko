---
date: 2026-01-01T15:30:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: eDNA-Explorer/tronko
topic: "BWA Implementation Changes and Crash Analysis"
tags: [research, codebase, bwa, crash, memory, bounds-checking]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: BWA Implementation Changes and Crash Analysis

**Date**: 2026-01-01T15:30:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: eDNA-Explorer/tronko

## Research Question

Does the codebase have changes made to the BWAMap implementation that could be causing crashes?

## Summary

**Yes, significant changes have been made to BWA-related code on the `experimental` branch.** The changes include:

1. **Memory layout optimization**: Changed `posteriornc` from 2D (`double**`) to 1D (`type_of_PP*`) array with `PP_IDX()` macro
2. **Bounds checking additions**: Added `MAX_NUM_BWA_MATCHES` bounds checks to prevent buffer overflows (was causing segfaults)
3. **Bug fix**: Fixed assignment-instead-of-comparison typo (`=-1` → `==-1`)
4. **Crash debugging infrastructure**: New comprehensive crash debugging system with BWA-specific context tracking
5. **Early termination optimization**: New parameters for tree pruning during score assignment

**Key finding**: The bounds checking fix in `tronko-assign.c` addresses a known buffer overflow bug where `leaf_iter` could exceed `MAX_NUM_BWA_MATCHES` (10), causing crashes on large datasets. However, `MAX_NUM_BWA_MATCHES` was also reduced from 5000 (tronko-build) to 10 (tronko-assign), which may cause excessive match dropping.

## Detailed Findings

### 1. BWA Source File Changes

**File**: `tronko-assign/bwa_source_files/bwamem.c`

```diff
-#ifdef HAVE_PTHREAD
 #include <pthread.h>
-#endif
```

**Impact**: pthread is now always included unconditionally. This is a minor change unlikely to cause crashes.

**File**: `tronko-assign/bwa_source_files_include.h`

```diff
+// Forward declaration of main_mem from fastmap.c
+int main_mem(char* databaseFile, int number_of_seqs, int number_of_threads,
+             bwaMatches* bwa_results, int concordant, int numberOfTrees,
+             int startline, int paired, int start, int end,
+             int max_query_length, int max_readname_length, int max_acc_name);
```

**Impact**: Cleaner header organization, no runtime impact.

### 2. Memory Layout Changes (CRITICAL)

**File**: `tronko-assign/global.h:87`

```diff
-	double **posteriornc;
+	type_of_PP *posteriornc;  /* 1D array: access with PP_IDX(pos, nuc) */
```

**New macro**: `PP_IDX(pos, nuc) = ((pos) * 4 + (nuc))`

**Files affected**:
- `assignment.c:143-227` - All `posteriornc[i][j]` → `posteriornc[PP_IDX(i, j)]`
- `allocatetreememory.c:4-23` - Changed from nested `malloc` loops to single `calloc`
- `tronko-assign.c:63-75` - Updated `store_PPs_Arr()` function

**Impact**: This is a significant memory layout change that could cause crashes if:
- Old binary files are used with new code (format mismatch)
- Any code path was missed during the conversion
- Index calculations are incorrect

### 3. Bounds Checking Fix (CRITICAL)

**File**: `tronko-assign/tronko-assign.c:320-480`

Multiple instances of bounds checking added:

```c
if (leaf_iter < MAX_NUM_BWA_MATCHES) {
    results->leaf_coordinates[leaf_iter][0] = ...;
    // ... other operations ...
    leaf_iter++;
} else {
    dropped_matches_count++;
}
```

**Root cause this fixes**: The original code had no bounds check on `leaf_iter`, which could increment beyond `MAX_NUM_BWA_MATCHES` (10), causing buffer overflow when accessing:
- `results->leaf_coordinates[leaf_iter]`
- `results->starts_forward[leaf_iter]`
- `results->cigars_forward[leaf_iter]`
- `results->nodeScores[leaf_iter]`

### 4. Bug Fix: Assignment vs Comparison

**File**: `tronko-assign/tronko-assign.c:462`

```diff
-						if (trees_search[k]=-1){
+						if (trees_search[k]==-1){
```

**Impact**: The original code was assigning `-1` to `trees_search[k]` instead of comparing. This would:
1. Always evaluate to true (since `-1` is non-zero)
2. Corrupt the `trees_search` array
3. Cause incorrect tree matching behavior

### 5. Crash Debug Infrastructure (NEW)

**New files**:
- `tronko-assign/crash_debug.c` (724+ lines)
- `tronko-assign/crash_debug.h` (185+ lines)

**BWA-specific tracking functions**:
```c
crash_set_bwa_context(int leaf_iter, int max_matches, int concordant, int discordant);
crash_set_bwa_bounds_violation(int leaf_iter, int max_matches, int dropped);
crash_clear_bwa_context();
```

**Fields in crash report**:
- `bwa_leaf_iter` - Current leaf iterator value
- `bwa_max_matches` - MAX_NUM_BWA_MATCHES constant
- `bwa_concordant_count` - Concordant matches found
- `bwa_discordant_count` - Discordant matches found
- `bwa_dropped_matches` - Matches dropped due to bounds

### 6. Early Termination Optimization (NEW)

**File**: `tronko-assign/assignment.c:24-76`

The `assignScores_Arr_paired()` function signature changed to include optimization parameters:

```diff
-void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
-    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
-    FILE* site_scores_file, char* readname);
+void assignScores_Arr_paired(int rootNum, int node, char *locQuery, int *positions,
+    type_of_PP ***scores, int alength, int search_number, int print_all_nodes,
+    FILE* site_scores_file, char* readname,
+    int early_termination, type_of_PP *best_score, int *strikes,
+    type_of_PP strike_box, int max_strikes,
+    int enable_pruning, type_of_PP pruning_threshold);
```

**Purpose**: Allows early termination of tree traversal when scores become too poor, and subtree pruning to skip unpromising branches.

### 7. Constants Mismatch (POTENTIAL ISSUE)

```c
// tronko-build/global.h:19
#define MAX_NUM_BWA_MATCHES 5000

// tronko-assign/global.h:43
#define MAX_NUM_BWA_MATCHES 10
```

**Impact**: The build tool uses 5000, but assign uses 10. This could cause:
- Excessive match dropping for datasets with many trees
- Different behavior than originally expected
- Warning spam in logs

## Potential Crash Causes

Based on the analysis, the most likely crash causes are:

### 1. Memory Layout Mismatch (if using old reference files)

If a reference file was built before the `posteriornc` layout change, loading it would cause memory corruption:
- Old format: `posteriornc[position][nucleotide]` as 2D array of pointers
- New format: `posteriornc[position * 4 + nucleotide]` as 1D array

**Symptom**: Segfault during reference loading or LCA calculation

### 2. Incomplete Bounds Check Coverage

While many bounds checks were added, verify all access patterns are covered:
- `tronko-assign.c:323-435` - Main match processing loop (mostly covered)
- `placement.c:449-853` - Score assignment calls (checked)

### 3. Early Termination State Corruption

The new early termination code introduces state variables:
```c
type_of_PP best_score = -9999999999999999;
int strikes = 0;
```

If these get corrupted or the threshold calculations overflow, it could cause premature termination or infinite loops.

### 4. Thread-Local Counter Issues

```c
static __thread int dropped_matches_count = 0;
```

Thread-local storage is used for tracking dropped matches. If `__thread` isn't properly supported on the target platform, this could cause data races.

## Recommendations

### Immediate Actions

1. **Verify reference file format**: Ensure reference files used were built with the same `posteriornc` layout
2. **Check MAX_NUM_BWA_MATCHES**: If processing large datasets, the value of 10 is likely too small
3. **Review crash logs**: Check `/tmp/tronko_assign_crash_*.crash` for specific crash context

### Testing

```bash
# Run with verbose logging
tronko-assign -r -f reference.txt -a query.fasta -s -o output.txt -v

# Check for bounds violations in output
grep "overflow" output.log
grep "dropped" output.log
```

### If Crashes Persist

1. Build with debug symbols: `make debug`
2. Run under valgrind: `valgrind --track-origins=yes ./tronko-assign ...`
3. Check the crash report generated in `/tmp/`

## Code References

- `tronko-assign/tronko-assign.c:320-480` - Bounds checking logic
- `tronko-assign/global.h:29-30` - `PP_IDX` macro definition
- `tronko-assign/global.h:87` - `posteriornc` type change
- `tronko-assign/assignment.c:24-76` - Modified score assignment function
- `tronko-assign/crash_debug.c:698-724` - BWA context tracking
- `tronko-assign/allocatetreememory.c:4-23` - Memory allocation changes

## Related Research

- `thoughts/shared/research/2026-01-01-16S-bacteria-segfault-analysis.md` - Original buffer overflow analysis
- `thoughts/shared/plans/2026-01-01-bwa-bounds-overflow-fix.md` - Fix implementation plan

## Open Questions

1. Should `MAX_NUM_BWA_MATCHES` be increased to match tronko-build's value (5000)?
2. Are there any code paths that still use the old 2D `posteriornc` access pattern?
3. Is the reference file format documented to ensure compatibility?
