---
date: 2026-01-01T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: eDNA-Explorer/tronko
topic: "16S_Bacteria Segfault Root Cause Analysis"
tags: [research, codebase, segfault, memory, bounds-checking]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: 16S_Bacteria Segfault Root Cause Analysis

**Date**: 2026-01-01T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: eDNA-Explorer/tronko

## Research Question

What causes the segmentation fault when processing the 16S_Bacteria dataset (17,868 trees, ~1.27M nodes), and what additional logging would help diagnose such issues?

## Summary

**Root Cause Identified**: The segfault is caused by a **buffer overflow** in the `leaf_iter` counter. Multiple arrays are allocated with size `MAX_NUM_BWA_MATCHES` (10), but `leaf_iter` can increment beyond 10 when processing reads that match many unique trees. With 17,868 trees in the 16S_Bacteria dataset, this boundary is much more likely to be exceeded than with the single_tree dataset.

**Secondary Issue**: Numerous memory allocations lack NULL checks, which could cause crashes if allocations fail under the ~13 GB memory pressure of this dataset.

## Detailed Findings

### Critical Bug: `leaf_iter` Buffer Overflow

#### The Problem

The variable `leaf_iter` tracks unique tree matches found for a read. It increments for each unique tree found across:
- Discordant forward matches
- Discordant reverse matches
- Concordant forward matches
- Concordant reverse matches

**Maximum potential value**: 4 × 10 = 40 unique trees

However, several arrays are sized only for `MAX_NUM_BWA_MATCHES` (10):

| Array | Allocation Size | Access Pattern |
|-------|-----------------|----------------|
| `nodeScores[i][j][k]` | `i < 10` | `nodeScores[leaf_iter][...]` |
| `starts_forward` | 10 | `starts_forward[leaf_iter]` |
| `starts_reverse` | 10 | `starts_reverse[leaf_iter]` |
| `cigars_forward` | 10 | `cigars_forward[leaf_iter]` |
| `cigars_reverse` | 10 | `cigars_reverse[leaf_iter]` |

#### Code References

**Allocation** (`allocateMemoryForResults.c:12-16`):
```c
results->nodeScores = (type_of_PP ***)malloc(MAX_NUM_BWA_MATCHES*(sizeof(type_of_PP **)));
for (i=0; i<MAX_NUM_BWA_MATCHES; i++){
    results->nodeScores[i] = (type_of_PP **)malloc(numberOfTrees*(sizeof(type_of_PP *)));
```

**Unbounded increment** (`tronko-assign.c:337, 370, 406, 443`):
```c
trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
leaf_iter++;  // NO BOUNDS CHECK!
```

**Out-of-bounds access** (`tronko-assign.c:315-316`):
```c
results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
```

#### Why This Affects Large Datasets

With 17,868 trees vs 1 tree:
- **single_tree**: Only 1 tree exists, so `leaf_iter` can never exceed 1
- **16S_Bacteria**: With 17,868 trees and diverse query sequences, reads frequently match more than 10 unique trees

### Secondary Issue: Missing NULL Checks

Almost all allocations lack NULL checks, which could cause crashes if memory allocation fails under ~13 GB pressure:

| File:Line | Allocation |
|-----------|------------|
| `readreference.c:559-561` | `numbaseArr`, `rootArr`, `numspecArr` |
| `readreference.c:592` | `treeArr` |
| `allocatetreememory.c:14` | Per-tree node array |
| `allocatetreememory.c:21` | Per-node `posteriornc` |
| `allocateMemoryForResults.c:12-16` | `nodeScores` (3D) |
| `allocateMemoryForResults.c:22-24` | `voteRoot` (2D) |
| `allocateMemoryForResults.c:75-78` | `minNodes`, `LCAnames` |

Only the taxonomy arrays use `calloc_check()` which properly exits on failure.

### Current Crash Debugging Infrastructure

The crash debugging system (`crash_debug.c`) is comprehensive but missing critical diagnostic data:

**Currently Captured**:
- Signal number, code, fault address
- Stack trace with symbols
- Register dump (RAX, RBX, RCX, RDX, RSI, RDI, RBP, RSP, RIP)
- Memory state (RSS, VM size)
- Processing stage, current file, current tree

**Missing for This Bug**:
- `leaf_iter` value at crash time
- `MAX_NUM_BWA_MATCHES` constant value
- Array allocation sizes vs access indices
- Per-read match count breakdown (concordant/discordant × forward/reverse)
- Bounds violation detection

## Recommended Additional Logging

### 1. Add Array Bounds Context to Crash Reports

Add to `crash_debug.h` and track in `g_app_context`:

```c
// Array bounds tracking
int current_leaf_iter;
int max_bwa_matches;  // Constant for reference
int concordant_match_count;
int discordant_match_count;
```

### 2. Add Bounds Checking with Diagnostic Output

Before incrementing `leaf_iter` (`tronko-assign.c:337`):

```c
if (leaf_iter >= MAX_NUM_BWA_MATCHES) {
    crash_set_context("BOUNDS VIOLATION: leaf_iter overflow");
    LOG_ERROR("leaf_iter=%d exceeds MAX_NUM_BWA_MATCHES=%d for read %s",
              leaf_iter, MAX_NUM_BWA_MATCHES, current_read_name);
    // Either cap or abort with diagnostic
}
leaf_iter++;
```

### 3. Add Per-Read Match Summary Logging

In `tronko-assign.c` after BWA matching:

```c
LOG_DEBUG("Read %s: %d concordant + %d discordant matches across %d unique trees",
          read_name, concordant_count, discordant_count, leaf_iter);
```

### 4. Track Allocation Sizes in Crash Context

After large allocations:

```c
crash_set_allocation_info("nodeScores", MAX_NUM_BWA_MATCHES, numberOfTrees, total_bytes);
```

### 5. Add Memory Allocation Wrapper with Logging

Create `checked_malloc()` that logs allocation attempts and failures:

```c
void* checked_malloc(size_t size, const char* name, const char* file, int line) {
    void* ptr = malloc(size);
    if (!ptr) {
        LOG_ERROR("malloc failed: %s, size=%zu at %s:%d", name, size, file, line);
        crash_flag_allocation_failure(name, size, file, line);
    }
    return ptr;
}
```

### 6. TSV Logging Enhancements

Add to the TSV memory log format (`tsv_memlog.h`):

```
phase   leaf_iter   max_matches   unique_trees   match_type
```

This would help correlate memory spikes with match patterns.

## Proposed Fix for Root Cause

### Option A: Cap `leaf_iter` at `MAX_NUM_BWA_MATCHES`

```c
// In tronko-assign.c before each leaf_iter++
if (leaf_iter < MAX_NUM_BWA_MATCHES) {
    trees_search[index1] = results->leaf_coordinates[leaf_iter][0];
    leaf_iter++;
} else {
    LOG_WARN("Exceeded MAX_NUM_BWA_MATCHES, dropping additional matches");
}
```

**Pros**: Minimal code change, safe
**Cons**: May lose valuable matches

### Option B: Increase `MAX_NUM_BWA_MATCHES`

Change from 10 to 40+ in `global.h:43`:

```c
#define MAX_NUM_BWA_MATCHES 50
```

**Pros**: Captures more matches
**Cons**: Increases memory usage per thread (~5× for nodeScores)

### Option C: Dynamic Allocation Based on Actual Matches

Allocate arrays after counting matches, sized to actual need.

**Pros**: Optimal memory usage
**Cons**: Significant code refactoring required

## Architecture Insights

The array dimension mismatch reveals a design inconsistency:

- `leaf_coordinates` is sized for `numberOfTrees` (17,868) - anticipating many unique trees
- `nodeScores`, `starts_*`, `cigars_*` are sized for `MAX_NUM_BWA_MATCHES` (10) - assuming few unique matches

These two assumptions are incompatible. The code should use consistent sizing based on the actual use case.

## Open Questions

1. **Is capping at 10 matches acceptable?** Need to understand biological implications of dropping matches beyond 10.

2. **Memory budget**: If increasing `MAX_NUM_BWA_MATCHES`, what's the memory impact at scale?

3. **Deduplication logic**: Should matches be deduplicated by tree+node, not just tree?

## Related Files

- Bug report: `thoughts/shared/bugs/2026-01-01-16S-bacteria-segfault.md`
- Benchmark results: `thoughts/shared/research/2026-01-01-tier1-benchmark-results.md`
- Crash reports: `/tmp/tronko_assign_crash_*.crash`

## Code References

- `tronko-assign/tronko-assign.c:298-443` - leaf_iter increment logic
- `tronko-assign/allocateMemoryForResults.c:12-48` - Array allocations with MAX_NUM_BWA_MATCHES
- `tronko-assign/global.h:43` - `MAX_NUM_BWA_MATCHES` definition
- `tronko-assign/crash_debug.c:29-46` - Current crash context structure
- `tronko-assign/placement.c:889-937` - nodeScores access pattern
