# BWA Bounds Overflow Fix - Implementation Plan

## Overview

Fix the buffer overflow bug in `tronko-assign` where `leaf_iter` can exceed `MAX_NUM_BWA_MATCHES` (10) causing segfaults when processing large datasets like 16S_Bacteria (17,868 trees). Add diagnostic logging to help troubleshoot similar issues in the future.

## Current State Analysis

### Root Cause
The variable `leaf_iter` tracks unique tree matches found for a read. It increments for each unique tree found across all match types (concordant forward, concordant reverse, discordant forward, discordant reverse).

**Maximum potential value**: 4 x 10 = 40 unique trees

However, these arrays are sized for only `MAX_NUM_BWA_MATCHES` (10):
- `results->nodeScores[leaf_iter]` - `allocateMemoryForResults.c:12`
- `results->starts_forward[leaf_iter]` - `allocateMemoryForResults.c:30`
- `results->starts_reverse[leaf_iter]` - `allocateMemoryForResults.c:32`
- `results->cigars_forward[leaf_iter]` - `allocateMemoryForResults.c:34`
- `results->cigars_reverse[leaf_iter]` - `allocateMemoryForResults.c:36`

### Key Discoveries
- `tronko-assign.c:298` - `leaf_iter` initialized to 0
- `tronko-assign.c:337,370,406,443` - `leaf_iter++` without bounds checking
- `tronko-assign.c:315-316,349-350,385-386,421-426` - Array accesses using `leaf_iter`
- Single_tree dataset: `leaf_iter` never exceeds 1 (only 1 tree)
- 16S_Bacteria: With 17,868 trees, reads frequently match >10 unique trees

### Missing Crash Context
The current crash debugging system captures:
- Signal info, stack trace, register dump
- Current file, current tree, processing stage

But NOT:
- `leaf_iter` value at crash time
- BWA match counts breakdown
- Array bounds violation detection

## Desired End State

After implementation:
1. Bounds checking prevents `leaf_iter` from exceeding `MAX_NUM_BWA_MATCHES`
2. Warning logged when matches are dropped due to exceeding limit
3. Crash reports include BWA-specific context (leaf_iter, match counts)
4. Per-read match summary available in verbose/debug mode

### Verification:
- Run against 16S_Bacteria dataset without segfault
- Warning messages appear for reads with >10 unique tree matches
- Crash reports (if any occur) include BWA bounds context

## What We're NOT Doing

- Increasing `MAX_NUM_BWA_MATCHES` (would increase memory usage ~5x)
- Dynamic allocation based on actual matches (significant refactoring)
- Changing the fundamental algorithm

## Implementation Approach

Use Option A from the research: Cap `leaf_iter` at `MAX_NUM_BWA_MATCHES` with diagnostic logging. This is the safest minimal change that prevents crashes while providing visibility into dropped matches.

---

## Phase 1: Add BWA Context to Crash Debug System

### Overview
Extend `crash_debug.h/c` to track BWA-specific context that will appear in crash reports.

### Changes Required:

#### 1. crash_debug.h - Add BWA context fields
**File**: `tronko-assign/crash_debug.h`
**Changes**: Add new fields to `crash_info_t` and new function declarations

After line 78 (after `corruption_type`), add:
```c
    // BWA bounds tracking (added for 16S_Bacteria segfault analysis)
    int bwa_leaf_iter;                // Current leaf_iter value
    int bwa_max_matches;              // MAX_NUM_BWA_MATCHES constant
    int bwa_concordant_count;         // Concordant matches for current read
    int bwa_discordant_count;         // Discordant matches for current read
    int bwa_unique_trees;             // Unique trees matched
    int bwa_dropped_matches;          // Matches dropped due to bounds
```

After line 133 (after `crash_clear_corruption_flags`), add:
```c
// BWA bounds tracking functions
void crash_set_bwa_context(int leaf_iter, int concordant_count, int discordant_count);
void crash_set_bwa_bounds_violation(int leaf_iter, int max_matches, int dropped);
void crash_clear_bwa_context(void);
```

#### 2. crash_debug.c - Add BWA context implementation
**File**: `tronko-assign/crash_debug.c`
**Changes**: Add BWA tracking to g_app_context and implement functions

After line 45 (after `pthread_mutex_t context_mutex`), add:
```c
    // BWA bounds tracking
    int bwa_leaf_iter;
    int bwa_max_matches;
    int bwa_concordant_count;
    int bwa_discordant_count;
    int bwa_unique_trees;
    int bwa_dropped_matches;
```

After the signal handler captures application context (~line 205), add:
```c
    // Capture BWA context
    crash_info.bwa_leaf_iter = g_app_context.bwa_leaf_iter;
    crash_info.bwa_max_matches = g_app_context.bwa_max_matches;
    crash_info.bwa_concordant_count = g_app_context.bwa_concordant_count;
    crash_info.bwa_discordant_count = g_app_context.bwa_discordant_count;
    crash_info.bwa_unique_trees = g_app_context.bwa_unique_trees;
    crash_info.bwa_dropped_matches = g_app_context.bwa_dropped_matches;
```

In `crash_generate_report` after application context section (~line 458), add:
```c
    // Write BWA bounds information if relevant
    if (crash_info->bwa_leaf_iter > 0 || crash_info->bwa_dropped_matches > 0) {
        fprintf(report_file, "\nBWA Bounds Context:\n");
        fprintf(report_file, "  leaf_iter: %d (max: %d)\n",
                crash_info->bwa_leaf_iter, crash_info->bwa_max_matches);
        fprintf(report_file, "  Concordant matches: %d\n", crash_info->bwa_concordant_count);
        fprintf(report_file, "  Discordant matches: %d\n", crash_info->bwa_discordant_count);
        fprintf(report_file, "  Unique trees: %d\n", crash_info->bwa_unique_trees);
        if (crash_info->bwa_dropped_matches > 0) {
            fprintf(report_file, "  *** MATCHES DROPPED: %d (bounds exceeded) ***\n",
                    crash_info->bwa_dropped_matches);
        }
    }
```

At end of file, add new functions:
```c
// BWA bounds tracking functions
void crash_set_bwa_context(int leaf_iter, int concordant_count, int discordant_count) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = leaf_iter;
    g_app_context.bwa_concordant_count = concordant_count;
    g_app_context.bwa_discordant_count = discordant_count;
    g_app_context.bwa_unique_trees = leaf_iter;  // leaf_iter tracks unique trees
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_bwa_bounds_violation(int leaf_iter, int max_matches, int dropped) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = leaf_iter;
    g_app_context.bwa_max_matches = max_matches;
    g_app_context.bwa_dropped_matches = dropped;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_clear_bwa_context(void) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = 0;
    g_app_context.bwa_max_matches = 0;
    g_app_context.bwa_concordant_count = 0;
    g_app_context.bwa_discordant_count = 0;
    g_app_context.bwa_unique_trees = 0;
    g_app_context.bwa_dropped_matches = 0;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}
```

### Success Criteria:

#### Automated Verification:
- [x] Clean build: `cd tronko-assign && make clean && make`
- [x] No compiler warnings related to new code

#### Manual Verification:
- [x] Code review confirms all new fields initialized properly

---

## Phase 2: Add Bounds Checking to leaf_iter

### Overview
Add bounds checking before each `leaf_iter++` to prevent buffer overflow. Log warnings when matches are dropped.

### Changes Required:

#### 1. tronko-assign.c - Add bounds checking macro and apply to all increments
**File**: `tronko-assign/tronko-assign.c`

Add after the includes (~line 28), a helper macro:
```c
// Bounds-checked leaf_iter increment with logging
// Returns 1 if increment was allowed, 0 if bounds would be exceeded
static int dropped_matches_count = 0;  // Track dropped matches per read

#define SAFE_LEAF_ITER_INCREMENT(leaf_iter, read_name) \
    do { \
        if ((leaf_iter) < MAX_NUM_BWA_MATCHES) { \
            (leaf_iter)++; \
        } else { \
            dropped_matches_count++; \
            if (dropped_matches_count == 1) { \
                LOG_WARN("BWA bounds: leaf_iter=%d >= MAX_NUM_BWA_MATCHES=%d for read %s, dropping additional matches", \
                         (leaf_iter), MAX_NUM_BWA_MATCHES, (read_name) ? (read_name) : "unknown"); \
            } \
        } \
    } while(0)
```

In `runAssignmentOnChunk_WithBWA`, after line 298 (`int leaf_iter=0;`), add:
```c
        dropped_matches_count = 0;  // Reset dropped counter for this read
```

Replace each `leaf_iter++;` with bounds-checked version:

**Line 337** (discordant concordant-mode):
```c
                if (found==0){
                    if (leaf_iter < MAX_NUM_BWA_MATCHES) {
                        trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
                        leaf_iter++;
                    } else {
                        dropped_matches_count++;
                    }
                }
```

**Line 370** (concordant matches):
```c
                    if (found==0){
                        if (leaf_iter < MAX_NUM_BWA_MATCHES) {
                            trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
                            leaf_iter++;
                        } else {
                            dropped_matches_count++;
                        }
                    }
```

**Line 406** (discordant non-concordant mode):
```c
                    if (found==0){
                        if (leaf_iter < MAX_NUM_BWA_MATCHES) {
                            trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
                            leaf_iter++;
                        } else {
                            dropped_matches_count++;
                        }
                    }
```

**Line 443** (concordant in non-concordant mode):
```c
                    if (found==0){
                        if (leaf_iter < MAX_NUM_BWA_MATCHES) {
                            trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
                            leaf_iter++;
                        } else {
                            dropped_matches_count++;
                        }
                    }
```

After each read processing loop (after `leaf_iter` usage complete, ~line 450), add:
```c
        // Update crash context and log if matches were dropped
        if (dropped_matches_count > 0) {
            crash_set_bwa_bounds_violation(leaf_iter, MAX_NUM_BWA_MATCHES, dropped_matches_count);
            LOG_WARN("Read %s: Dropped %d matches (leaf_iter capped at %d)",
                     paired ? pairedQueryMat->forward_name[lineNumber] : singleQueryMat->name[lineNumber],
                     dropped_matches_count, MAX_NUM_BWA_MATCHES);
        }
        crash_clear_bwa_context();
```

### Success Criteria:

#### Automated Verification:
- [x] Clean build: `cd tronko-assign && make clean && make`
- [ ] No segfault on 16S_Bacteria test: run with verbose mode and verify completion

#### Manual Verification:
- [ ] Warning messages appear for reads with >10 unique tree matches
- [ ] Output results match expected format

---

## Phase 3: Guard Array Accesses Before leaf_iter Increment

### Overview
The array accesses happen BEFORE the bounds-checked increment. We need to also guard those accesses.

### Changes Required:

#### 1. tronko-assign.c - Guard array accesses
**File**: `tronko-assign/tronko-assign.c`

For each section where arrays are accessed before increment:

**Lines 312-321** (discordant in concordant mode):
```c
                }else{
                    if (leaf_iter < MAX_NUM_BWA_MATCHES) {
                        results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
                        results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
                        if (use_leaf_portion==1){
                            results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
                            strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
                            if ( paired==1){
                                results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
                                strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
                            }
                        }
                    }
                }
```

Apply same pattern to:
- **Lines 345-355** (concordant matches)
- **Lines 380-391** (discordant in non-concordant mode)
- **Lines 418-427** (concordant in non-concordant mode)

### Success Criteria:

#### Automated Verification:
- [x] Clean build with no warnings
- [ ] Test with 16S_Bacteria dataset completes without crash

#### Manual Verification:
- [ ] Verify array bounds are never exceeded via debug logging

---

## Phase 4: Add TSV Log Enhancement

### Overview
Add BWA match context to TSV memory log format for post-hoc analysis.

### Changes Required:

#### 1. tsv_memlog.h - Add BWA-specific logging macro
**File**: `tronko-assign/tsv_memlog.h`

Add new macro for BWA context logging:
```c
#define TSV_LOG_BWA(fp, phase, leaf_iter, max_matches, concordant, discordant, dropped) \
    do { \
        if (fp) { \
            resource_stats_t __stats; \
            if (get_resource_stats(&__stats) == 0) { \
                fprintf(fp, "%.3f\t%s\t%zu\t%zu\t%zu\t%.2f\t%.2f\tleaf_iter=%d,max=%d,conc=%d,disc=%d,dropped=%d\n", \
                    __stats.wall_clock_sec, (phase), \
                    __stats.memory_rss_kb / 1024, \
                    __stats.memory_vm_size_kb / 1024, \
                    __stats.peak_rss_kb / 1024, \
                    __stats.cpu_user_sec, __stats.cpu_sys_sec, \
                    (leaf_iter), (max_matches), (concordant), (discordant), (dropped)); \
                fflush(fp); \
            } \
        } \
    } while(0)
```

### Success Criteria:

#### Automated Verification:
- [x] Clean build
- [ ] TSV log file contains BWA context when enabled

---

## Testing Strategy

### Unit Tests:
- Verify bounds checking prevents `leaf_iter > MAX_NUM_BWA_MATCHES`
- Verify dropped match counter increments correctly
- Verify crash context captures BWA state

### Integration Tests:
- Run single_tree dataset (baseline - should work as before)
- Run 16S_Bacteria dataset (should complete without segfault)
- Verify warning messages appear for high-match reads

### Manual Testing Steps:
1. Build with debug: `make clean && make debug`
2. Run single_tree test to verify no regression
3. Run 16S_Bacteria test with verbose mode: `-v 3`
4. Check for warning messages about dropped matches
5. If crash occurs, verify crash report includes BWA context

## Performance Considerations

- Bounds checking adds ~4 comparisons per unique tree match
- Logging only occurs when bounds exceeded (not on every match)
- No additional memory allocation
- Minimal impact on normal operation

## References

- Research: `thoughts/shared/research/2026-01-01-16S-bacteria-segfault-analysis.md`
- Bug report: `thoughts/shared/bugs/2026-01-01-16S-bacteria-segfault.md`
- Code: `tronko-assign/tronko-assign.c:298-450`
- Allocations: `tronko-assign/allocateMemoryForResults.c:12-48`
