# BWA Bounds Overflow Investigation

**Date**: 2026-01-01
**Status**: In Progress
**Branch**: experimental

## Summary

Implemented bounds checking for `leaf_iter` overflow in `tronko-assign.c`, but discovered the 16S_Bacteria crash has multiple root causes beyond the originally identified issue.

## Changes Made (Not Committed)

### 1. crash_debug.h - BWA Context Fields
Added new fields to `crash_info_t`:
```c
// BWA bounds tracking (added for 16S_Bacteria segfault analysis)
int bwa_leaf_iter;                // Current leaf_iter value
int bwa_max_matches;              // MAX_NUM_BWA_MATCHES constant
int bwa_concordant_count;         // Concordant matches for current read
int bwa_discordant_count;         // Discordant matches for current read
int bwa_unique_trees;             // Unique trees matched
int bwa_dropped_matches;          // Matches dropped due to bounds
```

Added function declarations:
```c
void crash_set_bwa_context(int leaf_iter, int concordant_count, int discordant_count);
void crash_set_bwa_bounds_violation(int leaf_iter, int max_matches, int dropped);
void crash_clear_bwa_context(void);
```

### 2. crash_debug.c - BWA Context Implementation
- Added corresponding fields to `g_app_context`
- Added BWA context capture in signal handler
- Added BWA context output in crash reports
- Implemented the three new tracking functions

### 3. tronko-assign.c - Bounds Checking
Added thread-local and global counters:
```c
static __thread int dropped_matches_count = 0;
static int g_overflow_read_count = 0;
static int g_total_dropped_matches = 0;
static int g_max_potential_matches = 0;
static pthread_mutex_t g_overflow_stats_mutex = PTHREAD_MUTEX_INITIALIZER;
```

Added bounds checking to all 4 locations where `leaf_iter` is incremented:
- Lines ~317-349: Discordant matches in concordant mode
- Lines ~356-390: Concordant matches in concordant mode
- Lines ~398-435: Discordant matches in non-concordant mode
- Lines ~441-479: Concordant matches in non-concordant mode

Pattern applied:
```c
if (leaf_iter < MAX_NUM_BWA_MATCHES) {
    // array access
    leaf_iter++;
} else {
    dropped_matches_count++;
}
```

Added diagnostic logging:
- Per-read warnings (first 100, then every 1000th overflow)
- End-of-run summary with total overflow count, dropped matches, max potential matches

### 4. tsv_memlog.h - BWA Logging Macro
Added `TSV_LOG_BWA` macro for structured BWA context logging.

### 5. fastmap.c - Unbounded Loop Fix (EXPERIMENTAL)
Added bounds check to while loops at lines 515 and 526:
```c
// Before:
while(strcmp(singleQueryMat->name[aux->startline+j],readname)!=0){

// After:
while(j < data->n_seqs && strcmp(singleQueryMat->name[aux->startline+j],readname)!=0){
```

## Test Results

| Test | Result |
|------|--------|
| single_tree dataset | ✅ PASS (165 results) |
| 16S_Bacteria dataset | ❌ CRASH |

## 16S_Bacteria Crash Analysis

### Timeline of Crashes

1. **Original crash** (before our changes): `fastmap.c:526` - unbounded `strcmp` loop
2. **After fastmap.c fix**: `bwa.c:394` - assertion failure in `bwa_idx2mem`

### Stack Trace (After fastmap.c Fix)
```
#6  bwa_idx2mem (bwa.c:394)
#7  ks_ksmall_mem_ars_hash (bwamem.c:398)
#8  mem_chain2aln (bwamem.c:708)
#9  mem_aln2sam (bwamem.c:922)
#10 __kb_delp_aux_chn (bwamem.c:187)
#11 sais_main (is.c:118)
```

### Crash Location
```c
// bwa.c:394
tmp += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
```
NULL pointer being passed to `strlen`.

### Memory Usage
- RSS: 16.2 GB at crash
- VM Size: 17.1 GB

### Hypothesis

The 16S_Bacteria dataset has multiple issues:

1. **Bounds overflow in fastmap.c** - FIXED (but reveals deeper issue)
2. **BWA index corruption or incompatibility** - The pre-built index may be:
   - Corrupted
   - Built with different BWA version
   - Incompatible with the reference file
3. **Memory pressure** - 16+ GB may cause issues with memory allocation

## Files Modified (Uncommitted)

```
tronko-assign/crash_debug.h
tronko-assign/crash_debug.c
tronko-assign/tronko-assign.c
tronko-assign/tsv_memlog.h
tronko-assign/bwa_source_files/fastmap.c
tronko-assign/bwa_source_files/bwamem.c (minor: pthread include)
```

## Next Steps

1. **Verify BWA index integrity** - Check if index files match the reference FASTA
2. **Test with freshly built BWA index** - Remove pre-built index and let tronko-assign rebuild
3. **Test with smaller 16S_Bacteria subset** - Determine if issue is size-related
4. **Investigate bwa_idx2mem** - Understand why `anns[i].name` is NULL
5. **Consider separating fixes**:
   - Commit leaf_iter bounds fix (validated on single_tree)
   - File separate ticket for fastmap.c and BWA issues

## Key Finding

The original plan assumed the segfault was caused by `leaf_iter` exceeding `MAX_NUM_BWA_MATCHES`. However, the 16S_Bacteria crash occurs **before** reaching the `leaf_iter` code - it happens during BWA's read processing in `fastmap.c` or during index operations in `bwa.c`.

The `leaf_iter` bounds fix is still valid and will protect against overflow when processing reads that match >10 unique trees, but we cannot test it with 16S_Bacteria until the BWA issues are resolved.

## Diagnostic Output Example (When Working)

When the bounds checking triggers, output will look like:
```
[WARN] Read example_read_123: 15 unique tree matches (capped at 10, dropped 5) [overflow #1]
...
[WARN] === BWA BOUNDS CAP SUMMARY ===
[WARN]   Reads hitting cap: 42
[WARN]   Total matches dropped: 187
[WARN]   Max potential matches seen: 23 (cap is 10)
[WARN]   Consider increasing MAX_NUM_BWA_MATCHES if accuracy is affected
```
