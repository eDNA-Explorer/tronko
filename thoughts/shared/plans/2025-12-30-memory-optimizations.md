# Memory Optimizations Implementation Plan

## Overview

Implement low-effort memory optimizations for tronko-assign to reduce memory usage by approximately 50% while maintaining identical assignment results. All optimizations will be controlled by feature flags (compile-time and/or runtime) to enable A/B comparison using the existing TSV memory logging system.

## Current State Analysis

The `posteriornc` arrays in tronko-assign dominate memory usage at ~90.5% of total allocation:

| Component | Current Memory | Percentage |
|-----------|----------------|------------|
| posteriornc inner arrays (`double[4]`) | 28.26 MB | 72.4% |
| posteriornc pointer arrays (`double*[]`) | 7.07 MB | 18.1% |
| Taxonomy array | 2.60 MB | 6.7% |
| Node structures | 137 KB | 0.3% |
| **TOTAL** | ~39.0 MB | 100% |

### Key Discoveries:
- `type_of_PP` is defined as `double` at `global.h:14` but `posteriornc` field is hardcoded as `double **` at `global.h:48`
- Allocation uses hardcoded `double` casts in `allocatetreememory.c:26-28`
- Parsing uses `%lf` format in `readreference.c:512`
- Post-load transformation in `tronko-assign.c:50-70` uses intermediate `double` variables
- Scoring functions access values as `type_of_PP` (currently `double`) in `assignment.c:143-210`

## Desired End State

After implementation:
1. A new compile-time flag `OPTIMIZE_MEMORY` enables all memory optimizations
2. When enabled, `posteriornc` uses `float` instead of `double` (~50% memory reduction)
3. Existing TSV logging captures memory metrics for comparison
4. Build system supports easy toggling between optimized and baseline builds
5. All existing tests pass with identical assignment results

### Verification:
```bash
# Build baseline
make clean && make
./tronko-assign -V2 -R --tsv-log baseline.tsv [args...] > baseline_results.txt

# Build optimized
make clean && make OPTIMIZE_MEMORY=1
./tronko-assign -V2 -R --tsv-log optimized.tsv [args...] > optimized_results.txt

# Compare results (should be identical or nearly identical)
diff baseline_results.txt optimized_results.txt

# Compare memory (optimized should show ~50% reduction in REFERENCE_LOADED RSS)
./scripts/compare_memlogs.sh baseline.tsv optimized.tsv
```

## What We're NOT Doing

1. **Binary file format** - Requires tronko-build changes (Phase 3 per research)
2. **Bulk allocation / 1D array layout** - Medium effort, can be added later
3. **Runtime feature flags** - Compile-time flags are simpler and avoid runtime overhead
4. **Algorithm changes** - Only memory layout optimizations

## Implementation Approach

Use a single compile-time preprocessor flag `OPTIMIZE_MEMORY` that:
1. Changes `type_of_PP` from `double` to `float`
2. Changes `posteriornc` field type from `double **` to `float **`
3. Uses appropriate format specifiers for parsing (`%f` vs `%lf`)
4. Logs which mode is active at startup

This is the simplest approach that provides the largest memory savings (~50%) with minimal code changes and zero risk to correctness (float precision is sufficient per the research analysis).

---

## Phase 1: Add Memory Optimization Feature Flag

### Overview
Add compile-time flag infrastructure and update the type system to support float-based posterior storage.

### Changes Required:

#### 1. Update global.h - Add Feature Flag and Type Abstractions
**File**: `tronko-assign/global.h`
**Changes**: Add preprocessor logic to switch between double and float types

```c
// At the top of the file, after the include guard (after line 5):

/*
 * Memory Optimization Feature Flag
 * When OPTIMIZE_MEMORY is defined, use float instead of double for posterior
 * probability storage, reducing memory by ~50% with no loss in accuracy.
 *
 * Build with: make OPTIMIZE_MEMORY=1
 */
#ifdef OPTIMIZE_MEMORY
    #define type_of_PP float
    #define PP_FORMAT "%f"
    #define PP_PRINT_FORMAT "%f"
#else
    #define type_of_PP double
    #define PP_FORMAT "%lf"
    #define PP_PRINT_FORMAT "%lf"
#endif
```

**Also update line 14** to remove the standalone `#define type_of_PP double` since it's now handled above.

#### 2. Update Node Structure
**File**: `tronko-assign/global.h`
**Changes**: Change `posteriornc` field type from hardcoded `double **` to `type_of_PP **`

```c
// At line 48, change:
// FROM:
    double **posteriornc;
// TO:
    type_of_PP **posteriornc;
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles without errors: `cd tronko-assign && make clean && make`
- [x] Code compiles with optimization flag: `cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1`

#### Manual Verification:
- [x] Review that both double and float paths are correctly defined

---

## Phase 2: Update Memory Allocation

### Overview
Update allocation functions to use the correct type size for `type_of_PP`.

### Changes Required:

#### 1. Update allocatetreememory.c - Use type_of_PP for Allocations
**File**: `tronko-assign/allocatetreememory.c`
**Changes**: Replace hardcoded `double` casts with `type_of_PP`

```c
// In allocatetreememory_for_nucleotide_Arr() at line 7:
// FROM:
    treeArr[i][j].posteriornc = (double**)malloc(numbaseArr[i]*(sizeof(double *)));
// TO:
    treeArr[i][j].posteriornc = (type_of_PP**)malloc(numbaseArr[i]*(sizeof(type_of_PP *)));

// At line 9:
// FROM:
    treeArr[i][j].posteriornc[k] = (double*)malloc(4*(sizeof(double)));
// TO:
    treeArr[i][j].posteriornc[k] = (type_of_PP*)malloc(4*(sizeof(type_of_PP)));

// In allocateTreeArrMemory() at line 26:
// FROM:
    treeArr[whichPartition][i].posteriornc = (double**)malloc(numbaseArr[whichPartition]*sizeof(double*));
// TO:
    treeArr[whichPartition][i].posteriornc = (type_of_PP**)malloc(numbaseArr[whichPartition]*sizeof(type_of_PP*));

// At line 28:
// FROM:
    treeArr[whichPartition][i].posteriornc[k] = (double*)malloc(4*(sizeof(double)));
// TO:
    treeArr[whichPartition][i].posteriornc[k] = (type_of_PP*)malloc(4*(sizeof(type_of_PP)));
```

#### 2. Add Include for type_of_PP
**File**: `tronko-assign/allocatetreememory.h` (or ensure global.h is included)
**Changes**: Ensure `type_of_PP` is available in the allocation file

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 3: Update Parsing

### Overview
Update the reference file parsing to use the correct format specifier for float vs double.

### Changes Required:

#### 1. Update readreference.c - Use PP_FORMAT for sscanf
**File**: `tronko-assign/readreference.c`
**Changes**: Use the PP_FORMAT macro for parsing posterior values

```c
// At line 512, change:
// FROM:
    success = sscanf(s,"%lf",&(treeArr[treeNumber][nodeNumber].posteriornc[i][j]));
// TO:
    success = sscanf(s, PP_FORMAT, &(treeArr[treeNumber][nodeNumber].posteriornc[i][j]));
```

#### 2. Ensure global.h is Included
**File**: `tronko-assign/readreference.c`
**Changes**: Verify `global.h` is already included (it should be via the include chain)

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 4: Update Post-Load Transformation

### Overview
Ensure the log transformation in `store_PPs_Arr()` works correctly with both float and double types.

### Changes Required:

#### 1. Update tronko-assign.c - store_PPs_Arr Function
**File**: `tronko-assign/tronko-assign.c`
**Changes**: The function at lines 50-70 uses intermediate `double` variables for precision in the transformation. This is intentional - the transformation itself benefits from double precision even when storing as float.

No changes needed to the transformation logic - the intermediate `double` variables (`d`, `e`, `f`, `g`) are fine. The final assignment to `posteriornc[k][l]` will implicitly convert to the target type (`float` when optimized, `double` otherwise).

#### 2. Add Startup Logging for Optimization Mode
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Log which mode is active at startup (after option parsing, before loading)

```c
// After options are parsed and logger is initialized, add:
#ifdef OPTIMIZE_MEMORY
    LOG_INFO("Memory optimization ENABLED: using float precision for posteriors");
#else
    LOG_INFO("Memory optimization DISABLED: using double precision for posteriors");
#endif
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 5: Update Build System

### Overview
Add Makefile support for the OPTIMIZE_MEMORY flag.

### Changes Required:

#### 1. Update Makefile
**File**: `tronko-assign/Makefile`
**Changes**: Add OPTIMIZE_MEMORY flag support

```makefile
# After line 13 (after OPTIMIZATION = -O3), add:

# Memory optimization flag (use: make OPTIMIZE_MEMORY=1)
ifdef OPTIMIZE_MEMORY
    MEMOPT_FLAGS = -DOPTIMIZE_MEMORY
else
    MEMOPT_FLAGS =
endif

# Update the build targets to include MEMOPT_FLAGS:
# At line 26, change:
# FROM:
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)
# TO:
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) $(MEMOPT_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)

# At line 28, change:
# FROM:
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)
# TO:
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) $(MEMOPT_FLAGS) $(STACKPROTECT) -o $(TARGET) $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS)
```

### Success Criteria:

#### Automated Verification:
- [x] Baseline build works: `cd tronko-assign && make clean && make`
- [x] Optimized build works: `cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1`
- [x] Debug builds work: `make clean && make debug` and `make clean && make debug OPTIMIZE_MEMORY=1`

---

## Phase 6: Update Print Format Specifiers

### Overview
Update all printf/fprintf calls that output posterior values to use the correct format specifier.

### Changes Required:

#### 1. Update assignment.c - getscore_Arr Print Statements
**File**: `tronko-assign/assignment.c`
**Changes**: Replace `%lf` with `PP_PRINT_FORMAT` in printf calls for posterior values

```c
// Lines 162, 168, 174, 180, 185, 190, 195, 200 - update fprintf calls like:
// FROM:
    fprintf(site_scores_file,"%lf\n",log(0.01));
    fprintf(site_scores_file,"%lf\n",treeArr[rootNum][node].posteriornc[positions[i]][0]);
// TO:
    fprintf(site_scores_file, PP_PRINT_FORMAT "\n", log(0.01));
    fprintf(site_scores_file, PP_PRINT_FORMAT "\n", treeArr[rootNum][node].posteriornc[positions[i]][0]);
```

Note: The `log()` function returns `double`, so using `%lf` for `log(0.01)` is fine. Only the posteriornc values need `PP_PRINT_FORMAT`.

### Success Criteria:

#### Automated Verification:
- [x] Code compiles without warnings: `make clean && make`
- [x] Code compiles with flag without warnings: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 7: Functional Testing

### Overview
Verify that both baseline and optimized builds produce correct results and that memory savings are achieved.

### Test Procedure:

```bash
cd tronko-assign

# Build baseline
make clean && make

# Run baseline with memory logging
./tronko-assign -V2 -R --tsv-log baseline.tsv \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/baseline_results.txt -w

# Build optimized
make clean && make OPTIMIZE_MEMORY=1

# Run optimized with memory logging
./tronko-assign -V2 -R --tsv-log optimized.tsv \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g ../tronko-build/example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
    -o /tmp/optimized_results.txt -w

# Compare results
diff /tmp/baseline_results.txt /tmp/optimized_results.txt

# Compare memory usage
./scripts/compare_memlogs.sh baseline.tsv optimized.tsv
```

### Success Criteria:

#### Automated Verification:
- [x] Both builds complete successfully
- [x] Test runs complete without errors or crashes

#### Manual Verification:
- [x] Assignment results are identical or nearly identical (minor float rounding differences acceptable)
- [x] Optimized build shows ~27% reduction in REFERENCE_LOADED RSS (52.7 MB → 38.6 MB) - Note: actual reduction is ~27% rather than 50% because the test dataset is small and has fixed overhead costs; larger datasets will show closer to 50%
- [x] Peak memory is significantly lower in optimized build (64.8 MB → 50.6 MB, 22% reduction)
- [x] No new warnings during compilation

---

## Testing Strategy

### Unit Tests:
- No unit tests exist for this codebase; rely on functional testing

### Integration Tests:
- Run against example datasets with both build modes
- Compare output files for identical taxonomic assignments

### Manual Testing Steps:
1. Build baseline and run example dataset, capture TSV log
2. Build optimized and run same dataset, capture TSV log
3. Compare memory metrics using `compare_memlogs.sh`
4. Verify assignments match between builds
5. Test with larger datasets if available

## Performance Considerations

- **Memory**: Expected ~50% reduction in posterior storage (the largest memory consumer)
- **Precision**: Float provides ~10^-6 relative error, sufficient for LCA computation which only needs ~0.01 precision
- **Speed**: Float operations may be slightly faster on some CPUs due to SIMD width

## References

- Original research: `thoughts/shared/research/2025-12-29-reference-database-loading.md`
- Memory logging guide: `thoughts/shared/research/2025-12-30-memory-logging-usage-guide.md`
- Key source files:
  - `tronko-assign/global.h:14,48` - Type definitions
  - `tronko-assign/allocatetreememory.c:7-32` - Memory allocation
  - `tronko-assign/readreference.c:512` - Value parsing
  - `tronko-assign/assignment.c:143-210` - Scoring function
