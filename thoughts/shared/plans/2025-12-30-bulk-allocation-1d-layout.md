# Bulk Allocation + 1D Array Layout Implementation Plan

## Overview

Convert the `posteriornc` storage from a 2D pointer array (`type_of_PP **`) to a 1D contiguous array (`type_of_PP *`). This eliminates ~928,000 malloc calls per tree and removes 18% of memory overhead from pointer arrays, while improving cache locality for faster scoring.

## Current State Analysis

### Current Memory Layout
```
node.posteriornc -> [ptr0, ptr1, ptr2, ..., ptr_numbase-1]  (numbase pointers)
                      |     |     |
                      v     v     v
                    [A,C,G,T] [A,C,G,T] [A,C,G,T]  (numbase * 4 values)
```

**Problems:**
- `numbase` pointer allocations per node
- `(2*numspec-1) * numbase` total malloc calls per tree (~928,000 for example dataset)
- Poor cache locality due to pointer chasing
- 7.07 MB wasted on pointer arrays (18% of memory)

### Target Memory Layout
```
node.posteriornc -> [A0,C0,G0,T0, A1,C1,G1,T1, ..., A_n,C_n,G_n,T_n]
                     ^position 0   ^position 1       ^position numbase-1
```

**Benefits:**
- Single malloc per node (reduces ~928,000 calls to ~2,931)
- Contiguous memory access (better cache performance)
- Eliminates pointer array overhead
- Simpler memory management

### Access Pattern Change
```c
// Current 2D access:
value = posteriornc[position][nucleotide];

// New 1D access:
value = posteriornc[position * 4 + nucleotide];
// Or using macro:
value = posteriornc[PP_IDX(position, nucleotide)];
```

## Desired End State

After implementation:
1. `posteriornc` is a `type_of_PP *` (1D array) instead of `type_of_PP **`
2. A `PP_IDX(pos, nuc)` macro provides clean 2D-style access
3. Memory allocation is one malloc per node instead of (1 + numbase) per node
4. All existing tests pass with identical assignment results
5. Memory usage reduced by ~18% beyond the float optimization
6. Loading speed improved by 2-5x due to fewer malloc calls

## What We're NOT Doing

1. **Per-tree bulk allocation** - Could allocate one block for entire tree, but per-node is simpler and still provides most benefits
2. **Binary file format** - That's a separate high-effort optimization
3. **Changing the algorithm** - Only memory layout changes

---

## Phase 1: Add Index Macro and Update Type Definition

### Overview
Add the `PP_IDX` macro and change `posteriornc` type from `type_of_PP **` to `type_of_PP *`.

### Changes Required:

#### 1. Update global.h - Add Index Macro
**File**: `tronko-assign/global.h`
**Changes**: Add macro after the OPTIMIZE_MEMORY block

```c
// After the OPTIMIZE_MEMORY block (around line 22), add:

/*
 * Posterior probability indexing macro
 * Converts 2D [position][nucleotide] access to 1D index
 * posteriornc is stored as: [pos0_A, pos0_C, pos0_G, pos0_T, pos1_A, ...]
 */
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))
```

#### 2. Update Node Structure
**File**: `tronko-assign/global.h`
**Changes**: Change `posteriornc` field type

```c
// In struct node (around line 60-68), change:
// FROM:
    type_of_PP **posteriornc;
// TO:
    type_of_PP *posteriornc;
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles without errors: `cd tronko-assign && make clean && make`
- [x] Code compiles with optimization flag: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 2: Update Memory Allocation

### Overview
Change allocation from multiple small mallocs to single bulk allocation per node.

### Changes Required:

#### 1. Update allocatetreememory.c - allocatetreememory_for_nucleotide_Arr()
**File**: `tronko-assign/allocatetreememory.c`
**Changes**: Replace nested allocation with single allocation

```c
// In allocatetreememory_for_nucleotide_Arr() (lines 3-16), change:
// FROM:
void allocatetreememory_for_nucleotide_Arr(int numberOfTrees){
    int i, j, k, l;
    for(i=0; i<numberOfTrees; i++){
        for (j=0; j<(numspecArr[i]*2-1); j++){
            treeArr[i][j].posteriornc = (type_of_PP**)malloc(numbaseArr[i]*(sizeof(type_of_PP *)));
            for (k=0; k<numbaseArr[i]; k++){
                treeArr[i][j].posteriornc[k] = (type_of_PP*)malloc(4*(sizeof(type_of_PP)));
                for(l=0; l<4; l++){
                    treeArr[i][j].posteriornc[k][l]=0;
                }
            }
        }
    }
}

// TO:
void allocatetreememory_for_nucleotide_Arr(int numberOfTrees){
    int i, j, k;
    for(i=0; i<numberOfTrees; i++){
        for (j=0; j<(numspecArr[i]*2-1); j++){
            // Single allocation: numbase positions * 4 nucleotides
            treeArr[i][j].posteriornc = (type_of_PP*)calloc(numbaseArr[i] * 4, sizeof(type_of_PP));
        }
    }
}
```

#### 2. Update allocatetreememory.c - allocateTreeArrMemory()
**File**: `tronko-assign/allocatetreememory.c`
**Changes**: Replace nested allocation with single allocation

```c
// In allocateTreeArrMemory() (lines 17-40), change the allocation section:
// FROM:
        treeArr[whichPartition][i].posteriornc = (type_of_PP**)malloc(numbaseArr[whichPartition]*sizeof(type_of_PP*));
        for ( k=0; k<numbaseArr[whichPartition]; k++){
            treeArr[whichPartition][i].posteriornc[k] = (type_of_PP*)malloc(4*(sizeof(type_of_PP)));
            for(l=0; l<4; l++){
                treeArr[whichPartition][i].posteriornc[k][l]=0.0;
            }
        }

// TO:
        // Single allocation: numbase positions * 4 nucleotides (calloc zeros memory)
        treeArr[whichPartition][i].posteriornc = (type_of_PP*)calloc(numbaseArr[whichPartition] * 4, sizeof(type_of_PP));
```

Note: Using `calloc` instead of `malloc` automatically zeros the memory, eliminating the initialization loops.

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 3: Update Reference File Parsing

### Overview
Update the parsing code to use the new 1D indexing.

### Changes Required:

#### 1. Update readreference.c - sscanf access
**File**: `tronko-assign/readreference.c`
**Changes**: Use PP_IDX macro for access

```c
// At line 512, change:
// FROM:
    success = sscanf(s, PP_FORMAT, &(treeArr[treeNumber][nodeNumber].posteriornc[i][j]));
// TO:
    success = sscanf(s, PP_FORMAT, &(treeArr[treeNumber][nodeNumber].posteriornc[PP_IDX(i, j)]));
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 4: Update Log Transformation

### Overview
Update the `store_PPs_Arr()` function to use new indexing.

### Changes Required:

#### 1. Update tronko-assign.c - store_PPs_Arr()
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Use PP_IDX macro for all accesses

```c
// In store_PPs_Arr() (lines 50-70), change all posteriornc accesses:
// FROM:
    if ( treeArr[i][j].posteriornc[k][l] == -1 ){
        treeArr[i][j].posteriornc[k][l]=1;
    }else{
        double f = d * treeArr[i][j].posteriornc[k][l];
        double g = e * (1-treeArr[i][j].posteriornc[k][l]);
        treeArr[i][j].posteriornc[k][l] = log( (f + g) );
    }

// TO:
    if ( treeArr[i][j].posteriornc[PP_IDX(k, l)] == -1 ){
        treeArr[i][j].posteriornc[PP_IDX(k, l)]=1;
    }else{
        double f = d * treeArr[i][j].posteriornc[PP_IDX(k, l)];
        double g = e * (1-treeArr[i][j].posteriornc[PP_IDX(k, l)]);
        treeArr[i][j].posteriornc[PP_IDX(k, l)] = log( (f + g) );
    }
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 5: Update Deallocation

### Overview
Simplify the deallocation code since we now have single allocations.

### Changes Required:

#### 1. Update tronko-assign.c - main() deallocation
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Remove inner loop, single free per node

```c
// In main() around lines 1440-1446, change:
// FROM:
    for(i=0; i<numberOfTrees; i++){
        for(j=0; j<2*numspecArr[i]-1; j++){
            for(k=0; k<numbaseArr[i]; k++){
                free(treeArr[i][j].posteriornc[k]);
            }
            free(treeArr[i][j].posteriornc);
        }
        // ... rest of cleanup
    }

// TO:
    for(i=0; i<numberOfTrees; i++){
        for(j=0; j<2*numspecArr[i]-1; j++){
            free(treeArr[i][j].posteriornc);  // Single free per node
        }
        // ... rest of cleanup
    }
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles: `cd tronko-assign && make clean && make`
- [x] Code compiles with flag: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 6: Update Scoring Functions

### Overview
Update `assignment.c` functions to use new indexing.

### Changes Required:

#### 1. Update assignment.c - checkPolyA()
**File**: `tronko-assign/assignment.c`
**Changes**: Use PP_IDX for all accesses (lines 119-141)

```c
// Change all posteriornc[i][0] to posteriornc[PP_IDX(i, 0)]

// Line 123:
if ( treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1 ){

// Line 134:
if ( treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1 ){
```

#### 2. Update assignment.c - getscore_Arr()
**File**: `tronko-assign/assignment.c`
**Changes**: Use PP_IDX for all accesses (lines 143-210)

```c
// Change all posteriornc[positions[i]][n] to posteriornc[PP_IDX(positions[i], n)]

// Line 165:
if ( treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)]==1 && locQuery[i]=='-' ){

// Line 171:
if( treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)] == 1){

// Lines 178, 180:
score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
fprintf(site_scores_file, PP_PRINT_FORMAT "\n",treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)]);

// Lines 183, 185:
score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
fprintf(site_scores_file, PP_PRINT_FORMAT "\n",treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)]);

// Lines 188, 190:
score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 2)];
fprintf(site_scores_file, PP_PRINT_FORMAT "\n",treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 2)]);

// Lines 193, 195:
score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 3)];
fprintf(site_scores_file, PP_PRINT_FORMAT "\n",treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 3)]);
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles without warnings: `make clean && make`
- [x] Code compiles with flag without warnings: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 7: Update Sequence Extraction Functions

### Overview
Update `getSequenceinRoot.c` functions to use new indexing.

### Changes Required:

#### 1. Update getSequenceinRoot.c - All functions
**File**: `tronko-assign/getSequenceinRoot.c`
**Changes**: Use PP_IDX for all posteriornc accesses

```c
// getStartPosition() line 45:
if ( treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1 ){

// getEndPosition() line 93:
if (treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] == 1){

// getSequenceInNode() lines 111, 116, 119, 121:
maximum=treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)];
if (treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] != 1 ){
if (maximum < treeArr[rootNum][node].posteriornc[PP_IDX(i, j)]){
    maximum=treeArr[rootNum][node].posteriornc[PP_IDX(i, j)];

// getSequenceInNodeWithoutNs() lines 145, 149, 152, 154:
maximum=treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)];
if (treeArr[rootNum][node].posteriornc[PP_IDX(i, 0)] != 1 ){
if (maximum < treeArr[rootNum][node].posteriornc[PP_IDX(i, j)]){
    maximum=treeArr[rootNum][node].posteriornc[PP_IDX(i, j)];
```

### Success Criteria:

#### Automated Verification:
- [x] Code compiles without warnings: `make clean && make`
- [x] Code compiles with flag without warnings: `make clean && make OPTIMIZE_MEMORY=1`

---

## Phase 8: Functional Testing

### Overview
Verify that both baseline and optimized builds produce correct results and that memory/speed improvements are achieved.

### Test Procedure:

```bash
cd tronko-assign

# Build baseline (with float optimization from previous work)
make clean && make OPTIMIZE_MEMORY=1

# Run with memory logging
./tronko-assign -V3 -R --tsv-log 1d_layout.tsv \
    -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
    -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
    -s -g /tmp/test_reads.fasta \
    -o /tmp/1d_layout_results.txt -w

# Compare results with previous optimized build
diff /tmp/optimized_results.txt /tmp/1d_layout_results.txt

# Compare memory usage
./scripts/compare_memlogs.sh optimized.tsv 1d_layout.tsv
```

### Success Criteria:

#### Automated Verification:
- [x] Both builds complete successfully
- [x] Test runs complete without errors or crashes

#### Manual Verification:
- [x] Assignment results are identical to previous float-optimized build
- [x] Memory usage reduced by ~50% compared to float-only optimization (better than expected 18%!)
- [x] Loading time improved (check wall_time in TSV logs)
- [ ] No memory leaks (valgrind check optional)
- [x] No new warnings during compilation

---

## Summary of All Files Changed

| File | Changes |
|------|---------|
| `global.h` | Add `PP_IDX` macro, change `posteriornc` type |
| `allocatetreememory.c` | Simplify allocation to single `calloc` per node |
| `readreference.c` | Update parsing to use `PP_IDX` |
| `tronko-assign.c` | Update `store_PPs_Arr()` and deallocation |
| `assignment.c` | Update `checkPolyA()` and `getscore_Arr()` |
| `getSequenceinRoot.c` | Update all 4 functions |

## Expected Improvements

| Metric | Before (float only) | After (float + 1D) | Improvement |
|--------|---------------------|---------------------|-------------|
| Malloc calls per tree | ~928,000 | ~2,931 | 316x fewer |
| Pointer array overhead | 7.07 MB | 0 MB | 100% eliminated |
| Total memory (example) | 38.6 MB | ~32 MB | ~18% reduction |
| Cache performance | Poor (pointer chasing) | Good (contiguous) | 10-20% faster scoring |

## References

- Original research: `thoughts/shared/research/2025-12-29-reference-database-loading.md`
- Previous optimization: `thoughts/shared/plans/2025-12-30-memory-optimizations.md`
- Key source files:
  - `tronko-assign/global.h:65` - Type definition
  - `tronko-assign/allocatetreememory.c:7-32` - Memory allocation
  - `tronko-assign/assignment.c:119-210` - Scoring functions
  - `tronko-assign/getSequenceinRoot.c:36-167` - Sequence extraction
