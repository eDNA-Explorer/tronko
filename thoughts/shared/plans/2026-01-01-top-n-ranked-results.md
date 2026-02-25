# Top-N Ranked Results Implementation Plan

## Overview

Add the ability for tronko-assign to output multiple ranked taxonomic assignments per query sequence instead of a single consensus result. This enables downstream analysis to consider alternative placements when the top result is uncertain.

## Current State Analysis

### How Results Are Currently Generated

1. **Score Calculation** (`placement.c:878-904`): Finds global maximum score across all BWA matches, trees, and nodes. Already computes scores for ALL nodes in `nodeScores[match][tree][node]`.

2. **Candidate Selection** (`placement.c:927-937`): Marks all nodes within Cinterval of maximum in `voteRoot[tree][node]` as candidates.

3. **Multi-Tree Consensus** (`tronko-assign.c:595-637`):
   - Computes LCA for each tree with candidates: `LCAs[i]`
   - Tracks which trees have candidates: `maxRoots[count]`
   - Climbs taxonomic levels until trees agree
   - **First tree wins** as tiebreaker when consensus reached

4. **Output** (`tronko-assign.c:643-767`): Single row per read with columns:
   ```
   Readname  Taxonomic_Path  Score  Forward_Mismatch  Reverse_Mismatch  Tree_Number  Node_Number
   ```

### Key Discovery: Infrastructure Already Exists

The code already computes per-tree results in `LCAs[]` and `maxRoots[]` arrays. Currently it discards all but the consensus. We simply need to output multiple tree results ranked by score.

## Desired End State

When invoked with `-N 3`:
```
Readname    Rank  Taxonomic_Path         Score      Fwd_Mismatch  Rev_Mismatch  Tree  Node
read_001    1     Ursus                  -1234.5    2             1             5     38    <- consensus (current behavior)
read_001    2     Ursus;maritimus        -1234.5    2             1             5     42    <- best tree's full taxonomy
read_001    3     Ursus;arctos           -1235.2    3             1             12    87    <- 2nd tree's full taxonomy
```

**Key Design Decision**:
- **Rank 1 = Consensus** (matches current output exactly for backward compatibility)
- **Rank 2+ = Per-tree results** sorted by best score, showing full species-level taxonomy

This means Rank 2 may be MORE specific than Rank 1 when trees disagree at species level.

**Verification**:
- `-N 1` produces identical output to current version (backward compatible)
- `-N 3` produces up to 3 rows per read: consensus + top-2 tree alternatives
- Existing pipelines continue to work unchanged

## What We're NOT Doing

1. **Not changing scoring algorithm** - Scores are already computed correctly
2. **Not changing Cinterval logic** - Only reporting what's already being calculated
3. **Not changing LCA calculation** - Per-tree LCAs remain the same
4. **Not adding new dependencies** - Pure C changes using existing data structures
5. **Not changing the output file format** - Same TSV, just with Rank column and potentially multiple rows

## Implementation Approach

**Consensus-First with Tree Alternatives**

For each read:
1. Compute consensus result using current algorithm (unchanged) → Rank 1
2. Track the best score per tree (already in `LCAs[]` / `maxRoots[]`)
3. Sort trees by their best node score
4. Output consensus as Rank 1, then top-(N-1) tree results as Rank 2+

This preserves backward compatibility while exposing the per-tree alternatives.

## Phase 1: Add CLI Option and Data Structures

### Overview
Add the `-N` / `--top-n` command-line option and necessary data structures to track per-tree scores.

### Changes Required:

#### 1. Options Structure
**File**: `tronko-assign/options.h`
**Changes**: Add `top_n` field to Options struct

```c
// Add after line ~238 in Options struct
int top_n;              // Number of results to output per read (default: 1)
```

#### 2. Options Parsing
**File**: `tronko-assign/options.c`
**Changes**: Add option parsing for `-N` / `--top-n`

In `long_options[]` array (after line ~48):
```c
{"top-n", required_argument, 0, 'N'},
```

In `usage[]` string (after line ~77):
```c
-N [INT], number of ranked results per read [default:1]\n\
```

In optstring (line ~109), add `N:`:
```c
c=getopt_long(argc,argv,"hpsrqw6yevUzP:75:f:u:t:m:d:o:x:g:1:2:a:c:n:3:4:C:L:N:V::l:RT",long_options, &option_index);
```

In switch statement (after line ~290):
```c
case 'N':
    success = sscanf(optarg, "%d", &(opt->top_n));
    if (!success || opt->top_n < 1 || opt->top_n > 10) {
        fprintf(stderr, "Invalid top-n value (1-10)\n");
        opt->top_n = 1;
    }
    break;
```

#### 3. Default Initialization
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Initialize default value

After line ~848:
```c
opt.top_n = 1;  // Default: single result (backward compatible)
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles without warnings: `cd tronko-assign && make clean && make`
- [ ] `./tronko-assign -h` shows new `-N` option in help text
- [ ] `./tronko-assign -N 3 ...` runs without error
- [ ] `./tronko-assign -N 0 ...` shows validation error and defaults to 1

#### Manual Verification:
- [ ] Help text is clear and matches existing style

---

## Phase 2: Track Per-Tree Best Scores

### Overview
Modify result construction to track the best score for each tree that has candidates, enabling ranking.

### Changes Required:

#### 1. Add Tree Score Tracking Structure
**File**: `tronko-assign/global.h`
**Changes**: Add structure to track per-tree results

After line ~188:
```c
typedef struct treeResult {
    int tree_number;        // Tree index
    int lca_node;           // LCA node within this tree
    type_of_PP best_score;  // Best node score in this tree
    int tax_index_0;        // Taxonomy species index
    int tax_index_1;        // Taxonomy level
} treeResult;
```

#### 2. Collect Per-Tree Scores in Result Construction
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Populate tree results array during LCA calculation

Replace lines 577-637 with enhanced version that tracks per-tree scores:

```c
// After line 576 (int numMinNodes=count;)
// Create array to hold per-tree results
treeResult tree_results[count > 0 ? count : 1];
int num_tree_results = 0;

// Collect per-tree LCAs and scores
for(i=0; i<count; i++){
    int tree_idx = maxRoots[i];
    LCAs[i] = getLCAofArray_Arr_Multiple(results->voteRoot[tree_idx], tree_idx, maxNumSpec, number_of_total_nodes);

    // Find best score in this tree for ranking
    type_of_PP best_tree_score = -9999999999999999;
    for(j=0; j<2*numspecArr[tree_idx]-1; j++){
        if(results->voteRoot[tree_idx][j] == 1){
            // Find the nodeScore for this node
            for(int m=0; m<leaf_iter; m++){
                if(results->leaf_coordinates[m][0] == tree_idx){
                    if(results->nodeScores[m][tree_idx][j] > best_tree_score){
                        best_tree_score = results->nodeScores[m][tree_idx][j];
                    }
                }
            }
        }
    }

    tree_results[num_tree_results].tree_number = tree_idx;
    tree_results[num_tree_results].lca_node = LCAs[i];
    tree_results[num_tree_results].best_score = best_tree_score;
    tree_results[num_tree_results].tax_index_0 = treeArr[tree_idx][LCAs[i]].taxIndex[0];
    tree_results[num_tree_results].tax_index_1 = treeArr[tree_idx][LCAs[i]].taxIndex[1];
    num_tree_results++;
}
```

#### 3. Sort Trees by Best Score
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Add comparison function and sort call

Before `runAssignmentOnChunk_WithBWA` function (around line 191):
```c
// Comparison function for sorting tree results by score (descending)
int compare_tree_results(const void *a, const void *b) {
    const treeResult *ta = (const treeResult *)a;
    const treeResult *tb = (const treeResult *)b;
    if (tb->best_score > ta->best_score) return 1;
    if (tb->best_score < ta->best_score) return -1;
    return 0;
}
```

After collecting tree results (in the new code block):
```c
// Sort by best score (descending - highest score first)
qsort(tree_results, num_tree_results, sizeof(treeResult), compare_tree_results);
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles without warnings: `make clean && make`
- [ ] Unit test with known input produces expected tree ordering

#### Manual Verification:
- [ ] Debug output shows trees sorted by score correctly

---

## Phase 3: Multi-Row Output Generation

### Overview
Modify output generation to produce multiple rows per read when N > 1.

### Changes Required:

#### 1. Pass top_n Through Thread Structure
**File**: `tronko-assign/global.h`
**Changes**: Add top_n to mystruct

In `mystruct` (after line ~158):
```c
int top_n;              // Number of results to output per read
```

**File**: `tronko-assign/tronko-assign.c`
**Changes**: Set top_n in thread initialization

After lines 1198-1199 (in the thread init loop):
```c
mstr[i].top_n = opt.top_n;
```

Also after lines 1465-1466 (paired-end thread init):
```c
mstr[i].top_n = opt.top_n;
```

#### 2. Modify taxonPath to Hold Multiple Results
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Allocate space for multiple results per read

In `runAssignmentOnChunk_WithBWA`, get top_n from mstr:
```c
int top_n = mstr->top_n;
```

Change the output buffer allocation to use a different approach - instead of storing multiple results in taxonPath, we'll store them as newline-separated results in the same string.

#### 3. Generate Multiple Output Rows
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Modify output formatting to include consensus as Rank 1, then alternatives

The key insight is that the current code already computes:
- `taxRoot`, `taxNode`, `taxIndex0`, `taxIndex1` for consensus (lines 629-632)
- `LCAs[]` and `maxRoots[]` for per-tree results

We need to:
1. Output current consensus logic as Rank 1 (preserve existing behavior)
2. Output sorted per-tree results as Rank 2+ (new)

**Pseudocode for the output section** (lines 643-767):

```c
char *readname = paired ? pairedQueryMat->forward_name[lineNumber] : singleQueryMat->name[lineNumber];
char *write_ptr = resultsPath;
int top_n = mstr->top_n;

if (count == 0) {
    // No matches - output unassigned (rank 1 only)
    sprintf(write_ptr, "%s\t1\tunassigned\t\t\t\t", readname);
} else {
    // === RANK 1: Consensus result (current behavior) ===
    char consensus_tax[MAXRESULTSNAME];
    int consensus_tree, consensus_node;

    if (count == 1) {
        // Single tree case - use that tree's LCA
        consensus_tree = maxRoot;
        consensus_node = LCA;
    } else {
        // Multi-tree case - use consensus (taxRoot/taxNode from existing logic)
        consensus_tree = taxRoot;
        consensus_node = taxNode;
    }

    // Build consensus taxonomy string (existing logic)
    build_taxonomy_path(consensus_tax, consensus_tree, consensus_node, taxIndex1);

    // Output Rank 1
    int written = sprintf(write_ptr, "%s\t1\t%s\t%lf\t%lf\t%lf\t%d\t%d",
        readname, consensus_tax, results->minimum[0],
        results->minimum[1], results->minimum[2],
        consensus_tree, consensus_node);
    write_ptr += written;

    // === RANK 2+: Per-tree alternatives (if N > 1) ===
    if (top_n > 1 && num_tree_results > 0) {
        // Sort tree_results by best_score descending
        qsort(tree_results, num_tree_results, sizeof(treeResult), compare_tree_results);

        int alternatives_to_output = (num_tree_results < top_n - 1) ? num_tree_results : (top_n - 1);

        for (int alt = 0; alt < alternatives_to_output; alt++) {
            treeResult *tr = &tree_results[alt];

            // Skip if this is the same as consensus (avoid duplicate)
            if (tr->tree_number == consensus_tree && tr->lca_node == consensus_node) {
                alternatives_to_output = (num_tree_results < top_n) ? num_tree_results : top_n;
                if (alt + 1 < num_tree_results) continue;
                else break;
            }

            char alt_tax[MAXRESULTSNAME];
            // Use full species-level taxonomy for alternatives (tax_index_1 from tree result)
            build_taxonomy_path(alt_tax, tr->tree_number, tr->lca_node, tr->tax_index_1);

            *write_ptr++ = '\n';
            written = sprintf(write_ptr, "%s\t%d\t%s\t%lf\t%lf\t%lf\t%d\t%d",
                readname, alt + 2, alt_tax, tr->best_score,
                results->minimum[1], results->minimum[2],
                tr->tree_number, tr->lca_node);
            write_ptr += written;
        }
    }
}

strcpy(results->taxonPath[iter], resultsPath);
```

Note: The actual implementation will integrate with the existing code structure, but this shows the logical flow.

#### 4. Update Output Header
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Add Rank column to header

Lines 1159 and 1420 - change:
```c
fprintf(results,"Readname\tRank\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles without warnings: `make clean && make`
- [ ] Output contains Rank column in header
- [ ] `-N 1` produces one row per read (backward compatible check)
- [ ] `-N 3` produces up to 3 rows per read

#### Manual Verification:
- [ ] Results are ranked correctly (best score = rank 1)
- [ ] Taxonomy paths are correct for each ranked result
- [ ] Existing pipelines still work with `-N 1` (default)

---

## Phase 4: Backward Compatibility and Edge Cases

### Overview
Handle edge cases and ensure full backward compatibility.

### Changes Required:

#### 1. Handle Fewer Results Than N Requested
Already handled in Phase 3: `output_count = min(num_tree_results, top_n)`

#### 2. Handle Duplicate Avoidance
When the consensus result is identical to the best-scoring tree's LCA (common in single-tree cases), we should not output the same result twice.

The Phase 3 pseudocode handles this with:
```c
// Skip if this is the same as consensus (avoid duplicate)
if (tr->tree_number == consensus_tree && tr->lca_node == consensus_node) {
    continue;  // Move to next alternative
}
```

This ensures that when N=3 is requested but only 2 unique results exist, we output 2 rows (not 3 with a duplicate).

#### 3. Ensure Buffer Size for Multiple Results
**File**: `tronko-assign/tronko-assign.c`
**Changes**: Increase resultsPath buffer to accommodate multiple results

```c
// Each result row is ~max_readname_length + max_lineTaxonomy + 120
// For N results, need N * (size + 1 for newline)
char* resultsPath = (char*)malloc(top_n * (max_readname_length + mstr->max_lineTaxonomy + 120) * sizeof(char));
```

Also update taxonPath allocation in the thread setup loops:
```c
for(k=0; k<end-start; k++){
    mstr[i].str->taxonPath[k] = malloc(top_n * (max_name_length + max_lineTaxonomy + 120) * sizeof(char));
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Compile and run test dataset with `-N 1` - compare output to baseline (should match exactly except Rank column)
- [ ] Run with `-N 5` on single-tree matches - should produce 1 row (can't exceed available trees)
- [ ] No memory leaks: `valgrind ./tronko-assign -N 3 ...`

#### Manual Verification:
- [ ] Edge cases work correctly:
  - [ ] Read with no matches: outputs "unassigned" with rank 1
  - [ ] Read matching only 1 tree: outputs 1 result regardless of N
  - [ ] Read matching 2 trees with N=5: outputs 2 results

---

## Testing Strategy

### Unit Tests

Create test script `test_top_n.sh`:
```bash
#!/bin/bash
# Test top-N functionality

TRONKO_ASSIGN="./tronko-assign"
REF_DB="path/to/test_reference.txt"
REF_FASTA="path/to/test_reference.fasta"
TEST_READS="path/to/test_reads.fasta"

# Test 1: N=1 (default) - backward compatibility
$TRONKO_ASSIGN -r -f $REF_DB -a $REF_FASTA -s -g $TEST_READS -o /tmp/n1.txt -N 1
# Count rows per read - should be 1

# Test 2: N=3
$TRONKO_ASSIGN -r -f $REF_DB -a $REF_FASTA -s -g $TEST_READS -o /tmp/n3.txt -N 3
# Verify ranks 1, 2, 3 present where multiple trees matched

# Test 3: Compare N=1 to baseline
# Remove Rank column and compare to previous version output
```

### Integration Tests

1. Run on example dataset:
```bash
cd tronko-assign
make clean && make
./tronko-assign -r -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_topn.txt -N 3
```

2. Verify output format and content

### Manual Testing Steps

1. [ ] Run with `-N 1` on a dataset and verify Rank 1 taxonomy matches current version exactly
2. [ ] Run with `-N 3` and verify:
   - [ ] Rank 1 = consensus result (same as current behavior)
   - [ ] Rank 2+ = per-tree alternatives sorted by score
   - [ ] No duplicate rows when consensus equals best tree's LCA
3. [ ] Test the polar bear / brown bear scenario:
   - [ ] If trees disagree at species but agree at genus, Rank 1 should be genus-level
   - [ ] Rank 2 should show the species from best-scoring tree
4. [ ] Check memory usage is reasonable (not N× increase)

## Performance Considerations

**Expected overhead: <1%**

- No additional score calculations (already computed for all nodes)
- Sorting is O(t log t) where t = number of trees with candidates (typically 1-10)
- Output writing is trivially more for N > 1 rows

Memory usage increase is negligible:
- Added `treeResult` array on stack: ~40 bytes × count (typically <20)
- Larger output buffer: ~N × 200 bytes (trivial)

## Migration Notes

### For Existing Pipelines

- **Default behavior unchanged**: `-N 1` is default, producing single-row output
- **New Rank column**: Existing parsers expecting column 2 = Taxonomic_Path need update
- **Suggested migration**: Update parsers to handle Rank column, or filter to rank=1

### Backward Compatibility Flag (Optional Future Work)

If needed, could add `--legacy-output` flag to suppress Rank column for full backward compatibility. Not recommended as it complicates codebase.

## References

- Research document: `thoughts/shared/research/2026-01-01-tie-handling-top-n-results.md`
- Score calculation: `tronko-assign/placement.c:878-942`
- Result construction: `tronko-assign/tronko-assign.c:552-767`
- Options parsing: `tronko-assign/options.c`
- Global structures: `tronko-assign/global.h:115-133`
