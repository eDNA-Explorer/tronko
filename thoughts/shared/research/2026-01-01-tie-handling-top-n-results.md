---
date: 2026-01-01T12:00:00-08:00
researcher: Claude
git_commit: 8fa99b63db2052b40a62e61b9726a6bb6960c389
branch: experimental
repository: tronko
topic: "Tie Handling and Top-N Results Feasibility"
tags: [research, codebase, tie-handling, multiple-results, LCA, scoring, algorithm]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Tie Handling and Top-N Results Feasibility

**Date**: 2026-01-01T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 8fa99b63db2052b40a62e61b9726a6bb6960c389
**Branch**: experimental
**Repository**: tronko

## Research Question

How does tronko currently handle ties in taxonomic assignment? When different species (e.g., polar bear vs brown bear) appear in different trees, how is the result selected? Is it possible to return top 2-3 guesses instead of just one, and what would be the performance implications?

## Summary

**Current Behavior**: Tronko uses a Cinterval-based voting system where all nodes scoring within a configurable threshold of the maximum are considered tied candidates. These candidates are aggregated via LCA (Lowest Common Ancestor) calculation. When candidates span multiple trees, a consensus taxonomy is found by climbing taxonomic levels until agreement.

**Key Finding**: Returning top-N results is **algorithmically feasible with minimal performance overhead** because scores for ALL nodes are already computed and stored in the `nodeScores[match][tree][node]` 3D array. The infrastructure to identify multiple candidates exists via `voteRoot[][]`. The main work would be modifying the output format and adding sorting logic.

**Polar Bear vs Brown Bear Issue**: This occurs when BWA finds matches in different reference trees. Currently, when trees disagree at species level, the algorithm climbs to genus (or higher) and uses the **first tree's result** as a tiebreaker - not necessarily the best-scoring one.

## Detailed Findings

### Current Tie-Handling Mechanism

#### Step 1: Score Calculation (Comprehensive)
Scores are calculated for **every node** in the phylogenetic tree, not just the best candidate.

Location: `tronko-assign/assignment.c:24-93`

```c
void assignScores_Arr_paired(int rootNum, int node, ...) {
    // Recursive traversal calculates score for EVERY node
    type_of_PP node_score = getscore_Arr(alength, node, rootNum, locQuery, positions, ...);
    scores[search_number][rootNum][node] += node_score;

    // Recurse to children
    assignScores_Arr_paired(rootNum, child0, ...);
    assignScores_Arr_paired(rootNum, child1, ...);
}
```

#### Step 2: Maximum Score Identification
A single pass finds the global maximum across all BWA matches and all trees.

Location: `tronko-assign/placement.c:878-904`

```c
type_of_PP maximum = -9999999999999999;
for (i=0; i<number_of_matches; i++) {
    for (j=...; j<...; j++) {
        for(k=0; k<2*numspecArr[j]-1; k++) {
            if (maximum < nodeScores[i][j][k]) {
                maximum = nodeScores[i][j][k];
                // Track best: match_number, minRoot, minNode
            }
        }
    }
}
```

#### Step 3: Cinterval-Based Candidate Selection
All nodes within `Cinterval` of the maximum are marked as "tied" candidates.

Location: `tronko-assign/placement.c:927-937`

```c
for(i=0; i<number_of_matches; i++) {
    for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++) {
        if (nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) &&
            nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval)) {
            voteRoot[leaf_coordinates[i][0]][k] = 1;  // Mark as candidate
        }
    }
}
```

**Configuration**: `Cinterval` defaults to 5 (log-probability units), set via `-c` option.

#### Step 4: LCA Aggregation
Multiple candidates within a single tree are resolved via LCA calculation.

Location: `tronko-assign/tronko-assign.c:140-151`

```c
int LCA_of_nodes(int whichRoot, int root_node, int* minNodes, int numMinNodes) {
    // Finds deepest common ancestor of all minNodes
    getKeysCount(whichRoot, root_node, minNodes, matching_nodes, ancestors, numMinNodes);
    int LCA = ancestors[0];  // First (deepest) common ancestor
    return LCA;
}
```

#### Step 5: Multi-Tree Consensus Resolution
When candidates span multiple trees, taxonomy is compared at progressively higher levels.

Location: `tronko-assign/tronko-assign.c:595-637`

```c
// Compute LCA for each tree
for(i=0; i<count; i++) {
    LCAs[i] = getLCAofArray_Arr_Multiple(results->voteRoot[maxRoots[i]], ...);
}

// Climb taxonomic levels until agreement
while(stop==0 && minLevel<=6) {
    // Compare taxonomy at minLevel across all trees
    if (all_trees_agree) {
        taxRoot = maxRoots[0];  // FIRST tree wins as tiebreaker
        taxNode = LCAs[0];
    } else {
        minLevel++;  // Move to higher taxonomic level
    }
}
```

### Why Polar Bear vs Brown Bear Ends Up at Genus Level

When a query matches both a polar bear in Tree 5 and a brown bear in Tree 12:

1. BWA returns matches in both trees
2. Both trees may have scores within Cinterval of the global maximum
3. Each tree computes its own LCA (polar bear node vs brown bear node)
4. At species level (level 0), the taxonomies differ: "Ursus maritimus" vs "Ursus arctos"
5. The consensus algorithm climbs to genus level (level 1): both are "Ursus"
6. Consensus reached at genus; **result is "Ursus" from Tree 5** (first in iteration order)

**The "first tree wins" behavior** (line 629: `taxRoot=maxRoots[0]`) means the result depends on tree iteration order, not on which species actually scored higher.

### Data Structures Already Supporting Multiple Results

| Structure | Location | Contents |
|-----------|----------|----------|
| `nodeScores[match][tree][node]` | `global.h:119` | Scores for ALL nodes (already computed) |
| `voteRoot[tree][node]` | `global.h:120` | Binary flags for nodes within Cinterval |
| `minNodes[]` | `resultsStruct` | Linear array of candidate node indices |
| `LCAs[]` | Local variable | Per-tree LCAs |
| `LCAnames[][]` | `resultsStruct` | Taxonomy names for candidates |

### Existing Option for All Scores

The `-7` / `--print-all-scores` option already dumps all node scores to a file:

Location: `tronko-assign/placement.c:884-907`

```c
if (print_all_nodes == 1) {
    fprintf(node_scores_file, "Tree_Number\tNode_Number\tScore\n");
    // ... writes all nodeScores to file
}
```

This demonstrates the infrastructure exists but is not exposed as structured top-N output.

## Feasibility Analysis: Returning Top-N Results

### What Would Need to Change

1. **Add sorting logic** (~30 lines): After populating `nodeScores`, sort to identify top-N by score
2. **Modify output format** (~50 lines): Add columns or separate rows for alternative assignments
3. **Add command-line option** (~10 lines): `-n N` or `--top-n N` to specify how many results
4. **Handle multi-tree alternatives** (~20 lines): Include results from different trees that didn't win consensus

### Performance Impact

| Operation | Current Cost | Additional Cost for Top-N |
|-----------|--------------|---------------------------|
| BWA alignment | O(query_len × ref_size) | None |
| Score calculation | O(nodes × positions) | None - already computed |
| Find maximum | O(matches × trees × nodes) | None |
| Mark candidates | O(matches × nodes) | None |
| LCA calculation | O(nodes) per tree | None |
| **New: Sort for top-N** | N/A | O(n log n) where n = candidate nodes |
| **New: Output N rows** | N/A | O(N) |

**Estimated overhead: <1%** of total runtime. The expensive operations (BWA alignment, score calculation) are unchanged. Sorting the already-computed scores is trivial compared to sequence alignment.

### Implementation Approach

**Option A: Top-N Individual Nodes**
Return the N highest-scoring individual nodes (may include nodes from different trees):
```
Readname    Rank    Taxonomic_Path           Score       Tree    Node
read_001    1       Ursus;maritimus          -1234.5     5       42
read_001    2       Ursus;arctos             -1235.2     12      87
read_001    3       Ursus                    -1236.1     5       38
```

**Option B: Top-N Tree Results**
Return the best LCA from each of the top-N scoring trees:
```
Readname    Tree_Rank    Taxonomic_Path    Best_Score    Tree
read_001    1            Ursus;maritimus   -1234.5       5
read_001    2            Ursus;arctos      -1235.2       12
```

**Option C: Alternative Placements Within Cinterval**
Return all distinct taxonomic assignments from nodes within Cinterval:
```
Readname    Taxonomic_Path      Score_Range          Trees
read_001    Ursus;maritimus     -1234.5 to -1236.0   5,8
read_001    Ursus;arctos        -1235.2 to -1237.1   12
```

### Recommended Implementation

**Option B (Top-N Tree Results)** is most aligned with the current architecture:
- Already computes per-tree LCAs in `LCAs[]` array
- Already tracks which trees have candidates in `maxRoots[]`
- Simply needs to output all tree results instead of just consensus

Code change locations:
- `tronko-assign/tronko-assign.c:595-767` - modify result construction
- `tronko-assign/options.c` - add `-n` option
- `tronko-assign/tronko-assign.c:1159,1420` - modify output header

## Code References

- `tronko-assign/assignment.c:24-93` - Recursive score calculation for all nodes
- `tronko-assign/assignment.c:172-239` - Per-node likelihood calculation
- `tronko-assign/placement.c:878-904` - Maximum score identification
- `tronko-assign/placement.c:927-937` - Cinterval-based candidate marking
- `tronko-assign/tronko-assign.c:552-576` - Vote counting across trees
- `tronko-assign/tronko-assign.c:588-637` - Multi-tree consensus resolution
- `tronko-assign/tronko-assign.c:140-151` - LCA calculation
- `tronko-assign/global.h:115-133` - `resultsStruct` definition
- `tronko-assign/allocateMemoryForResults.c:12-28` - Score array allocation
- `tronko-assign/options.c:67` - Cinterval option definition

## Architecture Insights

1. **Comprehensive Scoring**: Unlike greedy algorithms, tronko scores ALL tree nodes, making top-N results essentially "free" from a computation standpoint.

2. **Two-Phase Design**: The separation of (1) score calculation and (2) candidate selection via Cinterval enables flexibility in how results are aggregated.

3. **Deterministic Tiebreaking**: When exact ties occur, the first-encountered node/tree wins based on iteration order (lower indices). This is deterministic but arbitrary.

4. **LCA Aggregation Philosophy**: The algorithm prefers returning a less-specific but confident assignment (genus) over a more-specific but uncertain one (species). Top-N results would allow users to see the species-level alternatives.

## Open Questions

1. **Output format preference**: Should top-N results be on separate rows (easier to parse) or columns (compact)?

2. **Score normalization**: Should alternative scores be reported as absolute log-probabilities or relative to the best score (delta)?

3. **Cinterval interaction**: Should top-N be based on absolute score ranking or limited to nodes within Cinterval?

4. **Multi-tree handling**: When reporting top-N, should each tree's best result count as one entry, or should individual nodes compete across trees?
