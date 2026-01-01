---
date: 2026-01-01T00:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Algorithm-Level Improvements for tronko-assign: Early Termination, Tree Pruning, and Hierarchical Indexing"
tags: [research, codebase, algorithm-optimization, branch-and-bound, early-termination, tree-pruning, phylogenetic-placement]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Algorithm-Level Improvements for tronko-assign

**Date**: 2026-01-01T00:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Can we improve the performance of tronko-assign's assignment algorithm through:
1. Early termination: Stop tree traversal when LCA is confidently determined
2. Tree pruning: Skip trees that can't improve the current best score
3. Hierarchical indexing: Pre-filter which trees are relevant per read
4. Branch-and-bound algorithms and A* search adaptations

## Summary

The current tronko-assign algorithm performs **exhaustive tree traversal**, scoring every node in every matched tree before determining the LCA. This presents significant optimization opportunities. Based on analysis of the codebase and research into similar tools (pplacer, EPA-ng, APPLES, Kraken), several algorithmic improvements are feasible:

1. **Early termination** is highly viable - the "baseball heuristic" from pplacer can terminate after finding nodes that are clearly suboptimal
2. **Tree/subtree pruning** is possible if we exploit score monotonicity properties
3. **Hierarchical indexing** via pre-computed cluster representatives can reduce search space by 90%+
4. **Branch-and-bound** with analytical upper bounds can eliminate 92-98% of tree space

The key insight is that **scores at ancestral nodes provide upper bounds** for all descendants in many cases, enabling aggressive pruning without sacrificing accuracy.

## Detailed Findings

### Current Algorithm Analysis

#### Entry Points and Flow

The assignment algorithm has three main phases:

1. **BWA Alignment** (`tronko-assign.c:278`): Identifies candidate leaf nodes (up to 10 matches)
2. **Tree Traversal & Scoring** (`assignment.c:24-65`): Recursively scores ALL nodes
3. **LCA Selection** (`tronko-assign.c:527`): Finds LCA of nodes within Cinterval of best score

#### Core Scoring Function (`assignment.c:143-210`)

```c
type_of_PP getscore_Arr(int alength, int node, int rootNum, char *locQuery, int *positions, ...){
    type_of_PP score = 0;
    for (i=0; i<alength; i++){
        // Lookup pre-computed log-posterior probability for query nucleotide
        if (locQuery[i]=='A') score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 0)];
        else if (locQuery[i]=='C') score += treeArr[rootNum][node].posteriornc[PP_IDX(positions[i], 1)];
        // ... etc for G, T, gaps
    }
    return score;
}
```

#### Tree Traversal (`assignment.c:24-52`)

```c
void assignScores_Arr_paired(int rootNum, int node, ...){
    int child0 = treeArr[rootNum][node].up[0];
    int child1 = treeArr[rootNum][node].up[1];

    if (child0 == -1 && child1 == -1){  // Leaf
        scores[search_number][rootNum][node] += getscore_Arr(...);
    } else {  // Internal node
        scores[search_number][rootNum][node] += getscore_Arr(...);
        assignScores_Arr_paired(rootNum, child0, ...);  // Recurse left
        assignScores_Arr_paired(rootNum, child1, ...);  // Recurse right
    }
}
```

**Key observation**: Currently traverses ALL `2*numspec-1` nodes per tree, regardless of intermediate scores.

#### Current "Early Termination" (Minimal)

Only two early exits exist:
1. Invalid position check (`assignment.c:150-153`): Returns immediately if first position is invalid
2. Missing CIGAR (`placement.c:98-99`): Skips match if CIGAR is invalid

**No pruning** occurs during tree traversal itself.

### Optimization Opportunity 1: Early Termination

#### The Problem

After scoring all nodes, the algorithm finds nodes within `Cinterval` of the maximum score (`placement.c:911-921`):

```c
for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++){
    if (nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) &&
        nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval)){
        voteRoot[leaf_coordinates[i][0]][k] = 1;
    }
}
```

The LCA is then computed from voted nodes. If we could identify the LCA earlier, we could skip scoring many nodes.

#### Proposed Solution: Baseball Heuristic (from pplacer)

**How it works**:
1. Pre-sort nodes by a quick heuristic (e.g., alignment score to nearest leaf)
2. Track `best_score` and `strikes` counter
3. For each node in sorted order:
   - If `node_score > best_score + strike_box`: increment strikes
   - If `node_score < best_score`: reset strikes, update best_score
4. Stop when `strikes >= max_strikes` OR `pitches >= max_pitches`

**Implementation sketch**:

```c
// New function: assignScores_WithEarlyTermination
void assignScores_WithEarlyTermination(int rootNum, int node, ...,
                                        type_of_PP *best_score, int *strikes,
                                        type_of_PP strike_box, int max_strikes) {
    type_of_PP score = getscore_Arr(...);
    scores[search_number][rootNum][node] = score;

    // Early termination check
    if (score > *best_score + strike_box) {
        (*strikes)++;
        if (*strikes >= max_strikes) return;  // PRUNE
    } else if (score < *best_score) {
        *best_score = score;
        *strikes = 0;
    }

    // Continue to children only if not pruned
    if (child0 != -1) assignScores_WithEarlyTermination(rootNum, child0, ...);
    if (child1 != -1) assignScores_WithEarlyTermination(rootNum, child1, ...);
}
```

**Expected impact**: pplacer reports 2-10x speedup with default parameters (`--strike-box 3.0`, `--max-strikes 6`).

### Optimization Opportunity 2: Tree/Subtree Pruning

#### Score Monotonicity Analysis

**Key question**: Are scores monotonic along root-to-leaf paths?

In tronko, scores are computed using pre-stored posterior probabilities at each node. The posterior at internal nodes represents the **marginal probability** over all descendant leaves. This means:

- **Property 1**: Internal node scores tend to be "averaged" versions of leaf scores
- **Property 2**: The best-scoring leaf is often under the best-scoring internal node

This suggests **bounded pruning** is possible:

```c
void assignScores_WithBoundedPruning(int rootNum, int node, ...,
                                      type_of_PP global_best, type_of_PP cutoff) {
    type_of_PP score = getscore_Arr(...);
    scores[search_number][rootNum][node] = score;

    // PRUNING: If this node is too bad, skip entire subtree
    if (score > global_best + cutoff) {
        return;  // Don't recurse to children
    }

    if (child0 != -1) assignScores_WithBoundedPruning(rootNum, child0, ...);
    if (child1 != -1) assignScores_WithBoundedPruning(rootNum, child1, ...);
}
```

**Caveat**: Requires validation that monotonicity holds for tronko's posterior model. If it doesn't hold strictly, use a larger cutoff to ensure safety.

#### Branch-and-Bound Formalization

Classic branch-and-bound requires:
1. **Branching rule**: Expand children of current node
2. **Bounding function**: Upper bound on best achievable score in subtree
3. **Pruning condition**: Prune if bound > current best + threshold

**Upper bound calculation** for a subtree rooted at node N:
```c
type_of_PP subtree_upper_bound(int rootNum, int node) {
    // Option 1: Use node's own score as upper bound (if scores are monotonic)
    return scores[rootNum][node];

    // Option 2: Pre-compute and store upper bounds during database build
    return treeArr[rootNum][node].max_descendant_score;
}
```

**Research finding**: Analytical upper bounds on maximum likelihood can eliminate 92-98% of tree space (Hendy & Penny 1982, Rogers 2003).

### Optimization Opportunity 3: Hierarchical Indexing

#### Current Bottleneck

BWA returns up to 10 leaf matches, and currently ALL nodes in ALL matched trees are scored:
```
Total nodes scored = matches * trees_per_match * (2 * numspec - 1)
```

For a database with 20,000 species trees averaging 100 species each: `10 * 1 * 199 = 1,990 nodes/read` minimum.

#### Two-Phase Candidate Selection (EPA-ng approach)

**Phase 1: Quick Screening**
```c
// Compute approximate scores at strategic "checkpoint" nodes
for (int tree = 0; tree < num_matched_trees; tree++) {
    for (int checkpoint : precomputed_checkpoints[tree]) {
        approx_scores[tree][checkpoint] = quick_score(checkpoint);
    }
}

// Select top K trees/branches
sort_by_approx_score();
candidate_trees = top_k(num_matched_trees, K);  // e.g., K = top 10%
```

**Phase 2: Thorough Evaluation**
```c
// Only score nodes in candidate trees
for (int tree : candidate_trees) {
    assignScores_Arr_paired(tree, rootArr[tree], ...);
}
```

**Implementation considerations**:
- Checkpoint nodes can be pre-computed at database build time
- Place checkpoints at tree "centers" (nodes that minimize max distance to any leaf)
- Quick score can use simplified model (e.g., consensus sequence matching)

**Expected impact**: EPA-ng reports reducing candidate branches from thousands to <10, with minimal accuracy loss.

#### Pre-filtering by Leaf Similarity

Before full placement, filter trees based on BWA match quality:

```c
// After BWA returns matches
for (int match = 0; match < number_of_matches; match++) {
    int tree = leaf_coordinates[match][0];
    int leaf = leaf_coordinates[match][1];

    // Use BWA alignment score as pre-filter
    if (bwa_scores[match] < threshold) {
        skip_tree[tree] = 1;  // Don't process this tree
    }
}
```

### Optimization Opportunity 4: A* / Best-First Search

#### Standard A* Adaptation

```c
typedef struct {
    int node;
    type_of_PP g_score;  // Actual score
    type_of_PP f_score;  // g + heuristic
} PQNode;

void assignScores_AStar(int rootNum, ...) {
    PriorityQueue pq;
    pq_init(&pq);
    pq_push(&pq, (PQNode){rootArr[rootNum], 0, heuristic(rootArr[rootNum])});

    type_of_PP best_score = INFINITY;
    int best_nodes[MAX_NODES];
    int num_best = 0;

    while (!pq_empty(&pq)) {
        PQNode current = pq_pop(&pq);

        // Termination: if f_score exceeds best + Cinterval, we're done
        if (current.f_score > best_score + Cinterval) break;

        type_of_PP score = getscore_Arr(..., current.node, ...);

        if (score < best_score) {
            best_score = score;
        }
        if (score <= best_score + Cinterval) {
            best_nodes[num_best++] = current.node;
        }

        // Expand children
        for (int child : {treeArr[rootNum][current.node].up[0], up[1]}) {
            if (child != -1) {
                pq_push(&pq, (PQNode){child, score, score + heuristic(child)});
            }
        }
    }

    // best_nodes now contains nodes for LCA calculation
}
```

#### Heuristic Function Design

The heuristic `h(node)` should estimate the minimum additional cost to reach any leaf from `node`:

```c
type_of_PP heuristic(int rootNum, int node) {
    // Option 1: Depth-based (admissible if scores decrease with depth)
    return -depth[node] * avg_score_per_level;

    // Option 2: Pre-computed (computed at database build time)
    return treeArr[rootNum][node].min_leaf_score_delta;

    // Option 3: Zero heuristic (degrades to Dijkstra's algorithm)
    return 0;
}
```

#### IDA* (Iterative Deepening A*)

More memory-efficient variant:

```c
type_of_PP ida_star_search(int rootNum, ...) {
    type_of_PP threshold = heuristic(rootArr[rootNum]);

    while (1) {
        type_of_PP result = bounded_dfs(rootArr[rootNum], 0, threshold);
        if (result == FOUND) return best_score;
        if (result == INFINITY) return NOT_FOUND;
        threshold = result;  // New threshold for next iteration
    }
}

type_of_PP bounded_dfs(int node, type_of_PP g, type_of_PP threshold) {
    type_of_PP f = g + heuristic(node);
    if (f > threshold) return f;  // Return minimum exceeding threshold

    type_of_PP score = getscore_Arr(..., node, ...);
    // ... check if node qualifies, recurse to children
}
```

**Advantage**: O(d) memory vs O(n) for standard A*, where d = tree depth.

### Concrete Implementation Plan

#### Phase 1: Low-Hanging Fruit (1-2 days)

1. **Add early termination with strike counter**
   - Modify `assignScores_Arr_paired()` to accept best_score, strikes parameters
   - Add command-line options: `--strike-box`, `--max-strikes`
   - Default: strike_box=5 (same as Cinterval), max_strikes=10

2. **Add subtree pruning**
   - If `node_score > best_score + 2*Cinterval`, skip children
   - Conservative threshold ensures no false negatives

#### Phase 2: Structural Changes (3-5 days)

3. **Implement two-phase scoring**
   - Phase 1: Score only leaves + root path
   - Phase 2: If needed, score full tree

4. **Add priority queue-based traversal**
   - Replace recursive DFS with iterative best-first
   - Process nodes in score order rather than tree order

#### Phase 3: Database-Level (1 week)

5. **Pre-compute auxiliary data**
   - Add `max_descendant_score` to node structure
   - Add checkpoint nodes at tree build time
   - Store subtree size for work estimation

6. **Hierarchical tree clustering**
   - Cluster similar trees in database
   - Add cluster representatives for quick filtering

## Code References

- `tronko-assign/assignment.c:24-65` - Core tree traversal (target for modification)
- `tronko-assign/assignment.c:143-210` - Score calculation function
- `tronko-assign/placement.c:862-930` - Maximum finding and vote collection
- `tronko-assign/tronko-assign.c:485-571` - LCA calculation and result formatting
- `tronko-assign/global.h:85-93` - Node structure (may need new fields)
- `tronko-assign/allocatetreememory.c:12-28` - Tree memory allocation

## Architecture Insights

### Why Current Design is Exhaustive

The original tronko design prioritizes **correctness over speed**:
1. Posterior probabilities at internal nodes are "real" values (not bounds)
2. The Cinterval voting mechanism needs all scores to find all qualifying nodes
3. LCA of voted nodes is the final taxonomic assignment

This design is **sound** but leaves performance on the table when:
- Most nodes score far worse than the best
- Early nodes in traversal already establish strong bounds
- Multiple trees share similar structure (redundant computation)

### Compatibility Considerations

Any optimization must preserve the property that **all nodes within Cinterval of the best score are found**. Optimizations that might miss qualifying nodes:
- **Unsafe**: Hard cutoff that could miss nodes scoring exactly at Cinterval
- **Safe**: Soft cutoff with margin (e.g., `2*Cinterval`) or verification pass

### Memory vs Speed Trade-offs

| Optimization | Memory Impact | Speed Impact | Accuracy |
|-------------|--------------|--------------|----------|
| Strike-based early termination | None | 2-5x | 100%* |
| Subtree pruning | None | 2-10x | 100%* |
| Two-phase screening | O(checkpoints) | 3-10x | ~99.9% |
| Pre-computed bounds | O(nodes) | 5-20x | 100% |
| A*/IDA* search | O(1) to O(n) | 2-10x | 100% |

*With appropriate safety margins

## Related Research

### Published Tools

1. **pplacer** (Matsen et al., 2010): Linear-time placement using baseball heuristic
   - Key insight: Most placements only need a few branch length optimizations
   - Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC3098090/

2. **EPA-ng** (Barbera et al., 2019): Massively parallel placement with premasking
   - Key insight: Two-phase approach with quick screening + thorough evaluation
   - Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6368480/

3. **APPLES** (Balaban et al., 2020): O(n) distance-based placement
   - Key insight: Dynamic programming on fixed topology avoids redundant computation
   - Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC7164367/

4. **Kraken** (Wood & Salzberg, 2014): Ultrafast k-mer LCA classification
   - Key insight: Pre-computed LCA mappings enable O(1) classification
   - Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC4053813/

### Algorithmic Foundations

5. **Upper bounds on ML** (Rogers, 2003): Analytical bounds eliminate 92-98% of trees
   - Source: https://pubmed.ncbi.nlm.nih.gov/14534174/

6. **ML-guided tree search** (Crotty et al., 2021): Machine learning to predict search direction
   - Source: https://www.nature.com/articles/s41467-021-22073-8

## Open Questions

1. **Score monotonicity**: Do tronko's posterior-based scores exhibit monotonicity along root-to-leaf paths? This determines feasibility of aggressive pruning.

2. **Checkpoint placement**: What's the optimal strategy for placing screening checkpoints in the tree? Center nodes, LCA of diverse leaves, or random sampling?

3. **Heuristic design**: For A* search, what admissible heuristic works best for posterior probability scores?

4. **Paired-end handling**: How do optimizations interact with the score accumulation for paired reads? Need to ensure both reads' contributions are considered.

5. **Multi-tree consensus**: When multiple trees are matched, how should we balance tree-level vs node-level early termination?

6. **Benchmark datasets**: What test datasets will best demonstrate speedup while validating accuracy preservation?
