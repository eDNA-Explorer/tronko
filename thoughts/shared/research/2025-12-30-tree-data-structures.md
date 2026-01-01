---
date: 2025-12-30T11:00:00-07:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Tree Data Structures - Node Fields and Parent/Child Relationships"
tags: [research, codebase, tree-structures, node, phylogeny, tronko-assign]
status: complete
last_updated: 2025-12-30
last_updated_by: Claude
---

# Research: Tree Data Structures in tronko-assign

**Date**: 2025-12-30T11:00:00-07:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

What fields are stored per node? How are parent/child relationships represented?

## Summary

The tronko-assign module uses a flat array representation for binary phylogenetic trees. Each node stores parent/child indices (not pointers), depth, posterior probabilities for all alignment positions, and taxonomy indices. The naming convention is counterintuitive: `up[2]` stores **child** indices (toward leaves) while `down` stores the **parent** index (toward root).

## Detailed Findings

### Node Structure Definition

**File:** `tronko-assign/global.h:60-68`

```c
typedef struct node{
    int up[2];              // Child node indices (-1 for leaf, -2 uninitialized)
    int down;               // Parent node index (-1 for root, -2 uninitialized)
    int nd;                 // Node identifier
    int depth;              // Depth from root (0 at root, increases toward leaves)
    type_of_PP **posteriornc;  // 2D array: [numbase][4] for nucleotide posteriors
    char *name;             // Accession name (only allocated for leaf nodes)
    int taxIndex[2];        // Taxonomy coordinate: [species_index, rank_level]
}node;
```

### Field Details

| Field | Type | Purpose | Memory |
|-------|------|---------|--------|
| `up[2]` | `int[2]` | Child indices: `up[0]`=left, `up[1]`=right | 8 bytes |
| `down` | `int` | Parent index | 4 bytes |
| `nd` | `int` | Node identifier | 4 bytes |
| `depth` | `int` | Distance from root | 4 bytes |
| `posteriornc` | `type_of_PP**` | Log-likelihoods `[numbase][4]` | 8 + numbase*32 bytes |
| `name` | `char*` | Accession ID (leaf only) | 8 + MAX_NODENAME bytes |
| `taxIndex[2]` | `int[2]` | Taxonomy lookup coordinate | 8 bytes |

### Parent/Child Relationship Encoding

Despite the naming (`up` suggesting upward/parent, `down` suggesting downward/children), the convention is **inverted**:

```
                    Root (down = -1)
                   /              \
              up[0]               up[1]
               /                    \
          Child 0               Child 1
         /       \             /       \
     up[0]     up[1]      up[0]      up[1]
       ↓         ↓          ↓          ↓
    Leaf      Leaf       Leaf       Leaf
  (up[0]=-1, up[1]=-1)
```

**Sentinel Values:**
| Value | Meaning |
|-------|---------|
| `-2`  | Uninitialized (set during allocation) |
| `-1`  | Terminal: leaf node (no children) or root (no parent) |
| `>=0` | Valid node index within the tree array |

### Tree Storage Architecture

**Global Declaration** (`global.h:205`):
```c
extern node **treeArr;
```

The trees are stored as a 2D structure:
- `treeArr[partition]` - pointer to node array for partition/tree
- `treeArr[partition][node_index]` - specific node
- Each tree has `2*numspec - 1` nodes (standard binary tree)

**Supporting Arrays** (`global.h:44`):
```c
extern int *rootArr, *numspecArr, *numbaseArr;
```
- `rootArr[i]`: Root node index for tree `i`
- `numspecArr[i]`: Number of leaf nodes (species) in tree `i`
- `numbaseArr[i]`: Number of alignment positions in tree `i`

### Node Index Layout

Within each tree, nodes are indexed such that:
- **Internal nodes**: indices `0` to `numspec-2`
- **Leaf nodes**: indices `numspec-1` to `2*numspec-2`

Only leaf nodes have the `name` field allocated (`allocatetreememory.c:34-39`).

### Posterior Probability Storage

**Dimensions**: `posteriornc[numbase][4]`
- First index: alignment position (0 to numbase-1)
- Second index: nucleotide (0=A, 1=C, 2=G, 3=T)

**Data Type** (`global.h:14-22`):
```c
#ifdef OPTIMIZE_MEMORY
    #define type_of_PP float   // ~50% memory savings
#else
    #define type_of_PP double
#endif
```

**Usage in scoring** (`assignment.c:177-193`):
```c
if (locQuery[i]=='a' || locQuery[i]=='A'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][0];
}else if (locQuery[i]=='c' || locQuery[i]=='C'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][1];
}
// ... etc for G, T
```

### Taxonomy Index System

`taxIndex[2]` serves as a coordinate into the 4D `taxonomyArr`:
- `taxIndex[0]`: Species index (row in taxonomy table)
- `taxIndex[1]`: Taxonomic rank (0=species, 6=domain)

**Lookup Pattern**:
```c
taxonomyArr[tree][node.taxIndex[0]][node.taxIndex[1]]
// Returns the taxon name string
```

### Key Navigation Functions

**Leaf Detection**:
```c
if (up[0] == -1 && up[1] == -1) {
    // This is a leaf node
}
```

**Downward Traversal (root to leaves)** - `assignment.c:24-52`:
```c
void assignScores_Arr_paired(int rootNum, int node, ...) {
    int child0 = treeArr[rootNum][node].up[0];
    int child1 = treeArr[rootNum][node].up[1];
    if (child0 == -1 && child1 == -1) {
        // Leaf: compute score
    } else {
        // Internal: compute score, then recurse
        assignScores_Arr_paired(rootNum, child0, ...);
        assignScores_Arr_paired(rootNum, child1, ...);
    }
}
```

**Upward Traversal (LCA computation)** - `tronko-assign.c:98-107`:
```c
int getLCA_Arr(int node1, int node2, int whichRoot) {
    if (node1 == node2) return node1;
    if (treeArr[whichRoot][node1].depth > treeArr[whichRoot][node2].depth) {
        // Swap to ensure node2 is deeper
        int tmp = node1; node1 = node2; node2 = tmp;
    }
    node2 = treeArr[whichRoot][node2].down;  // Move up to parent
    return getLCA_Arr(node1, node2, whichRoot);
}
```

### leafMap Structure

**Definition** (`global.h:70-74`):
```c
typedef struct leafMap{
    char* name;   // Leaf node name (accession ID)
    int root;     // Which tree/partition this leaf belongs to
    int node;     // Node index within that tree
}leafMap;
```

Used as a hashmap for O(1) lookup from BWA alignment results to tree coordinates.

## Code References

- `tronko-assign/global.h:60-68` - Node struct definition
- `tronko-assign/global.h:70-74` - leafMap struct definition
- `tronko-assign/global.h:205` - treeArr declaration
- `tronko-assign/allocatetreememory.c:17-40` - Tree allocation function
- `tronko-assign/readreference.c:377-557` - Reference tree loading
- `tronko-assign/assignment.c:24-52` - Downward tree traversal (scoring)
- `tronko-assign/tronko-assign.c:98-107` - LCA calculation (upward traversal)
- `tronko-assign/tronko-assign.c:71-77` - Depth assignment function

## Architecture Insights

1. **Flat Array Representation**: Trees use index-based navigation rather than pointers, enabling efficient serialization/deserialization from the reference file.

2. **Counterintuitive Naming**: The `up`/`down` naming is reversed from typical conventions. This may be historical or reflect the original author's perspective (viewing the tree from leaves "up" toward root).

3. **Memory Optimization**: The `OPTIMIZE_MEMORY` flag allows trading precision for memory by using `float` instead of `double` for posteriors.

4. **Pre-computed Topology**: All parent/child relationships are serialized in the reference file, avoiding runtime tree construction.

5. **Separation of Concerns**:
   - `treeArr` holds structural/probability data
   - `taxonomyArr` holds taxonomy strings
   - `taxIndex` bridges the two

## Open Questions

1. Why is `nd` (node identifier) stored separately from the array index?
2. Could the `up`/`down` naming be refactored to `children`/`parent` for clarity?
3. Is there a memory savings opportunity by not allocating posteriors for leaf nodes (if they're never used)?
