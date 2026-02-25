---
date: 2025-12-31T00:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Reference Tree Database Format"
tags: [research, codebase, reference-tree, file-format, tronko-build, tronko-assign]
status: complete
last_updated: 2025-12-31
last_updated_by: Claude
---

# Research: Reference Tree Database Format

**Date**: 2025-12-31T00:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Understanding the `reference_tree.txt` file format - how it's structured, written by tronko-build, and read by tronko-assign.

## Summary

The reference tree is a custom binary phylogenetic database format that stores:
1. Tree topology (parent/child relationships)
2. Taxonomy mappings for each species
3. Posterior probability matrices for each node at each MSA position

The format is written by `tronko-build/printtree.c` and read by `tronko-assign/readreference.c`.

## Detailed Findings

### File Structure Overview

The file is a plain text format with the following sections:

```
[Global Header - 4 lines]
[Per-Tree Headers - 1 line per tree]
[Taxonomy Section - numspec lines per tree]
[Node Data Section - (2*numspec-1) nodes per tree, each with posteriors]
```

### Global Header (Lines 1-4)

```
Line 1: numberOfTrees     (int) - number of phylogenetic trees in database
Line 2: max_nodename      (int) - maximum length of node/accession names
Line 3: max_tax_name      (int) - maximum length of taxonomy field names
Line 4: max_lineTaxonomy  (int) - maximum line length for taxonomy entries
```

### Per-Tree Header

One line per tree, tab-separated:
```
numbase<TAB>root<TAB>numspec
```

| Field | Type | Description |
|-------|------|-------------|
| numbase | int | Number of base positions in MSA (alignment length) |
| root | int | Index of root node in tree |
| numspec | int | Number of species (leaf nodes) |

### Taxonomy Section

For each tree, `numspec` lines of semicolon-delimited taxonomy:
```
Species;Genus;Family;Order;Class;Phylum;Domain
```

Example:
```
Larus heuglini;Larus;Laridae;Charadriiformes;Aves;Chordata;Eukaryota
```

The taxonomy array is indexed as `taxonomyArr[tree][species][rank]` where rank is 0-6:
- 0: Species
- 1: Genus
- 2: Family
- 3: Order
- 4: Class
- 5: Phylum
- 6: Domain

### Node Data Section

For each tree, there are `2*numspec - 1` nodes (binary tree property).

Each node consists of:

**Node Header Line** (tab-separated):
```
treeNum<TAB>nodeNum<TAB>up[0]<TAB>up[1]<TAB>down<TAB>depth<TAB>taxIndex[0]<TAB>taxIndex[1]<TAB>[name]
```

| Field | Description |
|-------|-------------|
| treeNum | Tree index (0-based) |
| nodeNum | Node index within tree |
| up[0], up[1] | Child node indices (-1 for leaf nodes) |
| down | Parent node index (-1 for root) |
| depth | Depth in tree from root |
| taxIndex[0] | Species index into taxonomy array |
| taxIndex[1] | Rank index (0-6) |
| name | Accession name (only present for leaf nodes where up[0]==-1 && up[1]==-1) |

**Posterior Probability Lines** (`numbase` lines per node):
```
P(A)<TAB>P(C)<TAB>P(G)<TAB>P(T)
```

Posteriors are written with `%.17g` format (17 significant digits).

### Node Structure Definition

From `tronko-assign/global.h:60-68`:

```c
typedef struct node {
    int up[2];           // child node indices (-1 for leaf nodes)
    int down;            // parent node index (-1 for root)
    int nd;              // node descriptor
    int depth;           // depth in tree
    type_of_PP **posteriornc;  // [numbase][4] posterior probs (A,C,G,T)
    char *name;          // accession name (only for leaf nodes)
    int taxIndex[2];     // indices into taxonomy array
} node;
```

### Data Types

From `tronko-assign/global.h:14-21`:

```c
#ifdef OPTIMIZE_MEMORY
    #define type_of_PP float
    #define PP_FORMAT "%f"
#else
    #define type_of_PP double
    #define PP_FORMAT "%lf"
#endif
```

### Example: Single Tree Dataset

Header from `tronko-build/example_datasets/single_tree/reference_tree.txt`:
```
1          <- 1 tree
13         <- max node name length
28         <- max taxonomy name length
94         <- max taxonomy line length
316	0	1466   <- 316 bases, root=0, 1466 species
```

This means:
- 1466 taxonomy lines follow the header
- Then 2931 nodes (2 * 1466 - 1)
- Each node has 1 header line + 316 posterior probability lines
- Total: ~930,000 lines

### Node Type Detection

**Leaf nodes**: `up[0] == -1 && up[1] == -1`
- Have accession name in header line
- Represent actual sequences from the MSA

**Internal nodes**: `up[0] != -1 && up[1] != -1`
- No name field in header line
- Posterior probabilities are computed from children

**Root node**: `down == -1`
- Index stored in `rootArr[treeIndex]`

## Code References

| Purpose | File | Lines |
|---------|------|-------|
| Write reference tree | `tronko-build/printtree.c` | 31-83 |
| Read reference tree | `tronko-assign/readreference.c` | 377-556 |
| Node struct definition | `tronko-assign/global.h` | 60-68 |
| type_of_PP definition | `tronko-assign/global.h` | 14-21 |
| Tree memory allocation | `tronko-assign/allocatetreememory.c` | - |

## Architecture Insights

1. **Multi-tree support**: The format supports multiple phylogenetic trees (partitions), allowing tronko to handle different gene regions or taxonomic groups separately.

2. **Posterior probability storage**: The key innovation of tronko is storing posterior probabilities at ALL nodes (not just leaves), enabling LCA computation using fractional likelihoods.

3. **Memory optimization**: The `OPTIMIZE_MEMORY` flag allows switching from `double` to `float` for posterior probabilities, halving memory usage for large databases.

4. **Taxonomy indexing**: The two-level taxonomy index (`taxIndex[0]` for species, `taxIndex[1]` for rank) allows efficient lookup without string comparisons during assignment.

5. **File size**: For a database with N species and M base positions:
   - Nodes: 2N-1
   - Lines per node: M+1 (header + posteriors)
   - Total lines: ~(2N-1) * (M+1) + header + taxonomy

## Open Questions

1. How are the posterior probabilities computed during tronko-build? (See `likelihood.c`)
2. What determines the `taxIndex` values for internal nodes?
3. Is there compression support for the reference tree file?
