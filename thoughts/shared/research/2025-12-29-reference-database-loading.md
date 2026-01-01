---
date: 2025-12-29T00:00:00-08:00
researcher: Claude
git_commit: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
branch: main
repository: eDNA-Explorer/tronko
topic: "How tronko-assign loads reference database for assignment"
tags: [research, codebase, tronko-assign, reference-database, memory-loading]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
last_updated_note: "Added follow-up research on optimization opportunities"
---

# Research: Reference Database Loading in tronko-assign

**Date**: 2025-12-29
**Researcher**: Claude
**Git Commit**: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
**Branch**: main
**Repository**: eDNA-Explorer/tronko

## Research Question
How does tronko-assign load the reference database for assignment? Why does it need to load the entire database into memory? What files does it look for or need?

## Summary

tronko-assign loads the reference database through a multi-step process that reads a gzip-compressed file (`reference_tree.txt`) containing serialized phylogenetic trees with precomputed posterior nucleotide probabilities. **The entire database must be loaded into memory because the assignment algorithm requires random access to any node's likelihood values during phylogenetic placement** - the algorithm traverses the tree from matched leaf nodes up to the root, scoring each node along the path and computing LCA (Lowest Common Ancestor) across multiple tree partitions.

### Required Files
1. **Reference database** (`-f`): `reference_tree.txt` - gzip or plain text
2. **Reference FASTA** (`-a`): FASTA file for BWA indexing
3. **Query reads** (`-g` for single-end, `-1`/`-2` for paired-end): FASTA or FASTQ

## Detailed Findings

### 1. Database Loading Entry Point

The loading process begins in `main()` at `tronko-assign/tronko-assign.c:767-792`:

```c
// Validate file exists
struct stat st = {0};
if ( stat(opt.reference_file, &st) == -1 ){
    printf("Cannot find reference_tree.txt file. Exiting...\n");
    exit(-1);
}

// Open gzipped reference file
gzFile referenceTree = Z_NULL;
referenceTree = gzopen(opt.reference_file,"r");
assert(Z_NULL!=referenceTree);

// Read and parse the database
int* name_specs = (int*)malloc(3*sizeof(int));
numberOfTrees = readReferenceTree(referenceTree,name_specs);
gzclose(referenceTree);
```

### 2. Reference File Format

The `reference_tree.txt` file has a specific binary-like text format parsed by `readReferenceTree()` at `readreference.c:311-433`:

**Header Section:**
```
Line 1:  <numberOfTrees>          (e.g., 1)
Line 2:  <max_nodename>           (e.g., 13 - max accession name length)
Line 3:  <max_tax_name>           (e.g., 28 - max taxonomy name length)
Line 4:  <max_lineTaxonomy>       (e.g., 94 - max taxonomy line length)
```

**Per-Tree Metadata:**
```
<numbase>	<root>	<numspec>
```
Where:
- `numbase` = number of MSA alignment positions (e.g., 316)
- `root` = root node index (e.g., 0)
- `numspec` = number of species/leaf nodes (e.g., 1466)

**Taxonomy Section** (semicolon-delimited, 7 levels per species):
```
Larus heuglini;Larus;Laridae;Charadriiformes;Aves;Chordata;Eukaryota
```

**Node Data Section** (for each of `2*numspec-1` nodes):
```
<tree> <node> <up0> <up1> <down> <depth> <taxIdx0> <taxIdx1> [name]
```
Followed by `numbase` lines of 4 posterior probabilities (A, C, G, T).

### 3. Data Structures in Memory

The reference database populates these global structures defined in `global.h`:

#### Node Structure (`global.h:43-51`)
```c
typedef struct node{
    int up[2];           // Child node indices (-1 for leaves)
    int down;            // Parent node index
    int depth;           // Depth in tree (for LCA)
    double **posteriornc; // [numbase][4] likelihood values
    char *name;          // Leaf accession name (leaves only)
    int taxIndex[2];     // Indices into taxonomyArr
}node;
```

#### Global Arrays
| Variable | Type | Purpose |
|----------|------|---------|
| `treeArr` | `node **` | Array of phylogenetic trees |
| `taxonomyArr` | `char ****` | 4D taxonomy: `[tree][species][level][char]` |
| `numspecArr` | `int *` | Species count per tree |
| `numbaseArr` | `int *` | MSA length per tree |
| `rootArr` | `int *` | Root node index per tree |

### 4. Why Full Memory Loading is Required

The entire database must be in memory for three key reasons:

**A. Random Access During Tree Traversal**

The assignment algorithm (`assignment.c:143-210`) computes scores by traversing from leaf nodes to root:

```c
// Score computation requires random node access
if (locQuery[i]=='A'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][0];
}
```

Each query can match any leaf node, requiring immediate access to that node's likelihood values and its entire path to the root.

**B. LCA (Lowest Common Ancestor) Computation**

The core innovation of tronko is computing LCA using fractional likelihoods from *all* nodes, not just leaves. This requires:
- Accessing parent nodes via `treeArr[tree][node].down`
- Comparing depth values via `treeArr[tree][node].depth`
- Computing scores at every ancestral node

**C. Multi-Tree Voting**

When partitioned databases have multiple trees, the algorithm must:
- Score the same query against multiple trees simultaneously
- Compare results across trees for consensus
- Access any node in any tree at any time

### 5. Memory Allocation Pattern

The allocation happens in two phases:

**Phase 1: Allocate structural arrays** (`readreference.c:358-394`)
```c
numbaseArr = (int*)malloc(numberOfTrees*(sizeof(int)));
rootArr = (int*)malloc(numberOfTrees*(sizeof(int)));
numspecArr = (int*)malloc(numberOfTrees*(sizeof(int)));

treeArr = malloc(numberOfTrees*sizeof(struct node *));
for (i=0; i<numberOfTrees; i++){
    allocateTreeArrMemory(i, max_nodename);
}
```

**Phase 2: Per-node allocation** (`allocatetreememory.c:17-39`)
```c
// Binary tree: 2*numspec-1 nodes per tree
treeArr[partition] = malloc((numspecArr[partition]*2-1)*(sizeof(struct node)));

// Allocate posterior probability matrix per node
for (i=0; i<(numspecArr[partition]*2-1); i++){
    treeArr[partition][i].posteriornc = (double**)malloc(numbaseArr[partition]*sizeof(double*));
    for (k=0; k<numbaseArr[partition]; k++){
        treeArr[partition][i].posteriornc[k] = (double*)malloc(4*(sizeof(double)));
    }
}
```

### 6. Post-Loading Transformation

After loading, likelihoods are log-transformed (`tronko-assign.c:43-63`):

```c
void store_PPs_Arr(int numberOfRoots, double c){
    // Transform raw posteriors to log-likelihoods with scoring constant
    double f = d * treeArr[i][j].posteriornc[k][l];
    double g = e * (1-treeArr[i][j].posteriornc[k][l]);
    treeArr[i][j].posteriornc[k][l] = log( (f + g) );
}
```

This enables efficient score computation by summing log-probabilities.

## Code References

- `tronko-assign/tronko-assign.c:767-792` - Main loading entry point
- `tronko-assign/readreference.c:311-433` - `readReferenceTree()` parsing function
- `tronko-assign/allocatetreememory.c:17-39` - Node memory allocation
- `tronko-assign/allocatetreememory.c:54-71` - Taxonomy array allocation
- `tronko-assign/global.h:43-51` - Node structure definition
- `tronko-assign/global.h:27,38,183` - Global array declarations
- `tronko-assign/options.c:126-130` - Reference file option parsing
- `tronko-assign/assignment.c:143-210` - Score computation using likelihoods
- `tronko-assign/tronko-assign.c:43-63` - Log transformation of posteriors

## Architecture Insights

1. **Global State Pattern**: All tree data is stored in global variables for easy access across functions, typical of C scientific computing code.

2. **Lazy Name Allocation**: Only leaf nodes get `name` strings allocated, saving memory for internal nodes.

3. **Gzip Transparency**: Uses zlib's `gzopen/gzgets/gzclose` for transparent gzip handling.

4. **Two-Phase Parsing**: Header is read first to determine array sizes, then memory is allocated before reading node data - avoiding reallocation.

5. **Fixed Taxonomy Levels**: Hard-coded 7 taxonomic levels (domain through species) in the taxonomy array structure.

## Open Questions

1. **Memory Optimization**: Could streaming or memory-mapped files reduce peak memory usage for very large databases?

2. **Partial Loading**: Could a hierarchical index allow loading only relevant tree partitions?

3. **Parallelization**: The global state pattern limits multi-threaded loading; could this be refactored?

---

## Follow-up Research: Optimization Opportunities

### Memory Usage Analysis

For a typical database (1 tree, 1466 species, 316 base positions):

| Component | Memory | Percentage |
|-----------|--------|------------|
| **posteriornc inner arrays** | 28.26 MB | 72.4% |
| **posteriornc pointer arrays** | 7.07 MB | 18.1% |
| **Taxonomy array** | 2.60 MB | 6.7% |
| **Node structures** | 137 KB | 0.3% |
| **Leaf names** | 43 KB | 0.1% |
| **TOTAL** | ~39.0 MB | 100% |

**Memory Formula (simplified):**
```
Total ≈ T × (2S-1) × B × 40 bytes
```
Where T=trees, S=species, B=base positions.

The `posteriornc` arrays dominate memory usage at **90.5%**.

### Access Pattern Analysis

Critical finding: **Every query scores ALL nodes in matched tree partitions**, not just the path from leaf to root.

From `assignScores_Arr_paired()` at `assignment.c:24-65`:
- Pre-order traversal visits all `2*numspec-1` nodes
- Each node's `posteriornc[position][nucleotide]` is accessed once per aligned position
- Cannot skip nodes because LCA requires comparing scores across ALL nodes within the confidence interval

This means partial loading or lazy loading is **not feasible** without changing the algorithm.

### Optimization Opportunities

#### 1. **Reduce Precision: double → float (HIGH IMPACT)**

**Current state:**
- Posterior probabilities stored as `double` (8 bytes each)
- After log-transform, values range from -5.7 to 0.0
- Scores range from -1.5 (good) to -855 (poor) for 150bp reads

**Precision requirements:**
- LCA comparison uses `Cinterval` threshold (typically 1.0-10.0)
- Only needs ~0.01 relative precision for correct LCA computation
- `float` (23-bit mantissa) provides ~10^-6 relative error - **more than sufficient**

**Savings:**
```
Current: (2S-1) × B × 4 × 8 bytes = 28.26 MB
With float: (2S-1) × B × 4 × 4 bytes = 14.13 MB
Savings: ~50% reduction in posteriornc storage (14 MB for example database)
```

**Implementation:**
- Change `double **posteriornc` to `float **posteriornc` in `global.h:48`
- Change `type_of_PP` from `double` to `float` in `global.h:14`
- Update `sscanf("%lf", ...)` to `sscanf("%f", ...)` in `readreference.c`
- Update `getscore_Arr()` return type and accumulator

**Risk:** Low - precision analysis confirms float is adequate for LCA computation.

#### 2. **Bulk Memory Allocation (MEDIUM IMPACT)**

**Current state:**
Per node, allocates `1 + numbase` separate malloc calls:
```c
treeArr[i].posteriornc = malloc(numbase * sizeof(double*));  // 1 call
for (k=0; k<numbase; k++){
    treeArr[i].posteriornc[k] = malloc(4 * sizeof(double));  // numbase calls
}
```

For 2931 nodes × 317 calls = **~928,000 malloc calls** per tree.

**Optimization:**
Allocate one contiguous block per tree:
```c
// Single allocation for entire tree's posterior data
size_t total = (2*numspec-1) * numbase * 4 * sizeof(float);
float *block = malloc(total);

// Set up pointer arrays to reference into block
for (node=0; node < 2*numspec-1; node++) {
    treeArr[tree][node].posteriornc = block + (node * numbase * 4);
}
```

**Savings:**
- Reduce malloc calls from ~928,000 to ~1 per tree
- Better cache locality - contiguous memory access
- Faster loading (fewer system calls)
- Estimated 2-5x faster allocation time

#### 3. **Binary File Format (MEDIUM IMPACT)**

**Current state:**
- Text format with `sscanf()` per value
- For 2931 nodes × 316 positions × 4 nucleotides = **~3.7 million sscanf calls**

**Optimization:**
Create binary format for `reference_tree.txt.bin`:
```c
// Binary header
struct header { int trees, max_nodename, max_taxname, max_linetax; };

// Per-tree binary block
struct tree_meta { int numbase, root, numspec; };
// Followed by: float[nodes][numbase][4] as contiguous block
```

**Loading:**
```c
float *block = mmap(fd, total_size, PROT_READ, MAP_PRIVATE, 0, offset);
// Or: fread(block, sizeof(float), total_count, fp);
```

**Savings:**
- Eliminate all sscanf/strtok parsing overhead
- Potential 10-50x faster loading
- Could use memory-mapping for zero-copy loading

**Trade-off:** Requires tronko-build to output binary format (or conversion tool).

#### 4. **Eliminate Outer Pointer Arrays (MEDIUM IMPACT)**

**Current state:**
```c
double **posteriornc;  // [numbase][4]
// Requires: numbase pointers + numbase*4 doubles
```

**Optimization:**
Use 1D array with stride:
```c
float *posteriornc;  // [numbase * 4]
// Access: posteriornc[position * 4 + nucleotide]
```

**Savings:**
- Eliminate 7.07 MB of pointer arrays (18% of total)
- Better cache performance (no pointer chasing)
- Combined with float: 39 MB → ~14 MB total

#### 5. **Deferred Log Transform (LOW IMPACT)**

**Current state:**
`store_PPs_Arr()` transforms all values after loading - touches entire dataset.

**Optimization:**
Skip transformation during load; apply during first access:
```c
// Store raw, transform on-demand (or pre-transform in binary file)
```

**Trade-off:** Either adds per-access overhead OR requires modifying file format.

### Recommendations Summary

| Optimization | Memory Savings | Speed Improvement | Implementation Effort |
|--------------|----------------|-------------------|----------------------|
| float precision | ~50% | Minor | Low (type changes) |
| Bulk allocation | ~18% (indirect) | 2-5x load | Medium |
| Binary format | None | 10-50x load | High (both tools) |
| 1D array layout | ~18% | 10-20% score | Medium |
| Combined (all) | ~65% | 10-50x load | High |

### What Cannot Be Optimized (Algorithm Constraints)

1. **Cannot use partial loading:** Every node must be scored per query
2. **Cannot stream on-demand:** Random access pattern across all nodes
3. **Cannot skip internal nodes:** LCA requires scoring all nodes, not just leaf path
4. **Cannot reduce to path-only:** Core algorithm innovation is fractional LCA at ALL nodes

### Recommended Implementation Order

1. **Phase 1 (Quick win):** Change `double` → `float` (~50% memory, low risk)
2. **Phase 2 (Medium effort):** Bulk allocation + 1D array layout (faster loading)
3. **Phase 3 (Future):** Binary file format (requires tronko-build changes)
