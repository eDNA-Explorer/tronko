---
date: 2025-12-29T00:00:00-08:00
researcher: Claude
git_commit: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
branch: main
repository: tronko
topic: "tronko-assign Data Flow Overview"
tags: [research, codebase, tronko-assign, bioinformatics, phylogenetic-assignment]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
last_updated_note: "Added memory optimization analysis"
---

# Research: tronko-assign Data Flow Overview

**Date**: 2025-12-29
**Researcher**: Claude
**Git Commit**: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
**Branch**: main
**Repository**: tronko

## Research Question
What is the overall data flow of tronko-assign? Outline the steps at a high level.

## Summary

tronko-assign is a taxonomic assignment tool that assigns query DNA sequences to taxonomic classifications using phylogenetic trees with posterior probability profiles. The program takes query sequences (FASTA/FASTQ), aligns them against reference sequences in phylogenetic trees, scores each node based on posterior nucleotide probabilities, and outputs taxonomic assignments based on the Lowest Common Ancestor (LCA) of the best-scoring nodes.

## High-Level Data Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                          TRONKO-ASSIGN DATA FLOW                        │
└─────────────────────────────────────────────────────────────────────────┘

     ┌──────────────────┐      ┌──────────────────┐
     │  Reference DB    │      │  Query Sequences │
     │  (gzipped)       │      │  (FASTA/FASTQ)   │
     └────────┬─────────┘      └────────┬─────────┘
              │                         │
              ▼                         ▼
     ┌──────────────────┐      ┌──────────────────┐
     │ 1. Parse Trees   │      │ 2. Pre-scan      │
     │    & Taxonomy    │      │    Query Files   │
     └────────┬─────────┘      └────────┬─────────┘
              │                         │
              ▼                         │
     ┌──────────────────┐              │
     │ 3. Transform     │              │
     │    Posteriors    │              │
     │    (log-prob)    │              │
     └────────┬─────────┘              │
              │                         │
              ▼                         ▼
     ┌──────────────────────────────────────────┐
     │ 4. BWA Index Build (if needed)           │
     │    - Creates FM-index for leaf sequences │
     └────────────────────┬─────────────────────┘
                          │
                          ▼
     ┌──────────────────────────────────────────┐
     │ 5. Read Query Batches                    │
     │    - Reads N lines at a time             │
     │    - Handles single/paired reads         │
     └────────────────────┬─────────────────────┘
                          │
                          ▼
     ┌──────────────────────────────────────────┐
     │ 6. Multi-threaded Processing Loop        │
     │    ┌────────────────────────────────┐    │
     │    │  a. BWA Alignment              │    │
     │    │     - Find candidate leaves    │    │
     │    ├────────────────────────────────┤    │
     │    │  b. WFA2/NW Alignment          │    │
     │    │     - Align to leaf sequences  │    │
     │    ├────────────────────────────────┤    │
     │    │  c. Score Tree Traversal       │    │
     │    │     - Score all nodes in tree  │    │
     │    ├────────────────────────────────┤    │
     │    │  d. Confidence Interval        │    │
     │    │     - Select nodes within CI   │    │
     │    ├────────────────────────────────┤    │
     │    │  e. LCA Computation            │    │
     │    │     - Find common ancestor     │    │
     │    └────────────────────────────────┘    │
     └────────────────────┬─────────────────────┘
                          │
                          ▼
     ┌──────────────────────────────────────────┐
     │ 7. Write Results                         │
     │    - TSV: readname, taxonomy, scores     │
     └──────────────────────────────────────────┘
```

## Detailed Steps

### Step 1: Initialization and Option Parsing
**Location**: `tronko-assign.c:744-788`

The program starts in `main()`, initializes default options, and parses command-line arguments:
- Parses required inputs: reference DB file (`-f`), reference FASTA (`-a`), output file (`-o`)
- Determines read mode: single (`-s`) vs paired (`-p`)
- Sets algorithm options: WFA2 (default) or Needleman-Wunsch (`-w`)
- Configures parallelism: number of cores (`-C`), batch size (`-L`)

### Step 2: Reference Tree Loading
**Location**: `readreference.c:311-433`

Loads the reference database file (can be gzipped):
1. Reads header metadata: number of trees, max name lengths
2. For each tree, reads: number of bases, root index, number of species
3. Parses taxonomy for each species (7 levels: species→kingdom)
4. Loads tree structure: node relationships (parent/children), taxonomy indices
5. Loads posterior probability matrices: 4 values (A,C,G,T) per position per node

### Step 3: Posterior Probability Transformation
**Location**: `tronko-assign.c:43-63`

Transforms raw posterior probabilities into log-likelihood scores:
```
P'(base) = log((1-c) × P(base) + (c/3) × (1-P(base)))
```
Where `c` is the score constant (default: 0.01). Missing data positions are marked with value 1.

### Step 4: BWA Index Construction
**Location**: `tronko-assign.c:872` (single) / `1079` (paired)

If not skipped (`-6` flag), builds BWA FM-index from reference FASTA for fast seed alignment.

### Step 5: Query File Pre-scanning
**Location**: `readreference.c:454-477`

Scans query files to determine:
- Maximum read name length
- Maximum sequence length

This enables proper memory allocation for batch processing.

### Step 6: Batch Processing Loop
**Location**: `tronko-assign.c:870-1173`

Reads queries in batches (default: 50,000 lines) and processes in parallel:

#### 6a. BWA Seed Alignment
**Location**: `bwa_source_files/` integration

Uses BWA-MEM to find candidate leaf nodes:
- Maps query sequences against reference FASTA
- Returns matching leaf coordinates (tree index, node index)
- Handles concordant vs discordant paired-end matches

#### 6b. Sequence Alignment (WFA2 or Needleman-Wunsch)
**Location**: `placement.c:70-125` (WFA2) / `needleman_wunsch.c:26-145` (NW)

For each BWA hit, performs detailed pairwise alignment:
- Extracts leaf sequence from tree (`getSequenceinRoot.c`)
- Aligns query against leaf sequence
- Uses affine gap penalties (WFA2: mismatch=4, gap_open=6, gap_extend=2)
- Produces CIGAR string and aligned sequences

#### 6c. Tree Traversal and Scoring
**Location**: `assignment.c:24-210`

Recursively scores all nodes in the tree:
```c
for each aligned position:
    score += log_posterior[position][observed_base]
```
- Traverses from root through all nodes (not just leaves)
- Accumulates log-likelihood scores based on posterior probabilities
- Handles gaps, missing data, and ambiguous bases

#### 6d. Best Node Selection with Confidence Interval
**Location**: `placement.c:862-930`

Selects assignment based on scores:
1. Finds maximum score across all matches and nodes
2. Selects all nodes within confidence interval (default: 5) of maximum
3. Votes for nodes meeting the threshold

#### 6e. LCA Computation
**Location**: `tronko-assign.c:91-134`

Finds Lowest Common Ancestor of selected nodes:
- Uses depth-based algorithm to find common ancestor
- Works within single tree or across multiple trees
- Returns the node representing the most specific shared taxonomy

### Step 7: Taxonomy Extraction and Output
**Location**: `tronko-assign.c:556-683`, `970-974`

Constructs final output:
1. Extracts taxonomy path from LCA node (species→kingdom)
2. Formats as semicolon-delimited string
3. Writes TSV row: `readname\ttaxonomy\tscore\tfwd_mismatch\trev_mismatch\ttree\tnode`

## Key Data Structures

### Node Structure (`global.h:43-51`)
```c
typedef struct node {
    int up[2];           // Children indices (-1 for leaves)
    int down;            // Parent index
    int depth;           // Tree depth
    double **posteriornc; // [position][nucleotide] posteriors
    char *name;          // Leaf accession name
    int taxIndex[2];     // [taxonomy_row, taxonomy_level]
} node;
```

### Results Structure (`global.h:73-91`)
```c
typedef struct resultsStruct {
    type_of_PP ***nodeScores;  // [match][tree][node] scores
    int **voteRoot;            // [tree][node] vote flags
    char **taxonPath;          // Formatted result strings
    type_of_PP *minimum;       // [score, fwd_mm, rev_mm]
    // ... alignment structures
} resultsStruct;
```

## Input/Output Summary

### Inputs
| File | Flag | Description |
|------|------|-------------|
| Reference DB | `-f` | Gzipped tree/posterior file from tronko-build |
| Reference FASTA | `-a` | FASTA for BWA indexing |
| Query reads | `-g` (single) or `-1`/`-2` (paired) | FASTA/FASTQ sequences |

### Primary Output
| Column | Description |
|--------|-------------|
| Readname | Query sequence identifier |
| Taxonomic_Path | Semicolon-delimited taxonomy (kingdom→species) |
| Score | Log-likelihood assignment score |
| Forward_Mismatch | Forward read mismatch count |
| Reverse_Mismatch | Reverse read mismatch count |
| Tree_Number | Index of assigned tree |
| Node_Number | Index of LCA node |

### Optional Debug Outputs
- `scores_all_nodes.txt` - All node scores (`-7` flag)
- `site_scores.txt` - Per-site scoring details (`-7` flag)
- Node info file - Tree/node reference mapping (`-5` flag)

## Code References
- `tronko-assign/tronko-assign.c:744` - Main entry point
- `tronko-assign/options.c:70-217` - Command-line parsing
- `tronko-assign/readreference.c:311-433` - Reference tree loading
- `tronko-assign/placement.c:70-125` - WFA2 alignment
- `tronko-assign/assignment.c:24-210` - Score computation
- `tronko-assign/global.h:43-51` - Node data structure

## Architecture Insights

1. **Batch Processing**: Queries are processed in batches to manage memory, with configurable batch size via `-L`

2. **Two-Stage Alignment**: BWA provides fast initial filtering, WFA2/NW provides accurate scoring

3. **Log-Space Scoring**: All probability calculations use log-space to avoid numerical underflow

4. **Confidence-Based Assignment**: Rather than picking a single best node, selects all nodes within a confidence interval and finds their LCA for robust taxonomy assignment

5. **Multi-Tree Support**: Can work with multiple phylogenetic trees simultaneously, useful for handling multiple genetic markers or partitioned datasets

## Open Questions

1. How are the reference trees built? (See tronko-build component)
2. What is the optimal confidence interval (`-c`) for different use cases?
3. How does performance scale with number of trees and tree size?

---

## Follow-up Research: Memory Optimization Analysis

### Memory Bottleneck Summary

The **posterior probability matrices** (`posteriornc`) are the dominant memory consumer, accounting for 80-95% of total memory usage.

**Memory Formula per Tree:**
```
Memory = (2 × numspec - 1) × numbase × 4 × 8 bytes
```

**Example Calculation:**
- 100 trees, 1,000 species each, 1,500 base positions
- Per tree: `(2×1000-1) × 1500 × 4 × 8 = 95.9 MB`
- Total for posteriors: `100 × 95.9 MB ≈ 9.6 GB`

### Memory Allocation Breakdown

| Component | Size Formula | Typical Size | Shared/Per-Thread |
|-----------|-------------|--------------|-------------------|
| **Posterior matrices** | `Σ(2×numspec[i]-1) × numbase[i] × 32` | **9+ GB** | Shared |
| Taxonomy array | `Σ(numspec[i]) × 7 × max_tax_name` | 100-500 MB | Shared |
| BWA FM-index | ~4× reference FASTA size | 100-500 MB | Shared |
| nodeScores | `10 × numTrees × Σ(2×numspec[i]-1) × 8` | 10-100 MB | Per-thread |
| voteRoot | `numTrees × Σ(2×numspec[i]-1) × 4` | 5-50 MB | Per-thread |
| Query batch | `batchSize × (max_query_len + max_name_len)` | 50-200 MB | Shared |

### Key Observations

1. **All posteriors loaded upfront**: Every node's posterior matrix is loaded into RAM before any queries are processed (`allocatetreememory.c:17-40`)

2. **Sparse access pattern**: For any given query:
   - BWA matches ~1-10 leaves out of thousands
   - Only nodes along paths from matched leaves to root are scored
   - Most loaded posteriors are never accessed

3. **Sequential file format**: The reference DB is a text file read line-by-line (`readreference.c:311-433`), making random access impossible without format changes

4. **No memory-mapped I/O**: File is read entirely into heap memory, no `mmap()` usage

### Optimization Opportunities

#### 1. **Lazy Loading with Memory-Mapped Binary Format** (High Impact)

**Current state**: All posteriors loaded into heap at startup
**Proposed**: Binary format with mmap for on-demand page loading

```
┌─────────────────────────────────────────────────────────┐
│ NEW BINARY FORMAT                                       │
├─────────────────────────────────────────────────────────┤
│ Header (fixed size):                                    │
│   - numberOfTrees, metadata offsets                     │
├─────────────────────────────────────────────────────────┤
│ Tree Index Table:                                       │
│   - Per tree: numspec, numbase, root, file offset       │
├─────────────────────────────────────────────────────────┤
│ Node Index Table:                                       │
│   - Per node: file offset to posterior data             │
├─────────────────────────────────────────────────────────┤
│ Taxonomy Section (load entirely - small)                │
├─────────────────────────────────────────────────────────┤
│ Posterior Data Section (mmap, access on-demand):        │
│   - Packed binary: node posteriors contiguous           │
└─────────────────────────────────────────────────────────┘
```

**Benefits**:
- OS pages in only accessed posteriors
- Automatic eviction under memory pressure
- Near-zero startup time for posteriors

**Implementation**:
1. New tool: `tronko-build --binary-format` to generate `.tdb` files
2. Replace `readReferenceTree()` with mmap-based loader
3. Keep `posteriornc` pointers, but point into mmap region

#### 2. **Tree-Level Streaming** (Medium Impact)

**Current state**: All trees loaded before processing
**Proposed**: Load trees on-demand based on BWA matches

```c
// Pseudocode for streaming approach
for each query_batch:
    bwa_matches = run_bwa(queries)
    needed_trees = unique(bwa_matches.tree_indices)

    for tree_id in needed_trees:
        if not loaded[tree_id]:
            load_tree(tree_id)  // Load just this tree's posteriors

    process_queries(query_batch, bwa_matches)

    // Optionally evict trees not needed recently
    evict_unused_trees(LRU_threshold)
```

**Benefits**:
- Only loads trees that have BWA hits
- Works with current text format (with index file)
- 10-50% memory reduction depending on hit distribution

**Requirements**:
- Pre-built index file with byte offsets per tree
- Seek support (incompatible with gzip streaming)

#### 3. **Node-Level Lazy Loading** (Highest Impact, Most Complex)

**Current state**: All `2×numspec-1` nodes per tree loaded
**Proposed**: Only load leaf posteriors initially; compute internal node posteriors on-demand

**Key insight**: Internal node posteriors are derived from leaf posteriors during `tronko-build`. They could potentially be recomputed from leaves.

```
Access Pattern Analysis:
─────────────────────────────────────────────────────
For a query matching leaf L in tree T:
  1. BWA finds leaf L
  2. Align query to leaf L sequence (needs L.posteriornc)
  3. Score traversal from root → L (needs all ancestor posteriors)

Nodes accessed: ~log2(numspec) ancestors + 1 leaf
Nodes loaded: 2×numspec-1 (ALL nodes)
Waste ratio: ~(2×numspec-1) / log2(numspec) ≈ 200× for 1000 species
─────────────────────────────────────────────────────
```

**Two approaches**:

**A. Store only leaves, recompute ancestors**:
- Memory: `numspec × numbase × 32` (50% reduction)
- Tradeoff: CPU cost to recompute ancestral posteriors per query
- May be prohibitive for deep trees

**B. Store all, but load on-demand via mmap** (recommended):
- Keep full file, use memory mapping
- OS loads only accessed pages
- Best of both worlds

#### 4. **Quantized Posteriors** (Medium Impact)

**Current state**: 8-byte doubles for each probability
**Proposed**: Use 2-byte fixed-point or 1-byte quantized values

```c
// Current: 32 bytes per position per node
double posteriornc[numbase][4];  // 4 × 8 bytes

// Option A: 16-bit fixed point (8 bytes per position)
int16_t posteriornc[numbase][4]; // 4 × 2 bytes
// Convert: value = (double)raw / 32768.0

// Option B: 8-bit log-quantized (4 bytes per position)
uint8_t posteriornc[numbase][4]; // 4 × 1 byte
// Convert: value = log_table[raw]  // Pre-computed log values
```

**Benefits**:
- 2-4× reduction in posterior memory
- May improve cache performance
- Minimal accuracy loss (posteriors already log-transformed)

**Tradeoffs**:
- Requires file format change
- Small precision loss (likely negligible for classification)

#### 5. **Reduce Per-Thread Memory** (Low Impact)

**Current state**: Each thread allocates full `nodeScores[10][numTrees][numNodes]`
**Proposed**: Sparse score storage

```c
// Current: Dense array for all nodes
type_of_PP ***nodeScores;  // [10][numTrees][2*numspec-1]

// Proposed: Hash map for scored nodes only
typedef struct {
    int tree_id;
    int node_id;
    double score;
} NodeScore;
HashMap<(tree,node), score> nodeScores;  // Only stores computed scores
```

**Benefits**:
- Per-thread memory proportional to nodes scored, not total nodes
- Typical savings: 10-100× per thread

### Recommended Implementation Order

1. **Binary format with mmap** (Biggest win, foundational change)
   - Enables all other optimizations
   - Reduces startup time from minutes to seconds
   - Memory usage becomes demand-paged

2. **Quantized posteriors** (Easy win with binary format)
   - 2-4× memory reduction
   - Implement alongside binary format

3. **Sparse nodeScores** (Per-thread improvement)
   - Reduces memory scaling with thread count
   - Independent of file format changes

4. **Tree-level LRU caching** (For very large databases)
   - Useful when total trees exceed RAM
   - Can be layered on top of mmap

### Compatibility Considerations

| Change | Backward Compatible? | Migration Path |
|--------|---------------------|----------------|
| Binary mmap format | No | New `-f` file, keep text parser for old files |
| Quantized posteriors | No | Version flag in binary header |
| Sparse nodeScores | Yes | Internal refactor only |
| Tree streaming | Partial | Requires uncompressed or indexed-gzip |

### Memory vs. I/O Tradeoff

```
                    CURRENT                    PROPOSED (mmap)
                    ───────                    ───────────────
Startup time:       O(file_size)              O(1)
                    (read entire file)         (just open + mmap)

Query memory:       O(total_posteriors)        O(accessed_posteriors)
                    (all in RAM)               (OS pages on demand)

I/O pattern:        Sequential read once       Random reads per query
                    (good for HDD)             (good for SSD, bad for HDD)

Peak memory:        Full DB size               Depends on access pattern
                    (predictable)              (may spike, but self-limits)
```

**Recommendation**: mmap approach is ideal for SSD-backed systems. For HDD or network storage, consider tree-level prefetching based on BWA batch results.
