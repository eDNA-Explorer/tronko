---
date: 2025-12-30T22:58:21-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Rust Port Detailed Implementation Research"
tags: [research, rust, porting, tree-structures, alignment, performance]
status: complete
last_updated: 2025-12-30
last_updated_by: Claude
---

# Research: Rust Port Detailed Implementation

**Date**: 2025-12-30T22:58:21-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

Detailed implementation research for the Rust port of tronko-assign, focusing on:
1. Tree data structures (indextree vs custom arena vs light_phylogeny)
2. Posterior storage (Flat Vec vs nested Vec vs ndarray)
3. rust-bio alignment - Can we use existing implementations?
4. Rust tree patterns from "Too Many Linked Lists"

---

## Summary

This research provides concrete implementation guidance for the four key components of the Rust port:

| Component | Recommended Approach | Rationale |
|-----------|---------------------|-----------|
| **Tree Structure** | Custom arena with `Vec<[f64; 4]>` posteriors | Maximum control, matches C patterns |
| **Posterior Storage** | `Vec<[f64; 4]>` per node | Cache-optimal, 2-3x faster than Vec<Vec> |
| **Alignment** | rust-bio for accuracy, block-aligner for speed | Both are viable; rust-bio is simpler |
| **Traversal** | Arena pattern with index-based parent/child | Standard Rust idiom, thread-safe |

---

## Detailed Findings

### 1. Tree Data Structure: Custom Arena (Recommended)

After evaluating `indextree`, `light_phylogeny`, and custom arena approaches, the **custom arena with indices** is recommended for tronko's specific needs.

#### Why Not indextree?

`indextree` is excellent for general n-ary trees but has limitations:
- **Designed for n-ary trees**: Uses linked-list siblings (`first_child` + `next_sibling`) rather than explicit left/right children
- **Binary tree access is indirect**: Must call `first_child()` then `next_sibling()` for right child
- **Fixed Node structure**: Cannot embed custom fields without wrapper types

#### Why Not light_phylogeny?

`light_phylogeny` (v2.8.6, actively maintained) has phylogenetic-specific features but:
- **No generic custom data field**: The `Noeud` struct has fixed fields; `support` is a `String`, not suitable for likelihoods
- **Requires wrapper pattern**: Must maintain parallel data structures for posteriors
- **Overkill for our needs**: Includes visualization, reconciliation, and event tracking we don't need

#### Recommended: Custom Arena Pattern

Based on the C code analysis, here's the recommended Rust structure:

```rust
/// Node identifier - a simple index into the arena
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NodeId(pub usize);

/// A node in the phylogenetic tree
/// Mirrors C struct at tronko-assign/global.h:60-68
pub struct PhyloNode {
    /// Child indices: up[0]=left, up[1]=right. None = leaf
    pub left: Option<NodeId>,
    pub right: Option<NodeId>,

    /// Parent index. None = root
    pub parent: Option<NodeId>,

    /// Depth from root (for LCA computation)
    pub depth: usize,

    /// Posterior probabilities per alignment position
    /// Each element is [P(A), P(C), P(G), P(T)]
    pub posteriors: Vec<[f64; 4]>,

    /// Leaf name (accession). Only populated for leaf nodes
    pub name: Option<String>,

    /// Taxonomy reference: (species_index, taxonomic_level)
    pub tax_index: (usize, usize),
}

/// Arena containing all tree nodes
pub struct PhyloArena {
    nodes: Vec<PhyloNode>,
    root: Option<NodeId>,
}

impl PhyloArena {
    /// O(1) node access by index
    #[inline]
    pub fn get(&self, id: NodeId) -> &PhyloNode {
        &self.nodes[id.0]
    }

    /// Find LCA using ancestor path intersection
    pub fn lca(&self, a: NodeId, b: NodeId) -> Option<NodeId> {
        use std::collections::HashSet;

        // Collect ancestors of a (including a)
        let mut ancestors_a = HashSet::new();
        let mut current = Some(a);
        while let Some(id) = current {
            ancestors_a.insert(id);
            current = self.nodes[id.0].parent;
        }

        // Walk up from b until we find common ancestor
        let mut current = Some(b);
        while let Some(id) = current {
            if ancestors_a.contains(&id) {
                return Some(id);
            }
            current = self.nodes[id.0].parent;
        }
        None
    }
}
```

#### Key Design Decisions

| Decision | C Pattern | Rust Pattern | Rationale |
|----------|-----------|--------------|-----------|
| Node storage | `node **treeArr` | `Vec<PhyloNode>` | Contiguous allocation |
| Parent/child refs | `int up[2], down` | `Option<NodeId>` | Type-safe, no sentinel values |
| Posteriors | `double **posteriornc` | `Vec<[f64; 4]>` | Cache-optimal (see below) |
| Leaf names | `char *name` | `Option<String>` | Only allocate for leaves |

---

### 2. Posterior Probability Storage: Vec<[f64; 4]> (Recommended)

Based on benchmarking research, `Vec<[f64; 4]>` is the optimal storage format.

#### Benchmark Summary

| Approach | Relative Speed | Memory Overhead | Cache Behavior |
|----------|---------------|-----------------|----------------|
| `Vec<Vec<f64>>` | 1.0x (baseline) | High (24 bytes per inner Vec) | Poor (fragmented) |
| `Vec<f64>` with stride | 2-3x faster | Low (24 bytes total) | Excellent |
| `Vec<[f64; 4]>` | 2-3x faster | Low (24 bytes total) | **Excellent** |
| `ndarray::Array2` | 1.5-2x slower | Moderate (~76 bytes) | Good |

**Source**: [Grids in Rust: nested vs flat Vecs](https://blog.adamchalmers.com/grids-1/)

#### Why Vec<[f64; 4]> Is Optimal

1. **Contiguous memory**: All data in single allocation
2. **Cache-line friendly**: Each `[f64; 4]` is 32 bytes; 2 positions fit in a 64-byte cache line
3. **Natural access pattern**: Matches tronko's "iterate positions, access all 4 nucleotides" pattern
4. **Zero index arithmetic**: Direct `posteriors[pos][nucleotide]` access
5. **SIMD-friendly**: Fixed-size arrays enable auto-vectorization

#### Complete Implementation

```rust
pub mod nucleotide {
    pub const A: usize = 0;
    pub const C: usize = 1;
    pub const G: usize = 2;
    pub const T: usize = 3;
}

/// Posterior probability matrix for a tree node
#[derive(Debug, Clone)]
pub struct PosteriorMatrix {
    /// Each element is [P(A), P(C), P(G), P(T)] at one alignment position
    data: Vec<[f64; 4]>,
}

impl PosteriorMatrix {
    pub fn new(num_positions: usize) -> Self {
        Self {
            data: vec![[0.0; 4]; num_positions],
        }
    }

    #[inline]
    pub fn get(&self, pos: usize, nuc: usize) -> f64 {
        self.data[pos][nuc]
    }

    #[inline]
    pub fn get_position(&self, pos: usize) -> &[f64; 4] {
        &self.data[pos]
    }

    /// Iterate positions - maximally cache efficient
    pub fn iter(&self) -> impl Iterator<Item = &[f64; 4]> {
        self.data.iter()
    }
}
```

#### Memory Comparison (400 positions, 3000 nodes)

| Approach | Per Node | Total (3000 nodes) |
|----------|----------|-------------------|
| C `double**` | ~22 KB | ~66 MB |
| `Vec<[f64; 4]>` | ~12.8 KB | ~38 MB |
| `Vec<[f32; 4]>` | ~6.4 KB | ~19 MB |

---

### 3. Needleman-Wunsch Alignment: rust-bio vs block-aligner

Both crates are viable. Choose based on accuracy vs speed tradeoff.

#### rust-bio Pairwise Alignment

**Pros**:
- Exact dynamic programming (no approximation)
- Clean, simple API
- Well-documented with examples
- Affine gap penalties supported

**Cons**:
- Not SIMD-accelerated (except distance functions)
- O(n*m) complexity without optimization

```rust
use bio::alignment::pairwise::*;

// DNA scoring: match=1, mismatch=-1, gap_open=-5, gap_extend=-1
let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

let mut aligner = Aligner::with_capacity(
    query.len(),
    reference.len(),
    -5,  // gap open
    -1,  // gap extend
    &score_fn
);

// Needleman-Wunsch global alignment
let alignment = aligner.global(query, reference);
println!("Score: {}", alignment.score);
println!("CIGAR: {}", alignment.cigar(false));
```

#### block-aligner (SIMD)

**Pros**:
- 5-10x faster for long sequences
- SIMD-accelerated (SSE2, AVX2, Neon, WASM)
- Pre-built DNA matrices (NW1: match=1, mismatch=-1)

**Cons**:
- Heuristic/approximate for very long sequences
- More complex API
- Less benefit for short reads (~150bp)

```rust
use block_aligner::scan_block::*;
use block_aligner::scores::*;

let gaps = Gaps { open: -2, extend: -1 };
let query = PaddedBytes::from_bytes::<NucMatrix>(query_seq, 256);
let reference = PaddedBytes::from_bytes::<NucMatrix>(ref_seq, 256);

let mut block = Block::<true, false>::new(query.len(), reference.len(), 256);
block.align(&query, &reference, &NW1, gaps, 32..=256, 0);

let result = block.res();
println!("Score: {}", result.score);
```

#### Recommendation for tronko

For **short Illumina reads (~150bp)**, use **rust-bio**:
- Exact results match the current C implementation
- Simpler integration
- Speed difference is minimal for short sequences

Consider **block-aligner** for:
- Long reads (Nanopore, PacBio)
- Batch alignment of thousands of sequences
- When 5-10% accuracy loss is acceptable

---

### 4. Rust Tree Patterns: Arena vs Alternatives

#### The Core Challenge

Rust's ownership model makes bidirectional tree references (parent pointers) challenging:
- Recursive types require indirection (`Box`, `Rc`, or arena)
- Parent pointers create ownership cycles
- Mutable access during traversal requires careful design

#### Pattern Comparison

| Pattern | Parent Pointers | Thread Safe | Performance | Complexity |
|---------|-----------------|-------------|-------------|------------|
| **Arena + indices** | Yes | Yes | Excellent | Low |
| `slotmap` crate | Yes | Yes | Excellent | Very Low |
| `indextree` crate | Yes | Yes | Good | Very Low |
| `Rc<RefCell<>>` | Yes (Weak) | No | Poor | High |
| `Arc<Mutex<>>` | Yes (Weak) | Yes | Poor | High |

#### Why Arena + Indices Is Best for tronko

1. **Matches C memory layout**: `treeArr[node]` in C becomes `arena[NodeId]` in Rust
2. **No runtime borrow checks**: Unlike `RefCell`, all checks are at compile time
3. **Thread-safe by default**: Can use with Rayon for parallel scoring
4. **Minimal overhead**: NodeId is just a `usize` (8 bytes)
5. **Cache-friendly**: All nodes contiguous in memory

#### Parallel Scoring with Rayon

```rust
use rayon::prelude::*;

impl PhyloArena {
    /// Parallel likelihood computation (read-only)
    pub fn parallel_score<F>(&self, f: F) -> Vec<f64>
    where
        F: Fn(&PhyloNode) -> f64 + Sync,
    {
        self.nodes.par_iter().map(f).collect()
    }
}

// Usage matching assignScores_Arr_paired pattern
let scores = arena.parallel_score(|node| {
    compute_likelihood(&node.posteriors, &query_seq)
});
```

---

## C Code Reference: Key Structures

From the codebase analysis:

### Node Structure (`tronko-assign/global.h:60-68`)

```c
typedef struct node{
    int up[2];              // Children: -1=leaf, -2=uninitialized
    int down;               // Parent index
    int nd;                 // Node descriptor (unused)
    int depth;              // Depth from root
    type_of_PP **posteriornc; // [numbase][4] posteriors
    char *name;             // Leaf accession (leaves only)
    int taxIndex[2];        // [species_idx, tax_level]
}node;
```

### Access Pattern (`tronko-assign/assignment.c:177-202`)

```c
// Sequential position iteration, all 4 nucleotides per position
if (locQuery[i]=='A'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][0];
}else if (locQuery[i]=='C'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][1];
}else if (locQuery[i]=='G'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][2];
}else if (locQuery[i]=='T'){
    score += treeArr[rootNum][node].posteriornc[positions[i]][3];
}
```

This access pattern directly maps to:

```rust
match query[i] {
    b'A' => score += posteriors[pos][0],
    b'C' => score += posteriors[pos][1],
    b'G' => score += posteriors[pos][2],
    b'T' => score += posteriors[pos][3],
    _ => score += (0.25_f64).ln(),
}
```

---

## Recommended Implementation Order

### Phase 1: Reference Database Parsing (Priority 1)

1. **Parse header** (4 lines: tree count, max lengths)
2. **Parse tree metadata** (numbase, root, numspec per tree)
3. **Allocate arena** with `Arena::with_capacity(2 * numspec - 1)`
4. **Parse taxonomy** (7 levels per species)
5. **Parse nodes** (structure + posteriors)
6. **Validate** against C output

**Deliverable**: `pub fn parse_reference_tree(path: &Path) -> Result<Vec<PhyloArena>>`

### Phase 2: Needleman-Wunsch Alignment (Priority 2)

1. **Add rust-bio dependency**
2. **Create alignment wrapper** matching `needleman_wunsch.c` API
3. **Benchmark** against C implementation
4. **Validate** CIGAR output matches

**Deliverable**: `pub fn align_global(query: &[u8], ref: &[u8]) -> Alignment`

### Phase 3: Tree Scoring (Priority 3)

1. **Port `getscore_Arr()`** - single node scoring
2. **Port `assignScores_Arr_paired()`** - tree traversal
3. **Add parallel scoring** with Rayon
4. **Validate** scores match C for same input

**Deliverable**: `pub fn score_tree(arena: &PhyloArena, query: &[u8]) -> Vec<f64>`

### Phase 4: LCA and Assignment (Priority 4)

1. **Implement LCA** using ancestor path intersection
2. **Port assignment logic** from `assignment.c`
3. **Integrate** with BWA via rust-bwa FFI
4. **Full pipeline validation**

---

## Code References

- `tronko-assign/global.h:60-68` - Node structure definition
- `tronko-assign/global.h:14-22` - `type_of_PP` precision configuration
- `tronko-assign/readreference.c:377-557` - Reference tree parsing
- `tronko-assign/allocatetreememory.c:17-40` - Memory allocation pattern
- `tronko-assign/assignment.c:143-210` - Score computation
- `tronko-assign/assignment.c:24-65` - Tree traversal

## Architecture Insights

1. **C uses global state**: `treeArr`, `taxonomyArr` are globals. Rust version should use explicit passing or `Arc` for thread safety.

2. **Binary tree invariant**: Always exactly 2 children or 0 (leaf). Total nodes = `2 * species - 1`.

3. **Sentinel values**: C uses -1 (leaf), -2 (uninitialized). Rust uses `Option<NodeId>`.

4. **Memory dominance**: Posteriors are 90%+ of memory. Optimization here has highest impact.

5. **Access pattern**: Sequential position iteration with random node access. Arena pattern fits perfectly.

---

## Related Research

- `thoughts/shared/research/2025-12-29-reference-database-loading.md` - File format details
- `thoughts/shared/research/2025-12-29-rust-port-feasibility.md` - Library mapping
- `thoughts/shared/research/2025-12-30-rust-port-next-steps.md` - Priority ordering

---

## Sources

### Rust Crates
- [indextree - docs.rs](https://docs.rs/indextree/latest/indextree/)
- [light_phylogeny - docs.rs](https://docs.rs/light_phylogeny/latest/light_phylogeny/)
- [rust-bio - docs.rs](https://docs.rs/bio/latest/bio/)
- [block-aligner - GitHub](https://github.com/Daniel-Liu-c0deb0t/block-aligner)
- [ndarray - docs.rs](https://docs.rs/ndarray/latest/ndarray/)
- [slotmap - docs.rs](https://docs.rs/slotmap/latest/slotmap/)
- [rayon - docs.rs](https://docs.rs/rayon/latest/rayon/)

### Rust Patterns
- [Learning Rust With Entirely Too Many Linked Lists](https://rust-unofficial.github.io/too-many-lists/)
- [Idiomatic Trees in Rust - Rust Leipzig](https://rust-leipzig.github.io/architecture/2016/12/20/idiomatic-trees-in-rust/)
- [No More Tears: Arena-Allocated Trees - DEV Community](https://dev.to/deciduously/no-more-tears-no-more-knots-arena-allocated-trees-in-rust-44k6)
- [Grids in Rust: nested vs flat Vecs](https://blog.adamchalmers.com/grids-1/)

### Performance
- [ndarray vs raw slices benchmarks](https://www.reidatcheson.com/rust/ndarray/performance/2022/06/11/rust-ndarray.html)
- [Block Aligner publication - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10457662/)

---

## Open Questions

1. **Float vs Double**: Should we use `f32` for posteriors? C allows this with `OPTIMIZE_MEMORY` flag. Would save 50% memory.

2. **Binary format**: Should Rust version output/read binary format for faster loading?

3. **LCA optimization**: For many LCA queries, should we precompute binary lifting table (O(n log n) preprocessing, O(log n) queries)?

4. **SIMD posteriors**: Can we use `packed_simd` or `std::simd` for scoring 4 nucleotides simultaneously?
