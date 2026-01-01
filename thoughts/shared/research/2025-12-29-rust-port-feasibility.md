---
date: 2025-12-29T12:00:00-08:00
researcher: Claude
git_commit: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
branch: main
repository: tronko
topic: "Rust Port Feasibility Assessment"
tags: [research, rust, porting, bioinformatics, feasibility]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
---

# Research: Rust Port Feasibility Assessment

**Date**: 2025-12-29
**Researcher**: Claude
**Git Commit**: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
**Branch**: main
**Repository**: tronko

## Research Question

Based on the tronko-assign data flow analysis, what aspects of this project could be ported to Rust? Does Rust have existing libraries that do what is included in this project? What is the feasibility of porting this to Rust?

## Summary

**Porting tronko to Rust is feasible but will require significant effort.** The Rust bioinformatics ecosystem has matured considerably and offers high-quality replacements for most tronko functionality. The main challenges are:

1. **Tree data structures**: The hardest part - Rust's ownership model makes cyclic parent/child references challenging
2. **Embedded C libraries**: BWA and WFA2 have Rust bindings but no pure Rust replacements for all features
3. **External tools**: RAxML and FAMSA have no Rust equivalents; must continue as subprocess calls

**Key findings:**
- ~80% of tronko functionality has direct Rust library equivalents
- Performance should be comparable (within 5-10%) or better with SIMD optimizations
- Memory safety benefits are substantial for bioinformatics workloads
- Estimated effort: 3-6 months for experienced Rust developers

---

## Component-by-Component Analysis

### 1. FASTA/FASTQ Parsing

| C Implementation | Rust Replacement | Maturity |
|-----------------|------------------|----------|
| `readfasta.c` with zlib | **needletail** | Production-ready |

**needletail** (https://crates.io/crates/needletail):
- ~25x faster than Python equivalents
- Comparable speed to C `readfq` library
- Native gzip support
- Used in production at One Codex

**Alternative**: **seq_io** for zero-copy parsing with built-in parallel processing utilities.

### 2. Sequence Alignment

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| `needleman_wunsch.c` | **rust-bio::alignment::pairwise** | Full NW/SW implementations |
| WFA2 library | **libwfa** (FFI) or **rust_wfa** (pure Rust) | Both available |

**rust-bio pairwise alignment**:
- Global, semiglobal, and local alignment
- Banded alignment for long sequences
- Comparable speed to C++ SeqAn

**block-aligner** (https://crates.io/crates/block-aligner):
- SIMD-accelerated (SSE2, AVX2, ARM Neon)
- **5-10x faster** than traditional implementations
- Published in Bioinformatics 2023

### 3. BWA (Burrows-Wheeler Aligner)

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| `bwa_source_files/` embedded BWA | **rust-bwa** (FFI) | 10X Genomics maintains |

**rust-bwa** (https://github.com/10XGenomics/rust-bwa):
- Production-quality FFI bindings
- Returns `rust-htslib::Record` objects
- Supports paired-end alignment
- Same performance as C (uses C code via FFI)

**Alternative for FM-index only**: `rust-bio::data_structures::fmindex`

**Alternative for long reads**: **minimap2-rs** (successor to BWA for many use cases)

### 4. Phylogenetic Tree Handling

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| `getclade.c` Newick parsing | **phylotree** or **light_phylogeny** | Both support Newick |
| Tree traversal/LCA | **light_phylogeny** | Has LCA support |
| Node data structures | **indextree** or custom arena | See challenges below |

**phylotree** (https://crates.io/crates/phylotree):
- Newick parsing and writing
- Tree traversal (preorder, level-order)
- Python bindings available

**light_phylogeny** (https://crates.io/crates/light_phylogeny):
- **LCA (Lowest Common Ancestor) support** - critical for tronko
- Arena-based tree structure
- SVG visualization

### 5. Memory-Mapped I/O

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| Proposed mmap optimization | **memmap2** | Standard Rust mmap crate |

**memmap2** (https://crates.io/crates/memmap2):
- Cross-platform memory mapping
- Equivalent to POSIX `mmap()`
- The standard choice for Rust projects

### 6. Multi-threading

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| pthreads | **rayon** | Data parallelism with work-stealing |

**rayon** (https://crates.io/crates/rayon):
- Simple API: `.par_iter()` instead of `.iter()`
- Automatic load balancing
- **Compile-time data race prevention**

```rust
// Transform pthread-based batch processing to:
sequences.par_iter()
    .map(|seq| process_sequence(seq))
    .collect()
```

### 7. gzip Compression

| C Implementation | Rust Replacement | Notes |
|-----------------|------------------|-------|
| zlib via `gzopen()` | **flate2** with zlib-ng | 2-3x faster compression |

---

## Major Challenges

### Challenge 1: Tree Data Structures (HARDEST)

The tronko `node` structure has bidirectional parent/child references:

```c
typedef struct node {
    int up[1000];    // Children
    int down;        // Parent
    double **posteriornc; // [position][nucleotide] posteriors
    // ...
} node;
```

**Why this is hard in Rust:**
- Rust's ownership model disallows cyclic references
- A node referencing both parent and children creates ownership cycles
- This is famously described as "Learning Rust With Entirely Too Many Linked Lists"

**Solutions (in order of preference):**

1. **Arena allocation with indices** (Recommended):
   ```rust
   use indextree::{Arena, NodeId};

   struct PhyloNode {
       name: String,
       likelihood: Vec<f64>,
       posteriors: Vec<Vec<f64>>, // [position][nucleotide]
   }

   let arena: Arena<PhyloNode> = Arena::new();
   let root = arena.new_node(PhyloNode { ... });
   ```

2. **Vec + index-based approach**:
   ```rust
   struct Node {
       parent: Option<usize>,
       children: Vec<usize>,
       posteriors: Vec<Vec<f64>>,
   }
   struct Tree {
       nodes: Vec<Node>,
       root: usize,
   }
   ```

3. **Keep unsafe with raw pointers** (defeats purpose of porting)

### Challenge 2: No Rust ML Phylogenetic Inference

**RAxML** and **FastTree** have no Rust equivalents. Options:
- Continue calling as subprocesses (current approach)
- Create FFI bindings to RAxML-NG (C++)
- Wait for Rust implementations to emerge

### Challenge 3: FAMSA (Multiple Sequence Alignment)

FAMSA is heavily SIMD-optimized C++. No Rust equivalent exists.
- Must continue as subprocess call
- Could create FFI bindings (complex due to C++)

### Challenge 4: Likelihood Calculations

The GTR+gamma model in `likelihood.c` is mathematically complex:
- Eigenvalue decomposition
- Transition matrix construction
- Recursive tree traversal with underflow protection

**Porting considerations:**
- **nalgebra** crate for linear algebra (eigensystems)
- Careful translation of numerical code
- Extensive validation against original implementation

---

## Porting Strategy Options

### Option A: Incremental FFI Approach (Recommended)

**Phase 1 - FFI Foundation (2-4 weeks)**
- Create bindgen bindings for BWA, WFA2
- Set up Cargo workspace with mixed C/Rust build
- Keep all C code working, add Rust CLI wrapper

**Phase 2 - I/O Layer (2-4 weeks)**
- Replace FASTA parsing with needletail
- Replace gzip with flate2 + zlib-ng
- Add Rust CLI with clap

**Phase 3 - Data Structures (4-8 weeks)**
- Design arena-based tree structure
- Port tree parsing (Newick) using phylotree patterns
- **This is the hardest phase**

**Phase 4 - Core Algorithms (4-8 weeks)**
- Port likelihood calculations with nalgebra
- Port assignment logic
- Replace pthreads with Rayon

**Phase 5 - Remove C Dependencies (2-4 weeks)**
- Optionally replace rust-bwa with pure Rust FM-index
- Or keep FFI dependencies for production reliability

**Total: 3-6 months**

### Option B: Full c2rust Transpilation

1. Generate `compile_commands.json` from Makefile
2. Run `c2rust transpile` on entire codebase
3. Get tests passing with unsafe Rust
4. Incrementally refactor to safe Rust

**Pros:** Preserves functionality immediately
**Cons:** Results in non-idiomatic unsafe Rust; still requires significant refactoring

### Option C: Hybrid - Critical Paths Only

Port only the memory-critical paths to Rust:
- Posterior probability matrices (80-95% of memory)
- Batch processing loop
- Keep likelihood calculations in C via FFI

---

## Rust Library Mapping

| tronko Component | Rust Crate | Confidence |
|-----------------|------------|------------|
| FASTA/FASTQ parsing | needletail, seq_io | High |
| gzip handling | flate2 | High |
| Needleman-Wunsch alignment | rust-bio::alignment | High |
| WFA2 alignment | libwfa, rust_wfa | High |
| BWA alignment | rust-bwa | High |
| FM-index | rust-bio::data_structures::fmindex | High |
| Newick tree parsing | phylotree, light_phylogeny | High |
| LCA computation | light_phylogeny | Medium |
| Tree data structures | indextree, custom arena | Medium |
| Memory mapping | memmap2 | High |
| Parallelism | rayon | High |
| GTR likelihood | nalgebra + custom | Medium |
| RAxML | None (subprocess) | N/A |
| FAMSA | None (subprocess) | N/A |

---

## Performance Expectations

Based on benchmarks and real-world ports:

| Aspect | Expected Performance |
|--------|---------------------|
| FASTA parsing | Equal or faster |
| Sequence alignment | 5-10x faster with block-aligner SIMD |
| Tree traversal | Comparable |
| Overall throughput | Within 5-10% or better |
| Memory safety | Compile-time guarantees |
| Parallelization | Safer, similar speed |

**Key insight from rust-bio paper:**
> "Rust enables fearless parallelization that would be risky in C"

---

## Benefits of Porting

1. **Memory Safety**: No segmentation faults, buffer overflows
2. **Thread Safety**: Data races prevented at compile time
3. **Better Error Handling**: `Result<T, E>` vs C error codes
4. **Package Management**: Cargo vs manual dependency management
5. **Cross-Platform**: Single codebase for Linux/macOS/Windows
6. **Modern Tooling**: Built-in testing, documentation, benchmarking
7. **Active Ecosystem**: Growing bioinformatics community

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Tree structure complexity | Use arena allocation; study "Too Many Linked Lists" |
| Numerical precision differences | Extensive validation against C output |
| Performance regression | Benchmark continuously; use SIMD crates |
| Missing library functionality | Keep critical C code via FFI |
| Learning curve | Team Rust training; incremental approach |

---

## Recommendations

### Should You Port?

**Port if:**
- Memory safety is critical for your deployment
- You want to leverage Rust's parallelism safely
- You're planning significant new development
- Cross-platform distribution is important

**Don't port if:**
- Current C code is stable and meeting needs
- No resources for 3-6 month effort
- Team has no Rust experience and no time to learn

### Recommended Approach

1. **Start with tronko-assign** (more self-contained than tronko-build)
2. **Use incremental FFI approach** to maintain working code throughout
3. **Tackle tree structures early** - they're the biggest risk
4. **Keep RAxML/FAMSA as subprocesses** - no benefit to porting
5. **Leverage existing Rust crates** - don't reinvent needletail, phylotree, etc.

---

## Code References

- `tronko-assign/tronko-assign.c:744` - Main entry point
- `tronko-assign/readreference.c:311-433` - Reference tree loading (port to memmap2)
- `tronko-assign/placement.c:70-125` - WFA2 alignment (replace with libwfa)
- `tronko-assign/assignment.c:24-210` - Score computation (port to Rust)
- `tronko-assign/global.h:43-51` - Node structure (redesign with arena)
- `tronko-build/likelihood.c` - GTR+gamma calculations (port with nalgebra)

---

## Related Research

- `thoughts/shared/research/2025-12-29-tronko-assign-data-flow.md` - Data flow analysis this research builds upon

---

## Sources

### Rust Bioinformatics Libraries
- [rust-bio](https://rust-bio.github.io/) - Core bioinformatics library
- [needletail](https://github.com/onecodex/needletail) - Fast FASTA/FASTQ parsing
- [phylotree](https://github.com/lucblassel/phylotree-rs) - Phylogenetic tree handling
- [light_phylogeny](https://github.com/simonpenel/light_phylogeny) - LCA and tree visualization
- [rust-bwa](https://github.com/10XGenomics/rust-bwa) - BWA bindings by 10X Genomics
- [libwfa](https://lib.rs/crates/libwfa) - WFA2 Rust bindings
- [block-aligner](https://github.com/Daniel-Liu-c0deb0t/block-aligner) - SIMD alignment
- [minimap2-rs](https://github.com/jguhlin/minimap2-rs) - Minimap2 bindings
- [memmap2](https://github.com/RazrFalcon/memmap2-rs) - Memory mapping
- [rayon](https://github.com/rayon-rs/rayon) - Data parallelism

### C-to-Rust Porting
- [c2rust](https://github.com/immunant/c2rust) - Automated transpiler
- [rust-bindgen](https://github.com/rust-lang/rust-bindgen) - FFI binding generator
- [Learning Rust With Entirely Too Many Linked Lists](https://rust-unofficial.github.io/too-many-lists/)
- [Idiomatic Trees in Rust](https://rust-leipzig.github.io/architecture/2016/12/20/idiomatic-trees-in-rust/)

### Performance and Feasibility
- [Rust-Bio: a fast and safe bioinformatics library](https://academic.oup.com/bioinformatics/article/32/3/444/1743419) - Bioinformatics 2016
- [Why scientists are turning to Rust](https://www.nature.com/articles/d41586-020-03382-2) - Nature 2020
- [COMBINE-lab: Why use Rust for bioinformatics?](https://combine-lab.github.io/blog/2022/11/25/rust-for-bioinformatics-part-1.html)
- [Challenges Porting AV1 Decoder to Rust](https://www.infoq.com/news/2024/10/porting-av1-decoder-rust/) - InfoQ 2024

---

## Open Questions

1. What is the acceptable memory/performance trade-off for arena-based trees vs pointer-based?
2. Should BWA be kept as FFI or replaced with pure Rust FM-index (less mature)?
3. Is there interest in contributing Rust phylogenetic inference tools to the community?
4. Would a partial port (just the memory-critical paths) be sufficient?
