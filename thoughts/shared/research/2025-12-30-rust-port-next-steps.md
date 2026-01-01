---
date: 2025-12-30T00:00:00-08:00
researcher: Claude
git_commit: experimental
branch: experimental
repository: tronko
topic: "Rust Port Next Steps - After FASTA Parsing"
tags: [research, rust, porting, roadmap, reference-database]
status: active
last_updated: 2025-12-30
last_updated_by: Claude
---

# Research: Rust Port Next Steps

**Date**: 2025-12-30
**Status**: Active roadmap for incremental Rust port

## Current State

The Rust monorepo is set up with:
- `rust/crates/tronko-core` - Core library with FASTA/FASTQ parsing via needletail
- `rust/crates/tronko-rs` - CLI with `parse` and `validate` commands
- CI integration in `.github/workflows/build.yml`

**Progress: ~1-2% of tronko-assign functionality**

This document outlines the recommended next steps for porting.

---

## Recommended Next Step: Reference Database Parsing

### Why This Is the Right Next Step

1. **Low risk** - Pure I/O, like FASTA parsing
2. **Forces understanding** - Must learn tree data structures to parse them
3. **Enables validation** - Can verify Rust reads what C produces
4. **Foundation** - Required before any assignment logic can be ported
5. **Well-documented** - Existing research covers the format in detail

### What We Already Know

From `thoughts/shared/research/2025-12-29-reference-database-loading.md`:

**File Format (`reference_tree.txt`):**
```
Line 1:  <numberOfTrees>
Line 2:  <max_nodename>
Line 3:  <max_tax_name>
Line 4:  <max_lineTaxonomy>

Per tree:
<numbase>	<root>	<numspec>

Taxonomy lines (semicolon-delimited, 7 levels):
Larus heuglini;Larus;Laridae;Charadriiformes;Aves;Chordata;Eukaryota

Node data (for each of 2*numspec-1 nodes):
<tree> <node> <up0> <up1> <down> <depth> <taxIdx0> <taxIdx1> [name]
Followed by numbase lines of 4 posterior probabilities (A, C, G, T)
```

**Memory layout:**
- `posteriornc` arrays dominate at 90.5% of memory
- For 1466 species, 316 positions: ~39 MB total

### Research Needed

#### 1. Tree Data Structure Design

**Question**: What's the best Rust representation for the phylogenetic tree?

**Options to evaluate**:

| Approach | Crate | Pros | Cons |
|----------|-------|------|------|
| Arena + indices | `indextree` | Simple, no unsafe | Index-based access |
| Vec + indices | None (custom) | Full control | More boilerplate |
| Existing phylo | `light_phylogeny` | Has LCA support | May not fit our needs |

**Research tasks**:
- [ ] Read `indextree` documentation and examples
- [ ] Study `light_phylogeny` arena implementation
- [ ] Prototype small tree with posterior data attached
- [ ] Benchmark memory usage vs C implementation

**Key code to study**:
```bash
# C node structure
grep -A 20 "typedef struct node" tronko-assign/global.h

# How nodes are accessed during scoring
grep -B5 -A10 "posteriornc\[" tronko-assign/assignment.c
```

#### 2. Posterior Probability Storage

**Question**: How to store the `[numbase][4]` posterior matrices efficiently?

**Options**:

```rust
// Option A: Vec of Vecs (like C)
posteriornc: Vec<Vec<f64>>,  // or f32 per optimization research

// Option B: Flat Vec with stride (more cache-friendly)
posteriornc: Vec<f64>,  // access: posteriornc[pos * 4 + nucleotide]

// Option C: ndarray (numpy-like)
posteriornc: ndarray::Array2<f64>,
```

**Research tasks**:
- [ ] Benchmark flat vs nested Vec access patterns
- [ ] Evaluate `ndarray` crate for 2D arrays
- [ ] Test with actual access pattern from `getscore_Arr()`

#### 3. Gzip Handling

**Question**: How to handle gzip-compressed reference files?

The C code uses zlib's `gzopen/gzgets`. Options:
- `flate2` crate with `GzDecoder`
- `needletail` already handles this for FASTA (can we reuse patterns?)

**Research task**:
- [ ] Test `flate2::read::GzDecoder` with reference file

#### 4. Parsing Strategy

**Question**: Stream line-by-line or read entire file?

C approach: Line-by-line with `gzgets()`

Rust options:
- `BufReader::lines()` for streaming
- `read_to_string()` for small files
- Memory-mapped for large files (future optimization)

---

## Second Priority: Needleman-Wunsch Alignment

### Why Second

1. **Self-contained** - `needleman_wunsch.c` is ~300 lines, no dependencies
2. **Pure algorithm** - No I/O, no tree structures
3. **Learning exercise** - Good way to learn Rust patterns
4. **Replaceable** - `rust-bio` has implementations, but understanding helps

### Research Needed

#### 1. Compare with rust-bio

**Question**: Can we use `rust-bio::alignment::pairwise` directly?

```bash
# Study C implementation
head -100 tronko-assign/needleman_wunsch.c

# Key functions
grep "^int\|^void\|^double" tronko-assign/needleman_wunsch.c
```

**Research tasks**:
- [ ] Read rust-bio pairwise alignment API
- [ ] Compare scoring matrices used
- [ ] Check if gap penalties match
- [ ] Benchmark rust-bio vs C implementation

#### 2. SIMD Optimization

**Question**: Should we use `block-aligner` for SIMD acceleration?

The existing research notes block-aligner is 5-10x faster. Worth evaluating for batch alignment.

**Research task**:
- [ ] Evaluate `block-aligner` crate for our use case

---

## Third Priority: Tree Structures (Hardest)

### Why This Is Hard

From the feasibility research:
> "Rust's ownership model disallows cyclic references. A node referencing both parent and children creates ownership cycles."

This is the core challenge that will determine the port's success.

### Research Needed

#### 1. Study Existing Approaches

**Resources to read**:
- [ ] [Learning Rust With Entirely Too Many Linked Lists](https://rust-unofficial.github.io/too-many-lists/)
- [ ] [Idiomatic Trees in Rust](https://rust-leipzig.github.io/architecture/2016/12/20/idiomatic-trees-in-rust/)
- [ ] `light_phylogeny` source code (has arena-based phylogenetic trees)

#### 2. Prototype with Real Data

**Task**: Create minimal tree that can:
1. Parse a small Newick tree
2. Store posterior probabilities per node
3. Traverse from leaf to root (for LCA)
4. Access any node by index (for scoring)

---

## Fourth Priority: Scoring/Assignment Logic

### Depends On

- Reference database parsing (to load trees)
- Tree structures (to traverse and score)

### Research Needed

Once tree structures are designed:
- [ ] Port `getscore_Arr()` from `assignment.c`
- [ ] Port `assignScores_Arr_paired()` or single-end variant
- [ ] Validate output matches C for same input

---

## Not Porting (Keep as C/FFI or Subprocess)

Per the feasibility research:

| Component | Reason |
|-----------|--------|
| BWA | Use `rust-bwa` FFI bindings (maintained by 10X Genomics) |
| WFA2 | Use `libwfa` FFI bindings |
| RAxML | No Rust equivalent; keep as subprocess |
| FAMSA | No Rust equivalent; keep as subprocess |

---

## Immediate Action Items

### This Week

1. **Read** the `indextree` crate documentation
2. **Prototype** a simple arena-based tree with attached data
3. **Parse** the header of `reference_tree.txt` in Rust
4. **Validate** sequence count parsing matches C

### Code Exploration Commands

```bash
# Examine reference file format
head -20 tronko-build/example_datasets/single_tree/reference_tree.txt
zcat reference_tree.txt.gz | head -50  # if gzipped

# Study C reading code
grep -A 50 "readReferenceTree" tronko-assign/readreference.c

# Study C writing code (to understand format)
grep -A 50 "writeReferenceTree\|write.*reference" tronko-build/*.c

# Study node access patterns
grep -B2 -A5 "treeArr\[" tronko-assign/assignment.c
```

---

## Success Metrics

### For Reference Database Parsing

- [ ] Can parse `reference_tree.txt` header (tree count, max lengths)
- [ ] Can parse tree metadata (numbase, root, numspec)
- [ ] Can parse taxonomy lines
- [ ] Can parse node data (tree, node, up, down, depth, etc.)
- [ ] Can parse posterior probabilities
- [ ] Memory usage within 2x of C implementation
- [ ] Tree count matches C implementation output

### For Full Port

- [ ] Assignment output matches C for same input
- [ ] Performance within 20% of C (10% is target)
- [ ] All tests pass
- [ ] Clippy clean

---

## Related Research

- `2025-12-29-reference-database-loading.md` - Detailed format documentation
- `2025-12-29-rust-port-feasibility.md` - Library mapping and strategy
- `2025-12-29-tronko-assign-data-flow.md` - Overall data flow

---

## Open Questions

1. Should we support both gzipped and plain text reference files?
2. Should we add a binary format option for faster loading (per optimization research)?
3. How much should we diverge from C data structures to be more idiomatic Rust?
4. Should the tree structures be generic (reusable) or tronko-specific?
