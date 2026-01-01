---
date: 2025-12-31T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Reference database files and their usage in tronko-assign"
tags: [research, codebase, tronko-assign, bwa, fasta, taxonomy, reference-database]
status: complete
last_updated: 2025-12-31
last_updated_by: Claude
---

# Research: Reference Database Files and Their Usage in tronko-assign

**Date**: 2025-12-31T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

What is the purpose of the files in `example_datasets/16S_Bacteria/` beyond the reference_tree.txt/.trkb? Which files does tronko-assign utilize and why are they important?

## Summary

The `example_datasets/16S_Bacteria/` directory contains files for a complete tronko reference database:

| File | Used By | Purpose |
|------|---------|---------|
| `reference_tree.txt(.gz)` / `.trkb` | tronko-assign | Pre-computed phylogenetic likelihoods, tree topology, embedded taxonomy |
| `16S_Bacteria.fasta` | tronko-assign | Reference sequences for BWA alignment |
| `*.amb, *.ann, *.bwt, *.pac, *.sa` | tronko-assign (BWA) | Pre-built BWA index for fast alignment |
| `16S_Bacteria_taxonomy.txt` / `.tax.tsv` | **tronko-build only** | Input taxonomy for building reference database |

**Key insight**: tronko-assign requires BOTH the reference tree AND the FASTA file because:
1. The reference tree contains pre-computed likelihoods for phylogenetic placement
2. The FASTA file (+ BWA index) enables fast initial alignment to identify candidate leaf nodes

## Detailed Findings

### Reference Tree File (`reference_tree.txt` / `.trkb`)

The reference tree is a self-contained database created by tronko-build containing:

- **Tree topology**: Parent-child relationships (`up[0]`, `up[1]`, `down` fields)
- **Pre-computed posterior probabilities**: Fractional likelihoods for every node (not just leaves) - a 4-column matrix (A, C, G, T) for each alignment position
- **Embedded taxonomy**: 7-level taxonomy (domain→species) for each leaf, loaded at `readreference.c:578-591`
- **Leaf node names**: Sequence accession IDs that **must match** the FASTA headers

The tree format allows tronko-assign to perform LCA (Lowest Common Ancestor) calculations using the pre-computed likelihoods stored at every internal node.

### Reference FASTA File (`-a` option)

The FASTA file contains the actual nucleotide sequences for each leaf node. It serves a completely different purpose from the reference tree:

**Purpose**: Enable BWA (Burrows-Wheeler Aligner) to rapidly align query reads to identify candidate leaf nodes.

**Workflow**:
1. BWA aligns query sequences to the FASTA reference
2. BWA returns the names of reference sequences that match
3. A hash map (`fastmap.c:111-122`) maps sequence names to tree node indices
4. Phylogenetic placement then uses the likelihood data to find the LCA

**Why separate from tree file**: BWA requires raw sequence data in FASTA format to build its specialized index structures. The tree file stores phylogenetic/probabilistic data, not raw sequences.

**Critical requirement**: FASTA sequence headers must exactly match leaf node names in the tree file. The mapping code at `fastmap.c:215-228` looks up BWA results by name:
```c
leaf_map = hashmap_get(&map, read1);  // read1 = reference sequence name from BWA
aux->results[j-1].concordant_matches_roots[k] = leaf_map->root;
aux->results[j-1].concordant_matches_nodes[k] = leaf_map->node;
```

### BWA Index Files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`)

These are pre-built BWA index files, created by running `bwa index` on the FASTA file:

| Extension | Purpose |
|-----------|---------|
| `.bwt` | Burrows-Wheeler Transform - compressed suffix representation |
| `.sa` | Sampled Suffix Array - maps BWT positions to sequence positions |
| `.pac` | 2-bit packed sequence - used for alignment extension |
| `.ann` | Sequence annotations - names, lengths, offsets |
| `.amb` | Ambiguous bases - N positions in reference |

**Why pre-built matters**: Index building is computationally expensive. By including pre-built indices, users can skip the indexing step using `-6` flag:
```bash
tronko-assign -6 ...  # Skip BWA index building
```

Without `-6`, tronko-assign will rebuild the index at `tronko-assign.c:1039`:
```c
if (opt.skip_build==0){
    bwa_index(2,opt.fasta_file);
}
```

### Taxonomy Files (`*_taxonomy.txt`, `*.tax.tsv`)

**These are NOT used by tronko-assign at runtime.**

These files are input to **tronko-build** (via `-x` option), which embeds the taxonomy into the reference_tree.txt during database creation.

The `-x` option exists in tronko-assign's option parsing (`options.c:14`) but is **never used** - it's vestigial. All taxonomy data comes from the embedded data in the reference tree file.

## Code References

- `tronko-assign/options.c:46-52` - FASTA file option documentation
- `tronko-assign/options.c:69` - Skip BWA build option (`-6`)
- `tronko-assign/tronko-assign.c:1039` - BWA index building call
- `tronko-assign/bwa_source_files/fastmap.c:111-122` - Hash map creation for name→node mapping
- `tronko-assign/bwa_source_files/fastmap.c:787-789` - BWA index loading
- `tronko-assign/readreference.c:511-691` - Reference tree parsing
- `tronko-assign/readreference.c:578-591` - Taxonomy reading from embedded data
- `tronko-build/printtree.c:31-84` - Reference tree file creation

## Architecture Insights

### Two-Phase Assignment Process

1. **Alignment Phase** (BWA): Query reads are aligned to the reference FASTA to identify candidate leaf nodes. This is a fast initial filter using string matching.

2. **Phylogenetic Placement Phase** (Likelihood): Candidate matches are scored using the pre-computed posterior probabilities at tree nodes. The LCA is determined by climbing the tree and comparing likelihoods.

### Data Separation Rationale

The separation of FASTA and tree files reflects their different roles:
- **FASTA**: Raw sequence data for alignment algorithms (BWA's FM-index construction)
- **Tree**: Statistical/phylogenetic data for taxonomic inference (likelihood matrices)

This separation also allows flexibility - the same FASTA could theoretically be used with different tree builds, though in practice they should be consistent.

### Pre-computation Strategy

The key innovation of tronko is pre-computing fractional likelihoods at ALL nodes (not just leaves) during the build phase. This makes assignment fast - no likelihood calculations are needed at runtime, only lookups and comparisons.

## Open Questions

1. **Binary format performance**: How much faster is `.trkb` compared to `.txt.gz`? The binary format is documented at `readreference.c:701-981`.

2. **Multi-tree support**: The format supports multiple trees (`numberOfTrees` in header). When/why would multiple trees be used?

3. **Taxonomy file in example**: Why include `16S_Bacteria_taxonomy.txt` in the example_datasets if tronko-assign doesn't use it? (Answer: For users who want to rebuild the database)
