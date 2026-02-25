# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tronko is a phylogeny-based method for accurate community profiling of metabarcoding datasets. It consists of two C modules:

- **tronko-build**: Builds custom reference databases from phylogenetic trees, MSAs, and taxonomy files
- **tronko-assign**: Assigns species to query sequences using a tronko-build database

The key innovation is calculating LCA (Lowest Common Ancestor) using fractional likelihoods stored in all nodes of a phylogeny, not just leaf nodes.

## Build Commands

```bash
# Build tronko-build
cd tronko-build && make

# Build tronko-assign
cd tronko-assign && make

# Debug builds (with -g flag)
cd tronko-build && make debug
cd tronko-assign && make debug

# Clean builds
cd tronko-build && make clean
cd tronko-assign && make clean
```

Both use gcc with `-O3` optimization. tronko-assign links pthread, zlib, and rt libraries.

## Architecture

### tronko-build
Main source: `tronko-build/tronko-build.c`
- Reads phylogenetic trees (Newick), MSAs (FASTA), and taxonomy files
- Can operate in single-tree mode (`-l`) or partition mode (`-y`)
- Partitioning uses sum-of-pairs scoring or minimum leaf node thresholds
- Output: `reference_tree.txt` database file

Key modules:
- `likelihood.c` - Likelihood calculations for phylogenetic placement
- `readfasta.c` - FASTA file parsing
- `readreference.c` - Reference database I/O
- `opt.c` / `math.c` - Optimization and mathematical functions

### tronko-assign
Main source: `tronko-assign/tronko-assign.c`
- Uses BWA for alignment to leaf nodes (embedded in `bwa_source_files/`)
- Supports Wavefront Alignment (`WFA2/`) or Needleman-Wunsch alignment
- Handles paired-end (`-p`) or single-end (`-s`) reads
- Supports FASTA (default) or FASTQ (`-q`) input

Key modules:
- `placement.c` - Phylogenetic placement logic
- `assignment.c` - Species assignment
- `alignment.c`, `needleman_wunsch.c` - Sequence alignment
- `WFA2/` - Wavefront Alignment Algorithm v2 (third-party)

### External Dependencies
For partitioning in tronko-build (bundled in `bin/`):
- `raxmlHPC-PTHREADS` - Tree estimation
- `famsa` - Multiple sequence alignment
- `nw_reroot` - Newick utilities
- `fasta2phyml.pl` - Format conversion (also in `scripts/`)

## Testing

Example datasets are provided for testing builds:

```bash
# Single tree test
tronko-build -l -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta \
  -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt \
  -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
  -d tronko-build/example_datasets/single_tree

# Test assignment with single-end reads
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_results.txt -w
```

## File Naming Conventions

For multi-cluster builds, files must follow this naming:
- MSA: `[Number]_MSA.fasta`
- Taxonomy: `[Number]_taxonomy.txt`
- Tree: `RAxML_bestTree.[Number].reroot`

Taxonomy file format: `FASTA_header\tdomain;phylum;class;order;family;genus;species`
