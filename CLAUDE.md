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

## tronko-build CLI Options

```
tronko-build [OPTIONS] -d [OUTPUT DIRECTORY]

Required:
  -d [DIRECTORY]    Output directory for reference database

Mode (pick one):
  -l                Single-tree mode (use with -t, -m, -x)
  -y                Partition mode (multi-cluster, use with -e, -n)

Single-tree mode (-l):
  -t [FILE]         Rooted phylogenetic tree (Newick)
  -m [FILE]         Multiple sequence alignment (FASTA, can be gzipped)
  -x [FILE]         Taxonomy file

Partition mode (-y):
  -e [DIRECTORY]    Input directory with cluster files
  -n [INT]          Number of clusters in input directory
  -b [INT]          Restart from partition number (default: 0)
  -s                Use sum-of-pairs score for partitioning
  -u [FLOAT]        SP-score threshold (default: 0.5)
  -v                Use minimum leaf node count for partitioning
  -f [INT]          Minimum leaf nodes threshold (use with -v)
General:
  -a                Use FastTree instead of RAxML
  -c [INT]          FAMSA threads (0 = auto-detect, default: 0)
  -g                Don't flag missing data
  -i [STRING]       Prefix for output partition filenames
  -p                Two-step build (partition only, then exit)
  -r                Remove unused trees (use with -p)
  -E                Export final subtrees to exported_subtrees/ directory
  -h                Show help
```

## build-tronko-db.sh Options

End-to-end pipeline script wrapping AncestralClust + tronko-build + BWA indexing.

```
build-tronko-db.sh -f <input.fasta> -t <taxonomy.txt> -o <output_dir> [OPTIONS]

Required:
  -f    Input reference FASTA (unaligned, single-line)
  -t    Taxonomy file (accession<TAB>lineage)
  -o    Output directory

Options:
  -p    Primer/marker name (default: marker)
  -T    Threads for FAMSA/tree inference (default: 8)
  -s    SP-score threshold for tronko-build partitioning (default: 0.1)
  -F    Use FastTree instead of RAxML
  -E    Export subtrees for ablation studies
  -C    AncestralClust cutoff — max seqs before clustering (default: 25000)
  -B    AncestralClust bin size — target seqs per cluster (default: 20000)
  -P    AncestralClust descendants parameter (default: 75)
  -J    Parallel jobs for Step 2 cluster processing (default: 1)
```

The `-J` flag runs up to J cluster alignments/trees concurrently in Step 2.

## Testing

Example datasets are provided for testing builds:

```bash
# Single tree test
tronko-build -l -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta \
  -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt \
  -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
  -d tronko-build/example_datasets/single_tree

# Multi-cluster test (sequential)
tronko-build -y -e tronko-build/example_datasets/multiple_trees/multiple_MSA \
  -n 5 -d /tmp/test_multi -s -u 0.1

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
