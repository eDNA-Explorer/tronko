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
  -W [FLOAT]        Column gap mask threshold (default: 1.0 = no masking)
  -L, --legacy-sp   Use legacy SP normalization (divides by numspec, pre-fix behavior)
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

## tronko-assign Accuracy Tuning Options

These flags control how tronko-assign makes taxonomic assignments. All are optional and have safe defaults.

```
Accuracy Tuning:
  -c [FLOAT]                  LCA cutoff / Cinterval (default: 5)
  -u [FLOAT]                  Score constant for Jukes-Cantor correction (default: 0.01)
  --max-bwa-matches [INT]     Cap on candidate leaf alignments (default: 10)
  --best-leaf-threshold [FLOAT]   Best-leaf override score threshold (default: 0 = disabled)
  --best-leaf-max-votes [INT]     Max votes for best-leaf override (default: 0 = disabled)
  --adaptive-cinterval        Enable adaptive cinterval (default: disabled)
  --adaptive-gap-scale [FLOAT]    Scaling factor for adaptive cinterval (default: 0.5)

Aligner Selection:
  --aligner [STR]             'bwa' (default) or 'minimap2'
  --minimap2-kmer [INT]       minimap2 k-mer size (default: 15)
  --minimap2-window [INT]     minimap2 window size (default: 5)
```

## How Adaptive Cinterval Works

tronko-assign places a query by scoring every node in every candidate tree, then voting: all nodes within `Cinterval` log-likelihood units of the best score get a vote. The LCA of voted nodes determines the taxonomic assignment. A fixed `Cinterval` applies the same tolerance to every query, but different queries have fundamentally different score profiles:

- **Clear matches** (e.g. a haplotype that's in the database): the best leaf score is much higher than the second-best. A wide voting window adds irrelevant nodes and pushes the LCA up, losing species-level resolution.
- **Novel species** (not in the database): scores are spread across multiple distant leaves with no clear winner. A narrow window would pick one arbitrarily; a wider window correctly generalizes to genus level.

`--adaptive-cinterval` analyzes the gap between the top-1 and top-2 leaf scores before voting. When the gap is large (clear match), it shrinks the effective voting window for a more specific call. When scores are ambiguous, it keeps the window near the original `Cinterval` for a conservative call. The `--adaptive-gap-scale` parameter (0.0-1.0) controls how aggressively it shrinks: 0.0 = no adaptation, 1.0 = maximum shrinkage on clear matches.

This is opt-in and disabled by default. Without `--adaptive-cinterval`, behavior is identical to before.

## How Column Gap Masking Works

tronko-build computes posterior probabilities for every node at every alignment column. In diverse partitions (many species in one tree), alignments can have 60-80% gap fraction. Most of this comes from columns where the majority of sequences have gaps — these columns contribute noise to the likelihood calculation without carrying useful phylogenetic signal.

`-W <threshold>` masks columns where the gap fraction exceeds the threshold. Masked columns:
- Are skipped in the likelihood optimization (so model parameters are estimated from informative columns only)
- Get uniform posteriors (0.25, 0.25, 0.25, 0.25) in the output, which tronko-assign treats as uninformative

This preserves the tree's broad phylogenetic context (the tree topology and branch lengths are unchanged) while cleaning up the per-node posterior signal. The hypothesis is that this helps haplotype-level discrimination in diverse partitions without hurting species holdout performance.

Default `-W 1.0` means no masking (gap fraction can never exceed 100%). Useful values: `-W 0.5` masks columns with >50% gaps, `-W 0.7` masks only heavily gapped columns.

## File Naming Conventions

For multi-cluster builds, files must follow this naming:
- MSA: `[Number]_MSA.fasta`
- Taxonomy: `[Number]_taxonomy.txt`
- Tree: `RAxML_bestTree.[Number].reroot`

Taxonomy file format: `FASTA_header\tdomain;phylum;class;order;family;genus;species`
