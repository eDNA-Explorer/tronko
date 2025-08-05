# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Tronko is a phylogeny-based method for taxonomic classification of metabarcoding datasets. It consists of two main C programs:
- `tronko-build`: Creates reference databases from phylogenetic trees and sequence alignments
- `tronko-assign`: Assigns taxonomic classifications to query sequences using the reference databases

## Build Commands

### Building the binaries
```bash
# Build tronko-build
cd tronko-build
make

# Build tronko-assign  
cd ../tronko-assign
make

# Build both (clean builds)
cd tronko-build && make clean && make
cd ../tronko-assign && make clean && make

# Debug builds
cd tronko-build && make debug
cd ../tronko-assign && make debug
```

### Installation
```bash
# Copy binaries to PATH after building
cp tronko-build/tronko-build /usr/local/bin/
cp tronko-assign/tronko-assign /usr/local/bin/
```

## Testing Commands

### Running example datasets
```bash
# Single tree example with tronko-build
tronko-build -l \
  -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta \
  -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt \
  -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
  -d /full/path/to/output/directory

# Single-end read assignment example
tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o results.txt -w

# Paired-end read assignment example
tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -p -1 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta \
  -2 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta \
  -o results.txt -w
```

## Code Architecture

### Core Components

**tronko-build/**
- `tronko-build.c`: Main program for database building
- `readfasta.c`, `readreference.c`: Input file parsers
- `likelihood.c`: Node probability calculations  
- `getclade.c`: Taxonomic information management
- `printtree.c`: Output generation

**tronko-assign/**
- `tronko-assign.c`: Main program for sequence assignment
- `alignment.c`, `needleman_wunsch.c`: Sequence alignment modules
- `assignment.c`: Taxonomic classification logic
- `placement.c`: Phylogenetic tree placement
- BWA integration for initial sequence mapping
- WFA2 library for efficient alignment

### Memory Management
- Custom allocators in `allocatetreememory.c` and `allocateMemoryForResults.c`
- Manual memory management for C structures
- Thread-safe allocation for parallel processing

### External Dependencies
- **BWA**: Burrows-Wheeler Aligner (embedded source in `bwa_source_files/`)
- **WFA2**: Wavefront Alignment Algorithm (embedded in `WFA2/`)
- **HashMaps**: Custom implementation in `hashmap.c`

### Threading Model
- `tronko-build`: Single-threaded
- `tronko-assign`: Multi-threaded with `-C` parameter for thread count

## Key Configuration Files

### Input File Formats
- **MSA**: FASTA format multiple sequence alignment
- **Tree**: Newick format phylogenetic tree (must be rooted)
- **Taxonomy**: Tab-delimited format: `FASTA_header\tdomain;phylum;class;order;family;genus;species`

### Global Configuration
- Constants defined in `global.h`:
  - `MAXQUERYLENGTH 30000`: Maximum query sequence length
  - `MAX_NUMBEROFROOTS 20000`: Maximum number of reference trees
  - `STATESPACE 20`: Likelihood calculation categories

## Development Workflow

### Making Changes
1. Understand the two-phase workflow: database building → sequence assignment
2. Check existing documentation in `docs/` for detailed algorithm explanations
3. Both programs use custom data structures for trees and sequences
4. Memory allocation patterns follow manual C management
5. Thread safety required for any changes to `tronko-assign`

### Important Parameters
- **LCA cutoff** (`-c`): Score threshold for taxonomic assignment (default: 5)
- **Alignment method**: WFA2 (default) or Needleman-Wunsch (`-w`)
- **Score constant** (`-u`): Affects likelihood calculations (default: 0.01)

### File Dependencies
- Reference database format is binary and specific to tronko
- FASTA reference sequences required separately from database file
- Tree partitioning affects memory usage and performance