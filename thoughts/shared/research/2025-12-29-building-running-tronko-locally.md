---
date: 2025-12-29T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Building and Running Tronko Locally"
tags: [research, codebase, tronko-build, tronko-assign, build-system, testing]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
---

# Research: Building and Running Tronko Locally

**Date**: 2025-12-29
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question
What's involved with building and running tronko locally?

## Summary

Tronko consists of two C modules (tronko-build and tronko-assign) that can be built with standard `make` commands. Both use GCC with `-O3` optimization and require math, pthread, and zlib libraries. tronko-assign additionally requires the POSIX real-time library (`-lrt`). The project includes bundled executables for tree building/partitioning (famsa, raxmlHPC-PTHREADS, nw_reroot) and embeds BWA and WFA2 alignment libraries directly in the source. Example datasets are provided for testing single-tree and multi-cluster workflows.

## Detailed Findings

### System Requirements

#### Required Libraries
| Library | Flag | Purpose |
|---------|------|---------|
| Math library | `-lm` | Mathematical functions (log, exp, gamma) |
| POSIX threads | `-pthread` | Multi-threaded processing |
| zlib | `-lz` | Reading compressed FASTA/FASTQ files |
| POSIX real-time | `-lrt` | Timing functions (tronko-assign only) |

#### Compiler
- **GCC** with GNU C99 standard (`-std=gnu99`)
- Both modules compile all sources in a single compilation unit

### Building tronko-build

**Location**: `/home/jimjeffers/Work/tronko/tronko-build/`

**Build Commands**:
```bash
# Production build with -O3 optimization
cd tronko-build && make

# Debug build with symbols (-g flag)
cd tronko-build && make debug

# Clean build artifacts
cd tronko-build && make clean
```

**Source Files** (10 C files + 1 hashmap library):
| File | Purpose |
|------|---------|
| `tronko-build.c` | Main entry point and core logic |
| `getclade.c` | Newick tree parsing, clade extraction |
| `readfasta.c` | FASTA/MSA file parsing via zlib |
| `readreference.c` | Reference database I/O |
| `allocatetreememory.c` | Memory allocation for tree structures |
| `math.c` | Mathematical functions (gamma, eigenvectors) |
| `likelihood.c` | Phylogenetic likelihood calculations |
| `opt.c` | Numerical optimization |
| `options.c` | Command-line argument parsing |
| `printtree.c` | Tree output formatting |
| `hashmap.c` | Third-party hashmap implementation |

**Output**: `tronko-build` executable

### Building tronko-assign

**Location**: `/home/jimjeffers/Work/tronko/tronko-assign/`

**Build Commands**:
```bash
# Production build with -O3 optimization
cd tronko-assign && make

# Debug build with symbols
cd tronko-assign && make debug

# Clean build artifacts
cd tronko-assign && make clean
```

**Source Files** (~63 C files total):
- **Core tronko-assign** (13 files): Main logic, reference reading, assignment, placement, options parsing, logging, resource monitoring, crash debugging
- **Needleman-Wunsch alignment** (3 files): Classic global alignment from seq-align library
- **BWA library** (20 files): FM-index based read mapping (embedded)
- **WFA2 library** (31 files): Wavefront Alignment Algorithm v2 (embedded)
- **Hashmap** (1 file): Generic hashmap implementation

**Output**: `tronko-assign` executable

### Bundled External Tools

**Location**: `/home/jimjeffers/Work/tronko/bin/`

| Tool | Purpose | Used When |
|------|---------|-----------|
| `famsa` | Multiple sequence alignment | Partitioning workflow |
| `raxmlHPC-PTHREADS` | Phylogenetic tree estimation | Partitioning workflow |
| `nw_reroot` | Newick tree rerooting | After tree generation |
| `fasta2phyml.pl` | FASTA to PHYLIP conversion | Before RAxML |

These tools are only needed when using the partitioning option (`-y`) in tronko-build.

### Embedded Third-Party Libraries

#### WFA2 (Wavefront Alignment Algorithm v2)
- **Location**: `tronko-assign/WFA2/`
- **License**: MIT
- **Purpose**: Fast pairwise sequence alignment (default algorithm)
- **Files**: 68 source files

#### BWA (Burrows-Wheeler Aligner)
- **Location**: `tronko-assign/bwa_source_files/`
- **Purpose**: FM-index based initial read mapping to leaf nodes
- **Files**: 35 source files

#### Needleman-Wunsch
- **Location**: `tronko-assign/needleman_wunsch.c`
- **Source**: https://github.com/noporpoise/seq-align
- **License**: Public Domain
- **Purpose**: Global alignment (activated with `-w` flag)

### Testing with Example Datasets

#### Single Tree Test

**Build the reference database**:
```bash
./tronko-build/tronko-build -l \
  -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta \
  -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt \
  -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot \
  -d tronko-build/example_datasets/single_tree
```

**Run single-end assignment** (uses pre-built database):
```bash
./tronko-assign/tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/test_results.txt -w
```

**Run paired-end assignment**:
```bash
./tronko-assign/tronko-assign -r \
  -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -p -1 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta \
     -2 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta \
  -o /tmp/test_results.txt -w
```

#### Multiple Trees Test

**Build with multiple MSAs**:
```bash
./tronko-build/tronko-build \
  -d tronko-build/example_datasets/multiple_trees/multiple_MSA
```

**Note**: Multi-cluster builds require files to follow naming conventions:
- MSA: `[Number]_MSA.fasta`
- Taxonomy: `[Number]_taxonomy.txt`
- Tree: `RAxML_bestTree.[Number].reroot`

### Key Command-Line Options

#### tronko-build Options
| Flag | Description |
|------|-------------|
| `-l` | Single tree mode (use with -m, -x, -t) |
| `-m` | MSA FASTA file path |
| `-x` | Taxonomy file path |
| `-t` | Tree file path (Newick format) |
| `-d` | Output directory |
| `-y` | Partition mode (requires external tools in PATH) |
| `-a` | Use FastTree instead of RAxML (with -y) |

#### tronko-assign Options
| Flag | Description |
|------|-------------|
| `-f` | Reference database (reference_tree.txt) |
| `-a` | Reference FASTA for BWA index |
| `-s` | Single-end mode |
| `-g` | Single-end query FASTA |
| `-p` | Paired-end mode |
| `-1`, `-2` | Paired-end read files |
| `-q` | FASTQ input (default: FASTA) |
| `-r` | Print results header |
| `-o` | Output file path |
| `-w` | Use Needleman-Wunsch instead of WFA2 |

### Expected Test Results

The example datasets include expected results files for validation:
- `missingreads_singleend_150bp_2error_results.txt`
- `missingreads_pairedend_150bp_2error_results.txt`

Single tree test should assign all reads to: `Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge`

## Code References

- tronko-build Makefile: `tronko-build/Makefile`
- tronko-assign Makefile: `tronko-assign/Makefile`
- Global constants (build): `tronko-build/global.h:7-25`
- Global constants (assign): `tronko-assign/global.h:6-24`
- WFA2 integration: `tronko-assign/placement.c:70-116`
- BWA integration: `tronko-assign/tronko-assign.c:86-93`
- External tool invocations: `tronko-build/tronko-build.c:667-751`

## Architecture Documentation

### Build System Design
- Single compilation unit approach (all sources compiled together)
- No incremental builds (any change triggers full recompilation)
- Embedded libraries compiled directly rather than linked externally

### Alignment Pipeline (tronko-assign)
1. BWA-MEM maps reads to leaf nodes
2. WFA2 (default) or Needleman-Wunsch performs pairwise alignment
3. Placement algorithm traverses tree using likelihood calculations
4. Assignment score determines taxonomic classification

## Open Questions

- FastTree alternative to RAxML: What are the trade-offs?
- Thread count optimization for famsa (`--famsa-threads` option)
- Performance characteristics of WFA2 vs Needleman-Wunsch alignment
