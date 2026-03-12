# tronko

![Quick Tests](https://github.com/lpipes/tronko/workflows/Quick%20Tests/badge.svg)
![Comprehensive Tests](https://github.com/lpipes/tronko/workflows/Tests/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg)](https://doi.org/10.5281/zenodo.7407318)

A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets

In the tronko package there are two modules: `tronko-build` and `tronko-assign`. `tronko-build` is for building custom reference databases that tronko-assign uses as input. We have two reference databases currently available for download with `tronko-assign`. Cytochrome oxidase I (COI) which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GGWACWGGWTGAACWGTWTAYCCYCC` and reverse primer `TANACYTCnGGRTGNCCRAARAAYCA`. 16S which was custom built with <a href="https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools">CRUX</a> using forward primer `GTGCCAGCMGCCGCGGTAA` and reverse primer `GACTACHVGGGTATCTAATCC`.

Alignment-based and composition-based assignment methods calculate the lowest common ancestor (LCA) using data only in the leaf nodes of a phylogeny (A). The advantage of Tronko is that it stores fractional likelihoods in all nodes of a phylogeny and calculates the LCA based on all nodes in the tree (B).
<img src="https://github.com/lpipes/tronko/blob/main/Overview_Figure.jpg?raw=true">

# tronko-build
`tronko-build` is for building custom reference databases to be used with `tronko-assign`.

	tronko-build [OPTIONS] -d [OUTPUT DIRECTORY]
	
		-h, usage:
		-d [DIRECTORY], REQUIRED, full path to output directory
		-y, use a partition directory (you want to partition or you have multiple clusters)
		-l, use only single tree (do not partition)
		-t [FILE], compatible only with -l, rooted phylogenetic tree [FILE: Newick]
		-m [FILE], comptabile only with -l, multiple sequence alignment [FILE: FASTA] (can be gzipped)
		-x [FILE], taxonomy file [FILE: FASTA_header	domain;phylum;class;order;family;genus;species, use only with -l]
		-e [DIRECTORY], compatible only with -y, directory for reading multiple clusters
		-n [INT], compatible only with -y, number of partitions in read directory
		-b [INT], comptabile only with -y, restart partitions with partition number [default: 0]
		-s, compatible only with -y, partition using sum-of-pairs score [can't use with -f, use with -s]
		-u [FLOAT], compatible only with -y, minimum threshold for sum of pairs score [default: 0.5]
		-v, compatible only with -y, partition using minimum number of leaf nodes [can't use with -s, use with -f]
		-f [INT], don't partition less than the minimum number of leaf nodes [can't use with -s, use with -v, use only with -y]
		-g, don't flag missing data
		-c, [INT] Number of FAMSA threads to use (0 means use all threads) [default: 1]
		-p, break the db build into two steps
		-r, remove unused trees and copy trees from initial partition directory [can only be used with -p]
		-i, [STRING] set the prefix for output partitions in -d
		-a, use fasttree instead of RAxML

	Build-time options (passed to make):
		NUMCAT=[INT], Number of gamma rate categories for the nucleotide substitution model [default: 1]
			Example: make NUMCAT=4

# tronko-assign
`tronko-assign` is for species assignment of queries. It requires a `tronko-build` database.

	tronko-assign [OPTIONS] -r -f [TRONKO-BUILD DB FILE] -a [REF FASTA FILE] -o [OUTPUT FILE]
	
		-h, usage:
		-r, REQUIRED, use a reference
		-f [FILE], REQUIRED, path to reference database file, can be gzipped
		-a [FILE], REQUIRED, path to reference fasta file (for aligner index build)
		-o [FILE], REQUIRED, path to output file
		-p, use paired-end reads
		-s, use single reads
		-v, when using single reads, reverse-complement it
		-z, when using paired-end reads,  reverse-complement the second read
		-g [FILE], compatible only with -s, path to single-end reads file
		-1 [FILE], compatible only with -p, path to paired-end forward read file
		-2 [FILE], compatible only with -p, path to paired-end reverse read file
		-c [INT], LCA cut-off to use [default:5]
		-C [INT], number of cores [default:1] (use -C 1 for reproducible results)
		-L [INT], number of lines to read for assignment [default:50000]
		-P, print alignments to stdout
		-w, use Needleman-Wunsch Alignment Algorithm (default: WFA)
		-q, Query is FASTQ [default is FASTA]
		-e, Use only a portion of the reference sequences
		-n [INT], compatible only with -e, Padding (Number of bases) to use in the portion of the reference sequences
		-5 [FILE], Print tree number and leaf number and exit
		-6, Skip the bwa build if database already exists
		-u, Score constant [default: 0.01]
		-U, Print unassigned results
		-7, Print scores for all nodes [scores_all_nodes.txt]
		-3 [DIR], Save alignments to directory
		-4 [DIR], Print trees to directory
		-V [LEVEL], Enable verbose logging [0=ERROR, 1=WARN, 2=INFO, 3=DEBUG] [default: disabled]
		-l [FILE], Log file path [default: stderr only]
		-R, Enable resource monitoring (memory/CPU usage)
		-T, Enable timing information
		--tsv-log [FILE], Export memory stats to TSV file for analysis
		--aligner [STR], Aligner for leaf matching: "bwa" (default) or "minimap2"
		--minimap2-kmer [INT], minimap2 k-mer size [default: 15]
		--minimap2-window [INT], minimap2 minimizer window size [default: 5]
		--max-bwa-matches [INT], Maximum leaf matches per read [default: 10]
		--best-leaf-threshold [FLOAT], Score threshold for best-leaf override [default: 0, disabled]
		--best-leaf-max-votes [INT], Max total votes for best-leaf override [default: 0, disabled]
		--early-termination, Enable early termination of scoring
		--no-early-termination, Disable early termination (default)
		--strike-box [FLOAT], Strike box size as multiplier of Cinterval [default: 1.0]
		--max-strikes [INT], Max strikes before early termination [default: 6]
		--enable-pruning, Enable subtree pruning during scoring
		--disable-pruning, Disable subtree pruning (default)
		--pruning-factor [FLOAT], Pruning threshold as multiplier of Cinterval [default: 2.0]
		--trace-read [STR], Print diagnostic trace for a specific read name

## Tunable Parameters Reference

`tronko-assign` has several groups of tunable parameters that affect classification accuracy and performance. This section describes each parameter, what it does, and how parameters interact with each other. All defaults preserve the original tronko-assign behavior.

Use `--trace-read "READNAME"` (or `--trace-read "*"` for all reads) with any parameter combination to see detailed per-read diagnostics on stderr.

---

### Parameter Interaction Map

Parameters fall into **4 independent groups** that do not interact with each other, plus one build-time parameter. Within each group, parameters may interact and should be grid-searched together.

```
GROUP 1: Aligner Selection          GROUP 2: Candidate Cap
  --aligner                           --max-bwa-matches
  --minimap2-kmer
  --minimap2-window

GROUP 3: Voting & LCA               GROUP 4: Scoring Speedups
  --best-leaf-threshold                --early-termination
  --best-leaf-max-votes                --strike-box
  -c (Cinterval)                       --max-strikes
                                       --enable-pruning
                                       --pruning-factor

BUILD-TIME (tronko-build):
  NUMCAT (gamma rate categories)
```

**Cross-group interactions are minimal.** Group 1 determines *which* leaf candidates are found; Group 2 caps *how many* are kept; Group 3 controls *how votes are tallied and the best-leaf override*; Group 4 trades *scoring speed for completeness*. Each group can be benchmarked independently against a baseline, then promising values combined.

The one exception: `--max-bwa-matches` (Group 2) interacts weakly with Group 1 — if the aligner finds more candidates, a higher cap matters more. But the interaction is monotonic (more matches + higher cap = more information), so they don't need a full cross-product grid.

---

### Group 1: Aligner Selection (Independent)

These parameters control which aligner finds leaf candidates. They affect the *set of candidate leaves* that enter the scoring pipeline but do not affect scoring, voting, or LCA logic.

#### `--aligner bwa|minimap2`

Selects the aligner used to match query reads against leaf reference sequences.

| Aligner | Indexing | Strengths | Weaknesses |
|---------|---------|-----------|------------|
| `bwa` (default) | BWT/suffix array | Well-tested, stable | Misses candidates at >7% divergence |
| `minimap2` | Minimizer hash | Better sensitivity at high divergence, faster | Newer code path |

BWA uses exact-match seeds (MEMs) of 17+ bp, so reads with clustered mutations can fail to seed. minimap2 uses (w,k)-minimizers which tolerate mutations in non-sampled positions, finding candidates at up to ~15% divergence.

```bash
# Default (BWA)
tronko-assign -r -f ref.txt -a ref.fasta -s -g reads.fasta -o out.txt -w

# Use minimap2
tronko-assign --aligner minimap2 -r -f ref.txt -a ref.fasta -s -g reads.fasta -o out.txt -w
```

When `--aligner minimap2` is set, the BWA index build step (`bwa index`) is automatically skipped. The `-a` FASTA file is still required (minimap2 builds its index from it at runtime).

#### `--minimap2-kmer` and `--minimap2-window`

**Only relevant when `--aligner minimap2`.** Ignored when using BWA.

These control minimap2's minimizer indexing. A minimizer is the smallest k-mer hash within each sliding window of `w` consecutive k-mers. Together they determine the density and sensitivity of the index.

| Parameter | Default | Effect of increasing | Effect of decreasing |
|-----------|---------|---------------------|---------------------|
| `--minimap2-kmer` | 15 | More specific seeds, fewer spurious matches, may miss divergent hits | More sensitive to divergent sequences, more spurious matches |
| `--minimap2-window` | 5 | Sparser index, faster but less sensitive | Denser index, slower but more sensitive |

**These two interact with each other** and should be grid-searched together:

```bash
# Suggested grid for benchmarking (kmer x window):
--minimap2-kmer 11 --minimap2-window 3   # Most sensitive, slowest
--minimap2-kmer 11 --minimap2-window 5
--minimap2-kmer 13 --minimap2-window 3
--minimap2-kmer 13 --minimap2-window 5
--minimap2-kmer 15 --minimap2-window 3
--minimap2-kmer 15 --minimap2-window 5   # Default
--minimap2-kmer 15 --minimap2-window 10  # Fastest, least sensitive
```

For amplicon databases (references ~250bp, queries ~150bp), k=15 w=5 is a good starting point. Lowering k to 11-13 may help for highly divergent queries (>10% from nearest reference).

---

### Group 2: Candidate Cap (Independent)

#### `--max-bwa-matches`

Controls the maximum number of leaf matches (from either aligner) considered per read. After the aligner returns hits, only up to this many unique tree-leaf pairs are kept for downstream scoring.

| Value | Behavior |
|-------|----------|
| 10 (default) | Original behavior. May silently drop relevant candidates for reads that match many leaves. |
| 20-50 | More candidates scored, better accuracy for multi-tree databases, slightly slower. |
| 100+ | Diminishing returns; most reads match <20 unique leaves. |

```bash
tronko-assign --max-bwa-matches 30 [other options...]
```

**This parameter is independent** of all others. It only controls how many candidates pass through to scoring. Increasing it never hurts accuracy (just costs time), so it can be benchmarked as a simple sweep:

```bash
# Suggested sweep:
--max-bwa-matches 10   # Default
--max-bwa-matches 20
--max-bwa-matches 50
--max-bwa-matches 100
```

---

### Group 3: Voting & LCA (Grid-search together)

These parameters affect how phylogenetic scores are converted into taxonomic assignments. They interact because they operate on the same voting pipeline: scores → votes → LCA → final taxonomy.

#### `-c` (Cinterval / LCA cut-off)

The confidence interval for likelihood scoring. A node receives votes if its score is within `Cinterval` of the best score. Larger values spread votes across more nodes (more conservative, higher-level assignments); smaller values concentrate votes on the best-scoring nodes (more specific assignments).

```bash
-c 5    # Default — moderate confidence interval
-c 1    # Tight — only nodes very close to the best score get votes
-c 10   # Wide — many nodes get votes, assignments tend toward higher taxonomy
```

**Interacts with**: `--best-leaf-threshold` (the threshold is a score value that must be considered relative to Cinterval).

#### `--best-leaf-threshold` and `--best-leaf-max-votes`

Best-leaf override: when the single best-scoring leaf has a score exceeding `--best-leaf-threshold` and the total number of votes is below `--best-leaf-max-votes`, override the LCA result with the best leaf's taxonomy. This is a specificity optimization for reads that clearly belong to one species but whose LCA is pulled up by noise votes.

Both must be set to non-zero to enable the feature. Setting either to 0 (default) disables it entirely.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--best-leaf-threshold` | 0 (disabled) | Minimum score gap between best leaf and 2nd-best leaf. A leaf is "best" if its score is this much better. |
| `--best-leaf-max-votes` | 0 (disabled) | Only apply override when total votes are below this. Prevents overriding high-confidence LCA results. |

```bash
# Enable best-leaf override when gap > 5.0 and votes < 20
tronko-assign --best-leaf-threshold 5.0 --best-leaf-max-votes 20 [options...]
```

**These two interact tightly with each other** (both must be set) and with `-c` (which determines vote counts and score ranges). Grid-search all three together:

```bash
# Suggested grid (threshold x max-votes, at each Cinterval):
-c 5 --best-leaf-threshold 0 --best-leaf-max-votes 0       # Disabled (baseline)
-c 5 --best-leaf-threshold 3.0 --best-leaf-max-votes 10
-c 5 --best-leaf-threshold 5.0 --best-leaf-max-votes 20
-c 5 --best-leaf-threshold 10.0 --best-leaf-max-votes 50
```

#### Group 3 grid search strategy

Since these parameters interact, the recommended approach is:

1. **First**, sweep `-c` alone (e.g., 1, 3, 5, 10) with all others at defaults
2. **Then**, at the best `-c`, sweep `--best-leaf-threshold` + `--best-leaf-max-votes`

```bash
# Compact grid: 4 x 4 = 16 combinations per aligner
-c:            1, 3, 5, 10
--best-leaf:   off, (threshold=-2.0, max-votes=5), (threshold=5.0, max-votes=20), (threshold=10.0, max-votes=50)
```

---

### Group 4: Scoring Speedups (Independent)

These parameters trade scoring thoroughness for speed. They affect how exhaustively the phylogenetic likelihood is computed for each candidate leaf but do **not** change the voting or LCA logic. They are independent of Groups 1-3.

#### `--early-termination`, `--strike-box`, `--max-strikes`

When enabled, scoring for a candidate stops early if the score is clearly worse than the current best. A "strike" occurs when a position's contribution falls outside the strike box (a multiple of Cinterval). After `max-strikes` strikes, scoring terminates for that candidate.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--early-termination` | disabled | Enable the feature |
| `--strike-box` | 1.0 | Strike box size as multiplier of Cinterval. Larger = more lenient (fewer false terminations). |
| `--max-strikes` | 6 | Strikes before terminating. More = safer but slower. |

```bash
# Enable with defaults
tronko-assign --early-termination [options...]

# More aggressive (faster, may miss some candidates)
tronko-assign --early-termination --strike-box 0.5 --max-strikes 3 [options...]

# More conservative (safer, still faster than no termination)
tronko-assign --early-termination --strike-box 2.0 --max-strikes 10 [options...]
```

**These three interact with each other** and should be tuned together. They do NOT interact with Groups 1-3.

#### `--enable-pruning` and `--pruning-factor`

When enabled, entire subtrees are skipped during scoring if the partial score at the subtree root is worse than `pruning_factor * Cinterval` below the current best. This avoids scoring nodes that cannot possibly affect the final LCA.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--enable-pruning` | disabled | Enable subtree pruning |
| `--pruning-factor` | 2.0 | Pruning threshold multiplier. Larger = more aggressive pruning (faster but riskier). |

```bash
tronko-assign --enable-pruning --pruning-factor 2.0 [options...]
```

**Interacts with**: `--early-termination` (both are scoring shortcuts, can be combined). Does NOT interact with Groups 1-3.

#### Group 4 grid search strategy

These are speed optimizations. Benchmark them for **speed vs accuracy tradeoff**:

```bash
# Baseline (no speedups)
(no flags)

# Early termination sweep:
--early-termination --strike-box 0.5 --max-strikes 3
--early-termination --strike-box 1.0 --max-strikes 6
--early-termination --strike-box 2.0 --max-strikes 10

# Pruning sweep:
--enable-pruning --pruning-factor 1.5
--enable-pruning --pruning-factor 2.0
--enable-pruning --pruning-factor 3.0

# Combined:
--early-termination --strike-box 1.0 --max-strikes 6 --enable-pruning --pruning-factor 2.0
```

---

### Build-time Parameter: NUMCAT (Independent)

`NUMCAT` is set at **tronko-build compile time**, not at tronko-assign runtime. It controls the number of gamma rate categories in the nucleotide substitution model used when building the reference database. Changing it requires rebuilding the database.

| Value | Behavior |
|-------|----------|
| 1 (default) | Single rate category. Fastest build, original behavior. |
| 4 | Four discrete gamma categories. Slower build, may improve likelihood accuracy for rate-variable markers. |

```bash
# Build tronko-build with 4 gamma rate categories
cd tronko-build && make NUMCAT=4

# Then rebuild the database
tronko-build -l -m MSA.fasta -x taxonomy.txt -t tree.nwk -d output_dir/
```

**Completely independent** of all tronko-assign runtime parameters. A database built with NUMCAT=4 can be tested with any combination of tronko-assign parameters.

---

### Diagnostic Tool: `--trace-read`

Not a tuning parameter itself, but essential for understanding parameter effects. Prints detailed per-read diagnostics to stderr.

```bash
# Trace a specific read
tronko-assign --trace-read "GU572157.1_8_1" [other options...] 2>trace.log

# Trace ALL reads (verbose — use only on small test sets)
tronko-assign --trace-read "*" [other options...] 2>trace.log
```

Example trace output:
```
TRACE read=GU572157.1_0
  BWA hits: 10 unique trees (raw: 0 concordant + 10 discordant, max_bwa_matches=10)
    hit[0]: tree=0 node=2371 leaf=GU571681.1 numspec=1466
    hit[1]: tree=0 node=1903 leaf=DQ434794.1 numspec=1466
    ...
  SCORING: max_node=tree0:node1092 score=-29.98 is_leaf=0
  BEST_LEAF: tree0:node1903 score=-29.98 gap=0.00 tax=Uria aalge
  VOTES: 90 leaves + 60 internal = 150 total (Cinterval=5.0)
```

---

### Quick Reference: Python Benchmarking Sweep

For a systematic benchmark, here is a suggested parameter sweep organized by group. Run each group independently against a baseline, then combine winning values.

**Baseline command:**
```bash
tronko-assign -r -f ref.txt -a ref.fasta -s -g reads.fasta -o out.txt -w
```

**Group 1 (Aligner) — 7 runs:**
```bash
# BWA baseline (already covered by baseline)
--aligner minimap2 --minimap2-kmer 11 --minimap2-window 3
--aligner minimap2 --minimap2-kmer 13 --minimap2-window 5
--aligner minimap2 --minimap2-kmer 15 --minimap2-window 5   # minimap2 default
--aligner minimap2 --minimap2-kmer 15 --minimap2-window 10
--aligner minimap2 --minimap2-kmer 11 --minimap2-window 5
--aligner minimap2 --minimap2-kmer 13 --minimap2-window 3
```

**Group 2 (Candidate cap) — 4 runs:**
```bash
--max-bwa-matches 10   # default
--max-bwa-matches 20
--max-bwa-matches 50
--max-bwa-matches 100
```

**Group 3 (Voting/LCA) — sequential sweeps, ~8-16 runs total:**
```bash
# Step 1: Cinterval sweep (4 runs)
-c 1 / -c 3 / -c 5 / -c 10

# Step 2: At best -c, best-leaf sweep (4 runs)
(disabled) / --best-leaf-threshold 3 --best-leaf-max-votes 10 / 5,20 / 10,50
```

**Group 4 (Speed) — 7 runs:**
```bash
# Baseline (no speedups)
--early-termination --strike-box 0.5 --max-strikes 3
--early-termination --strike-box 1.0 --max-strikes 6
--early-termination --strike-box 2.0 --max-strikes 10
--enable-pruning --pruning-factor 1.5
--enable-pruning --pruning-factor 2.0
--early-termination --enable-pruning   # combined with defaults
```

**Build-time (NUMCAT) — 2 database builds:**
```bash
cd tronko-build && make              # NUMCAT=1 (default)
cd tronko-build && make NUMCAT=4     # 4 gamma categories
# Rebuild database with each, then test with best runtime params
```

## Verbose Logging and Performance Monitoring

`tronko-assign` includes comprehensive logging capabilities to monitor performance and troubleshoot issues:

### Basic Verbose Logging
```bash
# Enable INFO level logging to stderr
tronko-assign -V2 [other options...]

# Enable DEBUG level logging with maximum detail
tronko-assign -V3 [other options...]
```

### Logging to File
```bash
# Write logs to a file (also logs to stderr)
tronko-assign -V2 -l assignment.log [other options...]
```

### Resource Monitoring
```bash
# Monitor memory and CPU usage at each milestone
tronko-assign -V2 -R [other options...]
```

### Timing Information
```bash
# Include timing information between processing milestones
tronko-assign -V2 -T [other options...]
```

### Comprehensive Monitoring
```bash
# Enable all logging features for detailed performance analysis
tronko-assign -V3 -T -R -l comprehensive.log [other options...]
```

The logging system tracks 18 key milestones throughout the assignment process including:
- Program startup and option parsing
- Reference database loading 
- BWA index creation
- Memory allocation
- Batch processing start/completion
- Sequence alignment and placement
- Results writing and cleanup

Log levels:
- `-V0`: ERROR messages only
- `-V1`: WARN and ERROR messages  
- `-V2`: INFO, WARN, and ERROR messages (recommended)
- `-V3`: DEBUG, INFO, WARN, and ERROR messages (most verbose)

Tronko uses the <a href="https://github.com/smarco/WFA2-lib">Wavefront Alignment Algorithm (version 2)</a> or <a href="https://github.com/noporpoise/seq-align">Needleman-Wunsch Algorithm</a> for semi-global alignments. For initial leaf candidate matching, it supports two aligners: <a href="https://github.com/lh3/bwa">BWA-MEM</a> (default) and <a href="https://github.com/lh3/minimap2">minimap2</a> (via `--aligner minimap2`). Both are compiled into the binary; no external tools are needed. It uses <a href="https://github.com/DavidLeeds/hashmap">David Leeds' hashmap</a> for hashmap implementation in C. `tronko-assign` does not reverse complement your reads automatically. You must use options `-v` or `-z` to reverse complement your read for better alignment to the reference database. For more information on the direction of your reads based on your library prep, please refer to this helpful blog here: <a href="http://onetipperday.blogspot.com/2012/07/how-to-tell-which-library-type-to-use.html">http://onetipperday.blogspot.com/2012/07/how-to-tell-which-library-type-to-use.html</a>.

## Example output

The output file is a tab-delimited text file where only the forward readname is retained (if using paired-end reads). The output displays the taxonomic path for assignment, the score, the number of forward read mismatches with the `bwa` hit, the number of reverse read mismatches with the `bwa` hit, the tree number for the best assignment (0 if using 1 tree), and the node number the read (or reads in the case of paired-end reads) was assigned to. For single-end reads, the `Reverse_Mismatch` will always be 0 and the `Forward_Mismatch` is the number of read mismatches with the `bwa` hit.

```
Readname	Taxonomic_Path	Score	Forward_Mismatch	Reverse_Mismatch	Tree_Number	Node_Number
GU572157.1_8_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-54.258690	5.000000	4.00000	0	1095
GU572157.1_7_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-42.871226	1.000000	6.00000	0	1095
GU572157.1_6_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-59.952407	7.000000	3.00000	0	1098
GU572157.1_5_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-31.483761	2.000000	4.00000	0	1095
GU572157.1_4_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-31.483761	2.000000	3.00000	0	1095
GU572157.1_3_1	Eukaryota;Chordata;Aves;Charadriiformes;Alcidae;Uria;Uria aalge	-54.258690	2.000000	7.00000	0	1095
```

# INSTALLATION

1. Clone the [GitHub repo](https://github.com/lpipes/tronko), e.g. with `git clone https://github.com/lpipes/tronko.git`
2. Install build dependencies:
   - **Debian/Ubuntu**: `sudo apt-get install -y zlib1g-dev libzstd-dev`
   - **RHEL/CentOS/Fedora**: `sudo dnf install -y zlib-devel libzstd-devel`
   - **macOS (Homebrew)**: `brew install zlib zstd`
3. Run `make` in the `tronko-build` and `tronko-assign` directories.
4. Copy the `tronko-build` and `tronko-assign` binaries to your path.

`tronko-assign` uses [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads), [zlib](http://en.wikipedia.org/wiki/Zlib), and [zstd](https://github.com/facebook/zstd) as its dependencies.

## Supported Input Formats

`tronko-assign` accepts FASTA/FASTQ query files in the following formats:
- **Plain text** (`.fasta`, `.fastq`)
- **Gzip compressed** (`.fasta.gz`, `.fastq.gz`)
- **Zstandard compressed** (`.fasta.zst`, `.fastq.zst`)

Compression format is auto-detected from file magic bytes, not file extension. This allows for streaming decompression with minimal memory overhead (~256 KB for zstd buffers).

`tronko-build` only has dependencies if using the partition procedure for the initial trees. These are <a href="https://github.com/stamatak/standard-RAxML">`raxmlHPC-PTHREADS`</a>, <a href="https://github.com/refresh-bio/FAMSA">`famsa`</a>, `nw_reroot` from <a href="https://anaconda.org/bioconda/newick_utils/files">Newick utilties</a>, <a href="https://raw.githubusercontent.com/lpipes/tronko/main/scripts/fasta2phyml.pl">`fasta2phyml.pl`</a>, and <a href="https://ftp.gnu.org/gnu/sed/">`sed`</a>, which must be installed in your path. By default, `tronko-build` uses `raxmlHPC-PTHREADS` to estimate trees, if you choose the `-a` option to use `FastTree` instead, you must have `FastTree` installed in your path.

	cd tronko/tronko-build
	make
	../tronko-assign
	make

# SINGULARITY CONTAINER

To use the container download:

	singularity pull library://lpipes/tronko/tronko:1.0

To run `tronko-assign` with the container:

	singularity exec --bind <root-dir-to-bind-to-container> tronko_1.0.sif tronko-assign

To run `tronko-build` with the container:

	singularity exec --bind <root-dir-to-bind-to-container> tronko_1.0.sif tronko-build

# `tronko-assign` Usage
Tronko does not detect the correct orientation of the reads. If your reverse read needs to be reverse complemented use the option `-z`. The default options of Tronko assume that your reads are in FASTA format. If you want to assign reads in FASTQ format, use the option `-q`. You will also need a FASTA file (not gzipped) of all of your reference sequences in the reference database (use the option `-a`). When using BWA (default), `tronko-assign` will create a `bwa index` of the reference sequences with the extension of *.fasta.ann, etc. If you already have the `bwa index` files present in the same directory and naming scheme as your reference sequences, you can choose skip the `bwa index` build use `-6`. When using minimap2 (`--aligner minimap2`), the index is built at runtime from the FASTA file and `-6` is not needed. The reads (and reference database file) can be gzipped or not gzipped. Assigning paired-end reads in FASTA format:
```
tronko-assign -r -f [tronko-build REFERENCE DB FILE] -p -1 [FORWARD READS FASTA] -2 [REVERSE READS FASTA] -a [REFERENCE SEQUENCES FASTA] -o [OUTPUT FILE]
```
Assigning single-end reads in FASTA format:
```
tronko-assign -r -f [tronko-build REFERENCE DB FILE] -s -g [READS FASTA] -a [REFERENCE SEQUENCES FASTA] -o [OUTPUT FILE]
```

# `tronko-build` Usage

## `tronko-build` Simple Usage (using 1 phylogenetic tree)
```
tronko-build -l -t [Rooted Newick Tree] -m [Multiple Sequence Alignment FASTA] -x [TAXONOMY FILE] -d [FULL PATH TO OUTPUT DIRECTORY] 
```
The taxonomy file is a `.txt` file that has the following format:
```
FASTA_header\tdomain;phylum;class;order;family;genus;species
```
The tree file, MSA file, and the taxonomy file must all contain identical corresponding names. The MSA file should not contain any line breaks. To remove line breaks we recommend
```
sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' test.fasta
```
## `tronko-build` example datasets

An example dataset for a single tree is provided in `tronko-build/example_datasets/single_tree`. The dataset includes 1466 COI sequences from the order Charadriiformes. The MSA is `tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`, the taxonomy file is `tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt`, and the tree file is `tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot`. To build the reference database for `tronko-assign` for this dataset, run `tronko-build` with the following command (`-l` specifies using a single tree):
```
tronko-build -l -m tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta -x tronko-build/example_datasets/single_tree/Charadriiformes_taxonomy.txt -t tronko-build/example_datasets/single_tree/RAxML_bestTree.Charadriiformes.reroot -d tronko-build/example_datasets/single_tree
```
An example dataset for multiple trees is provided in `tronko-build/example_datasets/multiple_trees/one_MSA`. The dataset includes 99 COI sequences from different species. The MSA is `tronko-build/example_datasets/multiple_trees/one_MSA/1_MSA.fasta`, the taxonomy file is `tronko-build/example_datasets/multiple_trees/one_MSA/1_taxonomy.txt`, and the tree file is `tronko-build/example_datasets/multiple_trees/one_MSA/RAxML_bestTree.1.reroot`. To build the reference database and partition it using the Sum-of-Pairs score (see manuscript), run `tronko-build` with the following command [-d is REQUIRED and must contain the full path NOT relative path):
```
tronko-build -y -e tronko-build/example_datasets/multiple_trees/one_MSA -n 1 -d [FULL PATH TO OUTPUT DIRECTORY] -s
```

An example dataset to not partition multiple tree is provided in `tronko-build/example_datasets/multiple_trees/multiple_MSA`. The dataset include 99 COI sequences from different species. To build the reference database and not partition the database, set the `-f` parameter higher than the number of sequences in the dataset. Run `tronko-build` with the following command [Make sure `-d` is supplied with the full path to the output directory]:
```
tronko-build -y -e tronko-build/example_datasets/multiple_trees/multiple_MSA -n 5 -d outdir_multiple_MSA -v -f 500
```

A successful run will produce a `tronko-assign` reference database named `reference_tree.txt` in the output directory. The `reference_tree.txt` database file is what is used to run `tronko-assign`.

## `tronko-assign` example datasets

An example dataset for `tronko-assign` is provided in `example_datasets/single_tree`. This dataset contains single-end reads (`example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta`) and paired-end reads (forward read: `example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta` and reverse read: `example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta`). For the single-end reads, there are 164 150bp reads with 2% simulated error/polymorphisms from the following sequence (taxonomically classified on NCBI as Uria aalge):

```
>GU572157.1
CCTGGCTGGTAATCTAGCCCATGCCGGAGCTTCAGTGGATTTAGCAATCTTCTCCCTTCACTTAGCAGGTGTATCATCTATTCTAGGCGCTATCAACTTTATCACAACAGCCATCAACATAAAGCCTCCAGCCCTCTCACAATACCAAACCCCCCTATTCGTATGATCAGTACTTATCACTGCTGTCCTACTACTACTCTCACTCCCAGTACTTGCTGCTGGTATCACTATATTACTAACAGATCGAAACTTAAACACAACATTCTTTGATCCAGCTGGAGGTGGTGACCCAGTACTTTACCAACACCTCTTC
``` 

This particular sequence, `GU572157.1`, has been removed from the single tree example database. To obtain assignments, run `tronko-assign` with the single tree example database from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment and a default LCA cut-off of 5):

```
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta -o example_datasets/single_tree/missingreads_singleend_150bp_2error_results.txt -w
```

To obtain assignments, run `tronko-assign` with the single tree example database from `tronko-build` with the following command for the paired-end reads (using the Needleman-Wunsch alignment and a default LCA cut-off of 5):

```
tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta -p -1 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/single_tree/missingreads_pairedend_150bp_2error_results.txt -w
```

An example dataset for `tronko-assign` with multiple trees is provided in `example_datasets/multiple_trees`. This dataset contains single-end reads (`example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta`) and paired-end reads (forward read: `example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta` and reverse read: `example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta`). For the single-end reads, there are 84 150bp reads with 2% simulated error/polymorphisms from the following sequence (taxonomically classified on NCBI as Sagitta elegans):

```
>KP857443.1
TTTGAGCACTGTGGGACATAGGGGTGGAGCAGTGGATTTGGGTATCTTCTCTTTGCACCTGGCTGGCGTTAGAAGAATCTTGGGGAGAGCTAATTTTATTACCACTATCACCAATATAAAAGGGGAAGGTATGACTATAGAACTCATGCCTTTATTCGTGTGGGCGGTGCTCCTCACGGCTGTCTTACTTTTACTCTCTCTACCTGTATTAGCTGGGGCTATCACAATGTTAC
```

This particular sequence, `KP857443.1`, has been removed from the multiple trees example databases. To obtain assignments, run `tronko-assign` with the multiple trees example database with partitions from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment `-w` and a default LCA cut-off of 5 `-c 5`):

```
tronko-assign -r -f out_oneMSA/reference_tree.txt -s -g example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta -o example_datasets/multiple_trees/missingreads_singleend_150bp_2error_partition_results.txt -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To assign paired-end reads on the same reference database run `tronko-assign`:
```
tronko-assign -r -f out_oneMSA/reference_tree.txt -p -1 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_partition_results.txt -a -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To obtain assignments, run `tronko-assign` with the multiple trees example database with multiple MSAs from `tronko-build` with the following command for the single-end reads (using the Needleman-Wunsch alignment `-w` and a default LCA cut-off of 5 `-c 5`):

```
tronko-assign -r -f outdir_multiple_MSA/reference_tree.txt -s -g example_datasets/multiple_trees/missingreads_singleend_150bp_2error.fasta -o example_datasets/multiple_trees/missingreads_singleend_150bp_2error_multiple_MSA_results.txt -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

To assign paired-end reads on the same reference database run `tronko-assign`:
```
tronko-assign -r -f outdir_multiple_MSA/reference_tree.txt -p -1 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta -2 example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta -o example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_multiple_MSA_results.txt -a tronko-build/example_datasets/multiple_trees/one_MSA/1.fasta
```

## More on `tronko-build` Usage with multiple trees

`tronko-build` requires a multiple sequence alignment (FASTA format), rooted phylogenetic tree (Newick format), and a corresponding taxonomy file for each cluster build. All of the files should be in one directory and specify the directory with `-e` with each cluster being designated by a number. MSA files should be named `[Number]_MSA.fasta`, taxonomy files should be named `[Number]_taxonomy.txt`, and tree files should be named `RAxML_bestTree.[Number].reroot`. Example of the contents of a directory containing 3 clusters:
```
1_MSA.fasta
2_MSA.fasta
3_MSA.fasta
1_taxonomy.txt
2_taxonomy.txt
3_taxonomy.txt
RAxML_bestTree.1.reroot
RAxML_bestTree.2.reroot
RAxML_bestTree.3.reroot
```
Once you have the cluster files prepared and you have created an output directory (the program assumes the output directory already exists), an example command is
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n [NUMBER OF PARTITIONS] -d [OUTPUT DIRECTORY] -s
```
For the example with 3 clusters, an example command partitioning by sum-of-pairs score would be:
```
tronko-build -y -e [DIRECTORY CONTAINING MSA, TAX, and TREE FILES] -n 3 -d output -s
```
The reference database file will be output to `[OUTPUT DIRECTORY]/reference_tree.txt`. The `reference_tree.txt` file is the reference database file that `tronko-assign` requires for assignment.

To partition the reference database further (see manuscript for details) the following dependencies are needed to be installed in your path (no dependencies needed for `tronko-assign` and for building a reference database with a single tree or without partitioning the database): <a href="https://github.com/stamatak/standard-RAxML">`raxmlHPC-PTHREADS`</a>, <a href="https://github.com/refresh-bio/FAMSA">`famsa`</a>, `nw_reroot` from <a href="https://anaconda.org/bioconda/newick_utils/files">Newick utilties</a>, <a href="https://raw.githubusercontent.com/lpipes/tronko/main/scripts/fasta2phyml.pl">`fasta2phyml.pl`</a>, and <a href="https://ftp.gnu.org/gnu/sed/">`sed`</a>. Partitioning the database further is only needed when the underlying MSA is unreliable. An example command to create the reference database and partition a database that contains 100 initial clusters using the sum-of-squares approach:

```
tronko-build -y -e initial_clusters_directory -d outdir -n 100 -s
```

The `reference_tree.txt` file will be output to the `outdir` directory. An example command to create the reference database and partition a database that contains 100 initial clusters using a threshold for the number of leaf nodes (i.e., 500 leaf nodes):

```
tronko-build -y -e initial_clusters_directory -d outdir -n 100 -v -f 500
```

The `reference_tree.txt` file will be output to the `outdir` directory.

# `tronko-assign` using pre-built 16S and COI databases

First download the databases from Zenodo (<a href="https://doi.org/10.5281/zenodo.7407318"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7407318.svg" alt="DOI"></a>). To assign FASTQ (`-q`) paired-end reads (`-p`) using the 16S database (and reverse-complement the reverse read `-z`) with Needleman-Wunsch (`-w`), 16 cores (`-C 16`), and an LCA cut-off of 5 (`-c 5`):

```
tronko-assign -r -q -p -z -w -C 16 -c 5 -f 16S_tronko_build.txt.gz -a 16S.fasta -1 16S_TW-DR-1-S88_F_filt.fastq.gz -2 16S_TW-DR-1-S88_R_filt.fastq.gz -o 16S_TW-DR-1-S88_results.txt 
```  

With the CO1 database and same parameters:
```
tronko-assign -r -q -p -z -w -C 16 -c 5 -f CO1_tronko_build.txt.gz -a CO1.fasta -1 CO1_TW-DR-1-S88_F_filt.fastq.gz -2 CO1_TW-DR-1-S88_R_filt.fastq.gz -o CO1_TW-DR-1-S88_results.txt
```

# Testing

Tronko includes a comprehensive testing suite to ensure reliability and help with development.

## Quick Testing

### Unit Tests
Test core functionality including crash debugging, signal handling, and corruption detection:
```bash
cd tests
make -f Makefile.simple smoke
```

### Integration Tests  
Test binary functionality, command-line interfaces, and data validation:
```bash
cd tests
make -f Makefile.simple test_simple_workflows
./test_simple_workflows
```

### Functional Tests
Test end-to-end workflows with real datasets:
```bash
./test_with_example_data.sh
```

## Comprehensive Testing

### Full Test Suite
Run all tests (unit, integration, and functional):
```bash
cd tests
make -f Makefile.simple test
```

### Using the Test Runner
The test runner provides convenient options for different testing scenarios:
```bash
# Run all tests
./run_tests.sh

# Run only unit tests  
./run_tests.sh --unit-only

# Run only integration tests
./run_tests.sh --integration-only

# Generate code coverage report
./run_tests.sh --coverage

# Run with memory leak detection
./run_tests.sh --valgrind

# Verbose output for debugging
./run_tests.sh --verbose
```

## Docker Testing

Test in a consistent containerized environment:
```bash
# Build and run tests in Docker
docker compose up -d
docker compose exec tronko-dev bash -c "cd /app/tests && make -f Makefile.simple test"

# Or use the test runner in Docker
docker compose exec tronko-dev bash -c "cd /app && ./run_tests.sh"
```

## Continuous Integration

Tests run automatically on every push via GitHub Actions:
- **Quick Tests**: Fast validation on all branches (~2-3 minutes)
- **Comprehensive Tests**: Full validation on main/experimental/develop branches (~10-15 minutes)
- **Performance Tests**: Benchmarking on main/experimental branches

View test status: [![Quick Tests](https://github.com/lpipes/tronko/workflows/Quick%20Tests/badge.svg)](https://github.com/lpipes/tronko/actions)

## Test Categories

### Unit Tests
- **Framework**: Unity (lightweight C testing framework)
- **Coverage**: Crash debugging system, corruption detection, signal handling
- **Location**: `tests/unit/`

### Integration Tests
- **Coverage**: Binary existence, command-line interfaces, error handling
- **Location**: `tests/integration/`

### Functional Tests  
- **Coverage**: End-to-end workflows, real data processing, crash testing
- **Location**: `test_with_example_data.sh`

For detailed testing documentation, see [`tests/README.md`](tests/README.md).

# Reproducibility and Variance

`tronko-assign` uses BWA's multi-threaded alignment engine, which introduces non-deterministic ordering of candidate matches. This results in approximately **3% variance** in taxonomic assignments between runs with identical inputs when using multiple cores.

## Understanding the Variance

- **Magnitude**: ~3% of reads may receive different taxonomic assignments across runs
- **Nature**: Differences are unbiased (equal probability of more/less specific results)
- **Cause**: BWA's work-stealing thread scheduler finds matches in variable order
- **Scientific Impact**: Variance is typically smaller than biological variation in most studies

### What This Means

The variance affects individual read assignments, not overall community composition. For most metabarcoding studies:
- ✅ Community-level statistics (abundance, diversity) remain stable
- ✅ Biological replicates average out technical variance
- ✅ Differential abundance analysis is unaffected when using proper statistical methods
- ⚠️ Exact read-by-read reproducibility requires special modes (see below)

## Achieving Reproducible Results

### Option 1: Single-Threaded Mode (100% Deterministic)

Use the `-C 1` flag to run in deterministic single-threaded mode:

```bash
tronko-assign -C 1 -r -f reference.trkb -a reference.fasta -p \
  -1 forward.fasta -2 reverse.fasta -o results.txt
```

**Trade-offs**:
- ✅ Perfect reproducibility (identical results every run)
- ❌ 10-16x slower than multi-threaded mode
- 👍 **Use for**: Validation, testing, regression checks, publications requiring exact reproducibility

### Option 2: Standard Multi-Threaded Mode (Fast, ~3% Variance)

Use multiple cores for performance (default or explicit `-C N`):

```bash
tronko-assign -C 16 -r -f reference.trkb -a reference.fasta -p \
  -1 forward.fasta -2 reverse.fasta -o results.txt
```

**Trade-offs**:
- ✅ Fast execution (scales with CPU cores)
- ⚠️ ~3% variance between runs
- 👍 **Use for**: Standard workflows, exploratory analysis, production pipelines

### Option 3: Consensus Voting (Best of Both Worlds)

Run multiple times and aggregate results for both speed and reproducibility:

```bash
# Run 3 times in parallel
tronko-assign -C 16 [...] -o run1.txt &
tronko-assign -C 16 [...] -o run2.txt &
tronko-assign -C 16 [...] -o run3.txt &
wait

# Aggregate results (see consensus script in scripts/)
# (Consensus voting script to be implemented)
```

**Trade-offs**:
- ✅ Reduces variance from ~3% to ~0.3-0.6%
- ✅ Can parallelize across runs (3 runs ≈ 3x cost but parallelizable)
- 👍 **Use for**: Critical datasets, regulatory submissions, important publications

## Testing for Variance

Measure variance on your specific dataset using the test script:

```bash
cd tronko-assign

./scripts/test_determinism.sh \
  -r /path/to/reference.trkb \
  -a /path/to/reference.fasta \
  -1 /path/to/forward.fasta \
  -2 /path/to/reverse.fasta \
  -n 3 \
  -c 16 \
  -k
```

This will run 3 replicates and report variance statistics.

### Expected Results

- **Single-threaded** (`-c 1`): 0% variance (identical results)
- **Multi-threaded** (e.g., `-c 16`): 2-4% variance depending on dataset and core count
- **Higher core counts**: Slightly higher variance due to increased thread contention

## Recommendations by Use Case

| Use Case | Recommended Mode | Cores | Expected Variance |
|----------|------------------|-------|-------------------|
| Validation / Testing | Single-threaded | `-C 1` | 0% |
| Exploratory Analysis | Multi-threaded | `-C 4-16` | ~3% |
| Production Pipeline | Multi-threaded | `-C 8-16` | ~3% |
| Publication (standard) | Multi-threaded | `-C 4-16` | ~3% |
| Publication (exact reproducibility) | Single-threaded | `-C 1` | 0% |
| Regulatory Submission | Consensus voting | `-C 16`, 3 runs | <0.5% |

## Technical Details

For developers and advanced users:

- **Root cause**: BWA's `kt_for()` work-stealing scheduler in `bwa_source_files/kthread.c`
- **Propagation**: Variable match order → different tie-breaking in `placement.c:893`
- **Floating-point sensitivity**: Confidence interval comparisons (`placement.c:930`) sensitive to small score differences
- **Not affected by**: WFA2 alignment (deterministic), scoring algorithm (deterministic), result collection order (deterministic)

For detailed technical analysis, see:
- Research document: `thoughts/shared/research/2026-01-03-non-determinism-mitigation-strategies.md`
- Baseline measurement: `thoughts/shared/research/2026-01-03-non-determinism-baseline-measurement.md`

## FAQ

**Q: Is this a bug?**
A: No, this is inherent to BWA's multi-threaded design. It's a performance/reproducibility trade-off present in many bioinformatics tools.

**Q: Will this affect my biological conclusions?**
A: For most metabarcoding studies, no. The variance is small (~3%) and unbiased, similar to technical replicates.

**Q: Can I use tronko-assign in CI/CD pipelines?**
A: Yes, but use statistical comparison (e.g., expect <5% difference) rather than exact output matching. Or use `-C 1` for exact regression tests.

**Q: Does this affect accuracy?**
A: No, accuracy is identical. The algorithm is sound; only the specific reads assigned to borderline taxa vary slightly.

**Q: Will this be fixed in future versions?**
A: We are exploring options including BWA seed control and alternative aligners. For now, single-threaded mode provides perfect reproducibility when needed.

## Advanced Features

Tronko includes comprehensive logging and debugging capabilities:

### Performance Monitoring
```bash
# Enable performance logging with resource monitoring
tronko-assign -V2 -R -T -l performance.log [options...]
```

### Crash Debugging
Automatic crash detection with root cause analysis:
```bash
# Crash debugging is always enabled
# View crash reports in /tmp/tronko_assign_crash_*.crash
tronko-assign -V2 [options...]  # Enhanced context capture
```

### Documentation
- **[Performance Logging Guide](docs/performance-logging.md)**: Comprehensive performance monitoring, resource tracking, and optimization
- **[Crash Debugging Guide](docs/crash-debugging.md)**: Advanced crash detection, root cause analysis, and debugging workflows
- **[Complete Documentation Index](docs/index.md)**: Full documentation overview and quick reference
- **[Experiments Log](EXPERIMENTS_LOG.md)**: Benchmark results and optimization experiments (memory vs accuracy tradeoffs)

# tronko-assign Performance Optimizations

This section documents performance optimizations applied to `tronko-assign` and the methodology used to validate them. All optimizations produce **byte-identical output** to the original code — no accuracy or behavioral changes.

## Methodology

### Golden Output Testing

Optimizations are validated using a golden output regression framework. The original (unmodified) binary is run on multiple test configurations and the outputs are saved as golden reference files. After each optimization, the modified binary is run on the same inputs and the outputs are compared byte-for-byte.

**Test configurations** (all use `-C 1` for deterministic output):

| Test | Mode | Aligner | Reads |
|------|------|---------|-------|
| `single_tree_singleend_nw` | single-end | Needleman-Wunsch | 164 |
| `single_tree_singleend_wfa` | single-end | WFA | 164 |
| `single_tree_pairedend_nw` | paired-end | NW | 9 pairs |
| `single_tree_pairedend_wfa` | paired-end | WFA | 9 pairs |
| `benchmark_10k_singleend_nw` | single-end | NW | 10,000 |
| `benchmark_10k_singleend_wfa` | single-end | WFA | 10,000 |

The 10K-read benchmarks use a synthetic FASTA generated by replicating the 164 example reads. At this scale, per-read computation dominates total runtime (~1.76ms/read) with negligible startup overhead.

**Test scripts:**
- `tests/test_golden_outputs.sh [binary]` — runs all 6 configurations and diffs against golden files
- `tests/benchmark_optimization.sh [binary] [iterations]` — 3-run timing with correctness verification

### Benchmark Environment

- Dataset: Charadriiformes single tree (1,466 species, ~2,931 nodes), 10,000 single-end reads
- Hardware: Apple M-series, single thread (`-C 1`)
- Each measurement: average of 3 runs

## Results

### Combined Speedup

| Aligner | Before | After | Speedup | Per-read |
|---------|--------|-------|---------|----------|
| Needleman-Wunsch | 17.78s | 12.78s | **28.1% faster** | 1.78ms → 1.28ms |
| WFA | 14.54s | 9.48s | **34.8% faster** | 1.45ms → 0.95ms |

Memory usage unchanged (68.9 MB / 71.2 MB RSS). All 6 golden output tests pass.

### Optimization Details

#### 1. Scoring hot loop (`assignment.c: getscore_Arr`)

The inner scoring function is called for every node in every matched tree for every read (~29,310 calls per read). Three changes:

- **Precomputed log constants**: `log(0.01)` and `log(0.25)` were recomputed on every loop iteration. Replaced with `static const` values.
- **Lookup table for nucleotide dispatch**: The if/else chain mapping `A/a/C/c/G/g/T/t/-` to nucleotide indices was replaced with a 256-entry lookup table (`char → int`), eliminating branches in the inner loop.
- **Split fast/slow path on debug flag**: The `print_all_nodes` flag (almost always 0) was checked on every iteration. Restructured into two code paths so the common case has zero debug overhead.

This single optimization delivered the largest improvement (~24-28%).

#### 2. Output string formatting (`tronko-assign.c`)

Per-read output was built with 5 separate `asprintf()` heap allocations, each followed by `strcat()` and `free()`. Replaced with cursor-based `snprintf()` into the pre-allocated buffer — zero heap allocations in the output path.

#### 3. WFA aligner reuse (`placement.c`)

A `wavefront_aligner_t` was created and destroyed for every BWA match (up to 10 per read, doubled for paired-end). Moved aligner creation to once per `place_paired()` call and reuse across all matches. Per-match text buffers are still allocated/freed as required by the WFA2 API.

#### 4. Iterative tree traversal (`assignment.c: assignScores_Arr_paired`)

The recursive DFS traversal passed 17 parameters per call across ~2,931 stack frames. Replaced with an iterative traversal using an explicit stack array. The DFS visit order is preserved (push child1 then child0) so output is identical.

#### 5. Buffer clearing with memset (`placement.c`, `tronko-assign.c`)

Per-element clearing loops for `positions[]`, `locQuery[]`, `leaf_sequence[]`, `positionsInRoot[]`, `voteRoot[][]`, `nodeScores[][][]`, and `minNodes[]` were replaced with `memset()` calls.

## Build Flags

The Makefile was also updated to support building on macOS with the crash debugging system:

```makefile
# macOS requires _XOPEN_SOURCE for ucontext.h, plus _DARWIN_C_SOURCE to keep dladdr/Dl_info
ifeq ($(UNAME_S),Darwin)
    XOPEN_FLAGS = -D_XOPEN_SOURCE=600 -D_DARWIN_C_SOURCE
endif
```

## Reproducing

```bash
# Build
cd tronko-assign && make clean && make

# Verify correctness (byte-identical output to baseline)
bash tests/test_golden_outputs.sh

# Measure performance
bash tests/benchmark_optimization.sh
```

# ancestralclust Performance Optimizations

This section documents performance optimizations applied to `ancestralclust` (the initial clustering step used by `tronko-build`). All optimizations produce **byte-identical output** to the original code — no accuracy or behavioral changes.

## Methodology

### Golden Output Testing

Optimizations are validated using a golden output regression framework. The original binary is run with deterministic settings (`-c 1`, NW mode) and outputs are saved as reference files. After each optimization, the modified binary is run on the same inputs and the cluster FASTA files are compared byte-for-byte.

**Test configuration**: 1,466 Charadriiformes COI sequences, `-b 50 -r 100 -p 10 -c 1`, NW alignment mode. 14 output cluster files checked for sequence count and byte-identical content (28 checks total).

**Test script**: `ancestralclust/test_consistency.sh`

## Results

### Benchmark: 10,000 CO1 sequences, 6 threads, WFA mode, kseq=500

| Phase | Baseline | Optimized | Speedup |
|-------|----------|-----------|---------|
| Iter 1: Distance matrix | 86.4s | 85.2s | 1.3% |
| Iter 1: Assignment (9,500 seqs) | 22.6s | 21.5s | **4.7%** |
| Iter 2: Distance matrix | 94.9s | 94.2s | 0.8% |
| **Total** | ~209s | ~205s | **~2%** |

### Benchmark: 1,466 Charadriiformes, 6 threads, WFA mode, kseq=100

| Phase | Baseline | Optimized | Speedup |
|-------|----------|-----------|---------|
| Distance matrix | 3.72s | 3.45s | 8% |
| Tree + clustering | 4.07s | 3.78s | 8% |
| Assignment | 0.44s | 0.41s | 8% |

### Benchmark: 1,466 Charadriiformes, 1 thread, WFA mode, kseq=100

| Phase | Baseline | Optimized | Speedup |
|-------|----------|-----------|---------|
| Distance matrix | 0.180s | 0.090s | **2x** |

At scale (kseq=500), the distance matrix involves ~125,000 WFA pairwise alignments. Runtime is dominated by `wavefront_align()` inner loops — WFA2-lib is already heavily optimized (SIMD, O(ns) complexity). The structural optimizations reduce overhead around the WFA calls but cannot improve WFA itself.

### Golden Output Verification

All benchmarks verified with 28/28 golden output checks passing (byte-identical cluster files).

## Optimization Details

### 1. Nucleotide lookup table (`ancestralclust.c: encode_base`)

The `switch(c)` nucleotide encoder in the inner CIGAR-to-JC distance loop was replaced with a 256-entry static lookup table. Eliminates branch mispredictions in the hot loop of `compute_JC_from_CIGAR`.

### 2. Flat contiguous distance matrix

Two allocation sites changed from N separate `calloc()` calls (one per row) to a single `calloc(dim*dim)` with row-pointer setup. Improves cache locality for the pairwise distance matrix access pattern and reduces allocator overhead.

### 3. Balanced triangular thread partitioning

The upper-triangle distance matrix computation was partitioned by column ranges, causing thread imbalance (column j requires j comparisons). Replaced with work-balanced partitioning using `sqrt(2 * k * total_work / T)` cutpoints so each thread gets equal work.

### 4. Bubble sort → qsort

Branch length sorting used O(N²) bubble sort with parallel array swaps. Replaced with `qsort()` using a struct-of-pairs pattern for O(N log N) sorting.

### 5. O(N) tree diameter (`findLongestTipToTip`)

The original implementation found the longest tip-to-tip path using O(N²) all-pairs LCA distance comparisons. Replaced with a single O(N) iterative post-order DFS that tracks the deepest leaf in each subtree and finds the diameter through each internal node.

### 6. ML estimation malloc elimination

`estimatenucparameters` calls `estimatebranchlengths` (6×) and `getlike_gamma` (indirectly), each of which allocated and freed `templike_nc`, `UFCnc`, and `statevector` arrays. Hoisted allocation to the outer function, passed pre-allocated buffers through the call chain, and freed once at the end.

Additionally, `getlike_gamma` flattened `locloglike` from `double[numbase][NUMCAT]` to `double[numbase]` since `NUMCAT==1`.

### 7. strlen caching in `fillInMat_avg`

Pre-compute and cache `strlen()` results for cluster sequences in VLA arrays before the inner comparison loops, avoiding redundant length computations.

### 8. Build flags: NATIVE_ARCH and LTO

The Makefile supports optional architecture-specific tuning and link-time optimization:

```bash
# Standard build
make

# With native CPU tuning and link-time optimization
make NATIVE_ARCH=1 LTO=1
```

`NATIVE_ARCH=1` enables `-mcpu=native` on ARM64 or `-march=native -mtune=native` on x86. `LTO=1` enables `-flto` for cross-file inlining (particularly beneficial for WFA2 library calls).

## Reproducing

```bash
# Build
cd ancestralclust && make clean && make

# Verify correctness (28 byte-identical checks)
bash test_consistency.sh

# Benchmark (requires input FASTA)
bash benchmark.sh <input.fasta> [threads] [label]
```

# Performance

We performed a leave-one-species-out test comparing Tronko (with LCA cut-offs for the score of 0, 5, 10, 15, and 20 with Needleman-Wunsch alignment) to kraken2, metaphlan2, and MEGAN for 1,467 COI sequences from 253 species from the order Charadriiformes using 150bp x 2 paired-end sequences and 150bp and 300bp single-end sequences using 0, 1, and 2% error/polymorphism.
<img src="https://github.com/lpipes/tronko/blob/main/LSO.png?raw=true">
Using leave-one-species-out and simulating reads (both paired-end and single-end) with a 0-2% error (or polymorphism), Tronko detected the correct genus more accurately than the other methods even when using an aggressive cut-off (i.e., when cut-off=0) (D and G).

# Citation

Pipes L, and Nielsen R (2024) A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. eLife.
https://elifesciences.org/articles/85794.pdf
