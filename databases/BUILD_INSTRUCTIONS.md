# How to Build a Tronko Database

## Prerequisites

```bash
# Tools on PATH
tronko-build    # cd tronko-build && make
FastTree        # brew install fasttree
famsa           # brew install famsa
nw_reroot       # brew install newick-utils
bwa             # brew install bwa
tronko-convert  # cd tronko-convert && make (optional, for .trkb format)

# PASTA (for tree estimation)
git clone https://github.com/smirarab/pasta.git ~/pasta
python3 -m venv /tmp/pasta_env
source /tmp/pasta_env/bin/activate
pip install dendropy
export PYTHONPATH=~/pasta
```

## Step 1: Run PASTA

Produces a co-estimated alignment + tree from raw sequences. This is the most expensive step (~8-12 hours for 100K sequences).

```bash
source /tmp/pasta_env/bin/activate
export PYTHONPATH=~/pasta

python3 ~/pasta/run_pasta.py \
    --input input.fasta \
    --datatype dna \
    --iter-limit 3 \
    --max-subproblem-size 500 \
    --num-cpus 12 \
    --merger muscle \
    --job my_marker \
    --output-directory pasta_output/ \
    -o pasta_output/
```

**Key parameters:**
- `--max-subproblem-size`: Size of chunks MAFFT aligns at once. Larger = better alignment, slower. Tested: 200 (original), 500 (3.5% better scores).
- `--iter-limit`: PASTA iterations. 3 is sufficient — scores plateau after iteration 2.
- `--merger muscle`: Use muscle for sub-alignment merging. Default (opal) may not be installed.

**Output:** `pasta_output/my_marker.tre` — Newick tree of all input sequences.

## Step 2: Fix the PASTA tree

Strip single quotes that PASTA adds to leaf names:

```bash
sed "s/'//g" pasta_output/my_marker.tre > pasta_output/my_marker_fixed.tre
```

## Step 3: Decompose + per-partition alignment + tree building

Splits the PASTA tree into manageable partitions, aligns each with FAMSA, builds per-partition trees with FastTree.

```bash
python3 tronko-build/partition_and_build.py \
    --tree pasta_output/my_marker_fixed.tre \
    --fasta input.fasta \
    --taxonomy input_taxonomy.txt \
    --outdir partitions/ \
    --max-size 1000 \
    --min-size 3 \
    --strategy centroid \
    --threads 12 \
    --skip-tronko
```

**Key parameters:**
- `--max-size N`: Max leaves per partition. Controls the tradeoff between partition count and tree size. Values tested: 500, 1000, 2000, 5000.
- `--max-diam D`: Alternative to max-size — limits tree diameter per partition. Values tested: 10, 25.
- `--min-size 3`: Minimum leaves per partition (need at least 3 for a tree).
- `--strategy centroid`: Decomposition strategy. Centroid balances partition sizes.
- `--skip-tronko`: Only do decomposition + alignment, run tronko-build separately (required for norepartition).

**Output:** `partitions/` directory with per-partition MSA, taxonomy, and tree files.

## Step 4: Run tronko-build

There are two modes depending on whether you want SP-score partitioning:

### Option A: No SP partitioning (PASTA builds)

Use this when PASTA's tree decomposition already produced good partitions and you don't want further splitting. Benchmarks showed this outperforms SP partitioning for PASTA-based builds.

```bash
N_PARTITIONS=$(ls partitions/RAxML_bestTree.*.reroot | wc -l | tr -d ' ')

tronko-build -y \
    -e partitions/ \
    -n $N_PARTITIONS \
    -d my_database/ \
    -E -a -f 999999 \
    -c 12
```

### Option B: With SP partitioning (AncestralClust builds)

Use this when partitions come from AncestralClust clustering and need refinement via SP-score recursive splitting.

```bash
tronko-build -y \
    -e partitions/ \
    -n $N_PARTITIONS \
    -d my_database/ \
    -E -a -s -u 0.1 \
    -c 12
```

**Flag reference:**
- `-y`: Partition mode (reads multiple trees)
- `-e partitions/`: Directory containing partition files
- `-n $N_PARTITIONS`: Number of partitions
- `-d my_database/`: Output directory
- `-E`: Export subtrees for ablation studies
- `-a`: Use FastTree-style Newick parser (required for FastTree/FAMSA-built trees)
- `-c 12`: Threads for posterior computation
- **`-f 999999`**: **(Option A only)** Prevents tronko-build's internal 3-way variance-based repartitioning. Without this, every input partition is split into 3 subtrees regardless of other flags. Only use this when your partitions are already well-formed (e.g. from PASTA decomposition).
- **`-s -u 0.1`**: **(Option B only)** Enables SP-score recursive partitioning with threshold 0.1. Trees with SP score below the threshold are further split.

**Output:** `my_database/reference_tree.txt` + `exported_subtrees/` + `tree_list.txt`

## Step 5: Build marker.fasta + BWA index

Concatenate gap-free sequences from all partitions:

```bash
for i in $(seq 0 $((N_PARTITIONS - 1))); do
    if [ -f partitions/partition${i}_unaligned.fasta ]; then
        cat partitions/partition${i}_unaligned.fasta
    elif [ -f partitions/partition${i}_MSA.fasta ]; then
        awk '/^>/{print; next} {gsub(/-/,""); print}' partitions/partition${i}_MSA.fasta
    fi
done > my_database/marker.fasta

bwa index my_database/marker.fasta
```

## Step 6: Build marker_taxonomy.txt

```bash
for i in $(seq 0 $((N_PARTITIONS - 1))); do
    [ -f partitions/partition${i}_taxonomy.txt ] && cat partitions/partition${i}_taxonomy.txt
done > my_database/marker_taxonomy.txt
```

## Step 7: Convert to .trkb (optional)

Zstd-compressed binary format — 30x smaller, used by tronko-assign:

```bash
tronko-convert -i my_database/reference_tree.txt \
    -o my_database/reference_tree.trkb -c zstd
```

## Step 8: Copy partition files into database directory

The ablation benchmark system expects partition files alongside the database:

```bash
cp partitions/partition*_MSA.fasta my_database/
cp partitions/partition*_taxonomy.txt my_database/
cp partitions/RAxML_bestTree.*.reroot my_database/
```

## Step 9: Copy input files for provenance

```bash
cp input.fasta my_database/input.fasta
cp input_taxonomy.txt my_database/input_taxonomy.txt
```

## Step 10: Write metadata

### build_info.json

```bash
N_TREES=$(head -1 my_database/reference_tree.txt)
SIZE_MB=$(du -m my_database/reference_tree.txt | cut -f1)

cat > my_database/build_info.json << EOF
{
  "name": "my_database_name",
  "pipeline": "PASTA + partition_and_build.py + tronko-build (norepartition)",
  "marker": "12SV5",
  "input_sequences": $(grep -c '^>' input.fasta),
  "n_trees": $N_TREES,
  "n_input_partitions": $N_PARTITIONS,
  "db_size_mb": $SIZE_MB,
  "sp_partitioning": false,
  "internal_repartitioning": "none",
  "gamma": false,
  "decomposition": {
    "mode": "max_size",
    "value": 1000,
    "strategy": "centroid",
    "min_size": 3
  },
  "tronko_build_flags": "-y -e <parts> -n $N_PARTITIONS -d <out> -E -a -f 999999 -c 12",
  "date": "$(date +%Y-%m-%d)"
}
EOF
```

### pasta_stats.json

```bash
# Extract from PASTA output
SCORE=$(cat pasta_output/my_marker.score.txt)
cat > my_database/pasta_stats.json << EOF
{
  "pasta_version": "1.9.3",
  "iter_limit": 3,
  "max_subproblem_size": 500,
  "aligner": "mafft",
  "merger": "muscle",
  "tree_estimator": "fasttree",
  "num_cpus": 12,
  "input_sequences": $(grep -c '^>' input.fasta),
  "final_score": $SCORE
}
EOF
```

## Step 11: Verify

A complete database directory should contain:

```
my_database/
├── reference_tree.txt          # tronko database
├── reference_tree.trkb         # compressed binary (optional)
├── marker.fasta                # concatenated sequences for BWA
├── marker.fasta.{amb,ann,bwt,pac,sa}  # BWA index
├── marker_taxonomy.txt         # taxonomy for all sequences
├── tree_list.txt               # partition manifest (from tronko-build)
├── exported_subtrees/          # per-partition files (from -E flag)
├── partition*_MSA.fasta        # per-partition MSA (for ablation)
├── partition*_taxonomy.txt     # per-partition taxonomy (for ablation)
├── RAxML_bestTree.*.reroot     # per-partition trees (for ablation)
├── build_info.json             # build parameters
├── pasta_stats.json            # PASTA scores and config
├── input.fasta                 # original input sequences
└── input_taxonomy.txt          # original input taxonomy
```

## Using with tronko-assign

```bash
tronko-assign -r \
    -f my_database/reference_tree.txt \
    -a my_database/marker.fasta \
    -s -g query_reads.fasta \
    -o results.txt -w
```

## Using with benchmark sweep

```bash
poetry run python notebooks/tronko_db_sweep.py --db-only \
    --db-sweep-root /path/to/databases
```

## Optional: FastTree -gamma branch length reoptimization

After Step 3, before Step 4, reoptimize branch lengths under a gamma model:

```bash
for tree in partitions/RAxML_bestTree.*.reroot; do
    name=$(basename "$tree" | sed 's/RAxML_bestTree.//' | sed 's/.reroot//')
    msa="partitions/${name}_MSA.fasta"
    OMP_NUM_THREADS=1 FastTree -nt -gtr -gamma -nome -mllen \
        -intree "$tree" "$msa" > "${tree}.gamma" 2>/dev/null
    mv "${tree}.gamma" "$tree"
done
```

This doesn't change tree topology, only branch lengths. Better branch lengths may improve tronko-build's posterior probability computation.
