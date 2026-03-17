# Tronko Reference Databases

All tronko reference databases built during this project. Tracked with DVC for versioning large files.

**Marker gene**: 12SV5 (vertebrate 12S rRNA)
**Input sequences**: 101,088 (filtered) or 176,947 (unfiltered) vertebrate amplicon references from CRUXv2

---

## Database Inventory

### 1. `charadriiformes_single_tree/`

Small test dataset used for development and correctness verification.

| Property | Value |
|---|---|
| Marker | COI (Cytochrome Oxidase I) |
| Sequences | 1,413 |
| Alignment columns | 317 bp |
| Trees | 1 (single-tree mode) |
| DB size | 38 MB |
| Build mode | `tronko-build -l` (single tree, no partitioning) |

**Input**: Charadriiformes (shorebirds) from CRUXv2 COI database.

**Build command**:
```bash
tronko-build -l \
    -m Charadriiformes_MSA.fasta \
    -x Charadriiformes_taxonomy.txt \
    -t RAxML_bestTree.Charadriiformes.reroot \
    -d .
```

---

### 2. `vert12S_fasttree/`

vert12S database built with AncestralClust clustering + FastTree tree inference.

| Property | Value |
|---|---|
| Marker | 12SV5 |
| Input sequences | 101,088 |
| Tree inference | FastTree 2 |
| Trees | varies by SP threshold |
| DB size | 1.4 GB (text), 214 MB (gzipped), 68 MB (trkb/zstd) |
| Build time | ~45 min |

**Build command**:
```bash
./build-tronko-db.sh \
    -f 12SV5_species.fasta \
    -t 12SV5_species_taxonomy.txt \
    -o output_fasttree/ \
    -T 12 -s 0.1 -F
```

**Notes**: Contains all three formats (`.txt`, `.txt.gz`, `.trkb`) for format comparison experiments. The `.txt.bak` is a pre-patch_taxonomy version (left in original location).

---

### 3. `vert12S_veryfasttree/`

vert12S database built with AncestralClust clustering + VeryFastTree (parallelized FastTree).

| Property | Value |
|---|---|
| Marker | 12SV5 |
| Input sequences | 101,088 |
| Tree inference | VeryFastTree (multi-threaded) |
| DB size | 217 MB (gzipped) |

**Build command**:
```bash
./build-tronko-db.sh \
    -f 12SV5_species.fasta \
    -t 12SV5_species_taxonomy.txt \
    -o output_veryfasttree/ \
    -T 12 -s 0.1 -a    # -a enables VeryFastTree
```

---

### 4. `vert12S_pasta/`

vert12S database built using PASTA (SATe-enabled Phylogenetic Alignment and Tree Estimation) for simultaneous alignment and tree building, with SP-score partitioning.

| Property | Value |
|---|---|
| Marker | 12SV5 |
| Input sequences | 101,088 |
| Tree inference | PASTA (iterative alignment + FastTree) |
| DB size | 2.4 GB (text) |

**Build command**:
```bash
python3 partition_and_build.py \
    --tree pasta_output/vert12S_pasta1.fixed.tre \
    --fasta pasta_output/12SV5_species.fixed.fasta \
    --taxonomy pasta_output/12SV5_species_taxonomy.fixed.txt \
    --outdir pasta_sp_partitions \
    --max-diam 25.0 --min-size 3 \
    --threads 12 --sp-threshold 0.1 \
    --tronko-build tronko-build/tronko-build \
    --tronko-outdir pasta_tronko_db
```

**Notes**: Uses PASTA's guide tree instead of AncestralClust for clustering. PASTA produces a single large tree, which is then partitioned by SP-score threshold.

---

### 5. `12sv5_unfiltered/`

Raw unfiltered 12SV5 dataset (176,947 sequences) with a small test build.

| Property | Value |
|---|---|
| Marker | 12SV5 |
| Input sequences | 176,947 |
| Status | Partial build (3 partitions only) |

Contains the raw `12SV5_species.fasta` and `12SV5_species_taxonomy.txt` input files along with a few test partitions.

---

### 6-17. `sweep_*` directories

Parameter sweep across 4 AncestralClust clustering configs x 3 SP-score thresholds = 12 databases. All built from the same 101,088 filtered vert12S sequences.

#### Clustering Configurations

| Config | Description |
|---|---|
| `ac_default` | AncestralClust default parameters |
| `ac_more_bins` | More clustering bins (finer initial clusters) |
| `ac_fewer_bins` | Fewer clustering bins (coarser initial clusters, larger trees) |
| `ac_single` | Single cluster (no AncestralClust; one global FastTree, then SP-partitioning only) |

#### SP-Score Thresholds

The SP (Sum-of-Pairs) score threshold controls when tronko-build splits a tree into subtrees. Lower = fewer, larger trees. Higher = more, smaller trees.

| Threshold | Effect |
|---|---|
| `sp_0.05` | Fewer, larger trees (less aggressive partitioning) |
| `sp_0.10` | Default balance |
| `sp_0.20` | More, smaller trees (more aggressive partitioning) |

#### Sweep Results Summary

| Config | SP | Trees | Leaves/Tree (median) | DB Size | Build Time | Genus Purity | Family Purity |
|---|---|---|---|---|---|---|---|
| ac_default | 0.05 | 4,378 | 22 | 1,551 MB | 43 min | 0.60 | 0.77 |
| ac_default | 0.10 | 7,996 | 12 | 1,473 MB | 45 min | 0.65 | 0.80 |
| ac_default | 0.20 | 10,602 | 8 | 1,426 MB | 46 min | 0.66 | 0.81 |
| ac_more_bins | 0.05 | 4,347 | 22 | 1,555 MB | 30 min | 0.57 | 0.74 |
| ac_more_bins | 0.10 | 8,013 | 12 | 1,474 MB | 32 min | 0.62 | 0.78 |
| ac_more_bins | 0.20 | 10,575 | 8 | 1,429 MB | 33 min | 0.64 | 0.79 |
| ac_fewer_bins | 0.05 | 1,574 | — | 4,231 MB | 149 min | 0.41 | 0.52 |
| ac_fewer_bins | 0.10 | 2,902 | — | 4,200 MB | 150 min | 0.43 | 0.53 |
| ac_fewer_bins | 0.20 | 3,842 | — | 4,181 MB | 158 min | 0.43 | 0.53 |
| ac_single | 0.05 | 2,909 | — | 2,000 MB | — | — | — |
| ac_single | 0.10 | 7,005 | — | 1,600 MB | ~3.5 hr | — | — |
| ac_single | 0.20 | 10,605 | — | 1,400 MB | — | — | — |

**Key findings**:
- `ac_fewer_bins` produces 3x larger databases with worse purity — fewer initial clusters lead to larger, more taxonomically heterogeneous trees
- `ac_default` and `ac_more_bins` produce very similar results (AncestralClust is robust to bin count within this range)
- Higher SP thresholds improve genus/family purity at the cost of more trees (more partition overhead at assignment time)
- `ac_single` (no clustering, just SP-partitioning of one big tree) produces competitive results, suggesting AncestralClust may not be strictly necessary

**Build command** (all sweep DBs):
```bash
./build-tronko-db.sh \
    -f 12SV5_species.fasta \
    -t 12SV5_species_taxonomy.txt \
    -o dbs/<cluster_config>/<sp_threshold>/ \
    -T 12 -s <sp_threshold> -a
```

Sweep metadata (build times, tree statistics, purity metrics) in `sweep_metadata/`.

---

### `sweep_metadata/`

Logs and metrics from the parameter sweep:
- `sweep_summary.json` — Per-config statistics (trees, leaves, purity, build time)
- `sweep_metrics.json` — Extended metrics
- `sweep_report.html` — Visual report
- `sweep.log` — Build log
- `extended_analysis.json` — Additional analysis

---

## Database File Contents

Each complete database directory contains:

| File | Description |
|---|---|
| `reference_tree.txt` | Tronko reference database (text format) |
| `reference_tree.txt.gz` | Gzipped text format |
| `reference_tree.trkb` | Binary format with zstd compression (smallest) |
| `marker.fasta` | Concatenated leaf sequences (for BWA indexing) |
| `marker.fasta.{amb,ann,bwt,pac,sa}` | BWA index files |
| `final_partitions.txt` | Partition manifest (which leaves in which tree) |

Not all directories have all formats. The minimum needed for `tronko-assign` is one reference tree file + `marker.fasta` + BWA index.

## Using a Database with tronko-assign

```bash
tronko-assign -r \
    -f databases/<db>/reference_tree.txt.gz \
    -a databases/<db>/marker.fasta \
    -s -g query_reads.fasta \
    -o results.txt -w
```

## DVC Tracking

Large database files are tracked with [DVC](https://dvc.org/) and not stored in git.

```bash
# Pull databases (after cloning)
dvc pull

# Check status
dvc status
```
