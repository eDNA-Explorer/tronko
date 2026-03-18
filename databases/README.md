# Tronko Reference Databases

All tronko reference databases built during this project. Each database directory contains a `build_info.json` with exact build parameters and a `pasta_stats.json` (if PASTA-based) with alignment scores.

**Master index**: `metadata.json` — machine-readable inventory of all databases with completeness and ablation-readiness status.

**Marker gene**: 12SV5 (vertebrate 12S rRNA)
**Input sequences**: 101,088 (filtered) or 176,947 (unfiltered) vertebrate amplicon references from CRUXv2

---

## Database Inventory

### 1. `charadriiformes_single_tree/`

Small COI test dataset (1,413 seqs, 1 tree). Used for development and correctness verification. Not ablation-ready.

---

### 2-3. `vert12S_fasttree/`, `vert12S_veryfasttree/`

Early builds using AncestralClust + FastTree/VeryFastTree + SP partitioning. Superseded by PASTA-based approach.

---

### 4. `vert12S_pasta/`

**The original winning database (F1=0.696).** PASTA tree decomposed into 141 partitions (max_size=1000), FAMSA+FastTree per partition, tronko-build with no SP partitioning and no internal repartitioning. Built with an older tronko-build binary that did not trigger `createNewRoots`.

| Property | Value |
|---|---|
| Trees | 141 |
| DB size | 2.4 GB |
| PASTA | iter=3, sub=200 |
| Decomposition | centroid, max_size=1000 |
| Internal repartitioning | none |
| SP partitioning | none |

---

### 4b. `vert12S_pasta_sp_0.10/`

PASTA tree + diameter decomposition (max_diam=25) + SP 0.10 recursive partitioning. 7,804 trees. Performed worse than the no-SP version (4b). Ablation-ready with exported subtrees.

---

### 5. `12sv5_unfiltered/`

Unfiltered 12SV5 dataset (176,947 sequences). AncestralClust (fewer bins, bin=50K) + FastTree + SP 0.10 partitioning. 14,558 trees. Ablation-ready with exported subtrees.

---

### 6-17. `sweep_ac_*` (12 databases)

AncestralClust clustering sweep: 4 configs x 3 SP thresholds. All use SP-score partitioning. Built with `build-tronko-db.sh -F -E`. Ablation-ready via `exported_subtrees/`.

| Config | SP 0.05 | SP 0.10 | SP 0.20 |
|---|---|---|---|
| ac_default (bin=20K) | 4,376 trees | 7,992 trees | 10,538 trees |
| ac_more_bins (bin=10K) | 4,335 trees | 7,995 trees | 10,623 trees |
| ac_fewer_bins (bin=50K) | 1,574 trees | 2,908 trees | 3,848 trees |
| ac_single (no clustering) | 2,909 trees | 7,005 trees | 10,605 trees |

---

### 18+. `pasta_sweep/` (PASTA decomposition sweep)

Systematic sweep of PASTA tree decomposition parameters. All databases use **no SP partitioning** and **no internal repartitioning** (`-f 999999`).

#### Two PASTA trees tested

| PASTA Config | Subproblem Size | Final Score | Input |
|---|---|---|---|
| sub=200 (original) | 200 | -3,309,385 | 101K filtered |
| sub=500 | 500 | -3,168,989 (3.5% better) | 101K filtered |
| sub=500 unfiltered | 500 | (in progress) | 176K unfiltered |

#### Decomposition parameters swept

| Parameter | Values Tested |
|---|---|
| `--max-size` | 500, 1000, 2000, 5000 |
| `--max-diam` | 10, 25 |
| FastTree `-gamma` | false, true |
| `norepartition` (`-f 999999`) | all configs use this |

#### Pipeline

```
PASTA tree (101K or 176K seqs, 3 iterations)
  → centroid decomposition (max_size or max_diam)
    → per-partition FAMSA + FastTree
      → tronko-build -y -E -a -f 999999 (no repartitioning)
        → marker.fasta + BWA index
```

#### Current databases in `pasta_sweep/`

**sub=200 PASTA tree:**
- `maxsize{500,1000,2000,5000}_{nogamma,gamma}` — 8 databases
- `maxdiam{10,25}_{nogamma,gamma}` — 4 databases
- `maxsize1000_norepartition_{nogamma,gamma}` — 2 databases (identical pipeline, named explicitly)

**sub=500 PASTA tree:**
- `sub500_maxsize{500,1000,2000,5000}_{nogamma,gamma}` — 8 databases
- `sub500_maxsize1000_norepartition_{nogamma,gamma}` — 2 databases

**sub=500 unfiltered (176K, in progress):**
- `unfiltered_maxsize{500,1000,2000,5000}_{nogamma,gamma}` — pending
- `unfiltered_maxdiam{10,25}_{nogamma,gamma}` — pending

---

### `sweep_metadata/`

Logs and metrics from the AncestralClust parameter sweep (sweep_summary.json, sweep_metrics.json, etc).

---

## Database File Contents

Each complete database directory contains:

| File | Description |
|---|---|
| `reference_tree.txt` | Tronko reference database (text) |
| `marker.fasta` | Concatenated gap-free leaf sequences (for BWA) |
| `marker.fasta.{amb,ann,bwt,pac,sa}` | BWA index files |
| `marker_taxonomy.txt` | Taxonomy for all sequences in marker.fasta |
| `build_info.json` | Build parameters and provenance |
| `pasta_stats.json` | PASTA scores, config, runtime (if PASTA-based) |
| `exported_subtrees/` | Per-partition MSA + taxonomy + tree (for ablation) |
| `tree_list.txt` | Partition manifest |
| `input.fasta`, `input_taxonomy.txt` | Input files used for the build |

## Using a Database

```bash
tronko-assign -r \
    -f databases/<db>/reference_tree.txt \
    -a databases/<db>/marker.fasta \
    -s -g query_reads.fasta \
    -o results.txt -w
```

## Benchmarking

```bash
# Discover all databases and run ablation benchmarks
cd /path/to/assignment-tool-benchmarking
poetry run python notebooks/tronko_db_sweep.py --db-only \
    --db-sweep-root /path/to/tronko-fork/databases
```

## DVC Tracking

Large files tracked with DVC against `gs://edna-dvc-data/tronko-databases`.

```bash
poetry install        # sets up DVC
dvc pull              # all databases
dvc pull databases/pasta_sweep/maxsize1000_nogamma.dvc  # one database
```
