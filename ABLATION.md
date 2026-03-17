# Tronko Ablation Studies

Remove sequences from a tronko reference database without a full rebuild. Only recomputes posteriors — skips FAMSA alignment, tree building, and SP-score partitioning.

## Prerequisites

1. A completed master build with exported subtrees
2. Tools on PATH: `tronko-build`, `nw_prune`, `bwa`

## Step 1: Master Build with Subtree Export

Run the normal build pipeline with the `-E` flag added to the tronko-build invocation. This exports the final post-partition sub-trees to `exported_subtrees/`.

```bash
# Add -E to the tronko-build call inside build-tronko-db.sh,
# or run tronko-build manually after partitioning:
build-tronko-db.sh -f input.fasta -t taxonomy.txt -o master/ -T 8 -s 0.1 -F

# Then export subtrees (if -E wasn't in the pipeline):
tronko-build -y -e master/.cache/merged -n <N_CLUSTERS> -d master/ -s -u 0.1 -a -c 8 -E
```

This creates `master/exported_subtrees/` containing, for each cluster:
- `{N}_MSA.fasta` — unaligned FASTA
- `RAxML_bestTree.{N}.reroot` — Newick tree
- `{N}_taxonomy.txt` — taxonomy file

This only needs to be done once per marker/database.

## Step 2: Run Ablations

### From the command line

```bash
# Create a remove list (one accession per line)
cat > remove_list.txt << 'EOF'
gi|1395189203|gb|AC277927.1|_194159-194314
gi|1234567890|gb|XY123456.1|_1-500
EOF

# Run ablation
./ablate-tronko-db.sh \
    -m master/exported_subtrees/ \
    -r remove_list.txt \
    -o ablation_01/ \
    -T 8 \
    -p 12SV5
```

### From Python

```python
import subprocess
import tempfile
from pathlib import Path


def ablate_tronko_db(
    subtrees_dir: str,
    remove_accessions: list[str],
    output_dir: str,
    threads: int = 8,
    primer: str = "marker",
    ablate_script: str = "./ablate-tronko-db.sh",
) -> subprocess.CompletedProcess:
    """Remove sequences from a tronko reference database.

    Args:
        subtrees_dir: Path to exported_subtrees/ from master build.
        remove_accessions: List of accession strings to remove.
        output_dir: Where to write the ablated database.
        threads: Threads for tronko-build posterior computation.
        primer: Marker name for output FASTA (e.g. "12SV5").
        ablate_script: Path to ablate-tronko-db.sh.

    Returns:
        CompletedProcess with returncode, stdout, stderr.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, prefix="tronko_remove_"
    ) as f:
        f.write("\n".join(remove_accessions) + "\n")
        remove_list_path = f.name

    try:
        result = subprocess.run(
            [
                ablate_script,
                "-m", subtrees_dir,
                "-r", remove_list_path,
                "-o", output_dir,
                "-T", str(threads),
                "-p", primer,
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Ablation failed (exit {result.returncode}):\n{result.stderr}"
            )
        return result
    finally:
        Path(remove_list_path).unlink(missing_ok=True)
```

#### Batch ablation example

```python
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

MASTER_SUBTREES = "master/exported_subtrees"
ABLATIONS_DIR = "ablations"

# Define ablation experiments
experiments = {
    "no_mammals": ["acc1", "acc2", "acc3"],
    "no_birds": ["acc4", "acc5"],
    "leave_one_out_42": ["acc42"],
}

# Run sequentially (each ablation uses multiple threads internally)
for name, remove_list in experiments.items():
    output = f"{ABLATIONS_DIR}/{name}"
    print(f"Running ablation: {name} ({len(remove_list)} sequences)")
    ablate_tronko_db(
        subtrees_dir=MASTER_SUBTREES,
        remove_accessions=remove_list,
        output_dir=output,
        threads=8,
        primer="12SV5",
    )
    print(f"  Done: {output}/reference_tree.txt.gz")
```

#### Getting accessions from the master FASTA

```python
def get_accessions(fasta_path: str) -> list[str]:
    """Extract all accession IDs from a FASTA file."""
    accessions = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                accessions.append(line[1:].strip().split()[0])
    return accessions

# Example: leave-one-out for a specific species
all_accs = get_accessions("master/marker.fasta")
# Filter by taxonomy or other criteria to build remove lists
```

## Accession Format

The build pipeline replaces `:` with `_` in FASTA headers (to prevent VeryFastTree name truncation). Your remove list must use the post-rename format:

| Original FASTA header | Remove list entry |
|---|---|
| `>gi\|123\|gb\|AC277927.1\|:194-314` | `gi\|123\|gb\|AC277927.1\|_194-314` |
| `>Homo_sapiens_ABC123` | `Homo_sapiens_ABC123` |

To see the exact accessions in the exported subtrees:
```bash
grep '^>' master/exported_subtrees/*_MSA.fasta | sed 's/.*>//' | sort -u > all_accessions.txt
```

## Output

The ablated database has the same structure as a normal tronko-build output:

```
ablation_01/
├── reference_tree.txt.gz    # Ablated reference database
├── 12SV5.fasta              # Concatenated sequences (for BWA)
├── 12SV5.fasta.amb          # BWA index files
├── 12SV5.fasta.ann
├── 12SV5.fasta.bwt
├── 12SV5.fasta.pac
├── 12SV5.fasta.sa
└── .cache/ablated/          # Intermediate files (can be deleted)
```

Use with tronko-assign as normal:
```bash
tronko-assign -r \
    -f ablation_01/reference_tree.txt.gz \
    -a ablation_01/12SV5.fasta \
    -s -g query_reads.fasta \
    -o results.txt -w
```

## Performance

| Step | Full rebuild (101K seqs) | Ablation |
|---|---|---|
| FAMSA + VeryFastTree | ~30 min | Skipped |
| SP-score partitioning | ~15 min | Skipped |
| Posterior computation | ~5 min | ~5 min |
| BWA indexing | ~2 min | ~2 min |
| **Total** | **~50 min** | **~7 min** |

Actual times depend on dataset size and hardware. The key savings come from skipping alignment and tree building entirely.
