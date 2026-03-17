#!/bin/bash
# Ablate sequences from a tronko reference database without full rebuild.
#
# Requires exported subtrees from a master build (tronko-build -E).
# Only recomputes posteriors — skips FAMSA, tree building, and SP-score partitioning.
#
# Usage:
#   ablate-tronko-db.sh -m <exported_subtrees_dir> -r <remove_list.txt> -o <output_dir> \
#                        [-T threads] [-p primer_name]
#
# Options:
#   -m  Directory of exported subtrees (from tronko-build -E)
#   -r  File listing sequence accessions to remove (one per line)
#   -o  Output directory for ablated database
#   -T  Threads for tronko-build (default: 8)
#   -p  Primer/marker name for output FASTA (default: marker)

set -euo pipefail

THREADS=8
PRIMER="marker"

while getopts "m:r:o:T:p:" opt; do
    case $opt in
        m) SUBTREE_DIR="$OPTARG" ;;
        r) REMOVE_LIST="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        p) PRIMER="$OPTARG" ;;
        *) echo "Usage: $0 -m <subtrees_dir> -r <remove_list> -o <output_dir> [-T threads] [-p primer]" >&2; exit 1 ;;
    esac
done

# Validate required args
if [[ -z "${SUBTREE_DIR:-}" || -z "${REMOVE_LIST:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    echo "ERROR: -m, -r, and -o are required" >&2
    exit 1
fi

if [[ ! -d "$SUBTREE_DIR" ]]; then
    echo "ERROR: Subtree directory not found: $SUBTREE_DIR" >&2
    exit 1
fi

if [[ ! -f "$REMOVE_LIST" ]]; then
    echo "ERROR: Remove list not found: $REMOVE_LIST" >&2
    exit 1
fi

for tool in tronko-build nw_prune bwa; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found on PATH" >&2
        exit 1
    fi
done

TOTAL_REMOVE=$(wc -l < "$REMOVE_LIST" | tr -d ' ')
echo "================================================================"
echo "  tronko-build ablation"
echo "================================================================"
echo "  Subtrees:    $SUBTREE_DIR"
echo "  Remove list: $REMOVE_LIST ($TOTAL_REMOVE sequences)"
echo "  Output:      $OUTPUT_DIR"
echo "  Threads:     $THREADS"
echo "================================================================"

# Setup working directory
CACHE_DIR="$OUTPUT_DIR/.cache/ablated"
mkdir -p "$CACHE_DIR" "$OUTPUT_DIR"

# Sort remove list for efficient matching
SORTED_REMOVES="$CACHE_DIR/_removes_sorted.txt"
sort "$REMOVE_LIST" > "$SORTED_REMOVES"

# Step 1: Copy and prune each cluster
echo ""
echo "=== Step 1: Pruning subtrees ==="
removed_total=0
skipped_clusters=0
valid_cluster_count=0

for msa_file in "$SUBTREE_DIR"/*_MSA.fasta; do
    cluster_id=$(basename "$msa_file" _MSA.fasta)
    tree_file="$SUBTREE_DIR/RAxML_bestTree.${cluster_id}.reroot"
    tax_file="$SUBTREE_DIR/${cluster_id}_taxonomy.txt"

    if [[ ! -f "$tree_file" ]] || [[ ! -f "$tax_file" ]]; then
        echo "  Cluster $cluster_id: missing tree or taxonomy, skipping"
        skipped_clusters=$((skipped_clusters + 1))
        continue
    fi

    # Find which sequences from remove list exist in this cluster
    grep '^>' "$msa_file" | sed 's/^>//' | sort > "$CACHE_DIR/_cluster_accs.txt"
    cluster_removes="$CACHE_DIR/_cluster_${cluster_id}_removes.txt"
    comm -12 "$SORTED_REMOVES" "$CACHE_DIR/_cluster_accs.txt" > "$cluster_removes"

    n_remove=$(wc -l < "$cluster_removes" | tr -d ' ')
    n_total=$(wc -l < "$CACHE_DIR/_cluster_accs.txt" | tr -d ' ')
    n_remaining=$((n_total - n_remove))

    if [[ "$n_remove" -eq 0 ]]; then
        # No removals — copy as-is
        cp "$msa_file" "$CACHE_DIR/${valid_cluster_count}_MSA.fasta"
        cp "$tree_file" "$CACHE_DIR/RAxML_bestTree.${valid_cluster_count}.reroot"
        cp "$tax_file" "$CACHE_DIR/${valid_cluster_count}_taxonomy.txt"
        echo "  Cluster $cluster_id -> $valid_cluster_count: no removals ($n_total seqs)"
        valid_cluster_count=$((valid_cluster_count + 1))
        continue
    fi

    if [[ "$n_remaining" -lt 3 ]]; then
        echo "  Cluster $cluster_id: only $n_remaining seqs after pruning, dropping cluster"
        skipped_clusters=$((skipped_clusters + 1))
        removed_total=$((removed_total + n_remove))
        continue
    fi

    echo "  Cluster $cluster_id -> $valid_cluster_count: removing $n_remove/$n_total seqs ($n_remaining remaining)"

    # Prune tree
    nw_prune -f "$tree_file" "$cluster_removes" \
        > "$CACHE_DIR/RAxML_bestTree.${valid_cluster_count}.reroot" 2>/dev/null

    # Filter FASTA: remove matching sequences
    awk -v rmfile="$cluster_removes" '
        BEGIN { while ((getline line < rmfile) > 0) rm[line]=1 }
        /^>/ { name=substr($0,2); split(name,a," "); skip=rm[a[1]] ? 1 : 0 }
        !skip { print }
    ' "$msa_file" > "$CACHE_DIR/${valid_cluster_count}_MSA.fasta"

    # Filter taxonomy
    awk -F'\t' -v rmfile="$cluster_removes" '
        BEGIN { while ((getline line < rmfile) > 0) rm[line]=1 }
        !rm[$1]
    ' "$tax_file" > "$CACHE_DIR/${valid_cluster_count}_taxonomy.txt"

    removed_total=$((removed_total + n_remove))
    valid_cluster_count=$((valid_cluster_count + 1))
done

# Check for sequences not found in any cluster
found_removes=$(comm -12 "$SORTED_REMOVES" <(cat "$SUBTREE_DIR"/*_MSA.fasta | grep '^>' | sed 's/^>//' | sort -u) 2>/dev/null | wc -l | tr -d ' ')
not_found=$((TOTAL_REMOVE - found_removes))
if [[ "$not_found" -gt 0 ]]; then
    echo "  WARNING: $not_found sequence(s) from remove list not found in any cluster"
fi

# Clean up temp files
rm -f "$CACHE_DIR"/_cluster_*.txt "$CACHE_DIR"/_removes_sorted.txt "$CACHE_DIR"/_cluster_accs.txt

echo ""
echo "  Removed: $removed_total sequences"
echo "  Dropped clusters: $skipped_clusters"
echo "  Remaining clusters: $valid_cluster_count"

if [[ "$valid_cluster_count" -eq 0 ]]; then
    echo "ERROR: No valid clusters remaining after ablation" >&2
    exit 1
fi

# Step 2: Run tronko-build (posteriors only — no repartitioning)
echo ""
echo "=== Step 2: tronko-build (posteriors only) ==="
echo "Running: tronko-build -y -e $CACHE_DIR -n $valid_cluster_count -d $OUTPUT_DIR -v -f 999999 -a -c $THREADS"

time tronko-build -y \
    -e "$CACHE_DIR" \
    -n "$valid_cluster_count" \
    -d "$OUTPUT_DIR" \
    -v -f 999999 \
    -a \
    -c "$THREADS"

if [[ ! -f "$OUTPUT_DIR/reference_tree.txt" ]]; then
    echo "ERROR: tronko-build did not produce reference_tree.txt" >&2
    exit 1
fi

# Step 3: Finalize (concat FASTA + BWA index + gzip)
echo ""
echo "=== Step 3: Finalize ==="

# Build marker.fasta by stripping gaps from aligned MSA (BWA needs unaligned sequences)
for msa in "$CACHE_DIR"/*_MSA.fasta; do
    awk '/^>/ { print; next } { gsub(/-/, ""); print }' "$msa"
done > "$OUTPUT_DIR/${PRIMER}.fasta"

echo "Building BWA index..."
bwa index "$OUTPUT_DIR/${PRIMER}.fasta"

gzip "$OUTPUT_DIR/reference_tree.txt"

echo ""
echo "================================================================"
echo "  Ablation complete!"
echo "================================================================"
echo "  Sequences removed: $removed_total"
echo "  Clusters:          $valid_cluster_count"
echo "  Output files:"
ls -lh "$OUTPUT_DIR"/ | grep -v '^total' | sed 's/^/    /'
echo "================================================================"
