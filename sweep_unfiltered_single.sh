#!/bin/bash
# Unfiltered 176K sweep: ac_default (bin=20K, descendants=75)
# Inlines Steps 0-3 from build-tronko-db.sh, then runs tronko-build per threshold.
#
# Usage: nohup bash sweep_unfiltered_default.sh > databases/sweep_unfiltered_default.log 2>&1 &
set -euo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
export PATH="$REPO/bin:$REPO/tronko-build:$PATH"

# === Config ===
CONFIG="ac_single"
AC_BIN_SIZE=20000
AC_DESCENDANTS=75
AC_CUTOFF=999999
THREADS=15
PARALLEL_JOBS=1
THRESHOLDS=(2.00 2.50 2.75)

INPUT_FASTA="$REPO/databases/12sv5_unfiltered/input.fasta"
INPUT_TAXONOMY="$REPO/databases/12sv5_unfiltered/input_taxonomy.txt"
CACHE_DIR="$REPO/databases/sweep_clusters_unfiltered/$CONFIG"
AC_DIR="$CACHE_DIR/ancestralclust"
NEWICK_DIR="$CACHE_DIR/newick"
MERGED_DIR="$CACHE_DIR/merged"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }

log "=== Unfiltered sweep: $CONFIG (bin=$AC_BIN_SIZE, desc=$AC_DESCENDANTS) ==="

mkdir -p "$AC_DIR" "$NEWICK_DIR" "$MERGED_DIR"
NUM_SEQS=$(grep -c '^>' "$INPUT_FASTA")
USE_FASTTREE=1

# Helper
unwrap_fasta() {
    awk '/^>/ { if (seq) print seq; print; seq="" ; next } { seq = seq $0 } END { if (seq) print seq }' "$1" > "$2"
}

# =========================================================================
# STEP 0: Sanitize
# =========================================================================
CLEAN_FASTA="$CACHE_DIR/input_clean.fasta"
if [[ -f "$CLEAN_FASTA" ]]; then
    log "Sanitized FASTA cached"
else
    log "Sanitizing FASTA..."
    sed '/^>/!s/[RYWSMKHBVDrywsmkhbvd]/N/g' "$INPUT_FASTA" > "$CLEAN_FASTA"
fi

# =========================================================================
# STEP 1: AncestralClust
# =========================================================================
log "=== Step 1: AncestralClust ==="
if [[ -f "$AC_DIR/.step1_done" ]]; then
    NUM_CLUSTERS=$(ls "$AC_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
    log "CACHED: $NUM_CLUSTERS clusters"
    # Generate taxonomy for any clusters missing it
    for cluster_fasta in "$AC_DIR"/*.fasta; do
        cluster_id=$(basename "$cluster_fasta" .fasta)
        if [[ ! -f "$AC_DIR/${cluster_id}_taxonomy.txt" ]]; then
            log "  Generating missing taxonomy for cluster $cluster_id..."
            grep '^>' "$cluster_fasta" | sed 's/^>//' | cut -d' ' -f1 | sort > "$CACHE_DIR/accs_${cluster_id}.txt"
            sort "$INPUT_TAXONOMY" | join -t$'\t' "$CACHE_DIR/accs_${cluster_id}.txt" - > "$AC_DIR/${cluster_id}_taxonomy.txt"
            rm "$CACHE_DIR/accs_${cluster_id}.txt"
        fi
    done
elif [[ "$NUM_SEQS" -le "$AC_CUTOFF" ]]; then
    log "Below cutoff — single cluster"
    cp "$CLEAN_FASTA" "$AC_DIR/0.fasta"
    cp "$INPUT_TAXONOMY" "$AC_DIR/0_taxonomy.txt"
    NUM_CLUSTERS=1
    touch "$AC_DIR/.step1_done"
else
    NUM_BINS=$(( (NUM_SEQS + AC_BIN_SIZE - 1) / AC_BIN_SIZE ))
    NUM_SEEDS=$(( NUM_BINS * 30 ))
    [[ "$NUM_SEEDS" -gt 4000 ]] && NUM_SEEDS=4000
    NUM_LINES=$(wc -l < "$CLEAN_FASTA")
    AC_THREADS=4

    log "AncestralClust: $NUM_SEQS seqs -> $NUM_BINS bins, $NUM_SEEDS seeds"
    ancestralclust -f \
        -i "$CLEAN_FASTA" \
        -b "$NUM_BINS" \
        -r "$NUM_SEEDS" \
        -p "$AC_DESCENDANTS" \
        -l "$NUM_LINES" \
        -c "$AC_THREADS" \
        -d "$AC_DIR"

    log "Splitting taxonomy by cluster..."
    for cluster_fasta in "$AC_DIR"/*.fasta; do
        cluster_id=$(basename "$cluster_fasta" .fasta)
        grep '^>' "$cluster_fasta" | sed 's/^>//' | cut -d' ' -f1 | sort > "$CACHE_DIR/accs_${cluster_id}.txt"
        sort "$INPUT_TAXONOMY" | join -t$'\t' "$CACHE_DIR/accs_${cluster_id}.txt" - > "$AC_DIR/${cluster_id}_taxonomy.txt"
        rm "$CACHE_DIR/accs_${cluster_id}.txt"
        log "  Cluster $cluster_id: $(grep -c '^>' "$cluster_fasta") seqs"
    done
    NUM_CLUSTERS=$(ls "$AC_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
    touch "$AC_DIR/.step1_done"
fi
log "Step 1 complete: $NUM_CLUSTERS clusters"

# =========================================================================
# STEP 2: Per-cluster FAMSA + FastTree
# =========================================================================
log "=== Step 2: FAMSA + FastTree ==="

build_cluster_tree() {
    local cluster_fasta="$1"
    local cluster_id=$(basename "$cluster_fasta" .fasta)
    local outdir="$2"
    local cluster_threads="${3:-$THREADS}"
    local n_seqs=$(grep -c '^>' "$cluster_fasta")

    # Dedup
    local dedup_fasta="$outdir/${cluster_id}_dedup.fasta"
    awk '/^>/{gsub(/:/, "_"); split($0,a," "); h=a[1]; if(seen[h]++){p=0; next} p=1; print; next} p{print}' "$cluster_fasta" > "$dedup_fasta"
    local n_dedup=$(grep -c '^>' "$dedup_fasta")
    log "  Cluster $cluster_id ($n_seqs seqs, $n_dedup unique): aligning..."

    if [[ "$n_dedup" -lt 3 ]]; then
        log "  Cluster $cluster_id: only $n_dedup unique, skipping"
        rm -f "$dedup_fasta"
        return 1
    fi

    famsa -t "$cluster_threads" "$dedup_fasta" "$outdir/${cluster_id}_MSA.fasta" 2>/dev/null
    rm -f "$dedup_fasta"
    unwrap_fasta "$outdir/${cluster_id}_MSA.fasta" "$outdir/${cluster_id}_MSA_unwrapped.fasta"
    mv "$outdir/${cluster_id}_MSA_unwrapped.fasta" "$outdir/${cluster_id}_MSA.fasta"

    log "  Cluster $cluster_id: building tree..."
    set +e
    OMP_NUM_THREADS="$cluster_threads" FastTree -nt -gtr -nosupport "$outdir/${cluster_id}_MSA.fasta" \
        > "$outdir/RAxML_bestTree.${cluster_id}.unrooted" 2>"$outdir/${cluster_id}_tree.log"
    local tree_rc=$?
    set -e
    if [[ "$tree_rc" -ne 0 ]] || [[ ! -s "$outdir/RAxML_bestTree.${cluster_id}.unrooted" ]]; then
        log "  WARNING: FastTree failed for cluster $cluster_id"
        return 1
    fi

    # Strip support values
    python3 -c "
import re, sys
with open(sys.argv[1]) as f: data = f.read()
with open(sys.argv[1], 'w') as f: f.write(re.sub(r'\)([0-9][0-9.eE+-]*):', '):', data))
" "$outdir/RAxML_bestTree.${cluster_id}.unrooted"

    set +e
    nw_reroot "$outdir/RAxML_bestTree.${cluster_id}.unrooted" \
        > "$outdir/RAxML_bestTree.${cluster_id}.reroot" 2>/dev/null
    local reroot_rc=$?
    set -e
    if [[ "$reroot_rc" -ne 0 ]] || [[ ! -s "$outdir/RAxML_bestTree.${cluster_id}.reroot" ]]; then
        cp "$outdir/RAxML_bestTree.${cluster_id}.unrooted" "$outdir/RAxML_bestTree.${cluster_id}.reroot"
    fi
    log "  Cluster $cluster_id: done"
}

CLUSTER_THREADS=$(( THREADS / PARALLEL_JOBS ))
[[ "$CLUSTER_THREADS" -lt 1 ]] && CLUSTER_THREADS=1

for cluster_fasta in "$AC_DIR"/*.fasta; do
    cluster_id=$(basename "$cluster_fasta" .fasta)
    n_seqs=$(grep -c '^>' "$cluster_fasta")
    [[ "$n_seqs" -lt 3 ]] && continue

    if [[ -f "$NEWICK_DIR/RAxML_bestTree.${cluster_id}.reroot" ]] && [[ -s "$NEWICK_DIR/RAxML_bestTree.${cluster_id}.reroot" ]]; then
        log "  Cluster $cluster_id ($n_seqs seqs): CACHED"
        continue
    fi

    build_cluster_tree "$cluster_fasta" "$NEWICK_DIR" "$CLUSTER_THREADS"

    if [[ -f "$AC_DIR/${cluster_id}_taxonomy.txt" ]]; then
        awk -F'\t' '{gsub(/:/, "_", $1)} !seen[$1]++' "$AC_DIR/${cluster_id}_taxonomy.txt" > "$NEWICK_DIR/${cluster_id}_taxonomy.txt"
    fi
done
log "Step 2 complete"

# =========================================================================
# STEP 3: Merge/renumber
# =========================================================================
log "=== Step 3: Merge ==="
counter=0
for tree_file in "$NEWICK_DIR"/RAxML_bestTree.*.reroot; do
    orig_id=$(basename "$tree_file" | sed 's/RAxML_bestTree\.\(.*\)\.reroot/\1/')
    [[ ! -f "$NEWICK_DIR/${orig_id}_MSA.fasta" ]] && continue
    [[ ! -f "$NEWICK_DIR/${orig_id}_taxonomy.txt" ]] && continue

    cp "$tree_file" "$MERGED_DIR/RAxML_bestTree.${counter}.reroot"
    cp "$NEWICK_DIR/${orig_id}_MSA.fasta" "$MERGED_DIR/${counter}_MSA.fasta"
    cp "$NEWICK_DIR/${orig_id}_taxonomy.txt" "$MERGED_DIR/${counter}_taxonomy.txt"
    [[ -f "$AC_DIR/${orig_id}.fasta" ]] && cp "$AC_DIR/${orig_id}.fasta" "$MERGED_DIR/${counter}.fasta"
    counter=$((counter + 1))
done
N_CLUSTERS=$counter
log "Step 3 complete: $N_CLUSTERS clusters"

# =========================================================================
# STEP 4: tronko-build per threshold
# =========================================================================
for sp in "${THRESHOLDS[@]}"; do
    sp_label=$(printf "%.2f" "$sp")
    db_name="unfiltered_${CONFIG}_sp_${sp_label}"
    db_dir="$REPO/databases/$db_name"

    if [[ -f "$db_dir/reference_tree.txt" ]]; then
        log "CACHED: $db_name"
        continue
    fi

    log "--- $db_name (SP=$sp, $N_CLUSTERS clusters) ---"
    mkdir -p "$db_dir"

    t0=$(date +%s)
    "$REPO/tronko-build/tronko-build" -y \
        -e "$MERGED_DIR" \
        -n "$N_CLUSTERS" \
        -d "$db_dir" \
        -s -u "$sp" \
        -a -E \
        -c "$THREADS" \
        2>&1 | tee "$db_dir/build.log"
    t1=$(date +%s)

    if [[ ! -f "$db_dir/reference_tree.txt" ]]; then
        log "ERROR: tronko-build failed for $db_name"
        continue
    fi

    N_TREES=$(head -1 "$db_dir/reference_tree.txt")
    SIZE_MB=$(du -m "$db_dir/reference_tree.txt" | cut -f1)
    log "  $db_name: $N_TREES trees, ${SIZE_MB}MB, $((t1-t0))s"

    # marker.fasta from exported subtrees or merged clusters
    > "$db_dir/marker.fasta"
    if [[ -d "$db_dir/exported_subtrees" ]]; then
        for msa in "$db_dir/exported_subtrees"/*_MSA.fasta; do
            awk '/^>/{print; next} {gsub(/-/,""); print}' "$msa" >> "$db_dir/marker.fasta"
        done
    else
        for i in $(seq 0 $((N_CLUSTERS - 1))); do
            [[ -f "$MERGED_DIR/${i}_MSA.fasta" ]] && awk '/^>/{print; next} {gsub(/-/,""); print}' "$MERGED_DIR/${i}_MSA.fasta" >> "$db_dir/marker.fasta"
        done
    fi
    bwa index "$db_dir/marker.fasta" 2>/dev/null

    cat > "$db_dir/build_info.json" <<EOF
{
  "name": "$db_name",
  "pipeline": "AncestralClust + FastTree + corrected SP partitioning",
  "marker": "12SV5",
  "input_sequences": $NUM_SEQS,
  "clustering": "$CONFIG",
  "clustering_params": {"bin_size": $AC_BIN_SIZE, "descendants": $AC_DESCENDANTS},
  "sp_threshold": $sp,
  "sp_normalization": "corrected",
  "n_trees": $N_TREES,
  "n_input_clusters": $N_CLUSTERS,
  "db_size_mb": $SIZE_MB,
  "build_seconds": $((t1-t0)),
  "date": "$(date +%Y-%m-%d)"
}
EOF
done

log "=== $CONFIG sweep complete ==="
