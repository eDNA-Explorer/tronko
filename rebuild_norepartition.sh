#!/bin/bash
# Rebuild all PASTA sweep databases with norepartition (-f 999999)
# Runs sequentially with low thread count alongside the overnight PASTA run
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

TRONKO_BUILD="tronko-build/tronko-build"
SWEEP_DIR="databases/pasta_sweep"
THREADS=4

FILTERED_FASTA="databases/vert12S_pasta/input.fasta"
FILTERED_TAX="databases/vert12S_pasta/input_taxonomy.txt"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }

build_one() {
    local part_dir="$1"
    local db_dir="$2"
    local label="$3"
    local use_gamma="$4"

    local ref_tree="$db_dir/reference_tree.txt"
    if [[ -f "$ref_tree" ]]; then
        log "  CACHED: $label"
        return
    fi

    mkdir -p "$db_dir"
    local n_partitions=$(ls "$part_dir"/RAxML_bestTree.*.reroot 2>/dev/null | wc -l | tr -d ' ')

    # Gamma reoptimization
    local tronko_input="$part_dir"
    if [[ "$use_gamma" == "true" ]]; then
        local gamma_dir="$SWEEP_DIR/_cache/gamma_norepartition_$(basename "$db_dir")"
        if [[ ! -f "$gamma_dir/.gamma_done" ]]; then
            log "  Gamma reopt ($n_partitions trees)..."
            mkdir -p "$gamma_dir"
            for f in "$part_dir"/partition*_MSA.fasta "$part_dir"/partition*_taxonomy.txt; do
                local fname=$(basename "$f")
                [[ ! -e "$gamma_dir/$fname" ]] && ln -s "$(cd "$(dirname "$f")" && pwd)/$fname" "$gamma_dir/$fname"
            done
            local gdone=0
            for tree_file in "$part_dir"/RAxML_bestTree.*.reroot; do
                local fname=$(basename "$tree_file")
                local pname="${fname#RAxML_bestTree.}"
                pname="${pname%.reroot}"
                local msa="$part_dir/${pname}_MSA.fasta"
                if [[ -f "$msa" ]]; then
                    OMP_NUM_THREADS=1 FastTree -nt -gtr -gamma -nome -mllen \
                        -intree "$tree_file" "$msa" \
                        > "$gamma_dir/$fname" 2>/dev/null
                    [[ ! -s "$gamma_dir/$fname" ]] && cp "$tree_file" "$gamma_dir/$fname"
                else
                    cp "$tree_file" "$gamma_dir/$fname"
                fi
                gdone=$((gdone + 1))
                (( gdone % 50 == 0 )) && log "    Gamma: $gdone/$n_partitions"
            done
            touch "$gamma_dir/.gamma_done"
        fi
        tronko_input="$gamma_dir"
    fi

    # tronko-build with norepartition
    log "  tronko-build ($n_partitions partitions, norepartition)..."
    local t0=$(date +%s)
    "$TRONKO_BUILD" -y \
        -e "$tronko_input" \
        -n "$n_partitions" \
        -d "$db_dir" \
        -E -a -f 999999 \
        -c "$THREADS"
    local t1=$(date +%s)

    if [[ ! -f "$ref_tree" ]]; then
        log "  ERROR: tronko-build failed for $label"
        return 1
    fi

    local n_trees=$(head -1 "$ref_tree")
    log "  tronko-build: $n_trees trees ($((t1-t0))s)"

    # marker.fasta + taxonomy + BWA index
    local marker="$db_dir/marker.fasta"
    local marker_tax="$db_dir/marker_taxonomy.txt"
    > "$marker"
    > "$marker_tax"
    for i in $(seq 0 $((n_partitions - 1))); do
        if [[ -f "$part_dir/partition${i}_unaligned.fasta" ]]; then
            cat "$part_dir/partition${i}_unaligned.fasta" >> "$marker"
        elif [[ -f "$part_dir/partition${i}_MSA.fasta" ]]; then
            awk '/^>/{print; next} {gsub(/-/,""); print}' "$part_dir/partition${i}_MSA.fasta" >> "$marker"
        fi
        [[ -f "$part_dir/partition${i}_taxonomy.txt" ]] && cat "$part_dir/partition${i}_taxonomy.txt" >> "$marker_tax"
    done
    bwa index "$marker" 2>/dev/null

    # Convert to .trkb (zstd compressed binary)
    if command -v tronko-convert &>/dev/null; then
        log "  tronko-convert..."
        tronko-convert -i "$ref_tree" -o "$db_dir/reference_tree.trkb" -c zstd 2>/dev/null
    fi

    # Provenance
    cp "$FILTERED_FASTA" "$db_dir/input.fasta" 2>/dev/null || true
    cp "$FILTERED_TAX" "$db_dir/input_taxonomy.txt" 2>/dev/null || true

    local size_mb=$(du -m "$ref_tree" | cut -f1)
    cat > "$db_dir/build_info.json" <<EOF
{
  "name": "$label",
  "pipeline": "partition_and_build.py + tronko-build (norepartition)",
  "marker": "12SV5",
  "input_sequences": 101088,
  "n_trees": $n_trees,
  "n_input_partitions": $n_partitions,
  "db_size_mb": $size_mb,
  "sp_partitioning": false,
  "internal_repartitioning": "none",
  "gamma": $use_gamma,
  "tronko_build_flags": "-y -e <parts> -n $n_partitions -d <out> -E -a -f 999999 -c $THREADS",
  "build_seconds": $((t1-t0)),
  "date": "$(date +%Y-%m-%d)"
}
EOF

    log "  DONE: $label → $n_trees trees, ${size_mb} MB"
}

# =========================================================================
# Rebuild all configs with norepartition
# =========================================================================

# sub=200 (existing PASTA tree) decompositions
for decomp in maxsize500 maxsize1000 maxsize2000 maxsize5000 maxdiam10 maxdiam25; do
    part_dir="$SWEEP_DIR/_cache/partitions_${decomp}"
    [[ -d "$part_dir" ]] || continue
    for gamma in false true; do
        label="${decomp}_nogamma"
        [[ "$gamma" == "true" ]] && label="${decomp}_gamma"
        db_dir="$SWEEP_DIR/$label"
        log ""
        log "=== $label ==="
        build_one "$part_dir" "$db_dir" "$label" "$gamma"
    done
done

# sub=500 decompositions
for decomp in sub500_maxsize500 sub500_maxsize1000 sub500_maxsize2000 sub500_maxsize5000; do
    part_dir="$SWEEP_DIR/_cache/partitions_${decomp}"
    [[ -d "$part_dir" ]] || continue
    for gamma in false true; do
        label="${decomp}_nogamma"
        [[ "$gamma" == "true" ]] && label="${decomp}_gamma"
        db_dir="$SWEEP_DIR/$label"
        log ""
        log "=== $label ==="
        build_one "$part_dir" "$db_dir" "$label" "$gamma"
    done
done

# Also copy pasta_stats.json into new databases
ORIG_PASTA='{"pasta_version":"1.9.3","iter_limit":3,"max_subproblem_size":200,"aligner":"mafft","merger":"muscle","tree_estimator":"fasttree","num_cpus":4,"input_sequences":101088,"runtime_seconds":30429,"scores_per_iteration":[-5625060.805,-3418673.679,-3345837.877],"final_score":-3309385.351}'
SUB500_PASTA='{"pasta_version":"1.9.3","iter_limit":3,"max_subproblem_size":500,"aligner":"mafft","merger":"muscle","tree_estimator":"fasttree","num_cpus":8,"input_sequences":101088,"runtime_seconds":30112,"scores_per_iteration":[-4723503.39,-3235594.818,-3192804.36],"final_score":-3168988.702}'

for dir in "$SWEEP_DIR"/maxsize*/ "$SWEEP_DIR"/maxdiam*/; do
    [[ -d "$dir" ]] || continue
    echo "$ORIG_PASTA" | python3 -m json.tool > "$dir/pasta_stats.json"
done
for dir in "$SWEEP_DIR"/sub500_*/; do
    [[ -d "$dir" ]] || continue
    echo "$SUB500_PASTA" | python3 -m json.tool > "$dir/pasta_stats.json"
done

log ""
log "All rebuilds complete."
ls -d "$SWEEP_DIR"/*/ | grep -v _cache | wc -l | xargs -I{} echo "Total databases: {}"
