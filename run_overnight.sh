#!/bin/bash
# Overnight database build sweep
#
# Phase 1: Finish remaining decomposition configs with existing PASTA tree (filtered 101K)
# Phase 2: New PASTA run (iter=5, sub=500) on filtered 101K, then decompose
# Phase 3: PASTA run on unfiltered 176K, then decompose
#
# All databases written to databases/pasta_sweep/
#
# Usage:
#   nohup bash run_overnight.sh > databases/overnight.log 2>&1 &

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

# Activate PASTA env
source /tmp/pasta_test_env/bin/activate 2>/dev/null || true
export PYTHONPATH="$HOME/pasta:${PYTHONPATH:-}"

TRONKO_BUILD="$REPO_DIR/tronko-build/tronko-build"
PARTITION_SCRIPT="$REPO_DIR/tronko-build/partition_and_build.py"
PASTA_RUNNER="$HOME/pasta/run_pasta.py"
SWEEP_DIR="$REPO_DIR/databases/pasta_sweep"
THREADS=8

# Filtered 101K inputs
FILTERED_FASTA="$REPO_DIR/databases/vert12S_pasta/input.fasta"
FILTERED_TAX="$REPO_DIR/databases/vert12S_pasta/input_taxonomy.txt"
EXISTING_TREE="$REPO_DIR/databases/vert12S_pasta/input_tree.tre"

# Unfiltered 176K inputs
UNFILTERED_FASTA="$REPO_DIR/databases/12sv5_unfiltered/input.fasta"
UNFILTERED_TAX="$REPO_DIR/databases/12sv5_unfiltered/input_taxonomy.txt"

mkdir -p "$SWEEP_DIR"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') $*"
}

# =========================================================================
# Helper: build tronko DB from partitions
# =========================================================================
build_tronko_db() {
    local part_dir="$1"
    local db_dir="$2"
    local label="$3"
    local no_repartition="${4:-0}"
    local use_gamma="${5:-false}"
    local input_fasta="${6:-$FILTERED_FASTA}"
    local input_tax="${7:-$FILTERED_TAX}"

    local ref_tree="$db_dir/reference_tree.txt"

    if [[ -f "$ref_tree" ]]; then
        local n_trees=$(head -1 "$ref_tree")
        log "  CACHED: $label ($n_trees trees)"
        return 0
    fi

    mkdir -p "$db_dir"
    local n_partitions=$(ls "$part_dir"/RAxML_bestTree.*.reroot 2>/dev/null | wc -l | tr -d ' ')

    # Gamma reoptimization
    local tronko_input="$part_dir"
    if [[ "$use_gamma" == "true" ]]; then
        local gamma_dir="${db_dir}/_gamma_trees"
        if [[ ! -f "$gamma_dir/.gamma_done" ]]; then
            log "  Gamma reoptimization ($n_partitions trees)..."
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

    # tronko-build
    local extra_flags="-a"
    if [[ "$no_repartition" -eq 1 ]]; then
        extra_flags="-a -f 999999"
        log "  tronko-build ($n_partitions partitions, no repartition)..."
    else
        log "  tronko-build ($n_partitions partitions)..."
    fi

    local t0=$(date +%s)
    "$TRONKO_BUILD" -y \
        -e "$tronko_input" \
        -n "$n_partitions" \
        -d "$db_dir" \
        -E $extra_flags \
        -c "$THREADS"
    local t1=$(date +%s)

    if [[ ! -f "$ref_tree" ]]; then
        log "  ERROR: tronko-build failed for $label"
        return 1
    fi

    local n_trees=$(head -1 "$ref_tree")
    log "  tronko-build done: $n_trees trees ($((t1-t0))s)"

    # marker.fasta
    log "  Building marker.fasta..."
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
    cp "$input_fasta" "$db_dir/input.fasta" 2>/dev/null || true
    cp "$input_tax" "$db_dir/input_taxonomy.txt" 2>/dev/null || true

    # build_info.json
    local size_mb=$(du -m "$ref_tree" | cut -f1)
    cat > "$db_dir/build_info.json" <<EOF
{
  "name": "$label",
  "n_trees": $n_trees,
  "db_size_mb": $size_mb,
  "n_input_partitions": $n_partitions,
  "no_repartition": $( [[ "$no_repartition" -eq 1 ]] && echo "true" || echo "false" ),
  "gamma": $use_gamma,
  "tronko_build_flags": "-y -e <parts> -n $n_partitions -d <out> -E $extra_flags -c $THREADS",
  "build_seconds": $((t1-t0)),
  "date": "$(date +%Y-%m-%d)"
}
EOF

    log "  DONE: $label → $n_trees trees, ${size_mb} MB"
}

# =========================================================================
# Helper: run decomposition (partition_and_build.py --skip-tronko)
# =========================================================================
run_decomposition() {
    local tree="$1"
    local fasta="$2"
    local taxonomy="$3"
    local outdir="$4"
    local decomp_args="$5"  # e.g. "--max-size 1000"

    if [[ -f "$outdir/.decomp_done" ]]; then
        log "  Decomposition CACHED: $outdir"
        return 0
    fi

    log "  Decomposing ($decomp_args)..."
    python3 "$PARTITION_SCRIPT" \
        --tree "$tree" \
        --fasta "$fasta" \
        --taxonomy "$taxonomy" \
        --outdir "$outdir" \
        $decomp_args \
        --min-size 3 \
        --strategy centroid \
        --threads "$THREADS" \
        --famsa-threads 2 \
        --workers 4 \
        --skip-tronko

    touch "$outdir/.decomp_done"
}

# =========================================================================
# Helper: run PASTA
# =========================================================================
run_pasta() {
    local input_fasta="$1"
    local outdir="$2"
    local job_name="$3"
    local iter_limit="$4"
    local max_sub="$5"

    local done_marker="$outdir/.pasta_done"
    if [[ -f "$done_marker" ]]; then
        log "  PASTA CACHED: $outdir"
        return 0
    fi

    mkdir -p "$outdir"
    log "  Running PASTA (iter=$iter_limit, sub=$max_sub)..."
    local t0=$(date +%s)

    python3 "$PASTA_RUNNER" \
        --input "$input_fasta" \
        --datatype dna \
        --iter-limit "$iter_limit" \
        --max-subproblem-size "$max_sub" \
        --num-cpus "$THREADS" \
        --merger muscle \
        --job "$job_name" \
        --output-directory "$outdir" \
        -o "$outdir" \
        2>&1 | tee "$outdir/pasta_run.log"

    local t1=$(date +%s)

    # Find and fix the output tree (strip quotes)
    local tree_out="$outdir/${job_name}_final.tre"
    for f in "$outdir"/*.tre; do
        [[ "$f" == *temp* ]] && continue
        [[ "$f" == *topo* ]] && continue
        [[ "$f" == *final* ]] && continue
        # Copy and strip quotes
        sed "s/'//g" "$f" > "$tree_out"
        break
    done

    if [[ ! -f "$tree_out" ]]; then
        log "  ERROR: PASTA did not produce a tree"
        return 1
    fi

    echo "$((t1-t0))" > "$outdir/pasta_elapsed.txt"
    touch "$done_marker"
    log "  PASTA done ($((t1-t0))s): $tree_out"
}


# =========================================================================
# PHASE 1: Remaining configs with existing PASTA tree (filtered 101K)
# =========================================================================
log "============================================================"
log "PHASE 1: Decomposition sweep with existing PASTA tree"
log "============================================================"

# The norepartition config — replicates the original 141-tree winning DB
for gamma in false true; do
    label="maxsize1000_norepartition"
    [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
    part_dir="$SWEEP_DIR/_cache/partitions_maxsize1000"
    db_dir="$SWEEP_DIR/$label"

    log ""
    log "--- $label ---"

    run_decomposition "$EXISTING_TREE" "$FILTERED_FASTA" "$FILTERED_TAX" \
        "$part_dir" "--max-size 1000"

    build_tronko_db "$part_dir" "$db_dir" "$label" 1 "$gamma"
done

log ""
log "Phase 1 complete."


# =========================================================================
# PHASE 2: New PASTA (iter=5, sub=500) on filtered 101K
# =========================================================================
log ""
log "============================================================"
log "PHASE 2: PASTA iter=3 sub=500 on filtered 101K"
log "============================================================"

PASTA_DIR_FILTERED="$SWEEP_DIR/_cache/pasta_iter3_sub500_filtered"
PASTA_TREE_FILTERED="$PASTA_DIR_FILTERED/pasta_iter3_sub500_filtered_final.tre"

run_pasta "$FILTERED_FASTA" "$PASTA_DIR_FILTERED" "pasta_iter3_sub500_filtered" 3 500

if [[ -f "$PASTA_TREE_FILTERED" ]]; then
    # Decompose with same configs that worked well in Phase 1
    for max_size in 500 1000 2000 5000; do
        for gamma in false true; do
            label="sub500_maxsize${max_size}"
            [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
            part_dir="$SWEEP_DIR/_cache/partitions_sub500_maxsize${max_size}"
            db_dir="$SWEEP_DIR/$label"

            log ""
            log "--- $label ---"

            run_decomposition "$PASTA_TREE_FILTERED" "$FILTERED_FASTA" "$FILTERED_TAX" \
                "$part_dir" "--max-size $max_size"

            build_tronko_db "$part_dir" "$db_dir" "$label" 0 "$gamma"
        done
    done

    # maxdiam configs
    for max_diam in 10 25; do
        for gamma in false true; do
            label="sub500_maxdiam${max_diam}"
            [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
            part_dir="$SWEEP_DIR/_cache/partitions_sub500_maxdiam${max_diam}"
            db_dir="$SWEEP_DIR/$label"

            log ""
            log "--- $label ---"

            run_decomposition "$PASTA_TREE_FILTERED" "$FILTERED_FASTA" "$FILTERED_TAX" \
                "$part_dir" "--max-diam $max_diam"

            build_tronko_db "$part_dir" "$db_dir" "$label" 0 "$gamma"
        done
    done

    # Also try norepartition
    for gamma in false true; do
        label="sub500_maxsize1000_norepartition"
        [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
        part_dir="$SWEEP_DIR/_cache/partitions_sub500_maxsize1000"
        db_dir="$SWEEP_DIR/$label"

        log ""
        log "--- $label ---"
        build_tronko_db "$part_dir" "$db_dir" "$label" 1 "$gamma"
    done
else
    log "PASTA failed for filtered, skipping Phase 2 decompositions"
fi

log ""
log "Phase 2 complete."


# =========================================================================
# PHASE 3: PASTA on unfiltered 176K
# =========================================================================
log ""
log "============================================================"
log "PHASE 3: PASTA iter=3 sub=500 on unfiltered 176K"
log "============================================================"

PASTA_DIR_UNFILTERED="$SWEEP_DIR/_cache/pasta_iter3_sub500_unfiltered"
PASTA_TREE_UNFILTERED="$PASTA_DIR_UNFILTERED/pasta_iter3_sub500_unfiltered_final.tre"

run_pasta "$UNFILTERED_FASTA" "$PASTA_DIR_UNFILTERED" "pasta_iter3_sub500_unfiltered" 3 500

if [[ -f "$PASTA_TREE_UNFILTERED" ]]; then
    for max_size in 500 1000 2000 5000; do
        for gamma in false true; do
            label="unfiltered_maxsize${max_size}"
            [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
            part_dir="$SWEEP_DIR/_cache/partitions_unfiltered_maxsize${max_size}"
            db_dir="$SWEEP_DIR/$label"

            log ""
            log "--- $label ---"

            run_decomposition "$PASTA_TREE_UNFILTERED" "$UNFILTERED_FASTA" "$UNFILTERED_TAX" \
                "$part_dir" "--max-size $max_size"

            build_tronko_db "$part_dir" "$db_dir" "$label" 0 "$gamma" \
                "$UNFILTERED_FASTA" "$UNFILTERED_TAX"
        done
    done

    # maxdiam configs
    for max_diam in 10 25; do
        for gamma in false true; do
            label="unfiltered_maxdiam${max_diam}"
            [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
            part_dir="$SWEEP_DIR/_cache/partitions_unfiltered_maxdiam${max_diam}"
            db_dir="$SWEEP_DIR/$label"

            log ""
            log "--- $label ---"

            run_decomposition "$PASTA_TREE_UNFILTERED" "$UNFILTERED_FASTA" "$UNFILTERED_TAX" \
                "$part_dir" "--max-diam $max_diam"

            build_tronko_db "$part_dir" "$db_dir" "$label" 0 "$gamma" \
                "$UNFILTERED_FASTA" "$UNFILTERED_TAX"
        done
    done

    # norepartition variant
    for gamma in false true; do
        label="unfiltered_maxsize1000_norepartition"
        [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
        part_dir="$SWEEP_DIR/_cache/partitions_unfiltered_maxsize1000"
        db_dir="$SWEEP_DIR/$label"

        log ""
        log "--- $label ---"
        build_tronko_db "$part_dir" "$db_dir" "$label" 1 "$gamma" \
            "$UNFILTERED_FASTA" "$UNFILTERED_TAX"
    done
else
    log "PASTA failed for unfiltered, skipping Phase 3 decompositions"
fi

log ""
log "============================================================"
log "ALL PHASES COMPLETE"
log "============================================================"
log "Databases in: $SWEEP_DIR"
ls -d "$SWEEP_DIR"/*/  2>/dev/null | wc -l | xargs -I{} echo "Total database dirs: {}"
