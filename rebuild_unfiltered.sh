#!/bin/bash
# Sequential rebuild: AncestralClust unfiltered, then PASTA unfiltered sweep
# Runs one at a time to give each the full 64 cores / 492GB RAM
#
# Usage:
#   nohup bash rebuild_unfiltered.sh > databases/rebuild_unfiltered.log 2>&1 &

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

# Add bundled tools and compiled binaries to PATH
export PATH="$REPO_DIR/bin:$REPO_DIR/tronko-build:$REPO_DIR/tronko-assign:$REPO_DIR/tronko-convert:$PATH"

THREADS=60  # leave a few cores for OS overhead
PARALLEL_JOBS=1  # Step 4 -J parallelism causes segfaults; use J=1 for safety, all 60 threads for posteriors

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') $*"
}

# =========================================================================
# PART 1: Clean rebuild of AncestralClust unfiltered 12sv5
# =========================================================================
log "============================================================"
log "PART 1: AncestralClust unfiltered 12sv5 — CLEAN REBUILD"
log "============================================================"

UNFILT_DIR="$REPO_DIR/databases/12sv5_unfiltered"
INPUT_FASTA="$UNFILT_DIR/input.fasta"
INPUT_TAX="$UNFILT_DIR/input_taxonomy.txt"
OUTPUT_DIR="$UNFILT_DIR/tronko_db"

# Verify inputs exist
if [[ ! -f "$INPUT_FASTA" ]] || [[ ! -f "$INPUT_TAX" ]]; then
    log "ERROR: input.fasta or input_taxonomy.txt missing in $UNFILT_DIR"
    exit 1
fi

log "Input sequences: $(grep -c '^>' "$INPUT_FASTA")"
log "Cleaning old build artifacts..."

# Remove everything except input files and build cache
cd "$UNFILT_DIR"
find . -maxdepth 1 \
    ! -name '.' \
    ! -name 'input.fasta' \
    ! -name 'input_taxonomy.txt' \
    ! -name '.cache' \
    -exec rm -rf {} +
cd "$REPO_DIR"

log "Old artifacts removed. Starting clean build..."
log "Params: -T $THREADS -s 0.1 -F -B 50000 -P 30 -J $PARALLEL_JOBS"

t0=$(date +%s)

bash "$REPO_DIR/build-tronko-db.sh" \
    -f "$INPUT_FASTA" \
    -t "$INPUT_TAX" \
    -o "$UNFILT_DIR" \
    -T "$THREADS" \
    -s 0.1 \
    -F \
    -E \
    -B 50000 \
    -P 30 \
    -J "$PARALLEL_JOBS"

t1=$(date +%s)

if [[ -f "$UNFILT_DIR/reference_tree.txt" ]]; then
    n_trees=$(head -1 "$UNFILT_DIR/reference_tree.txt")
    size=$(du -h "$UNFILT_DIR/reference_tree.txt" | cut -f1)
    log "PART 1 COMPLETE: $n_trees trees, $size, $((t1-t0))s"
else
    log "PART 1 FAILED: no reference_tree.txt produced"
    log "Continuing to Part 2 anyway..."
fi

# =========================================================================
# PART 2: PASTA unfiltered sweep (Phase 3 from run_overnight.sh)
# =========================================================================
log ""
log "============================================================"
log "PART 2: PASTA unfiltered 176K sweep"
log "============================================================"

# Activate PASTA env
source /tmp/pasta_test_env/bin/activate 2>/dev/null || true
export PYTHONPATH="$HOME/pasta:${PYTHONPATH:-}"

TRONKO_BUILD="$REPO_DIR/tronko-build/tronko-build"
PARTITION_SCRIPT="$REPO_DIR/tronko-build/partition_and_build.py"
PASTA_RUNNER="$HOME/pasta/run_pasta.py"
SWEEP_DIR="$REPO_DIR/databases/pasta_sweep"
PASTA_THREADS=$THREADS

mkdir -p "$SWEEP_DIR"

# Clear failed PASTA cache for unfiltered
PASTA_CACHE="$SWEEP_DIR/_cache/pasta_iter3_sub500_unfiltered"
if [[ -d "$PASTA_CACHE" ]] && [[ ! -f "$PASTA_CACHE/.pasta_done" ]]; then
    log "Clearing failed PASTA cache: $PASTA_CACHE"
    rm -rf "$PASTA_CACHE"
fi

# --- Run PASTA ---
PASTA_DIR="$SWEEP_DIR/_cache/pasta_iter3_sub500_unfiltered"
PASTA_TREE="$PASTA_DIR/pasta_iter3_sub500_unfiltered_final.tre"

if [[ -f "$PASTA_DIR/.pasta_done" ]]; then
    log "PASTA CACHED: $PASTA_DIR"
else
    mkdir -p "$PASTA_DIR"
    log "Running PASTA (iter=3, sub=500) on 176K sequences..."
    t2=$(date +%s)

    python3 "$PASTA_RUNNER" \
        --input "$INPUT_FASTA" \
        --datatype dna \
        --iter-limit 3 \
        --max-subproblem-size 500 \
        --num-cpus "$PASTA_THREADS" \
        --merger muscle \
        --job "pasta_iter3_sub500_unfiltered" \
        --output-directory "$PASTA_DIR" \
        -o "$PASTA_DIR" \
        2>&1 | tee "$PASTA_DIR/pasta_run.log"

    t3=$(date +%s)

    # Find and fix the output tree (strip quotes)
    for f in "$PASTA_DIR"/*.tre; do
        [[ "$f" == *temp* ]] && continue
        [[ "$f" == *topo* ]] && continue
        [[ "$f" == *final* ]] && continue
        sed "s/'//g" "$f" > "$PASTA_TREE"
        break
    done

    if [[ -f "$PASTA_TREE" ]] && [[ -s "$PASTA_TREE" ]]; then
        echo "$((t3-t2))" > "$PASTA_DIR/pasta_elapsed.txt"
        touch "$PASTA_DIR/.pasta_done"
        log "PASTA done ($((t3-t2))s): $PASTA_TREE"
    else
        log "ERROR: PASTA did not produce a tree"
        exit 1
    fi
fi

# --- Decompose + build for each config ---
if [[ -f "$PASTA_TREE" ]]; then
    for max_size in 500 1000 2000 5000; do
        for gamma in false true; do
            label="unfiltered_maxsize${max_size}"
            [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
            part_dir="$SWEEP_DIR/_cache/partitions_unfiltered_maxsize${max_size}"
            db_dir="$SWEEP_DIR/$label"

            log ""
            log "--- $label ---"

            # Decomposition
            if [[ -f "$part_dir/.decomp_done" ]]; then
                log "  Decomposition CACHED: $part_dir"
            else
                log "  Decomposing (--max-size $max_size)..."
                python3 "$PARTITION_SCRIPT" \
                    --tree "$PASTA_TREE" \
                    --fasta "$INPUT_FASTA" \
                    --taxonomy "$INPUT_TAX" \
                    --outdir "$part_dir" \
                    --max-size "$max_size" \
                    --min-size 3 \
                    --strategy centroid \
                    --threads "$PASTA_THREADS" \
                    --famsa-threads 2 \
                    --workers 4 \
                    --skip-tronko
                touch "$part_dir/.decomp_done"
            fi

            # Check if DB already exists
            ref_tree="$db_dir/reference_tree.txt"
            if [[ -f "$ref_tree" ]]; then
                n_trees=$(head -1 "$ref_tree")
                log "  CACHED: $label ($n_trees trees)"
                continue
            fi

            mkdir -p "$db_dir"
            n_partitions=$(ls "$part_dir"/RAxML_bestTree.*.reroot 2>/dev/null | wc -l | tr -d ' ')

            # Gamma reoptimization
            tronko_input="$part_dir"
            if [[ "$gamma" == "true" ]]; then
                gamma_dir="${db_dir}/_gamma_trees"
                if [[ ! -f "$gamma_dir/.gamma_done" ]]; then
                    log "  Gamma reoptimization ($n_partitions trees)..."
                    mkdir -p "$gamma_dir"
                    for f in "$part_dir"/partition*_MSA.fasta "$part_dir"/partition*_taxonomy.txt; do
                        fname=$(basename "$f")
                        [[ ! -e "$gamma_dir/$fname" ]] && ln -s "$(cd "$(dirname "$f")" && pwd)/$fname" "$gamma_dir/$fname"
                    done
                    gdone=0
                    for tree_file in "$part_dir"/RAxML_bestTree.*.reroot; do
                        fname=$(basename "$tree_file")
                        pname="${fname#RAxML_bestTree.}"
                        pname="${pname%.reroot}"
                        msa="$part_dir/${pname}_MSA.fasta"
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
            log "  tronko-build ($n_partitions partitions)..."
            tb_t0=$(date +%s)
            "$TRONKO_BUILD" -y \
                -e "$tronko_input" \
                -n "$n_partitions" \
                -d "$db_dir" \
                -E -a \
                -c "$PASTA_THREADS"
            tb_t1=$(date +%s)

            if [[ -f "$ref_tree" ]]; then
                n_trees=$(head -1 "$ref_tree")
                log "  tronko-build done: $n_trees trees ($((tb_t1-tb_t0))s)"
            else
                log "  ERROR: tronko-build failed for $label"
                continue
            fi

            # marker.fasta
            log "  Building marker.fasta..."
            marker="$db_dir/marker.fasta"
            > "$marker"
            for i in $(seq 0 $((n_partitions - 1))); do
                if [[ -f "$part_dir/partition${i}_unaligned.fasta" ]]; then
                    cat "$part_dir/partition${i}_unaligned.fasta" >> "$marker"
                elif [[ -f "$part_dir/partition${i}_MSA.fasta" ]]; then
                    awk '/^>/{print; next} {gsub(/-/,""); print}' "$part_dir/partition${i}_MSA.fasta" >> "$marker"
                fi
            done
            bwa index "$marker" 2>/dev/null

            # build_info.json
            size_mb=$(du -m "$ref_tree" | cut -f1)
            cat > "$db_dir/build_info.json" <<EOF
{
  "name": "$label",
  "n_trees": $n_trees,
  "db_size_mb": $size_mb,
  "n_input_partitions": $n_partitions,
  "gamma": $gamma,
  "tronko_build_flags": "-y -e <parts> -n $n_partitions -d <out> -E -a -c $PASTA_THREADS",
  "build_seconds": $((tb_t1-tb_t0)),
  "date": "$(date +%Y-%m-%d)"
}
EOF
            log "  DONE: $label → $n_trees trees, ${size_mb} MB"
        done
    done

    # norepartition variants
    for gamma in false true; do
        label="unfiltered_maxsize1000_norepartition"
        [[ "$gamma" == "true" ]] && label="${label}_gamma" || label="${label}_nogamma"
        part_dir="$SWEEP_DIR/_cache/partitions_unfiltered_maxsize1000"
        db_dir="$SWEEP_DIR/$label"

        ref_tree="$db_dir/reference_tree.txt"
        if [[ -f "$ref_tree" ]]; then
            n_trees=$(head -1 "$ref_tree")
            log "  CACHED: $label ($n_trees trees)"
            continue
        fi

        log ""
        log "--- $label ---"
        mkdir -p "$db_dir"
        n_partitions=$(ls "$part_dir"/RAxML_bestTree.*.reroot 2>/dev/null | wc -l | tr -d ' ')

        tronko_input="$part_dir"
        if [[ "$gamma" == "true" ]]; then
            gamma_dir="${db_dir}/_gamma_trees"
            if [[ ! -f "$gamma_dir/.gamma_done" ]]; then
                log "  Gamma reoptimization ($n_partitions trees)..."
                mkdir -p "$gamma_dir"
                for f in "$part_dir"/partition*_MSA.fasta "$part_dir"/partition*_taxonomy.txt; do
                    fname=$(basename "$f")
                    [[ ! -e "$gamma_dir/$fname" ]] && ln -s "$(cd "$(dirname "$f")" && pwd)/$fname" "$gamma_dir/$fname"
                done
                gdone=0
                for tree_file in "$part_dir"/RAxML_bestTree.*.reroot; do
                    fname=$(basename "$tree_file")
                    pname="${fname#RAxML_bestTree.}"
                    pname="${pname%.reroot}"
                    msa="$part_dir/${pname}_MSA.fasta"
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

        log "  tronko-build ($n_partitions partitions, no repartition)..."
        tb_t0=$(date +%s)
        "$TRONKO_BUILD" -y \
            -e "$tronko_input" \
            -n "$n_partitions" \
            -d "$db_dir" \
            -E -a -f 999999 \
            -c "$PASTA_THREADS"
        tb_t1=$(date +%s)

        if [[ -f "$ref_tree" ]]; then
            n_trees=$(head -1 "$ref_tree")
            log "  tronko-build done: $n_trees trees ($((tb_t1-tb_t0))s)"
        else
            log "  ERROR: tronko-build failed for $label"
            continue
        fi

        marker="$db_dir/marker.fasta"
        > "$marker"
        for i in $(seq 0 $((n_partitions - 1))); do
            if [[ -f "$part_dir/partition${i}_unaligned.fasta" ]]; then
                cat "$part_dir/partition${i}_unaligned.fasta" >> "$marker"
            elif [[ -f "$part_dir/partition${i}_MSA.fasta" ]]; then
                awk '/^>/{print; next} {gsub(/-/,""); print}' "$part_dir/partition${i}_MSA.fasta" >> "$marker"
            fi
        done
        bwa index "$marker" 2>/dev/null

        size_mb=$(du -m "$ref_tree" | cut -f1)
        cat > "$db_dir/build_info.json" <<EOF
{
  "name": "$label",
  "n_trees": $n_trees,
  "db_size_mb": $size_mb,
  "n_input_partitions": $n_partitions,
  "no_repartition": true,
  "gamma": $gamma,
  "tronko_build_flags": "-y -e <parts> -n $n_partitions -d <out> -E -a -f 999999 -c $PASTA_THREADS",
  "build_seconds": $((tb_t1-tb_t0)),
  "date": "$(date +%Y-%m-%d)"
}
EOF
        log "  DONE: $label → $n_trees trees, ${size_mb} MB"
    done
else
    log "ERROR: No PASTA tree available, skipping Part 2"
fi

log ""
log "============================================================"
log "ALL BUILDS COMPLETE"
log "============================================================"
log "AncestralClust DB: $UNFILT_DIR/reference_tree.txt"
log "PASTA sweep DBs:   $SWEEP_DIR/unfiltered_*/"
ls -d "$SWEEP_DIR"/unfiltered_*/ 2>/dev/null | wc -l | xargs -I{} echo "Total PASTA sweep DBs: {}"
