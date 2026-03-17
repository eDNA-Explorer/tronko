#!/bin/bash
# Sweep PASTA decomposition parameters using an existing PASTA tree.
# Reuses the same tree for all configs — varies partition method, size, and gamma.
#
# Two decomposition modes:
#   max_size: cap number of leaves per partition
#   max_diam: cap tree diameter per partition
#
# Usage:
#   bash sweep_pasta_decomp.sh
#   bash sweep_pasta_decomp.sh --dry-run   # show what would run

set -euo pipefail

# =========================================================================
# Configuration
# =========================================================================
PASTA_TREE="databases/vert12S_pasta/input_tree.tre"
INPUT_FASTA="databases/vert12S_pasta/input.fasta"
INPUT_TAXONOMY="databases/vert12S_pasta/input_taxonomy.txt"
TRONKO_BUILD="tronko-build/tronko-build"
PARTITION_SCRIPT="tronko-build/partition_and_build.py"
SWEEP_DIR="databases/pasta_sweep"

THREADS=6
FAMSA_THREADS=2
WORKERS=3
MIN_SIZE=3
STRATEGY="centroid"

# Parameters to sweep
# Format: "mode:value[:norepartition]"
# norepartition = pass -f 999999 to tronko-build to skip internal 3-way splitting
DECOMP_CONFIGS=(
    "size:500"
    "size:1000"
    "size:2000"
    "size:5000"
    "diam:10"
    "diam:25"
    "size:1000:norepartition"
)
GAMMA_OPTIONS=(false true)

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
fi

# =========================================================================
# Validate inputs
# =========================================================================
for f in "$PASTA_TREE" "$INPUT_FASTA" "$INPUT_TAXONOMY" "$TRONKO_BUILD" "$PARTITION_SCRIPT"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
        exit 1
    fi
done

# Activate PASTA virtualenv if needed
if [[ -d "/tmp/pasta_test_env" ]]; then
    source /tmp/pasta_test_env/bin/activate
fi

mkdir -p "$SWEEP_DIR"

# =========================================================================
# Results tracking
# =========================================================================
RESULTS_FILE="$SWEEP_DIR/sweep_results.json"
echo "[" > "$RESULTS_FILE"
FIRST_RESULT=1

add_result() {
    local name="$1" status="$2" n_trees="$3" size_mb="$4" n_partitions="$5" elapsed="$6"
    if [[ "$FIRST_RESULT" -eq 0 ]]; then
        echo "," >> "$RESULTS_FILE"
    fi
    FIRST_RESULT=0
    cat >> "$RESULTS_FILE" <<JSONEOF
  {
    "name": "$name",
    "status": "$status",
    "n_trees": $n_trees,
    "size_mb": $size_mb,
    "n_partitions": $n_partitions,
    "elapsed_seconds": $elapsed
  }
JSONEOF
}

# =========================================================================
# Sweep loop
# =========================================================================
TOTAL=$(( ${#DECOMP_CONFIGS[@]} * ${#GAMMA_OPTIONS[@]} ))
CURRENT=0
SWEEP_START=$(date +%s)

for DECOMP_CONFIG in "${DECOMP_CONFIGS[@]}"; do
    # Parse mode:value[:norepartition]
    IFS=':' read -r DECOMP_MODE DECOMP_VALUE DECOMP_EXTRA <<< "$DECOMP_CONFIG"
    NO_REPARTITION=0
    if [[ "${DECOMP_EXTRA:-}" == "norepartition" ]]; then
        NO_REPARTITION=1
    fi

    for USE_GAMMA in "${GAMMA_OPTIONS[@]}"; do
        CURRENT=$((CURRENT + 1))

        if [[ "$DECOMP_MODE" == "size" ]]; then
            LABEL="maxsize${DECOMP_VALUE}"
        else
            LABEL="maxdiam${DECOMP_VALUE}"
        fi
        if [[ "$NO_REPARTITION" -eq 1 ]]; then
            LABEL="${LABEL}_norepartition"
        fi

        if [[ "$USE_GAMMA" == "true" ]]; then
            NAME="${LABEL}_gamma"
        else
            NAME="${LABEL}_nogamma"
        fi

        DB_DIR="$SWEEP_DIR/$NAME"
        PART_DIR="$SWEEP_DIR/_cache/partitions_${LABEL}"
        REF_TREE="$DB_DIR/reference_tree.txt"

        echo ""
        echo "================================================================"
        echo "[$CURRENT/$TOTAL] $NAME"
        echo "  ${DECOMP_MODE}=${DECOMP_VALUE}, gamma=$USE_GAMMA"
        echo "================================================================"

        # Check if already done
        if [[ -f "$REF_TREE" ]]; then
            N_TREES=$(head -1 "$REF_TREE")
            SIZE_MB=$(du -m "$REF_TREE" | cut -f1)
            echo "  CACHED: $N_TREES trees, ${SIZE_MB} MB"
            add_result "$NAME" "cached" "$N_TREES" "$SIZE_MB" "0" "0"
            continue
        fi

        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "  [DRY RUN] Would build: $NAME"
            add_result "$NAME" "dry_run" "0" "0" "0" "0"
            continue
        fi

        CONFIG_START=$(date +%s)
        mkdir -p "$DB_DIR"

        # -----------------------------------------------------------
        # Stage 2: Decomposition (cached across gamma variants)
        # -----------------------------------------------------------
        if [[ -f "$PART_DIR/.decomp_done" ]]; then
            echo "  Decomposition CACHED: $PART_DIR"
        else
            echo "  Decomposing tree (${DECOMP_MODE}=${DECOMP_VALUE})..."
            DECOMP_ARGS=()
            if [[ "$DECOMP_MODE" == "size" ]]; then
                DECOMP_ARGS+=(--max-size "$DECOMP_VALUE")
            else
                DECOMP_ARGS+=(--max-diam "$DECOMP_VALUE")
            fi

            python3 "$PARTITION_SCRIPT" \
                --tree "$PASTA_TREE" \
                --fasta "$INPUT_FASTA" \
                --taxonomy "$INPUT_TAXONOMY" \
                --outdir "$PART_DIR" \
                "${DECOMP_ARGS[@]}" \
                --min-size "$MIN_SIZE" \
                --strategy "$STRATEGY" \
                --threads "$THREADS" \
                --famsa-threads "$FAMSA_THREADS" \
                --workers "$WORKERS" \
                --skip-tronko

            touch "$PART_DIR/.decomp_done"
        fi

        # Count partitions
        N_PARTITIONS=$(ls "$PART_DIR"/RAxML_bestTree.*.reroot 2>/dev/null | wc -l | tr -d ' ')
        echo "  $N_PARTITIONS partitions"

        # -----------------------------------------------------------
        # Stage 3: Gamma reoptimization (if enabled)
        # -----------------------------------------------------------
        TRONKO_INPUT_DIR="$PART_DIR"

        if [[ "$USE_GAMMA" == "true" ]]; then
            GAMMA_DIR="$SWEEP_DIR/_cache/gamma_${LABEL}"
            if [[ -f "$GAMMA_DIR/.gamma_done" ]]; then
                echo "  Gamma reopt CACHED: $GAMMA_DIR"
            else
                echo "  Reoptimizing branch lengths with -gamma..."
                mkdir -p "$GAMMA_DIR"

                # Symlink MSA and taxonomy files
                for f in "$PART_DIR"/partition*_MSA.fasta "$PART_DIR"/partition*_taxonomy.txt; do
                    fname=$(basename "$f")
                    [[ ! -e "$GAMMA_DIR/$fname" ]] && ln -s "$(cd "$(dirname "$f")" && pwd)/$fname" "$GAMMA_DIR/$fname"
                done

                # Re-run FastTree with -gamma on each tree
                GAMMA_DONE=0
                GAMMA_TOTAL=$N_PARTITIONS
                for tree_file in "$PART_DIR"/RAxML_bestTree.*.reroot; do
                    fname=$(basename "$tree_file")
                    partition_name="${fname#RAxML_bestTree.}"
                    partition_name="${partition_name%.reroot}"
                    msa="$PART_DIR/${partition_name}_MSA.fasta"

                    if [[ -f "$msa" ]]; then
                        # FastTree -gamma -nome -mllen -intree reoptimizes branch lengths
                        OMP_NUM_THREADS=1 FastTree -nt -gtr -gamma \
                            -nome -mllen \
                            -intree "$tree_file" \
                            "$msa" \
                            > "$GAMMA_DIR/$fname" 2>/dev/null

                        # If gamma failed, fall back to original
                        if [[ ! -s "$GAMMA_DIR/$fname" ]]; then
                            cp "$tree_file" "$GAMMA_DIR/$fname"
                        fi
                    else
                        cp "$tree_file" "$GAMMA_DIR/$fname"
                    fi

                    GAMMA_DONE=$((GAMMA_DONE + 1))
                    if (( GAMMA_DONE % 20 == 0 )); then
                        echo "    Gamma: $GAMMA_DONE/$GAMMA_TOTAL"
                    fi
                done
                echo "    Gamma: $GAMMA_DONE/$GAMMA_TOTAL done"
                touch "$GAMMA_DIR/.gamma_done"
            fi
            TRONKO_INPUT_DIR="$GAMMA_DIR"
        fi

        # -----------------------------------------------------------
        # Stage 4: tronko-build
        # -----------------------------------------------------------
        TRONKO_EXTRA_FLAGS="-a"
        if [[ "$NO_REPARTITION" -eq 1 ]]; then
            TRONKO_EXTRA_FLAGS="-a -f 999999"
            echo "  Running tronko-build ($N_PARTITIONS partitions, no repartition)..."
        else
            echo "  Running tronko-build ($N_PARTITIONS partitions)..."
        fi
        "$TRONKO_BUILD" -y \
            -e "$TRONKO_INPUT_DIR" \
            -n "$N_PARTITIONS" \
            -d "$DB_DIR" \
            -E $TRONKO_EXTRA_FLAGS \
            -c "$THREADS"

        if [[ ! -f "$REF_TREE" ]]; then
            echo "  ERROR: tronko-build did not produce reference_tree.txt"
            add_result "$NAME" "failed" "0" "0" "$N_PARTITIONS" "0"
            continue
        fi

        N_TREES=$(head -1 "$REF_TREE")
        echo "  tronko-build done: $N_TREES trees"

        # -----------------------------------------------------------
        # Stage 5: Build marker.fasta + taxonomy + BWA index
        # -----------------------------------------------------------
        echo "  Building marker.fasta..."
        MARKER="$DB_DIR/marker.fasta"
        MARKER_TAX="$DB_DIR/marker_taxonomy.txt"

        # Concatenate gap-free FASTAs from all partitions
        > "$MARKER"
        > "$MARKER_TAX"
        for i in $(seq 0 $((N_PARTITIONS - 1))); do
            # Try unaligned first (no gaps), fall back to MSA with gap stripping
            if [[ -f "$PART_DIR/partition${i}_unaligned.fasta" ]]; then
                cat "$PART_DIR/partition${i}_unaligned.fasta" >> "$MARKER"
            elif [[ -f "$PART_DIR/partition${i}_MSA.fasta" ]]; then
                awk '/^>/{print; next} {gsub(/-/,""); print}' "$PART_DIR/partition${i}_MSA.fasta" >> "$MARKER"
            fi
            if [[ -f "$PART_DIR/partition${i}_taxonomy.txt" ]]; then
                cat "$PART_DIR/partition${i}_taxonomy.txt" >> "$MARKER_TAX"
            fi
        done

        N_SEQS=$(grep -c '^>' "$MARKER")
        echo "  marker.fasta: $N_SEQS sequences"

        echo "  BWA indexing..."
        bwa index "$MARKER" 2>/dev/null

        # Copy input files for provenance
        cp "$INPUT_FASTA" "$DB_DIR/input.fasta" 2>/dev/null || true
        cp "$INPUT_TAXONOMY" "$DB_DIR/input_taxonomy.txt" 2>/dev/null || true

        CONFIG_END=$(date +%s)
        ELAPSED=$((CONFIG_END - CONFIG_START))
        SIZE_MB=$(du -m "$REF_TREE" | cut -f1)
        echo "  DONE: $N_TREES trees, ${SIZE_MB} MB, ${ELAPSED}s"
        add_result "$NAME" "success" "$N_TREES" "$SIZE_MB" "$N_PARTITIONS" "$ELAPSED"
    done
done

echo "" >> "$RESULTS_FILE"
echo "]" >> "$RESULTS_FILE"

SWEEP_END=$(date +%s)
TOTAL_ELAPSED=$((SWEEP_END - SWEEP_START))

echo ""
echo "================================================================"
echo "Sweep complete: $TOTAL configs in ${TOTAL_ELAPSED}s"
echo "Results: $RESULTS_FILE"
echo "================================================================"
cat "$RESULTS_FILE"
