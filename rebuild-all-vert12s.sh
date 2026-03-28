#!/bin/bash
set -euo pipefail

# ============================================================
# Rebuild all PASTA and AncestralClust tronko databases
# for a given marker, for both species and LCA taxonomy variants.
#
# Usage: bash rebuild-all.sh <MARKER>
#   e.g. bash rebuild-all.sh 12S_MiFish_U
#
# Caching strategy:
#   PASTA:  One PASTA tree per (marker, variant), reused across
#           all partition configs (maxsize/maxdiam).
#   AC:     One AncestralClust clustering per (marker, variant, bin_config),
#           reused across all SP thresholds.
#
# Set VARIANTS="lca" or VARIANTS="species" to build only one.
# Set SKIP_PASTA=1 or SKIP_AC=1 to skip a pipeline.
# ============================================================

MARKER="${1:?Usage: bash rebuild-all.sh <MARKER>}"

# ── Tunables ─────────────────────────────────────────────────
THREADS="${THREADS:-8}"
FAMSA_THREADS="${FAMSA_THREADS:-4}"
PARALLEL_JOBS="${PARALLEL_JOBS:-1}"
VARIANTS="${VARIANTS:-lca species}"
SKIP_PASTA="${SKIP_PASTA:-0}"
SKIP_AC="${SKIP_AC:-0}"

TREE_BACKEND="${TREE_BACKEND:-veryfasttree}"
TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko-fork}"
PASTA_DIR="$TRONKO_DIR/pasta"
export PASTA_DIR
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"
export PASTA_TREE_BACKEND="$TREE_BACKEND"

# Java for OPAL merger (PASTA)
if [[ -d /opt/homebrew/opt/openjdk ]]; then
    export JAVA_HOME="/opt/homebrew/opt/openjdk"
elif [[ -d /usr/lib/jvm/java-17-openjdk-amd64 ]]; then
    export JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"
elif [[ -d /usr/lib/jvm/java-21-openjdk-amd64 ]]; then
    export JAVA_HOME="/usr/lib/jvm/java-21-openjdk-amd64"
fi
export PATH="${JAVA_HOME:-}/bin:$PATH"

DB_BASE="databases/${MARKER}"

# ── Input files ──────────────────────────────────────────────
RCRUX_BASE="$HOME/rcrux-py/databases/${MARKER}/unfiltered"

# ── PASTA partition configs ──────────────────────────────────
# Each entry: "label --partition-flag value"
# Partitions are cached and reused across gamma/nogamma/SP variants.
PASTA_SIZES=(500 1000 2000 5000)
PASTA_DIAMS=(10 15 20 25 30)

# SP thresholds for repartitioning sweep (applied on top of each partition config)
PASTA_SP_THRESHOLDS=(2.00 2.25 2.50 2.75)

# Whether to build gamma and/or nogamma variants
PASTA_GAMMA=(true false)  # true=gamma, false=nogamma

# ── AC configs ───────────────────────────────────────────────
# Each entry: "label BIN_SIZE DESCENDANTS"
AC_CONFIGS=(
    "ac_default 10000 75"
    "ac_fewer_bins 20000 75"
    "ac_more_bins 5000 75"
)

# SP thresholds to sweep for each AC config
AC_SP_THRESHOLDS=(0.05 0.10 0.20)

# ============================================================
# Helpers
# ============================================================

log() { echo "[$(date '+%H:%M:%S')] $*"; }

get_input_files() {
    local variant="$1"
    if [[ "$variant" == "species" ]]; then
        INPUT_FASTA="${RCRUX_BASE}/${MARKER}_species.fasta"
        INPUT_TAX="${RCRUX_BASE}/${MARKER}_species_taxonomy.txt"
    else
        INPUT_FASTA="${RCRUX_BASE}/${MARKER}_lca.fasta"
        INPUT_TAX="${RCRUX_BASE}/${MARKER}_lca_taxonomy.txt"
    fi
    for f in "$INPUT_FASTA" "$INPUT_TAX"; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: $f not found" >&2; exit 1
        fi
    done
}

# ── PASTA tree builder ───────────────────────────────────────
# Builds one PASTA tree per (marker, variant). Skips if rooted tree exists.
build_pasta_tree() {
    local variant="$1"
    local pasta_out="$2"
    local job_name="${MARKER}_${variant}_pasta"
    local rooted_tree="${pasta_out}/${job_name}.rooted.tre"

    if [[ -s "$rooted_tree" ]]; then
        log "PASTA tree exists for $variant — skipping"
        ROOTED_TREE="$rooted_tree"
        return
    fi

    log "Running PASTA for $variant..."
    pkill -9 -f "run_pasta.*${job_name}" 2>/dev/null || true
    sleep 1
    rm -rf "$HOME/.pasta/${job_name}"
    mkdir -p "$pasta_out"

    python3 "$PASTA_DIR/run_pasta.py" \
        --input="$INPUT_FASTA" \
        --datatype=dna \
        --num-cpus="$THREADS" \
        --iter-limit=3 \
        -o "$pasta_out" \
        -j "$job_name" || true  # PASTA exits non-zero on harmless "Refused to clean" warning

    # Find the output tree
    local pasta_tree
    pasta_tree=$(find "$pasta_out" -name "${job_name}*.tre" ! -name "*_temp_*" -size +1k 2>/dev/null | head -1)
    if [[ -z "$pasta_tree" ]]; then
        echo "ERROR: PASTA did not produce a tree file" >&2; exit 1
    fi
    log "PASTA tree: $pasta_tree"

    # Restore original names from safe names, always quoting for nw_reroot
    local name_trans="${pasta_out}/${job_name}_temp_name_translation.txt"
    python3 -c "
import re, sys

trans = {}
with open('$name_trans') as f:
    lines = f.read().strip().split('\n')
    i = 0
    while i < len(lines):
        safe = lines[i].strip()
        if i+1 < len(lines):
            orig = lines[i+1].strip()
            if safe and orig:
                trans[safe] = orig
        i += 3

with open('$pasta_tree') as f:
    data = f.read()

data = data.replace(\"'\", '')
data = re.sub(r'\)([0-9][0-9.eE+-]*):', '):', data)

def replace_name(m):
    safe = m.group(1)
    return trans.get(safe, safe)

data = re.sub(r'([A-Za-z0-9_.-]+)(?=:)', replace_name, data)
data = re.sub(r'[;\s]+$', '', data) + ';'

with open('$pasta_tree', 'w') as f:
    f.write(data + '\n')

print(f'Restored {len(trans)} names')
"

    # Midpoint root
    nw_reroot "$pasta_tree" > "$rooted_tree"
    log "Rooted tree: $rooted_tree"
    ROOTED_TREE="$rooted_tree"
}

# ── Run tronko-build on cached partitions with given flags ───
# Reused by gamma/nogamma/SP variants on the same partitions.
run_tronko_on_partitions() {
    local part_dir="$1"
    local build_dir="$2"
    local n_partitions="$3"
    shift 3
    # remaining args are extra tronko-build flags

    if [[ -f "${build_dir}/reference_tree.txt" ]]; then
        log "    Already exists — skipping"
        return 0
    fi

    mkdir -p "$build_dir"
    tronko-build -y \
        -e "$part_dir" \
        -n "$n_partitions" \
        -d "$build_dir" \
        -a -E \
        "$@"

    cp "$INPUT_FASTA" "$build_dir/input.fasta"
    cp "$INPUT_TAX" "$build_dir/input_taxonomy.txt"
}

# ============================================================
# PASTA pipeline
#
# Cache hierarchy:
#   1. PASTA tree:   one per (marker, variant)
#   2. Partitions:   one per (marker, variant, size/diam) — FAMSA + FastTree
#   3. tronko-build: one per (partitions, gamma, SP) — fast, ~1-2 min each
# ============================================================
run_pasta_pipeline() {
    if [[ "$SKIP_PASTA" == "1" ]]; then
        log "Skipping PASTA pipeline (SKIP_PASTA=1)"
        return
    fi

    log "========== PASTA PIPELINE =========="

    for variant in $VARIANTS; do
        get_input_files "$variant"
        local pasta_dir="${DB_BASE}/${variant}/pasta"
        local pasta_out="${pasta_dir}/tree"
        mkdir -p "$pasta_dir"

        local seq_count
        seq_count=$(grep -c '^>' "$INPUT_FASTA")
        log "PASTA $variant: $seq_count sequences"

        # Step 1: Build PASTA tree (cached per variant)
        build_pasta_tree "$variant" "$pasta_out"

        # Step 2: For each partition config, build partitions (cached),
        #         then run tronko-build for each gamma/SP variant.
        for size in "${PASTA_SIZES[@]}"; do
            local part_label="maxsize${size}"
            local config_dir="${pasta_dir}/${part_label}"
            local part_dir="${config_dir}/partitions"

            # Build partitions if not cached
            if [[ ! -d "$part_dir" ]] || [[ -z "$(ls "$part_dir"/*_MSA.fasta 2>/dev/null)" ]]; then
                log "  Partitioning $part_label..."
                python3 "$TRONKO_DIR/tronko-build/partition_and_build.py" \
                    --tree "$ROOTED_TREE" \
                    --fasta "$INPUT_FASTA" \
                    --taxonomy "$INPUT_TAX" \
                    --outdir "$part_dir" \
                    --tronko-outdir "/dev/null" \
                    --strategy centroid \
                    --threads "$THREADS" \
                    --famsa-threads "$FAMSA_THREADS" \
                    --tronko-build "$(which tronko-build)" \
                    --skip-tronko \
                    --max-size "$size"
            else
                log "  Partitions cached: $part_label"
            fi

            local n_parts
            n_parts=$(ls "$part_dir"/*_MSA.fasta 2>/dev/null | wc -l | tr -d ' ')

            # gamma/nogamma variants (norepartition: -f 999999)
            for gamma in "${PASTA_GAMMA[@]}"; do
                local suffix; [[ "$gamma" == "true" ]] && suffix="gamma" || suffix="nogamma"
                local build_dir="${config_dir}/${suffix}"
                log "    $variant/pasta/${part_label}/${suffix} ($n_parts partitions)"
                run_tronko_on_partitions "$part_dir" "$build_dir" "$n_parts" -f 999999
            done

            # SP repartitioning variants
            for sp in "${PASTA_SP_THRESHOLDS[@]}"; do
                local build_dir="${config_dir}/sp${sp}"
                log "    $variant/pasta/${part_label}/sp${sp} ($n_parts partitions)"
                run_tronko_on_partitions "$part_dir" "$build_dir" "$n_parts" -s -u "$sp"
            done
        done

        for diam in "${PASTA_DIAMS[@]}"; do
            local part_label="maxdiam${diam}"
            local config_dir="${pasta_dir}/${part_label}"
            local part_dir="${config_dir}/partitions"

            # Build partitions if not cached
            if [[ ! -d "$part_dir" ]] || [[ -z "$(ls "$part_dir"/*_MSA.fasta 2>/dev/null)" ]]; then
                log "  Partitioning $part_label..."
                python3 "$TRONKO_DIR/tronko-build/partition_and_build.py" \
                    --tree "$ROOTED_TREE" \
                    --fasta "$INPUT_FASTA" \
                    --taxonomy "$INPUT_TAX" \
                    --outdir "$part_dir" \
                    --tronko-outdir "/dev/null" \
                    --strategy centroid \
                    --threads "$THREADS" \
                    --famsa-threads "$FAMSA_THREADS" \
                    --tronko-build "$(which tronko-build)" \
                    --skip-tronko \
                    --max-diam "$diam"
            else
                log "  Partitions cached: $part_label"
            fi

            local n_parts
            n_parts=$(ls "$part_dir"/*_MSA.fasta 2>/dev/null | wc -l | tr -d ' ')

            for gamma in "${PASTA_GAMMA[@]}"; do
                local suffix; [[ "$gamma" == "true" ]] && suffix="gamma" || suffix="nogamma"
                local build_dir="${config_dir}/${suffix}"
                log "    $variant/pasta/${part_label}/${suffix} ($n_parts partitions)"
                run_tronko_on_partitions "$part_dir" "$build_dir" "$n_parts" -f 999999
            done

            for sp in "${PASTA_SP_THRESHOLDS[@]}"; do
                local build_dir="${config_dir}/sp${sp}"
                log "    $variant/pasta/${part_label}/sp${sp} ($n_parts partitions)"
                run_tronko_on_partitions "$part_dir" "$build_dir" "$n_parts" -s -u "$sp"
            done
        done
    done
}

# ============================================================
# AncestralClust pipeline
# ============================================================
run_ac_pipeline() {
    if [[ "$SKIP_AC" == "1" ]]; then
        log "Skipping AC pipeline (SKIP_AC=1)"
        return
    fi

    log "========== ANCESTRALCLUST PIPELINE =========="

    for variant in $VARIANTS; do
        get_input_files "$variant"
        local ac_base="${DB_BASE}/${variant}/ac"
        mkdir -p "$ac_base"

        local seq_count
        seq_count=$(grep -c '^>' "$INPUT_FASTA")
        log "AC $variant: $seq_count sequences"

        for ac_entry in "${AC_CONFIGS[@]}"; do
            local ac_label="${ac_entry%% *}"
            local ac_rest="${ac_entry#* }"
            local ac_bin_size="${ac_rest%% *}"
            local ac_descendants="${ac_rest#* }"
            # Strip "ac_" prefix for dir name (ac_default -> default)
            local config_name="${ac_label#ac_}"
            local config_dir="${ac_base}/${config_name}"

            log "  AC config: $config_name (bin=$ac_bin_size, desc=$ac_descendants)"

            for sp in "${AC_SP_THRESHOLDS[@]}"; do
                local build_dir="${config_dir}/sp${sp}"

                if [[ -f "$build_dir/reference_tree.txt" ]] || \
                   [[ -f "$build_dir/reference_tree.txt.gz" ]] || \
                   [[ -f "$build_dir/reference_tree.trkb" ]]; then
                    log "    $variant/ac/${config_name}/sp${sp} already exists — skipping"
                    continue
                fi

                log "    Building $variant/ac/${config_name}/sp${sp}..."
                bash "$TRONKO_DIR/build-tronko-db.sh" \
                    -f "$INPUT_FASTA" \
                    -t "$INPUT_TAX" \
                    -o "$build_dir" \
                    -p "$MARKER" \
                    -T "$THREADS" \
                    -s "$sp" \
                    -E \
                    -L \
                    -B "$ac_bin_size" \
                    -P "$ac_descendants" \
                    -J "$PARALLEL_JOBS" \
                    --cache-dir "${config_dir}/cache"

                cp "$INPUT_FASTA" "$build_dir/input.fasta"
                cp "$INPUT_TAX" "$build_dir/input_taxonomy.txt"
                log "    Done: $build_dir"
            done
        done
    done
}

# ============================================================
# Main
# ============================================================
log "============================================================"
log "Rebuilding all databases for $MARKER"
log "  Variants: $VARIANTS"
log "  Threads: $THREADS, FAMSA: $FAMSA_THREADS"
log "  PASTA sizes: ${PASTA_SIZES[*]}, diams: ${PASTA_DIAMS[*]}"
log "  PASTA variants per partition: gamma(${#PASTA_GAMMA[@]}) + SP(${#PASTA_SP_THRESHOLDS[@]})"
log "  AC configs: ${#AC_CONFIGS[@]} x ${#AC_SP_THRESHOLDS[@]} SP thresholds"
log "  Skip PASTA: $SKIP_PASTA, Skip AC: $SKIP_AC"
log "============================================================"

# Run per variant: all PASTA + all AC for each variant before moving to the next
ALL_VARIANTS="$VARIANTS"
for CURRENT_VARIANT in $ALL_VARIANTS; do
    log ">>>>>>>>>> Variant: $CURRENT_VARIANT <<<<<<<<<<<"
    VARIANTS="$CURRENT_VARIANT"
    run_pasta_pipeline
    run_ac_pipeline
done
VARIANTS="$ALL_VARIANTS"

log "============================================================"
log "All builds complete for $MARKER"
log "============================================================"
