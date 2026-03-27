#!/bin/bash
set -euo pipefail

# ============================================================
# Build PASTA-based tronko databases for a given marker
# Usage: bash pasta-3builds.sh <MARKER>
#   e.g. bash pasta-3builds.sh ITS2_Plants
#
# Builds both species and LCA taxonomy variants.
# Output layout:
#   databases/MARKER/species/{pasta_output,maxdiam25,maxsize1000,maxsize500}
#   databases/MARKER/lca/{pasta_output,maxdiam25,maxsize1000,maxsize500}
#
# Tree backend: set TREE_BACKEND=veryfasttree to use VeryFastTree
#               (default: veryfasttree)
# ============================================================

MARKER="${1:?Usage: bash pasta-3builds.sh <MARKER>}"

TREE_BACKEND="${TREE_BACKEND:-veryfasttree}"
THREADS=56
FAMSA_THREADS=16

# Paths to tronko-fork tools (uses bundled PASTA copy)
TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko}"
PASTA_DIR="$TRONKO_DIR/pasta"
export PASTA_DIR
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"
export PASTA_TREE_BACKEND="$TREE_BACKEND"

# Verify everything is on PATH
required_cmds=(mafft hmmbuild nw_reroot java famsa tronko-build)
if [[ "$TREE_BACKEND" == "veryfasttree" ]]; then
    required_cmds+=(VeryFastTree)
else
    required_cmds+=(FastTree)
fi
for cmd in "${required_cmds[@]}"; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "WARNING: $cmd not found on PATH"
    fi
done

# Java for OPAL merger
if [[ -d /opt/homebrew/opt/openjdk ]]; then
    export JAVA_HOME="/opt/homebrew/opt/openjdk"
elif [[ -d /usr/lib/jvm/java-17-openjdk-amd64 ]]; then
    export JAVA_HOME="/usr/lib/jvm/java-17-openjdk-amd64"
elif [[ -d /usr/lib/jvm/java-21-openjdk-amd64 ]]; then
    export JAVA_HOME="/usr/lib/jvm/java-21-openjdk-amd64"
fi
export PATH="${JAVA_HOME:-}/bin:$PATH"

echo "Tree backend: $TREE_BACKEND"
echo "PASTA dir:    $PASTA_DIR"
echo ""

# ── Helper: run PASTA for a marker ───────────────────────
run_pasta() {
    local input_fasta="$1"
    local job_name="$2"
    local output_dir="$3"

    echo "=== Running PASTA: $job_name ==="
    rm -rf "$HOME/.pasta/${job_name}"
    mkdir -p "$output_dir"

    python3 "$PASTA_DIR/run_pasta.py" \
        --input="$input_fasta" \
        --datatype=dna \
        --num-cpus="$THREADS" \
        --iter-limit=3 \
        -o "$output_dir" \
        -j "$job_name" || true  # PASTA exits non-zero on harmless "Refused to clean" warning

    # Find the output tree (non-empty, non-temp .tre file)
    local pasta_tree
    pasta_tree=$(find "$output_dir" -name "${job_name}*.tre" ! -name "*_temp_*" -size +1k 2>/dev/null | head -1)
    if [[ -z "$pasta_tree" ]]; then
        echo "ERROR: PASTA did not produce a tree file" >&2; exit 1
    fi
    echo "PASTA tree: $pasta_tree"

    # Post-process PASTA tree:
    # 1. Strip quotes and support values
    # 2. Restore original names from PASTA's safe-name translation
    # 3. Quote leaf names containing colons (coordinate ranges) for nw_reroot
    # 4. Midpoint root, then unquote
    local name_trans="${output_dir}/${job_name}_temp_name_translation.txt"
    python3 -c "
import re, sys

# Load name translation (safe -> original)
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

# Strip quotes and support values
data = data.replace(\"'\", '')
data = re.sub(r'\)([0-9][0-9.eE+-]*):', '):', data)

# Replace safe names with originals, always quoting for nw_reroot safety
def replace_name(m):
    safe = m.group(1)
    orig = trans.get(safe, safe)
    return \"'\" + orig + \"'\"

data = re.sub(r'([A-Za-z0-9_.-]+)(?=:)', replace_name, data)
# Ensure exactly one tree: strip trailing whitespace/semicolons, add one semicolon
data = re.sub(r'[;\s]+$', '', data) + ';'

with open('$pasta_tree', 'w') as f:
    f.write(data + '\n')

print(f'Restored {len(trans)} names, {sum(1 for v in trans.values() if \":\" in v)} quoted')
"

    # Midpoint root
    local rooted_tree="${output_dir}/${job_name}.rooted.tre"
    nw_reroot "$pasta_tree" > "$rooted_tree"

    # Strip quotes after rerooting — names now safe in Newick context
    sed -i '' "s/'//g" "$rooted_tree"
    echo "Rooted tree: $rooted_tree"

    # Return path via global
    _ROOTED_TREE="$rooted_tree"
}

# ── Helper: partition + tronko-build ─────────────────────
run_build() {
    local label="$1"
    local outdir="$2"
    local input_fasta="$3"
    local taxonomy="$4"
    local rooted_tree="$5"
    local partition_dir="${outdir}/partitions_${label}"
    shift 5

    echo ""
    echo "========================================"
    echo "=== Build: $label ==="
    echo "========================================"

    python3 "$TRONKO_DIR/tronko-build/partition_and_build.py" \
        --tree "$rooted_tree" \
        --fasta "$input_fasta" \
        --taxonomy "$taxonomy" \
        --outdir "$partition_dir" \
        --tronko-outdir "$outdir" \
        --strategy centroid \
        --threads "$THREADS" \
        --famsa-threads "$FAMSA_THREADS" \
        --tronko-build "$(which tronko-build)" \
        "$@"

    # Collect reference FASTA
    REF_FASTA="$outdir/marker.fasta"
    if [[ ! -f "$REF_FASTA" ]]; then
        cat ${partition_dir}/partition*_unaligned.fasta > "$REF_FASTA" 2>/dev/null || \
        cp "$input_fasta" "$REF_FASTA"
    fi

    # Copy input files into the database directory
    cp "$input_fasta" "$outdir/input.fasta"
    cp "$taxonomy" "$outdir/input_taxonomy.txt"

    echo "  Database: $outdir/reference_tree.txt"
    echo "  Reference FASTA: $REF_FASTA"
}

# ============================================================
# Marker: $MARKER
# ============================================================
echo "############################################################"
echo "# $MARKER"
echo "############################################################"

# ── Input files ──────────────────────────────────────────────
SPECIES_FASTA="$HOME/rcrux-py/databases/${MARKER}/unfiltered/${MARKER}_species.fasta"
SPECIES_TAX="$HOME/rcrux-py/databases/${MARKER}/unfiltered/${MARKER}_species_taxonomy.txt"
LCA_FASTA="$HOME/rcrux-py/databases/${MARKER}/unfiltered/${MARKER}_lca.fasta"
LCA_TAX="$HOME/rcrux-py/databases/${MARKER}/unfiltered/${MARKER}_lca_taxonomy.txt"

for f in "$SPECIES_FASTA" "$SPECIES_TAX" "$LCA_FASTA" "$LCA_TAX"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found" >&2; exit 1
    fi
done

DB_BASE="databases/${MARKER}"

# ── Migrate old flat layout to species/ subdirectory ─────────
SPECIES_DIR="${DB_BASE}/species"
LCA_DIR="${DB_BASE}/lca"
mkdir -p "$SPECIES_DIR" "$LCA_DIR"

for item in pasta_output maxdiam25 maxsize1000 maxsize500; do
    if [[ -d "${DB_BASE}/${item}" ]] && [[ ! -d "${SPECIES_DIR}/${item}" ]]; then
        echo "Migrating ${DB_BASE}/${item} -> ${SPECIES_DIR}/${item}"
        mv "${DB_BASE}/${item}" "${SPECIES_DIR}/${item}"
    fi
    if [[ -f "${DB_BASE}/${item}.dvc" ]] && [[ ! -f "${SPECIES_DIR}/${item}.dvc" ]]; then
        echo "Migrating ${DB_BASE}/${item}.dvc -> ${SPECIES_DIR}/${item}.dvc"
        mv "${DB_BASE}/${item}.dvc" "${SPECIES_DIR}/${item}.dvc"
    fi
done

# ── Build both variants ─────────────────────────────────────
for VARIANT in lca species; do
    echo ""
    echo "============================================================"
    echo "  Variant: $VARIANT"
    echo "============================================================"

    if [[ "$VARIANT" == "species" ]]; then
        INPUT_FASTA="$SPECIES_FASTA"
        INPUT_TAX="$SPECIES_TAX"
    else
        INPUT_FASTA="$LCA_FASTA"
        INPUT_TAX="$LCA_TAX"
    fi

    SEQ_COUNT=$(grep -c '^>' "$INPUT_FASTA")
    echo "  FASTA:    $INPUT_FASTA ($SEQ_COUNT sequences)"
    echo "  Taxonomy: $INPUT_TAX"
    echo ""

    VARIANT_DIR="${DB_BASE}/${VARIANT}"
    PASTA_OUT="${VARIANT_DIR}/pasta_output"
    JOB_NAME="${MARKER}_${VARIANT}_pasta"

    # Skip PASTA if rooted tree already exists (and is non-empty)
    PASTA_TREE_GLOB="${PASTA_OUT}/${JOB_NAME}.rooted.tre"
    if [[ -s "$PASTA_TREE_GLOB" ]]; then
        echo "PASTA tree already exists for $VARIANT — skipping PASTA"
        _ROOTED_TREE="$PASTA_TREE_GLOB"
    else
        run_pasta "$INPUT_FASTA" "$JOB_NAME" "$PASTA_OUT"
    fi
    ROOTED_TREE="$_ROOTED_TREE"

    # Build each partition strategy, skipping if already done
    for build_args in "maxdiam25 --max-diam 25" "maxsize1000 --max-size 1000" "maxsize500 --max-size 500"; do
        BUILD_NAME="${build_args%% *}"
        BUILD_FLAGS="${build_args#* }"
        BUILD_DIR="${VARIANT_DIR}/${BUILD_NAME}"

        if [[ -f "${BUILD_DIR}/reference_tree.txt" ]]; then
            echo "  Already exists: ${BUILD_DIR}/reference_tree.txt — skipping"
            continue
        fi

        run_build "${MARKER}_${VARIANT}_${BUILD_NAME}" "$BUILD_DIR" "$INPUT_FASTA" "$INPUT_TAX" "$ROOTED_TREE" $BUILD_FLAGS
    done
done

echo ""
echo "=== All builds complete ==="
echo ""
echo "$MARKER PASTA databases:"
for VARIANT in lca species; do
    echo "  ${VARIANT}:"
    echo "    1) ${DB_BASE}/${VARIANT}/maxdiam25/reference_tree.txt"
    echo "    2) ${DB_BASE}/${VARIANT}/maxsize1000/reference_tree.txt"
    echo "    3) ${DB_BASE}/${VARIANT}/maxsize500/reference_tree.txt"
done
