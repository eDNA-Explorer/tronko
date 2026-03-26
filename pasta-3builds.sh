#!/bin/bash
set -euo pipefail

# ============================================================
# Build PASTA-based tronko databases for multiple markers
# Run from ~/tronko-fork
#
# Tree backend: set TREE_BACKEND=veryfasttree to use VeryFastTree
#               (default: fasttree)
# ============================================================

TREE_BACKEND="${TREE_BACKEND:-veryfasttree}"
THREADS=128
FAMSA_THREADS=16

# Paths to tronko-fork tools (uses bundled PASTA copy)
TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko}"
PASTA_DIR="$TRONKO_DIR/pasta"
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
if [[ -d /usr/lib/jvm/java-17-openjdk-amd64 ]]; then
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
    mkdir -p "$output_dir"

    python3 "$PASTA_DIR/run_pasta.py" \
        --input="$input_fasta" \
        --datatype=dna \
        --num-cpus="$THREADS" \
        --iter-limit=3 \
        -o "$output_dir" \
        -j "$job_name"

    # Find the output tree
    local pasta_tree
    pasta_tree=$(ls ${output_dir}/${job_name}*.tre | head -1)
    echo "PASTA tree: $pasta_tree"

    # Fix: strip single quotes from leaf names
    sed -i "s/'//g" "$pasta_tree"

    # Fix: strip FastTree support values: )0.996:0.068 -> ):0.068
    python3 -c "
import re
with open('$pasta_tree') as f:
    data = f.read()
with open('$pasta_tree', 'w') as f:
    f.write(re.sub(r'\)([0-9][0-9.eE+-]*):', '):', data))
"

    # Midpoint root
    local rooted_tree="${output_dir}/${job_name}.rooted.tre"
    nw_reroot "$pasta_tree" > "$rooted_tree"
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
    local partition_dir="pasta_partitions_${label}"
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

    echo "  Database: $outdir/reference_tree.txt"
    echo "  Reference FASTA: $REF_FASTA"
}

# ============================================================
# Marker: 12S_MiFish_U filtered (~59K seqs)
# ============================================================
echo "############################################################"
echo "# 12S_MiFish_U"
echo "############################################################"

INPUT_MIFISH="$HOME/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species.fasta"
TAX_MIFISH="$HOME/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species_taxonomy.txt"

run_pasta "$INPUT_MIFISH" "12S_MiFish_pasta" "pasta_output_12S_MiFish"
ROOTED_MIFISH="$_ROOTED_TREE"

run_build "MiFish_maxdiam25"   "databases/MiFish_pasta_maxdiam25"   "$INPUT_MIFISH" "$TAX_MIFISH" "$ROOTED_MIFISH" --max-diam 25
run_build "MiFish_maxsize1000" "databases/MiFish_pasta_maxsize1000" "$INPUT_MIFISH" "$TAX_MIFISH" "$ROOTED_MIFISH" --max-size 1000
run_build "MiFish_maxsize500"  "databases/MiFish_pasta_maxsize500"  "$INPUT_MIFISH" "$TAX_MIFISH" "$ROOTED_MIFISH" --max-size 500

echo ""
echo "=== All builds complete ==="
echo ""
echo "12S_MiFish_U databases:"
echo "  1) databases/MiFish_pasta_maxdiam25/reference_tree.txt"
echo "  2) databases/MiFish_pasta_maxsize1000/reference_tree.txt"
echo "  3) databases/MiFish_pasta_maxsize500/reference_tree.txt"
