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
THREADS=32
FAMSA_THREADS=8

# Paths to tronko-fork tools (uses bundled PASTA copy)
TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko-fork}"
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
# Marker 1: 18S_Euk (187K seqs)
# ============================================================
echo "############################################################"
echo "# 18S_Euk"
echo "############################################################"

INPUT_18S="$HOME/rcrux-py/databases/18S_Euk/dc-megablast/18S_Euk_species.fasta"
TAX_18S="$HOME/rcrux-py/databases/18S_Euk/dc-megablast/18S_Euk_species_taxonomy.txt"

run_pasta "$INPUT_18S" "18S_pasta" "pasta_output_18S"
ROOTED_18S="$_ROOTED_TREE"

run_build "18S_maxdiam25"   "databases/18S_pasta_maxdiam25"   "$INPUT_18S" "$TAX_18S" "$ROOTED_18S" --max-diam 25
run_build "18S_maxsize1000" "databases/18S_pasta_maxsize1000" "$INPUT_18S" "$TAX_18S" "$ROOTED_18S" --max-size 1000
run_build "18S_maxsize500"  "databases/18S_pasta_maxsize500"  "$INPUT_18S" "$TAX_18S" "$ROOTED_18S" --max-size 500

# ============================================================
# Marker 2: vert12S / 12SV5 (101K seqs)
# ============================================================
echo ""
echo "############################################################"
echo "# vert12S (12SV5)"
echo "############################################################"

INPUT_12S="$HOME/rcrux-py/databases/12SV5/dc-megablast/12SV5_species.fasta"
TAX_12S="$HOME/rcrux-py/databases/12SV5/dc-megablast/12SV5_species_taxonomy.txt"

run_pasta "$INPUT_12S" "vert12S_pasta" "pasta_output_vert12S"
ROOTED_12S="$_ROOTED_TREE"

run_build "12S_maxdiam25"   "databases/vert12S_pasta_maxdiam25"   "$INPUT_12S" "$TAX_12S" "$ROOTED_12S" --max-diam 25
run_build "12S_maxsize1000" "databases/vert12S_pasta_maxsize1000" "$INPUT_12S" "$TAX_12S" "$ROOTED_12S" --max-size 1000
run_build "12S_maxsize500"  "databases/vert12S_pasta_maxsize500"  "$INPUT_12S" "$TAX_12S" "$ROOTED_12S" --max-size 500

# ============================================================
# Marker 3: vert12S unfiltered / 12SV5 (177K seqs)
# ============================================================
echo ""
echo "############################################################"
echo "# vert12S unfiltered (12SV5 dc-megablast-unfiltered)"
echo "############################################################"

INPUT_12S_UF="$HOME/rcrux-py/databases/12SV5/dc-megablast-unfiltered/12SV5_species.fasta"
TAX_12S_UF="$HOME/rcrux-py/databases/12SV5/dc-megablast-unfiltered/12SV5_species_taxonomy.txt"

run_pasta "$INPUT_12S_UF" "vert12S_unfiltered_pasta" "pasta_output_vert12S_unfiltered"
ROOTED_12S_UF="$_ROOTED_TREE"

run_build "12S_uf_maxdiam25"   "databases/vert12S_unfiltered_pasta_maxdiam25"   "$INPUT_12S_UF" "$TAX_12S_UF" "$ROOTED_12S_UF" --max-diam 25
run_build "12S_uf_maxsize1000" "databases/vert12S_unfiltered_pasta_maxsize1000" "$INPUT_12S_UF" "$TAX_12S_UF" "$ROOTED_12S_UF" --max-size 1000
run_build "12S_uf_maxsize500"  "databases/vert12S_unfiltered_pasta_maxsize500"  "$INPUT_12S_UF" "$TAX_12S_UF" "$ROOTED_12S_UF" --max-size 500

echo ""
echo "=== All builds complete ==="
echo ""
echo "18S databases:"
echo "  1) databases/18S_pasta_maxdiam25/reference_tree.txt"
echo "  2) databases/18S_pasta_maxsize1000/reference_tree.txt"
echo "  3) databases/18S_pasta_maxsize500/reference_tree.txt"
echo ""
echo "vert12S (filtered) databases:"
echo "  4) databases/vert12S_pasta_maxdiam25/reference_tree.txt"
echo "  5) databases/vert12S_pasta_maxsize1000/reference_tree.txt"
echo "  6) databases/vert12S_pasta_maxsize500/reference_tree.txt"
echo ""
echo "vert12S (unfiltered) databases:"
echo "  7) databases/vert12S_unfiltered_pasta_maxdiam25/reference_tree.txt"
echo "  8) databases/vert12S_unfiltered_pasta_maxsize1000/reference_tree.txt"
echo "  9) databases/vert12S_unfiltered_pasta_maxsize500/reference_tree.txt"
