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

# Replace safe names with originals, quoting names that contain colons
def replace_name(m):
    safe = m.group(1)
    orig = trans.get(safe, safe)
    if ':' in orig or '|' in orig:
        return \"'\" + orig + \"'\"
    return orig

data = re.sub(r'([A-Za-z0-9_]+)(?=:)', replace_name, data)
data = data.rstrip().rstrip(';') + ';'

with open('$pasta_tree', 'w') as f:
    f.write(data + '\n')

print(f'Restored {len(trans)} names, {sum(1 for v in trans.values() if \":\" in v)} quoted')
"

    # Midpoint root
    local rooted_tree="${output_dir}/${job_name}.rooted.tre"
    nw_reroot "$pasta_tree" > "$rooted_tree"

    # Strip quotes after rerooting — names now safe in Newick context
    sed -i "s/'//g" "$rooted_tree"
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
# Marker: 16Smamm (~99K seqs)
# ============================================================
echo ""
echo "############################################################"
echo "# 16Smamm"
echo "############################################################"

INPUT_16S="$HOME/rcrux-py/databases/16Smamm/filtered/16Smamm_species.fasta"
TAX_16S="$HOME/rcrux-py/databases/16Smamm/filtered/16Smamm_species_taxonomy.txt"

run_pasta "$INPUT_16S" "16Smamm_pasta" "pasta_output_16Smamm"
ROOTED_16S="$_ROOTED_TREE"

run_build "16S_maxdiam25"   "databases/16Smamm_pasta_maxdiam25"   "$INPUT_16S" "$TAX_16S" "$ROOTED_16S" --max-diam 25
run_build "16S_maxsize1000" "databases/16Smamm_pasta_maxsize1000" "$INPUT_16S" "$TAX_16S" "$ROOTED_16S" --max-size 1000
run_build "16S_maxsize500"  "databases/16Smamm_pasta_maxsize500"  "$INPUT_16S" "$TAX_16S" "$ROOTED_16S" --max-size 500

echo ""
echo "=== All builds complete ==="
echo ""
echo "16Smamm databases:"
echo "  1) databases/16Smamm_pasta_maxdiam25/reference_tree.txt"
echo "  2) databases/16Smamm_pasta_maxsize1000/reference_tree.txt"
echo "  3) databases/16Smamm_pasta_maxsize500/reference_tree.txt"
