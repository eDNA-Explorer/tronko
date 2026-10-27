#!/bin/bash
set -euo pipefail

# ============================================================
# Build all tronko databases (AC + PASTA) for a given marker
#
# Usage:
#   bash build-marker.sh <marker_name> <fasta> <taxonomy>
#
# Example:
#   bash build-marker.sh 12S_MiFish_U \
#     ~/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species.fasta \
#     ~/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species_taxonomy.txt
# ============================================================

MARKER="${1:?Usage: build-marker.sh <marker> <fasta> <taxonomy>}"
INPUT_FASTA="${2:?Missing input FASTA}"
INPUT_TAX="${3:?Missing input taxonomy}"

THREADS=128
FAMSA_THREADS=16
PARALLEL_JOBS=8
AC_BIN_SIZE=10000
AC_DESCENDANTS=75

TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko}"
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"
export PASTA_DIR="$TRONKO_DIR/pasta"
export PASTA_TREE_BACKEND="${TREE_BACKEND:-veryfasttree}"

DB_BASE="databases/${MARKER}"
mkdir -p "$DB_BASE"

echo "============================================================"
echo "Building all databases for: $MARKER"
echo "  Input FASTA:    $INPUT_FASTA ($(grep -c '^>' "$INPUT_FASTA") sequences)"
echo "  Input taxonomy: $INPUT_TAX"
echo "  Output base:    $DB_BASE"
echo "============================================================"
echo ""

# ── AncestralClust builds ────────────────────────────────────
echo "############################################################"
echo "# AncestralClust builds"
echo "############################################################"

for SP in 0.05 0.10 0.20; do
    LABEL="ac_sp${SP}"
    OUTDIR="${DB_BASE}/${LABEL}"

    echo ""
    echo "=== AC: $LABEL (SP threshold = $SP) ==="

    if [[ -f "$OUTDIR/reference_tree.txt" ]] || [[ -f "$OUTDIR/reference_tree.txt.gz" ]]; then
        echo "  Already exists — skipping"
        continue
    fi

    bash build-tronko-db.sh \
        -f "$INPUT_FASTA" \
        -t "$INPUT_TAX" \
        -o "$OUTDIR" \
        -p "$MARKER" \
        -T "$THREADS" \
        -s "$SP" \
        -F \
        -E \
        -L \
        -B "$AC_BIN_SIZE" \
        -P "$AC_DESCENDANTS" \
        -J "$PARALLEL_JOBS"

    echo "  Done: $OUTDIR/reference_tree.txt"
done

# ── PASTA builds ─────────────────────────────────────────────
echo ""
echo "############################################################"
echo "# PASTA builds"
echo "############################################################"

PASTA_OUT="${DB_BASE}/pasta_output"

# Run PASTA if tree doesn't exist yet
if [[ -f "${PASTA_OUT}/${MARKER}_pasta.rooted.tre" ]]; then
    echo "  PASTA tree already exists — skipping PASTA alignment"
else
    echo ""
    echo "=== Running PASTA ==="
    mkdir -p "$PASTA_OUT"

    python3 "$PASTA_DIR/run_pasta.py" \
        --input="$INPUT_FASTA" \
        --datatype=dna \
        --num-cpus="$THREADS" \
        --iter-limit=3 \
        -o "$PASTA_OUT" \
        -j "${MARKER}_pasta"

    # Find output tree
    PASTA_TREE=$(ls ${PASTA_OUT}/${MARKER}_pasta*.tre 2>/dev/null | head -1)
    NAME_TRANS="${PASTA_OUT}/${MARKER}_pasta_temp_name_translation.txt"

    # Restore original names, quote colons for nw_reroot
    python3 -c "
import re
trans = {}
with open('$NAME_TRANS') as f:
    lines = f.read().strip().split('\n')
    i = 0
    while i < len(lines):
        safe = lines[i].strip()
        if i+1 < len(lines):
            orig = lines[i+1].strip()
            if safe and orig:
                trans[safe] = orig
        i += 3

with open('$PASTA_TREE') as f:
    data = f.read()
data = data.replace(\"'\", '')
data = re.sub(r'\)([0-9][0-9.eE+-]*):', '):', data)

def replace_name(m):
    safe = m.group(1)
    orig = trans.get(safe, safe)
    if ':' in orig or '|' in orig:
        return \"'\" + orig + \"'\"
    return orig

data = re.sub(r'([A-Za-z0-9_]+)(?=:)', replace_name, data)
data = data.rstrip().rstrip(';') + ';'
with open('$PASTA_TREE', 'w') as f:
    f.write(data + '\n')
print(f'Restored {len(trans)} names')
"

    # Midpoint root
    ROOTED="${PASTA_OUT}/${MARKER}_pasta.rooted.tre"
    nw_reroot "$PASTA_TREE" > "$ROOTED"
    sed -i "s/'//g" "$ROOTED"
    echo "  Rooted tree: $ROOTED"
fi

ROOTED="${PASTA_OUT}/${MARKER}_pasta.rooted.tre"

# Partition + tronko-build for each config
for CONFIG in "maxdiam25 --max-diam 25" "maxsize1000 --max-size 1000" "maxsize500 --max-size 500"; do
    LABEL=$(echo $CONFIG | awk '{print $1}')
    ARGS=$(echo $CONFIG | cut -d' ' -f2-)
    OUTDIR="${DB_BASE}/${LABEL}"

    echo ""
    echo "=== PASTA partition: $LABEL ==="

    if [[ -f "$OUTDIR/reference_tree.txt" ]]; then
        echo "  Already exists — skipping"
        continue
    fi

    python3 "$TRONKO_DIR/tronko-build/partition_and_build.py" \
        --tree "$ROOTED" \
        --fasta "$INPUT_FASTA" \
        --taxonomy "$INPUT_TAX" \
        --outdir "${OUTDIR}/partitions" \
        --tronko-outdir "$OUTDIR" \
        --strategy centroid \
        --threads "$THREADS" \
        --famsa-threads "$FAMSA_THREADS" \
        --tronko-build "$(which tronko-build)" \
        $ARGS

    echo "  Done: $OUTDIR/reference_tree.txt"
done

echo ""
echo "============================================================"
echo "All builds complete for $MARKER"
echo "============================================================"
ls -lh ${DB_BASE}/*/reference_tree.txt 2>/dev/null
