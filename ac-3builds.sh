#!/bin/bash
set -euo pipefail

# ============================================================
# Build AncestralClust-based tronko databases for a given marker
# Usage: bash ac-3builds.sh <MARKER>
#   e.g. bash ac-3builds.sh ITS2_Plants
#
# Run from ~/tronko-fork
# ============================================================

MARKER="${1:?Usage: bash ac-3builds.sh <MARKER>}"

THREADS=64
FAMSA_THREADS=8
PARALLEL_JOBS=4
AC_BIN_SIZE=10000
AC_DESCENDANTS=75

TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko-fork}"
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"

# ── Input files ──────────────────────────────────────────────
INPUT_FASTA="$HOME/rcrux-py/databases/${MARKER}/filtered/${MARKER}_species.fasta"
INPUT_TAX="$HOME/rcrux-py/databases/${MARKER}/filtered/${MARKER}_species_taxonomy.txt"

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: $INPUT_FASTA not found" >&2; exit 1
fi
if [[ ! -f "$INPUT_TAX" ]]; then
    echo "ERROR: $INPUT_TAX not found" >&2; exit 1
fi

DB_BASE="databases/${MARKER}"

echo "============================================================"
echo "AncestralClust builds for $MARKER"
echo "  Input FASTA:    $INPUT_FASTA"
echo "  Input taxonomy: $INPUT_TAX"
echo "  Threads: $THREADS, FAMSA threads: $FAMSA_THREADS"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo "  AC bin size: $AC_BIN_SIZE, descendants: $AC_DESCENDANTS"
echo "============================================================"
echo ""

for SP in 0.05 0.10 0.20; do
    LABEL="ac_sp${SP}"
    OUTDIR="${DB_BASE}/${LABEL}"

    echo "############################################################"
    echo "# Building: $LABEL (SP threshold = $SP)"
    echo "############################################################"

    if [[ -f "$OUTDIR/reference_tree.txt" ]]; then
        echo "  Already exists: $OUTDIR/reference_tree.txt — skipping"
        echo ""
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

    # Copy input files into the database directory
    cp "$INPUT_FASTA" "$OUTDIR/input.fasta"
    cp "$INPUT_TAX" "$OUTDIR/input_taxonomy.txt"

    echo ""
    echo "  Done: $OUTDIR/reference_tree.txt"
    echo ""
done

echo ""
echo "=== All AncestralClust builds complete ==="
echo ""
echo "$MARKER AC databases:"
echo "  1) ${DB_BASE}/ac_sp0.05/reference_tree.txt"
echo "  2) ${DB_BASE}/ac_sp0.10/reference_tree.txt"
echo "  3) ${DB_BASE}/ac_sp0.20/reference_tree.txt"
