#!/bin/bash
set -euo pipefail

# ============================================================
# Build AncestralClust-based tronko databases for 12S_MiFish_U
# at three SP thresholds (0.05, 0.10, 0.20)
#
# Run from ~/tronko
# ============================================================

THREADS=128
FAMSA_THREADS=16
PARALLEL_JOBS=8
AC_BIN_SIZE=10000
AC_DESCENDANTS=75

TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko}"
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"

# ── Input files ──────────────────────────────────────────────
INPUT_FASTA="$HOME/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species.fasta"
INPUT_TAX="$HOME/rcrux-py/databases/12S_MiFish_U/filtered/12S_MiFish_U_species_taxonomy.txt"

DB_BASE="databases/12S_MiFish_U"

echo "============================================================"
echo "AncestralClust builds for 12S_MiFish_U"
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
        -p "12S_MiFish_U" \
        -T "$THREADS" \
        -s "$SP" \
        -F \
        -E \
        -L \
        -B "$AC_BIN_SIZE" \
        -P "$AC_DESCENDANTS" \
        -J "$PARALLEL_JOBS"

    echo ""
    echo "  Done: $OUTDIR/reference_tree.txt"
    echo ""
done

echo ""
echo "=== All AncestralClust builds complete ==="
echo ""
echo "12S_MiFish_U AC databases:"
echo "  1) ${DB_BASE}/ac_sp0.05/reference_tree.txt"
echo "  2) ${DB_BASE}/ac_sp0.10/reference_tree.txt"
echo "  3) ${DB_BASE}/ac_sp0.20/reference_tree.txt"
