#!/bin/bash
set -euo pipefail

# ============================================================
# Build AncestralClust-based tronko databases for a given marker
# Usage: bash ac-3builds.sh <MARKER> [THREADS]
#   e.g. bash ac-3builds.sh ITS2_Plants 64
#
# Builds both species and LCA taxonomy variants.
# Output layout:
#   databases/MARKER/species/{ac_sp0.05,ac_sp0.10,ac_sp0.20}
#   databases/MARKER/lca/{ac_sp0.05,ac_sp0.10,ac_sp0.20}
#
# Run from ~/tronko-fork
# ============================================================

MARKER="${1:?Usage: bash ac-3builds.sh <MARKER> [THREADS]}"

THREADS="${2:-64}"
FAMSA_THREADS=$(( THREADS / 8 ))
PARALLEL_JOBS=$(( THREADS / 8 ))
AC_BIN_SIZE=10000
AC_DESCENDANTS=75

TRONKO_DIR="${TRONKO_DIR:-$HOME/tronko}"
export PATH="$TRONKO_DIR/bin:$TRONKO_DIR/tronko-build:$PATH"

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

# ── Ensure directory structure ─────────────────────────────────
mkdir -p "${DB_BASE}/species/ac/default" "${DB_BASE}/lca/ac/default"

echo "============================================================"
echo "AncestralClust builds for $MARKER"
echo "  Species FASTA:    $SPECIES_FASTA"
echo "  Species taxonomy: $SPECIES_TAX"
echo "  LCA FASTA:        $LCA_FASTA"
echo "  LCA taxonomy:     $LCA_TAX"
echo "  Threads: $THREADS, FAMSA threads: $FAMSA_THREADS"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo "  AC bin size: $AC_BIN_SIZE, descendants: $AC_DESCENDANTS"
echo "============================================================"
echo ""

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

    AC_DIR="${DB_BASE}/${VARIANT}/ac/default"

    for SP in 0.05 0.10 0.20; do
        OUTDIR="${AC_DIR}/sp${SP}"

        echo "############################################################"
        echo "# Building: ${VARIANT}/ac/default/sp${SP}"
        echo "############################################################"

        if [[ -f "$OUTDIR/reference_tree.txt" ]] || [[ -f "$OUTDIR/reference_tree.txt.gz" ]] || [[ -f "$OUTDIR/reference_tree.trkb" ]]; then
            echo "  Already exists: $OUTDIR — skipping"
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
        echo "  Done: $OUTDIR"
        echo ""
    done
done

echo ""
echo "=== All AncestralClust builds complete ==="
echo ""
echo "$MARKER AC databases:"
for VARIANT in lca species; do
    echo "  ${VARIANT}:"
    echo "    1) ${DB_BASE}/${VARIANT}/ac/default/sp0.05"
    echo "    2) ${DB_BASE}/${VARIANT}/ac/default/sp0.10"
    echo "    3) ${DB_BASE}/${VARIANT}/ac/default/sp0.20"
done
