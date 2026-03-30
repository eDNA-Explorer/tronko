#!/bin/bash
set -euo pipefail

# ============================================================
# Build AncestralClust-based tronko databases for a given marker
# Usage: bash ac-6builds.sh <MARKER> [THREADS] [CONFIG] [SP_MODE]
#   e.g. bash ac-6builds.sh 16Smamm 64
#         bash ac-6builds.sh vert12S 64 more_bins
#         bash ac-6builds.sh vert12S 40 default corrected
#         bash ac-6builds.sh vert12S 40 all corrected
#
# CONFIG:  default | more_bins | fewer_bins | all (default: all)
# SP_MODE: legacy | corrected (default: legacy)
#   legacy:    SP thresholds 0.05, 0.10, 0.20 with -L flag
#   corrected: SP thresholds 2.00, 2.50, 2.75 without -L flag
#
# Output layout:
#   databases/MARKER/{lca,species}/ac/{config}/sp{threshold}/
#
# Run from ~/tronko
# ============================================================

MARKER="${1:?Usage: bash ac-6builds.sh <MARKER> [THREADS] [CONFIG] [SP_MODE]}"

THREADS="${2:-64}"
ONLY_CONFIG="${3:-all}"
SP_MODE="${4:-legacy}"

FAMSA_THREADS=$(( THREADS / 8 ))
PARALLEL_JOBS=$(( THREADS / 8 ))
AC_DESCENDANTS=75

AC_CONFIGS=(
    "default 10000"
    "more_bins 5000"
    "fewer_bins 20000"
)

if [[ "$SP_MODE" == "corrected" ]]; then
    SP_THRESHOLDS=(2.00 2.50 2.75)
    LEGACY_FLAG=""
else
    SP_THRESHOLDS=(0.05 0.10 0.20)
    LEGACY_FLAG="-L"
fi

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

echo "============================================================"
echo "AncestralClust builds for $MARKER"
echo "  Species FASTA:    $SPECIES_FASTA"
echo "  Species taxonomy: $SPECIES_TAX"
echo "  LCA FASTA:        $LCA_FASTA"
echo "  LCA taxonomy:     $LCA_TAX"
echo "  Threads: $THREADS, FAMSA threads: $FAMSA_THREADS"
echo "  Parallel jobs: $PARALLEL_JOBS"
echo "  Config: $ONLY_CONFIG"
echo "  SP mode: $SP_MODE (thresholds: ${SP_THRESHOLDS[*]})"
echo "  Legacy flag: ${LEGACY_FLAG:-(none)}"
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

    for AC_ENTRY in "${AC_CONFIGS[@]}"; do
        CONFIG_NAME="${AC_ENTRY%% *}"
        AC_BIN_SIZE="${AC_ENTRY#* }"

        # Skip if a specific config was requested and this isn't it
        if [[ "$ONLY_CONFIG" != "all" && "$CONFIG_NAME" != "$ONLY_CONFIG" ]]; then
            continue
        fi
        CONFIG_DIR="${DB_BASE}/${VARIANT}/ac/${CONFIG_NAME}"
        mkdir -p "$CONFIG_DIR"

        echo ""
        echo "  AC config: $CONFIG_NAME (bin_size=$AC_BIN_SIZE)"

        for SP in "${SP_THRESHOLDS[@]}"; do
            OUTDIR="${CONFIG_DIR}/sp${SP}"

            echo "############################################################"
            echo "# Building: ${VARIANT}/ac/${CONFIG_NAME}/sp${SP} ($SP_MODE)"
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
                -E \
                $LEGACY_FLAG \
                -B "$AC_BIN_SIZE" \
                -P "$AC_DESCENDANTS" \
                -J "$PARALLEL_JOBS" \
                --cache-dir "$CONFIG_DIR/cache"

            cp "$INPUT_FASTA" "$OUTDIR/input.fasta"
            cp "$INPUT_TAX" "$OUTDIR/input_taxonomy.txt"

            echo ""
            echo "  Done: $OUTDIR"
            echo ""
        done
    done
done

echo ""
echo "=== All AncestralClust builds complete ==="
echo ""
echo "$MARKER AC databases ($SP_MODE):"
for VARIANT in lca species; do
    echo "  ${VARIANT}:"
    for AC_ENTRY in "${AC_CONFIGS[@]}"; do
        CONFIG_NAME="${AC_ENTRY%% *}"
        if [[ "$ONLY_CONFIG" != "all" && "$CONFIG_NAME" != "$ONLY_CONFIG" ]]; then
            continue
        fi
        for SP in "${SP_THRESHOLDS[@]}"; do
            echo "    ${DB_BASE}/${VARIANT}/ac/${CONFIG_NAME}/sp${SP}"
        done
    done
done
