#!/usr/bin/env bash
# Build the CO1 (Leray) tronko database:
#   AncestralClust clustering + LCA taxonomy + legacy SP 0.10
#
# Safe to re-run verbatim after a crash — --cache-dir makes it resume.

set -euo pipefail

MARKER=CO1_mlCOIintF_Fol-degen-rev
SRC="$HOME/rcrux-py/databases/$MARKER/dedup"
BASE="$HOME/tronko-build-branch/databases/$MARKER/lca/ac/default"
OUT="$BASE/sp0.10"

# 2-column taxonomy: rcrux emits accession<TAB>taxid<TAB>lineage, but
# tronko-build treats EVERYTHING after the first tab as the lineage, which
# would silently shift every rank. Column stripped in advance.
TAX="$BASE/lca_taxonomy.2col.txt"

export PATH="$HOME/tronko-build-branch/bin:$HOME/tronko-build-branch/tronko-build:$HOME/tronko-build-branch/tronko-convert:$PATH"

mkdir -p "$OUT"
cd "$HOME/tronko-build-branch"

# Fail fast rather than burning hours on a bad input
[[ -f "$SRC/lca.fasta" ]] || { echo "missing $SRC/lca.fasta" >&2; exit 1; }
[[ -f "$TAX" ]]           || { echo "missing $TAX" >&2; exit 1; }
if [[ "$(head -1 "$TAX" | awk -F'\t' '{print NF}')" != "2" ]]; then
    echo "ERROR: $TAX is not 2-column — would corrupt taxonomy" >&2
    exit 1
fi

echo "START $(date -Is)" >> "$OUT/status.txt"

set +e
bash build-tronko-db.sh \
    -f "$SRC/lca.fasta" \
    -t "$TAX" \
    -o "$OUT" \
    -p "$MARKER" \
    -T 128 \
    -s 0.10 \
    -E \
    -L \
    -B 10000 \
    -P 75 \
    -J 16 \
    --cache-dir "$BASE/cache"
RC=$?
set -e

echo "EXIT $RC $(date -Is)" >> "$OUT/status.txt"
exit $RC
