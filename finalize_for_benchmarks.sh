#!/usr/bin/env bash
# Post-build finalize: make the CO1 tronko master DB consumable by
# ~/edna-explorer-data-pipelines (projects/assignment_benchmarks).
#
# The harness's ablate_database() hardcodes TronkoBuildConfig(primer_name="marker"),
# so it looks for marker.fasta / marker_taxonomy.txt regardless of what -p was used
# at build time. It also uses partition* files + tree_list.txt (NOT exported_subtrees/),
# and treats a missing input_taxonomy.txt as "database not built" — silently dropping
# tronko from --tools.
#
# Safe to re-run.

set -euo pipefail

MARKER=CO1_mlCOIintF_Fol-degen-rev
BASE="$HOME/tronko-build-branch/databases/$MARKER/lca/ac/default"
OUT="$BASE/sp0.10"
SRC="$HOME/rcrux-py/databases/$MARKER/dedup"

cd "$OUT"

echo "=== Aliasing ${MARKER}.* -> marker.* (harness hardcodes primer_name='marker')"
if [[ -f "${MARKER}.fasta" ]]; then
    ln -sf "${MARKER}.fasta" marker.fasta
    for ext in amb ann bwt pac sa; do
        [[ -f "${MARKER}.fasta.${ext}" ]] && ln -sf "${MARKER}.fasta.${ext}" "marker.fasta.${ext}"
    done
else
    echo "  WARNING: ${MARKER}.fasta not present yet — build may not have reached Step 5"
fi
[[ -f "${MARKER}_taxonomy.txt" ]] && ln -sf "${MARKER}_taxonomy.txt" marker_taxonomy.txt

echo "=== Ensuring build inputs are recorded"
[[ -f input_taxonomy.txt ]] || cp "$BASE/lca_taxonomy.2col.txt" input_taxonomy.txt
[[ -f input.fasta ]]        || cp "$SRC/lca.fasta" input.fasta

echo
echo "=== Required artifact check ==="
rc=0
check() {
    if compgen -G "$1" > /dev/null; then
        printf '  OK      %s\n' "$2"
    else
        printf '  MISSING %s\n' "$2"; rc=1
    fi
}

# Assign-time set
check "reference_tree.trkb"                  "reference_tree.trkb"
check "reference_tree.txt.gz"                "reference_tree.txt.gz (text form; ablation patches this)"
check "${MARKER}.fasta"                      "${MARKER}.fasta"
check "${MARKER}_taxonomy.txt"               "${MARKER}_taxonomy.txt"
for ext in amb ann bwt pac sa; do
    check "${MARKER}.fasta.${ext}"           "BWA index .${ext}"
done

# Ablation set
check "marker.fasta"                         "marker.fasta (alias)"
check "marker_taxonomy.txt"                  "marker_taxonomy.txt (alias)"
check "tree_list.txt"                        "tree_list.txt (absent => slow full-rebuild fallback)"
check "input_taxonomy.txt"                   "input_taxonomy.txt (absent => tronko silently dropped)"
check "input.fasta"                          "input.fasta"
check "partition*.fasta"                     "partition{N}.fasta"
check "partition*_MSA.fasta"                 "partition{N}_MSA.fasta"
check "partition*_taxonomy.txt"              "partition{N}_taxonomy.txt"
check "RAxML_bestTree.partition*.reroot"     "RAxML_bestTree.partition{N}.reroot"

echo
echo "  partitions: $(ls partition*.fasta 2>/dev/null | grep -v '_MSA\|_taxonomy' | wc -l)"
echo "  trees in tree_list.txt: $(wc -l < tree_list.txt 2>/dev/null || true)"
echo
[[ $rc -eq 0 ]] && echo "ALL REQUIRED ARTIFACTS PRESENT" || echo "INCOMPLETE — see MISSING above"
exit $rc
