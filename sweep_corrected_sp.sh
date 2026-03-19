#!/bin/bash
# Sweep corrected SP thresholds using cached AncestralClust clusters.
#
# The original sweep_ac_* databases were built with a buggy SP normalization
# (divided by numspec AND numpairs). This re-runs tronko-build with corrected
# SP thresholds on the corrected [-2, +3] scale.
#
# Reuses cached clusters from: databases/sweep_clusters/
# No need to re-run AncestralClust or FAMSA.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TRONKO_BUILD="$SCRIPT_DIR/tronko-build/tronko-build"
CLUSTER_BASE="$SCRIPT_DIR/databases/sweep_clusters"
DB_BASE="$SCRIPT_DIR/databases"
INPUT_FASTA="$SCRIPT_DIR/tronko-build/example_datasets/vert12S/combined_input.fasta"
INPUT_TAXONOMY="$SCRIPT_DIR/tronko-build/example_datasets/vert12S/12SV5_species_taxonomy.txt"

# Corrected SP thresholds
THRESHOLDS=(2.00 2.50 2.75)

FAMSA_THREADS=4

# Config definitions: config_name cluster_dir n_clusters
CONFIGS="ac_default ac_default 8
ac_more_bins ac_more_bins 21
ac_fewer_bins ac_fewer_bins 2
ac_single ac_single 1"

echo "================================================================"
echo "Corrected SP Threshold Sweep"
echo "================================================================"
echo "Thresholds: ${THRESHOLDS[*]}"
echo "tronko-build: $TRONKO_BUILD"
echo ""

echo "$CONFIGS" | while read -r config_name ac_dir n_clusters; do
    merged_dir="$CLUSTER_BASE/$ac_dir"

    if [[ ! -d "$merged_dir" ]]; then
        echo "SKIP $config_name: $merged_dir not found"
        continue
    fi

    for sp in "${THRESHOLDS[@]}"; do
        # Format threshold for directory name (e.g., 2.00 -> 2.00)
        sp_label=$(printf "%.2f" "$sp")
        db_name="sweep_${config_name}_sp_${sp_label}_corrected"
        db_dir="$DB_BASE/$db_name"

        if [[ -f "$db_dir/reference_tree.txt" ]]; then
            echo "CACHED: $db_name (reference_tree.txt exists)"
            continue
        fi

        echo ""
        echo "======================================================================"
        echo "DATABASE: $db_name"
        echo "  Config: $config_name ($n_clusters clusters)"
        echo "  SP threshold: $sp (corrected scale, range [-2, +3])"
        echo "======================================================================"

        mkdir -p "$db_dir"

        # Run tronko-build
        START=$(date +%s)
        echo "  Command: $TRONKO_BUILD -y -e $merged_dir -n $n_clusters -d $db_dir -s -u $sp -a -E -c $FAMSA_THREADS"

        "$TRONKO_BUILD" -y \
            -e "$merged_dir" \
            -n "$n_clusters" \
            -d "$db_dir" \
            -s -u "$sp" \
            -a -E \
            -c "$FAMSA_THREADS" \
            2>&1 | tee "$db_dir/build.log"

        END=$(date +%s)
        BUILD_TIME=$((END - START))
        echo "  tronko-build completed in ${BUILD_TIME}s"

        # Count trees
        N_TREES=0
        if [[ -d "$db_dir/exported_subtrees" ]]; then
            N_TREES=$(ls "$db_dir/exported_subtrees"/*_MSA.fasta 2>/dev/null | wc -l | tr -d ' ')
        fi

        # Build BWA index
        echo "  Building marker.fasta + BWA index..."
        if [[ -f "$db_dir/reference_tree.txt" ]]; then
            # Extract leaf sequences for BWA
            python3 -c "
import os, glob
db_dir = '$db_dir'
subtree_dir = os.path.join(db_dir, 'exported_subtrees')
if not os.path.isdir(subtree_dir):
    subtree_dir = db_dir
out_fasta = os.path.join(db_dir, 'marker.fasta')
seen = set()
with open(out_fasta, 'w') as out:
    for msa in sorted(glob.glob(os.path.join(subtree_dir, '*_MSA.fasta'))):
        with open(msa) as f:
            name = None
            seq_parts = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name and name not in seen:
                        seq = ''.join(seq_parts).replace('-', '')
                        if seq:
                            out.write(f'>{name}\n{seq}\n')
                            seen.add(name)
                    name = line[1:].split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if name and name not in seen:
                seq = ''.join(seq_parts).replace('-', '')
                if seq:
                    out.write(f'>{name}\n{seq}\n')
                    seen.add(name)
print(f'  marker.fasta: {len(seen)} sequences')
"
            bwa index "$db_dir/marker.fasta" 2>&1 | tail -1
        fi

        # Copy input files for provenance
        cp "$INPUT_FASTA" "$db_dir/input.fasta" 2>/dev/null || true
        cp "$INPUT_TAXONOMY" "$db_dir/input_taxonomy.txt" 2>/dev/null || true

        # Write build_info.json
        cat > "$db_dir/build_info.json" << BIEOF
{
  "name": "$db_name",
  "pipeline": "build-tronko-db.sh (AncestralClust + FastTree + SP partitioning)",
  "marker": "12SV5",
  "input_sequences": 101088,
  "clustering": "ancestralclust_${config_name#ac_}",
  "sp_threshold": $sp,
  "sp_normalization": "corrected",
  "sp_normalization_note": "SP score = raw_score / numpairs, range [-2, +3]. Bug fix: removed extra /numspec division.",
  "tree_tool": "FastTree",
  "n_trees": $N_TREES,
  "build_time_seconds": $BUILD_TIME,
  "ablation_ready": true,
  "tronko_build_flags": "-y -e <merged> -n $n_clusters -d <out> -s -u $sp -a -E -c $FAMSA_THREADS"
}
BIEOF

        echo "  Done: $db_name — $N_TREES trees in ${BUILD_TIME}s"
    done
done

echo ""
echo "================================================================"
echo "Sweep complete!"
echo "================================================================"
