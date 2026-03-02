#!/bin/bash
# Build a tronko reference database replicating the CruxV2 pipeline.
#
# Pipeline steps (based on CALeDNA/crux tronko/build/ scripts):
#   1. AncestralClust pre-clustering (~20K seqs per cluster)
#   2. Per-cluster: FAMSA alignment + RAxML tree + nw_reroot
#   3. Merge/renumber clusters into tronko-build input directory
#   4. tronko-build -y -s (with recursive SP-score partitioning)
#   5. BWA index the reference FASTA
#
# Usage:
#   build-tronko-db.sh -f <input.fasta> -t <taxonomy.txt> -o <output_dir> \
#                      [-p <primer_name>] [-T <threads>] [-s <sp_threshold>] [-F]
#
# Options:
#   -f  Input reference FASTA (unaligned, single-line)
#   -t  Taxonomy file (accession<TAB>lineage)
#   -o  Output directory (will contain reference_tree.txt, FASTA, BWA index)
#   -p  Primer/marker name (default: marker)
#   -T  Threads for RAxML/FAMSA (default: 8)
#   -s  SP-score threshold for tronko-build partitioning (default: 0.1)
#   -F  Use FastTree instead of RAxML (faster, slightly less accurate)
#   -C  AncestralClust cutoff — max seqs before clustering (default: 25000)
#   -B  AncestralClust bin size — target seqs per cluster (default: 20000)

set -euo pipefail

# Defaults
PRIMER="marker"
THREADS=8
SP_THRESHOLD=0.1
USE_FASTTREE=0
AC_CUTOFF=25000
AC_BIN_SIZE=20000

while getopts "f:t:o:p:T:s:FC:B:" opt; do
    case $opt in
        f) INPUT_FASTA="$OPTARG" ;;
        t) INPUT_TAXONOMY="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) PRIMER="$OPTARG" ;;
        T) THREADS="$OPTARG" ;;
        s) SP_THRESHOLD="$OPTARG" ;;
        F) USE_FASTTREE=1 ;;
        C) AC_CUTOFF="$OPTARG" ;;
        B) AC_BIN_SIZE="$OPTARG" ;;
        *) echo "Usage: $0 -f <fasta> -t <taxonomy> -o <outdir> [-p primer] [-T threads] [-s sp_threshold] [-F] [-C cutoff] [-B binsize]" >&2; exit 1 ;;
    esac
done

# Validate required args
if [[ -z "${INPUT_FASTA:-}" || -z "${INPUT_TAXONOMY:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    echo "ERROR: -f, -t, and -o are required" >&2
    exit 1
fi

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "ERROR: Input FASTA not found: $INPUT_FASTA" >&2
    exit 1
fi

if [[ ! -f "$INPUT_TAXONOMY" ]]; then
    echo "ERROR: Taxonomy file not found: $INPUT_TAXONOMY" >&2
    exit 1
fi

# Verify tools
for tool in ancestralclust famsa nw_reroot tronko-build bwa fasta2phyml.pl; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found on PATH" >&2
        exit 1
    fi
done

# Find raxml binary (could be raxmlHPC-PTHREADS, raxmlHPC-PTHREADS-SSE3, etc.)
RAXML_BIN=""
for name in raxmlHPC-PTHREADS-SSE3 raxmlHPC-PTHREADS-AVX2 raxmlHPC-PTHREADS raxmlHPC; do
    if command -v "$name" &>/dev/null; then
        RAXML_BIN="$name"
        break
    fi
done

if [[ "$USE_FASTTREE" -eq 0 ]]; then
    if [[ -z "$RAXML_BIN" ]]; then
        echo "ERROR: No raxmlHPC variant found on PATH (use -F for FastTree)" >&2
        exit 1
    fi
    echo "Using RAxML: $RAXML_BIN"
else
    if command -v VeryFastTree &>/dev/null; then
        TREE_BIN="VeryFastTree"
        TREE_THREAD_FLAG="-threads $THREADS -nosupport -fastexp 2"
        echo "Using VeryFastTree: $(command -v VeryFastTree)"
    elif command -v FastTree &>/dev/null; then
        TREE_BIN="FastTree"
        TREE_THREAD_FLAG=""
        echo "Using FastTree: $(command -v FastTree) (consider installing VeryFastTree for 2-7x speedup)"
    else
        echo "ERROR: Neither VeryFastTree nor FastTree found on PATH" >&2
        exit 1
    fi
fi

# Setup directories — use $OUTPUT_DIR/.cache for checkpoint persistence
CACHE_DIR="$OUTPUT_DIR/.cache"
AC_DIR="$CACHE_DIR/ancestralclust"
NEWICK_DIR="$CACHE_DIR/newick"
MERGED_DIR="$CACHE_DIR/merged"
mkdir -p "$AC_DIR" "$NEWICK_DIR" "$MERGED_DIR" "$OUTPUT_DIR"

# No cleanup trap — cache persists for resume on failure
# To force a clean rebuild, delete $OUTPUT_DIR/.cache/

NUM_SEQS=$(grep -c '^>' "$INPUT_FASTA")
echo "================================================================"
echo "  tronko-build CruxV2 pipeline"
echo "================================================================"
echo "  Input:      $INPUT_FASTA ($NUM_SEQS sequences)"
echo "  Taxonomy:   $INPUT_TAXONOMY"
echo "  Output:     $OUTPUT_DIR"
echo "  Threads:    $THREADS"
echo "  SP thresh:  $SP_THRESHOLD"
echo "  Tree tool:  $([ "$USE_FASTTREE" -eq 1 ] && echo "${TREE_BIN:-FastTree}" || echo "RAxML")"
echo "  AC cutoff:  $AC_CUTOFF (cluster if > this many seqs)"
echo "  AC binsize: $AC_BIN_SIZE (target seqs per cluster)"
echo "================================================================"

# Helper: unwrap multi-line FASTA to single-line
unwrap_fasta() {
    local infile="$1"
    local outfile="$2"
    awk '/^>/ { if (seq) print seq; print; seq="" ; next } { seq = seq $0 } END { if (seq) print seq }' "$infile" > "$outfile"
}

# =========================================================================
# STEP 1: AncestralClust pre-clustering
# =========================================================================
echo ""
echo "=== Step 1/5: AncestralClust pre-clustering ==="
STEP1_START=$(date +%s)

if [[ -f "$AC_DIR/.step1_done" ]]; then
    NUM_CLUSTERS=$(ls "$AC_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
    echo "  CACHED: Step 1 already complete ($NUM_CLUSTERS clusters). Skipping."
elif [[ "$NUM_SEQS" -le "$AC_CUTOFF" ]]; then
    echo "Only $NUM_SEQS sequences (below cutoff $AC_CUTOFF) — skipping clustering"
    # Single cluster: just copy
    cp "$INPUT_FASTA" "$AC_DIR/0.fasta"
    # Extract taxonomy for these sequences
    cp "$INPUT_TAXONOMY" "$AC_DIR/0_taxonomy.txt"
    NUM_CLUSTERS=1
    touch "$AC_DIR/.step1_done"
else
    # Calculate AncestralClust parameters (CruxV2 formula)
    NUM_BINS=$(( (NUM_SEQS + AC_BIN_SIZE - 1) / AC_BIN_SIZE ))
    NUM_SEEDS=$(( NUM_BINS * 30 ))
    if [[ "$NUM_SEEDS" -gt 4000 ]]; then
        NUM_SEEDS=4000
    fi

    echo "AncestralClust: $NUM_SEQS seqs -> $NUM_BINS bins, $NUM_SEEDS seeds"
    echo "Running: ancestralclust -f -i $INPUT_FASTA -b $NUM_BINS -r $NUM_SEEDS -p 75 -c $THREADS -d $AC_DIR"

    ancestralclust -f \
        -i "$INPUT_FASTA" \
        -b "$NUM_BINS" \
        -r "$NUM_SEEDS" \
        -p 75 \
        -c "$THREADS" \
        -d "$AC_DIR"

    # AncestralClust outputs {0..N-1}.fasta in the output directory.
    # We need to create matching taxonomy files for each cluster.
    echo "Splitting taxonomy by cluster..."
    for cluster_fasta in "$AC_DIR"/*.fasta; do
        cluster_id=$(basename "$cluster_fasta" .fasta)
        # Extract accessions from this cluster's FASTA
        grep '^>' "$cluster_fasta" | sed 's/^>//' | cut -d' ' -f1 | sort > "$CACHE_DIR/accs_${cluster_id}.txt"
        # Filter taxonomy to matching accessions
        sort "$INPUT_TAXONOMY" | join -t$'\t' "$CACHE_DIR/accs_${cluster_id}.txt" - > "$AC_DIR/${cluster_id}_taxonomy.txt"
        rm "$CACHE_DIR/accs_${cluster_id}.txt"
        n_seqs=$(grep -c '^>' "$cluster_fasta")
        n_tax=$(wc -l < "$AC_DIR/${cluster_id}_taxonomy.txt")
        echo "  Cluster $cluster_id: $n_seqs seqs, $n_tax taxonomy entries"
    done

    NUM_CLUSTERS=$(ls "$AC_DIR"/*.fasta 2>/dev/null | wc -l | tr -d ' ')
    touch "$AC_DIR/.step1_done"
fi

STEP1_END=$(date +%s)
echo "Step 1 complete: $NUM_CLUSTERS clusters in $((STEP1_END - STEP1_START))s"

# =========================================================================
# STEP 2: Per-cluster alignment + tree building
# =========================================================================
echo ""
echo "=== Step 2/5: Per-cluster FAMSA + tree building ==="
STEP2_START=$(date +%s)

build_cluster_tree() {
    local cluster_fasta="$1"
    local cluster_id=$(basename "$cluster_fasta" .fasta)
    local outdir="$2"

    local n_seqs
    n_seqs=$(grep -c '^>' "$cluster_fasta")

    # Deduplicate sequences by header (VeryFastTree requires unique names)
    local dedup_fasta="$outdir/${cluster_id}_dedup.fasta"
    awk '/^>/{h=$0; if(seen[h]++){next} p=1; print; next} p{print} /^>/{p=0}' "$cluster_fasta" > "$dedup_fasta"
    local n_dedup
    n_dedup=$(grep -c '^>' "$dedup_fasta")
    if [[ "$n_dedup" -ne "$n_seqs" ]]; then
        echo "  Cluster $cluster_id ($n_seqs seqs, $n_dedup unique): aligning..."
    else
        echo "  Cluster $cluster_id ($n_seqs seqs): aligning..."
    fi

    # Skip clusters that are too small after dedup
    if [[ "$n_dedup" -lt 3 ]]; then
        echo "  Cluster $cluster_id: only $n_dedup unique seqs, skipping"
        rm -f "$dedup_fasta"
        return 1
    fi

    # FAMSA alignment
    famsa -t "$THREADS" "$dedup_fasta" "$outdir/${cluster_id}_MSA.fasta" 2>/dev/null
    rm -f "$dedup_fasta"

    # Unwrap to single-line (tronko-build requirement)
    unwrap_fasta "$outdir/${cluster_id}_MSA.fasta" "$outdir/${cluster_id}_MSA_unwrapped.fasta"
    mv "$outdir/${cluster_id}_MSA_unwrapped.fasta" "$outdir/${cluster_id}_MSA.fasta"

    echo "  Cluster $cluster_id: building tree..."

    if [[ "$USE_FASTTREE" -eq 1 ]]; then
        # VeryFastTree or FastTree (outputs Newick to stdout)
        set +e
        $TREE_BIN -gtr -gamma -nt $TREE_THREAD_FLAG "$outdir/${cluster_id}_MSA.fasta" \
            > "$outdir/RAxML_bestTree.${cluster_id}.unrooted" 2>"$outdir/${cluster_id}_tree.log"
        local tree_rc=$?
        set -e
        if [[ "$tree_rc" -ne 0 ]] || [[ ! -s "$outdir/RAxML_bestTree.${cluster_id}.unrooted" ]]; then
            echo "  WARNING: $TREE_BIN failed for cluster $cluster_id (exit $tree_rc). Log:"
            tail -20 "$outdir/${cluster_id}_tree.log"
            return 1
        fi
    else
        # RAxML (CruxV2 uses raxmlHPC-PTHREADS with GTR+Gamma)
        local raxml_dir="$outdir/${cluster_id}_RAxML"
        mkdir -p "$raxml_dir"
        fasta2phyml.pl "$outdir/${cluster_id}_MSA.fasta"
        "$RAXML_BIN" --silent -m GTRGAMMA \
            -w "$(cd "$raxml_dir" && pwd)" \
            -n 1 -p 1234 -T "$THREADS" \
            -s "$(cd "$outdir" && pwd)/${cluster_id}_MSA.phymlAln" 2>/dev/null || true
        if [[ -f "$raxml_dir/RAxML_bestTree.1" ]]; then
            cp "$raxml_dir/RAxML_bestTree.1" "$outdir/RAxML_bestTree.${cluster_id}.unrooted"
        else
            echo "  WARNING: RAxML failed for cluster $cluster_id, falling back to ${TREE_BIN:-FastTree}"
            ${TREE_BIN:-FastTree} -gtr -gamma -nt ${TREE_THREAD_FLAG:-} "$outdir/${cluster_id}_MSA.fasta" \
                > "$outdir/RAxML_bestTree.${cluster_id}.unrooted" 2>/dev/null
        fi
    fi

    # Midpoint rooting
    set +e
    nw_reroot "$outdir/RAxML_bestTree.${cluster_id}.unrooted" \
        > "$outdir/RAxML_bestTree.${cluster_id}.reroot" 2>/dev/null
    local reroot_rc=$?
    set -e
    if [[ "$reroot_rc" -ne 0 ]] || [[ ! -s "$outdir/RAxML_bestTree.${cluster_id}.reroot" ]]; then
        echo "  WARNING: nw_reroot failed for cluster $cluster_id, using unrooted tree"
        cp "$outdir/RAxML_bestTree.${cluster_id}.unrooted" "$outdir/RAxML_bestTree.${cluster_id}.reroot"
    fi

    echo "  Cluster $cluster_id: done"
}

for cluster_fasta in "$AC_DIR"/*.fasta; do
    cluster_id=$(basename "$cluster_fasta" .fasta)

    # Skip empty or single-sequence clusters
    n_seqs=$(grep -c '^>' "$cluster_fasta")
    if [[ "$n_seqs" -lt 3 ]]; then
        echo "  Cluster $cluster_id: only $n_seqs seqs, skipping (need >=3 for tree)"
        continue
    fi

    # Skip if already cached
    if [[ -f "$NEWICK_DIR/RAxML_bestTree.${cluster_id}.reroot" ]] && [[ -s "$NEWICK_DIR/RAxML_bestTree.${cluster_id}.reroot" ]]; then
        echo "  Cluster $cluster_id ($n_seqs seqs): CACHED, skipping"
        continue
    fi

    set +e
    build_cluster_tree "$cluster_fasta" "$NEWICK_DIR"
    cluster_rc=$?
    set -e
    if [[ "$cluster_rc" -ne 0 ]]; then
        echo "  WARNING: Cluster $cluster_id failed, skipping"
        continue
    fi

    # Copy taxonomy (deduplicated to match FASTA dedup)
    if [[ -f "$AC_DIR/${cluster_id}_taxonomy.txt" ]]; then
        awk -F'\t' '!seen[$1]++' "$AC_DIR/${cluster_id}_taxonomy.txt" > "$NEWICK_DIR/${cluster_id}_taxonomy.txt"
    fi
done

STEP2_END=$(date +%s)
echo "Step 2 complete in $((STEP2_END - STEP2_START))s"

# =========================================================================
# STEP 3: Merge/renumber clusters for tronko-build
# =========================================================================
echo ""
echo "=== Step 3/5: Merge and renumber clusters ==="
STEP3_START=$(date +%s)

counter=0
for tree_file in "$NEWICK_DIR"/RAxML_bestTree.*.reroot; do
    # Extract original cluster ID from filename
    orig_id=$(basename "$tree_file" | sed 's/RAxML_bestTree\.\(.*\)\.reroot/\1/')

    # Check that all 3 files exist
    if [[ ! -f "$NEWICK_DIR/${orig_id}_MSA.fasta" ]] || [[ ! -f "$NEWICK_DIR/${orig_id}_taxonomy.txt" ]]; then
        echo "  Skipping cluster $orig_id (missing MSA or taxonomy)"
        continue
    fi

    cp "$tree_file" "$MERGED_DIR/RAxML_bestTree.${counter}.reroot"
    cp "$NEWICK_DIR/${orig_id}_MSA.fasta" "$MERGED_DIR/${counter}_MSA.fasta"
    cp "$NEWICK_DIR/${orig_id}_taxonomy.txt" "$MERGED_DIR/${counter}_taxonomy.txt"

    # Also keep original FASTA for tronko-build partition mode
    if [[ -f "$AC_DIR/${orig_id}.fasta" ]]; then
        cp "$AC_DIR/${orig_id}.fasta" "$MERGED_DIR/${counter}.fasta"
    fi

    counter=$((counter + 1))
done

NUM_FINAL_CLUSTERS=$counter
STEP3_END=$(date +%s)
echo "Step 3 complete: $NUM_FINAL_CLUSTERS clusters renumbered in $((STEP3_END - STEP3_START))s"

if [[ "$NUM_FINAL_CLUSTERS" -eq 0 ]]; then
    echo "ERROR: No valid clusters produced. Cannot build database." >&2
    exit 1
fi

# =========================================================================
# STEP 4: tronko-build
# =========================================================================
echo ""
echo "=== Step 4/5: tronko-build ==="
STEP4_START=$(date +%s)

TRONKO_FLAGS=""
if [[ "$USE_FASTTREE" -eq 1 ]]; then
    TRONKO_FLAGS="-a"
fi

if [[ "$NUM_FINAL_CLUSTERS" -gt 1 ]] || [[ "$NUM_FINAL_CLUSTERS" -eq 1 ]]; then
    echo "Running: tronko-build -y -e $MERGED_DIR -n $NUM_FINAL_CLUSTERS -d $OUTPUT_DIR -s -u $SP_THRESHOLD $TRONKO_FLAGS"
    time tronko-build -y \
        -e "$MERGED_DIR" \
        -n "$NUM_FINAL_CLUSTERS" \
        -d "$OUTPUT_DIR" \
        -s -u "$SP_THRESHOLD" \
        $TRONKO_FLAGS \
        -c "$THREADS"
fi

# Check for output
if [[ ! -f "$OUTPUT_DIR/reference_tree.txt" ]]; then
    echo "ERROR: tronko-build did not produce reference_tree.txt" >&2
    exit 1
fi

STEP4_END=$(date +%s)
echo "Step 4 complete in $((STEP4_END - STEP4_START))s"

# =========================================================================
# STEP 5: Finalize (concat FASTA + taxonomy, BWA index)
# =========================================================================
echo ""
echo "=== Step 5/5: Finalize (FASTA concat + BWA index) ==="
STEP5_START=$(date +%s)

# Concatenate all cluster FASTAs and taxonomy
cat "$MERGED_DIR"/*_MSA.fasta > "$OUTPUT_DIR/${PRIMER}.fasta" 2>/dev/null || \
    cat "$MERGED_DIR"/*.fasta > "$OUTPUT_DIR/${PRIMER}.fasta" 2>/dev/null || true
cat "$MERGED_DIR"/*_taxonomy.txt > "$OUTPUT_DIR/${PRIMER}_taxonomy.txt" 2>/dev/null || true

# BWA index
echo "Building BWA index..."
bwa index "$OUTPUT_DIR/${PRIMER}.fasta"

# Compress reference tree
gzip "$OUTPUT_DIR/reference_tree.txt"

STEP5_END=$(date +%s)
echo "Step 5 complete in $((STEP5_END - STEP5_START))s"

# =========================================================================
# Summary
# =========================================================================
TOTAL_TIME=$(( STEP5_END - STEP1_START ))
echo ""
echo "================================================================"
echo "  Build complete!"
echo "================================================================"
echo "  Total time:     ${TOTAL_TIME}s ($(( TOTAL_TIME / 60 ))m $(( TOTAL_TIME % 60 ))s)"
echo "  Clusters:       $NUM_FINAL_CLUSTERS"
echo "  Output files:"
ls -lh "$OUTPUT_DIR"/ | grep -v '^total' | sed 's/^/    /'
echo "================================================================"
