#!/usr/bin/env bash
# Build a tronko reference database from a tagged rcrux release and publish it,
# end to end.
#
#   1. preflight: clean tronko checkout; tronko-build and tronko-convert compiled
#      from it, so the recorded tronko commit is the code that built the tree
#   2. fetch: curated/lca.fasta + lca_taxonomy.txt from the rcrux DVC remote
#      at --rcrux-tag
#   3. build (the sequence edna-explorer-data-pipelines' build_from_scratch runs):
#      ancestralclust above --cluster-cutoff sequences, else one cluster;
#      per cluster FAMSA → VeryFastTree/FastTree → nw_reroot; tronko-build -y;
#      final partitions → <marker>.fasta + <marker>_taxonomy.txt;
#      tronko-convert → reference_tree.trkb; bwa index.
#      Every setting goes to params.json and all output to build.log.
#   4. publish: tronko_release.py publish (dvc push -r rcrux, commit on a
#      publish branch of the rcrux-py checkout); merging it tags
#      <marker>_tronko_<version>
#
# Usage (from a tronko checkout):
#   scripts/release_tronko.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0 \
#       --rcrux-checkout ~/rcrux-py [--work DIR] [--threads 32] \
#       [--sp-threshold 0.1 | --min-leaf 500] [--cluster-cutoff 25000] \
#       [--min-cluster-size 50] [--max-cluster-size 25000] [--min-descendants 75] \
#       [--bin-size 20000] [--tronko-version v2.0.0 | --bump minor] [--replace] [--push]
#
# Needs on PATH: famsa, VeryFastTree or FastTree, nw_reroot, bwa, and
# ancestralclust for inputs above --cluster-cutoff; plus AWS credentials for
# the rcrux DVC remote (JS2 S3).
set -euo pipefail

usage() { sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

MARKER="" RCRUX_TAG="" RCRUX_CHECKOUT="" WORK="" THREADS=""
SP_THRESHOLD="0.1" MIN_LEAF="" CUTOFF=25000 MIN_CLUSTER=50 MAX_CLUSTER=25000
MIN_DESCENDANTS=75 BIN_SIZE=20000 RCRUX_REPO=""
PUBLISH_OPTS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --marker) MARKER="$2"; shift 2 ;;
        --rcrux-tag) RCRUX_TAG="$2"; shift 2 ;;
        --rcrux-checkout) RCRUX_CHECKOUT="$2"; shift 2 ;;
        --rcrux-repo) RCRUX_REPO="$2"; shift 2 ;;
        --work) WORK="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --sp-threshold) SP_THRESHOLD="$2"; shift 2 ;;
        --min-leaf) MIN_LEAF="$2"; shift 2 ;;
        --cluster-cutoff) CUTOFF="$2"; shift 2 ;;
        --min-cluster-size) MIN_CLUSTER="$2"; shift 2 ;;
        --max-cluster-size) MAX_CLUSTER="$2"; shift 2 ;;
        --min-descendants) MIN_DESCENDANTS="$2"; shift 2 ;;
        --bin-size) BIN_SIZE="$2"; shift 2 ;;
        --tronko-version) PUBLISH_OPTS+=(--tronko-version "$2"); shift 2 ;;
        --bump) PUBLISH_OPTS+=(--bump "$2"); shift 2 ;;
        --replace) PUBLISH_OPTS+=(--replace); shift ;;
        --push) PUBLISH_OPTS+=(--push-branch); shift ;;
        -h|--help) usage 0 ;;
        *) echo "error: unknown option $1" >&2; usage 2 ;;
    esac
done

die() { echo "error: $*" >&2; exit 1; }
say() { printf '\n==> %s\n' "$*"; }
run() { printf '+ %s\n' "$*"; "$@"; }

if [ -z "$MARKER" ] || [ -z "$RCRUX_TAG" ] || [ -z "$RCRUX_CHECKOUT" ]; then
    die "--marker, --rcrux-tag and --rcrux-checkout are required"
fi
[ -d "$RCRUX_CHECKOUT/.git" ] || die "$RCRUX_CHECKOUT is not an rcrux-py checkout"
RCRUX_CHECKOUT="$(cd "$RCRUX_CHECKOUT" && pwd)"

REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-$PWD/tronko-release-$MARKER}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
if [ -z "$THREADS" ]; then
    THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
fi
if [ -n "${PYTHON:-}" ]; then
    PY=("$PYTHON")
elif command -v uv >/dev/null; then
    PY=(uv run --quiet --project "$REPO/scripts" python)
else
    PY=("$REPO/scripts/.venv/bin/python")
fi
TOOL=("${PY[@]}" "$REPO/scripts/tronko_release.py")

# ---------------------------------------------------------------------------
say "Preflight"
[ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] \
    || die "$REPO has uncommitted changes; the tree must be built from a commit"
TRONKO_SHA="$(git -C "$REPO" rev-parse HEAD)"
if [ -z "${TRONKO_BUILD:-}" ]; then
    run make -s -C "$REPO/tronko-build"
    TRONKO_BUILD="$REPO/tronko-build/tronko-build"
fi
if [ -z "${TRONKO_CONVERT:-}" ]; then
    run make -s -C "$REPO/tronko-convert"
    TRONKO_CONVERT="$REPO/tronko-convert/tronko-convert"
fi
for tool in famsa nw_reroot bwa git dvc; do
    command -v "$tool" >/dev/null || die "$tool not found on PATH"
done
if command -v VeryFastTree >/dev/null; then
    TREE_BUILDER="VeryFastTree"
elif command -v FastTree >/dev/null; then
    TREE_BUILDER="FastTree"
elif command -v fasttree >/dev/null; then
    TREE_BUILDER="fasttree"
else
    die "VeryFastTree or FastTree not found on PATH"
fi
echo "tronko $TRONKO_SHA; tree builder $TREE_BUILDER; $THREADS threads; work dir $WORK"

# ---------------------------------------------------------------------------
say "Fetch $RCRUX_TAG"
INPUT="$WORK/input"
FETCH_ARGS=(--marker "$MARKER" --rcrux-tag "$RCRUX_TAG" --out "$INPUT")
[ -z "$RCRUX_REPO" ] || FETCH_ARGS+=(--rcrux-repo "$RCRUX_REPO")
run "${TOOL[@]}" fetch "${FETCH_ARGS[@]}"

# ---------------------------------------------------------------------------
BUILD="$WORK/build"
OUT="$WORK/out"
LOG="$WORK/build.log"
rm -rf "$BUILD" "$OUT"
mkdir -p "$BUILD/clusters" "$BUILD/merged" "$BUILD/tronko-build" "$OUT"
: > "$LOG"

build() {
    local fasta="$INPUT/lca.fasta" taxonomy="$INPUT/lca_taxonomy.txt"
    N_SEQS="$(grep -c '^>' "$fasta")"
    echo "tronko build of $MARKER from $RCRUX_TAG: $N_SEQS sequences, tronko $TRONKO_SHA"

    say "Cluster"
    if [ "$N_SEQS" -gt "$CUTOFF" ]; then
        command -v ancestralclust >/dev/null || die "ancestralclust not found; needed above $CUTOFF sequences"
        local bins=$(( (N_SEQS + BIN_SIZE - 1) / BIN_SIZE ))
        local seeds=$(( bins * 30 < 4000 ? bins * 30 : 4000 ))
        echo "ancestralclust: $bins bins, $seeds seeds (bin size $BIN_SIZE)"
        # -u (Needleman-Wunsch) balances clusters better; its code path is not
        # thread-safe, hence -c 1 (as in build_from_scratch).
        run ancestralclust -f -u -i "$fasta" -b "$bins" -r "$seeds" -p "$MIN_DESCENDANTS" \
            -c 1 -d "$BUILD/clusters/" -m "$MIN_CLUSTER" -x "$MAX_CLUSTER"
    else
        echo "$N_SEQS sequences, at most $CUTOFF: one cluster"
        cp "$fasta" "$BUILD/clusters/0.fasta"
    fi

    say "Align, tree and root each cluster"
    N_CLUSTERS=0
    for cluster in "$BUILD"/clusters/*.fasta; do
        [ -e "$cluster" ] || die "clustering produced no clusters"
        local i="$N_CLUSTERS" work="$BUILD/clusters/$N_CLUSTERS"
        mkdir -p "$work"
        # The cluster's taxonomy rows, by FASTA header.
        awk -F '\t' 'NR == FNR { tax[$1] = $2; next }
            /^>/ { id = substr($1, 2); sub(/[ \t].*/, "", id); if (id in tax) print id "\t" tax[id] }' \
            "$taxonomy" "$cluster" > "$BUILD/merged/${i}_taxonomy.txt"
        run famsa -t "$THREADS" "$cluster" "$work/msa.raw.fasta"
        # tronko-build needs one line per sequence; FAMSA wraps them.
        awk '/^>/ { if (seq != "") print seq; print; seq = ""; next } { seq = seq $0 }
            END { if (seq != "") print seq }' "$work/msa.raw.fasta" > "$BUILD/merged/${i}_MSA.fasta"
        if [ "$TREE_BUILDER" = "VeryFastTree" ]; then
            VeryFastTree -gtr -gamma -nt -threads "$THREADS" -nosupport -fastexp 2 \
                "$BUILD/merged/${i}_MSA.fasta" > "$work/tree.unrooted"
        else
            "$TREE_BUILDER" -nt -gtr -gamma "$BUILD/merged/${i}_MSA.fasta" > "$work/tree.unrooted"
        fi
        if ! nw_reroot "$work/tree.unrooted" > "$BUILD/merged/RAxML_bestTree.${i}.reroot" \
                || [ ! -s "$BUILD/merged/RAxML_bestTree.${i}.reroot" ]; then
            echo "warning: nw_reroot failed for cluster $i; using the unrooted tree"
            cp "$work/tree.unrooted" "$BUILD/merged/RAxML_bestTree.${i}.reroot"
        fi
        N_CLUSTERS=$(( N_CLUSTERS + 1 ))
    done
    echo "$N_CLUSTERS cluster(s)"

    say "tronko-build"
    local partition=(-s -u "$SP_THRESHOLD")
    [ -z "$MIN_LEAF" ] || partition=(-v -f "$MIN_LEAF")
    run "$TRONKO_BUILD" -y -e "$BUILD/merged" -n "$N_CLUSTERS" -d "$BUILD/tronko-build" \
        "${partition[@]}" -a -c "$THREADS"

    say "Final partitions → $MARKER.fasta, ${MARKER}_taxonomy.txt"
    local finals="$BUILD/tronko-build/final_partitions.txt"
    [ -s "$finals" ] || die "tronko-build wrote no final_partitions.txt"
    : > "$OUT/$MARKER.fasta"
    : > "$OUT/${MARKER}_taxonomy.txt"
    while read -r p; do
        [ -n "$p" ] || continue
        cat "$BUILD/tronko-build/partition$p.fasta" >> "$OUT/$MARKER.fasta"
        cat "$BUILD/tronko-build/partition${p}_taxonomy.txt" >> "$OUT/${MARKER}_taxonomy.txt"
    done < "$finals"

    say "reference_tree.trkb and BWA index"
    run "$TRONKO_CONVERT" -i "$BUILD/tronko-build/reference_tree.txt" -o "$OUT/reference_tree.trkb" -v
    run bwa index "$OUT/$MARKER.fasta"
}
build 2>&1 | tee -a "$LOG"
# The build ran in a pipeline subshell; read back what it decided.
N_SEQS="$(grep -c '^>' "$INPUT/lca.fasta")"
N_CLUSTERS="$(find "$BUILD/merged" -name '*_MSA.fasta' | wc -l | tr -d ' ')"
CLUSTERED=false; [ "$N_SEQS" -le "$CUTOFF" ] || CLUSTERED=true
[ -s "$OUT/reference_tree.trkb" ] || die "build failed; see $LOG"

# Every setting the build used, written here so it is never left out.
PARAMS="$WORK/params.json"
"${PY[@]}" - "$PARAMS" <<PY
import json, sys
partition = {"mode": "min_leaf", "min_leaf": int("$MIN_LEAF")} if "$MIN_LEAF" else {
    "mode": "sum_of_pairs", "sp_threshold": float("$SP_THRESHOLD")}
json.dump({
    "script": "scripts/release_tronko.sh",
    "tronko_sha": "$TRONKO_SHA",
    "rcrux_tag": "$RCRUX_TAG",
    "input_tier": "lca",
    "input_sequences": $N_SEQS,
    "clusters": $N_CLUSTERS,
    "ancestralclust": {
        "used": $( [ "$CLUSTERED" = true ] && echo True || echo False ),
        "cutoff": $CUTOFF, "bin_size": $BIN_SIZE, "min_descendants": $MIN_DESCENDANTS,
        "min_cluster_size": $MIN_CLUSTER, "max_cluster_size": $MAX_CLUSTER,
        "flags": ["-f", "-u", "-c", "1"],
    },
    "famsa": {"threads": $THREADS},
    "tree_builder": "$TREE_BUILDER",
    "reroot": "nw_reroot",
    "tronko_build": {"flags": ["-y", "-a"], "partition": partition, "threads": $THREADS},
    "tronko_convert": "reference_tree.txt -> reference_tree.trkb",
    "bwa": "bwa index",
}, open(sys.argv[1], "w"), indent=2)
PY

# ---------------------------------------------------------------------------
say "Publish"
run "${TOOL[@]}" publish --marker "$MARKER" --build-dir "$OUT" --inputs "$INPUT/inputs.json" \
    --build-log "$LOG" --build-params "$PARAMS" --tronko-sha "$TRONKO_SHA" \
    --rcrux-checkout "$RCRUX_CHECKOUT" ${PUBLISH_OPTS[@]+"${PUBLISH_OPTS[@]}"}
