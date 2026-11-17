#!/usr/bin/env bash
# Build a tronko reference database from a tagged rcrux release and publish it.
#
# Steps (all by default; --steps picks some, e.g. to rebuild after editing a
# setting, or to publish a build already made):
#   fetch    curated/lca.fasta + lca_taxonomy.txt from the rcrux DVC remote at
#            --rcrux-tag
#   build    ancestralclust above CLUSTER_CUTOFF sequences, else one cluster;
#            per cluster FAMSA → VeryFastTree/FastTree → nw_reroot;
#            tronko-build -y; final partitions → <marker>.fasta + taxonomy;
#            tronko-convert → reference_tree.trkb; bwa index. Writes
#            params.json (every setting) and build.log (all output).
#            tronko-build and tronko-convert are compiled from this checkout,
#            which must be clean, so the recorded commit built the tree.
#   publish  tronko_release.py publish: dvc push -r rcrux, commit on a publish
#            branch of the rcrux-py checkout; merging it tags
#            <marker>_tronko_<version>
#
# Settings come from scripts/tronko_build.conf, then
# scripts/tronko_build.d/<marker>.conf if it exists, then each --config FILE
# and --set KEY=VALUE / option in the order given. --print-config shows the
# result and where it came from.
#
# Usage (from a tronko checkout):
#   scripts/release_tronko.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0 \
#       --rcrux-checkout ~/rcrux-py [--work DIR] [--steps fetch,build,publish] \
#       [--config FILE]... [--set KEY=VALUE]... [--threads N] \
#       [--sp-threshold 0.1 | --min-leaf 500] [--cluster-cutoff N] \
#       [--tronko-version v2.0.0 | --bump minor] [--replace] [--push] [--print-config]
#
# Needs on PATH: famsa, VeryFastTree or FastTree, nw_reroot (REROOT=true), bwa,
# ancestralclust above CLUSTER_CUTOFF, git, dvc; plus AWS credentials for the
# rcrux DVC remote (JS2 S3).

# Settings are assigned by name (printf -v from the config files).
# shellcheck disable=SC2153
set -euo pipefail

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit "${1:-0}"; }
die() { echo "error: $*" >&2; exit 1; }
say() { printf '\n==> %s\n' "$*"; }
run() { printf '+ %s\n' "$*"; "$@"; }

REPO="$(cd "$(dirname "$0")/.." && pwd)"

# Every setting a config file or --set may change; anything else is a typo.
SETTINGS="THREADS CLUSTER_CUTOFF ANCESTRALCLUST_BIN_SIZE ANCESTRALCLUST_SEEDS_PER_BIN \
ANCESTRALCLUST_MAX_SEEDS ANCESTRALCLUST_MIN_DESCENDANTS ANCESTRALCLUST_MIN_CLUSTER \
ANCESTRALCLUST_MAX_CLUSTER ANCESTRALCLUST_ARGS FAMSA_ARGS TREE_BUILDER VERYFASTTREE_ARGS \
FASTTREE_ARGS REROOT PARTITION SP_THRESHOLD MIN_LEAF TRONKO_BUILD_ARGS BWA_INDEX_ARGS"
for key in $SETTINGS; do printf -v "$key" '%s' ""; done
CONFIG_FILES=()

is_setting() { case " $SETTINGS " in *" $1 "*) return 0 ;; esac; return 1; }
set_setting() {  # KEY VALUE WHERE
    is_setting "$1" || die "unknown setting $1 ($3); see scripts/tronko_build.conf"
    printf -v "$1" '%s' "$2"
}
# Read KEY=VALUE lines; values are taken literally, never run as shell.
load_config() {
    local file="$1" line key value n=0
    [ -f "$file" ] || die "config file $file not found"
    while IFS= read -r line || [ -n "$line" ]; do
        n=$((n + 1))
        line="${line#"${line%%[![:space:]]*}"}"
        case "$line" in '' | '#'*) continue ;; esac
        case "$line" in *=*) ;; *) die "$file:$n: expected KEY=VALUE" ;; esac
        key="${line%%=*}"
        key="${key%"${key##*[![:space:]]}"}"
        value="${line#*=}"
        value="${value#"${value%%[![:space:]]*}"}"
        value="${value%"${value##*[![:space:]]}"}"
        case "$value" in \"*\" | \'*\') value="${value:1:${#value}-2}" ;; esac
        set_setting "$key" "$value" "$file:$n"
    done < "$file"
    CONFIG_FILES+=("$file")
}

# --- Options (overrides applied in order, after the config files) -----------
MARKER="" RCRUX_TAG="" RCRUX_CHECKOUT="" RCRUX_REPO="" WORK="" STEPS="fetch,build,publish"
PRINT_CONFIG=0
OVERRIDES=()  # "config:FILE" or "set:KEY=VALUE"
PUBLISH_OPTS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --marker) MARKER="$2"; shift 2 ;;
        --rcrux-tag) RCRUX_TAG="$2"; shift 2 ;;
        --rcrux-checkout) RCRUX_CHECKOUT="$2"; shift 2 ;;
        --rcrux-repo) RCRUX_REPO="$2"; shift 2 ;;
        --work) WORK="$2"; shift 2 ;;
        --steps) STEPS="$2"; shift 2 ;;
        --config) OVERRIDES+=("config:$2"); shift 2 ;;
        --set) OVERRIDES+=("set:$2"); shift 2 ;;
        --threads) OVERRIDES+=("set:THREADS=$2"); shift 2 ;;
        --cluster-cutoff) OVERRIDES+=("set:CLUSTER_CUTOFF=$2"); shift 2 ;;
        --sp-threshold) OVERRIDES+=("set:PARTITION=sum_of_pairs" "set:SP_THRESHOLD=$2"); shift 2 ;;
        --min-leaf) OVERRIDES+=("set:PARTITION=min_leaf" "set:MIN_LEAF=$2"); shift 2 ;;
        --tronko-version) PUBLISH_OPTS+=(--tronko-version "$2"); shift 2 ;;
        --bump) PUBLISH_OPTS+=(--bump "$2"); shift 2 ;;
        --replace) PUBLISH_OPTS+=(--replace); shift ;;
        --push) PUBLISH_OPTS+=(--push-branch); shift ;;
        --print-config) PRINT_CONFIG=1; shift ;;
        -h | --help) usage 0 ;;
        *) echo "error: unknown option $1" >&2; usage 2 ;;
    esac
done

load_config "$REPO/scripts/tronko_build.conf"
if [ -n "$MARKER" ] && [ -f "$REPO/scripts/tronko_build.d/$MARKER.conf" ]; then
    load_config "$REPO/scripts/tronko_build.d/$MARKER.conf"
fi
for override in ${OVERRIDES[@]+"${OVERRIDES[@]}"}; do
    case "$override" in
        config:*) load_config "${override#config:}" ;;
        set:*=*)
            pair="${override#set:}"
            set_setting "${pair%%=*}" "${pair#*=}" "--set"
            ;;
        *) die "--set needs KEY=VALUE, got ${override#set:}" ;;
    esac
done
[ -n "$THREADS" ] || THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
OVERRIDES_TEXT=""
[ "${#OVERRIDES[@]}" -eq 0 ] || OVERRIDES_TEXT="${OVERRIDES[*]}"

# --- Check the settings ---------------------------------------------------------
for key in THREADS CLUSTER_CUTOFF ANCESTRALCLUST_BIN_SIZE ANCESTRALCLUST_SEEDS_PER_BIN \
        ANCESTRALCLUST_MAX_SEEDS ANCESTRALCLUST_MIN_DESCENDANTS ANCESTRALCLUST_MIN_CLUSTER \
        ANCESTRALCLUST_MAX_CLUSTER MIN_LEAF; do
    case "${!key}" in '' | *[!0-9]*) die "$key must be a whole number, got '${!key}'" ;; esac
done
case "$SP_THRESHOLD" in '' | *[!0-9.]* | *.*.*) die "SP_THRESHOLD must be a number, got '$SP_THRESHOLD'" ;; esac
case "$PARTITION" in sum_of_pairs | min_leaf) ;; *) die "PARTITION must be sum_of_pairs or min_leaf" ;; esac
case "$TREE_BUILDER" in auto | VeryFastTree | FastTree | fasttree) ;;
    *) die "TREE_BUILDER must be auto, VeryFastTree, FastTree or fasttree" ;; esac
case "$REROOT" in true | false) ;; *) die "REROOT must be true or false" ;; esac
for step in $(echo "$STEPS" | tr ',' ' '); do
    case "$step" in fetch | build | publish) ;; *) die "unknown step $step (fetch, build, publish)" ;; esac
done
has_step() { case ",$STEPS," in *",$1,"*) return 0 ;; esac; return 1; }

if [ "$PRINT_CONFIG" = 1 ]; then
    echo "# from: ${CONFIG_FILES[*]}${OVERRIDES_TEXT:+ + $OVERRIDES_TEXT}"
    for key in $SETTINGS; do echo "$key=${!key}"; done
    exit 0
fi

[ -n "$MARKER" ] || die "--marker is required"
if has_step fetch && [ -z "$RCRUX_TAG" ]; then die "fetch needs --rcrux-tag"; fi
if has_step publish; then
    [ -n "$RCRUX_CHECKOUT" ] || die "publish needs --rcrux-checkout"
    [ -d "$RCRUX_CHECKOUT/.git" ] || die "$RCRUX_CHECKOUT is not an rcrux-py checkout"
    RCRUX_CHECKOUT="$(cd "$RCRUX_CHECKOUT" && pwd)"
fi

WORK="${WORK:-$PWD/tronko-release-$MARKER}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
INPUT="$WORK/input" BUILD="$WORK/build" OUT="$WORK/out"
LOG="$WORK/build.log" PARAMS="$WORK/params.json"
if [ -n "${PYTHON:-}" ]; then
    PY=("$PYTHON")
elif command -v uv >/dev/null; then
    PY=(uv run --quiet --project "$REPO/scripts" python)
else
    PY=("$REPO/scripts/.venv/bin/python")
fi
TOOL=("${PY[@]}" "$REPO/scripts/tronko_release.py")
echo "settings from ${CONFIG_FILES[*]}; steps: $STEPS; work dir $WORK"
if has_step build; then
    # Checked first, before anything is downloaded: the tree must come from a commit.
    [ -z "$(git -C "$REPO" status --porcelain --untracked-files=no)" ] \
        || die "$REPO has uncommitted changes; the tree must be built from a commit"
fi

# ---------------------------------------------------------------------------
if has_step fetch; then
    say "Fetch $RCRUX_TAG"
    FETCH_ARGS=(--marker "$MARKER" --rcrux-tag "$RCRUX_TAG" --out "$INPUT")
    [ -z "$RCRUX_REPO" ] || FETCH_ARGS+=(--rcrux-repo "$RCRUX_REPO")
    run "${TOOL[@]}" fetch "${FETCH_ARGS[@]}"
fi

# ---------------------------------------------------------------------------
build() {
    local fasta="$INPUT/lca.fasta" taxonomy="$INPUT/lca_taxonomy.txt"
    local famsa_args=() tree_args=() tronko_args=() bwa_args=() cluster_args=()
    read -r -a famsa_args <<< "$FAMSA_ARGS"
    read -r -a tronko_args <<< "$TRONKO_BUILD_ARGS"
    read -r -a bwa_args <<< "$BWA_INDEX_ARGS"
    read -r -a cluster_args <<< "$ANCESTRALCLUST_ARGS"
    if [ "$TREE_BUILDER" = "VeryFastTree" ]; then
        read -r -a tree_args <<< "$VERYFASTTREE_ARGS"
        tree_args+=(-threads "$THREADS")
    else
        read -r -a tree_args <<< "$FASTTREE_ARGS"
    fi
    local n_seqs n_clusters=0
    n_seqs="$(grep -c '^>' "$fasta")"
    echo "tronko build of $MARKER: $n_seqs sequences, tronko $TRONKO_SHA, tree builder $TREE_BUILDER"

    say "Cluster"
    if [ "$n_seqs" -gt "$CLUSTER_CUTOFF" ]; then
        command -v ancestralclust >/dev/null || die "ancestralclust not found; needed above $CLUSTER_CUTOFF sequences"
        local bins=$(( (n_seqs + ANCESTRALCLUST_BIN_SIZE - 1) / ANCESTRALCLUST_BIN_SIZE ))
        local seeds=$(( bins * ANCESTRALCLUST_SEEDS_PER_BIN ))
        [ "$seeds" -le "$ANCESTRALCLUST_MAX_SEEDS" ] || seeds="$ANCESTRALCLUST_MAX_SEEDS"
        echo "ancestralclust: $bins bins, $seeds seeds"
        run ancestralclust ${cluster_args[@]+"${cluster_args[@]}"} -i "$fasta" -b "$bins" -r "$seeds" \
            -p "$ANCESTRALCLUST_MIN_DESCENDANTS" -d "$BUILD/clusters/" \
            -m "$ANCESTRALCLUST_MIN_CLUSTER" -x "$ANCESTRALCLUST_MAX_CLUSTER"
    else
        echo "$n_seqs sequences, at most $CLUSTER_CUTOFF: one cluster"
        cp "$fasta" "$BUILD/clusters/0.fasta"
    fi

    say "Align, tree and root each cluster"
    for cluster in "$BUILD"/clusters/*.fasta; do
        [ -e "$cluster" ] || die "clustering produced no clusters"
        local i="$n_clusters" work="$BUILD/clusters/$n_clusters"
        local msa="$BUILD/merged/${i}_MSA.fasta" tree="$BUILD/merged/RAxML_bestTree.${i}.reroot"
        mkdir -p "$work"
        # The cluster's taxonomy rows, by FASTA header.
        awk -F '\t' 'NR == FNR { tax[$1] = $2; next }
            /^>/ { id = substr($1, 2); sub(/[ \t].*/, "", id); if (id in tax) print id "\t" tax[id] }' \
            "$taxonomy" "$cluster" > "$BUILD/merged/${i}_taxonomy.txt"
        run famsa -t "$THREADS" ${famsa_args[@]+"${famsa_args[@]}"} "$cluster" "$work/msa.raw.fasta"
        # tronko-build needs one line per sequence; FAMSA wraps them.
        awk '/^>/ { if (seq != "") print seq; print; seq = ""; next } { seq = seq $0 }
            END { if (seq != "") print seq }' "$work/msa.raw.fasta" > "$msa"
        "$TREE_BUILDER" ${tree_args[@]+"${tree_args[@]}"} "$msa" > "$work/tree.unrooted"
        if [ "$REROOT" = true ] && nw_reroot "$work/tree.unrooted" > "$tree" && [ -s "$tree" ]; then
            :
        else
            [ "$REROOT" = false ] || echo "warning: nw_reroot failed for cluster $i; using the unrooted tree"
            cp "$work/tree.unrooted" "$tree"
        fi
        n_clusters=$(( n_clusters + 1 ))
    done
    echo "$n_clusters cluster(s)"

    say "tronko-build"
    local partition=(-s -u "$SP_THRESHOLD")
    [ "$PARTITION" = sum_of_pairs ] || partition=(-v -f "$MIN_LEAF")
    run "$TRONKO_BUILD" -y -e "$BUILD/merged" -n "$n_clusters" -d "$BUILD/tronko-build" \
        "${partition[@]}" ${tronko_args[@]+"${tronko_args[@]}"} -c "$THREADS"

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
    run bwa index ${bwa_args[@]+"${bwa_args[@]}"} "$OUT/$MARKER.fasta"
}

if has_step build; then
    say "Preflight"
    [ -s "$INPUT/inputs.json" ] || die "no fetched inputs in $INPUT; run the fetch step"
    TRONKO_SHA="$(git -C "$REPO" rev-parse HEAD)"
    if [ -z "${TRONKO_BUILD:-}" ]; then
        run make -s -C "$REPO/tronko-build"
        TRONKO_BUILD="$REPO/tronko-build/tronko-build"
    fi
    if [ -z "${TRONKO_CONVERT:-}" ]; then
        run make -s -C "$REPO/tronko-convert"
        TRONKO_CONVERT="$REPO/tronko-convert/tronko-convert"
    fi
    for tool in famsa bwa; do command -v "$tool" >/dev/null || die "$tool not found on PATH"; done
    [ "$REROOT" = false ] || command -v nw_reroot >/dev/null || die "nw_reroot not found on PATH"
    if [ "$TREE_BUILDER" = auto ]; then
        for candidate in VeryFastTree FastTree fasttree ""; do
            [ -n "$candidate" ] || die "VeryFastTree or FastTree not found on PATH"
            if command -v "$candidate" >/dev/null; then TREE_BUILDER="$candidate"; break; fi
        done
    fi
    command -v "$TREE_BUILDER" >/dev/null || die "$TREE_BUILDER not found on PATH"

    rm -rf "$BUILD" "$OUT"
    mkdir -p "$BUILD/clusters" "$BUILD/merged" "$BUILD/tronko-build" "$OUT"
    : > "$LOG"
    build 2>&1 | tee -a "$LOG"
    [ -s "$OUT/reference_tree.trkb" ] || die "build failed; see $LOG"

    # Every effective setting and what the build found, so params.json says
    # exactly how this tree was made.
    for key in $SETTINGS; do export "TB_SETTING_$key=${!key}"; done
    TB_SETTINGS="$SETTINGS" TB_CONFIG_FILES="${CONFIG_FILES[*]}" \
        TB_OVERRIDES="$OVERRIDES_TEXT" TB_TRONKO_SHA="$TRONKO_SHA" \
        TB_INPUTS="$INPUT/inputs.json" TB_TREE_BUILDER="$TREE_BUILDER" \
        TB_INPUT="$INPUT/lca.fasta" TB_MERGED="$BUILD/merged" \
        "${PY[@]}" - "$PARAMS" "$REPO" <<'PY'
import json, os, sys
from pathlib import Path

env = os.environ
repo = Path(sys.argv[2])
n_seqs = sum(1 for line in open(env["TB_INPUT"], "rb") if line.startswith(b">"))
settings = {k: env[f"TB_SETTING_{k}"] for k in env["TB_SETTINGS"].split()}
config_files = [
    os.path.relpath(f, repo) if f.startswith(str(repo)) else f
    for f in env["TB_CONFIG_FILES"].split()
]
json.dump(
    {
        "script": "scripts/release_tronko.sh",
        "tronko_sha": env["TB_TRONKO_SHA"],
        "rcrux_tag": json.load(open(env["TB_INPUTS"]))["rcrux_tag"],
        "input_tier": "lca",
        "input_sequences": n_seqs,
        "clusters": len(list(Path(env["TB_MERGED"]).glob("*_MSA.fasta"))),
        "clustered": n_seqs > int(settings["CLUSTER_CUTOFF"]),
        "tree_builder": env["TB_TREE_BUILDER"],
        "config_files": config_files,
        "overrides": env["TB_OVERRIDES"].split(),
        "settings": settings,
    },
    open(sys.argv[1], "w"),
    indent=2,
)
PY
fi

# ---------------------------------------------------------------------------
if has_step publish; then
    say "Publish"
    if [ ! -s "$PARAMS" ] || [ ! -s "$OUT/reference_tree.trkb" ]; then
        die "no finished build in $WORK; run the build step"
    fi
    # The commit that built the tree, as the build recorded it.
    BUILT_SHA="$(sed -n 's/.*"tronko_sha": "\([^"]*\)".*/\1/p' "$PARAMS")"
    run "${TOOL[@]}" publish --marker "$MARKER" --build-dir "$OUT" --inputs "$INPUT/inputs.json" \
        --build-log "$LOG" --build-params "$PARAMS" --tronko-sha "$BUILT_SHA" \
        --rcrux-checkout "$RCRUX_CHECKOUT" ${PUBLISH_OPTS[@]+"${PUBLISH_OPTS[@]}"}
fi
