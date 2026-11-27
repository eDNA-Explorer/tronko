#!/usr/bin/env bash
# Build a tronko reference database from a tagged rcrux release and publish it.
#
# Steps (all by default; --steps picks some, e.g. to rebuild after editing a
# setting, or to publish a build already made):
#   fetch    curated/lca.fasta + lca_taxonomy.txt from the rcrux DVC remote at
#            --rcrux-tag
#   build    build-tronko-db.sh: ancestralclust above CLUSTER_CUTOFF sequences,
#            else one cluster; per cluster FAMSA → tree → nw_reroot;
#            tronko-build -y; final partitions → <marker>.fasta + taxonomy;
#            tronko-convert → reference_tree.trkb; bwa and bwa-mem3 index. Writes
#            params.json (every setting, and build-tronko-db.sh's own
#            build_metadata.json; published as build_params.json) and
#            build.log (all output).
#            tronko-build, tronko-convert and ancestralclust are compiled from
#            this checkout, which must be clean, so the recorded commit built
#            the tree.
#   publish  tronko_release.py publish: dvc push -r tronko, commit on a publish
#            branch of a second tronko checkout; merging it tags
#            <marker>_tronko_<version> in the tronko repository
#
# Settings come from scripts/tronko_build.conf, then
# scripts/tronko_build.d/<marker>.conf if it exists, then each --config FILE
# and --set KEY=VALUE / option in the order given. --print-config shows the
# result and where it came from.
#
# Usage (from a tronko checkout):
#   scripts/release_tronko.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0 \
#       --tronko-checkout ~/tronko-release [--work DIR] [--steps fetch,build,publish] \
#       [--config FILE]... [--set KEY=VALUE]... [--threads N] \
#       [--sp-threshold 0.1 | --min-leaf 500] [--cluster-cutoff N] \
#       [--tronko-version v2.0.0 | --bump minor] [--replace] [--push] [--print-config]
#
# Needs on PATH: famsa, the tree tool (VeryFastTree, FastTree or
# raxmlHPC-PTHREADS), nw_reroot, bwa, git, dvc; plus AWS credentials for the
# rcrux DVC remote (JS2 S3).

# Settings are assigned by name (printf -v from the config files).
# shellcheck disable=SC2153
set -euo pipefail

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit "${1:-0}"; }
say() { printf '\n==> %s\n' "$*"; }
run() { printf '+ %s\n' "$*"; "$@"; }

REPO="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=../lib/tronko-db.sh
. "$REPO/lib/tronko-db.sh"  # input checks; its die is replaced below
die() { echo "error: $*" >&2; exit 1; }

# Every setting a config file or --set may change; anything else is a typo.
SETTINGS="THREADS CLUSTER_CUTOFF ANCESTRALCLUST_BIN_SIZE ANCESTRALCLUST_MIN_DESCENDANTS \
ANCESTRALCLUST_ARGS CLUSTER_JOBS FAMSA_ARGS TREE_TOOL TREE_ARGS PARTITION SP_THRESHOLD MIN_LEAF \
SP_FORMULA TRONKO_BUILD_JOBS TRONKO_BUILD_ARGS BWA_INDEX_ARGS ASSIGNMENT_CHECK_COMMUNITIES \
ASSIGNMENT_CHECK_SPECIES ASSIGNMENT_CHECK_READS_PER_QUERY"
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
MARKER="" RCRUX_TAG="" TRONKO_CHECKOUT="" RCRUX_REPO="" RCRUX_POINTER="" FETCH_PRIMERS=() WORK="" STEPS="fetch,build,publish"
PRINT_CONFIG=0
OVERRIDES=()  # "config:FILE" or "set:KEY=VALUE"
PUBLISH_OPTS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --marker) MARKER="$2"; shift 2 ;;
        --rcrux-tag) RCRUX_TAG="$2"; shift 2 ;;
        --tronko-checkout) TRONKO_CHECKOUT="$2"; shift 2 ;;
        --rcrux-repo) RCRUX_REPO="$2"; shift 2 ;;
        --rcrux-pointer) RCRUX_POINTER="$2"; shift 2 ;;
        --forward-primer) FETCH_PRIMERS+=(--forward-primer "$2"); shift 2 ;;
        --reverse-primer) FETCH_PRIMERS+=(--reverse-primer "$2"); shift 2 ;;
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
# One tronko-build job per 8 threads; a THREADS that is not a number is reported below.
# An empty CLUSTER_JOBS stays empty: Step 2's thread shares set how many clusters run.
[ -n "$TRONKO_BUILD_JOBS" ] || [[ ! "$THREADS" =~ ^[0-9]+$ ]] \
    || TRONKO_BUILD_JOBS=$(( THREADS >= 16 ? THREADS / 8 : 1 ))
OVERRIDES_TEXT=""
[ "${#OVERRIDES[@]}" -eq 0 ] || OVERRIDES_TEXT="${OVERRIDES[*]}"

# --- Check the settings ---------------------------------------------------------
for key in THREADS CLUSTER_CUTOFF ANCESTRALCLUST_BIN_SIZE ANCESTRALCLUST_MIN_DESCENDANTS \
        MIN_LEAF TRONKO_BUILD_JOBS ASSIGNMENT_CHECK_COMMUNITIES ASSIGNMENT_CHECK_SPECIES ASSIGNMENT_CHECK_READS_PER_QUERY; do
    case "${!key}" in '' | *[!0-9]*) die "$key must be a whole number, got '${!key}'" ;; esac
done
case "$CLUSTER_JOBS" in *[!0-9]*) die "CLUSTER_JOBS must be a whole number, got '$CLUSTER_JOBS'" ;; esac
for key in THREADS ANCESTRALCLUST_BIN_SIZE TRONKO_BUILD_JOBS; do
    [ "${!key}" -ge 1 ] || die "$key must be at least 1"
done
[ -z "$CLUSTER_JOBS" ] || [ "$CLUSTER_JOBS" -ge 1 ] || die "CLUSTER_JOBS must be at least 1"
for key in ASSIGNMENT_CHECK_COMMUNITIES ASSIGNMENT_CHECK_SPECIES; do
    [ "$ASSIGNMENT_CHECK_READS_PER_QUERY" -eq 0 ] || [ "${!key}" -ge 1 ] \
        || die "$key must be at least 1; ASSIGNMENT_CHECK_READS_PER_QUERY=0 skips the check"
done
[[ "$SP_THRESHOLD" =~ ^[0-9]+(\.[0-9]+)?$ ]] || die "SP_THRESHOLD must be a number like 0.1, got '$SP_THRESHOLD'"
case "$PARTITION" in sum_of_pairs | min_leaf) ;; *) die "PARTITION must be sum_of_pairs or min_leaf" ;; esac
case "$TREE_TOOL" in veryfasttree | fasttree | raxml) ;;
    *) die "TREE_TOOL must be veryfasttree, fasttree or raxml" ;; esac
case "$SP_FORMULA" in legacy | corrected) ;; *) die "SP_FORMULA must be legacy or corrected" ;; esac
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
    [ -n "$TRONKO_CHECKOUT" ] || die "publish needs --tronko-checkout (a clean tronko checkout to commit the release on)"
    [ -d "$TRONKO_CHECKOUT/.git" ] || die "$TRONKO_CHECKOUT is not a git checkout"
    TRONKO_CHECKOUT="$(cd "$TRONKO_CHECKOUT" && pwd)"
fi

WORK="${WORK:-$PWD/tronko-release-$MARKER}"
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
INPUT="$WORK/input" OUT="$WORK/out"
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
    [ -z "$RCRUX_POINTER" ] || FETCH_ARGS+=(--rcrux-pointer "$RCRUX_POINTER")
    FETCH_ARGS+=(${FETCH_PRIMERS[@]+"${FETCH_PRIMERS[@]}"})
    run "${TOOL[@]}" fetch "${FETCH_ARGS[@]}"
fi

# ---------------------------------------------------------------------------
if has_step build; then
    say "Preflight"
    [ -s "$INPUT/inputs.json" ] || die "no fetched inputs in $INPUT; run the fetch step"
    TRONKO_SHA="$(git -C "$REPO" rev-parse HEAD)"
    # The check needs community-lib, vsearch 2.22.1 and tronko-assign. Its scoring goes
    # through polars, which crashes under amd64 emulation (an arm64 Mac): a DataFrame
    # here makes that stop now, not after the build.
    if [ "$ASSIGNMENT_CHECK_READS_PER_QUERY" -ne 0 ]; then
        skip="ASSIGNMENT_CHECK_READS_PER_QUERY=0 skips the assignment check"
        "${PY[@]}" -c 'import polars, community_lib.check; polars.DataFrame({"a": [1]})' > /dev/null 2>&1 \
            || die "the assignment check needs community-lib and a working polars (the release image on amd64 has both); $skip"
        case "$(vsearch --version 2>&1 | head -1)" in vsearch\ v2.22.1_*) ;;
            *) die "the assignment check needs vsearch 2.22.1 on PATH (the release image has it); $skip" ;; esac
        command -v tronko-assign > /dev/null || die "the assignment check needs tronko-assign on PATH; $skip"
    fi
    # BUILD_TRONKO_DB names a build script to run in place of this checkout's.
    BUILD_SCRIPT="${BUILD_TRONKO_DB:-}"
    if [ -z "$BUILD_SCRIPT" ]; then
        for tool in tronko-build tronko-convert ancestralclust; do
            run make -s -C "$REPO/$tool"
            PATH="$REPO/$tool:$PATH"
        done
        PATH="$PATH:$REPO/bin"  # fasta2phyml.pl
        # BWA-MEM3 is a pinned outside release, not this checkout's code.
        if ! command -v bwa-mem3 >/dev/null; then
            run make -s -C "$REPO/tronko-assign" bwa-mem3
            PATH="$REPO/tronko-assign:$PATH"
        fi
        BUILD_SCRIPT="$REPO/build-tronko-db.sh"
    fi

    # The names fetch saved the tier's files under, and the tag it came from.
    input_field() { sed -n "s/.*\"$1\": \"\([^\"]*\)\".*/\1/p" "$INPUT/inputs.json"; }
    FASTA="$INPUT/$(input_field fasta)" SOURCE_TAXONOMY="$INPUT/$(input_field taxonomy)"
    [ -f "$FASTA" ] && [ -f "$SOURCE_TAXONOMY" ] || die "$INPUT does not hold the files inputs.json names"

    # rCRUX writes accession, taxid, lineage; the build takes it as it is, with the
    # taxdump, and writes the taxid files. The checks below read accession, lineage.
    TAXONOMY="$WORK/taxonomy.2col.txt"
    TAXONOMY_COLUMNS="$(awk -F '\t' '{ print NF }' "$SOURCE_TAXONOMY" | sort -u)"
    case "$TAXONOMY_COLUMNS" in
        2) cp "$SOURCE_TAXONOMY" "$TAXONOMY" ;;
        3) cut -f1,3 "$SOURCE_TAXONOMY" > "$TAXONOMY" ;;
        *) die "$SOURCE_TAXONOMY must have 2 or 3 tab-separated columns on every line" ;;
    esac
    # Checked now, so a bad input stops in seconds and not after the build.
    check_taxonomy_ranks "$TAXONOMY" 2 || die "every lineage needs $TAX_RANKS ranks"
    check_accession_charset "$TAXONOMY" || die "accessions cannot contain ':' or spaces"
    [ -z "$(cut -f1 "$TAXONOMY" | sort | uniq -d | head -1)" ] || die "the taxonomy repeats an accession"
    check_accession_sets "$FASTA" "$TAXONOMY" || die "the FASTA and the taxonomy cover different accessions"

    BUILD_ARGS=(-f "$FASTA" -t "$SOURCE_TAXONOMY" -o "$OUT" -p "$MARKER" -T "$THREADS"
        --rcrux-tag "$(input_field rcrux_tag)"
        -C "$CLUSTER_CUTOFF" -B "$ANCESTRALCLUST_BIN_SIZE" -P "$ANCESTRALCLUST_MIN_DESCENDANTS"
        ${CLUSTER_JOBS:+-J "$CLUSTER_JOBS"} --tree-tool "$TREE_TOOL" --tb-jobs "$TRONKO_BUILD_JOBS"
        --ancestralclust-args "$ANCESTRALCLUST_ARGS" --famsa-args "$FAMSA_ARGS"
        --tree-args "$TREE_ARGS" --tronko-build-args "$TRONKO_BUILD_ARGS"
        --bwa-args "$BWA_INDEX_ARGS")
    if [ "$PARTITION" = sum_of_pairs ]; then
        BUILD_ARGS+=(-s "$SP_THRESHOLD")
    else
        BUILD_ARGS+=(--min-leaf "$MIN_LEAF")
    fi
    [ "$SP_FORMULA" = corrected ] || BUILD_ARGS+=(-L)
    [ "$TAXONOMY_COLUMNS" = 2 ] || BUILD_ARGS+=(--taxdump "$INPUT/taxdump")

    say "Build"
    rm -rf "$OUT"
    : > "$LOG"
    run bash "$BUILD_SCRIPT" "${BUILD_ARGS[@]}" 2>&1 | tee -a "$LOG"
    [ -s "$OUT/reference_tree.trkb" ] || die "build failed; see $LOG"
    # Sets and counts: a count alone hides one sequence dropped and another repeated.
    n_in="$(grep -c '^>' "$FASTA")" n_out="$(grep -c '^>' "$OUT/$MARKER.fasta")"
    [ "$n_in" = "$n_out" ] || die "only $n_out of $n_in sequences reached $OUT/$MARKER.fasta"
    check_accession_sets "$OUT/$MARKER.fasta" "$TAXONOMY" || die "$OUT/$MARKER.fasta does not hold the input's accessions"
    if [ "$TAXONOMY_COLUMNS" = 3 ]; then
        for name in taxids lineage_taxids; do
            [ -s "$OUT/${MARKER}_$name.txt" ] || die "the build wrote no ${MARKER}_$name.txt"
        done
    fi

    # Every effective setting, and what the build recorded about itself, so
    # params.json says exactly how this tree was made.
    for key in $SETTINGS; do export "TB_SETTING_$key=${!key}"; done
    TB_SETTINGS="$SETTINGS" TB_CONFIG_FILES="${CONFIG_FILES[*]}" \
        TB_OVERRIDES="$OVERRIDES_TEXT" TB_TRONKO_SHA="$TRONKO_SHA" \
        TB_INPUTS="$INPUT/inputs.json" TB_BUILD="$OUT/build_metadata.json" \
        "${PY[@]}" - "$PARAMS" "$REPO" <<'PY'
import json, os, sys
from pathlib import Path

env = os.environ
repo = Path(sys.argv[2])
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
        "config_files": config_files,
        "overrides": env["TB_OVERRIDES"].split(),
        "settings": {k: env[f"TB_SETTING_{k}"] for k in env["TB_SETTINGS"].split()},
        "build": json.load(open(env["TB_BUILD"])),
    },
    open(sys.argv[1], "w"),
    indent=2,
)
PY

    if [ "$ASSIGNMENT_CHECK_READS_PER_QUERY" -ne 0 ]; then
        say "Assignment check"
        # The holdout references are rebuilt from the files a release publishes, laid
        # out as a pull of it, so a published release is known to be enough for them.
        rm -rf "$WORK/assignment_check" "$WORK/holdout_master"
        run "${TOOL[@]}" holdout-master --marker "$MARKER" --build-dir "$OUT" --out "$WORK/holdout_master"
        run "${PY[@]}" "$REPO/scripts/assignment_check.py" --build-dir "$WORK/holdout_master" --marker "$MARKER" \
            --inputs "$INPUT/inputs.json" --threads "$THREADS" --work "$WORK/assignment_check" \
            --communities "$ASSIGNMENT_CHECK_COMMUNITIES" --richness "$ASSIGNMENT_CHECK_SPECIES" \
            --reads-per-query "$ASSIGNMENT_CHECK_READS_PER_QUERY" --out "$OUT/assignment_check.json" 2>&1 | tee -a "$LOG"
    fi
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
        --tronko-checkout "$TRONKO_CHECKOUT" ${PUBLISH_OPTS[@]+"${PUBLISH_OPTS[@]}"}
fi
