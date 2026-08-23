#!/bin/bash
# Shared helpers for the tronko database build and verification scripts.
#
# Sourced by build-rcrux-marker-db.sh, verify-tronko-db.sh and (optionally)
# repair-orphaned-partitions.sh. Nothing here executes on its own.
#
# Every function returns non-zero on failure and prints to stderr; callers are
# expected to run under `set -e` or check explicitly.

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------
TRONKO_DB_LIB_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TRONKO_REPO_DIR=$(dirname "$TRONKO_DB_LIB_DIR")
RCRUX_REPO="${RCRUX_REPO:-$HOME/rcrux-py}"

# BWA index suffixes, in the order bwa writes them.
BWA_EXTS=(amb ann bwt pac sa)

# The taxonomy always carries these 7 ranks, semicolon-delimited.
TAX_RANKS=7

log()  { printf '%s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# status.txt
#
# watch_build.sh polls for a line matching '^EXIT' to know the build finished,
# so the START/EXIT words are a contract -- do not rename them.
# ---------------------------------------------------------------------------
db_status() {
    local db="$1" word="$2"
    printf '%s %s\n' "$word" "$(date -Is)" >> "$db/status.txt"
}

db_status_has() {
    local db="$1" word="$2"
    [ -f "$db/status.txt" ] && grep -q "^${word}\b" "$db/status.txt"
}

# ---------------------------------------------------------------------------
# Primer resolution
#
# The primer name determines the DB FASTA/taxonomy filenames. Prefer the value
# recorded at build time; fall back to the same heuristic the orphan scripts
# use (a .fasta that is not a partition, an alias, or the recorded input).
# ---------------------------------------------------------------------------
db_primer() {
    local db="$1" primer=""
    if [ -f "$db/build_metadata.json" ]; then
        primer=$(python3 -c "
import json,sys
try:
    m=json.load(open('$db/build_metadata.json')).get('marker','')
    print(m if m and m!='marker' else '')
except Exception:
    print('')
" 2>/dev/null)
    fi
    if [ -n "$primer" ] && [ -f "$db/${primer}.fasta" ]; then
        printf '%s' "$primer"; return 0
    fi
    local f
    f=$(find "$db" -maxdepth 1 -name '*.fasta' \
            ! -name 'partition*' ! -name 'marker.fasta' ! -name 'input.fasta' \
            -printf '%f\n' 2>/dev/null | sort | head -1)
    [ -n "$f" ] || return 1
    printf '%s' "${f%.fasta}"
}

# ---------------------------------------------------------------------------
# Harness FASTA/taxonomy resolution
#
# Mirrors TronkoDatabaseFiles.from_directory() in the benchmark harness
# (assignment_benchmarks/.../tronko/database.py). Every production caller omits
# marker=, so the sorted-glob fallback is the path actually taken. It skips BWA
# suffixes and partition* but NOT input.fasta -- so a primer name sorting after
# "input" silently yields the raw unaligned rCRUX input, with a BWA index that
# does not match and leaf coordinates unrelated to reference_tree.trkb.
#
# Prints "<fasta_basename> <taxonomy_basename>".
# ---------------------------------------------------------------------------
harness_resolution() {
    local db="$1"
    python3 - "$db" <<'PY'
import sys
from pathlib import Path
d = Path(sys.argv[1])
BWA = {".fasta." + e for e in ("amb", "ann", "bwt", "pac", "sa")}
fa = tx = None
for p in sorted(d.glob("*.fasta")):
    if any(p.name.endswith(s) for s in BWA):
        continue
    if p.name.startswith("partition"):
        continue
    fa = p.name
    break
for p in sorted(d.glob("*_taxonomy.txt")):
    if p.name.startswith("partition"):
        continue
    tx = p.name
    break
print(f"{fa or 'NONE'} {tx or 'NONE'}")
PY
}

# ---------------------------------------------------------------------------
# marker.* aliases
#
# The harness hardcodes primer_name="marker" for ablation, so it needs
# marker.fasta / marker_taxonomy.txt / marker.fasta.{amb,ann,bwt,pac,sa}
# regardless of what -p was used.
#
# MUST be the last mutation of a database. repair-orphaned-partitions.sh
# rewrites the DB FASTA with `mv` and re-runs `bwa index`, so any alias made
# before a repair points at replaced content. dvc add also materialises
# symlinks into real files, at which point a stale alias is a silent stale
# COPY -- which is exactly what happened to 16Smamm this session.
# ---------------------------------------------------------------------------
alias_marker_files() {
    local db="$1" primer="$2" e
    [ -f "$db/${primer}.fasta" ] || { echo "alias: missing $db/${primer}.fasta" >&2; return 1; }
    rm -f "$db/marker.fasta" "$db/marker_taxonomy.txt"
    ln -sf "${primer}.fasta"          "$db/marker.fasta"
    ln -sf "${primer}_taxonomy.txt"   "$db/marker_taxonomy.txt"
    for e in "${BWA_EXTS[@]}"; do
        rm -f "$db/marker.fasta.$e"
        ln -sf "${primer}.fasta.$e"   "$db/marker.fasta.$e"
    done
}

# ---------------------------------------------------------------------------
# rCRUX marker resolution
#
# Sets: RC_MARKER_DIR, RC_TIER_DIR, RC_FASTA, RC_TAX, RC_TAG
# Only markers with a dedup tier AND a version tag are buildable; the older
# markers use the legacy raw/unfiltered/filtered layout with primer-prefixed
# filenames and have no tag.
# ---------------------------------------------------------------------------
resolve_marker() {
    local m="$1"
    case "$m" in
        */*|"") echo "marker must be a bare directory name, not a path: '$m'" >&2; return 1 ;;
    esac
    RC_MARKER_DIR="$RCRUX_REPO/databases/$m"
    [ -d "$RC_MARKER_DIR" ] || { echo "no such marker directory: $RC_MARKER_DIR" >&2; return 1; }
    if [ ! -f "$RC_MARKER_DIR/dedup.dvc" ]; then
        {
            echo "marker '$m' has no dedup.dvc (legacy raw/unfiltered/filtered layout)."
            echo "Buildable markers:"
            local d
            for d in "$RCRUX_REPO"/databases/*/dedup.dvc; do
                [ -f "$d" ] && echo "  $(basename "$(dirname "$d")")"
            done
        } >&2
        return 1
    fi
    RC_TIER_DIR="$RC_MARKER_DIR/dedup"
    RC_FASTA="$RC_TIER_DIR/lca.fasta"
    RC_TAX="$RC_TIER_DIR/lca_taxonomy.txt"
    RC_TAG=$(git -C "$RCRUX_REPO" tag -l "${m}_v*" --sort=-v:refname 2>/dev/null | head -1)
    [ -n "$RC_TAG" ] || { echo "marker '$m' has no ${m}_v* tag in $RCRUX_REPO" >&2; return 1; }
    return 0
}

# ---------------------------------------------------------------------------
# Input preflight
#
# rCRUX emits accession/taxid/lineage. tronko-build takes EVERYTHING after the
# first tab as the lineage, so a 3-column file shifts all 7 ranks and still
# exits 0. Guard the derived file, not just the source.
#
# The ':' check guards build-tronko-db.sh:447, where gawk rebuilds $0 with
# OFS=" " if the gsub on $1 ever fires -- silently turning a tab-separated
# record into a space-separated one and losing the lineage.
# ---------------------------------------------------------------------------
check_taxonomy_columns() {
    local f="$1" want="$2"
    local n
    n=$(awk -F'\t' '{print NF; exit}' "$f")
    [ "$n" = "$want" ] || { echo "$f: expected $want columns, found $n" >&2; return 1; }
    local bad
    bad=$(awk -F'\t' -v w="$want" 'NF!=w{c++} END{print c+0}' "$f")
    [ "$bad" = "0" ] || { echo "$f: $bad rows do not have $want columns" >&2; return 1; }
}

check_taxonomy_ranks() {
    local f="$1" col="$2"
    local bad
    bad=$(awk -F'\t' -v c="$col" -v r="$TAX_RANKS" \
          '{n=split($c,a,";"); if(n!=r) bad++} END{print bad+0}' "$f")
    [ "$bad" = "0" ] || { echo "$f: $bad rows do not have $TAX_RANKS ranks" >&2; return 1; }
}

check_accession_charset() {
    local f="$1"
    local bad
    bad=$(cut -f1 "$f" | grep -cE '[:[:space:]]' || true)
    [ "$bad" = "0" ] || { echo "$f: $bad accessions contain ':' or whitespace" >&2; return 1; }
}

# Accession sets must match in BOTH directions. Counting alone hides the case
# where one side drops N and duplicates N -- which is exactly how an 11%
# coverage hole in vert12S went unnoticed behind matching record counts.
check_accession_sets() {
    local fasta="$1" tax="$2" tmp
    tmp=$(mktemp -d); trap 'rm -rf "$tmp"' RETURN
    grep '^>' "$fasta" | sed 's/^>//' | awk '{print $1}' | sort -u > "$tmp/f"
    cut -f1 "$tax" | sort -u > "$tmp/t"
    local only_f only_t
    only_f=$(comm -23 "$tmp/f" "$tmp/t" | wc -l)
    only_t=$(comm -13 "$tmp/f" "$tmp/t" | wc -l)
    if [ "$only_f" != "0" ] || [ "$only_t" != "0" ]; then
        echo "accession sets differ: $only_f only in FASTA, $only_t only in taxonomy" >&2
        return 1
    fi
}

# ---------------------------------------------------------------------------
# 2-column taxonomy derivation, with a source-hash sidecar so a stale file from
# a previous rCRUX version is never silently reused.
# ---------------------------------------------------------------------------
make_2col_taxonomy() {
    local src="$1" out="$2"
    local src_md5 sidecar="${2}.src.md5"
    src_md5=$(md5sum "$src" | cut -d' ' -f1)
    if [ -f "$out" ] && [ -f "$sidecar" ] && [ "$(cat "$sidecar")" = "$src_md5" ]; then
        log "  2-column taxonomy up to date ($(wc -l < "$out") rows)"
    else
        cut -f1,3 "$src" > "$out"
        printf '%s\n' "$src_md5" > "$sidecar"
        log "  derived 2-column taxonomy: $(wc -l < "$out") rows"
    fi
    check_taxonomy_columns "$out" 2 || return 1
    check_taxonomy_ranks   "$out" 2 || return 1
    check_accession_charset "$out"  || return 1
    local a b
    a=$(wc -l < "$src"); b=$(wc -l < "$out")
    [ "$a" = "$b" ] || { echo "2-col row count $b != source $a" >&2; return 1; }
}

# ---------------------------------------------------------------------------
# Toolchain identity. 16Smamm shipped with orphaned partitions purely because
# its tronko-build binary predated the createNewRoots fix, so record what built
# a database and refuse to build with a binary older than its sources.
# ---------------------------------------------------------------------------
binary_fingerprint() {
    local b
    for b in "$TRONKO_REPO_DIR/tronko-build/tronko-build" \
             "$TRONKO_REPO_DIR/tronko-convert/tronko-convert" \
             "$TRONKO_REPO_DIR/tronko-assign/tronko-assign" \
             "$TRONKO_REPO_DIR/bin/ancestralclust"; do
        [ -x "$b" ] && printf '%s %s %s\n' \
            "$(basename "$b")" "$(sha256sum "$b" | cut -c1-16)" \
            "$(date -Is -r "$b")"
    done
}

assert_binary_fresh() {
    local bin="$TRONKO_REPO_DIR/tronko-build/tronko-build"
    [ -x "$bin" ] || { echo "tronko-build not built: $bin" >&2; return 1; }
    # Headers matter as much as .c files here: MAX_NUMBEROFROOTS and the tree
    # struct layout live in global.h, so a header-only change still requires a
    # rebuild before the binary reflects it.
    local newer
    newer=$(find "$TRONKO_REPO_DIR/tronko-build" -maxdepth 1 \( -name '*.c' -o -name '*.h' \) \
            -newer "$bin" -printf '%f ' 2>/dev/null)
    if [ -n "$newer" ]; then
        echo "tronko-build is older than its sources ($newer) -- run 'make' in tronko-build/ first" >&2
        return 1
    fi
}
