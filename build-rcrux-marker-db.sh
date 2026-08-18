#!/bin/bash
# Build a tronko LCA database from an rCRUX marker, given only the marker name.
#
# Replaces the per-marker scripts that were hand-copied and edited in eight
# places each -- a pattern that already produced one defect (the MiFish copy
# wrote to a different databases/ subdirectory than its siblings, and the output
# had to be relocated mid-build, invalidating every path in tree_list.txt).
#
# LCA only, by design. A tronko leaf holds exactly one lineage
# (taxonomyArr[i][j][k]; reference_tree.txt writes exactly numspec taxonomy
# lines per tree), and rCRUX accessions are content hashes, so a sequence shared
# across taxa cannot be split into multiple leaves. tronko-assign already
# computes an LCA over the nodes it votes for, so lca_taxonomy.txt encodes
# exactly the ambiguity tronko is built to propagate; species_taxonomy.txt would
# force an arbitrary single species onto sequences that genuinely span several.
# To relabel a finished database onto a different taxonomy, use
# remap-tronko-taxonomy.py -- do not add a variant mode here.
#
# Scope is build + verify. dvc add / dvc push / git commit stay deliberate.
#
# Usage: build-rcrux-marker-db.sh <RCRUX_MARKER_DIR> [options]
#
#   --sp N              SP threshold (default 0.10, legacy scale)
#   --threads N         FAMSA/tree threads (default 64)
#   --jobs N            parallel cluster jobs in step 2 (default 8)
#   --bin-size N        AncestralClust bin size (default 10000)
#   --descendants N     AncestralClust -p (default 75)
#   --primer-name NAME  output FASTA stem (default: the marker directory name)
#   --out-root DIR      databases/ root (default: this repo)
#   --force             discard any existing build AND cache, then rebuild
#   --allow-rcrux-drift proceed when the working tree differs from the tag
#   --no-verify         skip the post-build verification (not recommended)

set -uo pipefail

LIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/tronko-db.sh"
[ -f "$LIB" ] || { echo "ERROR: missing $LIB" >&2; exit 1; }
# shellcheck source=lib/tronko-db.sh
source "$LIB"

MARKER=""
SP=0.10
THREADS=64
JOBS=8
BIN_SIZE=10000
DESCENDANTS=75
PRIMER=""
OUT_ROOT="$TRONKO_REPO_DIR"
FORCE=0
ALLOW_DRIFT=0
DO_VERIFY=1

while [ $# -gt 0 ]; do
    case "$1" in
        --sp)                SP="$2"; shift 2 ;;
        --threads)           THREADS="$2"; shift 2 ;;
        --jobs)              JOBS="$2"; shift 2 ;;
        --bin-size)          BIN_SIZE="$2"; shift 2 ;;
        --descendants)       DESCENDANTS="$2"; shift 2 ;;
        --primer-name)       PRIMER="$2"; shift 2 ;;
        --out-root)          OUT_ROOT="$2"; shift 2 ;;
        --force)             FORCE=1; shift ;;
        --allow-rcrux-drift) ALLOW_DRIFT=1; shift ;;
        --no-verify)         DO_VERIFY=0; shift ;;
        -h|--help)           sed -n '2,32p' "$0"; exit 0 ;;
        -*)                  die "unknown option: $1" ;;
        *)                   [ -z "$MARKER" ] && MARKER="$1" || die "unexpected argument: $1"; shift ;;
    esac
done
[ -n "$MARKER" ] || die "usage: $0 <RCRUX_MARKER_DIR> [options]  (see --help)"

export PATH="$TRONKO_REPO_DIR/bin:$TRONKO_REPO_DIR/tronko-build:$TRONKO_REPO_DIR/tronko-convert:$TRONKO_REPO_DIR/tronko-assign:$PATH"

# ===========================================================================
# A. Resolve, validate, lock
# ===========================================================================
echo "=== A. resolve and lock ==="
resolve_marker "$MARKER" || exit 1
[ -n "$PRIMER" ] || PRIMER="$MARKER"

# The primer name decides the DB FASTA filename, and the harness picks the
# reference with a sorted glob that does NOT skip input.fasta. A primer sorting
# at or after "input" therefore silently selects the raw rCRUX input instead of
# the reference set. Refuse rather than build something that benchmarks wrong.
if [[ ! "$PRIMER" < "input" ]]; then
    die "primer name '$PRIMER' sorts at/after 'input', so the benchmark harness would
       load input.fasta instead of ${PRIMER}.fasta. Choose a name sorting before 'input'."
fi
case "$PRIMER" in
    marker|partition*) die "primer name '$PRIMER' collides with reserved filenames" ;;
esac

CONFIG="$OUT_ROOT/databases/$MARKER/lca/ac/default"
OUT="$CONFIG/sp$SP"
mkdir -p "$CONFIG" || die "cannot create $CONFIG"
CONFIG=$(cd "$CONFIG" && pwd)
OUT="$CONFIG/sp$SP"

# Keep the config .gitignore complete. Only the sp*/ output is DVC-tracked; the
# cache, the derived taxonomy, its sidecar and the lock are all build-local and
# must never be committed. (They are invisible here only because of a local
# .git/info/exclude, which does not travel.)
GI="$CONFIG/.gitignore"
touch "$GI"
for pat in "/sp$SP" /cache /lca_taxonomy.2col.txt /lca_taxonomy.2col.txt.src.md5 /.build.lock; do
    grep -qxF "$pat" "$GI" 2>/dev/null || printf '%s\n' "$pat" >> "$GI"
done

# The AncestralClust cache is shared by every sp* under this config directory,
# so two concurrent builds of the same marker would corrupt each other's
# clusters. Hold the lock for the whole run.
exec 9>"$CONFIG/.build.lock"
flock -n 9 || die "another build holds the lock on $CONFIG"

echo "  marker : $MARKER"
echo "  tag    : $RC_TAG"
echo "  primer : $PRIMER"
echo "  output : $OUT"

assert_binary_fresh || exit 1
echo "  toolchain:"
binary_fingerprint | sed 's/^/    /'
TRONKO_HEAD=$(git -C "$TRONKO_REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)
TRONKO_DIRTY=$(git -C "$TRONKO_REPO_DIR" status --porcelain 2>/dev/null | grep -qv '^??' && echo dirty || echo clean)
echo "    tronko repo: $TRONKO_HEAD ($TRONKO_DIRTY)"

# ===========================================================================
# B. Acquire and verify inputs
# ===========================================================================
echo "=== B. acquire and verify rCRUX inputs ==="
# Unconditional: 4 of the 7 buildable markers are pointer-only, and an EMPTY
# dedup/ directory exists for them, so "is lca.fasta missing" is too weak a test.
( cd "$RCRUX_REPO" && dvc pull "databases/$MARKER/dedup.dvc" ) >/dev/null 2>&1 \
    || warn "dvc pull returned non-zero (continuing if files are present)"
[ -s "$RC_FASTA" ] || die "no $RC_FASTA after dvc pull"
[ -s "$RC_TAX" ]   || die "no $RC_TAX after dvc pull"

# The tag is what provenance will claim; make sure the working tree still
# matches it, rather than asserting it in a comment as the old scripts did.
TAG_MD5=$(git -C "$RCRUX_REPO" show "$RC_TAG:databases/$MARKER/dedup.dvc" 2>/dev/null | awk '/md5:/{print $3; exit}')
WT_MD5=$(awk '/md5:/{print $3; exit}' "$RC_MARKER_DIR/dedup.dvc")
if [ -n "$TAG_MD5" ] && [ "$TAG_MD5" != "$WT_MD5" ]; then
    if [ "$ALLOW_DRIFT" = "1" ]; then
        warn "working tree dedup.dvc ($WT_MD5) != tag $RC_TAG ($TAG_MD5) -- proceeding"
    else
        die "working tree dedup.dvc ($WT_MD5) differs from tag $RC_TAG ($TAG_MD5).
       Provenance would claim a tag this build did not use. Re-run with --allow-rcrux-drift to override."
    fi
else
    echo "  dedup.dvc matches tag $RC_TAG"
fi

echo "  preflighting $(basename "$RC_TAX") ..."
check_taxonomy_columns  "$RC_TAX" 3 || die "source taxonomy is not 3-column"
check_taxonomy_ranks    "$RC_TAX" 3 || die "source taxonomy rank count wrong"
check_accession_charset "$RC_TAX"    || die "source taxonomy has unusable accessions"
DUPACC=$(cut -f1 "$RC_TAX" | sort | uniq -d | wc -l)
[ "$DUPACC" = "0" ] || die "source taxonomy has $DUPACC duplicated accessions (species tier?)"
check_accession_sets "$RC_FASTA" "$RC_TAX" || die "FASTA and taxonomy cover different accessions"
NSEQ=$(grep -c '^>' "$RC_FASTA")
echo "  $NSEQ sequences, accession sets agree, 7 ranks on every row"

# ===========================================================================
# C. Derive the 2-column taxonomy
# ===========================================================================
echo "=== C. derive 2-column taxonomy ==="
TAX2="$CONFIG/lca_taxonomy.2col.txt"
make_2col_taxonomy "$RC_TAX" "$TAX2" || die "2-column taxonomy failed its checks"

# ===========================================================================
# D. Resume decision
# ===========================================================================
echo "=== D. resume decision ==="
CACHE="$CONFIG/cache"
INPUT_SIG=$( { md5sum "$RC_FASTA" | cut -d' ' -f1; md5sum "$TAX2" | cut -d' ' -f1; } | md5sum | cut -d' ' -f1)
SIGFILE="$CACHE/.inputs.md5"

if [ "$FORCE" = "1" ]; then
    echo "  --force: discarding existing build and cache"
    rm -rf "$OUT" "$CACHE"
elif [ -d "$OUT" ]; then
    if db_status_has "$OUT" DONE; then
        echo "  existing completed build found"
        if [ "$DO_VERIFY" != "1" ]; then
            echo "  --no-verify: leaving it untouched (use --force to rebuild)"
            exit 0
        fi
        if "$TRONKO_REPO_DIR/verify-tronko-db.sh" "$OUT"; then
            echo "  already built and verified -- nothing to do (use --force to rebuild)"
            exit 0
        fi
        die "existing build is present but did not verify. Inspect it, or re-run with --force."
    fi
    if pgrep -f "tronko-build .*$OUT" >/dev/null 2>&1 || pgrep -f "ancestralclust .*$CACHE" >/dev/null 2>&1; then
        die "a build appears to be running against $OUT"
    fi
    # A partial build cannot be resumed in place: build-tronko-db.sh step 5
    # concatenates only what final_partitions.txt lists, so stale partition*
    # files from a previous run survive and are picked up later by the harness.
    echo "  incomplete build found -- discarding \$OUT, keeping cache"
    rm -rf "$OUT"
fi

# The cache short-circuits step 1 on a .step1_done marker that carries no
# reference to the input, so a cache built from different inputs would silently
# be reused under this build's provenance.
if [ -d "$CACHE" ] && [ -f "$SIGFILE" ] && [ "$(cat "$SIGFILE")" != "$INPUT_SIG" ]; then
    echo "  cache was built from different inputs -- discarding it"
    rm -rf "$CACHE"
fi
mkdir -p "$OUT" "$CACHE"
printf '%s\n' "$INPUT_SIG" > "$SIGFILE"
rm -rf "$CACHE/merged"   # renumbered each run; stale higher numbers confuse nothing but are noise

# ===========================================================================
# E. Build
# ===========================================================================
echo "=== E. build ==="
db_status "$OUT" START
BUILD_T0=$(date +%s)

bash "$TRONKO_REPO_DIR/build-tronko-db.sh" \
    -f "$RC_FASTA" \
    -t "$TAX2" \
    -o "$OUT" \
    -p "$PRIMER" \
    -T "$THREADS" \
    -s "$SP" \
    -E \
    -L \
    -B "$BIN_SIZE" \
    -P "$DESCENDANTS" \
    -J "$JOBS" \
    --cache-dir "$CACHE" \
    --rcrux-tag "$RC_TAG" 2>&1 | tee -a "$OUT/build.log"
rc=${PIPESTATUS[0]}
db_status "$OUT" "EXIT $rc"
[ "$rc" -eq 0 ] || die "build-tronko-db.sh failed (rc=$rc); see $OUT/build.log"
BUILD_SECS=$(( $(date +%s) - BUILD_T0 ))

# A failed cluster silently drops its whole block of sequences from the
# reference; build-tronko-db.sh only warns about it.
if grep -qiE 'cluster [0-9]+ (failed|FAILED)' "$OUT/build.log"; then
    die "build.log reports cluster failures -- the reference would be incomplete"
fi

# ===========================================================================
# F. Repair BEFORE anything is named
#
# repair-orphaned-partitions.sh rewrites the DB FASTA and re-runs bwa index, so
# any alias created before this point would reference replaced content.
# ===========================================================================
echo "=== F. orphan check ==="
if ! "$TRONKO_REPO_DIR/check-orphaned-partitions.sh" "$OUT"; then
    echo "  orphans found -- repairing"
    "$TRONKO_REPO_DIR/repair-orphaned-partitions.sh" "$OUT" || die "repair failed"
    db_status "$OUT" REPAIRED
    "$TRONKO_REPO_DIR/check-orphaned-partitions.sh" "$OUT" || die "orphans remain after repair"
fi

# ===========================================================================
# G. Record inputs and alias -- the last mutations
# ===========================================================================
echo "=== G. record inputs and alias ==="
cp "$RC_FASTA" "$OUT/input.fasta"
cp "$TAX2"     "$OUT/input_taxonomy.txt"
db_status "$OUT" INPUTS_RECORDED
alias_marker_files "$OUT" "$PRIMER" || die "aliasing failed"
db_status "$OUT" ALIASED
db_status "$OUT" DONE

# ===========================================================================
# H/I. Verify and report
# ===========================================================================
VERIFY_RC=0
if [ "$DO_VERIFY" = "1" ]; then
    echo "=== H. verify ==="
    "$TRONKO_REPO_DIR/verify-tronko-db.sh" "$OUT" || VERIFY_RC=1
fi

echo "=== I. summary ==="
echo "  marker        : $MARKER  (tag $RC_TAG)"
echo "  output        : $OUT"
echo "  build time    : ${BUILD_SECS}s"
echo "  input seqs    : $NSEQ"
[ -f "$OUT/tree_list.txt" ] && echo "  trees         : $(wc -l < "$OUT/tree_list.txt")"
[ -f "$OUT/${PRIMER}.fasta" ] && echo "  reference seqs: $(grep -c '^>' "$OUT/${PRIMER}.fasta")"
echo
echo "  Not done (deliberately): dvc add / dvc push / git commit."
echo "  The benchmark harness looks under TRONKO_DB_ROOT (default ~/tronko/databases)"
echo "  and its presets use short marker names, so to benchmark this build either:"
echo "    export TRONKO_DB_ROOT=$OUT_ROOT/databases"
echo "  and point the preset's tronko_db_dir at $MARKER/lca/ac/default/sp$SP"

exit $VERIFY_RC
