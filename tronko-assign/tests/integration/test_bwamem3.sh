#!/usr/bin/env bash
# Checks of tronko-assign's BWA-MEM3 run (bwamem3.c) that the golden tests cannot make: the
# command line and environment it gives BWA-MEM3, the shared-memory index's life, the index build,
# and each failure it must stop on. A stand-in bwa-mem3 (a script written here) logs how it was
# called, passes every command to the real BWA-MEM3, and in the failure cases changes what it
# returns:
#   version   reports version 0.13.0                    -> refused, "needs BWA-MEM3 0.14.0"
#   memfail   mem exits 3                                -> "bwa-mem3 mem failed on batch 1"
#   dropread  mem's SAM loses read 5's records           -> "gave no record for read 5"
#   reorder   mem's SAM has read 0's records last        -> "out of input order"
# Input: the first 20 reads of tests/data/assignment/single_4000.fasta on the example reference.
# Environment: TRONKO_ASSIGN_BIN, TRONKO_BWA_MEM3 (see golden_lib.sh).
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=golden_lib.sh
source "$SCRIPT_DIR/golden_lib.sh"
EX="$REPO_ROOT/tronko-build/example_datasets/single_tree"
FASTA="$EX/Charadriiformes.fasta"
REFERENCE="$REPO_ROOT/tests/data/assignment/reference_tree.trkb"
READS="$REPO_ROOT/tests/data/assignment/single_4000.fasta"
for f in "$FASTA" "$REFERENCE" "$READS"; do
	[[ -f $f ]] || { echo "SKIP: missing fixture file $f"; exit 77; }
done
mem3_index "$FASTA" || exit 2
head -40 "$READS" >"$TMP_DIR/reads.fasta"

REAL=$TRONKO_BWA_MEM3
FAKE=$TMP_DIR/bin/bwa-mem3
mkdir -p "$TMP_DIR/bin"
cat >"$FAKE" <<'EOF'
#!/usr/bin/env bash
{ printf 'ARGS'; printf ' %s' "$@"; echo; env | grep -E '^(BWAMEM3_|BWA_MEM3_|BWA3_)' | sed 's/^/ENV /'; } >>"$FAKE_LOG"
mode=${FAKE_MODE:-}
if [[ $1 == version && $mode == version ]]; then echo 0.13.0; exit 0; fi
if [[ $1 != mem ]]; then exec "$REAL_BWA_MEM3" "$@"; fi
[[ $mode == memfail ]] && exit 3
out=
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do [[ ${args[i]} == -o ]] && out=${args[i+1]}; done
"$REAL_BWA_MEM3" "$@" || exit $?
case $mode in
dropread) awk -F'\t' '/^@/ || $1 != "5"' "$out" >"$out.x" && mv "$out.x" "$out" ;;
reorder) { awk -F'\t' '/^@/ || $1 != "0"' "$out"; awk -F'\t' '!/^@/ && $1 == "0"' "$out"; } >"$out.x" && mv "$out.x" "$out" ;;
esac
EOF
chmod +x "$FAKE"
export REAL_BWA_MEM3=$REAL FAKE_LOG=$TMP_DIR/fake.log

pass() { echo "PASS: $1"; }
fail() { echo "FAIL: $1" >&2; FAILS=$((FAILS + 1)); }

# run_case <label> <mode> <extra tronko args...>: sets RC and LOG
run_case() {
	local label=$1 mode=$2 dir
	shift 2
	dir=$(mktemp -d "$TMP_DIR/run.XXXXXX")
	: >"$FAKE_LOG"
	(cd "$dir" && FAKE_MODE=$mode TRONKO_BWA_MEM3=$FAKE BWAMEM3_FORCE_TIER=scalar BWA_MEM3_COHORT_SLICES=6 BWA3_KS_DEDUP=on \
		"$ASSIGN" -r -f "$REFERENCE" -w --Cinterval 10 -s -g "$TMP_DIR/reads.fasta" -o "$dir/out.tsv" "$@") >"$dir/log" 2>&1
	RC=$?
	LOG=$dir/log
	OUT=$dir/out.tsv
}

# 1. A good run: the pinned command line, the thread count, a scrubbed environment, the index
# staged once and dropped, the output equal to the real BWA-MEM3's
run_case good "" -a "$FASTA" -6 --number-of-cores 3
if ((RC == 0)); then pass "good run exits 0"; else fail "good run exit $RC: $(tail -3 "$LOG")"; fi
grep -q -- '^ARGS mem --compat=bwa-mem -K 100000000 --dedup-reads off --dedup off --ks-dedup off --cohort-slices 0 -v 1 -t 3 -o ' "$FAKE_LOG" &&
	pass "mem command line pinned, -t = --number-of-cores" || fail "mem command line: $(grep '^ARGS mem' "$FAKE_LOG")"
grep -q '^ENV ' "$FAKE_LOG" && fail "BWA-MEM3 environment overrides reached the aligner: $(grep -m3 '^ENV ' "$FAKE_LOG")" ||
	pass "BWAMEM3_*, BWA_MEM3_*, BWA3_* removed from the aligner's environment"
[[ $(grep -c '^ARGS shm /' "$FAKE_LOG") == 1 && $(grep -c '^ARGS shm -l' "$FAKE_LOG") == 1 ]] &&
	pass "index staged once, listed once at exit" || fail "shm calls: $(grep '^ARGS shm' "$FAKE_LOG")"
staged=$(grep -m1 '^ARGS shm /' "$FAKE_LOG" | sed 's|.*/||')
if "$REAL" shm -l 2>/dev/null | cut -f1 | grep -qx "$staged"; then fail "$staged still staged after exit"; else pass "staged index dropped at exit"; fi
good_out=$OUT
(cd "$TMP_DIR" && TRONKO_BWA_MEM3=/nonexistent "$ASSIGN" -r -f "$REFERENCE" -w --Cinterval 10 -s -g "$TMP_DIR/reads.fasta" \
	-a "$FASTA" -6 --number-of-cores 3 --bwa-mem3 "$REAL" -o "$TMP_DIR/direct.tsv") >/dev/null 2>&1
cmp -s "$good_out" "$TMP_DIR/direct.tsv" && pass "--bwa-mem3 wins over TRONKO_BWA_MEM3; output equals the stand-in run's" ||
	fail "--bwa-mem3 run differs or failed"

# 2. --no-shm: no staging
run_case noshm "" -a "$FASTA" -6 --no-shm
((RC == 0)) && ! grep -q '^ARGS shm' "$FAKE_LOG" && pass "--no-shm stages nothing" || fail "--no-shm: exit $RC, $(grep '^ARGS shm' "$FAKE_LOG")"

# 3. Index: missing with -6 is refused; without -6 it is built beside the FASTA
mkdir -p "$TMP_DIR/noindex"
cp "$FASTA" "$TMP_DIR/noindex/ref.fasta"
run_case missing-index "" -a "$TMP_DIR/noindex/ref.fasta" -6
((RC != 0)) && grep -q 'BWA-MEM3 index file .*ref.fasta.amb is missing' "$LOG" && pass "-6 without an index refused" ||
	fail "-6 without an index: exit $RC, $(tail -2 "$LOG")"
run_case build-index "" -a "$TMP_DIR/noindex/ref.fasta"
((RC == 0)) && [[ -f $TMP_DIR/noindex/ref.fasta.bwt.2bit.64 ]] && grep -q '^ARGS index ' "$FAKE_LOG" && cmp -s "$OUT" "$good_out" &&
	pass "without -6 the index is built and the output is the same" || fail "index build: exit $RC, $(tail -2 "$LOG")"

# 4. Failures
expect_failure() {
	local mode=$1 pattern=$2
	run_case "$mode" "$mode" -a "$FASTA" -6
	if ((RC != 0)) && grep -q -- "$pattern" "$LOG"; then pass "$mode refused: $pattern"; else fail "$mode: exit $RC, $(tail -2 "$LOG")"; fi
}
expect_failure version "needs BWA-MEM3 0.14.0"
expect_failure memfail "bwa-mem3 mem failed on batch 1"
expect_failure dropread "gave no record for read 5"
expect_failure reorder "out of input order"

((FAILS == 0)) || { echo "bwa-mem3 checks: $FAILS failure(s)" >&2; exit 1; }
