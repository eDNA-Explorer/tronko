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
#   hang      mem waits; tronko-assign is sent SIGTERM   -> ends by SIGTERM, index dropped, temporary
#                                                           directory removed
#   hangindex index builds, then waits; SIGTERM          -> no index under the final names, the
#                                                           build's temporary files removed
# Input: the first 20 reads of tests/data/assignment/single_4000.fasta on the example reference;
# for -6's index build, the three-tree fixture's single-end reads against their golden.
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
if [[ $1 == index && $mode == hangindex ]]; then "$REAL_BWA_MEM3" "$@" >/dev/null 2>&1; exec sleep 120; fi
if [[ $1 != mem ]]; then exec "$REAL_BWA_MEM3" "$@"; fi
[[ $mode == memfail ]] && exit 3
[[ $mode == hang ]] && exec sleep 120
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

# 3. Index without -6: built beside the FASTA in every run
mkdir -p "$TMP_DIR/noindex"
cp "$FASTA" "$TMP_DIR/noindex/ref.fasta"
run_case build-index "" -a "$TMP_DIR/noindex/ref.fasta"
((RC == 0)) && [[ -f $TMP_DIR/noindex/ref.fasta.bwt.2bit.64 ]] && grep -q '^ARGS index ' "$FAKE_LOG" && cmp -s "$OUT" "$good_out" &&
	pass "without -6 the index is built and the output is the same" || fail "index build: exit $RC, $(tail -2 "$LOG")"

# 3b. -6 with no index: built in the run beside the FASTA (under <fasta>.tmp-<pid>, renamed into
# place), then found by the next run; one missing file makes it build again. Every output is the
# three-tree fixture's golden. A run in a read-only directory stops (not checkable as root).
MT=$TA/tests/data/multitree
MT_GOLDEN=$(golden_dir "$MT/goldens")/expected_mt_single_F.tsv
mem3_index "$MT/multitree.fasta" || exit 2
INRUN=$TMP_DIR/inrun
mkdir -p "$INRUN"
cp "$MT/multitree.fasta" "$INRUN/multitree.fasta"
# run_mt <label> <mode> <fasta>: the three-tree single-end case with -6; sets RC, LOG, OUT
run_mt() {
	local dir
	dir=$(mktemp -d "$TMP_DIR/run.XXXXXX")
	: >"$FAKE_LOG"
	(cd "$dir" && FAKE_MODE=$2 TRONKO_BWA_MEM3=$FAKE "$ASSIGN" -r -f "$MT/reference_tree.trkb" -a "$3" -w -6 --Cinterval 10 \
		--number-of-cores 4 -s -g "$MT/single_F.fasta" -o "$dir/out.tsv") >"$dir/log" 2>&1
	RC=$?
	LOG=$dir/log
	OUT=$dir/out.tsv
}
# index_ok: the four files beside the copy equal the ones mem3_index built, nothing else is there
index_ok() {
	local ext extra
	for ext in amb ann pac bwt.2bit.64; do cmp -s "$MT/multitree.fasta.$ext" "$INRUN/multitree.fasta.$ext" || return 1; done
	extra=$(find "$INRUN" -mindepth 1 ! -name multitree.fasta ! -name 'multitree.fasta.amb' ! -name 'multitree.fasta.ann' \
		! -name 'multitree.fasta.pac' ! -name 'multitree.fasta.bwt.2bit.64')
	[[ -z $extra ]]
}
run_mt inrun-build "" "$INRUN/multitree.fasta"
if ((RC == 0)) && grep -q '^ARGS index .*/multitree.fasta.tmp-[0-9]*$' "$FAKE_LOG" &&
	grep -q 'bwa-mem3: index built beside .*multitree.fasta in [0-9.]* s$' "$LOG" && index_ok && cmp -s "$OUT" "$MT_GOLDEN"; then
	pass "-6 with no index: built in the run under a temporary name, logged with its time, output = golden"
else
	fail "-6 with no index: exit $RC, $(grep '^ARGS index' "$FAKE_LOG"), $(grep 'bwa-mem3: index' "$LOG"), $(tail -2 "$LOG")"
fi
run_mt inrun-found "" "$INRUN/multitree.fasta"
if ((RC == 0)) && ! grep -q '^ARGS index' "$FAKE_LOG" && grep -q 'bwa-mem3: index found beside .*multitree.fasta$' "$LOG" &&
	index_ok && cmp -s "$OUT" "$MT_GOLDEN"; then
	pass "-6 second run: index found (logged), not rebuilt, output = golden"
else
	fail "-6 second run: exit $RC, $(grep '^ARGS index' "$FAKE_LOG"), $(grep 'bwa-mem3: index' "$LOG")"
fi
mv "$INRUN/multitree.fasta.ann" "$TMP_DIR/stale.ann"
run_mt inrun-partial "" "$INRUN/multitree.fasta"
if ((RC == 0)) && grep -q '^ARGS index ' "$FAKE_LOG" && grep -q 'bwa-mem3: index built beside' "$LOG" && index_ok &&
	cmp -s "$OUT" "$MT_GOLDEN"; then
	pass "-6 with one index file missing: rebuilt, output = golden"
else
	fail "-6 with one index file missing: exit $RC, $(grep 'bwa-mem3: index' "$LOG"), $(tail -2 "$LOG")"
fi
if (($(id -u) != 0)); then
	ro=$TMP_DIR/readonly
	mkdir -p "$ro"
	cp "$MT/multitree.fasta" "$ro/multitree.fasta"
	chmod a-w "$ro"
	run_mt inrun-readonly "" "$ro/multitree.fasta"
	chmod u+w "$ro"
	((RC != 0)) && grep -q 'index of .*multitree.fasta is missing and its directory is not writable' "$LOG" &&
		! grep -q '^ARGS index' "$FAKE_LOG" && pass "-6 with no index in a read-only directory stops with a message" ||
		fail "-6 read-only directory: exit $RC, $(tail -2 "$LOG")"
else
	echo "SKIP: read-only directory check (running as root)"
fi

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

# 5. Stopped by SIGTERM while the aligner runs: the process ends by the signal, its staged index is
# dropped and its temporary directory removed (bwamem3.c, stop_waiter)
dir=$(mktemp -d "$TMP_DIR/run.XXXXXX")
mkdir -p "$dir/tmp"
: >"$FAKE_LOG"
(cd "$dir" && FAKE_MODE=hang TRONKO_BWA_MEM3=$FAKE TMPDIR=$dir/tmp exec "$ASSIGN" -r -f "$REFERENCE" -w --Cinterval 10 \
	-s -g "$TMP_DIR/reads.fasta" -a "$FASTA" -6 -o "$dir/out.tsv") >"$dir/log" 2>&1 &
pid=$!
for ((i = 0; i < 300; i++)); do grep -q '^ARGS mem' "$FAKE_LOG" && break; sleep 0.1; done
staged=$(grep -m1 '^ARGS shm /' "$FAKE_LOG" | sed 's|.*/||')
kill -TERM "$pid"
wait "$pid"
RC=$?
left=$(find "$dir/tmp" -mindepth 1 | head -3)
if ((RC == 143)) && [[ -n $staged && -z $left ]] && ! "$REAL" shm -l 2>/dev/null | cut -f1 | grep -qx "$staged"; then
	pass "SIGTERM: ends by the signal, staged index dropped, temporary directory removed"
else
	fail "SIGTERM: exit $RC, staged '$staged', left: ${left:-none}, $(tail -2 "$dir/log")"
fi

# 6. Stopped by SIGTERM during -6's index build: no file under the final names (so the next run
# builds again), the build's temporary files removed, then a run builds and gives the golden
STOP=$TMP_DIR/stopbuild
mkdir -p "$STOP/tmp"
cp "$MT/multitree.fasta" "$STOP/multitree.fasta"
: >"$FAKE_LOG"
(cd "$STOP/tmp" && FAKE_MODE=hangindex TRONKO_BWA_MEM3=$FAKE TMPDIR=$STOP/tmp exec "$ASSIGN" -r -f "$MT/reference_tree.trkb" \
	-a "$STOP/multitree.fasta" -w -6 --Cinterval 10 -s -g "$MT/single_F.fasta" -o "$STOP/out.tsv") >"$TMP_DIR/stopbuild.log" 2>&1 &
pid=$!
for ((i = 0; i < 300; i++)); do [[ -f $STOP/multitree.fasta.tmp-$pid.bwt.2bit.64 ]] && break; sleep 0.1; done
built=$([[ -f $STOP/multitree.fasta.tmp-$pid.bwt.2bit.64 ]] && echo yes || echo no)
kill -TERM "$pid"
wait "$pid"
RC=$?
left=$(find "$STOP" -mindepth 1 ! -name multitree.fasta ! -name tmp | head -5)
if ((RC == 143)) && [[ $built == yes && -z $left ]]; then
	pass "SIGTERM during -6's build: ends by the signal, no index under the final names, temporary files removed"
else
	fail "SIGTERM during the build: exit $RC, build files seen: $built, left: ${left:-none}, $(tail -2 "$TMP_DIR/stopbuild.log")"
fi
run_mt stopbuild-after "" "$STOP/multitree.fasta"
((RC == 0)) && grep -q 'bwa-mem3: index built beside' "$LOG" && cmp -s "$OUT" "$MT_GOLDEN" &&
	pass "after a stopped build, the next -6 run builds the index and gives the golden" ||
	fail "after a stopped build: exit $RC, $(tail -2 "$LOG")"

((FAILS == 0)) || { echo "bwa-mem3 checks: $FAILS failure(s)" >&2; exit 1; }
