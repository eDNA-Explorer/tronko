# Shared by the golden parity scripts in this directory (sourced, not run).
#
# Each script runs tronko-assign on a fixture with the flags production uses and compares the
# output byte for byte with a golden: the output of production tronko-assign (high-perf,
# 71f6ec3) at --number-of-cores 1, built from the source as committed (MAX_NUM_BWA_MATCHES 10)
# and with production's match cap of 25. The goldens for cap 10 live in <dir>, those for cap 25
# in <dir>-cap25.
#
# Environment:
#   TRONKO_ASSIGN_BIN     binary under test (default tronko-assign/tronko-assign)
#   TRONKO_ASSIGN_CORES   thread counts that must reproduce every golden (default "1")
#   TRONKO_MATCH_CAP      MAX_NUM_BWA_MATCHES the binary was built with (default: read from
#                         tronko-assign/global.h, so a tree patched the way production patches it
#                         picks the cap-25 goldens by itself)
#   TRONKO_TESTS_XFAIL    case names whose mismatch is reported as XFAIL instead of failing
#   TRONKO_TESTS_RECORD=1 write each output as the golden instead of comparing (one thread count
#                         only); how the goldens are (re)made from a given commit
#
# The caller sets SCRIPT_DIR before sourcing. Sets REPO_ROOT, TA, ASSIGN, CORES, MATCH_CAP,
# RECORD, TMP_DIR (removed on exit) and FAILS.

REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
TA="$REPO_ROOT/tronko-assign"
ASSIGN=${TRONKO_ASSIGN_BIN:-"$TA/tronko-assign"}
CORES=${TRONKO_ASSIGN_CORES:-"1"}
MATCH_CAP=${TRONKO_MATCH_CAP:-$(sed -n 's/^#define MAX_NUM_BWA_MATCHES \([0-9][0-9]*\).*/\1/p' "$TA/global.h")}
RECORD=${TRONKO_TESTS_RECORD:-0}
XFAIL=" ${TRONKO_TESTS_XFAIL:-} "
FAILS=0
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/tronko-goldens.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

[[ -x "$ASSIGN" ]] || { echo "Missing executable: $ASSIGN (build tronko-assign first)" >&2; exit 1; }
case $MATCH_CAP in
10 | 25) ;;
*) echo "SKIP: no goldens for MAX_NUM_BWA_MATCHES ${MATCH_CAP:-?} (10 and 25 only)"; exit 77 ;;
esac
if [[ $RECORD == 1 && $CORES != 1 ]]; then
	echo "TRONKO_TESTS_RECORD=1 records at one thread only (TRONKO_ASSIGN_CORES=1)" >&2
	exit 2
fi

# golden_dir <cap-10 golden directory>: the directory for MATCH_CAP
golden_dir() {
	if [[ $MATCH_CAP == 25 ]]; then echo "${1%/}-cap25"; else echo "${1%/}"; fi
}

# check_output <case> <golden> <output> <label>: compare (or record); counts failures in FAILS
check_output() {
	local c=$1 golden=$2 out=$3 label=$4 rows
	if [[ $RECORD == 1 ]]; then
		mkdir -p "$(dirname "$golden")"
		cp "$out" "$golden"
		echo "RECORDED: $label -> $golden"
		return 0
	fi
	if [[ ! -f $golden ]]; then
		echo "FAIL: $label: no golden $golden" >&2
		FAILS=$((FAILS + 1))
		return 1
	fi
	if cmp -s "$golden" "$out"; then
		echo "PASS: $label"
		[[ $XFAIL == *" $c "* ]] && echo "  (XPASS: $c is listed in TRONKO_TESTS_XFAIL)"
		return 0
	fi
	rows=$(diff "$golden" "$out" | grep -c '^<')
	if [[ $XFAIL == *" $c "* ]]; then
		echo "XFAIL: $label: $rows golden rows missing or changed (listed in TRONKO_TESTS_XFAIL)"
		return 0
	fi
	echo "FAIL: $label: $rows golden rows missing or changed" >&2
	diff -u "$golden" "$out" | head -20 >&2 || true
	FAILS=$((FAILS + 1))
	return 1
}

# run_assign <log> <args...>: run tronko-assign from a fresh directory under TMP_DIR (the working
# directory receives --tsv-log files); a non-zero exit is a failure
run_assign() {
	local log=$1 dir
	shift
	dir=$(mktemp -d "$TMP_DIR/run.XXXXXX")
	if ! (cd "$dir" && "$ASSIGN" "$@") >"$log" 2>&1; then
		echo "FAIL: tronko-assign exited non-zero: $*" >&2
		tail -5 "$log" >&2
		return 1
	fi
}

# The mate-rescue fixture's output depends, before the fix of the rescue's score-matrix over-read,
# on the bytes of a freed heap copy of the -a path, so every input of that case is read through
# one fixed directory. The files are copied there (atomically, so concurrent runs are harmless)
# and every argument string is the same on every machine and checkout.
MATE_RESCUE_DIR=/tmp/tronko-assign-tests/mate-rescue
# stage_mate_rescue <file>...: copy the files into MATE_RESCUE_DIR unless an identical copy is there
stage_mate_rescue() {
	local f dst
	mkdir -p "$MATE_RESCUE_DIR" || return 1
	for f in "$@"; do
		dst=$MATE_RESCUE_DIR/$(basename "$f")
		cmp -s "$f" "$dst" && continue
		cp "$f" "$dst.tmp.$$" && mv -f "$dst.tmp.$$" "$dst" || return 1
	done
}
