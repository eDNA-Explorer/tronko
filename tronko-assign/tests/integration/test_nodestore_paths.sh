#!/usr/bin/env bash
# The node-store scoring path (nodestore.c) on the repository fixtures, against the goldens of
# test_assignment_production_parity.sh (made at --number-of-cores 1; see that script).
#
# With the production flags (-w -6 --Cinterval 10) tronko-assign relays each tree's posteriors out
# node-innermost at load and scores reads in blocks against node tiles. Every variant below must
# give the golden byte for byte:
#   - the default path, in the three read modes;
#   - TRONKO_NODESTORE_CHECK=1: before the per-node arrays are released, every stored value is
#     compared with them and every leaf string with the legacy reader (the run exits 3 and prints
#     result=FAIL on any difference);
#   - blocks of 1 and 7 reads, and a 1 MiB block budget that closes blocks before they are full;
#   - every scoring kernel the CPU supports (TRONKO_NS_KERNEL=plain|avx2|avx512);
#   - TRONKO_NODESTORE=0, the legacy loop in the same binary;
#   - reads without candidates at block boundaries, against the legacy loop of the same binary.
#
# Environment: TRONKO_ASSIGN_BIN (default tronko-assign/tronko-assign), TRONKO_BWA_MEM3 (mem3_lib.sh), TRONKO_ASSIGN_CORES
# (default 4), TRONKO_FIXTURE_DIR (default tests/data/assignment, the fixtures of
# test_assignment_production_parity.sh). Exit 77 if those fixtures are absent.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../../.." && pwd)
ASSIGN=${TRONKO_ASSIGN_BIN:-"$REPO_ROOT/tronko-assign/tronko-assign"}
CORES=${TRONKO_ASSIGN_CORES:-4}
DATA=${TRONKO_FIXTURE_DIR:-"$REPO_ROOT/tests/data/assignment"}
FASTA="$REPO_ROOT/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/tronko-nodestore.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

[[ -x "$ASSIGN" ]] || { echo "Missing executable: $ASSIGN" >&2; exit 1; }
for path in "$DATA/reference_tree.trkb" "$DATA/paired_2000_1.fasta" "$DATA/paired_2000_2.fasta" "$DATA/single_4000.fasta"; do
	[[ -f "$path" ]] || { echo "SKIP: missing fixture $path (set TRONKO_FIXTURE_DIR)"; exit 77; }
done
# shellcheck source=mem3_lib.sh
source "$SCRIPT_DIR/mem3_lib.sh"
mem3_index "$FASTA" || exit 1

COMMON=(-r -f "$DATA/reference_tree.trkb" -a "$FASTA" -w -6 --Cinterval 10 --number-of-cores "$CORES")
PAIRED=(-p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta")
UNPAIRED_R=(-s -v -g "$DATA/paired_2000_2.fasta")
SINGLE=(-s -g "$DATA/single_4000.fasta")
FAILS=0

# check <label> <golden> <env assignments...> -- <mode args...>
check() {
	local label=$1 golden=$2 envs=() out log
	shift 2
	while [[ $1 != -- ]]; do envs+=("$1"); shift; done
	shift
	out="$TMP_DIR/$label.tsv" log="$TMP_DIR/$label.log"
	if ! env "${envs[@]}" "$ASSIGN" "${COMMON[@]}" "$@" -o "$out" --tsv-log "$TMP_DIR/$label.tsvlog" >"$log" 2>&1; then
		echo "FAIL: $label: tronko-assign exited non-zero" >&2
		tail -5 "$log" >&2
		FAILS=$((FAILS + 1))
		return
	fi
	if ! cmp -s "$golden" "$out"; then
		echo "FAIL: $label: output differs from $golden" >&2
		diff "$golden" "$out" | head -10 >&2 || true
		FAILS=$((FAILS + 1))
		return
	fi
	if [[ " ${envs[*]} " != *" TRONKO_NODESTORE=0 "* ]] && ! grep -q 'RELAYOUT_DONE.*kernel=' "$TMP_DIR/$label.tsvlog"; then
		echo "FAIL: $label: the node-store path was not taken" >&2
		grep RELAYOUT "$TMP_DIR/$label.tsvlog" >&2 || true
		FAILS=$((FAILS + 1))
		return
	fi
	echo "PASS: $label"
}

check paired "$DATA/expected_paired_nw_trkb.tsv" -- "${PAIRED[@]}"
check unpaired_r "$DATA/expected_unpaired_r_nw_trkb.tsv" -- "${UNPAIRED_R[@]}"
check single "$DATA/expected_single_nw_trkb.tsv" -- "${SINGLE[@]}"

check selfcheck "$DATA/expected_paired_nw_trkb.tsv" TRONKO_NODESTORE_CHECK=1 -- "${PAIRED[@]}"
if ! grep -q '^nodestore check: .*leaf_mismatch=0 result=PASS$' "$TMP_DIR/selfcheck.log"; then
	echo "FAIL: selfcheck: no passing 'nodestore check' line" >&2
	grep 'nodestore check' "$TMP_DIR/selfcheck.log" >&2 || true
	FAILS=$((FAILS + 1))
else
	grep '^nodestore check:' "$TMP_DIR/selfcheck.log"
fi

check block1 "$DATA/expected_unpaired_r_nw_trkb.tsv" TRONKO_NS_BLOCK_READS=1 -- "${UNPAIRED_R[@]}"
check block7 "$DATA/expected_paired_nw_trkb.tsv" TRONKO_NS_BLOCK_READS=7 -- "${PAIRED[@]}"
check budget1mb "$DATA/expected_paired_nw_trkb.tsv" TRONKO_NS_BLOCK_MB=1 -- "${PAIRED[@]}"
if ! grep -q 'NODESTORE_STATS.*budget_closes=[1-9]' "$TMP_DIR/budget1mb.tsvlog"; then
	echo "FAIL: budget1mb: no block was closed by the budget" >&2
	FAILS=$((FAILS + 1))
fi

for kernel in plain avx2 avx512; do
	probe="$TMP_DIR/probe-$kernel.log"
	if TRONKO_NS_KERNEL=$kernel "$ASSIGN" "${COMMON[@]}" "${UNPAIRED_R[@]}" -o "$TMP_DIR/probe.tsv" >"$probe" 2>&1; then
		check "kernel-$kernel" "$DATA/expected_paired_nw_trkb.tsv" TRONKO_NS_KERNEL=$kernel -- "${PAIRED[@]}"
	elif grep -q 'not supported by this CPU' "$probe"; then
		echo "SKIP: kernel $kernel not supported by this CPU"
	else
		echo "FAIL: kernel $kernel probe run failed" >&2
		FAILS=$((FAILS + 1))
	fi
done

check legacy "$DATA/expected_paired_nw_trkb.tsv" TRONKO_NODESTORE=0 -- "${PAIRED[@]}"
if grep -q 'RELAYOUT' "$TMP_DIR/legacy.tsvlog"; then
	echo "FAIL: legacy: TRONKO_NODESTORE=0 still built the node store" >&2
	FAILS=$((FAILS + 1))
fi

# Reads without candidates at block ends and starts. Phase C of a block must use each read's own
# saved candidates: a read with none prints "unassigned" and clears no votes. Twelve pairs of 150 Ns
# are inserted at 0, 1, 6, 7, 8, 63, 64, 65, 127, 128, 1000 and 1999 (one of them still gets a hit
# whose alignment is empty, which takes getscore_Arr's positions[0] == -1 exit). No golden exists
# for these inputs, so the node-store path is compared with the legacy loop of the same binary.
NS=$(printf 'N%.0s' $(seq 150))
for f in paired_2000_1 paired_2000_2; do
	awk -v nn="$NS" 'BEGIN { n = split("0 1 6 7 8 63 64 65 127 128 1000 1999", a, " "); for (i = 1; i <= n; i++) ins[a[i]] = 1 }
		/^>/ { r++; if ((r - 1) in ins) { print ">nohit_" (r - 1); print nn } } { print }' "$DATA/$f.fasta" >"$TMP_DIR/nohit_$f.fasta"
done
NOHIT=(-p -z -1 "$TMP_DIR/nohit_paired_2000_1.fasta" -2 "$TMP_DIR/nohit_paired_2000_2.fasta")
if TRONKO_NODESTORE=0 "$ASSIGN" "${COMMON[@]}" "${NOHIT[@]}" -o "$TMP_DIR/nohit-legacy.tsv" >"$TMP_DIR/nohit-legacy.log" 2>&1; then
	for reads in 64 2 3 7; do
		check "nohit-block$reads" "$TMP_DIR/nohit-legacy.tsv" TRONKO_NS_BLOCK_READS=$reads -- "${NOHIT[@]}"
	done
	if [[ $(grep -c $'^nohit_[0-9]*\tunassigned$' "$TMP_DIR/nohit-legacy.tsv") -lt 6 ]]; then
		echo "FAIL: nohit: fewer than 6 of the inserted reads are unassigned; the test no longer covers reads without candidates" >&2
		FAILS=$((FAILS + 1))
	fi
else
	echo "FAIL: nohit: the legacy run failed" >&2
	FAILS=$((FAILS + 1))
fi

if ((FAILS)); then
	echo "test_nodestore_paths: $FAILS failure(s)" >&2
	exit 1
fi
echo "test_nodestore_paths: all variants byte-identical to the goldens"
