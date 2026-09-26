#!/bin/bash
# record_ksw_calls.sh <output.kswd.gz>: record tests/data/ksw_fixture_calls.kswd.gz again.
#
# The file holds ksw_extend2 calls with the outputs of BWA's original scalar code. This script
# takes the commit before the vectorised ksw_extend2 (the parent of the commit that added
# ksw_extend2_scalar; KSW_RECORD_COMMIT overrides it), applies ksw_record.patch (next to this
# script) to a copy of it, builds tronko-assign, runs the three fixture cases of
# tests/integration/test_assignment_production_parity.sh at one thread (paired, single, unpaired
# reverse, in that order) and writes every hundredth ksw_extend2 call of each run, inputs and
# outputs, into one gzip file. The decompressed output equals the committed file byte for byte:
#   cmp <(gzip -dc <output.kswd.gz>) <(gzip -dc tronko-assign/tests/data/ksw_fixture_calls.kswd.gz)
set -euo pipefail
out=$(realpath -m "${1:?usage: record_ksw_calls.sh <output.kswd.gz>}")
tools=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(git -C "$tools" rev-parse --show-toplevel)
ksw=tronko-assign/bwa_source_files/ksw.c
commit=${KSW_RECORD_COMMIT:-$(git -C "$root" log --format=%H -S ksw_extend2_scalar -- "$ksw" | tail -1)^}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
git -C "$root" archive "$commit" | tar x -C "$work"
(cd "$work" && patch -p1 --quiet < "$tools/ksw_record.patch")
make -C "$work/tronko-assign" > "$work/build.log" 2>&1 || { tail -20 "$work/build.log" >&2; exit 1; }
bin=$work/tronko-assign/tronko-assign
fx=$root/tests/data/assignment
common=(-r -f "$fx/reference_tree.trkb" -a "$root/tronko-build/example_datasets/single_tree/Charadriiformes.fasta" -w -6 --Cinterval 10 --number-of-cores 1)
record() { # record <case> <mode arguments...>
	local c=$1; shift
	KSW_DUMP=$work/$c.kswd "$bin" "${common[@]}" "$@" -o "$work/$c.tsv" > "$work/$c.log" 2>&1
	if cmp -s "$work/$c.tsv" "$fx/expected_${c}_nw_trkb.tsv"; then echo "$c: output equals its golden"; else echo "$c: output differs from its golden" >&2; exit 1; fi
}
record paired -p -z -1 "$fx/paired_2000_1.fasta" -2 "$fx/paired_2000_2.fasta"
record single -s -g "$fx/single_4000.fasta"
record unpaired_r -s -v -g "$fx/paired_2000_2.fasta"
cat "$work/paired.kswd" "$work/single.kswd" "$work/unpaired_r.kswd" | gzip -n -9 > "$out"
echo "recorded from $(git -C "$root" rev-parse --short "$commit"): $out"
