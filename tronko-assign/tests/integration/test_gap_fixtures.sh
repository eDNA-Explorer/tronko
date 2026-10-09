#!/usr/bin/env bash
# Parity of tronko-assign on the coverage-gap fixtures: inputs that make it run code production
# executes and the other fixtures never reach. Each case must reproduce its golden, made by
# production tronko-assign (71f6ec3) at --number-of-cores 1 (mr_paired: by the code with the
# mate-rescue fix, tests/README.md), at every thread count in TRONKO_ASSIGN_CORES.
#
#   mr_paired      620 pairs in which one mate cannot be seeded, so BWA runs mate rescue
#                  (ksw_align2, including the 16-bit kernel for the 40 mates of 260 bases)
#   pc_single, pc_unpaired_r, pc_paired
#                  production's own command line: reads as .fasta.zst (one frame, no checksum,
#                  as the pipeline writes them) and -R -T --tsv-log <file> -V2
#                  (zstd reader, logging, resource monitor, crash-handler set-up); the goldens
#                  equal those of tests/data/assignment
#   rc_single, rc_unpaired_r, rc_paired
#                  reads production's QC lets through: one N, no alignment at all (random
#                  sequence), 30 to 60 bases
#   pn_single, pn_unpaired_r, pn_paired
#                  read names as the pipeline writes them (12S_MiFish_U_paired_F_<idx>, ...)
#
# All cases use the example Charadriiformes reference (tronko-build/example_datasets/single_tree,
# tests/data/assignment/reference_tree.trkb). Inputs: tests/data/gaps/<fixture>/ under
# tronko-assign; goldens: tests/data/gaps/<fixture>/goldens[-cap25]/expected_<case>.tsv.
# Environment: see golden_lib.sh; TRONKO_GAP_CASES limits the cases run.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=golden_lib.sh
source "$SCRIPT_DIR/golden_lib.sh"
GAPS="$TA/tests/data/gaps"
DATA="$REPO_ROOT/tests/data/assignment"
EX="$REPO_ROOT/tronko-build/example_datasets/single_tree"
FASTA="$EX/Charadriiformes.fasta"
REFERENCE="$DATA/reference_tree.trkb"
CASES=${TRONKO_GAP_CASES:-"mr_paired pc_single pc_unpaired_r pc_paired rc_single rc_unpaired_r rc_paired pn_single pn_unpaired_r pn_paired"}

for f in "$FASTA" "$REFERENCE" "$GAPS/mate-rescue/mr_1.fasta" "$GAPS/prod-cmdline/single_4000.fasta.zst" \
	"$GAPS/read-content/rc_single.fasta" "$GAPS/prod-names/pn_paired_F.fasta"; do
	[[ -f $f ]] || { echo "SKIP: missing fixture file $f"; exit 77; }
done
mem3_index "$FASTA" || exit 2
stage_mate_rescue "$FASTA" "$FASTA".{amb,ann,pac,bwt.2bit.64} "$REFERENCE" "$GAPS/mate-rescue/mr_1.fasta" "$GAPS/mate-rescue/mr_2.fasta" ||
	{ echo "cannot stage the mate-rescue inputs in $MATE_RESCUE_DIR" >&2; exit 2; }

# case_args <case>: the argument list without --number-of-cores and the output option, in
# production's order for the pc_ cases
case_args() {
	local M=$MATE_RESCUE_DIR P=$GAPS/prod-cmdline Q=$GAPS/read-content N=$GAPS/prod-names
	local common=(-r -f "$REFERENCE" -a "$FASTA" -w -6 --Cinterval 10)
	local tel=(-R -T --tsv-log "tsvlog_$1.tsv" -V2)
	case $1 in
	mr_paired) ARGS=(-r -f "$M/reference_tree.trkb" -a "$M/Charadriiformes.fasta" -w -6 --Cinterval 10 -p -z -1 "$M/mr_1.fasta" -2 "$M/mr_2.fasta") ;;
	pc_single) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -s -w -g "$P/single_4000.fasta.zst" -6 --Cinterval 10 "${tel[@]}") ;;
	pc_unpaired_r) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -s -w -g "$P/paired_2000_2.fasta.zst" -v -6 --Cinterval 10 "${tel[@]}") ;;
	pc_paired) ARGS=(-r -f "$REFERENCE" -a "$FASTA" -p -z -w -1 "$P/paired_2000_1.fasta.zst" -2 "$P/paired_2000_2.fasta.zst" -6 --Cinterval 10 "${tel[@]}") ;;
	rc_single) ARGS=("${common[@]}" -s -g "$Q/rc_single.fasta") ;;
	rc_unpaired_r) ARGS=("${common[@]}" -s -v -g "$Q/rc_single.fasta") ;;
	rc_paired) ARGS=("${common[@]}" -p -z -1 "$Q/rc_1.fasta" -2 "$Q/rc_2.fasta") ;;
	pn_single) ARGS=("${common[@]}" -s -g "$N/pn_unpaired_F.fasta") ;;
	pn_unpaired_r) ARGS=("${common[@]}" -s -v -g "$N/pn_unpaired_R.fasta") ;;
	pn_paired) ARGS=("${common[@]}" -p -z -1 "$N/pn_paired_F.fasta" -2 "$N/pn_paired_R.fasta") ;;
	*) echo "unknown gap case: $1" >&2; return 1 ;;
	esac
}

fixture_of() {
	case $1 in mr_*) echo mate-rescue ;; pc_*) echo prod-cmdline ;; rc_*) echo read-content ;; pn_*) echo prod-names ;; esac
}

echo "gap fixtures: cap $MATCH_CAP, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"
for c in $CASES; do
	case_args "$c" || exit 2
	golden=$(golden_dir "$GAPS/$(fixture_of "$c")/goldens")/expected_$c.tsv
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		run_assign "$TMP_DIR/$c-c$n.log" "${ARGS[@]}" --number-of-cores "$n" -o "$out" || { FAILS=$((FAILS + 1)); continue; }
		check_output "$c" "$golden" "$out" "$c, --number-of-cores $n"
	done
done
((FAILS == 0)) || { echo "gap fixtures: $FAILS failure(s)" >&2; exit 1; }
