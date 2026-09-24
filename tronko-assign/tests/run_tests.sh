#!/usr/bin/env bash
# The tests of tronko-assign, as `make test` runs them (from tronko-assign/, or
# `make -C tronko-assign/tests test` from the repository root):
#   unit/test_nodestore                       node-store kernels and store against the legacy loop
#   unit/test_sam_seq                         SAM SEQ text for every read byte, both strands
#   unit/test_mate_rescue                     mate rescue with the score matrix against a guard page
#   ../../tests/integration/test_assignment_production_parity.sh   repository fixture, three read modes
#   integration/test_multibatch_parity.sh     repository fixture at -L 400 (ten batches)
#   integration/test_multitree_parity.sh      three-tree fixture (reads with several candidate trees)
#   integration/test_gap_fixtures.sh          mate rescue, production's command line, read content, read names
#   integration/test_nodestore_paths.sh       node-store paths on the fixtures: self-check, block shapes, kernels
#   integration/test_path_options.sh          path options against their buffer sizes
#   integration/test_slot_cap.sh              candidate slots at cap 2 under AddressSanitizer (~15 s)
# The golden tests compare tronko-assign's output byte for byte with goldens made at one thread;
# they, test_nodestore_paths.sh and test_path_options.sh need tronko-assign built (make).
# test_slot_cap.sh builds its own copy and needs -fsanitize=address. A test that finds no fixture
# or tool reports SKIP (exit 77), which is not a failure here. Environment of the golden tests:
# integration/golden_lib.sh.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
fail=0
run() {
	local name=$1 rc
	shift
	"$@"
	rc=$?
	case $rc in
	0) echo "== $name: PASS" ;;
	77) echo "== $name: SKIP" ;;
	*) echo "== $name: FAIL (exit $rc)"; fail=1 ;;
	esac
}
run unit/test_nodestore "$HERE/unit/test_nodestore"
run unit/test_sam_seq "$HERE/unit/test_sam_seq"
run unit/test_mate_rescue "$HERE/unit/test_mate_rescue"
run tests/integration/test_assignment_production_parity.sh bash "$ROOT/tests/integration/test_assignment_production_parity.sh"
run integration/test_multibatch_parity.sh bash "$HERE/integration/test_multibatch_parity.sh"
run integration/test_multitree_parity.sh bash "$HERE/integration/test_multitree_parity.sh"
run integration/test_gap_fixtures.sh bash "$HERE/integration/test_gap_fixtures.sh"
run integration/test_nodestore_paths.sh bash "$HERE/integration/test_nodestore_paths.sh"
run integration/test_path_options.sh bash "$HERE/integration/test_path_options.sh"
run integration/test_slot_cap.sh bash "$HERE/integration/test_slot_cap.sh"
exit $fail
