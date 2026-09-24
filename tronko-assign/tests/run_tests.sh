#!/usr/bin/env bash
# The tests of tronko-assign, as `make test` runs them (from tronko-assign/, or
# `make -C tronko-assign/tests test` from the repository root):
#   ../../tests/integration/test_assignment_production_parity.sh   repository fixture, three read modes
#   integration/test_multibatch_parity.sh     repository fixture at -L 400 (ten batches)
#   integration/test_multitree_parity.sh      three-tree fixture (reads with several candidate trees)
#   integration/test_gap_fixtures.sh          mate rescue, production's command line, read content, read names
# Every integration test compares tronko-assign's output byte for byte with goldens made by
# production tronko-assign at one thread; they need tronko-assign built (make). A test that finds
# no fixture reports SKIP (exit 77), which is not a failure here. Environment: golden_lib.sh.
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
run tests/integration/test_assignment_production_parity.sh bash "$ROOT/tests/integration/test_assignment_production_parity.sh"
run integration/test_multibatch_parity.sh bash "$HERE/integration/test_multibatch_parity.sh"
run integration/test_multitree_parity.sh bash "$HERE/integration/test_multitree_parity.sh"
run integration/test_gap_fixtures.sh bash "$HERE/integration/test_gap_fixtures.sh"
exit $fail
