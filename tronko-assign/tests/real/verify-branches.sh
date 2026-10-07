#!/usr/bin/env bash
# make -C tronko-assign/tests/real verify-branches BRANCHES="<rev> ...": test-real on each branch
# or commit, against the goldens and references of this checkout, with one summary at the end.
# Each revision's tronko-assign/ (and its carquet submodule) is exported from this repository and
# built as the pipeline builds it; the reference files are fetched once. A revision that fails to
# build or run is reported as ERROR and the next one still runs.
# Environment: BRANCHES (required), THREADS (default 16), and MARKERS, SETS, CASES, FORMATS,
# TRONKO_REAL_* (lib.sh). Writes <work>/summary.tsv, one row per revision and output.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
[[ -n ${BRANCHES:-} ]] || die "set BRANCHES, e.g. BRANCHES=\"main my-branch\""
export THREADS=${THREADS:-16}
WORK=${TRONKO_REAL_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/tronko-verify.XXXXXX")}
mkdir -p "$WORK"
log "branches: $BRANCHES; threads: $THREADS; work: $WORK"
fetch_all
verify_all
printf 'branch\tcommit\tmarker\tset\tcase\tformat\tthreads\tresult\n' >"$WORK/summary.tsv"
fail=0
for b in $BRANCHES; do
	key=$(echo "$b" | tr '/' '_')
	if ! commit=$(git -C "$REPO_ROOT" rev-parse --verify --quiet "$b^{commit}"); then
		printf '%s\t-\t-\t-\t-\t-\t-\tERROR (no such revision)\n' "$b" >>"$WORK/summary.tsv"; fail=1; continue
	fi
	log "== $b ($commit)"
	if (export_tree "$commit" "$WORK/$key/src") &&
		TRONKO_ASSIGN_SRC=$WORK/$key/src/tronko-assign TRONKO_REAL_WORK=$WORK/$key bash "$REAL_DIR/test-real.sh"; then
		:
	else
		fail=1
	fi
	if [[ -s $WORK/$key/results.tsv ]] && (($(wc -l <"$WORK/$key/results.tsv") > 1)); then
		tail -n +2 "$WORK/$key/results.tsv" | sed "s#^#$b\t${commit:0:12}\t#" >>"$WORK/summary.tsv"
	else
		printf '%s\t%s\t-\t-\t-\t-\t-\tERROR (see %s)\n' "$b" "${commit:0:12}" "$WORK/$key" >>"$WORK/summary.tsv"
	fi
done
echo
if command -v column >/dev/null; then column -t -s $'\t' "$WORK/summary.tsv"; else cat "$WORK/summary.tsv"; fi
echo
((fail == 0)) || die "not every branch reproduces the goldens ($WORK/summary.tsv)"
log "PASS: every output of every branch byte-identical to the goldens ($WORK/summary.tsv)"
