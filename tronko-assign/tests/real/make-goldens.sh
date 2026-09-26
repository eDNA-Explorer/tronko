#!/usr/bin/env bash
# make -C tronko-assign/tests/real goldens: make the real-reference goldens again from production's
# release (commit 71f6ec3) at one thread, built with the pipeline's recipe, and compare every file
# with the golden committed under data/real/ (the set's goldens.tar.zst).
# make -C tronko-assign/tests/real update-goldens: the same, and write every golden that is new or
# differs into the set's goldens.tar.zst (with its command line and time -v record under
# provenance/ in the archive, and a record in the set's provenance.txt), then rewrite the set's
# goldens.sha256. A set whose goldens are all identical is left untouched, so `git status` names
# exactly the sets that changed and `git diff` of goldens.sha256 the goldens; for each golden that
# differs the script prints the number of changed rows and the first of them.
#   1. fetch and check the reference files of manifests/<marker>.sha256 (as test-real)
#   2. export GOLDEN_COMMIT (default 71f6ec3) with its carquet submodule from this repository
#   3. build it as the pipeline builds it (cap 25), plain and parquet
#   4. run every case at --number-of-cores 1 with production's full command line
#   5. compare each output with the set's golden: IDENTICAL, DIFFER or NEW
# A new set needs only its read pairs, data/real/<marker>/<set>/<set>_{F,R}.fasta.zst, and a
# reads.json naming their source; a new marker also needs manifests/<marker>.sha256.
# Environment: GOLDEN_COMMIT, GOLDENS_UPDATE (1: write, as update-goldens does), and MARKERS, SETS,
# CASES, FORMATS, TRONKO_REAL_* (lib.sh). Measured at one thread: FWH dev 23 to 25 minutes per case.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
COMMIT=${GOLDEN_COMMIT:-71f6ec3}
UPDATE=${GOLDENS_UPDATE:-0}
WORK=${TRONKO_REAL_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/tronko-goldens.XXXXXX")}
mkdir -p "$WORK"
log "commit $COMMIT; markers: $MARKERS; sets: $SETS; cases: $CASES; formats: $FORMATS; update: $UPDATE; work: $WORK"
fetch_all
verify_all
full=$(git -C "$REPO_ROOT" rev-parse --verify "$COMMIT^{commit}") || die "commit $COMMIT is not in $REPO_ROOT (git fetch origin high-perf)"
export_tree "$full" "$WORK/src"
declare -A BIN
for f in $FORMATS; do
	build_pipeline "$WORK/src/tronko-assign" "$WORK/build-$f" "$f"
	BIN[$f]=$WORK/build-$f/tronko-assign
	log "built $COMMIT $f: $(sha_of "${BIN[$f]}")"
done
now=$(date -u +%FT%TZ)
cpu=$(lscpu 2>/dev/null | awk -F: '/^Model name/ {sub(/^ +/, "", $2); if ($2 != "-") print $2; exit}')
machine="$(uname -srm); ${cpu:-unknown CPU}; $(nproc) CPUs; $(gcc --version | head -1)"
fail=0 n=0 written=0
printf 'marker\tset\tcase\tformat\tsha256\tresult\n' >"$WORK/results.tsv"
for m in $MARKERS; do
	for s in $SETS; do
		sd=$(set_dir "$m" "$s") gd=$(goldens_dir "$m" "$s") stage=$WORK/stage/$m/$s
		changed=()
		if [[ $UPDATE == 1 ]]; then
			mkdir -p "$stage"
			[[ -n $gd ]] && cp -a "$gd/." "$stage/" && rm -f "$stage/.checked"
		fi
		for c in $CASES; do
			for f in $FORMATS; do
				d=$WORK/runs/$m-$s-$c-$f-t1
				run_case "${BIN[$f]}" "$m" "$s" "$c" 1 "$f" "$d"
				out=$(run_output "$d" "$f") rel=$(golden_rel "$m" "$c" "$f")
				if [[ -z $gd || ! -f $gd/$rel ]]; then r=NEW
				elif cmp -s "$gd/$rel" "$out"; then r=IDENTICAL
				else r=DIFFER; fi
				n=$((n + 1))
				printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$s" "$c" "$f" "$(sha_of "$out")" "$r" | tee -a "$WORK/results.tsv"
				[[ $r == IDENTICAL ]] && continue
				if [[ $r == DIFFER && $f == tsv ]]; then
					log "$m/$s $rel: $( (diff "$gd/$rel" "$out" || true) | grep -c '^>' || true) rows differ; the first:"
					(diff "$gd/$rel" "$out" || true) | grep '^[<>]' | head -6 >&2 || true
				fi
				if [[ $UPDATE != 1 ]]; then fail=1; continue; fi
				mkdir -p "$stage/$(dirname "$rel")" "$stage/provenance"
				cp "$out" "$stage/$rel"
				rm -f "$stage/provenance/${m}_$c.$f."*
				sed -e "s#$WORK/#<work>/#g" -e "s#$CACHE/#<cache>/#g" -e "s#$REPO_ROOT/#<repo>/#g" "$d/cmdline" >"$stage/provenance/${m}_$c.$f.cmd.txt"
				[[ -f $d/time.txt ]] && cp "$d/time.txt" "$stage/provenance/${m}_$c.$f.time.txt"
				changed+=("$rel ($r)")
				written=$((written + 1))
			done
		done
		if ((${#changed[@]})); then
			pack_goldens "$stage" "$m" "$s"
			{
				[[ -f $sd/provenance.txt ]] || echo "tronko-assign/tests/data/real/$m/$s: real-reference goldens of tronko-assign"
				echo
				echo "made by make-goldens.sh on $now: $(printf '%s, ' "${changed[@]}" | sed 's/, $//')"
				echo "  commit:  $COMMIT ($full)"
				echo "  build:   the pipeline's recipe, MAX_NUM_BWA_MATCHES 25, CC=$PIPELINE_CC"
				echo "  threads: --number-of-cores 1; production's full command line (provenance/<golden>.cmd.txt in goldens.tar.zst)"
				echo "  machine: $machine"
				for f in $FORMATS; do echo "  binary $f: $(sha_of "${BIN[$f]}")"; done
			} >>"$sd/provenance.txt"
		fi
	done
done
if [[ $UPDATE == 1 ]]; then
	log "$n goldens made from $COMMIT; $written written into $DATA_DIR ($WORK/results.tsv); review with git diff (goldens.sha256)"
else
	((fail == 0)) || die "goldens made from $COMMIT differ from the committed ones or are new ($WORK/results.tsv); make update-goldens writes them"
	log "PASS: $n goldens made again from $COMMIT, every one identical to the committed golden"
fi
