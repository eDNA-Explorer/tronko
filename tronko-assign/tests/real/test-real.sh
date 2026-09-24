#!/usr/bin/env bash
# make -C tronko-assign/tests/real test-real: this checkout's tronko-assign on production's real
# references, compared byte for byte with production's one-thread goldens in data/real/ (each set's
# goldens.tar.zst, unpacked into the cache and checked against its goldens.sha256).
#   1. fetch the reference files named in manifests/<marker>.sha256 (by URL, each checked against
#      its SHA-256; files already cached with the right sum are kept)
#   2. check every reference file against the manifest again; stop, naming the file, on any mismatch
#   3. build tronko-assign/ as the pipeline builds it (cap 25), plain and parquet
#   4. run every case with production's full command line at each of THREADS; cmp with the golden
# Environment: THREADS (default 16, production's thread count),
# TRONKO_ASSIGN_SRC (the tronko-assign/ source to build; default this checkout's), and MARKERS,
# SETS, CASES, FORMATS, TRONKO_REAL_* (lib.sh).
# Memory: MiFish about 16 GiB at 16 threads; FWH about 124 GiB at one thread, 146 GiB at 16.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
THREADS=${THREADS:-16}
SRC=${TRONKO_ASSIGN_SRC:-$REPO_ROOT/tronko-assign}
WORK=${TRONKO_REAL_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/tronko-test-real.XXXXXX")}
mkdir -p "$WORK"
log "source: $SRC; markers: $MARKERS; sets: $SETS; cases: $CASES; formats: $FORMATS; threads: $THREADS; work: $WORK"
fetch_all
verify_all
declare -A BIN
for f in $FORMATS; do
	build_pipeline "$SRC" "$WORK/build-$f" "$f"
	BIN[$f]=$WORK/build-$f/tronko-assign
	log "built $f: $(sha_of "${BIN[$f]}")"
done
fail=0 n=0
printf 'marker\tset\tcase\tformat\tthreads\tresult\n' >"$WORK/results.tsv"
for m in $MARKERS; do
	for s in $SETS; do
		gd=$(goldens_dir "$m" "$s")
		for c in $CASES; do
			for f in $FORMATS; do
				want=${gd:-/nonexistent}/$(golden_rel "$m" "$c" "$f")
				for t in $THREADS; do
					if [[ ! -f $want ]]; then
						r=NO_GOLDEN
					else
						d=$WORK/runs/$m-$s-$c-$f-t$t
						run_case "${BIN[$f]}" "$m" "$s" "$c" "$t" "$f" "$d"
						if cmp -s "$want" "$(run_output "$d" "$f")"; then r=IDENTICAL; else r=DIFFER; fi
					fi
					[[ $r == IDENTICAL ]] || fail=1
					n=$((n + 1))
					printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$m" "$s" "$c" "$f" "$t" "$r" | tee -a "$WORK/results.tsv"
				done
			done
		done
	done
done
((fail == 0)) || die "some outputs differ from the goldens or have none ($WORK/results.tsv)"
log "PASS: $n outputs byte-identical to production's one-thread goldens ($WORK/results.tsv)"
