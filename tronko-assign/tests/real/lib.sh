# Shared by test-real.sh, make-goldens.sh and verify-branches.sh: the reference manifest, fetching
# by URL with its checksum, the committed read sets and goldens, the pipeline's build recipe and
# production's command line. Sourced, not run.
#
# Environment:
#   MARKERS              markers to run (default: mifish); each has manifests/<marker>.sha256 and
#                        data/real/<marker>/<set>/
#   SETS                 read sets (default: dev5k); dev5k, dev, heldout
#   CASES                read modes (default: paired unpF unpR)
#   FORMATS              outputs (default: tsv parquet); parquet needs the carquet submodule
#   TRONKO_REAL_CACHE    where fetched reference files are kept (default: ~/.cache/tronko-test-real)
#   TRONKO_REAL_MIRROR   a local directory laid out like the bucket (gs://<bucket>/<path> is read
#                        from $TRONKO_REAL_MIRROR/<path>); without it files are fetched with
#                        `gcloud storage cp`, `gsutil cp`, or curl from storage.googleapis.com
#   TRONKO_REAL_WORK     scratch for builds and outputs (default: a new mktemp -d directory)
set -euo pipefail
REAL_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TESTS_DIR=$(cd "$REAL_DIR/.." && pwd)
REPO_ROOT=$(cd "$TESTS_DIR/../.." && pwd)
DATA_DIR=$TESTS_DIR/data/real
MARKERS=${MARKERS:-mifish}
SETS=${SETS:-dev5k}
CASES=${CASES:-paired unpF unpR}
FORMATS=${FORMATS:-tsv parquet}
CACHE=${TRONKO_REAL_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/tronko-test-real}
MIRROR=${TRONKO_REAL_MIRROR:-}
# The eDNA Explorer pipeline's build: the match cap patched to 25, then make with this CC;
# ENABLE_PARQUET=1 for the parquet build.
PIPELINE_CC="gcc -O3 -fcommon -Wno-error -Wno-implicit-function-declaration -Wno-incompatible-pointer-types -Wno-int-conversion"

die() { echo "real: $*" >&2; exit 1; }
log() { echo "real: $*" >&2; }

# manifest_lines <marker>: "<sha256> <url>" for each reference file of the marker
manifest_lines() {
	local f=$TESTS_DIR/manifests/$1.sha256
	[[ -f $f ]] || die "no manifest for marker '$1' ($f)"
	grep -v '^#' "$f" | awk 'NF == 2'
}
cache_path() { local u=${1#gs://}; echo "$CACHE/$u"; }
sha_of() { sha256sum "$1" | cut -d' ' -f1; }

# fetch_one <sha256> <url>: into the cache, unless a file with that checksum is there already
fetch_one() {
	local sha=$1 url=$2 dest part rel bucket
	dest=$(cache_path "$url")
	[[ -f $dest && $(sha_of "$dest") == "$sha" ]] && return 0
	mkdir -p "$(dirname "$dest")"
	part=$dest.part
	bucket=${url#gs://}; bucket=${bucket%%/*}; rel=${url#gs://$bucket/}
	if [[ -n $MIRROR ]]; then
		[[ -f $MIRROR/$rel ]] || die "not in the mirror: $url ($MIRROR/$rel)"
		cp "$MIRROR/$rel" "$part"
	elif command -v gcloud >/dev/null; then
		gcloud storage cp --quiet "$url" "$part" >&2 || die "download failed: $url"
	elif command -v gsutil >/dev/null; then
		gsutil -q cp "$url" "$part" || die "download failed: $url"
	else
		curl -fsSL -o "$part" "https://storage.googleapis.com/$bucket/$rel" || die "download failed: $url"
	fi
	local got; got=$(sha_of "$part")
	[[ $got == "$sha" ]] || die "checksum mismatch for $url: expected $sha, got $got (kept as $part)"
	mv "$part" "$dest"
	log "fetched $url"
}

# fetch_all: every reference file of the chosen markers
fetch_all() {
	local m sha url
	for m in $MARKERS; do
		while read -r sha url; do fetch_one "$sha" "$url"; done < <(manifest_lines "$m")
	done
}

# verify_all: every cached reference file against the manifest, and the committed read pairs of
# every chosen set, before a run; fails naming the file
verify_all() {
	local m s sha url f got bad=0
	for m in $MARKERS; do
		while read -r sha url; do
			f=$(cache_path "$url")
			if [[ ! -f $f ]]; then echo "real: missing: $url ($f)" >&2; bad=1; continue; fi
			got=$(sha_of "$f")
			[[ $got == "$sha" ]] || { echo "real: checksum mismatch: $f ($url): expected $sha, got $got" >&2; bad=1; }
		done < <(manifest_lines "$m")
		for s in $SETS; do
			for f in "$(set_dir "$m" "$s")/${s}_F.fasta.zst" "$(set_dir "$m" "$s")/${s}_R.fasta.zst"; do
				[[ -f $f ]] || { echo "real: no read set: $f" >&2; bad=1; }
			done
		done
	done
	((bad == 0)) || die "inputs are missing or do not match the manifests; nothing was run"
	log "every reference file matches manifests/{$(echo $MARKERS | tr ' ' ,)}.sha256; read sets present"
}

# ref_trkb / ref_fasta <marker>: the cached reference files named in the marker's manifest
ref_trkb() { cache_path "$(manifest_lines "$1" | awk '$2 ~ /\/reference_tree\.trkb$/ {print $2}')"; }
ref_fasta() { cache_path "$(manifest_lines "$1" | awk '$2 ~ /\.fasta$/ {print $2}')"; }
set_dir() { echo "$DATA_DIR/$1/$2"; }
# golden_rel <marker> <case> <tsv|parquet>: a golden's path inside the set's goldens.tar.zst
golden_rel() { if [[ $3 == tsv ]]; then echo "$1_$2.tsv"; else echo "parquet/$1_$2.parquet"; fi; }

# goldens_dir <marker> <set>: the set's goldens.tar.zst unpacked into the cache (once per archive),
# every member checked against the set's goldens.sha256; prints the directory, or nothing if the set
# has no goldens yet
goldens_dir() {
	local sd a d
	sd=$(set_dir "$1" "$2") a=$(set_dir "$1" "$2")/goldens.tar.zst
	[[ -f $a ]] || return 0
	d=$CACHE/goldens/$1/$2/$(sha_of "$a")
	if [[ ! -f $d/.checked ]]; then
		rm -rf "$d" && mkdir -p "$d"
		zstd -qdc "$a" | tar -x -C "$d" || die "cannot unpack $a"
		(cd "$d" && sha256sum --check --quiet "$sd/goldens.sha256") || die "the goldens of $1/$2 do not match goldens.sha256"
		touch "$d/.checked"
	fi
	echo "$d"
}

# pack_goldens <dir> <marker> <set>: every file under <dir> (goldens, parquet/, provenance/) into the
# set's goldens.tar.zst, and their SHA-256 into goldens.sha256; deterministic for a given zstd
# (sorted names, fixed times and owners), so an unchanged set packs to the same bytes
pack_goldens() {
	local src=$1 sd; sd=$(set_dir "$2" "$3")
	(cd "$src" && find . -type f ! -name .checked | sed 's#^\./##' | LC_ALL=C sort | xargs sha256sum) >"$sd/goldens.sha256"
	(cd "$src" && find . -type f ! -name .checked | sed 's#^\./##' | LC_ALL=C sort |
		tar --no-recursion --mtime=@0 --owner=0 --group=0 --numeric-owner --format=gnu -T - -cf -) |
		zstd -q -19 --no-check -c >"$sd/goldens.tar.zst"
}

# export_tree <commit> <dir>: <dir>/tronko-assign from the commit, with the commit's carquet submodule
export_tree() {
	local full sub
	full=$(git -C "$REPO_ROOT" rev-parse --verify "$1^{commit}") || die "commit $1 is not in $REPO_ROOT (git fetch origin)"
	mkdir -p "$2"
	git -C "$REPO_ROOT" archive "$full" tronko-assign | tar -x -C "$2"
	sub=$(git -C "$REPO_ROOT" rev-parse "$full:tronko-assign/carquet")
	git -C "$REPO_ROOT/tronko-assign/carquet" cat-file -e "$sub^{commit}" 2>/dev/null ||
		die "carquet $sub (of $1) is not available: run git submodule update --init tronko-assign/carquet, then git -C tronko-assign/carquet fetch"
	mkdir -p "$2/tronko-assign/carquet"
	git -C "$REPO_ROOT/tronko-assign/carquet" archive "$sub" | tar -x -C "$2/tronko-assign/carquet"
}

# build_pipeline <tronko-assign source dir> <out dir> <tsv|parquet>: a copy built as the pipeline builds
build_pipeline() {
	local src=$1 out=$2 fmt=$3
	mkdir -p "$out"
	cp -a "$src/." "$out/"
	(
		cd "$out"
		make clean >/dev/null 2>&1 || true
		sed -i 's/#define MAX_NUM_BWA_MATCHES 10/#define MAX_NUM_BWA_MATCHES 25/' global.h
		grep -q '#define MAX_NUM_BWA_MATCHES 25' global.h || die "the cap patch did not apply to $out/global.h"
		if [[ $fmt == parquet ]]; then make -j"$(nproc)" ENABLE_PARQUET=1 CC="$PIPELINE_CC"; else make -j"$(nproc)" CC="$PIPELINE_CC"; fi
	) >"$out.build.log" 2>&1 || die "build failed ($fmt), log: $out.build.log"
	[[ -x $out/tronko-assign ]] || die "no binary after the $fmt build ($out)"
}

# run_case <binary> <marker> <set> <case> <threads> <tsv|parquet> <run dir>: production's full
# command line (-R -T --tsv-log <file> -V2, .fasta.zst reads); writes out.tsv or out.parquet
run_case() {
	local bin=$1 m=$2 s=$3 c=$4 n=$5 fmt=$6 d=$7 sd
	local -a mode out argv tm=()
	sd=$(set_dir "$m" "$s")
	case $c in
	paired) mode=(-r -p -z -1 "$sd/${s}_F.fasta.zst" -2 "$sd/${s}_R.fasta.zst") ;;
	unpF) mode=(-r -s -g "$sd/${s}_F.fasta.zst") ;;
	unpR) mode=(-r -s -v -g "$sd/${s}_R.fasta.zst") ;;
	*) die "unknown case $c" ;;
	esac
	if [[ $fmt == tsv ]]; then out=(-o "$d/out.tsv"); else out=(--parquet "$d/out"); fi
	argv=("${mode[@]}" -f "$(ref_trkb "$m")" -a "$(ref_fasta "$m")" -w -6 --Cinterval 10
		--number-of-cores "$n" -R -T --tsv-log "$d/tsvlog.tsv" -V2 "${out[@]}")
	mkdir -p "$d"
	echo "tronko-assign ${argv[*]}" >"$d/cmdline"
	[[ -x /usr/bin/time ]] && tm=(/usr/bin/time -v -o "$d/time.txt")
	(cd "$d" && "${tm[@]}" "$bin" "${argv[@]}") >"$d/stdout" 2>"$d/stderr" || die "tronko-assign failed: $m $s $c $fmt at $n threads (see $d/stderr)"
}
# run_output <run dir> <tsv|parquet>: the output file run_case wrote
run_output() { if [[ $2 == tsv ]]; then echo "$1/out.tsv"; else echo "$1/out.parquet"; fi; }
