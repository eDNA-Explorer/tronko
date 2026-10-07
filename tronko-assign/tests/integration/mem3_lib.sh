# BWA-MEM3 for the golden tests (sourced by golden_lib.sh and the repository's production-parity
# script). tronko-assign runs BWA-MEM3 as a separate program and reads its index beside the -a
# FASTA (<fasta>.amb .ann .pac .bwt.2bit.64). The fixtures commit only the FASTA; the index is
# built from it here, at test time (well under a second per fixture reference), and rebuilt when
# the FASTA is newer than it.
#
# Sets and exports TRONKO_BWA_MEM3 (default: tronko-assign/bwa-mem3 of this checkout) so that
# every tronko-assign the tests run, wherever its binary lives, uses that BWA-MEM3.

MEM3_LIB_TA=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
export TRONKO_BWA_MEM3=${TRONKO_BWA_MEM3:-"$MEM3_LIB_TA/bwa-mem3"}
if [[ ! -x $TRONKO_BWA_MEM3 ]]; then
	echo "Missing BWA-MEM3: $TRONKO_BWA_MEM3 (run make bwa-mem3 in tronko-assign, or set TRONKO_BWA_MEM3)" >&2
	exit 1
fi

# mem3_index <fasta>...: build each FASTA's BWA-MEM3 index beside it unless it is there and newer
# than the FASTA. Built in a scratch directory and moved into place file by file, so a concurrent
# run never reads a half-written index.
mem3_index() {
	local fa tmp ext
	for fa in "$@"; do
		if [[ -f $fa.bwt.2bit.64 && $fa.bwt.2bit.64 -nt $fa && -f $fa.pac && -f $fa.ann && -f $fa.amb ]]; then
			continue
		fi
		tmp=$(mktemp -d "${TMPDIR:-/tmp}/mem3-index.XXXXXX") || return 1
		cp "$fa" "$tmp/ref.fasta" || { rm -rf "$tmp"; return 1; }
		if ! "$TRONKO_BWA_MEM3" index "$tmp/ref.fasta" >"$tmp/index.log" 2>&1; then
			echo "bwa-mem3 index failed for $fa:" >&2
			tail -5 "$tmp/index.log" >&2
			rm -rf "$tmp"
			return 1
		fi
		for ext in amb ann pac bwt.2bit.64; do
			cp "$tmp/ref.fasta.$ext" "$fa.$ext.tmp.$$" && mv -f "$fa.$ext.tmp.$$" "$fa.$ext" || { rm -rf "$tmp"; return 1; }
		done
		rm -rf "$tmp"
	done
}
