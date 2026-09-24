#!/usr/bin/env bash
# Records tronko-assign/tests/data/chain_flt_fixture.chf.gz: the input and the output of BWA's
# mem_chain_flt() for the 1st, 41st, 81st, ... call and for every call with 250 chains or more,
# of the calls the three fixture cases of tests/integration/test_assignment_production_parity.sh
# at the repository root make (single, unpaired reverse, paired; 4,000, 2,000 and 4,000 calls, one
# per read) at one thread, in that order. The record format is described in
# unit/test_chain_flt.c ("part 1").
#
# Usage, from the top of a checkout:
#   bash tronko-assign/tests/tools/record_chain_flt.sh <commit> <output.chf.gz>
# <commit> is a commit with BWA's original mem_chain_flt(): 71f6ec3 or any later commit up to this
# change's parent. The script exports that commit's tronko-assign/ and HEAD's fixtures into a
# temporary directory, wraps the mem_chain_flt() call in mem_align1_core() there with a recorder,
# builds it, runs the three cases, checks each output against its golden and writes
# <output.chf.gz>. The checkout itself is not touched.
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 <commit> <output.chf.gz>" >&2; exit 2; }
COMMIT=$1
OUT=$(realpath -m "$2")
TOP=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/record-chf.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

# the code from <commit>; the fixtures and the example reference with its BWA index from HEAD
git -C "$TOP" archive "$COMMIT" tronko-assign | tar -x -C "$WORK"
git -C "$TOP" archive HEAD tests/data/assignment tronko-build/example_datasets/single_tree | tar -x -C "$WORK"

SRC=$WORK/tronko-assign/bwa_source_files/bwamem.c
CALL='chn.n = mem_chain_flt(opt, chn.n, chn.a);'
[[ $(grep -cF "$CALL" "$SRC") -eq 1 ]] || { echo "mem_align1_core() not as expected in $COMMIT" >&2; exit 1; }
sed -i 's/chn\.n = mem_chain_flt(opt, chn\.n, chn\.a);/chn.n = chf_record(opt, chn.n, chn.a);/' "$SRC"
sed -i '0,/^mem_alnreg_v mem_align1_core(/s//static int chf_record(const mem_opt_t *opt, int n_chn, mem_chain_t *a);\n&/' "$SRC"
cat >>"$SRC" <<'EOF'

/* Recorder for tests/data/chain_flt_fixture.chf.gz (one thread only): writes the options, the
   input chains, the returned count and, per surviving chain, its input index, w, kept and first
   to the file named by TRONKO_CHF_DUMP. A surviving chain is found by its seeds pointer, which
   mem_chain_flt() moves with the chain. */
static int chf_record(const mem_opt_t *opt, int n_chn, mem_chain_t *a)
{
	static FILE *fp = NULL;
	static long calls = 0;
	int i, j, s, k, record;
	int32_t marker = 0x43484631, v[6];
	mem_seed_t **seeds;
	if (calls == 0 && getenv("TRONKO_CHF_DUMP")) fp = fopen(getenv("TRONKO_CHF_DUMP"), "wb");
	record = fp && n_chn > 0 && (calls % 40 == 0 || n_chn >= 250);
	calls++;
	if (!record) return mem_chain_flt(opt, n_chn, a);
	v[0] = opt->min_chain_weight; memcpy(&v[1], &opt->drop_ratio, 4); memcpy(&v[2], &opt->mask_level, 4);
	v[3] = opt->max_chain_gap; v[4] = opt->min_seed_len; v[5] = opt->max_chain_extend;
	fwrite(&marker, 4, 1, fp); fwrite(v, 4, 6, fp); fwrite(&n_chn, 4, 1, fp);
	seeds = malloc(n_chn * sizeof(*seeds));
	for (i = 0; i < n_chn; ++i) {
		int32_t h[3] = { a[i].n, a[i].rid, a[i].is_alt };
		fwrite(h, 4, 3, fp); fwrite(&a[i].frac_rep, 4, 1, fp); fwrite(&a[i].pos, 8, 1, fp);
		for (s = 0; s < a[i].n; ++s) {
			int32_t q[3] = { a[i].seeds[s].qbeg, a[i].seeds[s].len, a[i].seeds[s].score };
			fwrite(&a[i].seeds[s].rbeg, 8, 1, fp); fwrite(q, 4, 3, fp);
		}
		seeds[i] = a[i].seeds;
	}
	k = mem_chain_flt(opt, n_chn, a);
	fwrite(&k, 4, 1, fp);
	for (i = 0; i < k; ++i) {
		int32_t r[4];
		for (j = 0; j < n_chn && seeds[j] != a[i].seeds; ++j);
		r[0] = j; r[1] = (int32_t)a[i].w; r[2] = a[i].kept; r[3] = a[i].first;
		fwrite(r, 4, 4, fp);
	}
	free(seeds);
	return k;
}
EOF
make -C "$WORK/tronko-assign" >"$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log" >&2; exit 1; }

D=$WORK/tests/data/assignment
COMMON=(-r -f "$D/reference_tree.trkb" -a "$WORK/tronko-build/example_datasets/single_tree/Charadriiformes.fasta"
	-w -6 --Cinterval 10 --number-of-cores 1)
run() { # run <label> <tronko args...>
	local label=$1; shift
	TRONKO_CHF_DUMP="$WORK/$label.chf" "$WORK/tronko-assign/tronko-assign" "${COMMON[@]}" "$@" \
		-o "$WORK/$label.tsv" >/dev/null 2>"$WORK/$label.log"
	cmp -s "$WORK/$label.tsv" "$D/expected_$label.tsv" || { echo "$label: output differs from its golden" >&2; exit 1; }
}
run single_nw_trkb -s -g "$D/single_4000.fasta"
run unpaired_r_nw_trkb -s -v -g "$D/paired_2000_2.fasta"
run paired_nw_trkb -p -z -1 "$D/paired_2000_1.fasta" -2 "$D/paired_2000_2.fasta"

cat "$WORK"/single_nw_trkb.chf "$WORK"/unpaired_r_nw_trkb.chf "$WORK"/paired_nw_trkb.chf | gzip -9n >"$OUT"
echo "$OUT: $(stat -c %s "$OUT") bytes"
