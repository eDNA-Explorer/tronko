/* test_ksw_extend.c: the vectorised ksw_extend2 returns what BWA's scalar one returns.
 *
 * bwa_source_files/ksw.c keeps BWA's original seed-extension kernel as ksw_extend2_scalar()
 * and makes ksw_extend2() dispatch to a vectorised kernel (AVX2 if the CPU has it, SSE2
 * otherwise). The vectorised kernels are exact only if they return the same six values
 * (score, qle, tle, gtle, gscore, max_off) for every input. This test compares, for every call,
 * ksw_extend2(), ksw_extend2_sse2() and, where the CPU has AVX2, ksw_extend2_avx2() with the
 * expected values. Four parts, all deterministic, about a second in all:
 *   1. recorded calls: every record of tests/data/ksw_fixture_calls.kswd.gz, 2,115 calls with
 *      their inputs and the six outputs of the unmodified ksw_extend2 on the repository
 *      fixtures (every hundredth call of the paired, single and unpaired-reverse cases of
 *      tests/integration/test_assignment_production_parity.sh at one thread), recorded by
 *      tests/tools/record_ksw_calls.sh;
 *   2. random calls: 50,000 inputs from a fixed seed, each compared with ksw_extend2_scalar;
 *   3. targeted random calls for the rules most easily broken: h0 above the band (first-row
 *      cells beyond the band are non-zero and are read back when the band grows); gap
 *      extension above 1 with small z-drop thresholds and indels in both directions (the
 *      z-drop multiplies the diagonal offset by e_del or e_ins); tiny bands and empty targets;
 *   4. the 16-bit range: h0 + qlen * max(mat) within 8 of the bound (20000) above which
 *      ksw_extend2 runs the scalar code, on either side; and far above it (h0 of 32,768 to
 *      60,000), where a 16-bit kernel would overflow, so a missing fallback fails the test.
 *
 * Record format of a .kswd file (native byte order, little-endian on x86-64 and aarch64):
 * 20 int32 (0x4457534b "KSWD", call site, qlen, tlen, m, o_del, e_del, o_ins, e_ins, w,
 * end_bonus, zdrop, h0, score, qle, tle, gtle, gscore, max_off, 0 (unused)), then m * m bytes
 * of mat, qlen bytes of query, tlen bytes of target.
 *
 * Usage: test_ksw_extend <calls.kswd.gz>; exit status 0 when every comparison is equal.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "ksw.h"

typedef int (*kfn_t)(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat,
		int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0,
		int *qle, int *tle, int *gtle, int *gscore, int *max_off);

static const char *kname[4] = { "ksw_extend2", "scalar", "sse2", "avx2" };
static kfn_t kfn[4];
static int n_k;
static long n_cmp, n_bad;

static void call(kfn_t f, int out[6], int qlen, const uint8_t *q, int tlen, const uint8_t *t, const int8_t *mat,
		int od, int ed, int oi, int ei, int w, int eb, int zd, int h0)
{
	out[0] = f(qlen, q, tlen, t, 5, mat, od, ed, oi, ei, w, eb, zd, h0, &out[1], &out[2], &out[3], &out[4], &out[5]);
}

/* compare every kernel with ref[6]; print the first few differences */
static void check(const char *what, long id, const int ref[6], int qlen, const uint8_t *q, int tlen, const uint8_t *t,
		const int8_t *mat, int od, int ed, int oi, int ei, int w, int eb, int zd, int h0)
{
	int k, j, out[6];
	for (k = 0; k < n_k; ++k) {
		call(kfn[k], out, qlen, q, tlen, t, mat, od, ed, oi, ei, w, eb, zd, h0);
		++n_cmp;
		for (j = 0; j < 6; ++j) if (out[j] != ref[j]) break;
		if (j < 6 && n_bad++ < 10)
			fprintf(stderr, "FAIL %s %ld %s: qlen %d tlen %d w %d h0 %d o %d/%d e %d/%d zdrop %d: got %d %d %d %d %d %d want %d %d %d %d %d %d\n",
				what, id, kname[k], qlen, tlen, w, h0, od, oi, ed, ei, zd, out[0], out[1], out[2], out[3], out[4], out[5],
				ref[0], ref[1], ref[2], ref[3], ref[4], ref[5]);
	}
}

static long test_recorded(const char *fn)
{
	gzFile fp = gzopen(fn, "rb");
	int32_t h[20];
	long n = 0;
	if (fp == 0) { fprintf(stderr, "cannot open %s\n", fn); exit(2); }
	while (gzread(fp, h, 80) == 80) {
		int8_t mat[25];
		int ref[6], j;
		uint8_t *q, *t;
		if (h[0] != 0x4457534b || h[4] != 5 || h[2] < 0 || h[3] < 0) { fprintf(stderr, "bad record %ld\n", n); exit(2); }
		q = malloc(h[2] + 1); t = malloc(h[3] + 1);
		if (gzread(fp, mat, 25) != 25 || gzread(fp, q, h[2]) != h[2] || gzread(fp, t, h[3]) != h[3]) { fprintf(stderr, "truncated\n"); exit(2); }
		for (j = 0; j < 6; ++j) ref[j] = h[13 + j];
		check("recorded", n, ref, h[2], q, h[3], t, mat, h[5], h[6], h[7], h[8], h[9], h[10], h[11], h[12]);
		free(q); free(t);
		++n;
	}
	gzclose(fp);
	return n;
}

static uint64_t rs = 20260923;
static int rr(int lo, int hi)
{
	rs = rs * 6364136223846793005ULL + 1442695040888963407ULL;
	return lo + (int)((uint32_t)(rs >> 33) % (uint32_t)(hi - lo + 1));
}

/* mode 0: general; 1: h0 above the band; 2: gap extension 2 to 6; 3: tiny band or empty target;
 * 4: h0 + qlen * a within 8 of the 16-bit bound (20000) on either side, or far beyond it */
static void random_case(int mode, long id)
{
	static const char *mname[5] = { "random", "h0-above-band", "gap-extension>1", "tiny-band", "range-edge" };
	uint8_t q[600], t[900];
	int8_t mat[25];
	int qlen, tlen, i, j, k, a, b, od, ed, oi, ei, w, eb, zd, h0, sim, nrate, ref[6];
	qlen = mode == 1? rr(40, 300) : rr(1, 300);
	tlen = mode == 3 && rr(0, 3) == 0? rr(0, 2) : qlen + rr(-20, 150);
	if (tlen < 0) tlen = 0;
	nrate = rr(0, 3) == 0? rr(1, 15) : 0;
	for (j = 0; j < qlen; ++j) q[j] = rr(0, 99) < nrate? 4 : rr(0, 3);
	sim = mode == 2 || mode == 4? rr(70, 100) : rr(0, 100);
	for (i = j = 0; i < tlen; ++i) {
		int x = rr(0, 99);
		if (mode == 2 && j < qlen && x < 8) { t[i] = rr(0, 3); j += rr(1, 6); } /* query bases skipped: j runs ahead */
		else if (mode == 2 && x < 16) t[i] = rr(0, 3);                              /* target bases inserted: i runs ahead */
		else if (j < qlen && x < sim) t[i] = q[j++];
		else if (j < qlen && x < sim + 3) { t[i] = rr(0, 3); j += rr(0, 2); }
		else t[i] = rr(0, 99) < nrate? 4 : rr(0, 3);
	}
	a = rr(1, 2); b = rr(1, 6);
	for (i = k = 0; i < 4; ++i) { for (j = 0; j < 4; ++j) mat[k++] = i == j? a : -b; mat[k++] = -1; }
	for (j = 0; j < 5; ++j) mat[k++] = -1;
	od = rr(0, 12); oi = rr(0, 12);
	ed = mode == 2? rr(2, 6) : rr(1, 3); ei = mode == 2? rr(2, 6) : rr(1, 3);
	if ((mode == 0 || mode == 4) && rr(0, 2) == 0) od = oi = 6, ed = ei = 1;
	w = mode == 3? rr(0, 8) : rr(0, 3) == 0? rr(1, 40) : rr(50, 200);
	eb = rr(0, 10);
	zd = mode == 2? rr(1, 40) : rr(0, 5) == 0? rr(-3, 3) : rr(10, 150);
	if (mode == 1) { /* above the band: first-row cells beyond column w are non-zero */
		int wb = w < qlen? w : qlen;
		h0 = oi + ei * (wb + 1) + rr(1, 150);
	} else if (mode == 4) h0 = rr(0, 3)? 20000 - qlen * a + rr(-8, 8) : rr(32768, 60000);
	else h0 = rr(0, 2)? rr(1, 80) : rr(1, 300);
	ref[0] = ksw_extend2_scalar(qlen, q, tlen, t, 5, mat, od, ed, oi, ei, w, eb, zd, h0, &ref[1], &ref[2], &ref[3], &ref[4], &ref[5]);
	check(mname[mode], id, ref, qlen, q, tlen, t, mat, od, ed, oi, ei, w, eb, zd, h0);
}

int main(int argc, char *argv[])
{
	static const int n_mode[5] = { 50000, 5000, 5000, 5000, 2000 };
	long n_rec, n_rand = 0, i;
	int mode;
	if (argc < 2) { fprintf(stderr, "usage: test_ksw_extend <calls.kswd.gz>\n"); return 2; }
	kfn[n_k++] = ksw_extend2;
	kfn[n_k++] = ksw_extend2_scalar;
	kfn[n_k++] = ksw_extend2_sse2;
	if (ksw_simd_has_avx2()) kfn[n_k++] = ksw_extend2_avx2;
	n_rec = test_recorded(argv[1]);
	for (mode = 0; mode < 5; ++mode)
		for (i = 0; i < n_mode[mode]; ++i, ++n_rand) random_case(mode, i);
	printf("test_ksw_extend: %ld recorded calls, %ld random calls, %d kernels (%s), %ld comparisons, %ld differ\n",
		n_rec, n_rand, n_k, n_k == 4? "with avx2" : "without avx2", n_cmp, n_bad);
	return n_bad != 0 || n_rec == 0;
}
