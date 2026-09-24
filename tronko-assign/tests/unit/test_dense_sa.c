/*
 * test_dense_sa.c: the dense suffix-array sample of bwa_source_files/bwt.c.
 *
 * bwt_densify_sa() builds a uint32 suffix-array sample every 1 << shift indices, and bwt_sa() then
 * stops its LF walk at the first multiple of 1 << shift instead of the file sample's interval
 * (32). Tronko's output is unchanged only if bwt_sa() returns, for every index, the value the
 * original walk over the file sample returns. This test checks that for every index of:
 *
 *   - the repository's Charadriiformes BWA index (argv[1]; primary % 4 = 2, so walks through
 *     primary step into index 0 and read the $ sentinel), built at 1, 3, 4, 5 and 7 threads;
 *   - synthetic indexes whose suffix array is computed here by a plain sort, independently of
 *     BWA: random, two-letter, periodic and single-letter texts of 1 to 270,000 bases, file
 *     intervals 32, 16 and 8, dense intervals 2 to 16, 1 to 9 threads. Every residue of
 *     primary % 4 is required to occur.
 *
 * It also checks the build's failure paths (a corrupted file sample must leave the index on the
 * file sample; the skip conditions), the TRONKO_CHECK_SA mode (both lookups; a corrupted dense
 * entry must abort), and that the file sample is freed and set to null when not kept.
 *
 * Usage: test_dense_sa <Charadriiformes.fasta prefix>. Exit status 0 when every check passes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include "bwt.h"

int bwa_verbose = 1; // defined in bwa.c, which this test does not link; bwtindex.c refers to it

static long n_fail, n_compared;
#define FAIL(...) do { ++n_fail; fprintf(stderr, "FAIL: " __VA_ARGS__); fputc('\n', stderr); if (n_fail > 20) exit(1); } while (0)

static uint64_t rng = 0x243F6A8885A308D3ULL;
static uint64_t next_rand(void) // splitmix64
{
	uint64_t z = (rng += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

// the original lookup (file sample) and the dense lookup alone, on copies of the index struct
static bwtint_t sa_file(const bwt_t *b, bwtint_t k) { bwt_t o = *b; o.sa32 = 0; return bwt_sa(&o, k); }
static bwtint_t sa_dense(const bwt_t *b, bwtint_t k) { bwt_t d = *b; d.sa = 0; return bwt_sa(&d, k); }

// one LF step, as bwt.c's static bwt_invPsi(), from the public bwt_occ() and bwt_B0()
static bwtint_t lf(const bwt_t *b, bwtint_t k)
{
	bwtint_t x = k - (k > b->primary);
	int c = bwt_B0(b, x);
	x = b->L2[c] + bwt_occ(b, k, c);
	return k == b->primary? 0 : x;
}

static void drop_dense(bwt_t *b) // undo a build that kept the file sample
{
	if (b->sa32) munmap(b->sa32, ((b->seq_len >> b->sa32_shift) + 1) * sizeof(uint32_t));
	b->sa32 = 0;
	b->sa32_shift = 0;
}

/* Build at (shift, threads) keeping the file sample and compare:
 *   - every dense entry with the true suffix-array value (truth, when given) or the original
 *     lookup of its index;
 *   - the dense lookup with the original lookup (and with truth), for every index in [0, seq_len]
 *     when all is set, otherwise for 100,000 random indexes;
 *   - the checking lookup of TRONKO_CHECK_SA (both computed, abort on a difference), for every
 *     16th of those indexes.
 * Returns the number of those lookups whose dense walk ended at index 0 (the $ sentinel). The
 * dense sample is dropped again afterwards. */
static long check_all(bwt_t *b, const char *name, int shift, int threads, const bwtint_t *truth, int all)
{
	bwt_dense_stat_t st;
	bwtint_t i, k, n = b->seq_len + 1, mask = ((bwtint_t)1 << shift) - 1, n_look = all? n : 100000;
	long sentinel = 0, bad = 0;
	int rc = bwt_densify_sa(b, shift, threads, 1, &st);
	if (rc != BWT_DENSE_OK || b->sa32 == 0 || b->sa == 0) {
		FAIL("%s shift %d threads %d: build status %d", name, shift, threads, rc);
		return 0;
	}
	if (st.n_visited != n || st.n_written != st.n_entries || st.n_entries != (b->seq_len >> shift) + 1
			|| st.n_checked != b->n_sa - 1 || st.n_mismatch || st.n_unwritten || st.n_threads != threads)
		FAIL("%s shift %d threads %d: counts visited %llu written %llu entries %llu checked %llu threads %d", name, shift,
				threads, (unsigned long long)st.n_visited, (unsigned long long)st.n_written,
				(unsigned long long)st.n_entries, (unsigned long long)st.n_checked, st.n_threads);
	if (b->sa32[0] != UINT32_MAX) FAIL("%s: entry 0 is %u, not the sentinel", name, b->sa32[0]);
	for (i = 1; i < st.n_entries; ++i) {
		bwtint_t want = truth? truth[i << shift] : sa_file(b, i << shift);
		if (b->sa32[i] != want && ++bad < 5)
			FAIL("%s shift %d threads %d: entry %llu holds %u, SA is %llu", name, shift, threads, (unsigned long long)i,
					b->sa32[i], (unsigned long long)want);
	}
	for (i = 0; i < n_look; ++i) {
		bwtint_t want, d, j;
		k = all? i : next_rand() % n;
		want = sa_file(b, k), d = sa_dense(b, k), j = k;
		if (truth && want != truth[k] && ++bad < 5)
			FAIL("%s: original lookup of %llu is %llu, the suffix array says %llu", name, (unsigned long long)k,
					(unsigned long long)want, (unsigned long long)truth[k]);
		if (d != want && ++bad < 5)
			FAIL("%s shift %d threads %d: k %llu dense %llu original %llu", name, shift, threads,
					(unsigned long long)k, (unsigned long long)d, (unsigned long long)want);
		if (i % 16 == 0 && bwt_sa(b, k) != want && ++bad < 5)
			FAIL("%s: checking lookup of %llu differs", name, (unsigned long long)k);
		while (j & mask) j = lf(b, j);
		sentinel += j == 0;
	}
	if (bad) FAIL("%s shift %d threads %d: %ld differences", name, shift, threads, bad);
	n_compared += st.n_entries - 1 + n_look + (n_look + 15) / 16;
	drop_dense(b);
	return sentinel;
}

/* ---- the repository's BWA index */

static bwt_t *load_fixture(const char *prefix)
{
	char fn[4096];
	bwt_t *b;
	snprintf(fn, sizeof(fn), "%s.bwt", prefix);
	b = bwt_restore_bwt(fn);
	snprintf(fn, sizeof(fn), "%s.sa", prefix);
	bwt_restore_sa(fn, b);
	return b;
}

static void test_fixture(const char *prefix)
{
	static const int threads[] = { 1, 4, 7 };
	bwt_t *b = load_fixture(prefix);
	bwt_dense_stat_t st;
	bwtint_t saved, *sa, k;
	long sentinel = 0;
	int i, rc;
	pid_t pid;

	for (i = 0; i < 3; ++i) sentinel += check_all(b, "fixture", 2, threads[i], 0, i == 1); // every index at 4 threads
	sentinel += check_all(b, "fixture", 1, 3, 0, 0) + check_all(b, "fixture", 3, 5, 0, 0);
	if (b->primary % 4 != 0 && sentinel == 0) FAIL("fixture: primary %% 4 = %d but no walk read the sentinel", (int)(b->primary % 4));
	printf("fixture: seq_len %llu, primary %llu (%% 4 = %d), %ld dense walks through the $ sentinel\n",
			(unsigned long long)b->seq_len, (unsigned long long)b->primary, (int)(b->primary % 4), sentinel);

	// a corrupted file sample: at a segment start (index 1) and beyond the starts (the last index)
	for (i = 0; i < 2; ++i) {
		bwtint_t m = i? b->n_sa - 1 : 1;
		saved = b->sa[m];
		sa = b->sa;
		b->sa[m] ^= 1;
		rc = bwt_densify_sa(b, 2, 4, 0, &st);
		if ((rc != BWT_DENSE_CHECK_FAIL && rc != BWT_DENSE_BAD_STARTS) || b->sa32 != 0 || b->sa != sa)
			FAIL("fixture: corrupted file sample %llu: status %d, sa32 %p", (unsigned long long)m, rc, (void*)b->sa32);
		printf("fixture: file sample %llu corrupted: status %d (mismatches %llu, unwritten %llu), file sample kept\n",
				(unsigned long long)m, rc, (unsigned long long)st.n_mismatch, (unsigned long long)st.n_unwritten);
		b->sa[m] = saved;
	}

	// the skip conditions leave the index as it is
	{
		bwt_t c = *b;
		c.seq_len = UINT32_MAX; // the sentinel could collide: never builds (and touches nothing)
		if (bwt_densify_sa(&c, 2, 2, 0, &st) != BWT_DENSE_SKIPPED || c.sa32) FAIL("fixture: seq_len >= UINT32_MAX not skipped");
		c = *b; c.sa = 0;
		if (bwt_densify_sa(&c, 2, 2, 0, &st) != BWT_DENSE_SKIPPED || c.sa32) FAIL("fixture: no file sample not skipped");
		c = *b;
		if (bwt_densify_sa(&c, 5, 2, 0, &st) != BWT_DENSE_SKIPPED || c.sa32) FAIL("fixture: interval 32 not skipped");
	}

	// TRONKO_CHECK_SA mode aborts on a wrong dense entry (in a child process)
	rc = bwt_densify_sa(b, 2, 3, 1, &st);
	if (rc != BWT_DENSE_OK) FAIL("fixture: rebuild failed");
	for (k = 4; k < b->seq_len && b->sa32[k >> 2] == 0; k += 4) {}
	fflush(stdout);
	if ((pid = fork()) == 0) {
		struct rlimit no_core = { 0, 0 };
		setrlimit(RLIMIT_CORE, &no_core); // the abort is expected; leave no core file
		if (freopen("/dev/null", "w", stderr) == 0) _exit(3);
		b->sa32[k >> 2] ^= 1;
		bwt_sa(b, k);
		_exit(0);
	} else {
		int status = 0;
		waitpid(pid, &status, 0);
		if (!(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT)) FAIL("fixture: checking lookup did not abort on a wrong entry");
		else printf("fixture: TRONKO_CHECK_SA lookup aborts on a corrupted dense entry\n");
	}
	drop_dense(b);

	// normal mode (no TRONKO_CHECK_SA): the file sample is freed and nulled, lookups still equal
	// (against a second copy)
	{
		bwt_t *ref = load_fixture(prefix);
		rc = bwt_densify_sa(b, 2, 4, 0, &st);
		if (rc != BWT_DENSE_OK || b->sa != 0 || b->sa32 == 0) FAIL("fixture: normal build: status %d sa %p", rc, (void*)b->sa);
		for (k = 0; k <= b->seq_len; ++k)
			if (bwt_sa(b, k) != bwt_sa(ref, k)) { FAIL("fixture: normal-mode lookup %llu differs", (unsigned long long)k); break; }
		if (bwt_densify_sa(b, 2, 4, 0, &st) != BWT_DENSE_SKIPPED) FAIL("fixture: a second build was not skipped");
		n_compared += b->seq_len + 1;
		bwt_destroy(ref);
	}
	bwt_destroy(b); // frees sa32 (munmap) and the null sa
}

/* ---- synthetic indexes, suffix array by a plain sort */

static const uint8_t *g_text;
static int64_t g_len;
static int cmp_suffix(const void *pa, const void *pb)
{
	int64_t i = *(const int64_t*)pa, j = *(const int64_t*)pb;
	while (i < g_len && j < g_len) {
		if (g_text[i] != g_text[j]) return g_text[i] < g_text[j]? -1 : 1;
		++i, ++j;
	}
	return i == g_len? (j == g_len? 0 : -1) : 1; // the shorter suffix ($ is smallest) first
}

/* The BWA index of text[0, len) (codes 0-3) with a file sample every intv: the BWT as is_bwt()
 * lays it out ($ removed, primary = the index of the suffix at position 0), then BWA's own
 * bwt_bwtupdate_core() and bwt_cal_sa(). truth[k] = SA[k], with (bwtint_t)-1 for k = 0. */
static bwt_t *make_index(const uint8_t *text, int64_t len, int intv, bwtint_t **truth)
{
	int64_t i, j, *sa = malloc((len + 1) * sizeof(int64_t));
	bwt_t *b = calloc(1, sizeof(bwt_t));
	g_text = text, g_len = len;
	for (i = 0; i <= len; ++i) sa[i] = i;
	qsort(sa, len + 1, sizeof(int64_t), cmp_suffix);
	b->seq_len = len;
	b->bwt_size = (len + 15) >> 4;
	b->bwt = calloc(b->bwt_size + 1, 4);
	for (i = 0; i < len; ++i) ++b->L2[1 + text[i]];
	for (i = 2; i <= 4; ++i) b->L2[i] += b->L2[i-1];
	for (i = 0; i <= len; ++i) if (sa[i] == 0) b->primary = i;
	for (i = j = 0; i <= len; ++i) {
		if (sa[i] == 0) continue;
		b->bwt[j >> 4] |= (uint32_t)text[sa[i] - 1] << ((15 - (j & 15)) << 1);
		++j;
	}
	bwt_bwtupdate_core(b);
	bwt_gen_cnt_table(b);
	bwt_cal_sa(b, intv);
	*truth = malloc((len + 1) * sizeof(bwtint_t));
	for (i = 0; i <= len; ++i) (*truth)[i] = i? (bwtint_t)sa[i] : (bwtint_t)-1;
	for (i = 1; i < (int64_t)b->n_sa; ++i)
		if (b->sa[i] != (bwtint_t)sa[i * intv]) { FAIL("synthetic: bwt_cal_sa sample %lld differs", (long long)i); break; }
	free(sa);
	return b;
}

static void test_synthetic(void)
{
	static const int64_t lens[] = { 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 1000, 4097, 20000, 60001 };
	static const int intvs[] = { 32, 16, 8 };
	int n_lens = sizeof(lens) / sizeof(lens[0]), li, kind, residues = 0, n_idx = 0, big;
	long sentinel = 0;
	for (big = 0; big < 2; ++big)
	for (li = 0; li < (big? 1 : n_lens); ++li)
	for (kind = 0; kind < (big? 1 : 4); ++kind) {
		// big: 270,000 bases with a file sample every 8 has 33,750 file samples, more than twice the
		// 16,384 segments, so the segment starts are every other sample
		int64_t i, len = big? 270000 : lens[li];
		uint8_t *text;
		int ii;
		if ((kind == 3 && len > 1000) || (kind == 2 && len > 20000)) continue; // long repeats make the plain sort quadratic
		text = malloc(len + 1);
		for (i = 0; i < len; ++i) {
			uint64_t r = next_rand();
			text[i] = kind == 0? r & 3 : kind == 1? (r & 1) * 2 : kind == 2? ((i % 7) * 3 + (r % 97 == 0)) & 3 : 1;
		}
		for (ii = 0; ii < 3; ++ii) {
			int intv = intvs[ii], shift, threads;
			bwtint_t *truth;
			bwt_t *b;
			char name[96];
			if (big && intv != 8) continue;
			b = make_index(text, len, intv, &truth);
			++n_idx;
			for (shift = big? 2 : 1; (1 << shift) < intv && (!big || shift == 2); ++shift) {
				threads = 1 + (int)(next_rand() % 9);
				snprintf(name, sizeof(name), "synthetic len %lld kind %d intv %d", (long long)len, kind, intv);
				sentinel += check_all(b, name, shift, threads, truth, 1);
				if (shift == 2) residues |= 1 << (b->primary & 3);
			}
			bwt_destroy(b);
			free(truth);
		}
		free(text);
	}
	printf("synthetic: %d indexes, primary %% 4 residues seen 0x%x, %ld dense walks through the $ sentinel\n", n_idx, residues, sentinel);
	if (residues != 0xf) FAIL("synthetic: not every residue of primary %% 4 occurred (0x%x)", residues);
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		fprintf(stderr, "Usage: %s <prefix of the Charadriiformes BWA index>\n", argv[0]);
		return 2;
	}
	setvbuf(stdout, 0, _IOLBF, 0);
	test_fixture(argv[1]);
	test_synthetic();
	printf("test_dense_sa: %ld entries and lookups compared, %ld failures\n", n_compared, n_fail);
	return n_fail? 1 : 0;
}
