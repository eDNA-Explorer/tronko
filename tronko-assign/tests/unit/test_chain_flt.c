/*
 * Unit test for the chain filter of bwa_source_files/bwamem.c.
 *
 * BWA drops or keeps each chain of seeds in mem_chain_flt() before it extends the kept ones.
 * tronko-assign's copy marks the chains with mem_chain_flt_grouped(), which compares a chain with
 * one representative per group of kept chains sharing (chn_beg, chn_end, is_alt), instead of BWA's
 * loop over every kept chain. Both must set the same .kept and .first on every chain and return
 * the same chains in the same order, or BWA extends different chains and Tronko's output changes.
 *
 * The test compares the current code with a verbatim copy of BWA's original function
 * (base_mem_chain_flt() and base_mark() below, from tronko-assign at 71f6ec3) on:
 *   1. recorded chain sets: the input of mem_chain_flt() calls made by the repository's three
 *      fixture cases (every 40th call and every call with 250 chains or more), with the output the
 *      original binary produced (data/chain_flt_fixture.chf.gz, written by
 *      tools/record_chain_flt.sh); both functions must reproduce the recorded output;
 *   2. random chain sets through the whole function: weights computed from the seeds, the sort,
 *      the marking and the tail that caps and compacts; options away from BWA's defaults, among
 *      them a negative and a NaN drop_ratio (which take the retained pairwise loop), a
 *      min_chain_weight that drops some or all chains, and a small max_chain_extend;
 *   3. random chain lists with the weights set directly, through the marking step alone: ALT
 *      chains, weights of 2^24 and above with dense ties (where the float product rounds),
 *      non-default mask_level, drop_ratio, min_seed_len and max_chain_gap, more than 64 groups
 *      (heap storage) and more than 4,096 chains (heap storage), for mem_chain_flt_grouped() and
 *      for the retained mem_chain_flt_pairwise();
 *   4. constructed lists: a break at the very first kept chain, a group opened after others were
 *      covered, a group that joins a new member after all its members were covered, ALT chains on
 *      either side.
 * Compared per chain: .kept, .first, the kept-index list; after the whole function also the
 * returned count, the order of the surviving chains, .w, and every other field.
 *
 * Usage: test_chain_flt <chain_flt_fixture.chf.gz> [scale]   (scale multiplies parts 2 and 3)
 * Exit status 0 when every comparison agrees.
 */
#include "../../bwa_source_files/bwamem.c"
#include <zlib.h>

/* Defined in tronko-assign.c; BWA's tronko-specific code refers to them. */
queryMatPaired *pairedQueryMat;
queryMatSingle *singleQueryMat;

/* ---- the original, verbatim ------------------------------------------------------------------ */

/* mem_chain_flt() of BWA as vendored in tronko-assign at 71f6ec3, verbatim apart from its name. */
static int base_mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a)
{
	int i, k;
	kvec_t(int) chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
	if (n_chn == 0) return 0; // no need to filter
	// compute the weight of each chain and drop chains with small weight
	for (i = k = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		c->first = -1; c->kept = 0;
		c->w = mem_chain_weight(c);
		if (c->w < opt->min_chain_weight) free(c->seeds);
		else a[k++] = *c;
	}
	n_chn = k;
	ks_introsort(mem_flt, n_chn, a);
	// pairwise chain comparisons
	a[0].kept = 3;
	kv_push(int, chains, 0);
	for (i = 1; i < n_chn; ++i) {
		int large_ovlp = 0;
		for (k = 0; k < chains.n; ++k) {
			int j = chains.a[k];
			int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
			int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
			if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
				int li = chn_end(a[i]) - chn_beg(a[i]);
				int lj = chn_end(a[j]) - chn_beg(a[j]);
				int min_l = li < lj? li : lj;
				if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
					large_ovlp = 1;
					if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
					if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (k == chains.n) {
			kv_push(int, chains, i);
			a[i].kept = large_ovlp? 2 : 3;
		}
	}
	for (i = 0; i < chains.n; ++i) {
		mem_chain_t *c = &a[chains.a[i]];
		if (c->first >= 0) a[c->first].kept = 1;
	}
	free(chains.a);
	for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
		if (a[i].kept == 0 || a[i].kept == 3) continue;
		if (++k >= opt->max_chain_extend) break;
	}
	for (; i < n_chn; ++i)
		if (a[i].kept < 3) a[i].kept = 0;
	for (i = k = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &a[i];
		if (c->kept == 0) free(c->seeds);
		else a[k++] = a[i];
	}
	return k;
}

/* The original marking step alone, verbatim from the function above (from "a[0].kept = 3" to the
 * end of the pairwise loop), for chain lists whose weights the test sets directly. */
static void base_mark(const mem_opt_t *opt, int n_chn, mem_chain_t *a, mem_flt_idx_v *out)
{
	int i, k;
	kvec_t(int) chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
	// pairwise chain comparisons
	a[0].kept = 3;
	kv_push(int, chains, 0);
	for (i = 1; i < n_chn; ++i) {
		int large_ovlp = 0;
		for (k = 0; k < chains.n; ++k) {
			int j = chains.a[k];
			int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
			int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
			if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
				int li = chn_end(a[i]) - chn_beg(a[i]);
				int lj = chn_end(a[j]) - chn_beg(a[j]);
				int min_l = li < lj? li : lj;
				if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
					large_ovlp = 1;
					if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
					if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (k == chains.n) {
			kv_push(int, chains, i);
			a[i].kept = large_ovlp? 2 : 3;
		}
	}
	out->n = chains.n, out->m = chains.m, out->a = chains.a;
}

/* ---- helpers ----------------------------------------------------------------------------------- */

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static long n_cmp, n_bad;

static void fail(const char *part, long iter, const char *what, int i, long x, long y)
{
	if (++n_bad <= 10) fprintf(stderr, "MISMATCH %s call %ld: %s at %d: original %ld, new %ld\n", part, iter, what, i, x, y);
}

/* The marking step as mem_chain_flt() runs it, through the function named by `grouped`. */
static void new_mark(const mem_opt_t *opt, int n_chn, mem_chain_t *a, mem_flt_idx_v *out, int grouped)
{
	out->n = out->m = 0, out->a = 0;
	a[0].kept = 3;
	kv_push(int, *out, 0);
	if (grouped) mem_chain_flt_grouped(opt, n_chn, a, out);
	else mem_chain_flt_pairwise(opt, n_chn, a, out);
}

/* Compare the marks of two copies of one sorted chain list. */
static void cmp_marks(const char *part, long iter, int n, const mem_chain_t *x, const mem_flt_idx_v *cx, const mem_chain_t *y, const mem_flt_idx_v *cy)
{
	int i;
	++n_cmp;
	if (cx->n != cy->n) { fail(part, iter, "number of kept chains", -1, (long)cx->n, (long)cy->n); return; }
	for (i = 0; i < (int)cx->n; ++i)
		if (cx->a[i] != cy->a[i]) { fail(part, iter, "kept index", i, cx->a[i], cy->a[i]); return; }
	for (i = 0; i < n; ++i) {
		if (x[i].kept != y[i].kept) { fail(part, iter, ".kept", i, x[i].kept, y[i].kept); return; }
		if (x[i].first != y[i].first) { fail(part, iter, ".first", i, x[i].first, y[i].first); return; }
	}
}

/* Deep copy of a chain list: every chain gets its own seeds, as mem_chain() hands them over. */
static void copy_chains(int n, const mem_chain_t *src, mem_chain_t *dst)
{
	int i;
	for (i = 0; i < n; ++i) {
		dst[i] = src[i];
		dst[i].seeds = (mem_seed_t*)malloc(src[i].n * sizeof(mem_seed_t));
		memcpy(dst[i].seeds, src[i].seeds, src[i].n * sizeof(mem_seed_t));
	}
}

/* Compare the output of the whole function: count, order (by .pos, which the test sets to the
 * input index), and every field of every surviving chain. */
static void cmp_output(const char *part, long iter, int kx, const mem_chain_t *x, int ky, const mem_chain_t *y)
{
	int i;
	++n_cmp;
	if (kx != ky) { fail(part, iter, "returned count", -1, kx, ky); return; }
	for (i = 0; i < kx; ++i) {
		if (x[i].pos != y[i].pos) { fail(part, iter, "chain at output position", i, (long)x[i].pos, (long)y[i].pos); return; }
		if (x[i].w != y[i].w) { fail(part, iter, ".w", i, x[i].w, y[i].w); return; }
		if (x[i].kept != y[i].kept) { fail(part, iter, ".kept", i, x[i].kept, y[i].kept); return; }
		if (x[i].first != y[i].first) { fail(part, iter, ".first", i, x[i].first, y[i].first); return; }
		if (x[i].n != y[i].n || x[i].rid != y[i].rid || x[i].is_alt != y[i].is_alt || memcmp(&x[i].frac_rep, &y[i].frac_rep, sizeof(float)) != 0
				|| memcmp(x[i].seeds, y[i].seeds, x[i].n * sizeof(mem_seed_t)) != 0) {
			fail(part, iter, "other fields", i, 0, 1); return;
		}
	}
}

static void free_seeds(int n, mem_chain_t *a)
{
	int i;
	for (i = 0; i < n; ++i) free(a[i].seeds);
}

/* Run both whole functions on deep copies of `in` and compare; return the original's count and
 * leave its output in *xo (seeds freed) for the caller's own checks. */
static int run_both(const char *part, long iter, const mem_opt_t *opt, int n, const mem_chain_t *in, mem_chain_t *xo)
{
	mem_chain_t *x = (mem_chain_t*)malloc(n * sizeof(mem_chain_t)), *y = (mem_chain_t*)malloc(n * sizeof(mem_chain_t));
	int kx, ky;
	copy_chains(n, in, x);
	copy_chains(n, in, y);
	kx = base_mem_chain_flt(opt, n, x);
	ky = mem_chain_flt(opt, n, y);
	cmp_output(part, iter, kx, x, ky, y);
	if (xo) memcpy(xo, x, kx * sizeof(mem_chain_t));
	free_seeds(kx, x);
	free_seeds(ky, y);
	free(x); free(y);
	return kx;
}

/* ---- part 1: recorded chain sets ----------------------------------------------------------------
 * File: gzip of records. Record: int32 marker 0x43484631 ("CHF1"); options int32 min_chain_weight,
 * float drop_ratio, float mask_level, int32 max_chain_gap, int32 min_seed_len, int32
 * max_chain_extend; int32 n_chn; per chain int32 n, int32 rid, int32 is_alt, float frac_rep,
 * int64 pos, then n seeds of int64 rbeg, int32 qbeg, int32 len, int32 score; int32 k (the count
 * returned); per surviving chain int32 input index, uint32 w, int32 kept, int32 first. */

static int gz_get(gzFile fp, void *p, unsigned len) { return gzread(fp, p, len) == (int)len; }

static long part1(const char *path, long *n_chains_total)
{
	gzFile fp = gzopen(path, "rb");
	long n_rec = 0;
	int32_t marker;
	if (fp == 0) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
	while (gz_get(fp, &marker, 4)) {
		mem_opt_t *opt = mem_opt_init();
		int32_t n, k, i, v[6], rec_k;
		mem_chain_t *in, *xo;
		if (marker != 0x43484631) { fprintf(stderr, "%s: bad record %ld\n", path, n_rec); exit(2); }
		if (!gz_get(fp, v, sizeof(v)) || !gz_get(fp, &n, 4) || n <= 0) { fprintf(stderr, "%s: truncated record %ld\n", path, n_rec); exit(2); }
		opt->min_chain_weight = v[0]; memcpy(&opt->drop_ratio, &v[1], 4); memcpy(&opt->mask_level, &v[2], 4);
		opt->max_chain_gap = v[3]; opt->min_seed_len = v[4]; opt->max_chain_extend = v[5];
		in = (mem_chain_t*)calloc(n, sizeof(mem_chain_t));
		xo = (mem_chain_t*)calloc(n, sizeof(mem_chain_t));
		for (i = 0; i < n; ++i) {
			int32_t h[3], s;
			float fr;
			int64_t pos;
			if (!gz_get(fp, h, 12) || !gz_get(fp, &fr, 4) || !gz_get(fp, &pos, 8) || h[0] <= 0) { fprintf(stderr, "%s: truncated record %ld\n", path, n_rec); exit(2); }
			in[i].n = in[i].m = h[0], in[i].rid = h[1], in[i].is_alt = h[2], in[i].frac_rep = fr;
			in[i].pos = i; // identifies the chain in the output; mem_chain_flt() never reads .pos
			in[i].seeds = (mem_seed_t*)malloc(h[0] * sizeof(mem_seed_t));
			for (s = 0; s < h[0]; ++s) {
				mem_seed_t *p = &in[i].seeds[s];
				int32_t q[3];
				if (!gz_get(fp, &p->rbeg, 8) || !gz_get(fp, q, 12)) { fprintf(stderr, "%s: truncated record %ld\n", path, n_rec); exit(2); }
				p->qbeg = q[0], p->len = q[1], p->score = q[2];
			}
			(void)pos;
		}
		k = run_both("recorded", n_rec, opt, n, in, xo);
		if (!gz_get(fp, &rec_k, 4)) { fprintf(stderr, "%s: truncated record %ld\n", path, n_rec); exit(2); }
		++n_cmp;
		if (rec_k != k) fail("recorded (against the recording)", n_rec, "returned count", -1, rec_k, k);
		for (i = 0; i < rec_k; ++i) {
			int32_t r[4];
			if (!gz_get(fp, r, 16)) { fprintf(stderr, "%s: truncated record %ld\n", path, n_rec); exit(2); }
			if (i < k && (r[0] != xo[i].pos || (uint32_t)r[1] != xo[i].w || r[2] != xo[i].kept || r[3] != xo[i].first))
				fail("recorded (against the recording)", n_rec, "surviving chain", i, r[0], (long)xo[i].pos);
		}
		free_seeds(n, in);
		free(in); free(xo); free(opt);
		*n_chains_total += n;
		++n_rec;
	}
	gzclose(fp);
	return n_rec;
}

/* ---- random option sets -------------------------------------------------------------------------- */

static const float masks[] = { 0.5f, 0.0f, 0.3f, 0.95f, 1.0f };
static const float drops[] = { 0.5f, 0.0f, 0.99f, 1.0f, 0.3f };

static void random_marking_opt(mem_opt_t *opt)
{
	opt->min_seed_len = rnd() % 4 == 0? (int)(rnd() % 6) : 19;
	opt->max_chain_gap = rnd() % 6 == 0? 5 + (int)(rnd() % 60) : 10000;
	opt->mask_level = masks[rnd() % 5];
	opt->drop_ratio = drops[rnd() % 5];
}

/* ---- part 2: random chain sets through the whole function ---------------------------------------- */

static long part2(long iters)
{
	static mem_chain_t in[400];
	static mem_seed_t pool[400][4];
	long it;
	for (it = 0; it < iters; ++it) {
		mem_opt_t *opt = mem_opt_init();
		int i, s, n = 1 + (int)(rnd() % (rnd() % 10 == 0? 400 : 30));
		int big = rnd() % 10 == 0, alt = rnd() % 3 == 0, nb = 2 + (int)(rnd() % 6);
		random_marking_opt(opt);
		switch (rnd() % 8) { // include drop ratios that send mem_chain_flt() to the pairwise loop
			case 0: opt->drop_ratio = -0.25f; break;
			case 1: opt->drop_ratio = NAN; break;
			case 2: opt->drop_ratio = INFINITY; break;
			default: break;
		}
		if (rnd() % 8 == 0) opt->max_chain_extend = 1 + (int)(rnd() % 4);
		for (i = 0; i < n; ++i) {
			int ns = 1 + (int)(rnd() % 4), q, r;
			mem_chain_t *c = &in[i];
			c->n = c->m = ns, c->seeds = pool[i], c->rid = (int)(rnd() % 5), c->is_alt = alt? rnd() % 2 : 0;
			c->frac_rep = (float)(rnd() % 100) / 100.0f, c->pos = i, c->first = 7, c->kept = 1, c->w = 3;
			if (big) q = (int)(rnd() % nb) << 22, r = (int)(rnd() % nb) << 22;
			else q = (int)(rnd() % nb) * 17, r = (int)(rnd() % nb) * 23;
			for (s = 0; s < ns; ++s) {
				mem_seed_t *p = &c->seeds[s];
				p->qbeg = q, p->rbeg = r, p->score = (int)(rnd() % 100);
				p->len = big? (1 << 24) + (int)(rnd() % 64) : 19 + (int)(rnd() % nb) * 13;
				if (big) q += (int)(rnd() % nb) << 20, r += (int)(rnd() % 3) << 20;
				else q += (int)(rnd() % nb) * 11, r += (int)(rnd() % 3) * 11;
			}
		}
		if (rnd() % 4 == 0) { // drop some chains by weight; now and then all of them
			int w0 = mem_chain_weight(&in[rnd() % n]);
			opt->min_chain_weight = rnd() % 8 == 0? INT_MAX : w0;
		}
		run_both("random whole-function", it, opt, n, in, 0);
		free(opt);
	}
	return iters;
}

/* ---- part 3: random chain lists with direct weights, through the marking step ------------------- */

static void part3_one(const char *part, long it, const mem_opt_t *opt, int n, mem_chain_t *a, mem_chain_t *b, mem_chain_t *c)
{
	mem_flt_idx_v ca, cb, cc;
	int i;
	for (i = 0; i < n; ++i) a[i].first = -1, a[i].kept = 0;
	ks_introsort(mem_flt, n, a);
	memcpy(b, a, n * sizeof(mem_chain_t)); // the marking step reads the seeds and never frees them
	memcpy(c, a, n * sizeof(mem_chain_t));
	base_mark(opt, n, a, &ca);
	if (opt->drop_ratio >= 0.0f) {
		new_mark(opt, n, b, &cb, 1);
		cmp_marks(part, it, n, a, &ca, b, &cb);
		free(cb.a);
	}
	new_mark(opt, n, c, &cc, 0);
	cmp_marks(part, it, n, a, &ca, c, &cc);
	free(ca.a); free(cc.a);
}

static long part3(long iters)
{
	enum { NMAX = 6000 };
	static mem_chain_t a[NMAX], b[NMAX], c[NMAX];
	static mem_seed_t pool[NMAX][3];
	long it;
	for (it = 0; it < iters; ++it) {
		mem_opt_t *opt = mem_opt_init();
		int i, s, n, big = rnd() % 10 == 0, alt = rnd() % 3 == 0, nb = 2 + (int)(rnd() % 6), wide = 0;
		random_marking_opt(opt);
		if (it % 25000 == 24999) n = 4097 + (int)(rnd() % (NMAX - 4097)), wide = 1; // heap storage for next[]
		else if (it % 100 == 99) n = 200 + (int)(rnd() % 400), wide = 1; // more than 64 groups
		else n = 1 + (int)(rnd() % (rnd() % 10 == 0? 400 : 40));
		for (i = 0; i < n; ++i) {
			int ns = 1 + (int)(rnd() % 3), q = wide? (int)(rnd() % 4096) * 7 : (int)(rnd() % nb) * 17;
			uint32_t w;
			a[i].n = a[i].m = ns, a[i].seeds = pool[i], a[i].rid = i, a[i].is_alt = alt? rnd() % 2 : 0, a[i].pos = i;
			for (s = 0; s < ns; ++s) {
				pool[i][s].qbeg = q, pool[i][s].rbeg = 0, pool[i][s].score = 0;
				pool[i][s].len = 19 + (int)(rnd() % nb) * 13;
				q += (int)(rnd() % nb) * 11;
			}
			w = big? (uint32_t)((1u << 24) + rnd() % ((1u << 29) - (1u << 24))) : (uint32_t)(19 + rnd() % (rnd() % 2? 8 : 200));
			if (big && rnd() % 2) w = (uint32_t)((1u << 24) + rnd() % 64); // dense ties just above 2^24, where the product rounds
			a[i].w = w;
		}
		part3_one("random marking", it, opt, n, a, b, c);
		free(opt);
	}
	return iters;
}

/* ---- part 4: constructed lists ------------------------------------------------------------------
 * Each case lists chains as {qbeg, qend, weight, is_alt}, one seed per chain, heaviest first. */

typedef struct { int qb, qe, w, alt; } spec_t;

static void part4_case(const char *name, int n, const spec_t *sp, float drop_ratio)
{
	mem_chain_t a[16], b[16], c[16];
	mem_seed_t seeds[16];
	mem_opt_t *opt = mem_opt_init();
	int i;
	opt->drop_ratio = drop_ratio;
	for (i = 0; i < n; ++i) {
		seeds[i].qbeg = sp[i].qb, seeds[i].len = sp[i].qe - sp[i].qb, seeds[i].rbeg = 100 * i, seeds[i].score = 0;
		a[i].n = a[i].m = 1, a[i].seeds = &seeds[i], a[i].rid = 0, a[i].is_alt = sp[i].alt, a[i].pos = i, a[i].w = sp[i].w;
	}
	part3_one(name, 0, opt, n, a, b, c);
	free(opt);
}

static long part4(void)
{
	/* chain 1 breaks at chain 0, the very first kept chain */
	static const spec_t first_break[] = { {0, 100, 100, 0}, {0, 100, 20, 0}, {10, 90, 19, 0} };
	/* groups A=[0,100) and B=[200,300) are covered by a later chain that overlaps both; group C,
	 * opened after that, is then covered by its own later chains */
	static const spec_t late_group[] = { {0, 100, 90, 0}, {200, 300, 88, 0}, {0, 300, 80, 0}, {0, 100, 70, 0},
		{400, 500, 60, 0}, {400, 500, 59, 0}, {0, 500, 30, 0}, {400, 500, 10, 0} };
	/* a group whose members are all covered gains a new member, which must be covered next time */
	static const spec_t refill[] = { {0, 100, 100, 0}, {0, 100, 99, 0}, {0, 100, 98, 0}, {0, 100, 97, 0}, {0, 100, 10, 0} };
	/* ALT chains: a primary chain is not shadowed by a kept ALT chain; an ALT chain is */
	static const spec_t alt[] = { {0, 100, 100, 1}, {0, 100, 50, 0}, {0, 100, 40, 1}, {0, 100, 10, 0}, {0, 100, 5, 1} };
	part4_case("constructed first-break", 3, first_break, 0.5f);
	part4_case("constructed late-group", 8, late_group, 0.5f);
	part4_case("constructed late-group, drop 0.9", 8, late_group, 0.9f);
	part4_case("constructed refill", 5, refill, 0.5f);
	part4_case("constructed alt", 5, alt, 0.5f);
	part4_case("constructed alt, drop 1.0", 5, alt, 1.0f);
	return 6;
}

int main(int argc, char *argv[])
{
	long n_rec, n_chains = 0, scale = 1, n2, n3, n4;
	if (argc < 2) { fprintf(stderr, "Usage: %s <chain_flt_fixture.chf.gz> [scale]\n", argv[0]); return 2; }
	if (argc > 2) scale = atol(argv[2]);
	if (scale < 1) scale = 1;
	n_rec = part1(argv[1], &n_chains);
	if (n_rec == 0) { fprintf(stderr, "%s: no records\n", argv[1]); return 2; }
	n2 = part2(20000 * scale);
	n3 = part3(50000 * scale);
	n4 = part4();
	printf("test_chain_flt: %ld recorded calls (%ld chains), %ld random whole-function calls, %ld random marking calls, %ld constructed; %ld comparisons, %ld differ\n",
		n_rec, n_chains, n2, n3, n4, n_cmp, n_bad);
	return n_bad == 0? 0 : 1;
}
