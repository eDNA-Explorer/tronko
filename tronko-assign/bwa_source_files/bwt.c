/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <pthread.h>
#include <sched.h>
#include <sys/mman.h>
#include <time.h>
#include <alloca.h>
#include "utils.h"
#include "bwt.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	bwtint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
	bwtint_t isa, sa, i; // S(isa) = sa
	int intv_round = intv;

	kv_roundup32(intv_round);
	xassert(intv_round == intv, "SA sample interval is not a power of 2.");
	xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

	if (bwt->sa) free(bwt->sa);
	bwt->sa_intv = intv;
	bwt->n_sa = (bwt->seq_len + intv) / intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	// calculate SA value
	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % intv == 0) bwt->sa[isa/intv] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % intv == 0) bwt->sa[isa/intv] = sa;
	bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}

static inline bwtint_t bwt_sa_sampled(const bwt_t *bwt, bwtint_t k) // the original lookup, on the file sample
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}

bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k)
{
	if (bwt->sa32) {
		// Same LF walk as bwt_sa_sampled(), stopped at the first multiple of 1 << sa32_shift, which
		// comes at or before the first multiple of sa_intv. Both samples hold true SA values, so
		// the result is the same; index 0 ($) holds UINT32_MAX for the file sample's (bwtint_t)-1.
		bwtint_t sa = 0, k0 = k, mask = ((bwtint_t)1 << bwt->sa32_shift) - 1;
		uint32_t v;
		while (k & mask) {
			++sa;
			k = bwt_invPsi(bwt, k);
		}
		v = bwt->sa32[k >> bwt->sa32_shift];
		sa += v == UINT32_MAX? (bwtint_t)-1 : (bwtint_t)v;
		if (bwt->sa && sa != bwt_sa_sampled(bwt, k0)) { // TRONKO_CHECK_SA: the file sample was kept
			fprintf(stderr, "[E::%s] dense SA lookup of %llu gave %llu, the file sample %llu\n", __func__,
					(unsigned long long)k0, (unsigned long long)sa, (unsigned long long)bwt_sa_sampled(bwt, k0));
			abort();
		}
		return sa;
	}
	return bwt_sa_sampled(bwt, k);
}

/*********************************************************************
 * Dense suffix-array sample (tronko). bwt_densify_sa() writes SA[k] for every k that is a
 * multiple of 1 << shift into a uint32 array, by walking the text backwards from starting points
 * taken from the file sample: (isa, pos) = (m * sa_intv, sa[m]) for evenly spaced m >= 1, plus
 * (0, seq_len) for the $ suffix. Sorted by pos, start i owns the text positions from its pos down
 * to the previous start's pos + 1 (the lowest start down to 0), so the segments cover
 * [0, seq_len] once, and since positions and SA indices correspond one to one, every index is
 * visited once. A step is isa <- LF(isa), pos <- pos - 1, the walk of bwt_cal_sa().
 *
 * Self-check: the array is filled with UINT32_MAX first (no true value reaches it while
 * seq_len < UINT32_MAX) and scanned after the walks; every nonzero file sample is compared with
 * the walk's position; positions visited, entries written and samples compared are counted. Any
 * failure discards the dense sample.
 *
 * Nothing here allocates from the malloc heap: the arrays are mmap'd and the thread handles are
 * on the stack, so the main thread's heap history (which BWA's mate rescue reads past the end
 * of mem_opt_t) is unchanged. Threads claim work from shared counters; the caller runs as one of
 * the workers, so the build completes with however many threads could be created.
 *********************************************************************/

#define DSA_WALKS     16384        // text segments (plus the $ one)
#define DSA_IN_FLIGHT 8            // segments a thread advances round-robin, for memory-level parallelism
#define DSA_SLICE     (1ULL << 22) // entries per fill or scan work unit (16 MiB)

typedef struct { bwtint_t isa, pos, lo; } dsa_walk_t;

typedef struct {
	const bwt_t *bwt;
	uint32_t *sa32;
	int shift;
	const dsa_walk_t *walk;
	int64_t n_walk, n_slice;
	uint64_t n_ent;
	int64_t fill_next, fill_done, walk_next, walk_done, scan_next; // work counters, atomic
	uint64_t visited, written, checked, mismatch, unwritten;       // totals, atomic
} dsa_job_t;

static inline void dsa_prefetch(const bwt_t *bwt, bwtint_t k) // the occurrence block bwt_invPsi(bwt, k) reads
{
	const uint32_t *p = bwt_occ_intv(bwt, k - (k >= bwt->primary));
	__builtin_prefetch(p);
	__builtin_prefetch(p + 15);
}

static void dsa_wait(int64_t *counter, int64_t target)
{
	while (__atomic_load_n(counter, __ATOMIC_ACQUIRE) < target) sched_yield();
}

static void *dsa_worker(void *data)
{
	dsa_job_t *J = (dsa_job_t*)data;
	const bwt_t *bwt = J->bwt;
	const int shift = J->shift;
	const bwtint_t mask = ((bwtint_t)1 << shift) - 1, file_mask = (bwtint_t)bwt->sa_intv - 1;
	uint64_t visited = 0, written = 0, checked = 0, mismatch = 0, unwritten = 0;
	dsa_walk_t act[DSA_IN_FLIGHT];
	int na = 0, g, more = 1;
	int64_t i;

	// 1. fill with UINT32_MAX, so an entry no walk writes is visible afterwards
	while ((i = __atomic_fetch_add(&J->fill_next, 1, __ATOMIC_RELAXED)) < J->n_slice) {
		uint64_t beg = (uint64_t)i * DSA_SLICE, end = beg + DSA_SLICE < J->n_ent? beg + DSA_SLICE : J->n_ent;
		memset(J->sa32 + beg, 0xff, (end - beg) * sizeof(uint32_t));
		__atomic_fetch_add(&J->fill_done, 1, __ATOMIC_RELEASE);
	}
	dsa_wait(&J->fill_done, J->n_slice);

	// 2. walk the segments, DSA_IN_FLIGHT at a time
	for (;;) {
		while (more && na < DSA_IN_FLIGHT) {
			i = __atomic_fetch_add(&J->walk_next, 1, __ATOMIC_RELAXED);
			if (i >= J->n_walk) { more = 0; break; }
			act[na++] = J->walk[i];
			dsa_prefetch(bwt, J->walk[i].isa);
		}
		if (na == 0) break;
		for (g = 0; g < na;) {
			dsa_walk_t *x = &act[g];
			++visited;
			if (!(x->isa & mask)) {
				J->sa32[x->isa >> shift] = (uint32_t)x->pos;
				++written;
			}
			if (!(x->isa & file_mask) && x->isa != 0) { // index 0 holds -1 in the file, not seq_len
				++checked;
				if (bwt->sa[x->isa / bwt->sa_intv] != x->pos) ++mismatch;
			}
			if (x->pos == x->lo) { // segment done; tested before stepping, so pos never passes 0
				act[g] = act[--na];
				__atomic_fetch_add(&J->walk_done, 1, __ATOMIC_RELEASE);
				continue;
			}
			x->isa = bwt_invPsi(bwt, x->isa);
			--x->pos;
			dsa_prefetch(bwt, x->isa);
			++g;
		}
	}
	dsa_wait(&J->walk_done, J->n_walk);

	// 3. scan for entries no walk wrote (entry 0 is set by the caller)
	while ((i = __atomic_fetch_add(&J->scan_next, 1, __ATOMIC_RELAXED)) < J->n_slice) {
		uint64_t j, beg = (uint64_t)i * DSA_SLICE, end = beg + DSA_SLICE < J->n_ent? beg + DSA_SLICE : J->n_ent;
		for (j = beg? beg : 1; j < end; ++j)
			unwritten += J->sa32[j] == UINT32_MAX;
	}

	__atomic_fetch_add(&J->visited, visited, __ATOMIC_RELAXED);
	__atomic_fetch_add(&J->written, written, __ATOMIC_RELAXED);
	__atomic_fetch_add(&J->checked, checked, __ATOMIC_RELAXED);
	__atomic_fetch_add(&J->mismatch, mismatch, __ATOMIC_RELAXED);
	__atomic_fetch_add(&J->unwritten, unwritten, __ATOMIC_RELAXED);
	return 0;
}

static void dsa_sort_by_pos(dsa_walk_t *a, int64_t n) // heapsort: qsort may call malloc
{
	int64_t start = n / 2, end = n, root, child;
	dsa_walk_t t;
	for (;;) {
		if (start > 0) --start;
		else if (--end > 0) { t = a[0]; a[0] = a[end]; a[end] = t; }
		else break;
		for (root = start; (child = 2 * root + 1) < end; root = child) {
			if (child + 1 < end && a[child + 1].pos > a[child].pos) ++child;
			if (a[root].pos >= a[child].pos) break;
			t = a[root]; a[root] = a[child]; a[child] = t;
		}
	}
}

static double dsa_now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + t.tv_nsec * 1e-9;
}

int bwt_densify_sa(bwt_t *bwt, int shift, int n_threads, int keep_sa, bwt_dense_stat_t *st)
{
	double t0 = dsa_now();
	dsa_job_t J;
	dsa_walk_t *walk;
	pthread_t *tid;
	size_t walk_bytes, ent_bytes;
	int64_t i, n_file, n_walk, stride;
	int t, n_started = 0;

	memset(st, 0, sizeof(*st));
	st->status = BWT_DENSE_SKIPPED;
	if (bwt->sa == 0 || bwt->sa32 != 0 || shift < 1 || shift > 30 || bwt->sa_intv <= (1 << shift)
			|| bwt->seq_len >= UINT32_MAX)
		return st->status;
	if (n_threads < 1) n_threads = 1;
	if (n_threads > 1024) n_threads = 1024;

	memset(&J, 0, sizeof(J));
	J.bwt = bwt;
	J.shift = shift;
	J.n_ent = (bwt->seq_len >> shift) + 1;
	J.n_slice = (int64_t)((J.n_ent + DSA_SLICE - 1) / DSA_SLICE);
	n_file = (int64_t)bwt->n_sa - 1; // nonzero file samples
	n_walk = n_file < DSA_WALKS? n_file : DSA_WALKS;
	J.n_walk = n_walk + 1;
	st->n_entries = J.n_ent;
	st->n_walks = J.n_walk;

	ent_bytes = J.n_ent * sizeof(uint32_t);
	walk_bytes = J.n_walk * sizeof(dsa_walk_t);
	J.sa32 = (uint32_t*)mmap(0, ent_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	walk = (dsa_walk_t*)mmap(0, walk_bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (J.sa32 == MAP_FAILED || walk == MAP_FAILED) {
		if (J.sa32 != MAP_FAILED) munmap(J.sa32, ent_bytes);
		if (walk != MAP_FAILED) munmap(walk, walk_bytes);
		st->status = BWT_DENSE_NO_MEMORY;
		st->seconds = dsa_now() - t0;
		return st->status;
	}

	// starting points, sorted by text position; the $ suffix (pos seq_len) sorts last
	stride = n_walk? n_file / n_walk : 0;
	for (i = 0; i < n_walk; ++i) {
		bwtint_t m = 1 + (bwtint_t)i * stride;
		walk[i].isa = m * bwt->sa_intv;
		walk[i].pos = bwt->sa[m];
	}
	walk[n_walk].isa = 0;
	walk[n_walk].pos = bwt->seq_len;
	dsa_sort_by_pos(walk, J.n_walk);
	for (i = 0; i < J.n_walk; ++i) {
		if (i > 0 && walk[i].pos <= walk[i-1].pos) break; // a corrupt sample: segments would overlap
		walk[i].lo = i? walk[i-1].pos + 1 : 0;
	}
	if (i < J.n_walk || walk[J.n_walk - 1].isa != 0) {
		munmap(J.sa32, ent_bytes); munmap(walk, walk_bytes);
		st->status = BWT_DENSE_BAD_STARTS;
		st->seconds = dsa_now() - t0;
		return st->status;
	}
	J.walk = walk;

	tid = (pthread_t*)alloca(n_threads * sizeof(pthread_t));
	for (t = 1; t < n_threads; ++t)
		if (pthread_create(&tid[n_started], 0, dsa_worker, &J) == 0) ++n_started;
	dsa_worker(&J);
	for (t = 0; t < n_started; ++t) pthread_join(tid[t], 0);
	munmap(walk, walk_bytes);
	J.sa32[0] = UINT32_MAX; // the $ suffix; its walk wrote seq_len there

	st->n_threads = n_started + 1;
	st->n_visited = J.visited;
	st->n_written = J.written;
	st->n_checked = J.checked;
	st->n_mismatch = J.mismatch;
	st->n_unwritten = J.unwritten;
	if (J.mismatch != 0 || J.unwritten != 0 || J.visited != bwt->seq_len + 1 || J.written != J.n_ent
			|| J.checked != (uint64_t)n_file) {
		munmap(J.sa32, ent_bytes);
		st->status = BWT_DENSE_CHECK_FAIL;
		st->seconds = dsa_now() - t0;
		return st->status;
	}
	bwt->sa32 = J.sa32;
	bwt->sa32_shift = shift;
	if (!keep_sa) {
		free(bwt->sa);
		bwt->sa = 0;
	}
	st->status = BWT_DENSE_OK;
	st->seconds = dsa_now() - t0;
	return st->status;
}

static inline int __occ_aux(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
	bwtint_t _k, _l;
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		*ok = bwt_occ(bwt, k, c);
		*ol = bwt_occ(bwt, l, c);
	} else {
		bwtint_t m, n, i, j;
		uint32_t *p;
		if (k >= bwt->primary) --k;
		if (l >= bwt->primary) --l;
		n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
		p += sizeof(bwtint_t);
		// calculate *ok
		j = k >> 5 << 5;
		for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
			n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m = n;
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0) n -= ~k&31; // corrected for the masked bits
		*ok = n;
		// calculate *ol
		j = l >> 5 << 5;
		for (; i < j; i += 32, p += 2)
			m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
		m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
		if (c == 0) m -= ~l&31; // corrected for the masked bits
		*ol = m;
	}
}

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t x;
	uint32_t *p, tmp, *end;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(bwtint_t));
	p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
	end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p) x += __occ_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += __occ_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if (_l>>OCC_INTV_SHIFT != _k>>OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
		for (x = 0; p < endk; ++p) x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
	bwtint_t k, l, ok, ol;
	int i;
	k = 0; l = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) break; // no match
	}
	if (k > l) return 0; // no match
	if (sa_begin) *sa_begin = k;
	if (sa_end)   *sa_end = l;
	return l - k + 1;
}

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
	int i;
	bwtint_t k, l, ok, ol;
	k = *k0; l = *l0;
	for (i = len - 1; i >= 0; --i) {
		ubyte_t c = str[i];
		if (c > 3) return 0; // there is an N here. no match
		bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
		k = bwt->L2[c] + ok + 1;
		l = bwt->L2[c] + ol;
		if (k > l) return 0; // no match
	}
	*k0 = k; *l0 = l;
	return l - k + 1;
}

/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

static void bwt_reverse_intvs(bwtintv_v *p)
{
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}
// NOTE: $max_intv is not currently used in BWA-MEM
int bwt_smem1a(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			if (c >= 0 && ik.x[2] >= max_intv) bwt_extend(bwt, p, ok, 1);
			if (c < 0 || ik.x[2] < max_intv || ok[c].x[2] < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}
	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int bwt_smem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2])
{
	return bwt_smem1a(bwt, len, q, x, min_intv, 0, mem, tmpvec);
}

int bwt_seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem)
{
	int i, c;
	bwtintv_t ik, ok[4];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			if (ok[c].x[2] < max_intv && i - x >= min_len) {
				*mem = ok[c];
				mem->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	return len;
}

/*************************
 * Read/write BWT and SA *
 *************************/

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	err_fflush(fp);
	err_fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fflush(fp);
	err_fclose(fp);
}

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
	err_fread_noeof(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	err_fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size<<2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	err_fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	if (bwt->sa32) munmap(bwt->sa32, ((bwt->seq_len >> bwt->sa32_shift) + 1) * sizeof(uint32_t));
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}
