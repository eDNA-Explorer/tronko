/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

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

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#ifdef __aarch64__
#include "sse2neon.h"
#else
#include <emmintrin.h>
#include <immintrin.h>
#endif /* ifdef __aarch64__ */
#include "ksw.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

const kswr_t g_defr = { 0, -1, -1, -1, -1, -1, -1 };

struct _kswq_t {
	int qlen, slen;
	uint8_t shift, mdiff, max, size;
	__m128i *qp, *H0, *H1, *E, *Hmax;
};

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
	kswq_t *q;
	int slen, a, tmp, p;

	size = size > 1? 2 : 1;
	p = 8 * (3 - size); // # values per __m128i
	slen = (qlen + p - 1) / p; // segmented length
	q = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
	q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
	q->H0 = q->qp + slen * m;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen; q->qlen = qlen; q->size = size;
	// compute shift
	tmp = m * m;
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
		if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
		if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
	}
	q->max = q->mdiff;
	q->shift = 256 - q->shift; // NB: q->shift is uint8_t
	q->mdiff += q->shift; // this is the difference between the min and max scores
	// An example: p=8, qlen=19, slen=3 and segmentation:
	//  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
	if (size == 1) {
		int8_t *t = (int8_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
		}
	} else {
		int16_t *t = (int16_t*)q->qp;
		for (a = 0; a < m; ++a) {
			int i, k, nlen = slen * p;
			const int8_t *ma = mat + a * m;
			for (i = 0; i < slen; ++i)
				for (k = i; k < nlen; k += slen) // p iterations
					*t++ = (k >= qlen? 0 : ma[query[k]]);
		}
	}
	return q;
}

kswr_t ksw_u8(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
	kswr_t r;

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
    	(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi8(_o_del + _e_del);
	e_del = _mm_set1_epi8(_e_del);
	oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
	e_ins = _mm_set1_epi8(_e_ins);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, cmp, imax;
		__m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
		for (j = 0; LIKELY(j < slen); ++j) {
			/* SW cells are computed in the following order:
			 *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */
			// compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));
			h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
			e = _mm_load_si128(E + j); // e=E'(i,j)
			h = _mm_max_epu8(h, e);
			h = _mm_max_epu8(h, f); // h=H'(i,j)
			max = _mm_max_epu8(max, h); // set max
			_mm_store_si128(H1 + j, h); // save to H'(i,j)
			// now compute E'(i+1,j)
			e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
			t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
			e = _mm_max_epu8(e, t); // e=E'(i+1,j)
			_mm_store_si128(E + j, e); // save to E'(i+1,j)
			// now compute F'(i,j+1)
			f = _mm_subs_epu8(f, e_ins);
			t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
			f = _mm_max_epu8(f, t);
			// get H'(i-1,j) and prepare for the next j
			h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
		}
		// NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
		for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epu8(h, f); // h=H'(i,j)
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, oe_ins);
				f = _mm_subs_epu8(f, e_ins);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
			}
		}
end_loop16:
		//int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
		__max_16(imax, max); // imax is the maximum number in max
		if (imax >= minsc) { // write the b array; this condition adds branching unfornately
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i; // te is the end position on the target
			for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax + q->shift >= 255 || gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S; // swap H0 and H1
	}
	r.score = gmax + q->shift < 255? gmax : 255;
	r.te = te;
	if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
		int max = -1, tmp, low, high, qlen = slen * 16;
		uint8_t *t = (uint8_t*)Hmax;
		for (i = 0; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
			else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
		//printf("%d,%d\n", max, gmax);
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

kswr_t ksw_i16(kswq_t *q, int tlen, const uint8_t *target, int _o_del, int _e_del, int _o_ins, int _e_ins, int xtra) // the first gap costs -(_o+_e)
{
	int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
	uint64_t *b;
	__m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
	kswr_t r;

#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    	(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)

	// initialization
	r = g_defr;
	minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
	endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
	m_b = n_b = 0; b = 0;
	zero = _mm_set1_epi32(0);
	oe_del = _mm_set1_epi16(_o_del + _e_del);
	e_del = _mm_set1_epi16(_e_del);
	oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
	e_ins = _mm_set1_epi16(_e_ins);
	H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; ++i) {
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}
	// the core loop
	for (i = 0; i < tlen; ++i) {
		int j, k, imax;
		__m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
		h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
		h = _mm_slli_si128(h, 2);
		for (j = 0; LIKELY(j < slen); ++j) {
			h = _mm_adds_epi16(h, *S++);
			e = _mm_load_si128(E + j);
			h = _mm_max_epi16(h, e);
			h = _mm_max_epi16(h, f);
			max = _mm_max_epi16(max, h);
			_mm_store_si128(H1 + j, h);
			e = _mm_subs_epu16(e, e_del);
			t = _mm_subs_epu16(h, oe_del);
			e = _mm_max_epi16(e, t);
			_mm_store_si128(E + j, e);
			f = _mm_subs_epu16(f, e_ins);
			t = _mm_subs_epu16(h, oe_ins);
			f = _mm_max_epi16(f, t);
			h = _mm_load_si128(H0 + j);
		}
		for (k = 0; LIKELY(k < 16); ++k) {
			f = _mm_slli_si128(f, 2);
			for (j = 0; LIKELY(j < slen); ++j) {
				h = _mm_load_si128(H1 + j);
				h = _mm_max_epi16(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu16(h, oe_ins);
				f = _mm_subs_epu16(f, e_ins);
				if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
			}
		}
end_loop8:
		__max_8(imax, max);
		if (imax >= minsc) {
			if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
				if (n_b == m_b) {
					m_b = m_b? m_b<<1 : 8;
					b = (uint64_t*)realloc(b, 8 * m_b);
				}
				b[n_b++] = (uint64_t)imax<<32 | i;
			} else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
		}
		if (imax > gmax) {
			gmax = imax; te = i;
			for (j = 0; LIKELY(j < slen); ++j)
				_mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
			if (gmax >= endsc) break;
		}
		S = H1; H1 = H0; H0 = S;
	}
	r.score = gmax; r.te = te;
	{
		int max = -1, tmp, low, high, qlen = slen * 8;
		uint16_t *t = (uint16_t*)Hmax;
		for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
			if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
			else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp; 
		if (b) {
			i = (r.score + q->max - 1) / q->max;
			low = te - i; high = te + i;
			for (i = 0; i < n_b; ++i) {
				int e = (int32_t)b[i];
				if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
					r.score2 = b[i]>>32, r.te2 = e;
			}
		}
	}
	free(b);
	return r;
}

static inline void revseq(int l, uint8_t *s)
{
	int i, t;
	for (i = 0; i < l>>1; ++i)
		t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
}

kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int xtra, kswq_t **qry)
{
	int size;
	kswq_t *q;
	kswr_t r, rr;
	kswr_t (*func)(kswq_t*, int, const uint8_t*, int, int, int, int, int);

	q = (qry && *qry)? *qry : ksw_qinit((xtra&KSW_XBYTE)? 1 : 2, qlen, query, m, mat);
	if (qry && *qry == 0) *qry = q;
	func = q->size == 2? ksw_i16 : ksw_u8;
	size = q->size;
	r = func(q, tlen, target, o_del, e_del, o_ins, e_ins, xtra);
	if (qry == 0) free(q);
	if ((xtra&KSW_XSTART) == 0 || ((xtra&KSW_XSUBO) && r.score < (xtra&0xffff))) return r;
	revseq(r.qe + 1, query); revseq(r.te + 1, target); // +1 because qe/te points to the exact end, not the position after the end
	q = ksw_qinit(size, r.qe + 1, query, m, mat);
	rr = func(q, tlen, target, o_del, e_del, o_ins, e_ins, KSW_XSTOP | r.score);
	revseq(r.qe + 1, query); revseq(r.te + 1, target);
	free(q);
	if (r.score == rr.score)
		r.tb = r.te - rr.te, r.qb = r.qe - rr.qe;
	return r;
}

kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry)
{
	return ksw_align2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, xtra, qry);
}

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

int ksw_extend2_scalar(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	assert(h0 > 0);
	// allocate memory
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; LIKELY(i < tlen); ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = beg; LIKELY(j < end); ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
			e = e > t? e : t;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m == 0) break;
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
		beg = j;
		for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
		end = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

/*****************************************
 *** SW extension, vectorised (Tronko) ***
 *****************************************/

/* ksw_extend2() below returns, for every input, the same six results (score, qle, tle, gtle,
 * gscore, max_off) as ksw_extend2_scalar() above, which is BWA's ksw_extend2 unchanged.
 *
 * The kernels keep the scalar code's setup, row order, band, exits and band narrowing, and
 * vectorise only the cell loop of one row: eight (SSE2) or sixteen (AVX2) cells per vector, in
 * 16-bit integers. That is exact because gaps open from M, the diagonal term, never from H:
 * within row i, M(i,j) and E(i+1,j) depend only on row i-1, and the horizontal gap
 * F(i,j+1) = max(F(i,j) - e_ins, max(M(i,j) - oe_ins, 0)), F(i,beg) = 0, is a max-plus prefix
 * scan over values of row i-1, which shifts compute exactly in any grouping. H = max(M, E, F)
 * is then per cell. The rest is kept as follows:
 * - memory: a row writes H and E at [beg, end] and nothing else. Cells beyond the band keep
 *   what an earlier row left there, which a later row reads when the band grows back, so the
 *   last, partial vector of a row writes back the old contents of the cells at and past end,
 *   loaded before any write.
 * - ties: the row maximum goes to the last j holding it, the best score to the first row that
 *   reaches it, the whole-query score to the last row; the per-row code after the cells
 *   (KS_ROW_EPILOGUE) is the scalar code.
 * - integer range: every value stays in [min(mat) - max(oe_del, oe_ins), h0 + qlen * max(mat, 0)];
 *   ks_eligible() sends any call where that bound could leave 16 bits, or with parameters
 *   outside the range checked, to the scalar code.
 *
 * ksw_extend2() uses the AVX2 kernel when the CPU has AVX2 and the SSE2 kernel otherwise (every
 * x86-64 CPU; on aarch64 through sse2neon.h). The AVX2 kernel is compiled with a target
 * attribute, so the build flags do not change and one binary runs on any x86-64 CPU.
 * tronko-assign/tests/unit/test_ksw_extend.c compares every kernel with ksw_extend2_scalar()
 * on recorded and random calls. */

/* Unaligned 16-bit vector loads and stores. On aarch64, sse2neon.h implements
 * _mm_loadu_si128/_mm_storeu_si128 through int32 pointers, which UBSan reports as misaligned
 * for the 2-byte-aligned int16 arrays here; use the int16 NEON forms there instead. */
#if defined(__aarch64__)
#define KS_LOADU(p) vreinterpretq_s32_s16(vld1q_s16((const int16_t*)(p)))
#define KS_STOREU(p, x) vst1q_s16((int16_t*)(p), vreinterpretq_s16_s32(x))
#else
#define KS_LOADU(p) _mm_loadu_si128((const __m128i*)(p))
#define KS_STOREU(p, x) _mm_storeu_si128((__m128i*)(p), (x))
#endif

#define KS_PAD 32            /* int16 elements of slack before and after every row array */
#define KS_STACK_ELEMS 16384 /* 32 KiB on the stack; longer queries use malloc */

/* Whether the 16-bit kernels may run. Every value they store is in
 * [min(mat) - max(oe_del, oe_ins), h0 + qlen * max(mat, 0)]; the scan subtracts at most
 * 16 * e_ins more and the AVX2 carry at most 30000 from a non-negative value. With the limits
 * below, all of it stays inside int16. Anything else runs the scalar code, which then behaves
 * exactly as before (including its assert on h0 and its band arithmetic for w < 0). */
static int ks_eligible(int qlen, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int h0)
{
	int i, maxpos = 0;
	long bound;
	if (qlen < 1 || qlen > 16000 || m < 1 || m > 16 || w < 0) return 0;
	if (e_del < 1 || e_ins < 1 || e_del > 1000 || e_ins > 1000) return 0;
	if (o_del < 0 || o_ins < 0 || o_del > 4000 || o_ins > 4000) return 0;
	if (h0 < 1) return 0;
	for (i = 0; i < m * m; ++i) if (mat[i] > maxpos) maxpos = mat[i];
	bound = (long)h0 + (long)qlen * maxpos;
	return bound <= 20000;
}

typedef struct {
	int16_t *H;   /* H[j] mirrors the scalar code's eh[j].h, j in [0, qlen] */
	int16_t *E;   /* E[j] mirrors eh[j].e */
	int16_t *P;   /* query profile, row k at P + k * stride */
	int stride;
	void *heap;
} ks_buf_t;

/* Three zeroed arrays of stride elements plus a profile of m rows, with KS_PAD elements of
 * slack on each side, so that full-vector loads and stores near qlen stay inside the buffer. */
static void ks_buf_init(ks_buf_t *b, int16_t *stack, int qlen, int m, const int8_t *mat, const uint8_t *query, int fill_profile)
{
	int j, k, stride = ((qlen + 1 + 2 * KS_PAD) + 15) & ~15;
	size_t n = (size_t)(m + 2) * stride;
	int16_t *base;
	b->heap = 0;
	if (n <= KS_STACK_ELEMS) base = stack;
	else base = b->heap = malloc(n * sizeof(int16_t));
	memset(base, 0, n * sizeof(int16_t));
	b->stride = stride;
	b->H = base + KS_PAD;
	b->E = base + stride + KS_PAD;
	b->P = base + 2 * stride + KS_PAD;
	if (fill_profile) for (k = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		int16_t *pk = b->P + k * stride;
		for (j = 0; j < qlen; ++j) pk[j] = p[query[j]];
	}
}

/* The first row and the band clamp, as in ksw_extend2_scalar(). */
static int ks_setup(int16_t *H, int qlen, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int h0, int fill_first_row)
{
	int i, j, k, max, max_ins, max_del, oe_ins = o_ins + e_ins;
	if (fill_first_row) {
		H[0] = h0; H[1] = h0 > oe_ins? h0 - oe_ins : 0;
		for (j = 2; j <= qlen && H[j-1] > e_ins; ++j)
			H[j] = H[j-1] - e_ins;
	}
	k = m * m;
	for (i = 0, max = 0; i < k; ++i)
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del;
	return w;
}

/* The per-row code after the cells, as in ksw_extend2_scalar(): the whole-query score, the
 * zero-row exit, the best score, the z-drop and the band narrowing for the next row. jend is
 * where the scalar cell loop would have left j. */
#define KS_ROW_EPILOGUE(H, E) do { \
	if (jend == qlen) { \
		max_ie = gscore > h1? max_ie : i; \
		gscore = gscore > h1? gscore : h1; \
	} \
	if (mrow == 0) goto done; \
	if (mrow > max) { \
		max = mrow, max_i = i, max_j = mj; \
		max_off = max_off > abs(mj - i)? max_off : abs(mj - i); \
	} else if (zdrop > 0) { \
		if (i - max_i > mj - max_j) { \
			if (max - mrow - ((i - max_i) - (mj - max_j)) * e_del > zdrop) goto done; \
		} else { \
			if (max - mrow - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) goto done; \
		} \
	} \
	for (j = beg; j < end && H[j] == 0 && E[j] == 0; ++j); \
	beg = j; \
	for (j = end; j >= beg && H[j] == 0 && E[j] == 0; --j); \
	end = j + 2 < qlen? j + 2 : qlen; \
} while (0)

/* SSE2 kernel: 8 cells of one row per vector. */

static inline int ks_hmax8(__m128i x)
{
	x = _mm_max_epi16(x, _mm_shuffle_epi32(x, 0x4E));
	x = _mm_max_epi16(x, _mm_shuffle_epi32(x, 0xB1));
	x = _mm_max_epi16(x, _mm_srli_si128(x, 2)); /* only lane 0 is used: max(lane 0, lane 1) */
	return (int16_t)_mm_cvtsi128_si32(x);
}

static int ks_extend_sse2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	int16_t stack[KS_STACK_ELEMS] __attribute__((aligned(32)));
	ks_buf_t b;
	int16_t *H, *E;
	int i, j, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ie, gscore, max_off;
	const __m128i zero = _mm_setzero_si128(), ones = _mm_set1_epi16(-1);
	const __m128i v_oe_del = _mm_set1_epi16(oe_del), v_e_del = _mm_set1_epi16(e_del), v_oe_ins = _mm_set1_epi16(oe_ins);
	const __m128i ve1 = _mm_set1_epi16(e_ins), ve2 = _mm_set1_epi16(2 * e_ins), ve4 = _mm_set1_epi16(4 * e_ins);
	/* _mm_set_epi16 takes lane 7 first (sse2neon.h has no _mm_setr_epi16) */
	const __m128i vramp = _mm_set_epi16(8*e_ins, 7*e_ins, 6*e_ins, 5*e_ins, 4*e_ins, 3*e_ins, 2*e_ins, e_ins);
	const __m128i lane = _mm_set_epi16(7, 6, 5, 4, 3, 2, 1, 0);

	ks_buf_init(&b, stack, qlen, m, mat, query, 1);
	H = b.H, E = b.E;
	w = ks_setup(H, qlen, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, h0, 1);
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; i < tlen; ++i) {
		int h1, mrow, mj, jend;
		const int16_t *q = b.P + target[i] * b.stride;
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		if (beg < end) {
			__m128i Hcur = KS_LOADU(H + beg); /* H(i-1, j-1), read before any write */
			__m128i Sprev = zero, vmax = ones, vidx = ones;
			H[beg] = (int16_t)h1;
			for (j = beg; j < end; j += 8) {
				__m128i Hnext = KS_LOADU(H + j + 8); /* before H[j+8] is overwritten */
				__m128i Ev = KS_LOADU(E + j);
				__m128i Qv = KS_LOADU(q + j);
				__m128i M = _mm_andnot_si128(_mm_cmpeq_epi16(Hcur, zero), _mm_add_epi16(Hcur, Qv));
				__m128i S, c, F, Hv, Ev1, jv, ge;
				/* S(j) = F(i, j+1): inclusive max-plus scan of g = max(M - oe_ins, 0), decay e_ins */
				S = _mm_max_epi16(_mm_sub_epi16(M, v_oe_ins), zero);
				S = _mm_max_epi16(S, _mm_sub_epi16(_mm_slli_si128(S, 2), ve1));
				S = _mm_max_epi16(S, _mm_sub_epi16(_mm_slli_si128(S, 4), ve2));
				S = _mm_max_epi16(S, _mm_sub_epi16(_mm_slli_si128(S, 8), ve4));
				c = _mm_shuffle_epi32(_mm_shufflehi_epi16(Sprev, 0xFF), 0xFF); /* broadcast S(j0 - 1) */
				S = _mm_max_epi16(S, _mm_sub_epi16(c, vramp));
				/* F(i, j) = S(j - 1); F(i, beg) = 0 */
				F = _mm_or_si128(_mm_slli_si128(S, 2), _mm_srli_si128(Sprev, 14));
				Hv = _mm_max_epi16(_mm_max_epi16(M, Ev), F);
				Ev1 = _mm_max_epi16(_mm_sub_epi16(Ev, v_e_del), _mm_max_epi16(_mm_sub_epi16(M, v_oe_del), zero));
				jv = _mm_add_epi16(_mm_set1_epi16(j), lane);
				if (j + 8 <= end) {
					KS_STOREU(E + j, Ev1);
					KS_STOREU(H + j + 1, Hv); /* eh[j+1].h = H(i, j) */
				} else { /* last, partial vector: cells j >= end keep their old values */
					__m128i valid = _mm_cmpgt_epi16(_mm_set1_epi16(end), jv);
					/* old H[j+1 .. j+8] = [Hcur[1..7], Hnext[0]], both loaded before any write */
					__m128i oldH = _mm_or_si128(_mm_srli_si128(Hcur, 2), _mm_slli_si128(Hnext, 14));
					KS_STOREU(E + j, _mm_or_si128(_mm_and_si128(valid, Ev1), _mm_andnot_si128(valid, Ev)));
					KS_STOREU(H + j + 1, _mm_or_si128(_mm_and_si128(valid, Hv), _mm_andnot_si128(valid, oldH)));
					Hv = _mm_or_si128(_mm_and_si128(valid, Hv), _mm_andnot_si128(valid, ones));
				}
				/* per lane: running max and the last j reaching it (H >= running max) */
				ge = _mm_andnot_si128(_mm_cmpgt_epi16(vmax, Hv), ones);
				vidx = _mm_or_si128(_mm_and_si128(ge, jv), _mm_andnot_si128(ge, vidx));
				vmax = _mm_max_epi16(vmax, Hv);
				Sprev = S;
				Hcur = Hnext;
			}
			mrow = ks_hmax8(vmax);
			{
				__m128i sel = _mm_cmpeq_epi16(vmax, _mm_set1_epi16(mrow));
				mj = ks_hmax8(_mm_or_si128(_mm_and_si128(sel, vidx), _mm_andnot_si128(sel, ones)));
			}
			h1 = H[end];
			jend = end;
		} else {
			mrow = 0, mj = -1, jend = beg;
			H[end] = (int16_t)h1;
		}
		E[end] = 0;
		KS_ROW_EPILOGUE(H, E);
	}
done:
	if (b.heap) free(b.heap);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

/* AVX2 kernel: 16 cells of one row per vector (x86-64 only, chosen at run time). */

#if defined(__x86_64__) || defined(__i386__)
#define KS_AVX2 __attribute__((target("avx2")))

/* Query profile with a byte shuffle: row k is mat[k*m + query[j]], sign-extended to 16 bits.
 * Falls back to the scalar loop if any query byte is >= m (the scalar code's contract is
 * 0 <= query[j] < m; this keeps even out-of-contract inputs identical). */
static KS_AVX2 void ks_profile_avx2(int16_t *P, int stride, int qlen, int m, const int8_t *mat, const uint8_t *query)
{
	int j, k, n16 = qlen & ~15;
	__m128i mx = _mm_setzero_si128();
	for (j = 0; j < n16; j += 16) mx = _mm_max_epu8(mx, _mm_loadu_si128((const __m128i*)(query + j)));
	for (j = n16; j < qlen; ++j) if (query[j] >= m) n16 = -1;
	if (n16 >= 0) {
		uint8_t t[16];
		_mm_storeu_si128((__m128i*)t, mx);
		for (j = 0; j < 16; ++j) if (t[j] >= m) n16 = -1;
	}
	for (k = 0; k < m; ++k) {
		int16_t *pk = P + k * stride;
		j = 0;
		if (n16 > 0) {
			int8_t tb[16];
			__m128i tv;
			memset(tb, 0, sizeof(tb));
			memcpy(tb, &mat[k * m], m);
			tv = _mm_loadu_si128((const __m128i*)tb);
			for (; j < n16; j += 16)
				_mm256_storeu_si256((__m256i*)(pk + j), _mm256_cvtepi8_epi16(_mm_shuffle_epi8(tv, _mm_loadu_si128((const __m128i*)(query + j)))));
		}
		for (; j < qlen; ++j) pk[j] = mat[k * m + query[j]];
	}
}

/* The first row of ksw_extend2_scalar(): H[0] = h0 and, for 1 <= j <= qlen,
 * H[j] = max(h0 - oe_ins - (j - 1) * e_ins, 0) (the scalar loop stops at the first value <= 0
 * and leaves zeros, which is the same thing for e_ins >= 1). Saturating 16-bit arithmetic;
 * the last store may write up to 15 elements past qlen, into the padding. */
static KS_AVX2 void ks_first_row_avx2(int16_t *H, int qlen, int o_ins, int e_ins, int h0)
{
	int j, base = h0 - (o_ins + e_ins); /* value at j = 1 */
	const __m256i ramp = _mm256_setr_epi16(0, e_ins, 2*e_ins, 3*e_ins, 4*e_ins, 5*e_ins, 6*e_ins, 7*e_ins,
			8*e_ins, 9*e_ins, 10*e_ins, 11*e_ins, 12*e_ins, 13*e_ins, 14*e_ins, 15*e_ins);
	H[0] = h0;
	for (j = 1; j <= qlen && base > 0; j += 16, base -= 16 * e_ins)
		_mm256_storeu_si256((__m256i*)(H + j), _mm256_max_epi16(_mm256_subs_epi16(_mm256_set1_epi16(base), ramp), _mm256_setzero_si256()));
}

static inline KS_AVX2 int ks_hmax32(__m256i x)
{
	__m128i y = _mm_max_epi32(_mm256_castsi256_si128(x), _mm256_extracti128_si256(x, 1));
	y = _mm_max_epi32(y, _mm_shuffle_epi32(y, 0x4E));
	y = _mm_max_epi32(y, _mm_shuffle_epi32(y, 0xB1));
	return _mm_cvtsi128_si32(y);
}

static KS_AVX2 int ks_extend_avx2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	int16_t stack[KS_STACK_ELEMS] __attribute__((aligned(32)));
	ks_buf_t b;
	int16_t *H, *E;
	int i, j, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ie, gscore, max_off;
	const __m256i zero = _mm256_setzero_si256(), ones = _mm256_set1_epi16(-1);
	const __m256i v_oe_del = _mm256_set1_epi16(oe_del), v_e_del = _mm256_set1_epi16(e_del), v_oe_ins = _mm256_set1_epi16(oe_ins);
	const __m256i ve1 = _mm256_set1_epi16(e_ins), ve2 = _mm256_set1_epi16(2 * e_ins), ve4 = _mm256_set1_epi16(4 * e_ins);
	const __m256i vramp_hi = _mm256_setr_epi16(30000, 30000, 30000, 30000, 30000, 30000, 30000, 30000,
			e_ins, 2*e_ins, 3*e_ins, 4*e_ins, 5*e_ins, 6*e_ins, 7*e_ins, 8*e_ins);
	const __m256i vramp = _mm256_setr_epi16(e_ins, 2*e_ins, 3*e_ins, 4*e_ins, 5*e_ins, 6*e_ins, 7*e_ins, 8*e_ins,
			9*e_ins, 10*e_ins, 11*e_ins, 12*e_ins, 13*e_ins, 14*e_ins, 15*e_ins, 16*e_ins);
	const __m256i lane = _mm256_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

	ks_buf_init(&b, stack, qlen, m, mat, query, 0);
	H = b.H, E = b.E;
	ks_profile_avx2(b.P, b.stride, qlen, m, mat, query);
	ks_first_row_avx2(H, qlen, o_ins, e_ins, h0);
	w = ks_setup(H, qlen, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, h0, 0);
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; i < tlen; ++i) {
		int h1, mrow, mj, jend;
		const int16_t *q = b.P + target[i] * b.stride;
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		if (beg < end) {
			__m256i Hcur = _mm256_loadu_si256((const __m256i*)(H + beg));
			__m256i Sprev = zero, k0 = ones, k1 = ones; /* keys (H << 16 | j), -1 = none */
			H[beg] = (int16_t)h1;
			for (j = beg; j < end; j += 16) {
				__m256i Hnext = _mm256_loadu_si256((const __m256i*)(H + j + 16));
				__m256i Ev = _mm256_loadu_si256((const __m256i*)(E + j));
				__m256i Qv = _mm256_loadu_si256((const __m256i*)(q + j));
				__m256i M = _mm256_andnot_si256(_mm256_cmpeq_epi16(Hcur, zero), _mm256_add_epi16(Hcur, Qv));
				__m256i S, c, F, Hv, Ev1, jv;
				S = _mm256_max_epi16(_mm256_sub_epi16(M, v_oe_ins), zero);
				/* scan within each 128-bit half (in-lane byte shifts, zeros shifted in) */
				S = _mm256_max_epi16(S, _mm256_sub_epi16(_mm256_bslli_epi128(S, 2), ve1));
				S = _mm256_max_epi16(S, _mm256_sub_epi16(_mm256_bslli_epi128(S, 4), ve2));
				S = _mm256_max_epi16(S, _mm256_sub_epi16(_mm256_bslli_epi128(S, 8), ve4));
				/* upper half: S[7] - (l - 7) e_ins from the lower half; the lower half gets
				 * S[7] - 30000, below every S */
				c = _mm256_permute4x64_epi64(S, 0x55);
				c = _mm256_shufflehi_epi16(c, 0xFF);
				c = _mm256_unpackhi_epi64(c, c);
				S = _mm256_max_epi16(S, _mm256_sub_epi16(c, vramp_hi));
				/* previous vector: S(j0 - 1) - (l + 1) e_ins */
				c = _mm256_permute4x64_epi64(Sprev, 0xFF);
				c = _mm256_shufflehi_epi16(c, 0xFF);
				c = _mm256_unpackhi_epi64(c, c); /* broadcast S(j0 - 1) */
				S = _mm256_max_epi16(S, _mm256_sub_epi16(c, vramp));
				/* F = [Sprev[15], S[0..14]] */
				F = _mm256_alignr_epi8(S, _mm256_permute2x128_si256(Sprev, S, 0x21), 14);
				Hv = _mm256_max_epi16(_mm256_max_epi16(M, Ev), F);
				Ev1 = _mm256_max_epi16(_mm256_sub_epi16(Ev, v_e_del), _mm256_max_epi16(_mm256_sub_epi16(M, v_oe_del), zero));
				jv = _mm256_add_epi16(_mm256_set1_epi16(j), lane);
				if (j + 16 <= end) {
					_mm256_storeu_si256((__m256i*)(E + j), Ev1);
					_mm256_storeu_si256((__m256i*)(H + j + 1), Hv);
				} else { /* last, partial vector: cells j >= end keep their old values */
					__m256i valid = _mm256_cmpgt_epi16(_mm256_set1_epi16(end), jv);
					/* old H[j+1 .. j+16] = [Hcur[1..15], Hnext[0]], both loaded before any write */
					__m256i oldH = _mm256_alignr_epi8(_mm256_permute2x128_si256(Hcur, Hnext, 0x21), Hcur, 2);
					_mm256_storeu_si256((__m256i*)(E + j), _mm256_blendv_epi8(Ev, Ev1, valid));
					_mm256_storeu_si256((__m256i*)(H + j + 1), _mm256_blendv_epi8(oldH, Hv, valid));
					Hv = _mm256_blendv_epi8(ones, Hv, valid);
				}
				/* row maximum and the last j holding it, as one 32-bit key per cell: the largest
				 * key has the largest H and, among equal H, the largest j; H = -1 marks no cell */
				k0 = _mm256_max_epi32(k0, _mm256_unpacklo_epi16(jv, Hv));
				k1 = _mm256_max_epi32(k1, _mm256_unpackhi_epi16(jv, Hv));
				Sprev = S;
				Hcur = Hnext;
			}
			{
				int key = ks_hmax32(_mm256_max_epi32(k0, k1)); /* >= 0: the row has a cell */
				mrow = key >> 16, mj = key & 0xffff;
			}
			h1 = H[end];
			jend = end;
		} else {
			mrow = 0, mj = -1, jend = beg;
			H[end] = (int16_t)h1;
		}
		E[end] = 0;
		KS_ROW_EPILOGUE(H, E);
	}
done:
	if (b.heap) free(b.heap);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}
#endif

int ksw_extend2_sse2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	if (ks_eligible(qlen, m, mat, o_del, e_del, o_ins, e_ins, w, h0))
		return ks_extend_sse2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
	return ksw_extend2_scalar(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
}

int ksw_extend2_avx2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
#if defined(__x86_64__) || defined(__i386__)
	if (ks_eligible(qlen, m, mat, o_del, e_del, o_ins, e_ins, w, h0))
		return ks_extend_avx2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
#endif
	return ksw_extend2_sse2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
}

int ksw_simd_has_avx2(void)
{
#if defined(__x86_64__) || defined(__i386__)
	__builtin_cpu_init();
	return __builtin_cpu_supports("avx2");
#else
	return 0;
#endif
}

static int ks_use_avx2 = -1; /* -1 until the first call; every thread computes the same value */

int ksw_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	int avx2 = __atomic_load_n(&ks_use_avx2, __ATOMIC_RELAXED);
	if (avx2 < 0) {
		avx2 = ksw_simd_has_avx2();
		__atomic_store_n(&ks_use_avx2, avx2, __ATOMIC_RELAXED);
	}
	if (avx2) return ksw_extend2_avx2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
	return ksw_extend2_sse2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off);
}

int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off)
{
	return ksw_extend2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
}

/********************
 * Global alignment *
 ********************/

#define MINUS_INF -0x40000000

static inline uint32_t *push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, int op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

int ksw_global2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, score, n_col;
	uint8_t *z; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be a little more complex
	if (n_cigar_) *n_cigar_ = 0;
	// allocate memory
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	z = n_cigar_ && cigar_? malloc((long)n_col * tlen) : 0;
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = 0; eh[0].e = MINUS_INF;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(o_ins + e_ins * j), eh[j].e = MINUS_INF;
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = MINUS_INF; // everything is -inf outside the band
	// DP loop
	for (i = 0; LIKELY(i < tlen); ++i) { // target sequence is in the outer loop
		int32_t f = MINUS_INF, h1, beg, end, t;
		int8_t *q = &qp[target[i] * qlen];
		beg = i > w? i - w : 0;
		end = i + w + 1 < qlen? i + w + 1 : qlen; // only loop through [beg,end) of the query sequence
		h1 = beg == 0? -(o_del + e_del * (i + 1)) : MINUS_INF;
		if (n_cigar_ && cigar_) {
			uint8_t *zi = &z[(long)i * n_col];
			for (j = beg; LIKELY(j < end); ++j) {
				// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
				// Cells are computed in the following order:
				//   M(i,j)   = H(i-1,j-1) + S(i,j)
				//   H(i,j)   = max{M(i,j), E(i,j), F(i,j)}
				//   E(i+1,j) = max{M(i,j)-gapo, E(i,j)} - gape
				//   F(i,j+1) = max{M(i,j)-gapo, F(i,j)} - gape
				// We have to separate M(i,j); otherwise the direction may not be recorded correctly.
				// However, a CIGAR like "10M3I3D10M" allowed by local() is disallowed by global().
				// Such a CIGAR may occur, in theory, if mismatch_penalty > 2*gap_ext_penalty + 2*gap_open_penalty/k.
				// In practice, this should happen very rarely given a reasonable scoring system.
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				uint8_t d; // direction
				p->h = h1;
				m += q[j];
				d = m >= e? 0 : 1;
				h = m >= e? m : e;
				d = h >= f? d : 2;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				d |= e > t? 1<<2 : 0;
				e  = e > t? e    : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				d |= f > t? 2<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
				f  = f > t? f    : t;
				zi[j - beg] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
			}
		} else {
			for (j = beg; LIKELY(j < end); ++j) {
				eh_t *p = &eh[j];
				int32_t h, m = p->h, e = p->e;
				p->h = h1;
				m += q[j];
				h = m >= e? m : e;
				h = h >= f? h : f;
				h1 = h;
				t = m - oe_del;
				e -= e_del;
				e  = e > t? e : t;
				p->e = e;
				t = m - oe_ins;
				f -= e_ins;
				f  = f > t? f : t;
			}
		}
		eh[end].h = h1; eh[end].e = MINUS_INF;
	}
	score = eh[qlen].h;
	if (n_cigar_ && cigar_) { // backtrack
		int n_cigar = 0, m_cigar = 0, which = 0;
		uint32_t *cigar = 0, tmp;
		i = tlen - 1; k = (i + w + 1 < qlen? i + w + 1 : qlen) - 1; // (i,k) points to the last cell
		while (i >= 0 && k >= 0) {
			which = z[(long)i * n_col + (k - (i > w? i - w : 0))] >> (which<<1) & 3;
			if (which == 0)      cigar = push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --k;
			else if (which == 1) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i;
			else                 cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --k;
		}
		if (i >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 2, i + 1);
		if (k >= 0) cigar = push_cigar(&n_cigar, &m_cigar, cigar, 1, k + 1);
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
		*n_cigar_ = n_cigar, *cigar_ = cigar;
	}
	free(eh); free(qp); free(z);
	return score;
}

int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *n_cigar_, uint32_t **cigar_)
{
	return ksw_global2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, n_cigar_, cigar_);
}

/*******************************************
 * Main function (not compiled by default) *
 *******************************************/

#ifdef _KSW_MAIN

#include <unistd.h>
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, err_gzread)

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int main(int argc, char *argv[])
{
	int c, sa = 1, sb = 3, i, j, k, forward_only = 0, max_rseq = 0;
	int8_t mat[25];
	int gapo = 5, gape = 2, minsc = 0, xtra = KSW_XSTART;
	uint8_t *rseq = 0;
	gzFile fpt, fpq;
	kseq_t *kst, *ksq;

	// parse command line
	while ((c = getopt(argc, argv, "a:b:q:r:ft:1")) >= 0) {
		switch (c) {
			case 'a': sa = atoi(optarg); break;
			case 'b': sb = atoi(optarg); break;
			case 'q': gapo = atoi(optarg); break;
			case 'r': gape = atoi(optarg); break;
			case 't': minsc = atoi(optarg); break;
			case 'f': forward_only = 1; break;
			case '1': xtra |= KSW_XBYTE; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: ksw [-1] [-f] [-a%d] [-b%d] [-q%d] [-r%d] [-t%d] <target.fa> <query.fa>\n", sa, sb, gapo, gape, minsc);
		return 1;
	}
	if (minsc > 0xffff) minsc = 0xffff;
	xtra |= KSW_XSUBO | minsc;
	// initialize scoring matrix
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
	// open file
	fpt = xzopen(argv[optind],   "r"); kst = kseq_init(fpt);
	fpq = xzopen(argv[optind+1], "r"); ksq = kseq_init(fpq);
	// all-pair alignment
	while (kseq_read(ksq) > 0) {
		kswq_t *q[2] = {0, 0};
		kswr_t r;
		for (i = 0; i < (int)ksq->seq.l; ++i) ksq->seq.s[i] = seq_nt4_table[(int)ksq->seq.s[i]];
		if (!forward_only) { // reverse
			if ((int)ksq->seq.m > max_rseq) {
				max_rseq = ksq->seq.m;
				rseq = (uint8_t*)realloc(rseq, max_rseq);
			}
			for (i = 0, j = ksq->seq.l - 1; i < (int)ksq->seq.l; ++i, --j)
				rseq[j] = ksq->seq.s[i] == 4? 4 : 3 - ksq->seq.s[i];
		}
		gzrewind(fpt); kseq_rewind(kst);
		while (kseq_read(kst) > 0) {
			for (i = 0; i < (int)kst->seq.l; ++i) kst->seq.s[i] = seq_nt4_table[(int)kst->seq.s[i]];
			r = ksw_align(ksq->seq.l, (uint8_t*)ksq->seq.s, kst->seq.l, (uint8_t*)kst->seq.s, 5, mat, gapo, gape, xtra, &q[0]);
			if (r.score >= minsc)
				err_printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", kst->name.s, r.tb, r.te+1, ksq->name.s, r.qb, r.qe+1, r.score, r.score2, r.te2);
			if (rseq) {
				r = ksw_align(ksq->seq.l, rseq, kst->seq.l, (uint8_t*)kst->seq.s, 5, mat, gapo, gape, xtra, &q[1]);
				if (r.score >= minsc)
					err_printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n", kst->name.s, r.tb, r.te+1, ksq->name.s, (int)ksq->seq.l - r.qb, (int)ksq->seq.l - 1 - r.qe, r.score, r.score2, r.te2);
			}
		}
		free(q[0]); free(q[1]);
	}
	free(rseq);
	kseq_destroy(kst); err_gzclose(fpt);
	kseq_destroy(ksq); err_gzclose(fpq);
	return 0;
}
#endif
