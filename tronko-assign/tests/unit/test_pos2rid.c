/*
 * Unit test for the position-to-accession lookup of bwa_source_files/bntseq.c.
 *
 * bns_pos2rid() maps a forward-strand position to the index of the accession that contains it.
 * tronko-assign's copy reads a bucket table, the answer at the start of every 2^BNS_RID_BUCKET_SHIFT
 * bases built when the index is loaded, and steps forward from there, instead of binary-searching
 * every accession. It must return what BWA's binary search returns for every position, or seeds
 * are assigned to other accessions and Tronko's output changes. The test compares bns_pos2rid()
 * and bns_intv2rid() with verbatim copies of BWA's originals (base_pos2rid(), base_intv2rid(), from
 * tronko-assign at 71f6ec3) on:
 *   1. the repository's fixture index (tronko-build/example_datasets/single_tree/Charadriiformes):
 *      the table against the original at every bucket start; every position from -3 to l_pac + 3;
 *      intervals at every accession boundary (offset - 1, offset, offset + len - 1, offset + len)
 *      on both strands;
 *   2. synthetic references with zero-length accessions, accessions from 1 base to several buckets,
 *      and total lengths of exact multiples of a bucket: the table at every bucket start, every
 *      position from -3 to l_pac + 3, and random intervals on both strands;
 *   3. references where no table may be built (offsets that decrease, no bases, no accessions);
 *   4. the index loaded as tronko-assign loads it (bwa_idx_load, which must build the table) and
 *      then converted to BWA's shared-memory layout (bwa_idx2mem, as `bwa shm` does), which must
 *      keep BWA's own layout, without the table pointer, and fall back to the binary search.
 *
 * Usage: test_pos2rid <index prefix of the fixture FASTA>
 * Exit status 0 when every comparison agrees.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>
#include "../../bwa_source_files/bntseq.h"
#include "../../bwa_source_files/bwa.h"
#include "../../bwa_source_files/utils.h"
#include "../../global.h"

/* Defined in tronko-assign.c; BWA's tronko-specific code refers to them. */
queryMatPaired *pairedQueryMat;
queryMatSingle *singleQueryMat;

/* bns_pos2rid() and bns_intv2rid() of BWA as vendored in tronko-assign at 71f6ec3, verbatim apart
 * from their names. */
static int base_pos2rid(const bntseq_t *bns, int64_t pos_f)
{
	int left, mid, right;
	if (pos_f >= bns->l_pac) return -1;
	left = 0; mid = 0; right = bns->n_seqs;
	while (left < right) { // binary search
		mid = (left + right) >> 1;
		if (pos_f >= bns->anns[mid].offset) {
			if (mid == bns->n_seqs - 1) break;
			if (pos_f < bns->anns[mid+1].offset) break; // bracketed
			left = mid + 1;
		} else right = mid;
	}
	return mid;
}

static int base_intv2rid(const bntseq_t *bns, int64_t rb, int64_t re)
{
	int is_rev, rid_b, rid_e;
	if (rb < bns->l_pac && re > bns->l_pac) return -2;
	assert(rb <= re);
	rid_b = base_pos2rid(bns, bns_depos(bns, rb, &is_rev));
	rid_e = rb < re? base_pos2rid(bns, bns_depos(bns, re - 1, &is_rev)) : rid_b;
	return rid_b == rid_e? rid_b : -1;
}

static uint64_t rng_state = 0x2545F4914F6CDD1DULL;
static uint64_t rnd(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static long n_cmp, n_bad;

static void check(const char *part, const char *what, int64_t x, int64_t y, long expect, long got)
{
	++n_cmp;
	if (expect != got && ++n_bad <= 10)
		fprintf(stderr, "MISMATCH %s: %s(%lld, %lld): original %ld, new %ld\n", part, what, (long long)x, (long long)y, expect, got);
}

/* The table against the original at every bucket start. */
static void check_table(const char *part, const bntseq_t *bns)
{
	int64_t b, n_bucket = ((bns->l_pac - 1) >> BNS_RID_BUCKET_SHIFT) + 1;
	if (bns->rid_bucket == 0) { check(part, "table built", 0, 0, 1, 0); return; }
	for (b = 0; b < n_bucket; ++b)
		check(part, "rid_bucket", b, 0, base_pos2rid(bns, b << BNS_RID_BUCKET_SHIFT), bns->rid_bucket[b]);
}

/* Every position from -3 to l_pac + 3. */
static void check_positions(const char *part, const bntseq_t *bns)
{
	int64_t p;
	for (p = -3; p < bns->l_pac + 3; ++p)
		check(part, "bns_pos2rid", p, 0, base_pos2rid(bns, p), bns_pos2rid(bns, p));
}

static void check_interval(const char *part, const bntseq_t *bns, int64_t rb, int64_t re)
{
	if (rb < 0 || re < rb || re > bns->l_pac << 1) return;
	check(part, "bns_intv2rid", rb, re, base_intv2rid(bns, rb, re), bns_intv2rid(bns, rb, re));
}

/* Intervals that start or end at every accession boundary, on both strands. */
static void check_boundaries(const char *part, const bntseq_t *bns)
{
	int r, d, e;
	for (r = 0; r < bns->n_seqs; ++r) {
		int64_t off = bns->anns[r].offset, len = bns->anns[r].len;
		int64_t pts[4] = { off - 1, off, off + len - 1, off + len };
		for (d = 0; d < 4; ++d) {
			int64_t f = pts[d], rev = (bns->l_pac << 1) - 1 - pts[d];
			for (e = 0; e <= 20; e += 1 + 19 * (e > 0)) { // lengths 0, 1 and 21
				check_interval(part, bns, f, f + e);
				check_interval(part, bns, f - e, f);
				check_interval(part, bns, rev, rev + e);
				check_interval(part, bns, rev - e, rev);
			}
		}
	}
}

static void check_random_intervals(const char *part, const bntseq_t *bns, int n)
{
	int i;
	for (i = 0; i < n; ++i) {
		int64_t rb = (int64_t)(rnd() % (uint64_t)((bns->l_pac << 1) + 1)), re = rb + (int64_t)(rnd() % 600);
		check_interval(part, bns, rb, re);
	}
}

/* A synthetic reference: n accessions with the given lengths, contiguous as bwa index writes them. */
static bntseq_t *make_bns(int n, const int32_t *len)
{
	bntseq_t *bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	int64_t off = 0;
	int r;
	bns->n_seqs = n;
	bns->anns = (bntann1_t*)calloc(n > 0? n : 1, sizeof(bntann1_t));
	for (r = 0; r < n; ++r) bns->anns[r].offset = off, bns->anns[r].len = len[r], off += len[r];
	bns->l_pac = off;
	return bns;
}

static void free_bns(bntseq_t *bns)
{
	bns_rid_bucket_destroy(bns);
	free(bns->anns);
	free(bns);
}

static long part2(int n_refs)
{
	int32_t len[300];
	long n_pos = 0;
	int it, r;
	for (it = 0; it < n_refs; ++it) {
		int n = 1 + (int)(rnd() % 300);
		bntseq_t *bns;
		for (r = 0; r < n; ++r) {
			switch (rnd() % 6) {
				case 0: len[r] = 0; break;                                 // zero-length accession
				case 1: len[r] = 1 + (int32_t)(rnd() % 3); break;          // a few bases
				case 2: len[r] = 255 + (int32_t)(rnd() % 3); break;        // about one bucket
				case 3: len[r] = 1 + (int32_t)(rnd() % 1000); break;
				case 4: len[r] = 256 * (1 + (int32_t)(rnd() % 3)); break; // whole buckets
				default: len[r] = 3000 + (int32_t)(rnd() % 3000); break;   // several buckets
			}
		}
		if (it % 10 == 0) for (r = 0; r < n; ++r) len[r] = 256; // l_pac an exact multiple of a bucket
		bns = make_bns(n, len);
		if (bns->l_pac == 0) { free_bns(bns); continue; }
		bns_rid_bucket_build(bns);
		check_table("synthetic", bns);
		check_positions("synthetic", bns);
		check_random_intervals("synthetic", bns, 200);
		n_pos += bns->l_pac + 6;
		free_bns(bns);
	}
	return n_pos;
}

static void part3(void)
{
	static const int32_t zero[3] = { 0, 0, 0 }, some[4] = { 100, 300, 50, 700 };
	bntseq_t *bns;
	/* no bases at all: nothing to build */
	bns = make_bns(3, zero);
	bns_rid_bucket_build(bns);
	check("no bases", "table left out", 0, 0, 0, bns->rid_bucket != 0);
	check_positions("no bases", bns);
	free_bns(bns);
	/* no accessions */
	bns = make_bns(0, zero);
	bns->l_pac = 1000;
	bns_rid_bucket_build(bns);
	check("no accessions", "table left out", 0, 0, 0, bns->rid_bucket != 0);
	free_bns(bns);
	/* offsets that decrease: the forward scan would be wrong, so no table */
	bns = make_bns(4, some);
	bns->anns[2].offset = 10;
	bns_rid_bucket_build(bns);
	check("decreasing offsets", "table left out", 0, 0, 0, bns->rid_bucket != 0);
	check_positions("decreasing offsets", bns);
	free_bns(bns);
	/* built twice: the second call keeps the first table */
	bns = make_bns(4, some);
	bns_rid_bucket_build(bns);
	{
		int32_t *first = bns->rid_bucket;
		bns_rid_bucket_build(bns);
		check("built twice", "same table", 0, 0, 1, bns->rid_bucket == first);
	}
	check_positions("built twice", bns);
	free_bns(bns);
}

/* The fixture index as tronko-assign loads it, then in the shared-memory layout. */
static long part1_and_4(const char *prefix)
{
	bwaidx_t *idx = bwa_idx_load(prefix, BWA_IDX_ALL);
	int64_t expect_mem;
	long n_pos;
	int i;
	if (idx == 0) { fprintf(stderr, "cannot load the index %s\n", prefix); exit(2); }
	check("fixture", "table built at load", 0, 0, 1, idx->bns->rid_bucket != 0);
	check_table("fixture", idx->bns);
	check_positions("fixture", idx->bns);
	check_boundaries("fixture", idx->bns);
	check_random_intervals("fixture", idx->bns, 100000);
	n_pos = idx->bns->l_pac + 6;
	/* BWA's shared-memory layout: BWA's bwt_t (up to sa32, the dense suffix-array fields), the BWT,
	 * the SA, BWA's bntseq_t (up to rid_bucket), the holes, the annotations, their names and
	 * comments, the packed sequence */
	expect_mem = offsetof(bwt_t, sa32) + idx->bwt->bwt_size * 4 + idx->bwt->n_sa * sizeof(bwtint_t)
		+ offsetof(bntseq_t, rid_bucket) + idx->bns->n_holes * sizeof(bntamb1_t) + idx->bns->n_seqs * sizeof(bntann1_t)
		+ idx->bns->l_pac / 4 + 1;
	for (i = 0; i < idx->bns->n_seqs; ++i)
		expect_mem += strlen(idx->bns->anns[i].name) + strlen(idx->bns->anns[i].anno) + 2;
	bwa_idx2mem(idx); // asserts that the layout it wrote is the layout bwa_mem2idx() reads
	check("shared-memory layout", "segment size", 0, 0, (long)expect_mem, (long)idx->l_mem);
	check("shared-memory layout", "no table", 0, 0, 0, idx->bns->rid_bucket != 0);
	check_positions("shared-memory layout", idx->bns);
	bwa_idx_destroy(idx);
	return n_pos;
}

int main(int argc, char *argv[])
{
	long n_fix, n_syn;
	if (argc < 2) { fprintf(stderr, "Usage: %s <index prefix>\n", argv[0]); return 2; }
	bwa_verbose = 1;
	n_fix = part1_and_4(argv[1]);
	n_syn = part2(2000);
	part3();
	printf("test_pos2rid: fixture index %ld positions and its accession boundaries, %ld synthetic positions; %ld comparisons, %ld differ\n",
		n_fix, n_syn, n_cmp, n_bad);
	return n_bad == 0? 0 : 1;
}
