/* Unit test of the SAM SEQ column written by mem_aln2sam (bwa_source_files/bwamem.c).
 *
 * The read handed to mem_aln2sam is ASCII (worker1 aligns a converted copy, mem_align1), and the
 * SAM text is later parsed as a C string (fastmap.c: strtok_r on "\n", then sscanf). The property
 * checked: for every byte a read can hold, on both strands, for a primary record and a clipped
 * supplementary record,
 *   - the SAM text holds no NUL before its end (strlen == length) and no TAB or newline inside SEQ;
 *   - SEQ is the read's bases as A, C, G, T or N (reverse complemented on the reverse strand),
 *     every other byte printing as N, the same length as before;
 *   - parsing the text as fastmap.c does yields both records, with the fields it reads.
 * Before the fix, SEQ indexed "ACGTN"/"TGCAN" with the ASCII code, and on builds where that read a
 * 0 byte (aarch64: forward-strand T) the text ended inside the first record.
 *
 * Build and run: make test (from tronko-assign/). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../../bwa_source_files/bwamem.h"
#include "../../bwa_source_files/bntseq.h"
#include "../../bwa_source_files/kstring.h"
#include "../../global.h"

/* The fork's mem_aln2sam, which no header declares (bwamem.c). */
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, int concordant, int seq_index, int startline, int paired);

/* Tronko globals that mem_aln2sam reverse-complements a copy of the read in (reverse strand). */
queryMatPaired *pairedQueryMat;
queryMatSingle *singleQueryMat;

static int failures, checks;
#define CHECK(cond, ...) do { checks++; if (!(cond)) { failures++; if (failures <= 20) { fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } } } while (0)

static char expect_fwd(unsigned char c)
{
	switch (c) {
	case 'A': case 'a': return 'A';
	case 'C': case 'c': return 'C';
	case 'G': case 'g': return 'G';
	case 'T': case 't': return 'T';
	default: return 'N';
	}
}
static char expect_rev(unsigned char c)
{
	switch (expect_fwd(c)) {
	case 'A': return 'T';
	case 'C': return 'G';
	case 'G': return 'C';
	case 'T': return 'A';
	default: return 'N';
	}
}

/* One alignment of the read to accession rid: a primary (M over the whole read) or a supplementary
 * (hard-clipped by `clip` at both ends). The MD string follows the CIGAR, as mem_reg2aln stores it. */
static mem_aln_t make_aln(int rid, int is_rev, int l_seq, int clip, int supplementary)
{
	mem_aln_t a;
	memset(&a, 0, sizeof(a));
	a.rid = rid;
	a.pos = 99 + rid;
	a.is_rev = is_rev;
	a.mapq = 60;
	a.NM = 0;
	a.score = l_seq - 2 * clip;
	a.sub = -1;
	a.flag = supplementary? 0x800 : 0;
	a.n_cigar = clip? 3 : 1;
	a.cigar = calloc(a.n_cigar + 4, sizeof(uint32_t)); /* room for the MD string after the CIGAR */
	if (clip) {
		a.cigar[0] = (uint32_t)clip << 4 | 4;
		a.cigar[1] = (uint32_t)(l_seq - 2 * clip) << 4 | 0;
		a.cigar[2] = (uint32_t)clip << 4 | 4;
	} else a.cigar[0] = (uint32_t)l_seq << 4 | 0;
	strcpy((char *)(a.cigar + a.n_cigar), "5");
	return a;
}

/* Field k (0-based) of a TAB-separated record, copied into out. */
static void field(const char *rec, int k, char *out, size_t n)
{
	const char *p = rec, *e;
	while (k-- > 0 && p) { p = strchr(p, '\t'); if (p) p++; }
	if (!p) { out[0] = 0; return; }
	e = p;
	while (*e && *e != '\t' && *e != '\n') e++;
	if ((size_t)(e - p) >= n) e = p + n - 1;
	memcpy(out, p, e - p);
	out[e - p] = 0;
}

static void run_case(const mem_opt_t *opt, const bntseq_t *bns, const char *read, int is_rev, const char *label)
{
	int l = (int)strlen(read), clip = l >= 8? 2 : 0, i, nrec;
	char *seq = strdup(read), *copy = strdup(read);
	char name[] = "read_1", buf[4096], exp[4096], *rest, *tok;
	char *qmat[1] = { copy }, *qname[1] = { name };
	queryMatSingle qs = { qmat, qname };
	bseq1_t s;
	mem_aln_t list[2];
	kstring_t str = { 0, 0, 0 };

	singleQueryMat = &qs;
	memset(&s, 0, sizeof(s));
	s.name = name; s.seq = seq; s.l_seq = l; s.id = 0;
	list[0] = make_aln(0, is_rev, l, 0, 0);
	list[1] = make_aln(1, is_rev, l, clip, 1);
	/* as mem_reg2sam writes a read's records: one call per record into the same string */
	mem_aln2sam(opt, bns, &str, &s, 2, list, 0, 0, 0, 0, 0, 0);
	mem_aln2sam(opt, bns, &str, &s, 2, list, 1, 0, 0, 0, 0, 0);

	CHECK(strlen(str.s) == str.l, "%s: SAM text has a NUL at %zu of %zu", label, strlen(str.s), str.l);
	CHECK(strcmp(s.seq, read) == 0, "%s: s->seq modified", label);
	/* SEQ of each record */
	for (int r = 0; r < 2; r++) {
		const char *rec = str.s;
		int qb = r? clip : 0, qe = r? l - clip : l, k = 0;
		if (r) { rec = memchr(str.s, '\n', str.l); rec = rec? rec + 1 : str.s + str.l; }
		field(rec, 9, buf, sizeof(buf));
		if (!is_rev) for (i = qb; i < qe; i++) exp[k++] = expect_fwd((unsigned char)read[i]);
		else for (i = qe - 1; i >= qb; i--) exp[k++] = expect_rev((unsigned char)read[i]);
		exp[k] = 0;
		CHECK(strcmp(buf, exp) == 0, "%s: record %d SEQ '%s', expected '%s'", label, r, buf, exp);
	}
	/* the parse of fastmap.c */
	nrec = 0;
	rest = str.s;
	while ((tok = strtok_r(rest, "\n", &rest))) {
		char readname[64], read1[64], read2[64], cigar[64];
		int flag = -1, pos = -1, ok;
		ok = sscanf(tok, "%s %d %s %d %*d %s %s %*d %*d %*s %*s %*s %*s %*s %*s %*[^.,;]", readname, &flag, read1, &pos, cigar, read2);
		CHECK(ok >= 6, "%s: record %d parsed %d fields", label, nrec, ok);
		CHECK(strcmp(readname, "read_1") == 0, "%s: record %d QNAME %s", label, nrec, readname);
		CHECK(strcmp(read1, nrec? "acc_B" : "acc_A") == 0, "%s: record %d RNAME %s", label, nrec, read1);
		CHECK(pos == 100 + nrec, "%s: record %d POS %d", label, nrec, pos);
		CHECK(flag == ((is_rev? 0x10 : 0) | (nrec? 0x800 : 0)), "%s: record %d FLAG %d", label, nrec, flag);
		nrec++;
	}
	CHECK(nrec == 2, "%s: %d records parsed, expected 2", label, nrec);
	free(list[0].cigar); free(list[1].cigar);
	free(str.s); free(seq); free(copy);
}

int main(void)
{
	mem_opt_t *opt = mem_opt_init();
	bntann1_t anns[2];
	bntseq_t bns;
	char read[64], label[64];
	int c, strand;

	memset(anns, 0, sizeof(anns));
	anns[0].name = "acc_A"; anns[1].name = "acc_B";
	memset(&bns, 0, sizeof(bns));
	bns.anns = anns; bns.n_seqs = 2;

	/* every byte value a read string can hold, as a run of 12, on both strands */
	for (strand = 0; strand < 2; strand++) {
		for (c = 1; c < 256; c++) {
			memset(read, c, 12); read[12] = 0;
			snprintf(label, sizeof(label), "byte %d %s", c, strand? "reverse" : "forward");
			run_case(opt, &bns, read, strand, label);
		}
		/* mixed reads: each base and IUPAC code, and a gap */
		run_case(opt, &bns, "ACGTACGTNNACGTRYKMSWBDHVAC-GT", strand, strand? "mixed reverse" : "mixed forward");
		run_case(opt, &bns, "TTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGAAAAACCCCC", strand, strand? "runs reverse" : "runs forward");
	}
	free(opt);
	printf("test_sam_seq: %d checks, %d failures\n", checks, failures);
	return failures? 1 : 0;
}
