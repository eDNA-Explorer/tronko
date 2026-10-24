/* Unit test of BWA's mate rescue as tronko-assign calls it (bwa_source_files/bwamem_pair.c).
 *
 * mem_sam_pe passes the mate to mem_matesw as upper-case ASCII (worker1 aligns a converted copy
 * of each read), and mem_matesw hands its query to ksw_align2, whose ksw_qinit indexes the 5x5
 * score matrix opt->mat with each query code. With ASCII codes that read 65 to 109 bytes past
 * the matrix, beyond the end of mem_opt_t, into whatever the heap held there.
 *
 * Here mem_opt_t is placed so that it ends exactly at a page with no access rights, so any read
 * past it is a SIGSEGV. Each case runs in a child process:
 *   - guard: ksw_align2 on the guarded matrix with an ASCII query must be killed by SIGSEGV (the
 *     guard works) and with the same query as 2-bit codes must finish;
 *   - rescue: mem_matesw with an ASCII mate, in each orientation (FF, FR, RF, RR) and for a
 *     150-base mate (8-bit kernel) and a 260-base mate (16-bit kernel), must finish, find the
 *     mate where it was taken from the reference, with the full match score, as upstream BWA
 *     does for the same mate given as 2-bit codes.
 * Before the fix the forward orientations (FF, RR) read past the matrix (SIGSEGV here) and the
 * reverse ones (FR, RF) aligned a query of N only and never rescued.
 *
 * Build and run: make test (from tronko-assign/). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include "../../bwa_source_files/bwamem.h"
#include "../../bwa_source_files/bntseq.h"
#include "../../bwa_source_files/ksw.h"
#include "../../global.h"

queryMatPaired *pairedQueryMat;   /* Tronko globals bwamem.c refers to */
queryMatSingle *singleQueryMat;

int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *ma);

#define L 4000                    /* reference length */
static uint8_t ref2[L];           /* reference as 2-bit codes */
static uint8_t *pac;
static bntseq_t bns;
static bntann1_t ann;
static mem_opt_t *opt;            /* ends at a PROT_NONE page */

static void *guarded(size_t size)
{
	long pg = sysconf(_SC_PAGESIZE);
	uint8_t *p = mmap(0, 2 * pg, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (p == MAP_FAILED || mprotect(p + pg, pg, PROT_NONE) != 0) { perror("guard page"); exit(2); }
	return p + pg - size;
}

static void setup(void)
{
	uint32_t x = 12345;
	int i;
	mem_opt_t *o = mem_opt_init();
	for (i = 0; i < L; i++) { x = x * 1103515245u + 12345u; ref2[i] = (x >> 16) & 3; }
	pac = calloc(L / 4 + 1, 1);
	for (i = 0; i < L; i++) pac[i >> 2] |= ref2[i] << ((~i & 3) << 1);
	memset(&ann, 0, sizeof(ann)); ann.offset = 0; ann.len = L; ann.name = "ref";
	memset(&bns, 0, sizeof(bns)); bns.l_pac = L; bns.n_seqs = 1; bns.anns = &ann;
	opt = guarded(sizeof(mem_opt_t));
	*opt = *o;
	free(o);
}

/* Runs f(arg) in a child; returns the child's exit status, or 128 + signal. */
static int in_child(int (*f)(int), int arg)
{
	pid_t pid = fork();
	int st;
	if (pid == 0) {
		struct rlimit no_core = { 0, 0 };
		setrlimit(RLIMIT_CORE, &no_core);   /* the guard check crashes on purpose */
		_exit(f(arg));
	}
	if (pid < 0 || waitpid(pid, &st, 0) < 0) { perror("fork"); exit(2); }
	return WIFSIGNALED(st)? 128 + WTERMSIG(st) : WEXITSTATUS(st);
}

static int guard_ksw(int ascii)
{
	uint8_t q[100], t[100];
	int i;
	for (i = 0; i < 100; i++) { t[i] = ref2[500 + i]; q[i] = ascii? "ACGT"[ref2[500 + i]] : ref2[500 + i]; }
	kswr_t r = ksw_align2(100, q, 100, t, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
	return r.score == 100? 0 : 1;
}

/* orientation r as in mem_matesw: 0 FF, 1 FR, 2 RF, 3 RR; case index = r * 2 + (long mate) */
static int rescue(int c)
{
	int r = c >> 1, l_ms = (c & 1)? 260 : 150, is_rev = (r >> 1 != (r & 1)), is_larger = !(r >> 1), i;
	int64_t pos = is_larger? 2400 : 1000;   /* where the mate comes from, forward strand */
	char *ms = malloc(l_ms);
	mem_pestat_t pes[4];
	mem_alnreg_t a;
	mem_alnreg_v ma = { 0, 0, 0 };
	int64_t rb, re;
	/* the mate as tronko-assign passes it: ASCII; reverse complemented in the FR and RF cases */
	for (i = 0; i < l_ms; i++) {
		int b = ref2[pos + i];
		if (is_rev) ms[l_ms - 1 - i] = "TGCA"[b]; else ms[i] = "ACGT"[b];
	}
	memset(pes, 0, sizeof(pes));
	for (i = 0; i < 4; i++) pes[i].failed = (i != r);
	pes[r].low = 100; pes[r].high = 1500;
	memset(&a, 0, sizeof(a));
	a.rid = 0; a.rb = 2000 - (is_larger? 400 : 0) + (is_larger? 0 : 400); a.re = a.rb + 100;
	a.qb = 0; a.qe = 100; a.score = 100;
	/* expected position of the rescued mate: forward-strand [pos, pos+l_ms), or its image on the reverse strand */
	rb = is_rev? (L << 1) - (pos + l_ms) : pos;
	re = is_rev? (L << 1) - pos : pos + l_ms;
	int n = mem_matesw(opt, &bns, pac, pes, &a, l_ms, (uint8_t *)ms, &ma);
	int ok = n == 1 && ma.n == 1 && ma.a[0].score == l_ms && ma.a[0].rb == rb && ma.a[0].re == re
		&& ma.a[0].qb == 0 && ma.a[0].qe == l_ms;
	if (!ok) fprintf(stderr, "rescue r=%d l_ms=%d: n=%d ma.n=%zu score=%d rb=%ld re=%ld (expected %ld %ld)\n", r, l_ms, n,
		ma.n, ma.n? ma.a[0].score : -1, ma.n? (long)ma.a[0].rb : -1L, ma.n? (long)ma.a[0].re : -1L, (long)rb, (long)re);
	return ok? 0 : 1;
}

int main(void)
{
	static const char *orient[4] = { "FF", "FR", "RF", "RR" };
	int fail = 0, c, st;
	setup();
	st = in_child(guard_ksw, 1);
	if (st != 128 + SIGSEGV) { fprintf(stderr, "FAIL guard: ksw_align2 with an ASCII query was not stopped by the guard page (status %d)\n", st); fail++; }
	st = in_child(guard_ksw, 0);
	if (st != 0) { fprintf(stderr, "FAIL guard: ksw_align2 with a 2-bit query, status %d\n", st); fail++; }
	for (c = 0; c < 8; c++) {
		st = in_child(rescue, c);
		if (st != 0) {
			fprintf(stderr, "FAIL rescue %s, %d-base mate: %s\n", orient[c >> 1], (c & 1)? 260 : 150,
				st == 128 + SIGSEGV? "read past the score matrix (SIGSEGV at the guard page)" : "wrong or missing rescue");
			fail++;
		}
	}
	printf("test_mate_rescue: 10 checks, %d failures\n", fail);
	return fail? 1 : 0;
}
