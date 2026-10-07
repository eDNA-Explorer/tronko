/* The candidate search of tronko-assign: BWA-MEM3 v0.14.0, run unchanged as a separate program,
 * once per batch (design (a) of the aligner plan: run per batch, index in shared memory).
 *
 * For each batch, bwamem3_run_batch
 *   1. writes the reads exactly as Tronko used to hand them to its vendored BWA (the sequences
 *      of singleQueryMat/pairedQueryMat as they stand: read 2 and -v reads already reverse-
 *      complemented by Tronko, reverse-complement truncation included) to FASTA files in the
 *      process's temporary directory, each read named by its index in the batch;
 *   2. runs `bwa-mem3 mem` once on them with the options of mem3_mem_options below, at Tronko's
 *      thread count;
 *   3. reads the SAM, groups the records by read (the vendored BWA's seqs[i].sam: read 1 and
 *      read 2 of a pair are two reads), checks that every read has records and that they come in
 *      input order, and puts back each read's name as the vendored BWA wrote it (the name after
 *      trim_readno), so the SAM parse sees the same text in the fields it reads;
 *   4. reverse-complements Tronko's copy of a read when the vendored BWA's SAM writer would have
 *      (rule H8: single reads after an odd number of records with flag 0x10; for a pair, read 1
 *      after an odd number of read-1 records with 0x10, read 2 after an odd number of read-1
 *      records with 0x20);
 *   5. runs the SAM parse (sam_parse.c).
 *
 * The index is BWA-MEM3's, beside the -a FASTA (<fasta>.amb .ann .pac .bwt.2bit.64). The process
 * reaches it through symbolic links named tronko-<pid>-<basename> in its temporary directory, so
 * the shared-memory segment it stages (BWA-MEM3 keys segments by the prefix's basename) cannot be
 * confused with a segment another process staged from a different file of the same name, nor a
 * stale one. */
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "bwamem3.h"
#include "sam_parse.h"

extern char **environ;

/* `bwa-mem3 mem` options besides -t and the files. Each is pinned for a reason:
 *   --compat=bwa-mem      records identical to upstream BWA 0.7.19 (= the fork with F1 removed)
 *   -K 100000000          the batch is one BWA chunk: one insert-size estimate per batch, and
 *                         tie-break ids from 0 at the batch's first read (each run starts
 *                         n_processed at 0), as the vendored BWA did per call
 *   --dedup-reads off     the three adaptive duplicate modes default to 'auto', which measures
 *   --dedup off           at run time, latches and re-probes; BWA-MEM3 states all arms give the
 *   --ks-dedup off        same records. Pinned so every run takes one code path; 'off' is the
 *                         path closest to upstream BWA
 *   --cohort-slices 0     read the batch in one piece (the default reads a first batch above
 *                         16,000,000 bases as a ramp of slices, stated output-neutral)
 *   -v 1                  errors only (BWA-MEM3 still prints its insert-size lines and a
 *                         start-up report; bwamem3_run_batch passes on the insert-size lines)
 * Every BWA-MEM3 environment override (BWAMEM3_*, BWA_MEM3_*, BWA3_*) is removed from the
 * aligner's environment, since some of them win over these options. */
static const char *const mem3_mem_options[] = {
	"--compat=bwa-mem", "-K", "100000000",
	"--dedup-reads", "off", "--dedup", "off", "--ks-dedup", "off",
	"--cohort-slices", "0", "-v", "1",
};
#define MEM3_CHUNK_BASES 100000000L
static const char *const mem3_index_suffixes[] = { "amb", "ann", "pac", "bwt.2bit.64" };

static struct {
	int ready;
	char bin[PATH_MAX];
	char fasta[PATH_MAX];      /* -a, as given */
	char dir[PATH_MAX];        /* the process's temporary directory */
	char prefix[PATH_MAX];     /* dir/tronko-<pid>-<basename>: the index prefix bwa-mem3 is given */
	char threads[16];
	int staged;                /* this process staged prefix in shared memory */
	char **envp;               /* environ without BWA-MEM3's overrides */
	long batch;
} A;

static void die(const char *fmt, const char *a, const char *b)
{
	fprintf(stderr, "tronko-assign: bwa-mem3: ");
	fprintf(stderr, fmt, a ? a : "", b ? b : "");
	fprintf(stderr, "\n");
	exit(1);
}

static void build_envp(void)
{
	size_t n = 0, k = 0;
	char **e;
	for (e = environ; *e; e++) n++;
	A.envp = malloc((n + 1) * sizeof(char *));
	for (e = environ; *e; e++) {
		if (strncmp(*e, "BWAMEM3_", 8) == 0 || strncmp(*e, "BWA_MEM3_", 9) == 0 || strncmp(*e, "BWA3_", 5) == 0)
			continue;
		A.envp[k++] = *e;
	}
	A.envp[k] = 0;
}

/* Runs argv (argv[0] searched on PATH when it has no slash) with A.envp. Its standard output goes
 * to out_path (truncated) or, when out_path is NULL, to our standard error, so Tronko's standard
 * output stays clean; its standard error goes to err_path (truncated) or, when NULL, to ours.
 * Returns the exit status, or -1 if it could not run or was killed. */
static int run_to(char *const argv[], const char *out_path, const char *err_path)
{
	posix_spawn_file_actions_t fa;
	pid_t pid;
	int status, rc;
	posix_spawn_file_actions_init(&fa);
	if (err_path)
		posix_spawn_file_actions_addopen(&fa, STDERR_FILENO, err_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	if (out_path)
		posix_spawn_file_actions_addopen(&fa, STDOUT_FILENO, out_path, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	else
		posix_spawn_file_actions_adddup2(&fa, STDERR_FILENO, STDOUT_FILENO);
	fflush(stdout);
	fflush(stderr);
	rc = posix_spawnp(&pid, argv[0], &fa, 0, argv, A.envp);
	posix_spawn_file_actions_destroy(&fa);
	if (rc != 0) {
		fprintf(stderr, "tronko-assign: bwa-mem3: cannot run %s: %s\n", argv[0], strerror(rc));
		return -1;
	}
	while (waitpid(pid, &status, 0) < 0)
		if (errno != EINTR) return -1;
	return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

static int run(char *const argv[], const char *out_path)
{
	return run_to(argv, out_path, 0);
}

static char *slurp(const char *path, size_t *len)
{
	FILE *f = fopen(path, "rb");
	char *s;
	long n;
	if (!f) return 0;
	fseek(f, 0, SEEK_END);
	n = ftell(f);
	fseek(f, 0, SEEK_SET);
	s = malloc(n + 1);
	if (n > 0 && fread(s, 1, n, f) != (size_t)n) { fclose(f); free(s); return 0; }
	fclose(f);
	s[n] = 0;
	if (len) *len = n;
	return s;
}

static void find_binary(const char *bin_option)
{
	const char *env = getenv("TRONKO_BWA_MEM3");
	char self[PATH_MAX];
	ssize_t n;
	if (bin_option && bin_option[0]) {
		snprintf(A.bin, sizeof A.bin, "%s", bin_option);
		return;
	}
	if (env && env[0]) {
		snprintf(A.bin, sizeof A.bin, "%s", env);
		return;
	}
	n = readlink("/proc/self/exe", self, sizeof self - 1);
	if (n > 0) {
		char *slash;
		self[n] = 0;
		slash = strrchr(self, '/');
		if (slash) {
			slash[1] = 0;
			if (strlen(self) + strlen("bwa-mem3") < sizeof A.bin) {
				snprintf(A.bin, sizeof A.bin, "%sbwa-mem3", self);
				if (access(A.bin, X_OK) == 0) return;
			}
		}
	}
	snprintf(A.bin, sizeof A.bin, "bwa-mem3");
}

static void check_version(void)
{
	char out[PATH_MAX], *s, *nl;
	char *argv[] = { A.bin, "version", 0 };
	snprintf(out, sizeof out, "%s/version.txt", A.dir);
	if (run(argv, out) != 0) die("`%s version` failed; set --bwa-mem3 or TRONKO_BWA_MEM3 to BWA-MEM3 " BWAMEM3_VERSION "%s", A.bin, "");
	s = slurp(out, 0);
	if (!s) die("cannot read %s%s", out, "");
	nl = strchr(s, '\n');
	if (nl) *nl = 0;
	if (strcmp(s, BWAMEM3_VERSION) != 0)
		die("%s reports version '%s'; tronko-assign needs BWA-MEM3 " BWAMEM3_VERSION, A.bin, s);
	free(s);
	unlink(out);
}

static void index_file(char *dst, size_t size, const char *prefix, const char *suffix)
{
	if ((size_t)snprintf(dst, size, "%s.%s", prefix, suffix) >= size) die("path too long: %s.%s", prefix, suffix);
}

void bwamem3_init(const char *bin_option, const char *fasta, int threads, int skip_build, int use_shm)
{
	char tmpl[PATH_MAX], real[PATH_MAX], src[PATH_MAX], dst[PATH_MAX];
	const char *tmp = getenv("TMPDIR"), *base;
	size_t i;
	if (A.ready) return;
	memset(&A, 0, sizeof A);
	build_envp();
	snprintf(A.fasta, sizeof A.fasta, "%s", fasta);
	snprintf(A.threads, sizeof A.threads, "%d", threads > 1 ? threads : 1);
	snprintf(tmpl, sizeof tmpl, "%s/tronko-assign.XXXXXX", tmp && tmp[0] ? tmp : "/tmp");
	if (strlen(tmpl) > PATH_MAX - 256) die("TMPDIR is too long: %s%s", tmp, "");
	if (!mkdtemp(tmpl)) die("cannot make a temporary directory %s: %s", tmpl, strerror(errno));
	snprintf(A.dir, sizeof A.dir, "%s", tmpl);
	A.ready = 1;
	atexit(bwamem3_finish);

	find_binary(bin_option);
	check_version();

	if (!skip_build) {
		char *argv[] = { A.bin, "index", A.fasta, 0 };
		fprintf(stderr, "tronko-assign: building the BWA-MEM3 index of %s\n", A.fasta);
		if (run(argv, 0) != 0) die("`%s index %s` failed", A.bin, A.fasta);
	}
	/* the index files beside the -a path as given (as the vendored BWA read them), made absolute */
	if (fasta[0] == '/') snprintf(real, sizeof real, "%s", fasta);
	else {
		char cwd[PATH_MAX];
		if (!getcwd(cwd, sizeof cwd)) die("cannot read the working directory: %s%s", strerror(errno), "");
		if ((size_t)snprintf(real, sizeof real, "%s/%s", cwd, fasta) >= sizeof real) die("path too long: %s/%s", cwd, fasta);
	}
	base = strrchr(real, '/');
	base = base ? base + 1 : real;
	if ((size_t)snprintf(A.prefix, sizeof A.prefix, "%s/tronko-%ld-%s", A.dir, (long)getpid(), base) >= sizeof A.prefix)
		die("path too long: %s/%s", A.dir, base);
	for (i = 0; i < sizeof mem3_index_suffixes / sizeof *mem3_index_suffixes; i++) {
		index_file(src, sizeof src, real, mem3_index_suffixes[i]);
		index_file(dst, sizeof dst, A.prefix, mem3_index_suffixes[i]);
		if (access(src, R_OK) != 0)
			die("the BWA-MEM3 index file %s is missing (build it with `bwa-mem3 index <fasta>`, or run without -6)%s", src, "");
		if (symlink(src, dst) != 0) die("cannot link %s: %s", dst, strerror(errno));
	}
	index_file(src, sizeof src, real, "alt");
	if (access(src, R_OK) == 0) {
		index_file(dst, sizeof dst, A.prefix, "alt");
		if (symlink(src, dst) != 0) die("cannot link %s: %s", dst, strerror(errno));
	}
	if (use_shm) {
		char *argv[] = { A.bin, "shm", A.prefix, 0 };
		if (run(argv, 0) != 0) die("`%s shm %s` failed (is /dev/shm large enough?)", A.bin, A.prefix);
		A.staged = 1;
	}
}

void bwamem3_finish(void)
{
	char path[PATH_MAX];
	size_t i;
	if (!A.ready) return;
	A.ready = 0;
	if (A.staged) {
		/* BWA-MEM3 can drop only every staged index at once (shm -d), so drop ours only when it
		 * is the only one listed; otherwise leave it and say so. */
		char *list_argv[] = { A.bin, "shm", "-l", 0 }, *drop_argv[] = { A.bin, "shm", "-d", 0 };
		char *list, *line, *save = 0;
		const char *name = strrchr(A.prefix, '/') + 1;
		int others = 0, ours = 0;
		snprintf(path, sizeof path, "%s/shm-list.txt", A.dir);
		if (run(list_argv, path) == 0 && (list = slurp(path, 0)) != 0) {
			for (line = strtok_r(list, "\n", &save); line; line = strtok_r(0, "\n", &save)) {
				size_t l = strcspn(line, "\t");
				if (l == strlen(name) && strncmp(line, name, l) == 0) ours = 1;
				else if (l > 0) others = 1;
			}
			free(list);
			if (ours && !others) run(drop_argv, 0);
			else if (ours) fprintf(stderr, "tronko-assign: bwa-mem3: left %s staged in shared memory: other indices are staged too\n", name);
		}
		unlink(path);
		A.staged = 0;
	}
	for (i = 0; i < sizeof mem3_index_suffixes / sizeof *mem3_index_suffixes; i++) {
		snprintf(path, sizeof path, "%s.%s", A.prefix, mem3_index_suffixes[i]);
		unlink(path);
	}
	snprintf(path, sizeof path, "%s.alt", A.prefix); unlink(path);
	snprintf(path, sizeof path, "%s/batch.sam", A.dir); unlink(path);
	snprintf(path, sizeof path, "%s/batch.log", A.dir); unlink(path);
	snprintf(path, sizeof path, "%s/batch_1.fa", A.dir); unlink(path);
	snprintf(path, sizeof path, "%s/batch_2.fa", A.dir); unlink(path);
	snprintf(path, sizeof path, "%s/batch.fa", A.dir); unlink(path);
	rmdir(A.dir);
}

/* The read name the vendored BWA gave the aligner: Tronko's name less a trailing "/<digit>"
 * (BWA's trim_readno). */
static size_t trimmed_name_length(const char *name)
{
	size_t l = strlen(name);
	if (l > 2 && name[l-2] == '/' && isdigit((unsigned char)name[l-1])) l -= 2;
	return l;
}

static char complement_base(char b)
{ /* as the vendored BWA's complement(): other bytes are kept */
	switch (b) {
	case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C';
	case 'a': return 't'; case 't': return 'a'; case 'c': return 'g'; case 'g': return 'c';
	default: return b;
	}
}

static void reverse_complement_in_place(char *s)
{ /* as the vendored BWA's reverse_complement() */
	int len = strlen(s), i;
	for (i = 0; i < len / 2; i++) {
		char t = s[i];
		s[i] = complement_base(s[len - 1 - i]);
		s[len - 1 - i] = complement_base(t);
	}
	if (len % 2 == 1) s[len / 2] = complement_base(s[len / 2]);
}

static long write_fasta(const char *path, int n, char **seqs)
{
	FILE *f = fopen(path, "w");
	long bases = 0;
	int i;
	if (!f) die("cannot write %s: %s", path, strerror(errno));
	for (i = 0; i < n; i++) {
		fprintf(f, ">%d\n%s\n", i, seqs[i]);
		bases += strlen(seqs[i]);
	}
	if (fclose(f) != 0) die("cannot write %s: %s", path, strerror(errno));
	return bases;
}

void bwamem3_run_batch(int n, int paired, bwaMatches *results, int concordant, int ntree,
                       int max_readname_length, int max_acc_name)
{
	char sam_path[PATH_MAX], log_path[PATH_MAX], f1[PATH_MAX], f2[PATH_MAX], batch_no[32];
	char *argv[32], *sam, *p, **text, *log;
	int rc;
	size_t *len, *cap, sam_len;
	int n_seqs = paired ? 2 * n : n, argc = 0, i, last = -1;
	long bases;
	int *p10, *p20;
	size_t k;

	if (!A.ready) die("bwamem3_run_batch before bwamem3_init%s%s", 0, 0);
	if (n <= 0) return;
	A.batch++;
	snprintf(sam_path, sizeof sam_path, "%s/batch.sam", A.dir);
	if (paired) {
		snprintf(f1, sizeof f1, "%s/batch_1.fa", A.dir);
		snprintf(f2, sizeof f2, "%s/batch_2.fa", A.dir);
		bases = write_fasta(f1, n, pairedQueryMat->query1Mat);
		bases += write_fasta(f2, n, pairedQueryMat->query2Mat);
	} else {
		snprintf(f1, sizeof f1, "%s/batch.fa", A.dir);
		bases = write_fasta(f1, n, singleQueryMat->queryMat);
	}
	if (bases >= MEM3_CHUNK_BASES) {
		snprintf(batch_no, sizeof batch_no, "%ld", bases);
		die("a batch of %s bases does not fit one BWA-MEM3 chunk; lower -L%s", batch_no, "");
	}

	argv[argc++] = A.bin;
	argv[argc++] = "mem";
	for (k = 0; k < sizeof mem3_mem_options / sizeof *mem3_mem_options; k++) argv[argc++] = (char *)mem3_mem_options[k];
	argv[argc++] = "-t";
	argv[argc++] = A.threads;
	argv[argc++] = "-o";
	argv[argc++] = sam_path;
	argv[argc++] = A.prefix;
	argv[argc++] = f1;
	if (paired) argv[argc++] = f2;
	argv[argc] = 0;
	/* BWA-MEM3's messages go to a file: its insert-size lines ("[PE]", the counterpart of the
	 * vendored BWA's mem_pestat lines) are passed on, everything on failure */
	snprintf(log_path, sizeof log_path, "%s/batch.log", A.dir);
	rc = run_to(argv, 0, log_path);
	if ((log = slurp(log_path, 0)) != 0) {
		char *line, *save = 0;
		if (rc != 0) fputs(log, stderr);
		else for (line = strtok_r(log, "\n", &save); line; line = strtok_r(0, "\n", &save))
			if (strstr(line, "[PE]")) fprintf(stderr, "%s\n", line);
		free(log);
	}
	if (rc != 0) {
		snprintf(batch_no, sizeof batch_no, "%ld", A.batch);
		die("bwa-mem3 mem failed on batch %s (%s)", batch_no, A.bin);
	}

	/* Group the records by read, as the vendored BWA's seqs[i].sam */
	sam = slurp(sam_path, &sam_len);
	if (!sam) die("cannot read %s%s", sam_path, "");
	text = calloc(n_seqs, sizeof(char *));
	len = calloc(n_seqs, sizeof(size_t));
	cap = calloc(n_seqs, sizeof(size_t));
	p10 = calloc(n_seqs, sizeof(int));
	p20 = calloc(n_seqs, sizeof(int));
	for (p = sam; *p; ) {
		char *eol = strchr(p, '\n'), *tab, *end;
		const char *name;
		size_t rest, nl;
		long q;
		int flag, s;
		if (!eol) eol = p + strlen(p);
		if (*p == '@' && last < 0) { p = *eol ? eol + 1 : eol; continue; }
		tab = memchr(p, '\t', eol - p);
		q = strtol(p, &end, 10);
		if (!tab || end != tab || q < 0 || q >= n) die("unexpected SAM record in %s: %.60s", sam_path, p);
		flag = (int)strtol(tab + 1, &end, 10);
		if (end == tab + 1 || *end != '\t') die("unexpected SAM record in %s: %.60s", sam_path, p);
		s = paired ? 2 * (int)q + ((flag & 0x80) ? 1 : 0) : (int)q;
		if (s < last) die("SAM records out of input order in %s: %.60s", sam_path, p);
		last = s;
		if (flag & 0x10) p10[s]++;
		if (flag & 0x20) p20[s]++;
		name = paired ? pairedQueryMat->forward_name[q] : singleQueryMat->name[q];
		nl = trimmed_name_length(name);
		rest = eol - tab;          /* from the tab after QNAME to the end of the line */
		if (len[s] + nl + rest + 2 > cap[s]) {
			cap[s] = 2 * (len[s] + nl + rest + 2);
			text[s] = realloc(text[s], cap[s]);
		}
		memcpy(text[s] + len[s], name, nl);
		memcpy(text[s] + len[s] + nl, tab, rest);
		len[s] += nl + rest;
		text[s][len[s]++] = '\n';
		text[s][len[s]] = 0;
		p = *eol ? eol + 1 : eol;
	}
	free(sam);
	for (i = 0; i < n_seqs; i++)
		if (!text[i]) {
			snprintf(batch_no, sizeof batch_no, "%d", i);
			die("bwa-mem3 gave no record for read %s of the batch (%s)", batch_no, sam_path);
		}

	/* Rule H8: Tronko's copy of the read, reoriented as the vendored BWA's SAM writer left it */
	for (i = 0; i < n; i++) {
		if (!paired) {
			if (p10[i] & 1) reverse_complement_in_place(singleQueryMat->queryMat[i]);
		} else {
			if (p10[2*i] & 1) reverse_complement_in_place(pairedQueryMat->query1Mat[i]);
			if (p20[2*i] & 1) reverse_complement_in_place(pairedQueryMat->query2Mat[i]);
		}
	}

	{
		sam_parse_ctx_t ctx = { results, concordant, ntree, 0, paired, max_readname_length, max_acc_name };
		sam_parse_batch(&ctx, n_seqs, text);
	}
	for (i = 0; i < n_seqs; i++) free(text[i]);
	free(text); free(len); free(cap); free(p10); free(p20);
}
