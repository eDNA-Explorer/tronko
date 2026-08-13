#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "hashmap.h"
#include "hashmap_base.h"
#include "minimap2_src/minimap.h"
#include "minimap2_wrapper.h"

extern queryMatPaired *pairedQueryMat;
extern queryMatSingle *singleQueryMat;
extern node **treeArr;
extern int *numspecArr;

static mm_idx_t *g_mm2_idx = NULL;
static char g_mm2_db_path[MAXFILENAME] = "";
static int g_mm2_kmer = 0;
static int g_mm2_window = 0;
static pthread_mutex_t g_mm2_index_mutex = PTHREAD_MUTEX_INITIALIZER;
static int g_mm2_atexit_registered = 0;

static char complement_base(char base)
{
	switch (base) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		default: return base;
	}
}

static void reverse_complement_in_place(char *sequence)
{
	size_t left = 0;
	size_t right = strlen(sequence);
	if (right == 0) return;
	right--;
	while (left < right) {
		char left_base = sequence[left];
		sequence[left] = complement_base(sequence[right]);
		sequence[right] = complement_base(left_base);
		left++;
		right--;
	}
	if (left == right) sequence[left] = complement_base(sequence[left]);
}

/* Release the cached index at normal program exit (registered via atexit on
 * first build). For a CLI run the OS would reclaim it anyway, but this keeps
 * leak checkers clean and is correct if the wrapper is ever reused. */
static void mm2_destroy_index_atexit(void)
{
	pthread_mutex_lock(&g_mm2_index_mutex);
	if (g_mm2_idx != NULL) {
		mm_idx_destroy(g_mm2_idx);
		g_mm2_idx = NULL;
	}
	pthread_mutex_unlock(&g_mm2_index_mutex);
}

static void cigar_to_string(const uint32_t *cigar, int n_cigar, char *buf, int buf_size)
{
	int pos = 0;
	int i;

	for (i = 0; i < n_cigar && pos < buf_size - 12; i++) {
		int op = cigar[i] & 0xf;
		int len = cigar[i] >> 4;
		char op_char = (op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH)
			? 'M'
			: MM_CIGAR_STR[op];
		pos += snprintf(buf + pos, buf_size - pos, "%d%c", len, op_char);
	}
	if (pos >= buf_size) pos = buf_size - 1;
	buf[pos] = '\0';
}

static mm_idx_t *build_index_locked(const char *db_path, int kmer, int window)
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt_dummy;
	mm_idx_reader_t *reader;
	mm_idx_t *idx;

	mm_verbose = 1;
	mm_set_opt(0, &iopt, &mopt_dummy);
	iopt.k = kmer;
	iopt.w = window;
	/* mm_idx_reader_read() is a *batched* API: one call returns at most one
	   mini-batch of reference sequence, and callers are expected to loop until
	   it returns NULL. Tronko keeps a single global index and calls it once, so
	   everything past the first mini-batch was silently dropped -- measured at
	   97k of 218k sequences for V16S-U and 167k of 257k for vert12S. Reads whose
	   true reference landed past the cutoff came back unassigned with no error
	   and no warning. BWA was unaffected because it loads a prebuilt on-disk
	   index spanning the whole reference, which is why this hid for so long.
	   Size one mini-batch to cover the entire reference instead.
	   CAUTION: mm_idx_gen() takes mini_batch_size as an int. A value above
	   INT_MAX wraps negative and the index degenerates to a single sequence
	   (observed with 8e9). 2e9 stays inside int32 and covers the largest
	   reference we ship (CO1_Metazoa, ~300 Mbases). */
	iopt.batch_size = 0x7FFFFFFFFFFFFFFFULL;
	iopt.mini_batch_size = 2000000000;

	reader = mm_idx_reader_open(db_path, &iopt, NULL);
	if (reader == NULL) {
		fprintf(stderr, "[minimap2] ERROR: cannot open reference file '%s'\n", db_path);
		return NULL;
	}

	idx = mm_idx_reader_read(reader, 1);
	if (idx != NULL) {
		/* Cheap guard: at EOF this returns NULL immediately. A non-NULL batch
		   means the reference outgrew one mini-batch and the index is short. */
		mm_idx_t *overflow = mm_idx_reader_read(reader, 1);
		if (overflow != NULL) {
			fprintf(stderr,
				"[minimap2] ERROR: reference '%s' exceeds one index batch; "
				"%u sequences were NOT indexed. Alignments against them will "
				"silently fail.\n", db_path, overflow->n_seq);
			mm_idx_destroy(overflow);
		}
	}
	mm_idx_reader_close(reader);
	if (idx == NULL) {
		fprintf(stderr, "[minimap2] ERROR: failed to build index from '%s'\n", db_path);
		return NULL;
	}

	return idx;
}

/* Serialized lazy build of the shared index: every worker takes the mutex,
 * the first to find no valid cached index builds it, the rest reuse it.
 * Cache-hit hold time is a strcmp + 3 int compares, so contention is
 * negligible; the index is rebuilt only if the path/k/w change. */
static mm_idx_t *get_or_build_index(const char *db_path, int kmer, int window)
{
	mm_idx_t *idx;

	pthread_mutex_lock(&g_mm2_index_mutex);
	if (g_mm2_idx != NULL &&
	    strcmp(g_mm2_db_path, db_path) == 0 &&
	    g_mm2_kmer == kmer &&
	    g_mm2_window == window) {
		idx = g_mm2_idx;
		pthread_mutex_unlock(&g_mm2_index_mutex);
		return idx;
	}

	if (g_mm2_idx != NULL) {
		mm_idx_destroy(g_mm2_idx);
		g_mm2_idx = NULL;
		g_mm2_db_path[0] = '\0';
		g_mm2_kmer = 0;
		g_mm2_window = 0;
	}

	idx = build_index_locked(db_path, kmer, window);
	if (idx != NULL) {
		g_mm2_idx = idx;
		strncpy(g_mm2_db_path, db_path, sizeof(g_mm2_db_path) - 1);
		g_mm2_db_path[sizeof(g_mm2_db_path) - 1] = '\0';
		g_mm2_kmer = kmer;
		g_mm2_window = window;
		if (!g_mm2_atexit_registered) {
			atexit(mm2_destroy_index_atexit);
			g_mm2_atexit_registered = 1;
		}
	}
	pthread_mutex_unlock(&g_mm2_index_mutex);
	return idx;
}

/* Keep only best-scoring records. minimap2 flags a hit "secondary" when it overlaps
   another on the *query*, not when it is worse: 58% of secondaries score exactly equal
   to their primary. Those co-equal hits are what the LCA is meant to consume, so filter
   on score rather than on id != parent. Frees dropped records' extras so the caller's
   cleanup loop stays correct after n_regs shrinks. */
static int mm2_keep_top_scoring(mm_reg1_t *regs, int n)
{
	int i, k = 0, best = 0;
	for (i = 0; i < n; ++i) if (regs[i].score > best) best = regs[i].score;
	for (i = 0; i < n; ++i) {
		if (regs[i].score == best) regs[k++] = regs[i];
		else if (regs[i].p) free(regs[i].p);
	}
	return k;
}

static int add_match(bwaMatches *result, int root, int node_id,
                     const char *cigar_str, int start_pos,
                     int is_concordant, int is_forward,
                     int max_leaf_matches)
{
	int *roots = is_concordant ? result->concordant_matches_roots : result->discordant_matches_roots;
	int *nodes = is_concordant ? result->concordant_matches_nodes : result->discordant_matches_nodes;
	int k, l;

	for (k = 0; k < max_leaf_matches; k++) {
		if (roots[k] == -1) break;
	}
	for (l = 0; l < k; l++) {
		if (roots[l] == root && nodes[l] == node_id) return k;
	}
	if (k >= max_leaf_matches) return k;

	roots[k] = root;
	nodes[k] = node_id;

	if (result->use_portion == 1) {
		if (is_forward) {
			strncpy(result->cigars_forward[k], cigar_str, MAX_CIGAR - 1);
			result->cigars_forward[k][MAX_CIGAR - 1] = '\0';
			result->starts_forward[k] = start_pos;
		} else {
			strncpy(result->cigars_reverse[k], cigar_str, MAX_CIGAR - 1);
			result->cigars_reverse[k][MAX_CIGAR - 1] = '\0';
			result->starts_reverse[k] = start_pos;
		}
	}

	return k + 1;
}

static void mark_unmatched(bwaMatches *result)
{
	result->concordant_matches_roots[0] = -1;
	result->concordant_matches_nodes[0] = -1;
	result->discordant_matches_roots[0] = -1;
	result->discordant_matches_nodes[0] = -1;
}

void run_minimap2(int start, int end, bwaMatches *bwa_results, int concordant,
                  int numberOfTrees, char *databaseFile, int paired,
                  int max_query_length, int max_readname_length,
                  int max_acc_name, int max_leaf_matches,
                  int mm2_kmer, int mm2_window)
{
	int i, j;
	int n_reads = end - start;
	mm_idx_t *mi = get_or_build_index(databaseFile, mm2_kmer, mm2_window);
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	HASHMAP(char, leafMap) map;
	mm_tbuf_t *tbuf;

	(void)max_query_length;
	(void)max_readname_length;
	(void)max_acc_name;

	/* Reorienting a reverse-strand query invalidates the alignment start and CIGAR
	   recorded alongside it, which -e (leaf-portion mode) is the only consumer of.
	   Refuse the combination rather than score against a mislocated reference window. */
	if (n_reads > 0 && bwa_results[0].use_portion == 1) {
		fprintf(stderr, "[minimap2] ERROR: -e (leaf-portion mode) is not supported with "
		        "--aligner minimap2. Re-run with --aligner bwa, or without -e.\n");
		exit(EXIT_FAILURE);
	}

	if (mi == NULL) {
		for (i = 0; i < n_reads; i++) mark_unmatched(&bwa_results[i]);
		return;
	}

	mm_set_opt(0, &iopt, &mopt);
	mm_set_opt("sr", &iopt, &mopt);
	mopt.flag |= MM_F_CIGAR;
	mopt.best_n = max_leaf_matches;
	mm_mapopt_update(&mopt, mi);

	hashmap_init(&map, hashmap_hash_string, strcmp);
	for (i = 0; i < numberOfTrees; i++) {
		for (j = numspecArr[i] - 1; j < 2 * numspecArr[i] - 1; j++) {
			struct leafMap *lm = malloc(sizeof(struct leafMap));
			lm->name = treeArr[i][j].name;
			lm->root = i;
			lm->node = j;
			/* hashmap_put keeps the first entry and returns -EEXIST (without
			   storing lm) if another tree already has a leaf with this name,
			   as can happen when an accession appears in multiple AncestralClust
			   per-cluster trees. Free lm in that case so duplicates don't leak.
			   (Matches the first-wins resolution of the BWA fastmap.c path.) */
			if (hashmap_put(&map, lm->name, lm) != 0) {
				free(lm);
			}
		}
	}

	tbuf = mm_tbuf_init();
	for (i = 0; i < n_reads; i++) {
		char *seq1 = NULL;
		char *seq2 = NULL;
		int len1 = 0;
		int len2 = 0;
		int n_regs1 = 0;
		int n_regs2 = 0;
		mm_reg1_t *regs1 = NULL;
		mm_reg1_t *regs2 = NULL;

		if (paired != 0) {
			seq1 = pairedQueryMat->query1Mat[start + i];
			seq2 = pairedQueryMat->query2Mat[start + i];
			if (seq1) len1 = strlen(seq1);
			if (seq2) len2 = strlen(seq2);
		} else {
			seq1 = singleQueryMat->queryMat[start + i];
			if (seq1) len1 = strlen(seq1);
		}

		if (len1 == 0 && len2 == 0) {
			mark_unmatched(&bwa_results[i]);
			continue;
		}

		if (paired != 0 && len1 > 0 && len2 > 0) {
			int qlens[2] = { len1, len2 };
			const char *seqs[2] = { seq1, seq2 };
			int n_regs[2] = { 0, 0 };
			mm_reg1_t *regs[2] = { NULL, NULL };
			const char *qname = pairedQueryMat->forward_name[start + i];

			/*
			 * BWA parity: pass Tronko's stored mate sequences as-is. The -z option
			 * already reverse-complements read 2 before this point when callers want
			 * the reverse read pre-oriented for Tronko placement.
			 */
			mm_map_frag(mi, 2, qlens, seqs, n_regs, regs, tbuf, &mopt, qname);
			n_regs1 = n_regs[0];
			n_regs2 = n_regs[1];
			regs1 = regs[0];
			regs2 = regs[1];
		} else {
			if (len1 > 0) regs1 = mm_map(mi, len1, seq1, &n_regs1, tbuf, &mopt, NULL);
			if (paired != 0 && len2 > 0) regs2 = mm_map(mi, len2, seq2, &n_regs2, tbuf, &mopt, NULL);
		}

		n_regs1 = mm2_keep_top_scoring(regs1, n_regs1);
		n_regs2 = mm2_keep_top_scoring(regs2, n_regs2);

		/* Tronko stores paired ASVs in one common orientation and placement scores
		 * both directly against reference-forward tree sequences. Orient the pair
		 * from the strongest minimap2 hit before scoring; orienting mates separately
		 * incorrectly puts a stored pair onto opposite strands. */
		{
			int best_score = -1;
			int reverse_pair = 0;
			for (j = 0; j < n_regs1; j++) {
				if (regs1[j].score > best_score) {
					best_score = regs1[j].score;
					reverse_pair = regs1[j].rev;
				}
			}
			for (j = 0; j < n_regs2; j++) {
				if (regs2[j].score > best_score) {
					best_score = regs2[j].score;
					reverse_pair = regs2[j].rev;
				}
			}
			if (reverse_pair) {
				if (seq1) reverse_complement_in_place(seq1);
				if (seq2) reverse_complement_in_place(seq2);
			}
		}

		if (paired != 0) {
			char cigar_buf[MAX_CIGAR];
			for (j = 0; j < n_regs1; j++) {
				mm_reg1_t *r1 = &regs1[j];
				const char *leaf_name1;
				struct leafMap *lm1;
				int found_concordant = 0;
				int r2_idx = -1;
				int k;

				if (r1->rid < 0 || (uint32_t)r1->rid >= mi->n_seq) continue;
				leaf_name1 = mi->seq[r1->rid].name;
				lm1 = hashmap_get(&map, leaf_name1);
				if (lm1 == NULL) continue;

				cigar_buf[0] = '\0';
				if (r1->p) cigar_to_string(r1->p->cigar, r1->p->n_cigar, cigar_buf, MAX_CIGAR);

				for (k = 0; k < n_regs2; k++) {
					mm_reg1_t *r2 = &regs2[k];
					const char *leaf_name2;
					struct leafMap *lm2;

					if (r2->rid < 0 || (uint32_t)r2->rid >= mi->n_seq) continue;
					leaf_name2 = mi->seq[r2->rid].name;
					lm2 = hashmap_get(&map, leaf_name2);
					if (lm2 == NULL) continue;
					if (lm2->root == lm1->root) {
						found_concordant = 1;
						r2_idx = k;
						break;
					}
				}

				if (found_concordant || concordant == 0) {
					int is_conc = found_concordant;
					add_match(&bwa_results[i], lm1->root, lm1->node, cigar_buf, r1->rs + 1,
					          is_conc, 1, max_leaf_matches);

					if (is_conc && r2_idx >= 0 && bwa_results[i].use_portion == 1) {
						mm_reg1_t *r2 = &regs2[r2_idx];
						char cigar_buf2[MAX_CIGAR];
						int slot;

						cigar_buf2[0] = '\0';
						if (r2->p) cigar_to_string(r2->p->cigar, r2->p->n_cigar, cigar_buf2, MAX_CIGAR);
						for (slot = 0; slot < max_leaf_matches; slot++) {
							if (bwa_results[i].concordant_matches_roots[slot] == lm1->root &&
							    bwa_results[i].concordant_matches_nodes[slot] == lm1->node) {
								strncpy(bwa_results[i].cigars_reverse[slot], cigar_buf2, MAX_CIGAR - 1);
								bwa_results[i].cigars_reverse[slot][MAX_CIGAR - 1] = '\0';
								bwa_results[i].starts_reverse[slot] = r2->rs + 1;
								break;
							}
						}
					}
				}

				if (!found_concordant) {
					add_match(&bwa_results[i], lm1->root, lm1->node, cigar_buf, r1->rs + 1,
					          0, 1, max_leaf_matches);
				}
			}

			for (j = 0; j < n_regs2; j++) {
				mm_reg1_t *r2 = &regs2[j];
				const char *leaf_name2;
				struct leafMap *lm2;
				char cigar_buf[MAX_CIGAR];
				int already = 0;
				int k;

				if (r2->rid < 0 || (uint32_t)r2->rid >= mi->n_seq) continue;
				leaf_name2 = mi->seq[r2->rid].name;
				lm2 = hashmap_get(&map, leaf_name2);
				if (lm2 == NULL) continue;

				for (k = 0; k < max_leaf_matches; k++) {
					if (bwa_results[i].concordant_matches_roots[k] == -1) break;
					if (bwa_results[i].concordant_matches_roots[k] == lm2->root &&
					    bwa_results[i].concordant_matches_nodes[k] == lm2->node) {
						already = 1;
						break;
					}
				}
				if (already) continue;

				cigar_buf[0] = '\0';
				if (r2->p) cigar_to_string(r2->p->cigar, r2->p->n_cigar, cigar_buf, MAX_CIGAR);
				add_match(&bwa_results[i], lm2->root, lm2->node, cigar_buf, r2->rs + 1,
				          0, 0, max_leaf_matches);
			}
		} else {
			char cigar_buf[MAX_CIGAR];
			for (j = 0; j < n_regs1; j++) {
				mm_reg1_t *r1 = &regs1[j];
				const char *leaf_name;
				struct leafMap *lm;

				if (r1->rid < 0 || (uint32_t)r1->rid >= mi->n_seq) continue;
				leaf_name = mi->seq[r1->rid].name;
				lm = hashmap_get(&map, leaf_name);
				if (lm == NULL) continue;

				cigar_buf[0] = '\0';
				if (r1->p) cigar_to_string(r1->p->cigar, r1->p->n_cigar, cigar_buf, MAX_CIGAR);
				add_match(&bwa_results[i], lm->root, lm->node, cigar_buf, r1->rs + 1,
				          0, 1, max_leaf_matches);
			}
		}

		if (regs1) {
			for (j = 0; j < n_regs1; j++) {
				if (regs1[j].p) free(regs1[j].p);
			}
			free(regs1);
		}
		if (regs2) {
			for (j = 0; j < n_regs2; j++) {
				if (regs2[j].p) free(regs2[j].p);
			}
			free(regs2);
		}
	}

	mm_tbuf_destroy(tbuf);

	struct leafMap *blob;
	hashmap_foreach_data(blob, &map) {
		free(blob);
	}
	hashmap_cleanup(&map);
}
