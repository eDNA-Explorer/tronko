/*
 * minimap2_wrapper.c
 *
 * Provides run_minimap2() with the same contract as run_bwa().
 * Uses minimap2's C library API to map reads and populate bwaMatches structs.
 *
 * Key differences from BWA path:
 * - No SAM text intermediary — results are structured mm_reg1_t arrays
 * - Index built at runtime from reference FASTA (cached across calls)
 * - CIGAR converted from integer array to string for use_portion compatibility
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "minimap2_wrapper.h"
#include "minimap2_src/minimap.h"
#include "hashmap.h"
#include "hashmap_base.h"

extern queryMatPaired *pairedQueryMat;
extern queryMatSingle *singleQueryMat;
extern node **treeArr;
extern int *numspecArr;

/* Cached minimap2 index — built once, reused across calls */
static mm_idx_t *g_mm2_idx = NULL;
static char g_mm2_db_path[MAXFILENAME] = "";
static int g_mm2_kmer = 0;
static int g_mm2_window = 0;

/*
 * Convert minimap2's integer CIGAR array to a CIGAR string.
 * minimap2 encodes each CIGAR op as: (length << 4) | op
 * where op is 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
 */
static void cigar_to_string(const uint32_t *cigar, int n_cigar, char *buf, int buf_size)
{
	int pos = 0;
	int i;
	for (i = 0; i < n_cigar && pos < buf_size - 12; ++i) {
		int op = cigar[i] & 0xf;
		int len = cigar[i] >> 4;
		/* Map =/X to M for compatibility with BWA-style CIGAR parsing */
		char op_char;
		if (op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH)
			op_char = 'M';
		else
			op_char = MM_CIGAR_STR[op];
		pos += snprintf(buf + pos, buf_size - pos, "%d%c", len, op_char);
	}
	if (pos >= buf_size) pos = buf_size - 1;
	buf[pos] = '\0';
}

/*
 * Build or retrieve the cached minimap2 index for the given reference FASTA.
 */
static mm_idx_t *get_or_build_index(const char *db_path, int kmer, int window)
{
	/* Return cached index if parameters match */
	if (g_mm2_idx != NULL &&
	    strcmp(g_mm2_db_path, db_path) == 0 &&
	    g_mm2_kmer == kmer &&
	    g_mm2_window == window) {
		return g_mm2_idx;
	}

	/* Destroy old index if it exists */
	if (g_mm2_idx != NULL) {
		mm_idx_destroy(g_mm2_idx);
		g_mm2_idx = NULL;
	}

	/* Suppress minimap2 stderr messages */
	mm_verbose = 1;

	/* Build index from reference FASTA */
	mm_idxopt_t iopt;
	mm_mapopt_t mopt_dummy;
	mm_set_opt(0, &iopt, &mopt_dummy);
	iopt.k = kmer;
	iopt.w = window;

	mm_idx_reader_t *reader = mm_idx_reader_open(db_path, &iopt, NULL);
	if (reader == NULL) {
		fprintf(stderr, "[minimap2] ERROR: cannot open reference file '%s'\n", db_path);
		return NULL;
	}

	g_mm2_idx = mm_idx_reader_read(reader, 1);
	mm_idx_reader_close(reader);

	if (g_mm2_idx == NULL) {
		fprintf(stderr, "[minimap2] ERROR: failed to build index from '%s'\n", db_path);
		return NULL;
	}

	/* Cache parameters */
	strncpy(g_mm2_db_path, db_path, sizeof(g_mm2_db_path) - 1);
	g_mm2_db_path[sizeof(g_mm2_db_path) - 1] = '\0';
	g_mm2_kmer = kmer;
	g_mm2_window = window;

	return g_mm2_idx;
}

/*
 * Add a match to the concordant or discordant arrays of a bwaMatches entry.
 * Returns the new count (k) after insertion, or the same k if duplicate/full.
 */
static int add_match(bwaMatches *result, int root, int node_id,
                     const char *cigar_str, int start_pos,
                     int is_concordant, int is_forward,
                     int max_bwa_matches)
{
	int *roots, *nodes;
	int k, l;

	if (is_concordant) {
		roots = result->concordant_matches_roots;
		nodes = result->concordant_matches_nodes;
	} else {
		roots = result->discordant_matches_roots;
		nodes = result->discordant_matches_nodes;
	}

	/* Find current count */
	for (k = 0; k < max_bwa_matches; k++) {
		if (roots[k] == -1) break;
	}

	/* Check for duplicate */
	for (l = 0; l < k; l++) {
		if (roots[l] == root && nodes[l] == node_id) {
			return k; /* Already present */
		}
	}

	/* Check capacity */
	if (k >= max_bwa_matches) return k;

	/* Add the match */
	roots[k] = root;
	nodes[k] = node_id;

	/* Store CIGAR and start position for use_portion mode */
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

void run_minimap2(int start, int end, bwaMatches *bwa_results, int concordant,
                  int numberOfTrees, char *databaseFile, int paired,
                  int max_query_length, int max_readname_length,
                  int max_acc_name, int max_bwa_matches,
                  int mm2_kmer, int mm2_window)
{
	int i, j;

	/* Build or retrieve cached index */
	mm_idx_t *mi = get_or_build_index(databaseFile, mm2_kmer, mm2_window);
	if (mi == NULL) {
		/* Mark all reads as unmatched */
		for (i = 0; i < end - start; i++) {
			bwa_results[i].concordant_matches_roots[0] = -1;
			bwa_results[i].concordant_matches_nodes[0] = -1;
			bwa_results[i].discordant_matches_roots[0] = -1;
			bwa_results[i].discordant_matches_nodes[0] = -1;
		}
		return;
	}

	/* Set up mapping parameters */
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	mm_set_opt(0, &iopt, &mopt);   /* Initialize defaults */
	mm_set_opt("sr", &iopt, &mopt); /* Apply short-read preset */
	mopt.flag |= MM_F_CIGAR;       /* Generate CIGAR for use_portion */
	mopt.best_n = max_bwa_matches;  /* Report up to max_bwa_matches hits */
	mm_mapopt_update(&mopt, mi);

	/* Build leaf name → (tree_id, node_id) hashmap (same as fastmap.c:112-123) */
	HASHMAP(char, leafMap) map;
	hashmap_init(&map, hashmap_hash_string, strcmp);
	for (i = 0; i < numberOfTrees; i++) {
		for (j = numspecArr[i] - 1; j < 2 * numspecArr[i] - 1; j++) {
			struct leafMap *lm = malloc(sizeof(struct leafMap));
			lm->name = treeArr[i][j].name;
			lm->root = i;
			lm->node = j;
			hashmap_put(&map, lm->name, lm);
		}
	}

	/* Thread-local buffer for minimap2 */
	mm_tbuf_t *tbuf = mm_tbuf_init();

	/* Map each read */
	int n_reads = end - start;
	for (i = 0; i < n_reads; i++) {
		char *seq1 = NULL, *seq2 = NULL;
		int len1 = 0, len2 = 0;

		if (paired != 0) {
			seq1 = pairedQueryMat->query1Mat[start + i];
			seq2 = pairedQueryMat->query2Mat[start + i];
			if (seq1) len1 = strlen(seq1);
			if (seq2) len2 = strlen(seq2);
		} else {
			seq1 = singleQueryMat->queryMat[start + i];
			if (seq1) len1 = strlen(seq1);
		}

		/* Skip reads with no sequence */
		if (len1 == 0 && len2 == 0) {
			bwa_results[i].concordant_matches_roots[0] = -1;
			bwa_results[i].concordant_matches_nodes[0] = -1;
			bwa_results[i].discordant_matches_roots[0] = -1;
			bwa_results[i].discordant_matches_nodes[0] = -1;
			continue;
		}

		/* Map forward read */
		int n_regs1 = 0;
		mm_reg1_t *regs1 = NULL;
		if (len1 > 0) {
			regs1 = mm_map(mi, len1, seq1, &n_regs1, tbuf, &mopt, NULL);
		}

		/* Map reverse read (paired-end only) */
		int n_regs2 = 0;
		mm_reg1_t *regs2 = NULL;
		if (paired != 0 && len2 > 0) {
			regs2 = mm_map(mi, len2, seq2, &n_regs2, tbuf, &mopt, NULL);
		}

		if (paired != 0) {
			/*
			 * Paired-end logic:
			 * For each R1 hit, check if any R2 hit maps to a leaf in the same tree.
			 * If yes → concordant. If no → discordant.
			 * Mirrors the BWA SAM parsing logic in fastmap.c:196-271.
			 */
			char cigar_buf[MAX_CIGAR];

			for (j = 0; j < n_regs1; j++) {
				mm_reg1_t *r1 = &regs1[j];
				if (r1->rid < 0 || (uint32_t)r1->rid >= mi->n_seq) continue;

				const char *leaf_name1 = mi->seq[r1->rid].name;
				struct leafMap *lm1 = hashmap_get(&map, leaf_name1);
				if (lm1 == NULL) continue;

				/* Convert CIGAR for use_portion */
				cigar_buf[0] = '\0';
				if (r1->p) {
					cigar_to_string(r1->p->cigar, r1->p->n_cigar, cigar_buf, MAX_CIGAR);
				}
				/* minimap2 uses 0-based start; BWA SAM uses 1-based */
				int start_pos1 = r1->rs + 1;

				/* Check if any R2 hit is in the same tree → concordant */
				int found_concordant = 0;
				int r2_idx = -1;
				int k;
				for (k = 0; k < n_regs2; k++) {
					mm_reg1_t *r2 = &regs2[k];
					if (r2->rid < 0 || (uint32_t)r2->rid >= mi->n_seq) continue;
					const char *leaf_name2 = mi->seq[r2->rid].name;
					struct leafMap *lm2 = hashmap_get(&map, leaf_name2);
					if (lm2 == NULL) continue;
					if (lm2->root == lm1->root) {
						found_concordant = 1;
						r2_idx = k;
						break;
					}
				}

				if (found_concordant || concordant == 0) {
					int is_conc = found_concordant;
					add_match(&bwa_results[i], lm1->root, lm1->node,
					          cigar_buf, start_pos1,
					          is_conc, 1, max_bwa_matches);

					/* Add reverse CIGAR if concordant and we have an R2 match */
					if (is_conc && r2_idx >= 0 && bwa_results[i].use_portion == 1) {
						mm_reg1_t *r2 = &regs2[r2_idx];
						char cigar_buf2[MAX_CIGAR];
						cigar_buf2[0] = '\0';
						if (r2->p) {
							cigar_to_string(r2->p->cigar, r2->p->n_cigar, cigar_buf2, MAX_CIGAR);
						}
						/* Find the slot we just added to and store reverse CIGAR */
						int slot;
						for (slot = 0; slot < max_bwa_matches; slot++) {
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
					add_match(&bwa_results[i], lm1->root, lm1->node,
					          cigar_buf, start_pos1,
					          0, 1, max_bwa_matches);
				}
			}

			/* Also add any R2 hits that weren't matched as concordant */
			for (j = 0; j < n_regs2; j++) {
				mm_reg1_t *r2 = &regs2[j];
				if (r2->rid < 0 || (uint32_t)r2->rid >= mi->n_seq) continue;

				const char *leaf_name2 = mi->seq[r2->rid].name;
				struct leafMap *lm2 = hashmap_get(&map, leaf_name2);
				if (lm2 == NULL) continue;

				cigar_buf[0] = '\0';
				if (r2->p) {
					cigar_to_string(r2->p->cigar, r2->p->n_cigar, cigar_buf, MAX_CIGAR);
				}
				int start_pos2 = r2->rs + 1;

				/* Check if already added as concordant */
				int already = 0;
				int k;
				for (k = 0; k < max_bwa_matches; k++) {
					if (bwa_results[i].concordant_matches_roots[k] == -1) break;
					if (bwa_results[i].concordant_matches_roots[k] == lm2->root &&
					    bwa_results[i].concordant_matches_nodes[k] == lm2->node) {
						already = 1;
						break;
					}
				}
				if (!already) {
					add_match(&bwa_results[i], lm2->root, lm2->node,
					          cigar_buf, start_pos2,
					          0, 0, max_bwa_matches);
				}
			}
		} else {
			/*
			 * Single-end: all hits go to discordant matches
			 * (BWA single-end has no mate, so SAM parser classifies as discordant)
			 */
			char cigar_buf[MAX_CIGAR];

			for (j = 0; j < n_regs1; j++) {
				mm_reg1_t *r1 = &regs1[j];
				if (r1->rid < 0 || (uint32_t)r1->rid >= mi->n_seq) continue;

				const char *leaf_name = mi->seq[r1->rid].name;
				struct leafMap *lm = hashmap_get(&map, leaf_name);
				if (lm == NULL) continue;

				cigar_buf[0] = '\0';
				if (r1->p) {
					cigar_to_string(r1->p->cigar, r1->p->n_cigar, cigar_buf, MAX_CIGAR);
				}
				int start_pos = r1->rs + 1;

				add_match(&bwa_results[i], lm->root, lm->node,
				          cigar_buf, start_pos,
				          0, 1, max_bwa_matches);
			}
		}

		/* Free minimap2 results */
		if (regs1) {
			for (j = 0; j < n_regs1; j++)
				if (regs1[j].p) free(regs1[j].p);
			free(regs1);
		}
		if (regs2) {
			for (j = 0; j < n_regs2; j++)
				if (regs2[j].p) free(regs2[j].p);
			free(regs2);
		}
	}

	mm_tbuf_destroy(tbuf);

	/* Cleanup hashmap */
	struct leafMap *blob;
	hashmap_foreach_data(blob, &map) {
		free(blob);
	}
	hashmap_cleanup(&map);
}
