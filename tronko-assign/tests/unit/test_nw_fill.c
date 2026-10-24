/*
 test_nw_fill.c: the two-pass Needleman-Wunsch fill in alignment.c computes the same three score
 matrices as the original single-pass fill, cell for cell.

 tronko-assign aligns every candidate leaf to its read with needleman_wunsch_align(). The
 matrices it fills decide the traceback, and with it the alignment strings, the positions and
 the mismatch counts that placement scores. alignment.c keeps the original fill,
 alignment_fill_matrices(), and adds alignment_fill_matrices_fast(), which aligner_align() uses
 for tronko's scoring. Both are static, so this file includes alignment.c itself.

 For every alignment the test runs aligner_align() (the fast fill), copies the three matrices,
 runs the original fill on the same aligner and compares all (len_a+1)*(len_b+1) cells of each
 matrix with memcmp. Inputs:
  1. alignments recorded from the repository fixtures: every fifth needleman_wunsch_align() call
     of the three production-parity cases (single, unpaired reverse, paired) on the
     Charadriiformes example reference (data/nw_fixture_alignments.txt.gz, leaf TAB read;
     tools/record_nw_alignments.sh records it again);
  2. random alignments from a fixed seed, lengths 0 to 300 (up to 1,174 for leaves): ACGT,
     both cases with N and '-', reads derived from the leaf by substitutions and indels, all
     ASCII bytes 1 to 127, and reads with 129 to 255 distinct byte values (bytes 1 to 255,
     so 0x80 to 0xFC too); every length pair 0..8 x 0..8; nine scorings: tronko's
     (match 2, mismatch -1, gap open -3, extend -1, no start or end gap penalty,
     case-insensitive), both end-gap settings, both start-gap settings, case-sensitive,
     a wildcard and a swap table, another penalty scale, and two that must take the
     original fill (no_mismatches, no_gaps_in_a), as must Smith-Waterman.
 It also checks that the fast fill is the one aligner_align() ran (the profile was built), and
 hashes the alignment strings and score of needleman_wunsch_align2() for the fixture and the
 ASCII random alignments; the hashes must equal the ones the original seq-align code produces
 (alignment.c, alignment_scoring.c and needleman_wunsch.c as at 71f6ec3, unchanged up to the
 parent of the two-pass fill; tests/README.md shows how to recompute them with
 -DNW_TEST_REFERENCE).

 Usage: test_nw_fill <nw_fixture_alignments.txt.gz> [random alignments, default 200000]
 Exit status 0 when nothing differs.
*/
#include "alignment.c"
#include "needleman_wunsch.h"

#include <stdint.h>
#include <zlib.h>

// needleman_wunsch_align2() of the original seq-align code (71f6ec3) on this file's inputs (see
// above): the fixture alignments, and the random ones after the first 20,000 and after 200,000
#define EXPECTED_FIXTURE_HASH 0xcb40bf124bff058dULL
static const long checkpoint_at[2] = { 20000, 200000 };
static const uint64_t checkpoint_expected[2] = { 0x7a4f15eddcf82f66ULL, 0x2517715214af8866ULL };

static uint64_t rng_state = 0x6e772d74776f2d70ULL; // fixed seed

static uint64_t rng_next(void)
{
  uint64_t z = (rng_state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static size_t rng_below(size_t n) { return n ? (size_t)(rng_next() % n) : 0; }

static uint64_t fnv(uint64_t h, const void *p, size_t n)
{
  const unsigned char *c = p;
  for(size_t i = 0; i < n; i++) h = (h ^ c[i]) * 0x100000001b3ULL;
  return h;
}

enum { SC_TRONKO, SC_END_PENALTY, SC_START_END_PENALTY, SC_START_PENALTY, SC_CASE,
       SC_TABLES, SC_SCALE, SC_FAST_COUNT, // the fast fill must run for these
       SC_NO_MISMATCHES = SC_FAST_COUNT, SC_NO_GAPS_A, SC_COUNT };

/* Zeroed, with 8 KiB of zeroed memory before it. Where char is signed (x86-64), scoring_lookup()
   indexes its bitsets and tables with negative values for bytes above 0x7F, reading up to 4 KiB
   before the struct and within it: undefined behaviour of the original code (never reached by
   DNA reads), identical for both fills, which call it with the same arguments. The padding keeps
   those reads inside this allocation and their values fixed. */
#define SCORING_PAD 8192
static scoring_t *make_scoring(int v)
{
  char *block = calloc(1, SCORING_PAD + sizeof(scoring_t));
  scoring_t *s = (scoring_t *)(block + SCORING_PAD);
  bool no_start = !(v == SC_START_END_PENALTY || v == SC_START_PENALTY);
  bool no_end = !(v == SC_END_PENALTY || v == SC_START_END_PENALTY);
  if(v == SC_SCALE)
    scoring_init(s, 5, -4, -10, -2, no_start, no_end, false, false, false, false);
  else
    scoring_init(s, 2, -1, -3, -1, no_start, no_end,
                 v == SC_NO_GAPS_A, false, v == SC_NO_MISMATCHES, v == SC_CASE);
  if(v == SC_TABLES)
  {
    scoring_add_wildcard(s, 'n', -2);
    scoring_add_mutation(s, 'a', 'g', 1);
    scoring_add_mutation(s, 'g', 'a', 1);
    scoring_add_mutation(s, 'c', 't', 0);
  }
  return s;
}

static scoring_t *scorings[SC_COUNT];
static score_t *copy_m, *copy_a, *copy_b;
static size_t copy_cap;
static long fills, fills_differ, fast_not_used, fallback_misused, profile_rows_max;

// Fill with aligner_align() (the fast fill for the fast scorings), then with the original fill
// on the same aligner, and compare every cell of the three matrices.
static void compare_fill(aligner_t *al, const char *a, size_t la, const char *b, size_t lb,
                         int v, char is_sw, const char *what)
{
#ifdef NW_TEST_REFERENCE
  (void)al; (void)a; (void)la; (void)b; (void)lb; (void)v; (void)is_sw; (void)what;
#else
  const scoring_t *sc = scorings[v];
  const int expect_fast = !is_sw && v < SC_FAST_COUNT;
  al->prof_rows = -1;
  aligner_align(al, a, b, la, lb, sc, is_sw);
  if(expect_fast && la > 0 && lb > 0 && al->prof_rows <= 0) fast_not_used++;
  if(!expect_fast && al->prof_rows != -1) fallback_misused++;
  if(al->prof_rows > profile_rows_max) profile_rows_max = al->prof_rows;

  const size_t n = al->score_width * al->score_height, bytes = n * sizeof(score_t);
  if(n > copy_cap)
  {
    copy_cap = n;
    copy_m = realloc(copy_m, bytes);
    copy_a = realloc(copy_a, bytes);
    copy_b = realloc(copy_b, bytes);
  }
  memcpy(copy_m, al->match_scores, bytes);
  memcpy(copy_a, al->gap_a_scores, bytes);
  memcpy(copy_b, al->gap_b_scores, bytes);
  alignment_fill_matrices(al, is_sw);
  fills++;
  if(memcmp(copy_m, al->match_scores, bytes) || memcmp(copy_a, al->gap_a_scores, bytes) ||
     memcmp(copy_b, al->gap_b_scores, bytes))
  {
    if(fills_differ++ < 10)
    {
      size_t k;
      for(k = 0; k < n; k++)
        if(copy_m[k] != al->match_scores[k] || copy_a[k] != al->gap_a_scores[k] ||
           copy_b[k] != al->gap_b_scores[k]) break;
      fprintf(stderr, "DIFFER %s scoring %d is_sw %d: len_a %zu len_b %zu, first cell i=%zu j=%zu:"
              " fast %d %d %d, original %d %d %d\n", what, v, is_sw, la, lb,
              k % al->score_width, k / al->score_width, copy_m[k], copy_a[k], copy_b[k],
              al->match_scores[k], al->gap_a_scores[k], al->gap_b_scores[k]);
    }
  }
#endif
}

// needleman_wunsch_align2() as tronko calls it, hashed: score, then both alignment strings
static uint64_t hash_alignment(uint64_t h, nw_aligner_t *nw, alignment_t *aln,
                               const char *a, size_t la, const char *b, size_t lb, int v)
{
  needleman_wunsch_align2(a, b, la, lb, scorings[v], nw, aln);
  unsigned char s[4];
  uint32_t u = (uint32_t)aln->score;
  for(int k = 0; k < 4; k++) s[k] = (unsigned char)(u >> (8*k));
  h = fnv(h, s, 4);
  h = fnv(h, aln->result_a, strlen(aln->result_a) + 1);
  return fnv(h, aln->result_b, strlen(aln->result_b) + 1);
}

static const char dna[] = "ACGT", mixed[] = "ACGTacgtNn-";

static void random_seq(char *s, size_t n, const char *alphabet, size_t k)
{
  for(size_t i = 0; i < n; i++) s[i] = alphabet[rng_below(k)];
  s[n] = '\0';
}

static void random_bytes(char *s, size_t n, int lo, int hi)
{
  for(size_t i = 0; i < n; i++) s[i] = (char)(lo + (int)rng_below((size_t)(hi - lo + 1)));
  s[n] = '\0';
}

// a read from the leaf: a window with substitutions, insertions, deletions, case and N changes
static size_t derived_read(char *r, size_t cap, const char *leaf, size_t la)
{
  size_t start = rng_below(la), len = 1 + rng_below(la - start), n = 0;
  for(size_t i = start; i < start + len && n < cap; i++)
  {
    size_t x = rng_below(1000);
    if(x < 15) continue;                                   // deletion
    char c = leaf[i];
    if(x < 60) c = dna[rng_below(4)];                      // substitution
    else if(x < 70) c = 'N';
    else if(x < 90) c = (char)tolower((unsigned char)c);
    r[n++] = c;
    if(x >= 990 && n < cap) r[n++] = dna[rng_below(4)];    // insertion
  }
  r[n] = '\0';
  return n;
}

static size_t random_length(void)
{
  size_t x = rng_below(10);
  return x < 6 ? rng_below(41) : x < 9 ? rng_below(151) : rng_below(301);
}

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <nw_fixture_alignments.txt.gz> [random alignments]\n", argv[0]);
    return 2;
  }
  const long n_random = argc > 2 ? atol(argv[2]) : 200000;
  for(int v = 0; v < SC_COUNT; v++) scorings[v] = make_scoring(v);

  nw_aligner_t *nw = needleman_wunsch_new();
  alignment_t *aln = alignment_create(256);
  aligner_t *al = needleman_wunsch_new();
  const size_t cap = 1400;
  char *a = malloc(cap + 1), *b = malloc(cap + 1);
  uint64_t h_fix = 0xcbf29ce484222325ULL, h_rand = 0xcbf29ce484222325ULL, h_check[2] = {0, 0};

  // 1. recorded fixture alignments, tronko's scoring
  gzFile gz = gzopen(argv[1], "rb");
  if(!gz) { fprintf(stderr, "cannot open %s\n", argv[1]); return 2; }
  static char line[8192];
  long n_fixture = 0;
  while(gzgets(gz, line, sizeof(line)))
  {
    char *tab = strchr(line, '\t'), *end = line + strcspn(line, "\r\n");
    if(!tab || tab > end) { fprintf(stderr, "bad fixture line %ld\n", n_fixture + 1); return 2; }
    *tab = '\0'; *end = '\0';
    const char *leaf = line, *read = tab + 1;
    compare_fill(al, leaf, strlen(leaf), read, strlen(read), SC_TRONKO, 0, "fixture");
    h_fix = hash_alignment(h_fix, nw, aln, leaf, strlen(leaf), read, strlen(read), SC_TRONKO);
    n_fixture++;
  }
  gzclose(gz);
  if(n_fixture == 0) { fprintf(stderr, "no fixture alignments in %s\n", argv[1]); return 2; }

  // 2a. every length pair 0..8 x 0..8, every scoring, both fill kinds
  long n_small = 0;
  for(size_t la = 0; la <= 8; la++)
    for(size_t lb = 0; lb <= 8; lb++)
      for(int v = 0; v < SC_COUNT; v++)
        for(int rep = 0; rep < 4; rep++)
        {
          random_seq(a, la, rep < 2 ? dna : mixed, rep < 2 ? 4 : sizeof(mixed) - 1);
          random_seq(b, lb, rep < 2 ? dna : mixed, rep < 2 ? 4 : sizeof(mixed) - 1);
          compare_fill(al, a, la, b, lb, v, 0, "small");
          compare_fill(al, a, la, b, lb, v, 1, "small-sw");
          n_small++;
        }

  // 2b. random alignments
  long n_high = 0;
  for(long t = 0; t < n_random; t++)
  {
    int kind = (int)rng_below(100), v = (int)rng_below(10);
    v = v < 4 ? SC_TRONKO : v - 4 + 1; // 40 % tronko's scoring, 10 % each of six others
    if(v >= SC_COUNT) v = SC_TRONKO;
    size_t la = random_length(), lb = random_length();
    int ascii = 1;
    if(kind < 30) { random_seq(a, la, dna, 4); random_seq(b, lb, dna, 4); }
    else if(kind < 50) { random_seq(a, la, mixed, sizeof(mixed) - 1); random_seq(b, lb, mixed, sizeof(mixed) - 1); }
    else if(kind < 75)
    {
      la = kind == 74 ? 100 + rng_below(1075) : 1 + random_length(); // up to FWH's 1,174 columns
      random_seq(a, la, dna, 4);
      lb = derived_read(b, 300, a, la);
    }
    else if(kind < 90) { random_bytes(a, la, 1, 127); random_bytes(b, lb, 1, 127); }
    else
    {
      // a read with k distinct byte values, 129 <= k <= 255, from 1..255
      unsigned char vals[255];
      size_t k = 129 + rng_below(127), i;
      for(i = 0; i < 255; i++) vals[i] = (unsigned char)(i + 1);
      for(i = 0; i < k; i++) { size_t j = i + rng_below(255 - i); unsigned char x = vals[i]; vals[i] = vals[j]; vals[j] = x; }
      lb = k + rng_below(200);
      for(i = 0; i < lb; i++) b[i] = (char)(i < k ? vals[i] : vals[rng_below(k)]);
      for(i = 0; i < lb; i++) { size_t j = i + rng_below(lb - i); char x = b[i]; b[i] = b[j]; b[j] = x; }
      b[lb] = '\0';
      la = 1 + random_length();
      if(kind < 95) random_seq(a, la, dna, 4); else random_bytes(a, la, 1, 255);
      ascii = 0;
      n_high++;
    }
    compare_fill(al, a, la, b, lb, v, 0, ascii ? "random" : "random-high-bytes");
    // the hash is platform independent only for ASCII (char is signed on x86-64, not on aarch64)
    if(ascii && la > 0 && lb > 0 && v < SC_FAST_COUNT)
      h_rand = hash_alignment(h_rand, nw, aln, a, la, b, lb, v);
    for(int c = 0; c < 2; c++) if(t + 1 == checkpoint_at[c]) h_check[c] = h_rand;
  }

#ifdef NW_TEST_REFERENCE
  printf("test_nw_fill reference: fixture 0x%016llxULL, random after %ld 0x%016llxULL, after %ld"
         " 0x%016llxULL (%ld fixture, %ld random alignments)\n", (unsigned long long)h_fix,
         checkpoint_at[0], (unsigned long long)h_check[0], checkpoint_at[1],
         (unsigned long long)h_check[1], n_fixture, n_random);
  return 0;
#else
  int fail = fills_differ || fast_not_used || fallback_misused;
  int fix_ok = h_fix == EXPECTED_FIXTURE_HASH, rand_ok = 1, rand_checked = 0;
  for(int c = 0; c < 2; c++)
    if(n_random >= checkpoint_at[c])
    {
      rand_checked = c + 1;
      if(h_check[c] != checkpoint_expected[c])
      {
        rand_ok = 0;
        fprintf(stderr, "end-to-end hash after %ld random alignments: 0x%016llx, original 0x%016llx\n",
                checkpoint_at[c], (unsigned long long)h_check[c],
                (unsigned long long)checkpoint_expected[c]);
      }
    }
  printf("test_nw_fill: %ld fixture, %ld small, %ld random alignments (%ld with 129+ distinct"
         " bytes); %ld fills compared, %ld differ; fast fill not used %ld, fallback misused %ld;"
         " end-to-end hash vs the original code: fixture %s, random %s\n",
         n_fixture, n_small, n_random, n_high, fills, fills_differ, fast_not_used,
         fallback_misused, fix_ok ? "equal" : "DIFFERENT",
         !rand_checked ? "not checked (fewer than 20000)" : !rand_ok ? "DIFFERENT"
         : rand_checked == 2 ? "equal at 20000 and 200000" : "equal at 20000");
  if(!fix_ok)
    fprintf(stderr, "end-to-end hash of the fixture alignments: 0x%016llx, original 0x%016llx\n",
            (unsigned long long)h_fix, (unsigned long long)EXPECTED_FIXTURE_HASH);
  needleman_wunsch_free(nw);
  needleman_wunsch_free(al);
  alignment_free(aln);
  for(int v = 0; v < SC_COUNT; v++) free((char *)scorings[v] - SCORING_PAD);
  free(a); free(b); free(copy_m); free(copy_a); free(copy_b);
  return fail || !fix_ok || !rand_ok;
#endif
}
