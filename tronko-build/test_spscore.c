/*
 * test_spscore.c — golden test for calculateSPArr
 *
 * Constructs synthetic masterArr objects with known MSA data and
 * verifies calculateSPArr produces the expected output.
 *
 * Compile:
 *   gcc -O3 -o test_spscore test_spscore.c -lm
 *
 * Run:
 *   ./test_spscore
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Minimal definitions to compile calculateSPArr standalone */
#define STATESPACE 20
#define MAX_NODE_CHILDREN 1000
#define MAX_NODENAME 30
#define type_of_PP double

typedef struct node {
    int up[MAX_NODE_CHILDREN];
    int down;
    int nd;
    int depth;
    double bl;
    double *like;
    double **likenc;
    double *posterior;
    double **posteriornc;
    int s;
    int numsites;
    int spec;
    int mrca;
    char *name;
    int taxIndex[2];
    type_of_PP score;
    int SP_fail;
} node;

typedef struct masterArr {
    char index[10];
    node *tree;
    int **msa;
    char ***taxonomy;
    int numspec;
    int numNodes;
    int treeCapacity;
    int root;
    int numbase;
    char **names;
    char *filename;
} masterArr;

/* ---- calculateSPArr (optimized version under test) ---- */
double calculateSPArr(struct masterArr *m){
    int i,j,k;
    int *partition = (int*)malloc(sizeof(int)*m->numspec);
    for(i=0; i<m->numspec; i++){
        partition[i] = m->numspec-1+i;
    }
    long long raw_score = 0;
    long long numpairs = 0;
    for( j=0; j<m->numspec; j++){
        int *row_j = m->msa[partition[j]-m->numspec+1];
        for( k=j+1; k<m->numspec; k++){
            int *row_k = m->msa[partition[k]-m->numspec+1];
            long long pair_score = 0;
            for (i=0; i<m->numbase; i++){
                int a = row_j[i];
                int b = row_k[i];
                int eq  = (a == b);
                int gap = (a == 4) | (b == 4);
                pair_score += eq * (5 - gap) - 2 + gap;
            }
            raw_score += pair_score;
            numpairs += m->numbase;
        }
    }
    free(partition);
    return (double)raw_score / (double)numpairs;
}
/* ---------------------------------------------------------- */

/*
 * Build a minimal masterArr for testing.
 * seqs: array of numspec sequences, each of length numbase (encoded as ints:
 *       0=A, 1=C, 2=G, 3=T, 4=gap).
 * The tree array needs leaf nodes at indices [numspec-1 .. 2*numspec-2].
 */
static struct masterArr *make_master(int numspec, int numbase, int seqs[][32]) {
    struct masterArr *m = calloc(1, sizeof(masterArr));
    m->numspec = numspec;
    m->numbase = numbase;
    /* Leaves live at [numspec-1 .. 2*numspec-2] in tronko convention */
    m->numNodes = 2 * numspec - 1;
    m->tree = calloc(m->numNodes, sizeof(node));
    m->msa  = malloc(numspec * sizeof(int*));
    for (int i = 0; i < numspec; i++) {
        m->msa[i] = malloc(numbase * sizeof(int));
        for (int b = 0; b < numbase; b++)
            m->msa[i][b] = seqs[i][b];
    }
    return m;
}

static void free_master(struct masterArr *m) {
    for (int i = 0; i < m->numspec; i++) free(m->msa[i]);
    free(m->msa);
    free(m->tree);
    free(m);
}

/* Compute expected SPscore using exact long long arithmetic (reference) */
static double expected_spscore(int numspec, int numbase, int seqs[][32]) {
    long long raw = 0;
    long long numpairs = 0;
    for (int j = 0; j < numspec; j++) {
        for (int k = j+1; k < numspec; k++) {
            for (int i = 0; i < numbase; i++) {
                if (seqs[j][i] == seqs[k][i])                      raw += 3;
                else if (seqs[j][i] != 4 && seqs[k][i] != 4)      raw -= 2;
                else                                                raw -= 1;
            }
            numpairs += numbase;
        }
    }
    return (double)raw / (double)numpairs;
}

static int run_test(const char *name, int numspec, int numbase, int seqs[][32]) {
    struct masterArr *m = make_master(numspec, numbase, seqs);
    double got = calculateSPArr(m);
    double want = expected_spscore(numspec, numbase, seqs);
    free_master(m);

    int pass = (fabs(got - want) < 1e-12);
    printf("%-40s SPscore=%.15f  %s\n", name, got, pass ? "PASS" : "FAIL");
    if (!pass)
        printf("  EXPECTED %.15f  DIFF %.3e\n", want, fabs(got - want));
    return pass;
}

int main(void) {
    int all_pass = 1;

    /* Test 1: all identical — score should be +3.0 */
    {
        int seqs[4][32] = {
            {0,1,2,3, 0,1,2,3},
            {0,1,2,3, 0,1,2,3},
            {0,1,2,3, 0,1,2,3},
            {0,1,2,3, 0,1,2,3},
        };
        all_pass &= run_test("all_identical(4seq,8base)", 4, 8, seqs);
    }

    /* Test 2: all mismatched (no gaps) — score should be -2.0 */
    {
        int seqs[3][32] = {
            {0,0,0,0},
            {1,1,1,1},
            {2,2,2,2},
        };
        all_pass &= run_test("all_mismatch_no_gap(3seq,4base)", 3, 4, seqs);
    }

    /* Test 3: all gaps — score should be -1.0 */
    {
        int seqs[3][32] = {
            {4,4,4,4},
            {4,4,4,4},
            {4,4,4,4},
        };
        all_pass &= run_test("all_gaps(3seq,4base)", 3, 4, seqs);
    }

    /* Test 4: mixed matches, mismatches, gaps */
    {
        int seqs[3][32] = {
            {0,1,4,2},
            {0,2,4,3},
            {1,1,0,2},
        };
        all_pass &= run_test("mixed(3seq,4base)", 3, 4, seqs);
    }

    /* Test 5: two sequences, single base — match */
    {
        int seqs[2][32] = {{0},{0}};
        all_pass &= run_test("2seq_1base_match", 2, 1, seqs);
    }

    /* Test 6: two sequences, single base — mismatch */
    {
        int seqs[2][32] = {{0},{1}};
        all_pass &= run_test("2seq_1base_mismatch", 2, 1, seqs);
    }

    /* Test 7: larger case — 10 seqs, 20 bases, random-ish */
    {
        int seqs[10][32] = {
            {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3},
            {0,0,2,3,1,1,2,3,0,1,2,4,0,1,2,3,0,1,2,3},
            {1,1,2,3,0,1,3,3,0,1,2,3,0,2,2,3,0,1,2,3},
            {0,1,2,2,0,1,2,3,4,1,2,3,0,1,2,3,0,1,2,3},
            {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3},
            {2,1,2,3,0,1,2,0,0,1,2,3,0,1,2,3,0,1,2,3},
            {0,1,3,3,0,1,2,3,0,1,4,3,0,1,2,3,0,1,2,3},
            {0,1,2,3,0,2,2,3,0,1,2,3,1,1,2,3,0,1,2,3},
            {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,1,1,2,3},
            {0,1,2,1,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,0},
        };
        all_pass &= run_test("10seq_20base_mixed", 10, 20, seqs);
    }

    /* Benchmark on synthetic data matching typical partition sizes */
    printf("\n=== Benchmark ===\n");
    int bench_sizes[] = {100, 500, 1000, 5000, 0};
    int bench_numbase = 1800;
    for (int s = 0; bench_sizes[s] != 0; s++) {
        int ns = bench_sizes[s];
        /* Build a synthetic masterArr with random sequences */
        struct masterArr *bm = calloc(1, sizeof(masterArr));
        bm->numspec = ns;
        bm->numbase = bench_numbase;
        bm->numNodes = 2*ns-1;
        bm->tree = calloc(bm->numNodes, sizeof(node));
        bm->msa = malloc(ns * sizeof(int*));
        srand(42);
        for (int i = 0; i < ns; i++) {
            bm->msa[i] = malloc(bench_numbase * sizeof(int));
            for (int b = 0; b < bench_numbase; b++)
                bm->msa[i][b] = rand() % 5; /* 0-3 = ACGT, 4 = gap */
        }
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        double score = calculateSPArr(bm);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
        printf("  numspec=%5d  numbase=%d  SPscore=%.6f  time=%.3fs\n",
               ns, bench_numbase, score, elapsed);
        for (int i = 0; i < ns; i++) free(bm->msa[i]);
        free(bm->msa); free(bm->tree); free(bm);
    }

    printf("\n%s\n", all_pass ? "ALL TESTS PASSED" : "SOME TESTS FAILED");
    return all_pass ? 0 : 1;
}
