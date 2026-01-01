#include "format_text.h"
#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define BUFFER_SIZE 4096
#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))

// Allocate taxonomy array for a tree
static char ***alloc_taxonomy(int numspec, int max_tax_name) {
    char ***tax = calloc(numspec, sizeof(char **));
    if (!tax) return NULL;

    for (int s = 0; s < numspec; s++) {
        tax[s] = calloc(TAXONOMY_LEVELS, sizeof(char *));
        if (!tax[s]) return NULL;  // Simplified error handling
        for (int l = 0; l < TAXONOMY_LEVELS; l++) {
            tax[s][l] = calloc(max_tax_name + 1, sizeof(char));
            if (!tax[s][l]) return NULL;
        }
    }
    return tax;
}

// Allocate nodes array for a tree
static tronko_node_t *alloc_nodes(int num_nodes, int numbase, int max_nodename) {
    (void)max_nodename;  // Used for name allocation separately
    tronko_node_t *nodes = calloc(num_nodes, sizeof(tronko_node_t));
    if (!nodes) return NULL;

    for (int n = 0; n < num_nodes; n++) {
        nodes[n].up[0] = -2;  // Uninitialized sentinel
        nodes[n].up[1] = -2;
        nodes[n].down = -2;
        nodes[n].depth = -2;
        nodes[n].posteriors = calloc(numbase * 4, sizeof(float));
        if (!nodes[n].posteriors) return NULL;
    }
    return nodes;
}

tronko_db_t *load_text(const char *filename, int verbose) {
    gzFile fp = err_xzopen(filename, "r");
    char buffer[BUFFER_SIZE];
    char *s;

    // Read header: 4 lines
    int numberOfTrees, max_nodename, max_tax_name, max_lineTaxonomy;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &numberOfTrees) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_nodename) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_tax_name) != 1) goto error;

    if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error;
    s = strtok(buffer, "\n");
    if (sscanf(s, "%d", &max_lineTaxonomy) != 1) goto error;

    if (verbose) {
        fprintf(stderr, "  Header: %d trees, max_nodename=%d, max_tax=%d\n",
                numberOfTrees, max_nodename, max_tax_name);
    }

    // Allocate database
    tronko_db_t *db = alloc_db(numberOfTrees);
    if (!db) goto error;

    db->max_nodename = max_nodename;
    db->max_tax_name = max_tax_name;
    db->max_line_taxonomy = max_lineTaxonomy;

    // Read per-tree metadata
    for (int i = 0; i < numberOfTrees; i++) {
        if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;
        s = strtok(buffer, "\n");
        int numbase, root, numspec;
        if (sscanf(s, "%d\t%d\t%d", &numbase, &root, &numspec) != 3) goto error_db;

        db->trees[i].numbase = numbase;
        db->trees[i].root = root;
        db->trees[i].numspec = numspec;
        db->trees[i].num_nodes = 2 * numspec - 1;

        if (verbose) {
            fprintf(stderr, "  Tree %d: numbase=%d, root=%d, numspec=%d\n",
                    i, numbase, root, numspec);
        }
    }

    // Allocate and read taxonomy
    for (int i = 0; i < numberOfTrees; i++) {
        db->trees[i].taxonomy = alloc_taxonomy(db->trees[i].numspec, max_tax_name);
        if (!db->trees[i].taxonomy) goto error_db;

        for (int j = 0; j < db->trees[i].numspec; j++) {
            if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;
            s = strtok(buffer, ";\n");
            if (s) strcpy(db->trees[i].taxonomy[j][0], s);

            for (int k = 1; k < TAXONOMY_LEVELS; k++) {
                s = strtok(NULL, ";\n");
                if (s) strcpy(db->trees[i].taxonomy[j][k], s);
            }
        }
    }

    // Allocate nodes
    for (int i = 0; i < numberOfTrees; i++) {
        db->trees[i].nodes = alloc_nodes(db->trees[i].num_nodes,
                                          db->trees[i].numbase,
                                          max_nodename);
        if (!db->trees[i].nodes) goto error_db;

        // Allocate names only for leaf nodes (indices numspec-1 to 2*numspec-2)
        int numspec = db->trees[i].numspec;
        for (int n = numspec - 1; n < 2 * numspec - 1; n++) {
            db->trees[i].nodes[n].name = calloc(max_nodename + 1, sizeof(char));
            if (!db->trees[i].nodes[n].name) goto error_db;
        }
    }

    // Read node structures and posteriors
    // Total nodes across all trees
    int total_nodes = 0;
    for (int i = 0; i < numberOfTrees; i++) {
        total_nodes += db->trees[i].num_nodes;
    }

    if (verbose) {
        fprintf(stderr, "  Reading %d total nodes...\n", total_nodes);
    }

    int nodes_read = 0;
    while (gzgets(fp, buffer, BUFFER_SIZE)) {
        int treeNumber, nodeNumber;
        int up0, up1, down, depth, taxIdx0, taxIdx1;
        char acc_name[BUFFER_SIZE];
        acc_name[0] = '\0';

        // Parse node header line
        s = strtok(buffer, "\n");
        int parsed = sscanf(s, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s",
                           &treeNumber, &nodeNumber,
                           &up0, &up1, &down, &depth,
                           &taxIdx0, &taxIdx1, acc_name);

        if (parsed < 8) {
            // Might be end of file or parse error
            if (nodes_read > 0) break;
            goto error_db;
        }

        if (treeNumber < 0 || treeNumber >= numberOfTrees) goto error_db;
        if (nodeNumber < 0 || nodeNumber >= db->trees[treeNumber].num_nodes) goto error_db;

        tronko_node_t *node = &db->trees[treeNumber].nodes[nodeNumber];
        node->up[0] = up0;
        node->up[1] = up1;
        node->down = down;
        node->depth = depth;
        node->taxIndex[0] = taxIdx0;
        node->taxIndex[1] = taxIdx1;

        // Copy name for leaf nodes (when both children are -1)
        if (up0 == -1 && up1 == -1 && acc_name[0] != '\0') {
            if (node->name) {
                strcpy(node->name, acc_name);
            }
        }

        // Read posterior probabilities for this node
        int numbase = db->trees[treeNumber].numbase;
        for (int pos = 0; pos < numbase; pos++) {
            if (!gzgets(fp, buffer, BUFFER_SIZE)) goto error_db;

            double p0, p1, p2, p3;
            if (sscanf(buffer, "%lf\t%lf\t%lf\t%lf", &p0, &p1, &p2, &p3) != 4) {
                goto error_db;
            }

            // Store as float (conversion happens here)
            node->posteriors[PP_IDX(pos, 0)] = (float)p0;
            node->posteriors[PP_IDX(pos, 1)] = (float)p1;
            node->posteriors[PP_IDX(pos, 2)] = (float)p2;
            node->posteriors[PP_IDX(pos, 3)] = (float)p3;
        }

        nodes_read++;

        if (verbose && nodes_read % 1000 == 0) {
            fprintf(stderr, "\r  Read %d/%d nodes...", nodes_read, total_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\r  Read %d nodes total    \n", nodes_read);
    }

    err_gzclose(fp);
    return db;

error_db:
    free_db(db);
error:
    gzclose(fp);
    return NULL;
}

int write_text(tronko_db_t *db, const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "w");

    // Write header
    fprintf(fp, "%d\n", db->num_trees);
    fprintf(fp, "%d\n", db->max_nodename);
    fprintf(fp, "%d\n", db->max_tax_name);
    fprintf(fp, "%d\n", db->max_line_taxonomy);

    // Write per-tree metadata
    for (int i = 0; i < db->num_trees; i++) {
        fprintf(fp, "%d\t%d\t%d\n",
                db->trees[i].numbase, db->trees[i].root, db->trees[i].numspec);
    }

    // Write taxonomy
    for (int i = 0; i < db->num_trees; i++) {
        for (int j = 0; j < db->trees[i].numspec; j++) {
            for (int k = 0; k < TAXONOMY_LEVELS; k++) {
                if (k == TAXONOMY_LEVELS - 1) {
                    fprintf(fp, "%s\n", db->trees[i].taxonomy[j][k]);
                } else {
                    fprintf(fp, "%s;", db->trees[i].taxonomy[j][k]);
                }
            }
        }
    }

    // Write nodes and posteriors
    int nodes_written = 0;
    for (int i = 0; i < db->num_trees; i++) {
        for (int j = 0; j < db->trees[i].num_nodes; j++) {
            tronko_node_t *node = &db->trees[i].nodes[j];

            // Node header line
            fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
                    i, j, node->up[0], node->up[1],
                    node->down, node->depth,
                    node->taxIndex[0], node->taxIndex[1]);

            // Name (only for leaf nodes)
            if (node->up[0] == -1 && node->up[1] == -1 && node->name) {
                fprintf(fp, "%s\n", node->name);
            } else {
                fprintf(fp, "\n");
            }

            // Posteriors (match %.17g format from printtree.c)
            int numbase = db->trees[i].numbase;
            for (int pos = 0; pos < numbase; pos++) {
                fprintf(fp, "%.17g\t%.17g\t%.17g\t%.17g\n",
                        (double)node->posteriors[PP_IDX(pos, 0)],
                        (double)node->posteriors[PP_IDX(pos, 1)],
                        (double)node->posteriors[PP_IDX(pos, 2)],
                        (double)node->posteriors[PP_IDX(pos, 3)]);
            }

            nodes_written++;
            if (verbose && nodes_written % 1000 == 0) {
                fprintf(stderr, "\r  Wrote %d nodes...", nodes_written);
            }
        }
    }

    if (verbose) {
        fprintf(stderr, "\r  Wrote %d nodes total    \n", nodes_written);
    }

    err_fclose(fp);
    return 0;
}
