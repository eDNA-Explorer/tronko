#include "format_common.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int detect_format(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return FORMAT_UNKNOWN;

    uint8_t magic[4];
    size_t n = fread(magic, 1, 4, fp);
    fclose(fp);

    if (n < 4) return FORMAT_UNKNOWN;

    // Check for binary format magic
    if (magic[0] == TRONKO_MAGIC_0 && magic[1] == TRONKO_MAGIC_1 &&
        magic[2] == TRONKO_MAGIC_2 && magic[3] == TRONKO_MAGIC_3) {
        return FORMAT_BINARY;
    }

    // Check for gzip magic (text files may be gzipped)
    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        return FORMAT_TEXT;  // Gzipped text
    }

    // Assume text (first bytes should be ASCII digits)
    if (magic[0] >= '0' && magic[0] <= '9') {
        return FORMAT_TEXT;
    }

    return FORMAT_UNKNOWN;
}

tronko_db_t *alloc_db(int num_trees) {
    tronko_db_t *db = calloc(1, sizeof(tronko_db_t));
    if (!db) return NULL;

    db->num_trees = num_trees;
    db->trees = calloc(num_trees, sizeof(tronko_tree_t));
    if (!db->trees) {
        free(db);
        return NULL;
    }

    return db;
}

void free_db(tronko_db_t *db) {
    if (!db) return;

    for (int t = 0; t < db->num_trees; t++) {
        tronko_tree_t *tree = &db->trees[t];

        // Free taxonomy
        if (tree->taxonomy) {
            for (int s = 0; s < tree->numspec; s++) {
                if (tree->taxonomy[s]) {
                    for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                        free(tree->taxonomy[s][l]);
                    }
                    free(tree->taxonomy[s]);
                }
            }
            free(tree->taxonomy);
        }

        // Free nodes
        if (tree->nodes) {
            for (int n = 0; n < tree->num_nodes; n++) {
                free(tree->nodes[n].name);
                free(tree->nodes[n].posteriors);
            }
            free(tree->nodes);
        }
    }

    free(db->trees);
    free(db);
}
