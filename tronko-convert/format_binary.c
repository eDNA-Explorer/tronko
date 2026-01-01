#include "format_binary.h"
#include "utils.h"
#include "crc32.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PP_IDX(pos, nuc) ((pos) * 4 + (nuc))

// Binary file header size constants
#define FILE_HEADER_SIZE 64
#define GLOBAL_META_SIZE 16
#define TREE_META_SIZE 12
#define NODE_RECORD_SIZE 32

int write_binary(tronko_db_t *db, const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "wb");

    if (verbose) {
        fprintf(stderr, "  Writing binary format v%d.%d\n",
                TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR);
    }

    // Calculate section offsets
    uint64_t header_end = FILE_HEADER_SIZE + GLOBAL_META_SIZE;
    uint64_t tree_meta_size = db->num_trees * TREE_META_SIZE;
    uint64_t taxonomy_offset = header_end + tree_meta_size;

    // Calculate taxonomy section size
    uint64_t taxonomy_size = db->num_trees * 8;  // Header per tree (size + reserved)
    for (int t = 0; t < db->num_trees; t++) {
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                // Length prefix (2 bytes) + string + null
                int len = strlen(db->trees[t].taxonomy[s][l]) + 1;
                taxonomy_size += 2 + len;
            }
        }
    }

    uint64_t node_offset = taxonomy_offset + taxonomy_size;

    // Calculate node section size
    uint64_t node_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        node_size += 4;  // Node count header
        node_size += db->trees[t].num_nodes * NODE_RECORD_SIZE;
        // Name table for leaf nodes
        for (int n = db->trees[t].numspec - 1; n < db->trees[t].num_nodes; n++) {
            if (db->trees[t].nodes[n].name) {
                int len = strlen(db->trees[t].nodes[n].name) + 1;
                node_size += 2 + len;  // Length prefix + string + null
            }
        }
    }

    uint64_t posterior_offset = node_offset + node_size;

    // Calculate posterior section size
    uint64_t posterior_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        posterior_size += (uint64_t)db->trees[t].num_nodes * db->trees[t].numbase * 4 * sizeof(float);
    }

    uint64_t total_size = posterior_offset + posterior_size + 8;  // +8 for footer

    if (verbose) {
        fprintf(stderr, "  Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset);
        fprintf(stderr, "  Total file size: %lu bytes (%.2f MB)\n",
                (unsigned long)total_size, total_size / (1024.0 * 1024.0));
    }

    // === Write File Header (64 bytes) ===
    uint8_t header_bytes[8] = {
        TRONKO_MAGIC_0, TRONKO_MAGIC_1, TRONKO_MAGIC_2, TRONKO_MAGIC_3,
        TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR,
        0x01,  // Little-endian
        0x01   // Float precision
    };
    err_fwrite(header_bytes, 1, 8, fp);

    // Calculate header CRC (of bytes 0-7)
    uint32_t header_crc = tronko_crc32(0, header_bytes, 8);
    write_u32(fp, header_crc);

    write_u32(fp, 0);  // Reserved (alignment)
    write_u64(fp, taxonomy_offset);
    write_u64(fp, node_offset);
    write_u64(fp, posterior_offset);
    write_u64(fp, total_size);

    // Padding to 64 bytes
    uint8_t reserved[16] = {0};
    err_fwrite(reserved, 1, 16, fp);

    // === Write Global Metadata (16 bytes) ===
    write_i32(fp, db->num_trees);
    write_i32(fp, db->max_nodename);
    write_i32(fp, db->max_tax_name);
    write_i32(fp, db->max_line_taxonomy);

    // === Write Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < db->num_trees; i++) {
        write_i32(fp, db->trees[i].numbase);
        write_i32(fp, db->trees[i].root);
        write_i32(fp, db->trees[i].numspec);
    }

    // Verify we're at taxonomy_offset
    if ((uint64_t)err_ftell(fp) != taxonomy_offset) {
        fprintf(stderr, "Error: Taxonomy offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)taxonomy_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Taxonomy Section ===
    for (int t = 0; t < db->num_trees; t++) {
        // Calculate and write size for this tree's taxonomy
        uint32_t tree_tax_size = 0;
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                tree_tax_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
        write_u32(fp, tree_tax_size);
        write_u32(fp, 0);  // Reserved

        // Write taxonomy strings
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                const char *str = db->trees[t].taxonomy[s][l];
                uint16_t len = strlen(str) + 1;  // Include null terminator
                write_u16(fp, len);
                err_fwrite(str, 1, len, fp);
            }
        }
    }

    // Verify we're at node_offset
    if ((uint64_t)err_ftell(fp) != node_offset) {
        fprintf(stderr, "Error: Node offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)node_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Node Section ===
    for (int t = 0; t < db->num_trees; t++) {
        write_u32(fp, db->trees[t].num_nodes);

        // Track name offsets (relative to start of name table for this tree)
        uint32_t name_offset = 0;

        // First pass: write node records with calculated name offsets
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            write_i32(fp, node->up[0]);
            write_i32(fp, node->up[1]);
            write_i32(fp, node->down);
            write_i32(fp, node->depth);
            write_i32(fp, node->taxIndex[0]);
            write_i32(fp, node->taxIndex[1]);

            // Name offset (0 for internal nodes)
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                write_u32(fp, name_offset + 1);  // +1 so 0 means "no name"
                name_offset += 2 + strlen(node->name) + 1;
            } else {
                write_u32(fp, 0);
            }
            write_u32(fp, 0);  // Reserved
        }

        // Second pass: write name table
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                uint16_t len = strlen(node->name) + 1;
                write_u16(fp, len);
                err_fwrite(node->name, 1, len, fp);
            }
        }
    }

    // Verify we're at posterior_offset
    if ((uint64_t)err_ftell(fp) != posterior_offset) {
        fprintf(stderr, "Error: Posterior offset mismatch (expected %lu, got %ld)\n",
                (unsigned long)posterior_offset, err_ftell(fp));
        err_fclose(fp);
        return -1;
    }

    // === Write Posterior Section ===
    if (verbose) {
        fprintf(stderr, "  Writing posteriors...\n");
    }

    for (int t = 0; t < db->num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            // Bulk write all posteriors for this node
            err_fwrite(db->trees[t].nodes[n].posteriors,
                      sizeof(float),
                      db->trees[t].numbase * 4,
                      fp);
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: wrote %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // === Write Footer ===
    // CRC of all preceding data - need to close and reopen in read mode
    long data_end = err_ftell(fp);
    err_fclose(fp);

    // Reopen in read mode to calculate CRC
    FILE *rfp = err_xopen(filename, "rb");
    uint32_t data_crc = tronko_crc32_file(rfp, 0, data_end);
    err_fclose(rfp);

    // Reopen in append mode to write footer
    fp = err_xopen(filename, "ab");
    write_u32(fp, data_crc);
    write_u32(fp, TRONKO_FOOTER_MAGIC);

    // Verify total size
    if ((uint64_t)err_ftell(fp) != total_size) {
        fprintf(stderr, "Warning: Total size mismatch (expected %lu, got %ld)\n",
                (unsigned long)total_size, err_ftell(fp));
    }

    err_fclose(fp);
    return 0;
}

tronko_db_t *load_binary(const char *filename, int verbose) {
    FILE *fp = err_xopen(filename, "rb");

    // Read and validate magic
    uint8_t magic[4];
    err_fread(magic, 1, 4, fp);
    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        fprintf(stderr, "Error: Invalid magic number\n");
        err_fclose(fp);
        return NULL;
    }

    // Read version
    uint8_t version_major = fgetc(fp);
    uint8_t version_minor = fgetc(fp);
    if (version_major > TRONKO_VERSION_MAJOR) {
        fprintf(stderr, "Error: Unsupported format version %d.%d (max supported: %d.%d)\n",
                version_major, version_minor, TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR);
        err_fclose(fp);
        return NULL;
    }

    if (verbose) {
        fprintf(stderr, "  Binary format v%d.%d\n", version_major, version_minor);
    }

    // Read endianness and precision flags
    uint8_t endianness = fgetc(fp);
    uint8_t precision = fgetc(fp);
    if (endianness != 0x01) {
        fprintf(stderr, "Error: Only little-endian format is supported\n");
        err_fclose(fp);
        return NULL;
    }
    if (precision != 0x01) {
        fprintf(stderr, "Error: Only float precision is supported\n");
        err_fclose(fp);
        return NULL;
    }

    // Skip header CRC and reserved
    read_u32(fp);  // header_crc (could validate)
    read_u32(fp);  // reserved

    // Read section offsets
    uint64_t taxonomy_offset = read_u64(fp);
    uint64_t node_offset = read_u64(fp);
    uint64_t posterior_offset = read_u64(fp);
    uint64_t total_size = read_u64(fp);

    // Skip reserved bytes to reach global metadata
    err_fseek(fp, FILE_HEADER_SIZE, SEEK_SET);

    // Read global metadata
    int32_t num_trees = read_i32(fp);
    int32_t max_nodename = read_i32(fp);
    int32_t max_tax_name = read_i32(fp);
    int32_t max_line_taxonomy = read_i32(fp);

    if (verbose) {
        fprintf(stderr, "  %d trees, max_nodename=%d, max_tax=%d\n",
                num_trees, max_nodename, max_tax_name);
        fprintf(stderr, "  Offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu, total=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset, (unsigned long)total_size);
    }

    // Allocate database
    tronko_db_t *db = alloc_db(num_trees);
    if (!db) {
        err_fclose(fp);
        return NULL;
    }

    db->max_nodename = max_nodename;
    db->max_tax_name = max_tax_name;
    db->max_line_taxonomy = max_line_taxonomy;

    // Read tree metadata
    for (int i = 0; i < num_trees; i++) {
        db->trees[i].numbase = read_i32(fp);
        db->trees[i].root = read_i32(fp);
        db->trees[i].numspec = read_i32(fp);
        db->trees[i].num_nodes = 2 * db->trees[i].numspec - 1;

        if (verbose) {
            fprintf(stderr, "  Tree %d: numbase=%d, root=%d, numspec=%d\n",
                    i, db->trees[i].numbase, db->trees[i].root, db->trees[i].numspec);
        }
    }

    // === Read Taxonomy Section ===
    err_fseek(fp, taxonomy_offset, SEEK_SET);

    for (int t = 0; t < num_trees; t++) {
        uint32_t tree_tax_size = read_u32(fp);
        read_u32(fp);  // reserved
        (void)tree_tax_size;  // We don't need this, read strings directly

        // Allocate taxonomy
        db->trees[t].taxonomy = calloc(db->trees[t].numspec, sizeof(char **));
        if (!db->trees[t].taxonomy) goto error;

        for (int s = 0; s < db->trees[t].numspec; s++) {
            db->trees[t].taxonomy[s] = calloc(TAXONOMY_LEVELS, sizeof(char *));
            if (!db->trees[t].taxonomy[s]) goto error;

            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                uint16_t len = read_u16(fp);
                db->trees[t].taxonomy[s][l] = calloc(len, sizeof(char));
                if (!db->trees[t].taxonomy[s][l]) goto error;
                err_fread(db->trees[t].taxonomy[s][l], 1, len, fp);
            }
        }
    }

    // === Read Node Section ===
    err_fseek(fp, node_offset, SEEK_SET);

    for (int t = 0; t < num_trees; t++) {
        uint32_t num_nodes = read_u32(fp);
        if ((int32_t)num_nodes != db->trees[t].num_nodes) {
            fprintf(stderr, "Error: Node count mismatch for tree %d\n", t);
            goto error;
        }

        // Allocate nodes
        db->trees[t].nodes = calloc(num_nodes, sizeof(tronko_node_t));
        if (!db->trees[t].nodes) goto error;

        // Store name offsets for later
        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) goto error;

        // Read node records
        for (int n = 0; n < (int)num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            node->up[0] = read_i32(fp);
            node->up[1] = read_i32(fp);
            node->down = read_i32(fp);
            node->depth = read_i32(fp);
            node->taxIndex[0] = read_i32(fp);
            node->taxIndex[1] = read_i32(fp);
            name_offsets[n] = read_u32(fp);
            read_u32(fp);  // reserved
        }

        // Read name table
        long name_table_start = err_ftell(fp);
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                err_fseek(fp, name_table_start + name_offsets[n] - 1, SEEK_SET);
                uint16_t len = read_u16(fp);
                db->trees[t].nodes[n].name = calloc(len, sizeof(char));
                if (!db->trees[t].nodes[n].name) {
                    free(name_offsets);
                    goto error;
                }
                err_fread(db->trees[t].nodes[n].name, 1, len, fp);
            }
        }

        free(name_offsets);
    }

    // === Read Posterior Section ===
    err_fseek(fp, posterior_offset, SEEK_SET);

    if (verbose) {
        fprintf(stderr, "  Reading posteriors...\n");
    }

    for (int t = 0; t < num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            db->trees[t].nodes[n].posteriors = calloc(db->trees[t].numbase * 4, sizeof(float));
            if (!db->trees[t].nodes[n].posteriors) goto error;

            // Bulk read all posteriors for this node
            err_fread(db->trees[t].nodes[n].posteriors,
                     sizeof(float),
                     db->trees[t].numbase * 4,
                     fp);
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: read %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // Optionally validate footer
    uint32_t file_crc = read_u32(fp);
    uint32_t footer_magic = read_u32(fp);
    (void)file_crc;  // Could validate

    if (footer_magic != TRONKO_FOOTER_MAGIC) {
        fprintf(stderr, "Warning: Invalid footer magic (file may be corrupted)\n");
    }

    err_fclose(fp);
    return db;

error:
    free_db(db);
    err_fclose(fp);
    return NULL;
}
