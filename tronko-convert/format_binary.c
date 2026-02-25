#include "format_binary.h"
#include "utils.h"
#include "crc32.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>

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

// === Zstd streaming reader implementation ===

typedef struct {
    FILE *file;
    ZSTD_DCtx *dctx;
    uint8_t *in_buf;
    uint8_t *out_buf;
    size_t in_buf_size;
    size_t out_buf_size;
    size_t in_pos;
    size_t in_end;
    size_t out_pos;
    size_t out_end;
    int eof;
} zstd_reader_t;

static zstd_reader_t *zstd_reader_open(const char *filename) {
    zstd_reader_t *r = calloc(1, sizeof(zstd_reader_t));
    if (!r) return NULL;

    r->file = fopen(filename, "rb");
    if (!r->file) { free(r); return NULL; }

    r->dctx = ZSTD_createDCtx();
    if (!r->dctx) { fclose(r->file); free(r); return NULL; }

    r->in_buf_size = ZSTD_DStreamInSize();
    r->out_buf_size = ZSTD_DStreamOutSize();
    r->in_buf = malloc(r->in_buf_size);
    r->out_buf = malloc(r->out_buf_size);

    if (!r->in_buf || !r->out_buf) {
        if (r->in_buf) free(r->in_buf);
        if (r->out_buf) free(r->out_buf);
        ZSTD_freeDCtx(r->dctx);
        fclose(r->file);
        free(r);
        return NULL;
    }

    return r;
}

static int zstd_reader_refill(zstd_reader_t *r) {
    if (r->eof) return 0;

    if (r->in_pos > 0 && r->in_pos < r->in_end) {
        memmove(r->in_buf, r->in_buf + r->in_pos, r->in_end - r->in_pos);
        r->in_end -= r->in_pos;
        r->in_pos = 0;
    } else {
        r->in_pos = 0;
        r->in_end = 0;
    }

    size_t space = r->in_buf_size - r->in_end;
    if (space > 0) {
        size_t n = fread(r->in_buf + r->in_end, 1, space, r->file);
        r->in_end += n;
        if (n == 0 && r->in_end == 0) {
            r->eof = 1;
            return 0;
        }
    }

    ZSTD_inBuffer in = { r->in_buf, r->in_end, r->in_pos };
    ZSTD_outBuffer out = { r->out_buf, r->out_buf_size, 0 };

    size_t ret = ZSTD_decompressStream(r->dctx, &out, &in);
    if (ZSTD_isError(ret)) {
        fprintf(stderr, "ZSTD decompression error: %s\n", ZSTD_getErrorName(ret));
        return -1;
    }

    r->in_pos = in.pos;
    r->out_pos = 0;
    r->out_end = out.pos;

    return (int)out.pos;
}

static size_t zstd_reader_read(zstd_reader_t *r, void *buf, size_t size) {
    size_t total = 0;
    uint8_t *dst = buf;

    while (total < size) {
        if (r->out_pos >= r->out_end) {
            int ret = zstd_reader_refill(r);
            if (ret <= 0) break;
        }

        size_t avail = r->out_end - r->out_pos;
        size_t need = size - total;
        size_t copy = (avail < need) ? avail : need;

        memcpy(dst + total, r->out_buf + r->out_pos, copy);
        r->out_pos += copy;
        total += copy;
    }

    return total;
}

static void zstd_reader_close(zstd_reader_t *r) {
    if (!r) return;
    ZSTD_freeDCtx(r->dctx);
    free(r->in_buf);
    free(r->out_buf);
    fclose(r->file);
    free(r);
}

static uint16_t zstd_read_u16(zstd_reader_t *r) {
    uint8_t b[2];
    if (zstd_reader_read(r, b, 2) != 2) return 0;
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static uint32_t zstd_read_u32(zstd_reader_t *r) {
    uint8_t b[4];
    if (zstd_reader_read(r, b, 4) != 4) return 0;
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static int32_t zstd_read_i32(zstd_reader_t *r) {
    return (int32_t)zstd_read_u32(r);
}

static uint64_t zstd_read_u64(zstd_reader_t *r) {
    uint8_t b[8];
    if (zstd_reader_read(r, b, 8) != 8) return 0;
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= ((uint64_t)b[i] << (i * 8));
    }
    return v;
}

tronko_db_t *load_binary_zstd(const char *filename, int verbose) {
    zstd_reader_t *r = zstd_reader_open(filename);
    if (!r) {
        fprintf(stderr, "Error: Cannot open zstd file: %s\n", filename);
        return NULL;
    }

    // Read and validate magic
    uint8_t magic[4];
    if (zstd_reader_read(r, magic, 4) != 4) {
        fprintf(stderr, "Error: Failed to read magic number (zstd)\n");
        zstd_reader_close(r);
        return NULL;
    }

    if (magic[0] != TRONKO_MAGIC_0 || magic[1] != TRONKO_MAGIC_1 ||
        magic[2] != TRONKO_MAGIC_2 || magic[3] != TRONKO_MAGIC_3) {
        fprintf(stderr, "Error: Invalid magic number (zstd)\n");
        zstd_reader_close(r);
        return NULL;
    }

    // Read version
    uint8_t version_bytes[2];
    zstd_reader_read(r, version_bytes, 2);
    uint8_t version_major = version_bytes[0];
    uint8_t version_minor = version_bytes[1];

    if (version_major > TRONKO_VERSION_MAJOR) {
        fprintf(stderr, "Error: Unsupported format version %d.%d\n", version_major, version_minor);
        zstd_reader_close(r);
        return NULL;
    }

    if (verbose) {
        fprintf(stderr, "  Zstd binary format v%d.%d\n", version_major, version_minor);
    }

    // Read flags
    uint8_t flags[2];
    zstd_reader_read(r, flags, 2);
    if (flags[0] != 0x01 || flags[1] != 0x01) {
        fprintf(stderr, "Error: Unsupported format flags\n");
        zstd_reader_close(r);
        return NULL;
    }

    // Skip header CRC and reserved
    zstd_read_u32(r);
    zstd_read_u32(r);

    // Read section offsets (for logging)
    uint64_t taxonomy_offset = zstd_read_u64(r);
    uint64_t node_offset = zstd_read_u64(r);
    uint64_t posterior_offset = zstd_read_u64(r);
    uint64_t total_size = zstd_read_u64(r);

    // Skip reserved bytes
    uint8_t skip[16];
    zstd_reader_read(r, skip, 16);

    // Read global metadata
    int32_t num_trees = zstd_read_i32(r);
    int32_t max_nodename = zstd_read_i32(r);
    int32_t max_tax_name = zstd_read_i32(r);
    int32_t max_line_taxonomy = zstd_read_i32(r);

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
        zstd_reader_close(r);
        return NULL;
    }

    db->max_nodename = max_nodename;
    db->max_tax_name = max_tax_name;
    db->max_line_taxonomy = max_line_taxonomy;

    // Read tree metadata
    for (int i = 0; i < num_trees; i++) {
        db->trees[i].numbase = zstd_read_i32(r);
        db->trees[i].root = zstd_read_i32(r);
        db->trees[i].numspec = zstd_read_i32(r);
        db->trees[i].num_nodes = 2 * db->trees[i].numspec - 1;

        if (verbose) {
            fprintf(stderr, "  Tree %d: numbase=%d, root=%d, numspec=%d\n",
                    i, db->trees[i].numbase, db->trees[i].root, db->trees[i].numspec);
        }
    }

    // Read taxonomy
    for (int t = 0; t < num_trees; t++) {
        uint32_t tree_tax_size = zstd_read_u32(r);
        zstd_read_u32(r);  // reserved
        (void)tree_tax_size;

        db->trees[t].taxonomy = calloc(db->trees[t].numspec, sizeof(char **));
        if (!db->trees[t].taxonomy) goto error;

        for (int s = 0; s < db->trees[t].numspec; s++) {
            db->trees[t].taxonomy[s] = calloc(TAXONOMY_LEVELS, sizeof(char *));
            if (!db->trees[t].taxonomy[s]) goto error;

            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                uint16_t len = zstd_read_u16(r);
                db->trees[t].taxonomy[s][l] = calloc(len, sizeof(char));
                if (!db->trees[t].taxonomy[s][l]) goto error;
                zstd_reader_read(r, db->trees[t].taxonomy[s][l], len);
            }
        }
    }

    // Read nodes
    for (int t = 0; t < num_trees; t++) {
        uint32_t num_nodes = zstd_read_u32(r);
        if ((int32_t)num_nodes != db->trees[t].num_nodes) {
            fprintf(stderr, "Error: Node count mismatch for tree %d\n", t);
            goto error;
        }

        db->trees[t].nodes = calloc(num_nodes, sizeof(tronko_node_t));
        if (!db->trees[t].nodes) goto error;

        uint32_t *name_offsets = calloc(num_nodes, sizeof(uint32_t));
        if (!name_offsets) goto error;

        // Read node records
        for (int n = 0; n < (int)num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            node->up[0] = zstd_read_i32(r);
            node->up[1] = zstd_read_i32(r);
            node->down = zstd_read_i32(r);
            node->depth = zstd_read_i32(r);
            node->taxIndex[0] = zstd_read_i32(r);
            node->taxIndex[1] = zstd_read_i32(r);
            name_offsets[n] = zstd_read_u32(r);
            zstd_read_u32(r);  // reserved
        }

        // Read name table (sequentially for streaming)
        for (int n = 0; n < (int)num_nodes; n++) {
            if (name_offsets[n] > 0) {
                uint16_t len = zstd_read_u16(r);
                db->trees[t].nodes[n].name = calloc(len, sizeof(char));
                if (!db->trees[t].nodes[n].name) {
                    free(name_offsets);
                    goto error;
                }
                zstd_reader_read(r, db->trees[t].nodes[n].name, len);
            }
        }

        free(name_offsets);
    }

    // Read posteriors
    if (verbose) {
        fprintf(stderr, "  Reading posteriors (zstd)...\n");
    }

    for (int t = 0; t < num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            db->trees[t].nodes[n].posteriors = calloc(db->trees[t].numbase * 4, sizeof(float));
            if (!db->trees[t].nodes[n].posteriors) goto error;

            zstd_reader_read(r, db->trees[t].nodes[n].posteriors,
                           db->trees[t].numbase * 4 * sizeof(float));
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: read %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    zstd_reader_close(r);
    return db;

error:
    free_db(db);
    zstd_reader_close(r);
    return NULL;
}

// === Zstd streaming writer implementation ===

typedef struct {
    FILE *fp;
    ZSTD_CCtx *cctx;
    uint8_t *out_buf;
    size_t out_buf_size;
} zstd_writer_t;

static zstd_writer_t *zstd_writer_open(const char *filename, int level) {
    zstd_writer_t *w = calloc(1, sizeof(zstd_writer_t));
    if (!w) return NULL;

    w->fp = fopen(filename, "wb");
    if (!w->fp) { free(w); return NULL; }

    w->cctx = ZSTD_createCCtx();
    if (!w->cctx) { fclose(w->fp); free(w); return NULL; }

    ZSTD_CCtx_setParameter(w->cctx, ZSTD_c_compressionLevel, level);
    ZSTD_CCtx_setParameter(w->cctx, ZSTD_c_checksumFlag, 1);

    w->out_buf_size = ZSTD_CStreamOutSize();
    w->out_buf = malloc(w->out_buf_size);
    if (!w->out_buf) {
        ZSTD_freeCCtx(w->cctx);
        fclose(w->fp);
        free(w);
        return NULL;
    }

    return w;
}

static int zstd_writer_write(zstd_writer_t *w, const void *data, size_t size) {
    ZSTD_inBuffer in = { data, size, 0 };

    while (in.pos < in.size) {
        ZSTD_outBuffer out = { w->out_buf, w->out_buf_size, 0 };
        size_t ret = ZSTD_compressStream2(w->cctx, &out, &in, ZSTD_e_continue);
        if (ZSTD_isError(ret)) {
            fprintf(stderr, "ZSTD compression error: %s\n", ZSTD_getErrorName(ret));
            return -1;
        }
        if (out.pos > 0) {
            if (fwrite(w->out_buf, 1, out.pos, w->fp) != out.pos) return -1;
        }
    }
    return 0;
}

static int zstd_writer_close(zstd_writer_t *w) {
    if (!w) return -1;

    // Flush remaining data
    ZSTD_inBuffer in = { NULL, 0, 0 };
    size_t ret;
    do {
        ZSTD_outBuffer out = { w->out_buf, w->out_buf_size, 0 };
        ret = ZSTD_compressStream2(w->cctx, &out, &in, ZSTD_e_end);
        if (ZSTD_isError(ret)) {
            fprintf(stderr, "ZSTD flush error: %s\n", ZSTD_getErrorName(ret));
            break;
        }
        if (out.pos > 0) {
            fwrite(w->out_buf, 1, out.pos, w->fp);
        }
    } while (ret > 0);

    ZSTD_freeCCtx(w->cctx);
    free(w->out_buf);
    fclose(w->fp);
    free(w);
    return 0;
}

// Little-endian write helpers for zstd stream
static int zstd_write_u16(zstd_writer_t *w, uint16_t v) {
    uint8_t b[2] = { v & 0xFF, (v >> 8) & 0xFF };
    return zstd_writer_write(w, b, 2);
}

static int zstd_write_u32(zstd_writer_t *w, uint32_t v) {
    uint8_t b[4] = { v & 0xFF, (v >> 8) & 0xFF, (v >> 16) & 0xFF, (v >> 24) & 0xFF };
    return zstd_writer_write(w, b, 4);
}

static int zstd_write_i32(zstd_writer_t *w, int32_t v) {
    return zstd_write_u32(w, (uint32_t)v);
}

static int zstd_write_u64(zstd_writer_t *w, uint64_t v) {
    uint8_t b[8];
    for (int i = 0; i < 8; i++) b[i] = (v >> (i * 8)) & 0xFF;
    return zstd_writer_write(w, b, 8);
}

int write_binary_zstd(tronko_db_t *db, const char *filename, int compression_level, int verbose) {
    zstd_writer_t *w = zstd_writer_open(filename, compression_level);
    if (!w) {
        fprintf(stderr, "Error: Cannot open '%s' for zstd writing\n", filename);
        return -1;
    }

    if (verbose) {
        fprintf(stderr, "  Writing zstd-compressed binary format v%d.%d (level %d)\n",
                TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR, compression_level);
    }

    // Calculate section offsets (same as write_binary)
    uint64_t header_end = FILE_HEADER_SIZE + GLOBAL_META_SIZE;
    uint64_t tree_meta_size = db->num_trees * TREE_META_SIZE;
    uint64_t taxonomy_offset = header_end + tree_meta_size;

    // Calculate taxonomy section size
    uint64_t taxonomy_size = db->num_trees * 8;
    for (int t = 0; t < db->num_trees; t++) {
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                taxonomy_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
    }

    uint64_t node_offset = taxonomy_offset + taxonomy_size;

    // Calculate node section size
    uint64_t node_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        node_size += 4;
        node_size += db->trees[t].num_nodes * NODE_RECORD_SIZE;
        for (int n = db->trees[t].numspec - 1; n < db->trees[t].num_nodes; n++) {
            if (db->trees[t].nodes[n].name) {
                node_size += 2 + strlen(db->trees[t].nodes[n].name) + 1;
            }
        }
    }

    uint64_t posterior_offset = node_offset + node_size;

    // Calculate posterior section size
    uint64_t posterior_size = 0;
    for (int t = 0; t < db->num_trees; t++) {
        posterior_size += (uint64_t)db->trees[t].num_nodes * db->trees[t].numbase * 4 * sizeof(float);
    }

    uint64_t total_size = posterior_offset + posterior_size + 8;

    if (verbose) {
        fprintf(stderr, "  Section offsets: taxonomy=%lu, nodes=%lu, posteriors=%lu\n",
                (unsigned long)taxonomy_offset, (unsigned long)node_offset,
                (unsigned long)posterior_offset);
    }

    // === Write File Header (64 bytes) ===
    uint8_t header_bytes[8] = {
        TRONKO_MAGIC_0, TRONKO_MAGIC_1, TRONKO_MAGIC_2, TRONKO_MAGIC_3,
        TRONKO_VERSION_MAJOR, TRONKO_VERSION_MINOR,
        0x01,  // Little-endian
        0x01   // Float precision
    };
    zstd_writer_write(w, header_bytes, 8);

    // Header CRC (of bytes 0-7)
    uint32_t header_crc = tronko_crc32(0, header_bytes, 8);
    zstd_write_u32(w, header_crc);
    zstd_write_u32(w, 0);  // Reserved

    zstd_write_u64(w, taxonomy_offset);
    zstd_write_u64(w, node_offset);
    zstd_write_u64(w, posterior_offset);
    zstd_write_u64(w, total_size);

    // Padding to 64 bytes
    uint8_t reserved[16] = {0};
    zstd_writer_write(w, reserved, 16);

    // === Write Global Metadata (16 bytes) ===
    zstd_write_i32(w, db->num_trees);
    zstd_write_i32(w, db->max_nodename);
    zstd_write_i32(w, db->max_tax_name);
    zstd_write_i32(w, db->max_line_taxonomy);

    // === Write Tree Metadata (12 bytes per tree) ===
    for (int i = 0; i < db->num_trees; i++) {
        zstd_write_i32(w, db->trees[i].numbase);
        zstd_write_i32(w, db->trees[i].root);
        zstd_write_i32(w, db->trees[i].numspec);
    }

    // === Write Taxonomy Section ===
    for (int t = 0; t < db->num_trees; t++) {
        uint32_t tree_tax_size = 0;
        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                tree_tax_size += 2 + strlen(db->trees[t].taxonomy[s][l]) + 1;
            }
        }
        zstd_write_u32(w, tree_tax_size);
        zstd_write_u32(w, 0);  // Reserved

        for (int s = 0; s < db->trees[t].numspec; s++) {
            for (int l = 0; l < TAXONOMY_LEVELS; l++) {
                const char *str = db->trees[t].taxonomy[s][l];
                uint16_t len = strlen(str) + 1;
                zstd_write_u16(w, len);
                zstd_writer_write(w, str, len);
            }
        }
    }

    // === Write Node Section ===
    for (int t = 0; t < db->num_trees; t++) {
        zstd_write_u32(w, db->trees[t].num_nodes);

        uint32_t name_offset = 0;

        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            zstd_write_i32(w, node->up[0]);
            zstd_write_i32(w, node->up[1]);
            zstd_write_i32(w, node->down);
            zstd_write_i32(w, node->depth);
            zstd_write_i32(w, node->taxIndex[0]);
            zstd_write_i32(w, node->taxIndex[1]);

            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                zstd_write_u32(w, name_offset + 1);
                name_offset += 2 + strlen(node->name) + 1;
            } else {
                zstd_write_u32(w, 0);
            }
            zstd_write_u32(w, 0);  // Reserved
        }

        // Write name table
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            tronko_node_t *node = &db->trees[t].nodes[n];
            if (node->up[0] == -1 && node->up[1] == -1 && node->name && node->name[0]) {
                uint16_t len = strlen(node->name) + 1;
                zstd_write_u16(w, len);
                zstd_writer_write(w, node->name, len);
            }
        }
    }

    // === Write Posterior Section ===
    if (verbose) {
        fprintf(stderr, "  Writing posteriors...\n");
    }

    for (int t = 0; t < db->num_trees; t++) {
        for (int n = 0; n < db->trees[t].num_nodes; n++) {
            zstd_writer_write(w, db->trees[t].nodes[n].posteriors,
                             db->trees[t].numbase * 4 * sizeof(float));
        }

        if (verbose) {
            fprintf(stderr, "\r  Tree %d: wrote %d nodes", t, db->trees[t].num_nodes);
        }
    }

    if (verbose) {
        fprintf(stderr, "\n");
    }

    // === Write Footer ===
    // For zstd, we skip CRC calculation since zstd has its own checksumming
    zstd_write_u32(w, 0);  // data_crc placeholder
    zstd_write_u32(w, TRONKO_FOOTER_MAGIC);

    zstd_writer_close(w);

    if (verbose) {
        fprintf(stderr, "  Zstd compression complete\n");
    }

    return 0;
}
