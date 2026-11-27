/*
 * parquet_writer.c
 * Parquet output support for tronko-assign using Carquet library
 */
#ifdef ENABLE_PARQUET

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "global.h"
#include "hashmap.h"
#include "parquet_writer.h"
#include <carquet/carquet.h>

/* One row of a release's <marker>_lineage_taxids.txt; the map's key is path, its data &taxid. */
struct taxid_entry {
    int64_t taxid;
    char path[];
};

struct parquet_writer {
    carquet_writer_t *cw;
    carquet_schema_t *schema;
    char filename[MAXFILENAME];
    int has_taxids;                 /* a Taxid column, from --taxids */
    HASHMAP(char, int64_t) taxids;  /* Taxonomic_Path -> NCBI taxid */
};

/* Read <marker>_lineage_taxids.txt (path, rank, taxid per line, tab-separated; a path
 * with no taxid has an empty one) into pw->taxids. Returns 0, or -1 with err set. */
static int load_taxids(parquet_writer_t *pw, const char *path, char *err) {
    FILE *f = fopen(path, "r");
    if (!f) {
        snprintf(err, 256, "Cannot open taxids file %.200s", path);
        return -1;
    }
    hashmap_init(&pw->taxids, hashmap_hash_string, strcmp);
    pw->has_taxids = 1;
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    long lineno = 0;
    int status = 0;
    while (status == 0 && (len = getline(&line, &cap, f)) > 0) {
        lineno++;
        if (line[len - 1] == '\n') line[--len] = '\0';
        char *rank = strchr(line, '\t');
        char *id = rank ? strchr(rank + 1, '\t') : NULL;
        if (!id || strchr(id + 1, '\t')) {
            snprintf(err, 256, "%.150s line %ld: expected path, rank and taxid", path, lineno);
            status = -1;
            break;
        }
        *rank = '\0';
        if (*++id == '\0') continue;  /* a path with no taxid */
        char *end;
        long long taxid = strtoll(id, &end, 10);
        if (*end != '\0' || taxid <= 0) {
            snprintf(err, 256, "%.150s line %ld: taxid '%.40s' is not a positive integer", path, lineno, id);
            status = -1;
            break;
        }
        struct taxid_entry *e = malloc(sizeof *e + strlen(line) + 1);
        if (!e) {
            snprintf(err, 256, "Out of memory reading %.200s", path);
            status = -1;
            break;
        }
        e->taxid = taxid;
        strcpy(e->path, line);
        if (hashmap_put(&pw->taxids, e->path, &e->taxid) != 0) {
            snprintf(err, 256, "%.150s line %ld: path listed twice", path, lineno);
            free(e);
            status = -1;
        }
    }
    free(line);
    fclose(f);
    return status;
}

static void free_taxids(parquet_writer_t *pw) {
    if (!pw->has_taxids) return;
    int64_t *taxid;
    hashmap_foreach_data(taxid, &pw->taxids) {
        free(taxid);  /* taxid is the first member of its entry */
    }
    hashmap_cleanup(&pw->taxids);
}

/* Create schema matching tronko-assign output columns */
static carquet_schema_t* create_assignment_schema(int has_taxids, carquet_error_t *err) {
    carquet_schema_t *schema = carquet_schema_create(err);
    if (!schema) return NULL;

    /* All columns are required (non-nullable) */
    if (carquet_schema_add_column(schema, "Readname",
            CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Taxonomic_Path",
            CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Score",
            CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Forward_Mismatch",
            CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Reverse_Mismatch",
            CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Tree_Number",
            CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    if (carquet_schema_add_column(schema, "Node_Number",
            CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    /* Null where the path has no taxid: unassigned reads, or an all-NA path */
    if (has_taxids && carquet_schema_add_column(schema, "Taxid",
            CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_OPTIONAL, 0) != CARQUET_OK) {
        carquet_schema_free(schema);
        return NULL;
    }

    return schema;
}

parquet_writer_t* parquet_writer_create(const char *filename, const char *taxids_path, char *err) {
    parquet_writer_t *pw = calloc(1, sizeof(parquet_writer_t));
    if (!pw) {
        snprintf(err, 256, "Failed to allocate parquet_writer");
        return NULL;
    }

    strncpy(pw->filename, filename, MAXFILENAME - 1);
    pw->filename[MAXFILENAME - 1] = '\0';

    if (taxids_path && load_taxids(pw, taxids_path, err) != 0) {
        free_taxids(pw);
        free(pw);
        return NULL;
    }

    carquet_error_t cerr = CARQUET_ERROR_INIT;

    /* Create schema */
    pw->schema = create_assignment_schema(pw->has_taxids, &cerr);
    if (!pw->schema) {
        snprintf(err, 256, "Failed to create schema: %s", cerr.message);
        free_taxids(pw);
        free(pw);
        return NULL;
    }

    /* Configure writer options */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.compression_level = 3;  /* Balance speed/ratio */
    opts.row_group_size = 64 * 1024 * 1024;  /* 64 MB row groups */
    opts.write_statistics = true;

    /* Create writer */
    pw->cw = carquet_writer_create(filename, pw->schema, &opts, &cerr);
    if (!pw->cw) {
        snprintf(err, 256, "Failed to create writer: %s", cerr.message);
        carquet_schema_free(pw->schema);
        free_taxids(pw);
        free(pw);
        return NULL;
    }

    return pw;
}

int parquet_writer_write_batch(parquet_writer_t *writer,
                                const struct assignmentResult *results,
                                int count) {
    if (!writer || !results || count <= 0) return -1;

    /* Allocate temporary arrays for columnar data */
    carquet_byte_array_t *readnames = malloc(count * sizeof(carquet_byte_array_t));
    carquet_byte_array_t *taxon_paths = malloc(count * sizeof(carquet_byte_array_t));
    double *scores = malloc(count * sizeof(double));
    double *fwd_mismatches = malloc(count * sizeof(double));
    double *rev_mismatches = malloc(count * sizeof(double));
    int32_t *tree_numbers = malloc(count * sizeof(int32_t));
    int32_t *node_numbers = malloc(count * sizeof(int32_t));

    if (!readnames || !taxon_paths || !scores || !fwd_mismatches ||
        !rev_mismatches || !tree_numbers || !node_numbers) {
        free(readnames); free(taxon_paths); free(scores);
        free(fwd_mismatches); free(rev_mismatches);
        free(tree_numbers); free(node_numbers);
        return -1;
    }

    /* Convert row-oriented to columnar */
    for (int i = 0; i < count; i++) {
        readnames[i].length = results[i].readname ? strlen(results[i].readname) : 0;
        readnames[i].data = (uint8_t*)results[i].readname;

        taxon_paths[i].length = results[i].taxonomic_path ? strlen(results[i].taxonomic_path) : 0;
        taxon_paths[i].data = (uint8_t*)results[i].taxonomic_path;

        scores[i] = results[i].score;
        fwd_mismatches[i] = results[i].forward_mismatch;
        rev_mismatches[i] = results[i].reverse_mismatch;
        tree_numbers[i] = results[i].tree_number;
        node_numbers[i] = results[i].node_number;
    }

    /* Write each column */
    int failed = 0;
    if (carquet_writer_write_batch(writer->cw, 0, readnames, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 1, taxon_paths, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 2, scores, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 3, fwd_mismatches, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 4, rev_mismatches, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 5, tree_numbers, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (carquet_writer_write_batch(writer->cw, 6, node_numbers, count, NULL, NULL) != CARQUET_OK) failed = 1;
    if (writer->has_taxids) {
        int64_t *taxids = calloc(count, sizeof(int64_t));
        int16_t *present = calloc(count, sizeof(int16_t));
        if (!taxids || !present) {
            failed = 1;
        } else {
            for (int i = 0; i < count; i++) {
                int64_t *taxid = results[i].taxonomic_path
                    ? hashmap_get(&writer->taxids, results[i].taxonomic_path) : NULL;
                if (taxid) {
                    taxids[i] = *taxid;
                    present[i] = 1;
                }
            }
            if (carquet_writer_write_batch(writer->cw, 7, taxids, count, present, NULL) != CARQUET_OK) failed = 1;
        }
        free(taxids); free(present);
    }

    free(readnames); free(taxon_paths); free(scores);
    free(fwd_mismatches); free(rev_mismatches);
    free(tree_numbers); free(node_numbers);

    return failed ? -1 : 0;
}

int parquet_writer_close(parquet_writer_t *writer) {
    if (!writer) return -1;

    carquet_status_t status = carquet_writer_close(writer->cw);
    carquet_schema_free(writer->schema);
    free_taxids(writer);
    free(writer);

    return (status == CARQUET_OK) ? 0 : -1;
}

#endif /* ENABLE_PARQUET */
