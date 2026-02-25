#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "format_common.h"
#include "format_text.h"
#include "format_binary.h"
#include "utils.h"

static void print_usage(const char *prog) {
    fprintf(stderr, "tronko-convert: Convert tronko reference database formats\n\n");
    fprintf(stderr, "Usage: %s [options]\n", prog);
    fprintf(stderr, "  -i <file>   Input file (required)\n");
    fprintf(stderr, "  -o <file>   Output file (required)\n");
    fprintf(stderr, "  -t          Output as text (default: zstd-compressed binary)\n");
    fprintf(stderr, "  -u          Output as uncompressed binary\n");
    fprintf(stderr, "  -c <level>  Zstd compression level 1-19 (default: 19, max compression)\n");
    fprintf(stderr, "              Lower = faster compression, larger file\n");
    fprintf(stderr, "              Higher = slower compression, smaller file\n");
    fprintf(stderr, "  -v          Verbose output\n");
    fprintf(stderr, "  -h          Show this help\n\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb\n", prog);
    fprintf(stderr, "  %s -i reference_tree.trkb -o reference_tree.txt -t\n", prog);
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb -u  # uncompressed\n", prog);
    fprintf(stderr, "  %s -i reference_tree.txt -o reference_tree.trkb -c 3  # fast compression\n", prog);
}

int main(int argc, char *argv[]) {
    char *input_file = NULL;
    char *output_file = NULL;
    int output_text = 0;
    int output_uncompressed = 0;
    int compression_level = ZSTD_COMPRESSION_LEVEL;
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "i:o:tuc:vh")) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'o': output_file = optarg; break;
            case 't': output_text = 1; break;
            case 'u': output_uncompressed = 1; break;
            case 'c':
                compression_level = atoi(optarg);
                if (compression_level < 1 || compression_level > 19) {
                    fprintf(stderr, "Error: Compression level must be 1-19\n");
                    return 1;
                }
                break;
            case 'v': verbose = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    if (!input_file || !output_file) {
        fprintf(stderr, "Error: Both -i and -o are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    // Detect input format
    int input_format = detect_format(input_file);
    if (input_format < 0) {
        fprintf(stderr, "Error: Cannot open or detect format of '%s'\n", input_file);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Input format: %s\n",
                input_format == FORMAT_BINARY ? "binary (.trkb)" :
                input_format == FORMAT_BINARY_ZSTD ? "zstd binary (.trkb)" : "text");
        fprintf(stderr, "Output format: %s\n",
                output_text ? "text" :
                output_uncompressed ? "uncompressed binary (.trkb)" :
                "zstd binary (.trkb)");
    }

    // Load database
    tronko_db_t *db = NULL;
    if (input_format == FORMAT_BINARY) {
        db = load_binary(input_file, verbose);
    } else if (input_format == FORMAT_BINARY_ZSTD) {
        db = load_binary_zstd(input_file, verbose);
    } else {
        db = load_text(input_file, verbose);
    }

    if (!db) {
        fprintf(stderr, "Error: Failed to load '%s'\n", input_file);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Loaded: %d trees\n", db->num_trees);
        for (int i = 0; i < db->num_trees; i++) {
            fprintf(stderr, "  Tree %d: %d species, %d positions, %d nodes\n",
                    i, db->trees[i].numspec, db->trees[i].numbase, db->trees[i].num_nodes);
        }
    }

    // Write output
    int result;
    if (output_text) {
        result = write_text(db, output_file, verbose);
    } else if (output_uncompressed) {
        result = write_binary(db, output_file, verbose);
    } else {
        result = write_binary_zstd(db, output_file, compression_level, verbose);
    }

    if (result != 0) {
        fprintf(stderr, "Error: Failed to write '%s'\n", output_file);
        free_db(db);
        return 1;
    }

    if (verbose) {
        fprintf(stderr, "Successfully wrote '%s'\n", output_file);
    }

    free_db(db);
    return 0;
}
