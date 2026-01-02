#ifndef FORMAT_COMMON_H
#define FORMAT_COMMON_H

#include <stdint.h>

// Binary format magic and version
#define TRONKO_MAGIC_0 0x89
#define TRONKO_MAGIC_1 'T'
#define TRONKO_MAGIC_2 'R'
#define TRONKO_MAGIC_3 'K'
#define TRONKO_VERSION_MAJOR 1
#define TRONKO_VERSION_MINOR 0
#define TRONKO_FOOTER_MAGIC 0x454E4421  // "END!"

// Format identifiers
#define FORMAT_UNKNOWN        -1
#define FORMAT_TEXT            1
#define FORMAT_BINARY          2
#define FORMAT_BINARY_GZIPPED  3  // Gzipped binary format (.trkb.gz)
#define FORMAT_BINARY_ZSTD     4  // Zstd-compressed binary format

// Default compression level (max compression - write-once, read-many use case)
#define ZSTD_COMPRESSION_LEVEL 19

// Taxonomy levels (fixed at 7: domain, phylum, class, order, family, genus, species)
#define TAXONOMY_LEVELS 7

// Node structure (mirrors tronko-assign/global.h:67-75)
typedef struct {
    int32_t up[2];        // Child indices (-1 for leaf nodes)
    int32_t down;         // Parent index
    int32_t depth;        // Tree depth
    int32_t taxIndex[2];  // [species_index, taxonomy_level]
    char *name;           // Node name (NULL for internal nodes)
    float *posteriors;    // [numbase * 4] contiguous array (A, C, G, T per position)
} tronko_node_t;

// Tree structure
typedef struct {
    int32_t numbase;      // MSA alignment length (number of positions)
    int32_t root;         // Root node index
    int32_t numspec;      // Number of species (leaf nodes)
    int32_t num_nodes;    // Total nodes: 2*numspec - 1
    char ***taxonomy;     // [numspec][TAXONOMY_LEVELS] taxonomy strings
    tronko_node_t *nodes; // [num_nodes] node array
} tronko_tree_t;

// Database structure
typedef struct {
    int32_t num_trees;
    int32_t max_nodename;      // Max node name length
    int32_t max_tax_name;      // Max taxonomy name length
    int32_t max_line_taxonomy; // Max full taxonomy line length
    tronko_tree_t *trees;
} tronko_db_t;

// Format detection
int detect_format(const char *filename);

// Memory management
tronko_db_t *alloc_db(int num_trees);
void free_db(tronko_db_t *db);

#endif
