#ifndef DISCARD_LOGGER_H
#define DISCARD_LOGGER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Logging levels
#define LOG_NONE 0   // No logging
#define LOG_BASIC 1  // Basic info (stage, reason)
#define LOG_DETAILED 2  // Detailed info (includes values, thresholds)
#define LOG_VERBOSE 3   // Most verbose (includes additional context)

// Filtering stages
#define STAGE_BWA_ALIGNMENT "BWA_Alignment"
#define STAGE_ALIGNMENT_QUALITY "Alignment_Quality"
#define STAGE_TREE_PLACEMENT "Tree_Placement"
#define STAGE_LCA_DETERMINATION "LCA_Determination"

// Common reasons
#define REASON_NO_BWA_HIT "No_BWA_hit"
#define REASON_LOW_MAPQ "Low_mapping_quality"
#define REASON_LOW_ALIGNMENT_SCORE "Low_alignment_score"
#define REASON_SHORT_ALIGNMENT "Short_alignment_length"
#define REASON_EXCESSIVE_GAPS "Excessive_gaps"
#define REASON_AMBIGUOUS_PLACEMENT "Ambiguous_placement"
#define REASON_LOW_LIKELIHOOD "Low_likelihood_value"
#define REASON_INCONSISTENT_PAIRS "Inconsistent_paired_placements"
#define REASON_BELOW_LCA_CUTOFF "Below_LCA_cutoff"
#define REASON_HIGH_TAXONOMIC_LEVEL "High_taxonomic_level"

// Discard logger statistics structure
typedef struct {
    int total_reads_processed;
    int total_reads_discarded;
    int discarded_bwa_alignment;
    int discarded_alignment_quality;
    int discarded_tree_placement;
    int discarded_lca_determination;
} discard_stats_t;

/**
 * Initialize the discard logger
 * @param log_file_path Path to the log file to create/open
 * @param log_level Verbosity level (0-3)
 * @return 1 if successful, 0 if failed
 */
int init_discard_logger(const char* log_file_path, int log_level);

/**
 * Log a discarded read with detailed information
 * @param read_name Name/ID of the read being discarded
 * @param stage Processing stage where read was discarded
 * @param reason Specific reason for discard
 * @param actual_value Measured value that caused discard
 * @param threshold_value Threshold value that was not met
 * @param additional_info Any additional context (can be NULL)
 */
void log_discarded_read(const char* read_name, const char* stage, const char* reason,
                        double actual_value, double threshold_value,
                        const char* additional_info);

/**
 * Increment the processed reads counter
 * @param count Number of reads to add to count (usually 1)
 */
void increment_processed_reads(int count);

/**
 * Write summary statistics for discarded reads to the log file
 * and optionally to an output stream
 * @param output_stream Additional stream to write summary (can be NULL)
 */
void summarize_discarded_reads(FILE* output_stream);

/**
 * Get current discard statistics
 * @return Struct containing discard statistics
 */
discard_stats_t get_discard_stats();

/**
 * Clean up and close the logging system
 */
void close_discard_logger();

#endif /* DISCARD_LOGGER_H */