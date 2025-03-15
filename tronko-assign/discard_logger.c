#include "discard_logger.h"

// Static variables
static FILE* log_file = NULL;
static int current_log_level = LOG_NONE;
static discard_stats_t stats = {0};
static int is_initialized = 0;

// Helper function to get timestamp
static char* get_timestamp() {
    static char timestamp_buf[32];
    time_t now = time(NULL);
    struct tm* tm_info = localtime(&now);
    strftime(timestamp_buf, sizeof(timestamp_buf), "%Y-%m-%dT%H:%M:%S", tm_info);
    return timestamp_buf;
}

// Helper function to get stage-specific counter
static int* get_stage_counter(const char* stage) {
    if (strcmp(stage, STAGE_BWA_ALIGNMENT) == 0) {
        return &stats.discarded_bwa_alignment;
    } else if (strcmp(stage, STAGE_ALIGNMENT_QUALITY) == 0) {
        return &stats.discarded_alignment_quality;
    } else if (strcmp(stage, STAGE_TREE_PLACEMENT) == 0) {
        return &stats.discarded_tree_placement;
    } else if (strcmp(stage, STAGE_LCA_DETERMINATION) == 0) {
        return &stats.discarded_lca_determination;
    }
    return NULL;
}

int init_discard_logger(const char* log_file_path, int log_level) {
    // If already initialized, close previous log
    if (is_initialized) {
        close_discard_logger();
    }
    
    // Initialize stats
    memset(&stats, 0, sizeof(discard_stats_t));
    
    // Set log level
    current_log_level = log_level;
    
    // If log level is LOG_NONE, just set initialized flag but don't open file
    if (log_level == LOG_NONE) {
        is_initialized = 1;
        return 1;
    }
    
    // Open log file
    log_file = fopen(log_file_path, "w");
    if (log_file == NULL) {
        fprintf(stderr, "Error: Could not open discard log file '%s'\n", log_file_path);
        return 0;
    }
    
    // Write header
    fprintf(log_file, "# Tronko Discarded Reads Log\n");
    fprintf(log_file, "# Created: %s\n", get_timestamp());
    fprintf(log_file, "# Log Level: %d\n\n", log_level);
    
    // Write column headers
    fprintf(log_file, "ReadName\tStage\tReason\tActualValue\tThresholdValue\tAdditionalInfo\tTimestamp\n");
    
    // Flush to ensure header is written
    fflush(log_file);
    
    is_initialized = 1;
    return 1;
}

void log_discarded_read(const char* read_name, const char* stage, const char* reason,
                      double actual_value, double threshold_value, const char* additional_info) {
    // Check if initialized and logging is enabled
    if (!is_initialized || current_log_level == LOG_NONE || log_file == NULL) {
        return;
    }
    
    // Increment discard counter
    stats.total_reads_discarded++;
    
    // Increment stage-specific counter
    int* stage_counter = get_stage_counter(stage);
    if (stage_counter) {
        (*stage_counter)++;
    }
    
    // Basic logging (level 1+)
    if (current_log_level >= LOG_BASIC) {
        fprintf(log_file, "%s\t%s\t%s", 
                read_name ? read_name : "unknown", 
                stage ? stage : "unknown", 
                reason ? reason : "unknown");
        
        // Detailed logging (level 2+): include values and thresholds
        if (current_log_level >= LOG_DETAILED) {
            fprintf(log_file, "\t%.2f\t%.2f", actual_value, threshold_value);
            
            // Verbose logging (level 3+): include additional info
            if (current_log_level >= LOG_VERBOSE) {
                fprintf(log_file, "\t%s", additional_info ? additional_info : "NULL");
            } else {
                fprintf(log_file, "\tNULL");
            }
        } else {
            fprintf(log_file, "\tNULL\tNULL\tNULL");
        }
        
        // Add timestamp
        fprintf(log_file, "\t%s\n", get_timestamp());
    }
    
    // Periodically flush to prevent data loss (every 100 entries)
    if (stats.total_reads_discarded % 100 == 0) {
        fflush(log_file);
    }
}

void increment_processed_reads(int count) {
    if (is_initialized) {
        stats.total_reads_processed += count;
    }
}

void summarize_discarded_reads(FILE* output_stream) {
    if (!is_initialized) {
        return;
    }
    
    // Calculate percentages
    double discard_percent = 0.0;
    double bwa_percent = 0.0;
    double align_percent = 0.0;
    double tree_percent = 0.0;
    double lca_percent = 0.0;
    double success_percent = 0.0;
    
    if (stats.total_reads_processed > 0) {
        discard_percent = (double)stats.total_reads_discarded / stats.total_reads_processed * 100.0;
        bwa_percent = (double)stats.discarded_bwa_alignment / stats.total_reads_processed * 100.0;
        align_percent = (double)stats.discarded_alignment_quality / stats.total_reads_processed * 100.0;
        tree_percent = (double)stats.discarded_tree_placement / stats.total_reads_processed * 100.0;
        lca_percent = (double)stats.discarded_lca_determination / stats.total_reads_processed * 100.0;
        success_percent = 100.0 - discard_percent;
    }
    
    // Write to log file if available
    if (log_file != NULL && current_log_level > LOG_NONE) {
        fprintf(log_file, "\n# Summary Statistics\n");
        fprintf(log_file, "Total_Reads_Processed: %d\n", stats.total_reads_processed);
        fprintf(log_file, "Total_Reads_Discarded: %d (%.1f%%)\n", stats.total_reads_discarded, discard_percent);
        fprintf(log_file, "BWA_Alignment: %d (%.1f%%)\n", stats.discarded_bwa_alignment, bwa_percent);
        fprintf(log_file, "Alignment_Quality: %d (%.1f%%)\n", stats.discarded_alignment_quality, align_percent);
        fprintf(log_file, "Tree_Placement: %d (%.1f%%)\n", stats.discarded_tree_placement, tree_percent);
        fprintf(log_file, "LCA_Determination: %d (%.1f%%)\n", stats.discarded_lca_determination, lca_percent);
        fprintf(log_file, "Reads_Successfully_Assigned: %d (%.1f%%)\n", 
                stats.total_reads_processed - stats.total_reads_discarded, success_percent);
        fflush(log_file);
    }
    
    // Also write to provided output stream if any
    if (output_stream != NULL) {
        fprintf(output_stream, "\n# Discard Statistics\n");
        fprintf(output_stream, "Total Reads Processed: %d\n", stats.total_reads_processed);
        fprintf(output_stream, "Total Reads Discarded: %d (%.1f%%)\n", stats.total_reads_discarded, discard_percent);
        fprintf(output_stream, "  - BWA Alignment: %d (%.1f%%)\n", stats.discarded_bwa_alignment, bwa_percent);
        fprintf(output_stream, "  - Alignment Quality: %d (%.1f%%)\n", stats.discarded_alignment_quality, align_percent);
        fprintf(output_stream, "  - Tree Placement: %d (%.1f%%)\n", stats.discarded_tree_placement, tree_percent);
        fprintf(output_stream, "  - LCA Determination: %d (%.1f%%)\n", stats.discarded_lca_determination, lca_percent);
        fprintf(output_stream, "Reads Successfully Assigned: %d (%.1f%%)\n", 
                stats.total_reads_processed - stats.total_reads_discarded, success_percent);
    }
}

discard_stats_t get_discard_stats() {
    return stats;
}

void close_discard_logger() {
    if (is_initialized && log_file != NULL) {
        // Write final summary
        summarize_discarded_reads(NULL);
        
        // Close file
        fclose(log_file);
        log_file = NULL;
    }
    
    is_initialized = 0;
}