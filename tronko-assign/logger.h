/*
 * logger.h - Logging infrastructure for tronko-assign
 * Provides configurable logging with timestamps and resource monitoring
 */

#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdio.h>
#include <time.h>
#include <pthread.h>

// Log levels
typedef enum {
    LOG_DEBUG = 0,
    LOG_INFO = 1,
    LOG_WARN = 2,
    LOG_ERROR = 3,
    LOG_NONE = 4
} log_level_t;

// Milestone identifiers
typedef enum {
    MILESTONE_STARTUP = 0,
    MILESTONE_OPTIONS_PARSED,
    MILESTONE_REFERENCE_LOADED,
    MILESTONE_MEMORY_ALLOCATED,
    MILESTONE_BWA_INDEX_BUILT,
    MILESTONE_READ_SPECS_DETECTED,
    MILESTONE_THREADS_INITIALIZED,
    MILESTONE_BATCH_START,
    MILESTONE_BATCH_LOADED,
    MILESTONE_BWA_ALIGNMENT_COMPLETE,
    MILESTONE_DETAILED_ALIGNMENT_COMPLETE,
    MILESTONE_PLACEMENT_COMPLETE,
    MILESTONE_LCA_COMPLETE,
    MILESTONE_RESULTS_WRITTEN,
    MILESTONE_BATCH_COMPLETE,
    MILESTONE_CLEANUP_START,
    MILESTONE_CLEANUP_COMPLETE,
    MILESTONE_PROGRAM_END
} milestone_t;

// Logger configuration structure
typedef struct {
    log_level_t level;
    FILE* log_file;
    int log_to_stderr;
    int log_to_file;
    int enable_resource_monitoring;
    int enable_timing;
    pthread_mutex_t log_mutex;
} logger_config_t;

// Function declarations
void logger_init(log_level_t level, const char* log_filename, int log_to_stderr);
void logger_cleanup(void);
void logger_set_level(log_level_t level);
void logger_enable_resource_monitoring(int enable);
void logger_enable_timing(int enable);

// Main logging functions
void log_message(log_level_t level, const char* format, ...);
void log_milestone(milestone_t milestone, const char* additional_info);
void log_milestone_with_timing(milestone_t milestone, const char* additional_info);

// Convenience macros
#define LOG_DEBUG(fmt, ...) log_message(LOG_DEBUG, fmt, ##__VA_ARGS__)
#define LOG_INFO(fmt, ...) log_message(LOG_INFO, fmt, ##__VA_ARGS__)
#define LOG_WARN(fmt, ...) log_message(LOG_WARN, fmt, ##__VA_ARGS__)
#define LOG_ERROR(fmt, ...) log_message(LOG_ERROR, fmt, ##__VA_ARGS__)

// Milestone logging macros
#define LOG_MILESTONE(milestone) log_milestone(milestone, NULL)
#define LOG_MILESTONE_INFO(milestone, info) log_milestone(milestone, info)
#define LOG_MILESTONE_TIMED(milestone) log_milestone_with_timing(milestone, NULL)
#define LOG_MILESTONE_TIMED_INFO(milestone, info) log_milestone_with_timing(milestone, info)

// Thread-safe logging for multi-threaded sections
void log_thread_start(int thread_id, int start_read, int end_read);
void log_thread_complete(int thread_id, int reads_processed);

// Utility functions
const char* get_milestone_name(milestone_t milestone);
const char* get_log_level_name(log_level_t level);
void log_current_timestamp(FILE* file);

// Global logger instance
extern logger_config_t g_logger;

#endif /* _LOGGER_H_ */