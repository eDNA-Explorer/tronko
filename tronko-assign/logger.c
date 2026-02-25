/*
 * logger.c - Implementation of logging infrastructure for tronko-assign
 */

#include "logger.h"
#include "resource_monitor.h"
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/time.h>
#include <errno.h>

// Global logger instance
logger_config_t g_logger = {
    .level = LOG_INFO,
    .log_file = NULL,
    .log_to_stderr = 1,
    .log_to_file = 0,
    .enable_resource_monitoring = 0,
    .enable_timing = 0,
    .log_mutex = PTHREAD_MUTEX_INITIALIZER
};

// Static variables for timing
static struct timeval program_start_time;
static struct timeval last_milestone_time;
static int timing_initialized = 0;

// Milestone names
static const char* milestone_names[] = {
    "STARTUP",
    "OPTIONS_PARSED", 
    "REFERENCE_LOADED",
    "MEMORY_ALLOCATED",
    "BWA_INDEX_BUILT",
    "READ_SPECS_DETECTED",
    "THREADS_INITIALIZED",
    "BATCH_START",
    "BATCH_LOADED",
    "BWA_ALIGNMENT_COMPLETE",
    "DETAILED_ALIGNMENT_COMPLETE",
    "PLACEMENT_COMPLETE",
    "LCA_COMPLETE",
    "RESULTS_WRITTEN",
    "BATCH_COMPLETE",
    "CLEANUP_START",
    "CLEANUP_COMPLETE",
    "PROGRAM_END"
};

// Log level names
static const char* log_level_names[] = {
    "DEBUG", "INFO", "WARN", "ERROR", "NONE"
};

void logger_init(log_level_t level, const char* log_filename, int log_to_stderr) {
    pthread_mutex_lock(&g_logger.log_mutex);
    
    g_logger.level = level;
    g_logger.log_to_stderr = log_to_stderr;
    g_logger.log_to_file = 0;
    g_logger.log_file = NULL;
    
    // Attempt to open log file if specified
    if (log_filename != NULL && strlen(log_filename) > 0) {
        // Basic path validation - reject obviously invalid paths
        if (strlen(log_filename) > 1000) {
            fprintf(stderr, "Warning: Log file path too long, logging to stderr only\n");
        } else {
            g_logger.log_file = fopen(log_filename, "w");
            if (g_logger.log_file != NULL) {
                g_logger.log_to_file = 1;
                // Test that we can actually write to the file
                if (fprintf(g_logger.log_file, "# Log started\n") < 0) {
                    fprintf(stderr, "Warning: Cannot write to log file %s, logging to stderr only\n", log_filename);
                    fclose(g_logger.log_file);
                    g_logger.log_file = NULL;
                    g_logger.log_to_file = 0;
                } else {
                    fflush(g_logger.log_file);
                }
            } else {
                fprintf(stderr, "Warning: Could not open log file %s (errno=%d), logging to stderr only\n", 
                       log_filename, errno);
            }
        }
    }
    
    // Initialize timing with error checking
    if (gettimeofday(&program_start_time, NULL) != 0) {
        fprintf(stderr, "Warning: Failed to initialize logging timestamp\n");
        // Use a default time to prevent crashes
        memset(&program_start_time, 0, sizeof(program_start_time));
    }
    last_milestone_time = program_start_time;
    timing_initialized = 1;
    
    pthread_mutex_unlock(&g_logger.log_mutex);
    
    // Log initialization status (this is safe because mutex is unlocked and stderr is always available)
    LOG_INFO("Logger initialized - Level: %s, File: %s, Stderr: %s", 
             get_log_level_name(level),
             (g_logger.log_to_file && log_filename) ? log_filename : "none",
             log_to_stderr ? "yes" : "no");
}

void logger_cleanup(void) {
    pthread_mutex_lock(&g_logger.log_mutex);
    
    if (g_logger.log_file) {
        fclose(g_logger.log_file);
        g_logger.log_file = NULL;
        g_logger.log_to_file = 0;
    }
    
    pthread_mutex_unlock(&g_logger.log_mutex);
}

void logger_set_level(log_level_t level) {
    pthread_mutex_lock(&g_logger.log_mutex);
    g_logger.level = level;
    pthread_mutex_unlock(&g_logger.log_mutex);
}

void logger_enable_resource_monitoring(int enable) {
    pthread_mutex_lock(&g_logger.log_mutex);
    g_logger.enable_resource_monitoring = enable;
    pthread_mutex_unlock(&g_logger.log_mutex);
}

void logger_enable_timing(int enable) {
    pthread_mutex_lock(&g_logger.log_mutex);
    g_logger.enable_timing = enable;
    pthread_mutex_unlock(&g_logger.log_mutex);
}

void log_current_timestamp(FILE* file) {
    struct timeval tv;
    struct tm* tm_info;
    char timestamp[64];
    
    gettimeofday(&tv, NULL);
    tm_info = localtime(&tv.tv_sec);
    
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
    fprintf(file, "%s.%03ld", timestamp, tv.tv_usec / 1000);
}

static double get_elapsed_time(struct timeval start, struct timeval end) {
    return (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
}

void log_message(log_level_t level, const char* format, ...) {
    if (level < g_logger.level) {
        return;
    }
    
    pthread_mutex_lock(&g_logger.log_mutex);
    
    va_list args;
    char message[1024];
    
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    
    // Log to stderr if enabled
    if (g_logger.log_to_stderr) {
        fprintf(stderr, "[");
        log_current_timestamp(stderr);
        fprintf(stderr, "] [%s] %s\n", get_log_level_name(level), message);
        fflush(stderr);
    }
    
    // Log to file if enabled
    if (g_logger.log_to_file && g_logger.log_file) {
        fprintf(g_logger.log_file, "[");
        log_current_timestamp(g_logger.log_file);
        fprintf(g_logger.log_file, "] [%s] %s\n", get_log_level_name(level), message);
        fflush(g_logger.log_file);
    }
    
    pthread_mutex_unlock(&g_logger.log_mutex);
}

void log_milestone(milestone_t milestone, const char* additional_info) {
    pthread_mutex_lock(&g_logger.log_mutex);
    
    char message[256];
    if (additional_info) {
        snprintf(message, sizeof(message), "MILESTONE: %s - %s", 
                get_milestone_name(milestone), additional_info);
    } else {
        snprintf(message, sizeof(message), "MILESTONE: %s", 
                get_milestone_name(milestone));
    }
    
    // Add resource monitoring if enabled
    if (g_logger.enable_resource_monitoring) {
        resource_stats_t stats;
        if (get_resource_stats(&stats) == 0) {
            char resource_info[128];
            snprintf(resource_info, sizeof(resource_info), 
                    " [Memory: %ld KB, CPU: %.1f%%]", 
                    stats.memory_rss_kb, stats.cpu_percent);
            strncat(message, resource_info, sizeof(message) - strlen(message) - 1);
        }
    }
    
    pthread_mutex_unlock(&g_logger.log_mutex);
    
    log_message(LOG_INFO, "%s", message);
}

void log_milestone_with_timing(milestone_t milestone, const char* additional_info) {
    if (!timing_initialized) {
        log_milestone(milestone, additional_info);
        return;
    }
    
    struct timeval current_time;
    gettimeofday(&current_time, NULL);
    
    double elapsed_total = get_elapsed_time(program_start_time, current_time);
    double elapsed_since_last = get_elapsed_time(last_milestone_time, current_time);
    
    char timing_info[128];
    snprintf(timing_info, sizeof(timing_info), 
            "%.3fs total, %.3fs since last", elapsed_total, elapsed_since_last);
    
    char combined_info[384];
    if (additional_info) {
        snprintf(combined_info, sizeof(combined_info), "%s [%s]", 
                additional_info, timing_info);
    } else {
        snprintf(combined_info, sizeof(combined_info), "[%s]", timing_info);
    }
    
    log_milestone(milestone, combined_info);
    
    last_milestone_time = current_time;
}

void log_thread_start(int thread_id, int start_read, int end_read) {
    LOG_DEBUG("Thread %d starting - processing reads %d to %d", 
              thread_id, start_read, end_read);
}

void log_thread_complete(int thread_id, int reads_processed) {
    LOG_DEBUG("Thread %d completed - processed %d reads", 
              thread_id, reads_processed);
}

const char* get_milestone_name(milestone_t milestone) {
    if (milestone >= 0 && milestone < (sizeof(milestone_names) / sizeof(milestone_names[0]))) {
        return milestone_names[milestone];
    }
    return "UNKNOWN_MILESTONE";
}

const char* get_log_level_name(log_level_t level) {
    if (level >= 0 && level < (sizeof(log_level_names) / sizeof(log_level_names[0]))) {
        return log_level_names[level];
    }
    return "UNKNOWN_LEVEL";
}