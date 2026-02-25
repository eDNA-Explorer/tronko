/*
 * resource_monitor.c - Implementation of resource monitoring for tronko-assign
 */

#include "resource_monitor.h"
#include "logger.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <errno.h>

// Static variables for baseline measurements
static resource_stats_t baseline_stats;
static int monitoring_initialized = 0;
static struct timeval process_start_time;

// Internal function to parse /proc/self/status
static int parse_proc_status(resource_stats_t* stats) {
    FILE* status_file = fopen("/proc/self/status", "r");
    if (!status_file) {
        // /proc might not be available (container, chroot, permissions, etc.)
        // This is not critical, continue with values from getrusage()
        return -1;
    }
    
    char line[256];
    int lines_read = 0;
    const int MAX_LINES_TO_READ = 150; // Prevent potential infinite loops or hangs
    
    while (fgets(line, sizeof(line), status_file) && lines_read < MAX_LINES_TO_READ) {
        lines_read++;
        
        if (strncmp(line, "VmRSS:", 6) == 0) {
            sscanf(line, "VmRSS: %ld kB", &stats->memory_rss_kb);
        } else if (strncmp(line, "VmSize:", 7) == 0) {
            sscanf(line, "VmSize: %ld kB", &stats->memory_vm_size_kb);
        } else if (strncmp(line, "VmPeak:", 7) == 0) {
            sscanf(line, "VmPeak: %ld kB", &stats->memory_vm_peak_kb);
        } else if (strncmp(line, "VmHWM:", 6) == 0) {
            sscanf(line, "VmHWM: %ld kB", &stats->memory_vm_rss_peak_kb);
        } else if (strncmp(line, "voluntary_ctxt_switches:", 24) == 0) {
            sscanf(line, "voluntary_ctxt_switches: %ld", &stats->voluntary_ctx_switches);
        } else if (strncmp(line, "nonvoluntary_ctxt_switches:", 27) == 0) {
            sscanf(line, "nonvoluntary_ctxt_switches: %ld", &stats->involuntary_ctx_switches);
        }
    }
    
    fclose(status_file);
    
    // Log if we hit the safety limit
    if (lines_read >= MAX_LINES_TO_READ) {
        LOG_WARN("Reached maximum line limit while reading /proc/self/status");
    }
    
    return 0;
}

// Internal function to parse /proc/self/io
static int parse_proc_io(resource_stats_t* stats) {
    FILE* io_file = fopen("/proc/self/io", "r");
    if (!io_file) {
        // /proc/self/io may not be available, not an error
        stats->io_read_bytes = 0;
        stats->io_write_bytes = 0;
        return 0;
    }
    
    char line[256];
    int lines_read = 0;
    const int MAX_LINES_TO_READ = 20; // /proc/self/io is usually just a few lines
    
    while (fgets(line, sizeof(line), io_file) && lines_read < MAX_LINES_TO_READ) {
        lines_read++;
        
        if (strncmp(line, "read_bytes:", 11) == 0) {
            sscanf(line, "read_bytes: %ld", &stats->io_read_bytes);
        } else if (strncmp(line, "write_bytes:", 12) == 0) {
            sscanf(line, "write_bytes: %ld", &stats->io_write_bytes);
        }
    }
    
    fclose(io_file);
    
    // Log if we hit the safety limit (unlikely for /proc/self/io)
    if (lines_read >= MAX_LINES_TO_READ) {
        LOG_WARN("Reached maximum line limit while reading /proc/self/io");
    }
    
    return 0;
}

int init_resource_monitoring(void) {
    if (monitoring_initialized) {
        return 0;
    }
    
    // Initialize timing
    if (gettimeofday(&process_start_time, NULL) != 0) {
        LOG_WARN("Failed to get process start time for resource monitoring");
        // Continue anyway, timing will be less accurate
        memset(&process_start_time, 0, sizeof(process_start_time));
    }
    
    // Get baseline resource statistics - always succeeds now with defensive code
    get_resource_stats(&baseline_stats);
    
    monitoring_initialized = 1;
    LOG_DEBUG("Resource monitoring initialized");
    return 0;
}

void cleanup_resource_monitoring(void) {
    monitoring_initialized = 0;
}

int get_resource_stats(resource_stats_t* stats) {
    if (!stats) {
        return -1;
    }
    
    // Initialize all fields to 0
    memset(stats, 0, sizeof(resource_stats_t));
    
    // Get resource usage from getrusage
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        stats->user_time_sec = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1000000.0;
        stats->system_time_sec = usage.ru_stime.tv_sec + usage.ru_stime.tv_usec / 1000000.0;
        stats->minor_page_faults = usage.ru_minflt;
        stats->major_page_faults = usage.ru_majflt;
        
        // Convert maxrss to KB (on Linux it's already in KB, on macOS it's in bytes)
#ifdef __APPLE__
        stats->memory_rss_kb = usage.ru_maxrss / 1024;
#else
        stats->memory_rss_kb = usage.ru_maxrss;
#endif
    }
    
    // Get additional memory information from /proc/self/status (Linux only)
    parse_proc_status(stats);
    
    // Get I/O information from /proc/self/io (Linux only)
    parse_proc_io(stats);
    
    // Calculate wall time
    if (monitoring_initialized) {
        struct timeval current_time;
        gettimeofday(&current_time, NULL);
        stats->wall_time_sec = (current_time.tv_sec - process_start_time.tv_sec) + 
                              (current_time.tv_usec - process_start_time.tv_usec) / 1000000.0;
        
        // Calculate approximate CPU percentage
        if (stats->wall_time_sec > 0) {
            stats->cpu_percent = ((stats->user_time_sec + stats->system_time_sec) / stats->wall_time_sec) * 100.0;
        }
    }
    
    return 0;
}

void print_resource_stats(const resource_stats_t* stats, const char* label) {
    printf("=== Resource Stats: %s ===\n", label ? label : "Current");
    printf("Memory RSS: %ld KB (%.1f MB)\n", stats->memory_rss_kb, stats->memory_rss_kb / 1024.0);
    printf("Memory VM Size: %ld KB (%.1f MB)\n", stats->memory_vm_size_kb, stats->memory_vm_size_kb / 1024.0);
    printf("Memory VM Peak: %ld KB (%.1f MB)\n", stats->memory_vm_peak_kb, stats->memory_vm_peak_kb / 1024.0);
    printf("CPU Time: %.3fs user, %.3fs system\n", stats->user_time_sec, stats->system_time_sec);
    printf("CPU Percent: %.1f%%\n", stats->cpu_percent);
    printf("Wall Time: %.3fs\n", stats->wall_time_sec);
    printf("Page Faults: %ld minor, %ld major\n", stats->minor_page_faults, stats->major_page_faults);
    printf("Context Switches: %ld voluntary, %ld involuntary\n", stats->voluntary_ctx_switches, stats->involuntary_ctx_switches);
    printf("I/O: %ld bytes read, %ld bytes written\n", stats->io_read_bytes, stats->io_write_bytes);
    printf("=====================================\n");
}

void log_resource_stats(const resource_stats_t* stats, const char* label) {
    LOG_INFO("Resource Stats [%s]: Memory RSS=%.1fMB, VM=%.1fMB, CPU=%.1f%%, Wall=%.3fs", 
             label ? label : "Current",
             stats->memory_rss_kb / 1024.0,
             stats->memory_vm_size_kb / 1024.0,
             stats->cpu_percent,
             stats->wall_time_sec);
}

double get_memory_usage_mb(const resource_stats_t* stats) {
    return stats->memory_rss_kb / 1024.0;
}

double get_cpu_usage_percent(const resource_stats_t* stats) {
    return stats->cpu_percent;
}

int start_operation_monitoring(operation_monitor_t* monitor, const char* operation_name) {
    if (!monitor) {
        return -1;
    }
    
    strncpy(monitor->operation_name, operation_name ? operation_name : "Unknown", 
            sizeof(monitor->operation_name) - 1);
    monitor->operation_name[sizeof(monitor->operation_name) - 1] = '\0';
    
    return get_resource_stats(&monitor->start_stats);
}

int end_operation_monitoring(operation_monitor_t* monitor) {
    if (!monitor) {
        return -1;
    }
    
    return get_resource_stats(&monitor->end_stats);
}

void log_operation_stats(const operation_monitor_t* monitor) {
    if (!monitor) {
        return;
    }
    
    long memory_delta = monitor->end_stats.memory_rss_kb - monitor->start_stats.memory_rss_kb;
    double time_delta = monitor->end_stats.wall_time_sec - monitor->start_stats.wall_time_sec;
    double cpu_time_delta = (monitor->end_stats.user_time_sec + monitor->end_stats.system_time_sec) -
                           (monitor->start_stats.user_time_sec + monitor->start_stats.system_time_sec);
    
    LOG_INFO("Operation [%s]: Time=%.3fs, Memory Δ=%+ldKB, CPU Time=%.3fs", 
             monitor->operation_name, time_delta, memory_delta, cpu_time_delta);
}

void log_memory_allocation(const char* structure_name, size_t size_bytes) {
    LOG_DEBUG("Memory allocated: %s - %zu bytes (%.1f KB)", 
              structure_name, size_bytes, size_bytes / 1024.0);
}

void log_memory_deallocation(const char* structure_name, size_t size_bytes) {
    LOG_DEBUG("Memory deallocated: %s - %zu bytes (%.1f KB)", 
              structure_name, size_bytes, size_bytes / 1024.0);
}

void log_current_resource_usage(const char* context) {
    resource_stats_t stats;
    if (get_resource_stats(&stats) == 0) {
        log_resource_stats(&stats, context);
    } else {
        LOG_WARN("Failed to get resource stats for context: %s", context);
    }
}

void log_resource_delta(const resource_stats_t* before, const resource_stats_t* after, const char* operation) {
    if (!before || !after) {
        return;
    }
    
    long memory_delta = after->memory_rss_kb - before->memory_rss_kb;
    double time_delta = after->wall_time_sec - before->wall_time_sec;
    double cpu_delta = (after->user_time_sec + after->system_time_sec) - 
                      (before->user_time_sec + before->system_time_sec);
    long minor_faults_delta = after->minor_page_faults - before->minor_page_faults;
    long major_faults_delta = after->major_page_faults - before->major_page_faults;
    
    LOG_INFO("Resource Delta [%s]: Time=%.3fs, Memory Δ=%+ldKB, CPU=%.3fs, Page Faults=%ld+%ld", 
             operation ? operation : "Operation",
             time_delta, memory_delta, cpu_delta, minor_faults_delta, major_faults_delta);
}