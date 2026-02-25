/*
 * resource_monitor.h - Memory and CPU resource monitoring for tronko-assign
 * Provides functions to track system resource usage at runtime
 */

#ifndef _RESOURCE_MONITOR_H_
#define _RESOURCE_MONITOR_H_

#include <sys/types.h>
#include <sys/resource.h>

// Resource statistics structure
typedef struct {
    // Memory statistics (in KB)
    long memory_rss_kb;        // Resident Set Size (physical memory currently used)
    long memory_vm_size_kb;    // Virtual Memory Size (total virtual memory used)
    long memory_vm_peak_kb;    // Peak Virtual Memory Size
    long memory_vm_rss_peak_kb; // Peak Resident Set Size
    
    // CPU statistics
    double cpu_percent;        // Current CPU percentage (approximation)
    double user_time_sec;      // User CPU time in seconds
    double system_time_sec;    // System CPU time in seconds
    
    // I/O statistics
    long io_read_bytes;        // Bytes read from storage
    long io_write_bytes;       // Bytes written to storage
    
    // Context switches
    long voluntary_ctx_switches;   // Voluntary context switches
    long involuntary_ctx_switches; // Involuntary context switches
    
    // Page faults
    long minor_page_faults;    // Minor page faults (no I/O required)
    long major_page_faults;    // Major page faults (I/O required)
    
    // Timing
    double wall_time_sec;      // Wall clock time since process start
} resource_stats_t;

// Function declarations
int get_resource_stats(resource_stats_t* stats);
int init_resource_monitoring(void);
void cleanup_resource_monitoring(void);

// Utility functions
void print_resource_stats(const resource_stats_t* stats, const char* label);
void log_resource_stats(const resource_stats_t* stats, const char* label);
double get_memory_usage_mb(const resource_stats_t* stats);
double get_cpu_usage_percent(const resource_stats_t* stats);

// Resource monitoring for specific operations
typedef struct {
    resource_stats_t start_stats;
    resource_stats_t end_stats;
    char operation_name[64];
} operation_monitor_t;

int start_operation_monitoring(operation_monitor_t* monitor, const char* operation_name);
int end_operation_monitoring(operation_monitor_t* monitor);
void log_operation_stats(const operation_monitor_t* monitor);

// Memory usage tracking for specific data structures
void log_memory_allocation(const char* structure_name, size_t size_bytes);
void log_memory_deallocation(const char* structure_name, size_t size_bytes);

// High-level monitoring functions
void log_current_resource_usage(const char* context);
void log_resource_delta(const resource_stats_t* before, const resource_stats_t* after, const char* operation);

#endif /* _RESOURCE_MONITOR_H_ */