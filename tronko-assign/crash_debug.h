/*
 * crash_debug.h - Comprehensive crash debugging system for tronko-assign
 * Provides signal handling, stack traces, and crash forensics
 */

#ifndef _CRASH_DEBUG_H_
#define _CRASH_DEBUG_H_

#include <signal.h>
#include <ucontext.h>
#include <sys/types.h>
#include <time.h>
#include <stdio.h>
#include <stddef.h>

#ifdef __APPLE__
#include <execinfo.h>
#include <mach/mach.h>
#include <mach/thread_info.h>
#endif

#ifdef __linux__
#include <execinfo.h>
#include <sys/ucontext.h>
#endif

// Maximum stack trace depth
#define MAX_STACK_TRACE_DEPTH 128
#define MAX_CRASH_MESSAGE_LEN 4096
#define MAX_SYMBOL_NAME_LEN 256

// Crash information structure
typedef struct {
    int signal_number;
    int signal_code;
    void* fault_address;
    pid_t process_id;
    pid_t thread_id;
    time_t crash_time;
    char crash_message[MAX_CRASH_MESSAGE_LEN];
    
    // Stack trace information
    void* stack_trace[MAX_STACK_TRACE_DEPTH];
    int stack_depth;
    char** symbol_names;
    
    // Register state (platform specific)
#ifdef __APPLE__
    mcontext_t machine_context;
#endif
#ifdef __linux__
    mcontext_t machine_context;
#endif
    
    // Memory state at crash
    size_t memory_rss_kb;
    size_t memory_vm_kb;
    double cpu_percent;
    
    // Additional context
    char program_name[256];
    char working_directory[1024];
    char command_line[2048];
    
    // Application-specific context
    char current_operation[256];      // What tronko was doing when it crashed
    char current_file[512];           // File being processed
    int current_file_line;            // Line number in current file
    char current_read_name[256];      // Current read/sequence name
    int current_batch;                // Batch number
    int current_read_index;           // Read index within batch
    int current_tree;                 // Tree being processed
    char processing_stage[128];       // BWA, alignment, placement, LCA, etc.
} crash_info_t;

// Crash handler configuration
typedef struct {
    int enable_stack_trace;
    int enable_register_dump;
    int enable_memory_dump;
    int enable_core_dump;
    int enable_crash_log;
    char crash_log_directory[1024];
    char crash_log_prefix[256];
    void (*user_crash_handler)(const crash_info_t* crash_info);
} crash_config_t;

// Function declarations
int crash_debug_init(const crash_config_t* config);
void crash_debug_cleanup(void);
void crash_debug_set_user_handler(void (*handler)(const crash_info_t* crash_info));
void crash_debug_enable_core_dumps(void);
void crash_debug_set_crash_log_dir(const char* directory);

// Signal handling functions
void crash_signal_handler(int sig, siginfo_t* info, void* context);
void crash_setup_signal_handlers(void);
void crash_restore_default_handlers(void);

// Stack trace functions
int crash_capture_stack_trace(void** trace, int max_depth);
char** crash_resolve_symbols(void* const* trace, int depth);
void crash_print_stack_trace(const void* const* trace, int depth, FILE* output);
void crash_log_stack_trace(const void* const* trace, int depth);

// Register and memory state functions
void crash_capture_register_state(const ucontext_t* context, crash_info_t* crash_info);
void crash_capture_memory_state(crash_info_t* crash_info);
void crash_print_register_state(const crash_info_t* crash_info, FILE* output);
void crash_print_memory_state(const crash_info_t* crash_info, FILE* output);

// Crash report generation
void crash_generate_report(const crash_info_t* crash_info, const char* filename);
void crash_write_crash_log(const crash_info_t* crash_info);
char* crash_format_crash_report(const crash_info_t* crash_info);

// Application context tracking functions
void crash_set_context(const char* operation, const char* details);
void crash_set_current_file(const char* filename);
void crash_set_current_file_line(const char* filename, int line_number);
void crash_set_current_read(const char* read_name, int batch, int read_index);
void crash_set_current_tree(int tree_number);
void crash_set_processing_stage(const char* stage);
void crash_clear_context(void);

// Utility functions
const char* crash_signal_name(int signal);
const char* crash_signal_description(int signal, int code);
void crash_get_program_info(crash_info_t* crash_info);
void crash_safe_string_copy(char* dest, const char* src, size_t dest_size);

// Memory debugging helpers
void crash_check_heap_corruption(void);
void crash_dump_memory_region(void* addr, size_t size, FILE* output);
int crash_is_valid_pointer(void* ptr);

// Thread-safe crash handling
extern volatile sig_atomic_t g_crash_in_progress;
extern crash_config_t g_crash_config;

// Convenience macros
#define CRASH_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            crash_info_t info = {0}; \
            snprintf(info.crash_message, sizeof(info.crash_message), \
                     "ASSERTION FAILED: %s at %s:%d in %s", \
                     message, __FILE__, __LINE__, __func__); \
            crash_generate_report(&info, NULL); \
            abort(); \
        } \
    } while(0)

#define CRASH_CHECK_POINTER(ptr, name) \
    do { \
        if (!crash_is_valid_pointer(ptr)) { \
            crash_info_t info = {0}; \
            snprintf(info.crash_message, sizeof(info.crash_message), \
                     "INVALID POINTER: %s (%p) at %s:%d in %s", \
                     name, ptr, __FILE__, __LINE__, __func__); \
            crash_generate_report(&info, NULL); \
        } \
    } while(0)

#endif /* _CRASH_DEBUG_H_ */