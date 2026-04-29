/*
 * crash_debug.c - Implementation of comprehensive crash debugging system
 * Provides signal handling, stack traces, and crash forensics
 */

/*
 * crash_debug.c - Implementation of comprehensive crash debugging system
 * Provides signal handling, stack traces, and crash forensics
 *
 * Include dlfcn.h FIRST before _XOPEN_SOURCE restricts visibility.
 * On macOS, _XOPEN_SOURCE hides Dl_info/dladdr from dlfcn.h.
 * On Linux, _GNU_SOURCE exposes them.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <dlfcn.h>

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include "crash_debug.h"
#include "logger.h"
#include "resource_monitor.h"
#include "symbol_resolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <errno.h>

// Global state
volatile sig_atomic_t g_crash_in_progress = 0;
crash_config_t g_crash_config = {0};
static struct sigaction g_old_handlers[32];
static int g_handlers_installed = 0;

// Global application context (thread-safe with mutex)
static struct {
    char current_operation[256];
    char current_file[512];
    int current_file_line;
    char current_read_name[256];
    int current_batch;
    int current_read_index;
    int current_tree;
    char processing_stage[128];
    
    // Data corruption tracking
    char corrupted_file[512];
    int corrupted_line;
    char corruption_type[128];

    pthread_mutex_t context_mutex;

    // BWA bounds tracking
    int bwa_leaf_iter;
    int bwa_max_matches;
    int bwa_concordant_count;
    int bwa_discordant_count;
    int bwa_unique_trees;
    int bwa_dropped_matches;
} g_app_context = {0};

// Signal information mapping
static struct {
    int signal;
    const char* name;
    const char* description;
} signal_info[] = {
    {SIGSEGV, "SIGSEGV", "Segmentation fault"},
    {SIGABRT, "SIGABRT", "Abort signal"},
    {SIGBUS, "SIGBUS", "Bus error"},
    {SIGFPE, "SIGFPE", "Floating point exception"},
    {SIGILL, "SIGILL", "Illegal instruction"},
    {SIGPIPE, "SIGPIPE", "Broken pipe"},
    {SIGTERM, "SIGTERM", "Termination signal"},
    {SIGINT, "SIGINT", "Interrupt signal"},
    {0, NULL, NULL}
};

// Initialize crash debugging system
int crash_debug_init(const crash_config_t* config) {
    if (config) {
        memcpy(&g_crash_config, config, sizeof(crash_config_t));
    } else {
        // Default configuration
        g_crash_config.enable_stack_trace = 1;
        g_crash_config.enable_register_dump = 1;
        g_crash_config.enable_memory_dump = 1;
        g_crash_config.enable_core_dump = 1;
        g_crash_config.enable_crash_log = 1;
        strcpy(g_crash_config.crash_log_directory, "/tmp");
        strcpy(g_crash_config.crash_log_prefix, "tronko_crash");
        g_crash_config.user_crash_handler = NULL;
    }
    
    // Create crash log directory if it doesn't exist
    struct stat st = {0};
    if (stat(g_crash_config.crash_log_directory, &st) == -1) {
        if (mkdir(g_crash_config.crash_log_directory, 0755) != 0) {
            LOG_WARN("Failed to create crash log directory: %s", 
                     g_crash_config.crash_log_directory);
        }
    }
    
    // Enable core dumps if requested
    if (g_crash_config.enable_core_dump) {
        crash_debug_enable_core_dumps();
    }
    
    // Initialize application context mutex
    if (pthread_mutex_init(&g_app_context.context_mutex, NULL) != 0) {
        LOG_WARN("Failed to initialize context mutex");
    }
    
    // Install signal handlers
    crash_setup_signal_handlers();
    
    LOG_INFO("Crash debugging system initialized");
    return 0;
}

// Cleanup crash debugging system
void crash_debug_cleanup(void) {
    if (g_handlers_installed) {
        crash_restore_default_handlers();
    }
    
    // Cleanup context mutex
    pthread_mutex_destroy(&g_app_context.context_mutex);
    
    LOG_DEBUG("Crash debugging system cleaned up");
}

// Enable core dumps
void crash_debug_enable_core_dumps(void) {
    struct rlimit core_limit;
    core_limit.rlim_cur = RLIM_INFINITY;
    core_limit.rlim_max = RLIM_INFINITY;
    
    if (setrlimit(RLIMIT_CORE, &core_limit) != 0) {
        LOG_WARN("Failed to enable core dumps: %s", strerror(errno));
    } else {
        LOG_DEBUG("Core dumps enabled");
    }
}

// Set up signal handlers for crash detection
void crash_setup_signal_handlers(void) {
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    
    sa.sa_sigaction = crash_signal_handler;
    sa.sa_flags = SA_SIGINFO | SA_RESTART | SA_RESETHAND;
    sigemptyset(&sa.sa_mask);
    
    // Install handlers for critical signals
    int signals[] = {SIGSEGV, SIGABRT, SIGBUS, SIGFPE, SIGILL, 0};
    
    for (int i = 0; signals[i] != 0; i++) {
        if (sigaction(signals[i], &sa, &g_old_handlers[signals[i]]) != 0) {
            LOG_ERROR("Failed to install handler for signal %d: %s", 
                      signals[i], strerror(errno));
        } else {
            LOG_DEBUG("Installed crash handler for signal %s", 
                      crash_signal_name(signals[i]));
        }
    }
    
    g_handlers_installed = 1;
}

// Restore default signal handlers
void crash_restore_default_handlers(void) {
    int signals[] = {SIGSEGV, SIGABRT, SIGBUS, SIGFPE, SIGILL, 0};
    
    for (int i = 0; signals[i] != 0; i++) {
        sigaction(signals[i], &g_old_handlers[signals[i]], NULL);
    }
    
    g_handlers_installed = 0;
    LOG_DEBUG("Restored default signal handlers");
}

// Main crash signal handler
void crash_signal_handler(int sig, siginfo_t* info, void* context) {
    // Prevent recursive crashes
    if (g_crash_in_progress) {
        _exit(1);
    }
    g_crash_in_progress = 1;
    
    // Prepare crash information
    crash_info_t crash_info;
    memset(&crash_info, 0, sizeof(crash_info));
    
    crash_info.signal_number = sig;
    crash_info.signal_code = info ? info->si_code : 0;
    crash_info.fault_address = info ? info->si_addr : NULL;
    crash_info.process_id = getpid();
    crash_info.thread_id = getpid(); // pthread_self() for threads
    crash_info.crash_time = time(NULL);
    
    // Get program information
    crash_get_program_info(&crash_info);
    
    // Capture application context (thread-safe)
    pthread_mutex_lock(&g_app_context.context_mutex);
    crash_safe_string_copy(crash_info.current_operation, g_app_context.current_operation, sizeof(crash_info.current_operation));
    crash_safe_string_copy(crash_info.current_file, g_app_context.current_file, sizeof(crash_info.current_file));
    crash_safe_string_copy(crash_info.current_read_name, g_app_context.current_read_name, sizeof(crash_info.current_read_name));
    crash_safe_string_copy(crash_info.processing_stage, g_app_context.processing_stage, sizeof(crash_info.processing_stage));
    crash_info.current_file_line = g_app_context.current_file_line;
    crash_info.current_batch = g_app_context.current_batch;
    crash_info.current_read_index = g_app_context.current_read_index;
    crash_info.current_tree = g_app_context.current_tree;
    
    // Capture corruption tracking
    crash_safe_string_copy(crash_info.corrupted_file, g_app_context.corrupted_file, sizeof(crash_info.corrupted_file));
    crash_info.corrupted_line = g_app_context.corrupted_line;
    crash_safe_string_copy(crash_info.corruption_type, g_app_context.corruption_type, sizeof(crash_info.corruption_type));

    // Capture BWA context
    crash_info.bwa_leaf_iter = g_app_context.bwa_leaf_iter;
    crash_info.bwa_max_matches = g_app_context.bwa_max_matches;
    crash_info.bwa_concordant_count = g_app_context.bwa_concordant_count;
    crash_info.bwa_discordant_count = g_app_context.bwa_discordant_count;
    crash_info.bwa_unique_trees = g_app_context.bwa_unique_trees;
    crash_info.bwa_dropped_matches = g_app_context.bwa_dropped_matches;

    pthread_mutex_unlock(&g_app_context.context_mutex);
    
    // Create crash message
    snprintf(crash_info.crash_message, sizeof(crash_info.crash_message),
             "CRASH: %s (%s) at address %p in process %d",
             crash_signal_name(sig),
             crash_signal_description(sig, crash_info.signal_code),
             crash_info.fault_address,
             crash_info.process_id);
    
    // Capture stack trace
    if (g_crash_config.enable_stack_trace) {
        crash_info.stack_depth = crash_capture_stack_trace(
            crash_info.stack_trace, MAX_STACK_TRACE_DEPTH);
        crash_info.symbol_names = crash_resolve_symbols(
            crash_info.stack_trace, crash_info.stack_depth);
    }
    
    // Capture register state
    if (g_crash_config.enable_register_dump && context) {
        crash_capture_register_state((ucontext_t*)context, &crash_info);
    }
    
    // Capture memory state
    if (g_crash_config.enable_memory_dump) {
        crash_capture_memory_state(&crash_info);
    }
    
    // Generate crash report
    if (g_crash_config.enable_crash_log) {
        crash_write_crash_log(&crash_info);
    }
    
    // Call user handler if set
    if (g_crash_config.user_crash_handler) {
        g_crash_config.user_crash_handler(&crash_info);
    }
    
    // Log to stderr for immediate visibility
    fprintf(stderr, "\n=== CRASH DETECTED ===\n");
    fprintf(stderr, "%s\n", crash_info.crash_message);
    if (crash_info.stack_depth > 0) {
        fprintf(stderr, "\nEnhanced stack trace:\n");
        print_enhanced_stack_trace(crash_info.stack_trace, crash_info.stack_depth, stderr);
    }
    fprintf(stderr, "=== END CRASH REPORT ===\n");
    
    // Clean up and exit
    if (crash_info.symbol_names) {
        free(crash_info.symbol_names);
    }
    
    // Re-raise signal to generate core dump if enabled
    if (g_crash_config.enable_core_dump) {
        signal(crash_info.signal_number, SIG_DFL);
        raise(crash_info.signal_number);
    }
    
    _exit(1);
}

// Capture stack trace
int crash_capture_stack_trace(void** trace, int max_depth) {
    return backtrace(trace, max_depth);
}

// Resolve symbol names from stack trace
char** crash_resolve_symbols(void* const* trace, int depth) {
    return backtrace_symbols(trace, depth);
}

// Print stack trace to file
void crash_print_stack_trace(const void* const* trace, int depth, FILE* output) {
    char** symbols = backtrace_symbols((void* const*)trace, depth);
    if (!symbols) {
        fprintf(output, "Failed to resolve stack trace symbols\n");
        return;
    }
    
    for (int i = 0; i < depth; i++) {
        fprintf(output, "#%d  %s\n", i, symbols[i]);
        
        // Try to get more detailed information using dladdr
        Dl_info info;
        if (dladdr(trace[i], &info)) {
            if (info.dli_sname) {
                fprintf(output, "    Symbol: %s\n", info.dli_sname);
            }
            if (info.dli_fname) {
                fprintf(output, "    Object: %s\n", info.dli_fname);
            }
        }
    }
    
    free(symbols);
}

// Log stack trace using the logger
void crash_log_stack_trace(const void* const* trace, int depth) {
    LOG_ERROR("Stack trace (%d frames):", depth);
    
    char** symbols = backtrace_symbols((void* const*)trace, depth);
    if (!symbols) {
        LOG_ERROR("Failed to resolve stack trace symbols");
        return;
    }
    
    for (int i = 0; i < depth; i++) {
        LOG_ERROR("#%d  %s", i, symbols[i]);
    }
    
    free(symbols);
}

// Capture register state (platform specific)
void crash_capture_register_state(const ucontext_t* context, crash_info_t* crash_info) {
    if (!context) return;
    
#ifdef __APPLE__
    crash_info->machine_context = context->uc_mcontext;
#endif
#ifdef __linux__
    crash_info->machine_context = context->uc_mcontext;
#endif
}

// Capture memory state
void crash_capture_memory_state(crash_info_t* crash_info) {
    resource_stats_t stats;
    if (get_resource_stats(&stats) == 0) {
        crash_info->memory_rss_kb = stats.memory_rss_kb;
        crash_info->memory_vm_kb = stats.memory_vm_size_kb;
        crash_info->cpu_percent = stats.cpu_percent;
    }
}

// Print register state
void crash_print_register_state(const crash_info_t* crash_info, FILE* output) {
    fprintf(output, "\nRegister state:\n");
    
#ifdef __APPLE__
    #ifdef __x86_64__
    const _STRUCT_X86_THREAD_STATE64* regs = &crash_info->machine_context->__ss;
    fprintf(output, "RAX: 0x%016llx  RBX: 0x%016llx\n", regs->__rax, regs->__rbx);
    fprintf(output, "RCX: 0x%016llx  RDX: 0x%016llx\n", regs->__rcx, regs->__rdx);
    fprintf(output, "RSI: 0x%016llx  RDI: 0x%016llx\n", regs->__rsi, regs->__rdi);
    fprintf(output, "RBP: 0x%016llx  RSP: 0x%016llx\n", regs->__rbp, regs->__rsp);
    fprintf(output, "RIP: 0x%016llx  RFLAGS: 0x%016llx\n", regs->__rip, regs->__rflags);
    #endif
    #ifdef __arm64__
    const _STRUCT_ARM_THREAD_STATE64* regs = &crash_info->machine_context->__ss;
    fprintf(output, "X0:  0x%016llx  X1:  0x%016llx\n", regs->__x[0], regs->__x[1]);
    fprintf(output, "X2:  0x%016llx  X3:  0x%016llx\n", regs->__x[2], regs->__x[3]);
    fprintf(output, "LR:  0x%016llx  SP:  0x%016llx\n", regs->__lr, regs->__sp);
    fprintf(output, "PC:  0x%016llx  CPSR: 0x%08x\n", regs->__pc, regs->__cpsr);
    #endif
#endif

#ifdef __linux__
    #ifdef __x86_64__
    fprintf(output, "RAX: 0x%016lx  RBX: 0x%016lx\n", 
            crash_info->machine_context.gregs[REG_RAX],
            crash_info->machine_context.gregs[REG_RBX]);
    fprintf(output, "RCX: 0x%016lx  RDX: 0x%016lx\n",
            crash_info->machine_context.gregs[REG_RCX],
            crash_info->machine_context.gregs[REG_RDX]);
    fprintf(output, "RSI: 0x%016lx  RDI: 0x%016lx\n",
            crash_info->machine_context.gregs[REG_RSI],
            crash_info->machine_context.gregs[REG_RDI]);
    fprintf(output, "RBP: 0x%016lx  RSP: 0x%016lx\n",
            crash_info->machine_context.gregs[REG_RBP],
            crash_info->machine_context.gregs[REG_RSP]);
    fprintf(output, "RIP: 0x%016lx\n",
            crash_info->machine_context.gregs[REG_RIP]);
    #endif
#endif
}

// Print memory state
void crash_print_memory_state(const crash_info_t* crash_info, FILE* output) {
    fprintf(output, "\nMemory state:\n");
    fprintf(output, "RSS: %zu KB\n", crash_info->memory_rss_kb);
    fprintf(output, "VM Size: %zu KB\n", crash_info->memory_vm_kb);
    fprintf(output, "CPU Usage: %.2f%%\n", crash_info->cpu_percent);
}

// Generate comprehensive crash report
void crash_generate_report(const crash_info_t* crash_info, const char* filename) {
    char report_filename[1024];
    
    if (filename) {
        strncpy(report_filename, filename, sizeof(report_filename) - 1);
    } else {
        snprintf(report_filename, sizeof(report_filename),
                 "%s/%s_%d_%ld.crash",
                 g_crash_config.crash_log_directory,
                 g_crash_config.crash_log_prefix,
                 crash_info->process_id,
                 crash_info->crash_time);
    }
    
    FILE* report_file = fopen(report_filename, "w");
    if (!report_file) {
        LOG_ERROR("Failed to create crash report file: %s", report_filename);
        return;
    }
    
    // Write crash report header
    fprintf(report_file, "=== TRONKO CRASH REPORT ===\n");
    fprintf(report_file, "Crash Time: %s", ctime(&crash_info->crash_time));
    fprintf(report_file, "Process ID: %d\n", crash_info->process_id);
    fprintf(report_file, "Program: %s\n", crash_info->program_name);
    fprintf(report_file, "Working Directory: %s\n", crash_info->working_directory);
    fprintf(report_file, "Command Line: %s\n", crash_info->command_line);
    fprintf(report_file, "\n%s\n\n", crash_info->crash_message);
    
    // Write application context
    fprintf(report_file, "Application Context:\n");
    if (strlen(crash_info->current_operation) > 0) {
        fprintf(report_file, "  Current Operation: %s\n", crash_info->current_operation);
    }
    if (strlen(crash_info->processing_stage) > 0) {
        fprintf(report_file, "  Processing Stage: %s\n", crash_info->processing_stage);
    }
    if (strlen(crash_info->current_file) > 0) {
        if (crash_info->current_file_line > 0) {
            fprintf(report_file, "  Current File: %s (line %d)\n", crash_info->current_file, crash_info->current_file_line);
        } else {
            fprintf(report_file, "  Current File: %s\n", crash_info->current_file);
        }
    }
    if (strlen(crash_info->current_read_name) > 0) {
        fprintf(report_file, "  Current Read: %s (batch %d, index %d)\n", 
                crash_info->current_read_name, crash_info->current_batch, crash_info->current_read_index);
    }
    if (crash_info->current_tree >= 0) {
        fprintf(report_file, "  Current Tree: %d\n", crash_info->current_tree);
    }
    
    // Write corruption information if available
    if (strlen(crash_info->corrupted_file) > 0) {
        fprintf(report_file, "  *** DATA CORRUPTION DETECTED ***\n");
        if (crash_info->corrupted_line > 0) {
            fprintf(report_file, "  Corrupted File: %s (line %d)\n", crash_info->corrupted_file, crash_info->corrupted_line);
        } else {
            fprintf(report_file, "  Corrupted File: %s\n", crash_info->corrupted_file);
        }
        if (strlen(crash_info->corruption_type) > 0) {
            fprintf(report_file, "  Corruption Type: %s\n", crash_info->corruption_type);
        }
    }

    // Write BWA bounds information if relevant
    if (crash_info->bwa_leaf_iter > 0 || crash_info->bwa_dropped_matches > 0) {
        fprintf(report_file, "\nBWA Bounds Context:\n");
        fprintf(report_file, "  leaf_iter: %d (max: %d)\n",
                crash_info->bwa_leaf_iter, crash_info->bwa_max_matches);
        fprintf(report_file, "  Concordant matches: %d\n", crash_info->bwa_concordant_count);
        fprintf(report_file, "  Discordant matches: %d\n", crash_info->bwa_discordant_count);
        fprintf(report_file, "  Unique trees: %d\n", crash_info->bwa_unique_trees);
        if (crash_info->bwa_dropped_matches > 0) {
            fprintf(report_file, "  *** MATCHES DROPPED: %d (bounds exceeded) ***\n",
                    crash_info->bwa_dropped_matches);
        }
    }

    fprintf(report_file, "\n");

    // Write stack trace
    if (crash_info->stack_depth > 0) {
        fprintf(report_file, "Stack Trace:\n");
        crash_print_stack_trace((const void* const*)crash_info->stack_trace, crash_info->stack_depth, report_file);
        fprintf(report_file, "\n");
    }
    
    // Write register state
    if (g_crash_config.enable_register_dump) {
        crash_print_register_state(crash_info, report_file);
        fprintf(report_file, "\n");
    }
    
    // Write memory state
    if (g_crash_config.enable_memory_dump) {
        crash_print_memory_state(crash_info, report_file);
        fprintf(report_file, "\n");
    }
    
    fprintf(report_file, "=== END CRASH REPORT ===\n");
    fclose(report_file);
    
    LOG_ERROR("Crash report written to: %s", report_filename);
}

// Write crash log using the logger
void crash_write_crash_log(const crash_info_t* crash_info) {
    LOG_ERROR("=== CRASH DETECTED ===");
    LOG_ERROR("Process: %d, Signal: %s", crash_info->process_id, 
              crash_signal_name(crash_info->signal_number));
    LOG_ERROR("Message: %s", crash_info->crash_message);
    
    if (crash_info->stack_depth > 0) {
        crash_log_stack_trace((const void* const*)crash_info->stack_trace, crash_info->stack_depth);
    }
    
    // Generate detailed report file
    crash_generate_report(crash_info, NULL);
}

// Get signal name
const char* crash_signal_name(int signal) {
    for (int i = 0; signal_info[i].name; i++) {
        if (signal_info[i].signal == signal) {
            return signal_info[i].name;
        }
    }
    return "UNKNOWN";
}

// Get signal description
const char* crash_signal_description(int signal, int code) {
    for (int i = 0; signal_info[i].description; i++) {
        if (signal_info[i].signal == signal) {
            return signal_info[i].description;
        }
    }
    return "Unknown signal";
}

// Get program information
void crash_get_program_info(crash_info_t* crash_info) {
    // Get program name
    if (readlink("/proc/self/exe", crash_info->program_name, 
                 sizeof(crash_info->program_name) - 1) == -1) {
        strcpy(crash_info->program_name, "unknown");
    }
    
    // Get working directory
    if (!getcwd(crash_info->working_directory, sizeof(crash_info->working_directory))) {
        strcpy(crash_info->working_directory, "unknown");
    }
    
    // Try to get command line (platform specific)
    strcpy(crash_info->command_line, "unknown");
}

// Check if pointer is valid (basic check)
int crash_is_valid_pointer(void* ptr) {
    if (!ptr) return 0;
    
    // Basic alignment check for most architectures
    if ((uintptr_t)ptr % sizeof(void*) != 0) return 0;
    
    // Try to read one byte (this is not foolproof)
    volatile char test;
    volatile char* test_ptr = (volatile char*)ptr;
    
    // This is a very basic check - in production you might want
    // to use more sophisticated methods
    return 1;
}

// Set user crash handler
void crash_debug_set_user_handler(void (*handler)(const crash_info_t* crash_info)) {
    g_crash_config.user_crash_handler = handler;
}

// Set crash log directory
void crash_debug_set_crash_log_dir(const char* directory) {
    if (directory) {
        crash_safe_string_copy(g_crash_config.crash_log_directory, 
                               directory, sizeof(g_crash_config.crash_log_directory));
    }
}

// Safe string copy
void crash_safe_string_copy(char* dest, const char* src, size_t dest_size) {
    if (!dest || !src || dest_size == 0) return;
    
    strncpy(dest, src, dest_size - 1);
    dest[dest_size - 1] = '\0';
}

// Application context tracking functions
void crash_set_context(const char* operation, const char* details) {
    if (!operation) return;
    
    pthread_mutex_lock(&g_app_context.context_mutex);
    crash_safe_string_copy(g_app_context.current_operation, operation, sizeof(g_app_context.current_operation));
    if (details) {
        // Append details to operation
        size_t op_len = strlen(g_app_context.current_operation);
        if (op_len < sizeof(g_app_context.current_operation) - 3) {
            strcat(g_app_context.current_operation, ": ");
            strncat(g_app_context.current_operation, details, 
                    sizeof(g_app_context.current_operation) - op_len - 3);
        }
    }
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_current_file(const char* filename) {
    if (!filename) return;
    
    pthread_mutex_lock(&g_app_context.context_mutex);
    crash_safe_string_copy(g_app_context.current_file, filename, sizeof(g_app_context.current_file));
    g_app_context.current_file_line = 0;  // Reset line number
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_current_file_line(const char* filename, int line_number) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    if (filename) {
        crash_safe_string_copy(g_app_context.current_file, filename, sizeof(g_app_context.current_file));
    }
    g_app_context.current_file_line = line_number;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_current_read(const char* read_name, int batch, int read_index) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    if (read_name) {
        crash_safe_string_copy(g_app_context.current_read_name, read_name, sizeof(g_app_context.current_read_name));
    }
    g_app_context.current_batch = batch;
    g_app_context.current_read_index = read_index;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_current_tree(int tree_number) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.current_tree = tree_number;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_processing_stage(const char* stage) {
    if (!stage) return;
    
    pthread_mutex_lock(&g_app_context.context_mutex);
    crash_safe_string_copy(g_app_context.processing_stage, stage, sizeof(g_app_context.processing_stage));
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_clear_context(void) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    memset(g_app_context.current_operation, 0, sizeof(g_app_context.current_operation));
    memset(g_app_context.current_file, 0, sizeof(g_app_context.current_file));
    memset(g_app_context.current_read_name, 0, sizeof(g_app_context.current_read_name));
    memset(g_app_context.processing_stage, 0, sizeof(g_app_context.processing_stage));
    g_app_context.current_file_line = -1;
    g_app_context.current_batch = -1;
    g_app_context.current_read_index = -1;
    g_app_context.current_tree = -1;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

// Data corruption tracking functions
void crash_flag_corruption(const char* filename, int line_number, const char* corruption_type) {
    if (!filename || !corruption_type) return;
    
    pthread_mutex_lock(&g_app_context.context_mutex);
    crash_safe_string_copy(g_app_context.corrupted_file, filename, sizeof(g_app_context.corrupted_file));
    g_app_context.corrupted_line = line_number;
    crash_safe_string_copy(g_app_context.corruption_type, corruption_type, sizeof(g_app_context.corruption_type));
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_clear_corruption_flags(void) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    memset(g_app_context.corrupted_file, 0, sizeof(g_app_context.corrupted_file));
    g_app_context.corrupted_line = -1;
    memset(g_app_context.corruption_type, 0, sizeof(g_app_context.corruption_type));
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

// BWA bounds tracking functions
void crash_set_bwa_context(int leaf_iter, int concordant_count, int discordant_count) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = leaf_iter;
    g_app_context.bwa_concordant_count = concordant_count;
    g_app_context.bwa_discordant_count = discordant_count;
    g_app_context.bwa_unique_trees = leaf_iter;  // leaf_iter tracks unique trees
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_set_bwa_bounds_violation(int leaf_iter, int max_matches, int dropped) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = leaf_iter;
    g_app_context.bwa_max_matches = max_matches;
    g_app_context.bwa_dropped_matches = dropped;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}

void crash_clear_bwa_context(void) {
    pthread_mutex_lock(&g_app_context.context_mutex);
    g_app_context.bwa_leaf_iter = 0;
    g_app_context.bwa_max_matches = 0;
    g_app_context.bwa_concordant_count = 0;
    g_app_context.bwa_discordant_count = 0;
    g_app_context.bwa_unique_trees = 0;
    g_app_context.bwa_dropped_matches = 0;
    pthread_mutex_unlock(&g_app_context.context_mutex);
}