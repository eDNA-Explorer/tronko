/*
 * symbol_resolver.c - Implementation of enhanced symbol resolution
 * Provides addr2line integration and debug symbol management
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "symbol_resolver.h"
#include "logger.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include <sys/stat.h>
#include <errno.h>

// Global configuration and state
symbol_resolver_config_t g_symbol_config = {0};
static symbol_cache_entry_t* g_symbol_cache = NULL;
static int g_cache_size = 0;
static const int MAX_CACHE_SIZE = 1000;

// Initialize symbol resolver
int symbol_resolver_init(const symbol_resolver_config_t* config) {
    if (config) {
        memcpy(&g_symbol_config, config, sizeof(symbol_resolver_config_t));
    } else {
        // Default configuration
        g_symbol_config.enable_addr2line = 1;
        g_symbol_config.enable_source_lookup = 1;
        g_symbol_config.cache_symbols = 1;
        strcpy(g_symbol_config.debug_info_path, "/usr/lib/debug");
        strcpy(g_symbol_config.addr2line_path, "addr2line");
    }
    
    // Verify addr2line availability
    if (g_symbol_config.enable_addr2line) {
        char cmd[512];
        snprintf(cmd, sizeof(cmd), "which %s > /dev/null 2>&1", g_symbol_config.addr2line_path);
        if (system(cmd) != 0) {
            LOG_WARN("addr2line not found, disabling enhanced symbol resolution");
            g_symbol_config.enable_addr2line = 0;
        } else {
            LOG_DEBUG("addr2line found: %s", g_symbol_config.addr2line_path);
        }
    }
    
    LOG_INFO("Symbol resolver initialized (addr2line: %s, cache: %s)",
             g_symbol_config.enable_addr2line ? "enabled" : "disabled",
             g_symbol_config.cache_symbols ? "enabled" : "disabled");
    
    return 0;
}

// Cleanup symbol resolver
void symbol_resolver_cleanup(void) {
    if (g_symbol_config.cache_symbols) {
        symbol_cache_clear();
    }
    LOG_DEBUG("Symbol resolver cleaned up");
}

// Resolve symbol information for a single address
int resolve_symbol_info(void* address, symbol_info_t* info) {
    if (!info) return -1;
    
    memset(info, 0, sizeof(symbol_info_t));
    info->address = address;
    
    // Check cache first
    if (g_symbol_config.cache_symbols && symbol_cache_lookup(address, info)) {
        return 0;
    }
    
    // Use dladdr for basic symbol information
    Dl_info dl_info;
    if (dladdr(address, &dl_info)) {
        if (dl_info.dli_sname) {
            strncpy(info->symbol_name, dl_info.dli_sname, sizeof(info->symbol_name) - 1);
        }
        if (dl_info.dli_fname) {
            strncpy(info->object_file, dl_info.dli_fname, sizeof(info->object_file) - 1);
        }
        
        // Calculate offset
        if (dl_info.dli_saddr) {
            info->offset = (uintptr_t)address - (uintptr_t)dl_info.dli_saddr;
        }
        
        info->is_resolved = 1;
    }
    
    // Try addr2line for source file and line number
    if (g_symbol_config.enable_addr2line && strlen(info->object_file) > 0) {
        addr2line_resolve(address, info);
    }
    
    // Demangle C++ symbols if possible
    if (strlen(info->symbol_name) > 0) {
        const char* demangled = demangle_symbol_name(info->symbol_name);
        if (demangled && demangled != info->symbol_name) {
            strncpy(info->symbol_name, demangled, sizeof(info->symbol_name) - 1);
        }
    }
    
    // Cache the result
    if (g_symbol_config.cache_symbols && info->is_resolved) {
        symbol_cache_store(address, info);
    }
    
    return info->is_resolved ? 0 : -1;
}

// Resolve symbols for entire stack trace
int resolve_stack_trace_symbols(void* const* addresses, int count, symbol_info_t* symbols) {
    if (!addresses || !symbols || count <= 0) return -1;
    
    int resolved_count = 0;
    
    for (int i = 0; i < count; i++) {
        if (resolve_symbol_info(addresses[i], &symbols[i]) == 0) {
            resolved_count++;
        }
    }
    
    return resolved_count;
}

// Format symbol information as string
char* format_symbol_info(const symbol_info_t* info) {
    if (!info) return NULL;
    
    static char buffer[2048];
    char address_str[32];
    
    snprintf(address_str, sizeof(address_str), "%p", info->address);
    
    if (info->is_resolved) {
        if (strlen(info->source_file) > 0 && info->line_number > 0) {
            snprintf(buffer, sizeof(buffer),
                     "%s: %s at %s:%d",
                     address_str,
                     strlen(info->symbol_name) > 0 ? info->symbol_name : "??",
                     info->source_file,
                     info->line_number);
        } else if (strlen(info->symbol_name) > 0) {
            snprintf(buffer, sizeof(buffer),
                     "%s: %s+0x%lx (%s)",
                     address_str,
                     info->symbol_name,
                     info->offset,
                     strlen(info->object_file) > 0 ? info->object_file : "??");
        } else {
            snprintf(buffer, sizeof(buffer),
                     "%s: ?? (%s)",
                     address_str,
                     strlen(info->object_file) > 0 ? info->object_file : "??");
        }
    } else {
        snprintf(buffer, sizeof(buffer), "%s: ??", address_str);
    }
    
    return buffer;
}

// Use addr2line to resolve address to source file and line
int addr2line_resolve(void* address, symbol_info_t* info) {
    if (!info || strlen(info->object_file) == 0) return -1;
    
    char cmd[1024];
    char result[512];
    FILE* pipe;
    
    // Calculate relative address for addr2line
    uintptr_t relative_addr = (uintptr_t)address;
    
    // Build addr2line command
    snprintf(cmd, sizeof(cmd), "%s -C -f -e \"%s\" 0x%lx 2>/dev/null",
             g_symbol_config.addr2line_path, info->object_file, relative_addr);
    
    pipe = popen(cmd, "r");
    if (!pipe) {
        LOG_DEBUG("Failed to run addr2line: %s", strerror(errno));
        return -1;
    }
    
    // Read function name (first line)
    if (fgets(result, sizeof(result), pipe)) {
        result[strcspn(result, "\n")] = 0; // Remove newline
        if (strcmp(result, "??") != 0 && strlen(info->symbol_name) == 0) {
            strncpy(info->symbol_name, result, sizeof(info->symbol_name) - 1);
        }
    }
    
    // Read source file and line (second line)
    if (fgets(result, sizeof(result), pipe)) {
        result[strcspn(result, "\n")] = 0; // Remove newline
        
        char* colon = strrchr(result, ':');
        if (colon && strcmp(result, "??:0") != 0) {
            *colon = '\0';
            strncpy(info->source_file, result, sizeof(info->source_file) - 1);
            info->line_number = atoi(colon + 1);
        }
    }
    
    pclose(pipe);
    return 0;
}

// Batch resolve multiple addresses using addr2line
int addr2line_batch_resolve(void* const* addresses, int count, symbol_info_t* symbols) {
    if (!addresses || !symbols || count <= 0) return -1;
    
    // Group addresses by object file for efficient batch processing
    for (int i = 0; i < count; i++) {
        if (symbols[i].is_resolved && strlen(symbols[i].object_file) > 0) {
            addr2line_resolve(addresses[i], &symbols[i]);
        }
    }
    
    return 0;
}

// Get source line from file
char* addr2line_get_source_line(const char* filename, int line_number) {
    if (!filename || line_number <= 0) return NULL;
    
    static char line_buffer[512];
    FILE* file = fopen(filename, "r");
    if (!file) return NULL;
    
    int current_line = 1;
    while (fgets(line_buffer, sizeof(line_buffer), file)) {
        if (current_line == line_number) {
            fclose(file);
            // Remove trailing newline
            line_buffer[strcspn(line_buffer, "\n")] = 0;
            return line_buffer;
        }
        current_line++;
    }
    
    fclose(file);
    return NULL;
}

// Symbol cache implementation
int symbol_cache_lookup(void* address, symbol_info_t* info) {
    symbol_cache_entry_t* entry = g_symbol_cache;
    
    while (entry) {
        if (entry->address == address) {
            memcpy(info, &entry->info, sizeof(symbol_info_t));
            return 1;
        }
        entry = entry->next;
    }
    
    return 0;
}

void symbol_cache_store(void* address, const symbol_info_t* info) {
    if (g_cache_size >= MAX_CACHE_SIZE) {
        // Simple cache eviction - remove oldest entry
        symbol_cache_entry_t* old = g_symbol_cache;
        if (old) {
            g_symbol_cache = old->next;
            free(old);
            g_cache_size--;
        }
    }
    
    symbol_cache_entry_t* entry = malloc(sizeof(symbol_cache_entry_t));
    if (!entry) return;
    
    entry->address = address;
    memcpy(&entry->info, info, sizeof(symbol_info_t));
    entry->next = g_symbol_cache;
    g_symbol_cache = entry;
    g_cache_size++;
}

void symbol_cache_clear(void) {
    symbol_cache_entry_t* entry = g_symbol_cache;
    
    while (entry) {
        symbol_cache_entry_t* next = entry->next;
        free(entry);
        entry = next;
    }
    
    g_symbol_cache = NULL;
    g_cache_size = 0;
    LOG_DEBUG("Symbol cache cleared");
}

// Print enhanced stack trace
void print_enhanced_stack_trace(void* const* addresses, int count, FILE* output) {
    if (!addresses || !output || count <= 0) return;
    
    symbol_info_t* symbols = malloc(count * sizeof(symbol_info_t));
    if (!symbols) {
        fprintf(output, "Failed to allocate memory for symbol resolution\n");
        return;
    }
    
    int resolved = resolve_stack_trace_symbols(addresses, count, symbols);
    
    fprintf(output, "Stack trace (%d frames, %d resolved):\n", count, resolved);
    
    for (int i = 0; i < count; i++) {
        const char* formatted = format_symbol_info(&symbols[i]);
        fprintf(output, "#%-2d %s\n", i, formatted);
        
        // Show source code line if available
        if (g_symbol_config.enable_source_lookup && 
            strlen(symbols[i].source_file) > 0 && symbols[i].line_number > 0) {
            const char* source_line = addr2line_get_source_line(
                symbols[i].source_file, symbols[i].line_number);
            if (source_line) {
                fprintf(output, "     %s\n", source_line);
            }
        }
    }
    
    free(symbols);
}

// Log enhanced stack trace
void log_enhanced_stack_trace(void* const* addresses, int count) {
    if (!addresses || count <= 0) return;
    
    symbol_info_t* symbols = malloc(count * sizeof(symbol_info_t));
    if (!symbols) {
        LOG_ERROR("Failed to allocate memory for symbol resolution");
        return;
    }
    
    int resolved = resolve_stack_trace_symbols(addresses, count, symbols);
    LOG_ERROR("Enhanced stack trace (%d frames, %d resolved):", count, resolved);
    
    for (int i = 0; i < count; i++) {
        const char* formatted = format_symbol_info(&symbols[i]);
        LOG_ERROR("#%-2d %s", i, formatted);
    }
    
    free(symbols);
}

// Simple C++ name demangling (basic implementation)
const char* demangle_symbol_name(const char* mangled_name) {
    // This is a simplified version - in practice you might want to use
    // a proper demangling library like libbfd or call c++filt
    if (!mangled_name || mangled_name[0] != '_' || mangled_name[1] != 'Z') {
        return mangled_name;
    }
    
    // For now, just return the mangled name
    // TODO: Implement proper C++ demangling or use external tools
    return mangled_name;
}

// Check if this is a debug build
int is_debug_build(void) {
#ifdef DEBUG
    return 1;
#else
    return 0;
#endif
}

// Print compilation information
void print_compilation_info(void) {
    printf("Compilation info:\n");
    printf("  Debug build: %s\n", is_debug_build() ? "Yes" : "No");
    printf("  Compiler: %s\n", 
#ifdef __clang__
           "Clang"
#elif defined(__GNUC__)
           "GCC"
#else
           "Unknown"
#endif
    );
    printf("  Architecture: %s\n",
#ifdef __x86_64__
           "x86_64"
#elif defined(__aarch64__)
           "ARM64"
#elif defined(__arm__)
           "ARM"
#else
           "Unknown"
#endif
    );
}

// Find debug information for binary
int find_debug_info(const char* binary_path, char* debug_path, size_t debug_path_size) {
    if (!binary_path || !debug_path) return -1;
    
    struct stat st;
    
    // Try standard debug info locations
    const char* debug_locations[] = {
        "%s.debug",                    // Same directory with .debug suffix
        "%s/.debug/%s",               // .debug subdirectory
        "/usr/lib/debug%s",           // System debug directory
        "/usr/lib/debug%s.debug",     // System debug with suffix
        NULL
    };
    
    for (int i = 0; debug_locations[i]; i++) {
        snprintf(debug_path, debug_path_size, debug_locations[i], 
                 binary_path, strrchr(binary_path, '/') ? strrchr(binary_path, '/') + 1 : binary_path);
        
        if (stat(debug_path, &st) == 0) {
            return 0;
        }
    }
    
    return -1;
}

// Verify debug symbols are available
int verify_debug_symbols(const char* binary_path) {
    if (!binary_path) return 0;
    
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "file \"%s\" | grep -q 'not stripped'", binary_path);
    
    return system(cmd) == 0;
}