/*
 * symbol_resolver.h - Enhanced symbol resolution for crash debugging
 * Provides addr2line integration and debug symbol management
 */

#ifndef _SYMBOL_RESOLVER_H_
#define _SYMBOL_RESOLVER_H_

#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

// Symbol information structure
typedef struct {
    void* address;
    char symbol_name[256];
    char source_file[512];
    char object_file[512];
    int line_number;
    int column_number;
    uintptr_t offset;
    int is_resolved;
} symbol_info_t;

// Symbol resolver configuration
typedef struct {
    int enable_addr2line;
    int enable_source_lookup;
    int cache_symbols;
    char debug_info_path[1024];
    char addr2line_path[256];
} symbol_resolver_config_t;

// Function declarations
int symbol_resolver_init(const symbol_resolver_config_t* config);
void symbol_resolver_cleanup(void);

// Symbol resolution functions
int resolve_symbol_info(void* address, symbol_info_t* info);
int resolve_stack_trace_symbols(void* const* addresses, int count, symbol_info_t* symbols);
char* format_symbol_info(const symbol_info_t* info);

// addr2line integration
int addr2line_resolve(void* address, symbol_info_t* info);
int addr2line_batch_resolve(void* const* addresses, int count, symbol_info_t* symbols);
char* addr2line_get_source_line(const char* filename, int line_number);

// Debug information management
int find_debug_info(const char* binary_path, char* debug_path, size_t debug_path_size);
int verify_debug_symbols(const char* binary_path);
void print_debug_info_status(void);

// Symbol caching
typedef struct symbol_cache_entry {
    void* address;
    symbol_info_t info;
    struct symbol_cache_entry* next;
} symbol_cache_entry_t;

int symbol_cache_lookup(void* address, symbol_info_t* info);
void symbol_cache_store(void* address, const symbol_info_t* info);
void symbol_cache_clear(void);

// Enhanced stack trace formatting
void print_enhanced_stack_trace(void* const* addresses, int count, FILE* output);
void log_enhanced_stack_trace(void* const* addresses, int count);
char* format_stack_trace_html(void* const* addresses, int count);

// Binary analysis helpers
int get_binary_load_address(const char* binary_path, void** load_address);
int get_symbol_from_binary(const char* binary_path, void* address, char* symbol_name, size_t name_size);

// Utility functions
const char* demangle_symbol_name(const char* mangled_name);
int is_debug_build(void);
void print_compilation_info(void);

// Global configuration
extern symbol_resolver_config_t g_symbol_config;

#endif /* _SYMBOL_RESOLVER_H_ */