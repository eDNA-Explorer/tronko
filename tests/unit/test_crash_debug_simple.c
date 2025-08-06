/*
 * Simple unit tests for crash_debug.c
 * Tests core crash debugging functionality in isolation
 */

#include "unity.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>

// Include just the crash debug header and minimal dependencies
#define _GNU_SOURCE
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>

// Declare only the functions we need to test
int crash_debug_init(void);
void crash_flag_corruption(const char* filename, int line_number, const char* corruption_type);
void crash_clear_corruption_flags(void);
void crash_set_current_file(const char* filename);
void crash_set_processing_stage(const char* stage);
const char* crash_signal_name(int signal);
const char* crash_signal_description(int signal, int code);

// Test setup and teardown
void setUp(void) {
    // Initialize crash debugging for tests
    crash_debug_init();
    crash_clear_corruption_flags();
}

void tearDown(void) {
    // Clean up after each test
    crash_clear_corruption_flags();
}

// Test Cases

void test_crash_debug_initialization(void) {
    // Test that crash debug system initializes correctly
    int result = crash_debug_init();
    TEST_ASSERT_EQUAL_INT(0, result);
    
    // Test double initialization
    result = crash_debug_init();
    TEST_ASSERT_EQUAL_INT(0, result);
}

void test_corruption_flagging(void) {
    // Test corruption flagging functions don't crash
    crash_flag_corruption("corrupt_file.fasta", 123, "Test corruption type");
    
    // Multiple corruption flags
    crash_flag_corruption("file1.fasta", 10, "First corruption");
    crash_flag_corruption("file2.fasta", 20, "Second corruption");
    
    TEST_PASS_MESSAGE("Corruption flagging executed without error");
}

void test_corruption_clearing(void) {
    // Flag some corruption
    crash_flag_corruption("corrupt_file.fasta", 123, "Test corruption");
    
    // Clear corruption flags
    crash_clear_corruption_flags();
    
    TEST_PASS_MESSAGE("Corruption clearing executed without error");
}

void test_context_setting(void) {
    // Test setting various context values
    crash_set_current_file("test_file.fasta");
    crash_set_processing_stage("Testing stage");
    
    TEST_PASS_MESSAGE("Context setting executed without error");
}

void test_null_parameter_handling(void) {
    // Test that functions handle NULL parameters gracefully
    crash_flag_corruption(NULL, 123, "Test corruption");
    crash_flag_corruption("file.fasta", 123, NULL);
    crash_set_current_file(NULL);
    crash_set_processing_stage(NULL);
    
    TEST_PASS_MESSAGE("NULL parameter handling completed without crash");
}

void test_signal_name_mapping(void) {
    // Test signal name and description functions
    const char* sig_name = crash_signal_name(SIGSEGV);
    const char* sig_desc = crash_signal_description(SIGSEGV, 0);
    
    TEST_ASSERT_NOT_NULL(sig_name);
    TEST_ASSERT_NOT_NULL(sig_desc);
    TEST_ASSERT_EQUAL_STRING("SIGSEGV", sig_name);
}

void test_signal_handling_setup(void) {
    // Test that signal handlers can be installed without error
    // This is mostly testing that the initialization doesn't crash
    int result = crash_debug_init();
    TEST_ASSERT_EQUAL_INT(0, result);
    
    // Test some signal name lookups
    TEST_ASSERT_EQUAL_STRING("SIGABRT", crash_signal_name(SIGABRT));
    TEST_ASSERT_EQUAL_STRING("SIGBUS", crash_signal_name(SIGBUS));
    TEST_ASSERT_EQUAL_STRING("SIGFPE", crash_signal_name(SIGFPE));
    TEST_ASSERT_EQUAL_STRING("SIGILL", crash_signal_name(SIGILL));
}

// Test runner - specific to this test file
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_crash_debug_initialization);
    RUN_TEST(test_corruption_flagging);
    RUN_TEST(test_corruption_clearing);
    RUN_TEST(test_context_setting);
    RUN_TEST(test_null_parameter_handling);
    RUN_TEST(test_signal_name_mapping);
    RUN_TEST(test_signal_handling_setup);
    
    return UNITY_END();
}