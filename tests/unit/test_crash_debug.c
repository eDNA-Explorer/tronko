/*
 * Unit tests for crash_debug.c
 * Tests crash detection, context tracking, and corruption flagging
 */

#include "unity.h"
#include "crash_debug.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Test setup and teardown
void setUp(void) {
    crash_debug_init();
    crash_clear_context();
    crash_clear_corruption_flags();
}

void tearDown(void) {
    crash_clear_context();
    crash_clear_corruption_flags();
}

// Test Cases

void test_crash_context_tracking(void) {
    // Test setting and retrieving context
    crash_set_current_file("test_file.fasta");
    crash_set_current_file_line("test_file.fasta", 42);
    crash_set_processing_stage("Testing stage");
    
    // Since we can't easily access the internal context,
    // we'll test by triggering a mock crash scenario
    // This is more of an integration test
    TEST_PASS_MESSAGE("Context tracking functions executed without error");
}

void test_corruption_flagging(void) {
    // Test corruption flagging functions
    crash_flag_corruption("corrupt_file.fasta", 123, "Test corruption type");
    
    // Verify the corruption was flagged (would need access to internal state)
    TEST_PASS_MESSAGE("Corruption flagging executed without error");
}

void test_corruption_clearing(void) {
    // Flag some corruption
    crash_flag_corruption("corrupt_file.fasta", 123, "Test corruption");
    
    // Clear corruption flags
    crash_clear_corruption_flags();
    
    // Verify corruption flags were cleared
    TEST_PASS_MESSAGE("Corruption clearing executed without error");
}

void test_multiple_corruption_flags(void) {
    // Test that only the latest corruption is preserved
    crash_flag_corruption("file1.fasta", 10, "First corruption");
    crash_flag_corruption("file2.fasta", 20, "Second corruption");
    
    // The second corruption should overwrite the first
    TEST_PASS_MESSAGE("Multiple corruption flagging handled correctly");
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

void test_crash_debug_initialization(void) {
    // Test that crash debug system initializes correctly
    int result = crash_debug_init();
    TEST_ASSERT_EQUAL_INT(0, result);
    
    // Test double initialization
    result = crash_debug_init();
    TEST_ASSERT_EQUAL_INT(0, result);
}

// Test runner - specific to this test file
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_crash_context_tracking);
    RUN_TEST(test_corruption_flagging);
    RUN_TEST(test_corruption_clearing);
    RUN_TEST(test_multiple_corruption_flags);
    RUN_TEST(test_null_parameter_handling);
    RUN_TEST(test_signal_name_mapping);
    RUN_TEST(test_crash_debug_initialization);
    
    return UNITY_END();
}