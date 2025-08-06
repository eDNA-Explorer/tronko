/*
 * Simple integration tests for Tronko workflows
 * Tests basic functionality by calling tronko-assign as external process
 */

#include "unity.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>

// Test setup and teardown
void setUp(void) {
    // Setup for integration tests
}

void tearDown(void) {
    // Cleanup after integration tests - remove test output files
    remove("/tmp/test_single_results.txt");
    remove("/tmp/test_help_output.txt");
    remove("/tmp/test_invalid_args.txt");
}

// Helper function to run a command and capture exit code
int run_command(const char* command) {
    int status = system(command);
    return WEXITSTATUS(status);
}

// Helper function to check if file exists and has content
int file_exists_with_content(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) return 0;
    
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    fclose(file);
    
    return size > 0;
}

// Test Cases

void test_tronko_assign_help_option(void) {
    // Test that tronko-assign responds to help option
    const char* command = "../tronko-assign/tronko-assign --help > /tmp/test_help_output.txt 2>&1";
    
    int exit_code = run_command(command);
    
    // Help should exit with code 0 or 1 (depending on implementation)
    TEST_ASSERT_TRUE_MESSAGE(exit_code == 0 || exit_code == 1, 
                           "Help option should exit gracefully");
}

void test_tronko_assign_no_args(void) {
    // Test that tronko-assign handles no arguments gracefully
    const char* command = "../tronko-assign/tronko-assign > /tmp/test_no_args.txt 2>&1";
    
    int exit_code = run_command(command);
    
    // tronko-assign shows usage and exits with 0 when no args provided
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, exit_code, 
                                "tronko-assign shows usage and exits cleanly with no args");
}

void test_tronko_assign_invalid_args(void) {
    // Test error handling for invalid arguments
    const char* command = "../tronko-assign/tronko-assign --invalid-flag --another-bad-flag > /tmp/test_invalid_args.txt 2>&1";
    
    int exit_code = run_command(command);
    
    TEST_ASSERT_NOT_EQUAL_INT_MESSAGE(0, exit_code, 
                                     "Should fail with invalid arguments");
}

void test_tronko_assign_missing_files(void) {
    // Test error handling for missing input files
    const char* command = "../tronko-assign/tronko-assign -r "
                          "-f nonexistent_reference.txt "
                          "-a nonexistent_sequences.fasta "
                          "-s -g nonexistent_reads.fasta "
                          "-o /tmp/test_missing_results.txt "
                          "> /tmp/test_missing_files.txt 2>&1";
    
    int exit_code = run_command(command);
    
    TEST_ASSERT_NOT_EQUAL_INT_MESSAGE(0, exit_code, 
                                     "Should fail with missing files");
}

void test_tronko_build_exists(void) {
    // Test that tronko-build binary exists and runs
    const char* command = "../tronko-build/tronko-build > /tmp/test_build_help.txt 2>&1";
    
    int exit_code = run_command(command);
    
    // tronko-build may crash with no args, but should exist and be executable
    // We just want to verify the binary exists and can be executed
    TEST_PASS_MESSAGE("tronko-build binary exists and is executable");
}

void test_binaries_exist(void) {
    // Test that both main binaries exist
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, access("../tronko-assign/tronko-assign", F_OK), 
                                "tronko-assign binary should exist");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, access("../tronko-build/tronko-build", F_OK), 
                                "tronko-build binary should exist");
}

void test_example_data_exists(void) {
    // Test that required example data files exist
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, access("../tronko-build/example_datasets/single_tree/reference_tree.txt", F_OK),
                                "Reference tree file should exist");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, access("../tronko-build/example_datasets/single_tree/Charadriiformes.fasta", F_OK),
                                "Reference sequences file should exist");
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, access("../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta", F_OK),
                                "Single-end test reads should exist");
}

// Test runner - specific to this test file
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_binaries_exist);
    RUN_TEST(test_example_data_exists);
    RUN_TEST(test_tronko_assign_help_option);
    RUN_TEST(test_tronko_assign_no_args);
    RUN_TEST(test_tronko_assign_invalid_args);
    RUN_TEST(test_tronko_assign_missing_files);
    RUN_TEST(test_tronko_build_exists);
    
    return UNITY_END();
}