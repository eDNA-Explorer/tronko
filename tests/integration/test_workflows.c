/*
 * Integration tests for Tronko workflows
 * Tests end-to-end functionality with real data
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
    // Cleanup after integration tests
}

// Helper function to run tronko-assign with specific parameters
int run_tronko_assign(const char* args, char* output_buffer, size_t buffer_size) {
    char command[1024];
    snprintf(command, sizeof(command), "../tronko-assign/tronko-assign %s 2>&1", args);
    
    FILE* pipe = popen(command, "r");
    if (!pipe) return -1;
    
    size_t bytes_read = 0;
    if (output_buffer && buffer_size > 0) {
        bytes_read = fread(output_buffer, 1, buffer_size - 1, pipe);
        output_buffer[bytes_read] = '\0';
    }
    
    int status = pclose(pipe);
    return WEXITSTATUS(status);
}

// Test Cases

void test_single_end_workflow_success(void) {
    // Test successful single-end read processing
    const char* args = "-r "
                      "-f ../tronko-build/example_datasets/single_tree/reference_tree.txt "
                      "-a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta "
                      "-s "
                      "-g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta "
                      "-o /tmp/test_single_results.txt "
                      "-w";
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, exit_code, "Single-end workflow should succeed");
    
    // Check that output file was created
    FILE* results = fopen("/tmp/test_single_results.txt", "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(results, "Results file should be created");
    if (results) fclose(results);
}

void test_paired_end_workflow_success(void) {
    // Test successful paired-end read processing
    const char* args = "-r "
                      "-f ../tronko-build/example_datasets/single_tree/reference_tree.txt "
                      "-a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta "
                      "-p "
                      "-1 ../example_datasets/single_tree/missingreads_pairedend_150bp_2error_read1.fasta "
                      "-2 ../example_datasets/single_tree/missingreads_pairedend_150bp_2error_read2.fasta "
                      "-o /tmp/test_paired_results.txt "
                      "-w";
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, exit_code, "Paired-end workflow should succeed");
    
    // Check that output file was created
    FILE* results = fopen("/tmp/test_paired_results.txt", "r");
    TEST_ASSERT_NOT_NULL_MESSAGE(results, "Results file should be created");
    if (results) fclose(results);
}

void test_corrupted_data_handling(void) {
    // Test that corrupted data is handled gracefully (may crash, but should be detected)
    const char* args = "-r "
                      "-f ../tronko-build/example_datasets/single_tree/reference_tree.txt "
                      "-a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta "
                      "-p "
                      "-1 ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta "
                      "-2 ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta "
                      "-o /tmp/test_corrupt_results.txt "
                      "-w -V 1";  // Use minimal logging to reduce noise
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    // This test may fail due to corruption - we're testing that it fails gracefully
    if (exit_code != 0) {
        // Check if crash debugging caught the issue
        TEST_ASSERT_TRUE_MESSAGE(strstr(output, "CRASH") != NULL || 
                               strstr(output, "corruption") != NULL ||
                               strstr(output, "error") != NULL,
                               "Should detect corruption or provide error information");
    }
}

void test_missing_file_error_handling(void) {
    // Test error handling for missing input files
    const char* args = "-r "
                      "-f nonexistent_reference.txt "
                      "-a nonexistent_sequences.fasta "
                      "-s "
                      "-g nonexistent_reads.fasta "
                      "-o /tmp/test_missing_results.txt";
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    TEST_ASSERT_NOT_EQUAL_INT_MESSAGE(0, exit_code, "Should fail with missing files");
}

void test_invalid_arguments_handling(void) {
    // Test error handling for invalid command line arguments
    const char* args = "--invalid-flag --another-invalid-flag";
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    TEST_ASSERT_NOT_EQUAL_INT_MESSAGE(0, exit_code, "Should fail with invalid arguments");
}

void test_memory_usage_monitoring(void) {
    // Test with resource monitoring enabled
    const char* args = "-r "
                      "-f ../tronko-build/example_datasets/single_tree/reference_tree.txt "
                      "-a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta "
                      "-s "
                      "-g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta "
                      "-o /tmp/test_memory_results.txt "
                      "-w -R -V 2";  // Enable resource monitoring and verbose logging
    
    char output[4096];
    int exit_code = run_tronko_assign(args, output, sizeof(output));
    
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, exit_code, "Workflow with monitoring should succeed");
    
    // Check that resource information is present in output
    TEST_ASSERT_TRUE_MESSAGE(strstr(output, "Memory") != NULL || 
                           strstr(output, "RSS") != NULL,
                           "Should include memory usage information");
}

// Test runner - specific to this test file
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_single_end_workflow_success);
    RUN_TEST(test_paired_end_workflow_success);
    RUN_TEST(test_corrupted_data_handling);
    RUN_TEST(test_missing_file_error_handling);
    RUN_TEST(test_invalid_arguments_handling);
    RUN_TEST(test_memory_usage_monitoring);
    
    return UNITY_END();
}