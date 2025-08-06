/*
 * Unit tests for readreference.c
 * Tests file reading, parsing, and corruption detection
 */

#include "unity.h"
#include "readreference.h"
#include "crash_debug.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

// Test helper functions
void create_test_fasta_file(const char* filename, const char* content) {
    FILE* file = fopen(filename, "w");
    TEST_ASSERT_NOT_NULL_MESSAGE(file, "Failed to create test file");
    fprintf(file, "%s", content);
    fclose(file);
}

void cleanup_test_file(const char* filename) {
    remove(filename);
}

// Test Cases

void test_corruption_detection_whitespace_sequence(void) {
    const char* test_content = 
        ">valid_header\n"
        "ATCGATCGATCG\n"
        ">corrupt_sequence\n"
        "  ATCGATCGATCG\n";  // Starts with whitespace
    
    const char* test_file = "/tmp/test_corrupt_whitespace.fasta";
    create_test_fasta_file(test_file, test_content);
    
    // TODO: Call readreference function and verify corruption detection
    // This would require refactoring readreference functions to be more testable
    
    cleanup_test_file(test_file);
}

void test_corruption_detection_short_sequence(void) {
    const char* test_content = 
        ">valid_header\n"
        "ATCGATCGATCGATCGATCGATCG\n"
        ">short_sequence\n"
        "ATG\n";  // Very short sequence
    
    const char* test_file = "/tmp/test_corrupt_short.fasta";
    create_test_fasta_file(test_file, test_content);
    
    // TODO: Test short sequence detection
    
    cleanup_test_file(test_file);
}

void test_valid_fasta_parsing(void) {
    const char* test_content = 
        ">sequence_1\n"
        "ATCGATCGATCGATCGATCGATCG\n"
        ">sequence_2\n"
        "GCTAGCTAGCTAGCTAGCTAGCTA\n";
    
    const char* test_file = "/tmp/test_valid.fasta";
    create_test_fasta_file(test_file, test_content);
    
    // TODO: Test valid FASTA parsing
    
    cleanup_test_file(test_file);
}

void test_empty_file_handling(void) {
    const char* test_file = "/tmp/test_empty.fasta";
    create_test_fasta_file(test_file, "");
    
    // TODO: Test empty file handling
    
    cleanup_test_file(test_file);
}

void test_malformed_header_detection(void) {
    const char* test_content = 
        ">\n"  // Empty header
        "ATCGATCGATCGATCGATCGATCG\n";
    
    const char* test_file = "/tmp/test_empty_header.fasta";
    create_test_fasta_file(test_file, test_content);
    
    // TODO: Test empty header detection
    
    cleanup_test_file(test_file);
}

// Test runner - specific to this test file
int main(void) {
    UNITY_BEGIN();
    
    RUN_TEST(test_corruption_detection_whitespace_sequence);
    RUN_TEST(test_corruption_detection_short_sequence);
    RUN_TEST(test_valid_fasta_parsing);
    RUN_TEST(test_empty_file_handling);
    RUN_TEST(test_malformed_header_detection);
    
    return UNITY_END();
}