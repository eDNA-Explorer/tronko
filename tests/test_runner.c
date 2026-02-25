/* Test Runner for Tronko Tests
 * This file provides the main() function for Unity-based tests
 */

#include "unity.h"

// External test function declarations
// These will be defined in individual test files
extern void setUp(void);
extern void tearDown(void);

// Test function declarations - add your test functions here
// Example: extern void test_function_name(void);

int main(void) {
    UNITY_BEGIN();
    
    // Add your test function calls here
    // Example: RUN_TEST(test_function_name);
    
    return UNITY_END();
}