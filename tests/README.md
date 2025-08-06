# Tronko Testing Suite

This directory contains comprehensive tests for the Tronko phylogenetic classification software.

## Overview

The testing suite includes:
- **Unit tests**: Test individual functions and modules
- **Integration tests**: Test end-to-end workflows
- **Performance tests**: Benchmark performance and memory usage
- **Crash tests**: Verify crash debugging functionality

## Quick Start

### Run All Tests
```bash
# From tronko root directory
./run_tests.sh
```

### Run Specific Test Types
```bash
./run_tests.sh --unit-only          # Only unit tests
./run_tests.sh --integration-only   # Only integration tests
./run_tests.sh --coverage          # Generate coverage report
./run_tests.sh --valgrind          # Memory leak detection
```

## Test Framework

We use **Unity** - a lightweight C testing framework well-suited for embedded/scientific software:
- Simple assertion macros
- Automatic test discovery
- No external dependencies
- Good performance

## Test Structure

```
tests/
├── unit/                    # Unit tests
│   ├── test_readreference.c # File I/O and parsing tests
│   ├── test_crash_debug.c   # Crash debugging tests
│   └── test_assignment.c    # Classification algorithm tests
├── integration/             # Integration tests
│   └── test_workflows.c     # End-to-end workflow tests
├── fixtures/                # Test data files
├── unity/                   # Unity testing framework
├── Makefile                 # Build system for tests
├── CMakeLists.txt          # Alternative CMake build
└── README.md               # This file
```

## Writing Tests

### Unit Test Example
```c
#include "unity.h"
#include "module_to_test.h"

void setUp(void) {
    // Initialize before each test
}

void tearDown(void) {
    // Cleanup after each test
}

void test_function_should_return_expected_value(void) {
    // Arrange
    int input = 5;
    int expected = 10;
    
    // Act
    int result = function_under_test(input);
    
    // Assert
    TEST_ASSERT_EQUAL_INT(expected, result);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_function_should_return_expected_value);
    return UNITY_END();
}
```

### Common Assertions
```c
TEST_ASSERT_EQUAL_INT(expected, actual);
TEST_ASSERT_EQUAL_STRING(expected, actual);
TEST_ASSERT_NOT_NULL(pointer);
TEST_ASSERT_TRUE(condition);
TEST_ASSERT_FALSE(condition);
TEST_ASSERT_EQUAL_FLOAT(expected, actual, delta);
```

## Test Categories

### 1. Unit Tests

**File I/O Tests** (`test_readreference.c`):
- FASTA/FASTQ parsing
- Corruption detection
- Error handling
- Memory management

**Crash Debugging Tests** (`test_crash_debug.c`):
- Signal handling
- Context tracking
- Corruption flagging
- Report generation

**Algorithm Tests** (`test_assignment.c`):
- LCA calculations
- Placement algorithms
- Scoring functions

### 2. Integration Tests

**Workflow Tests** (`test_workflows.c`):
- Single-end read processing
- Paired-end read processing
- Error handling
- Memory monitoring

### 3. Performance Tests

**Benchmarking**:
- Processing speed
- Memory usage
- Scalability testing

## CI/CD Integration

### GitHub Actions
The `.github/workflows/test.yml` file defines our CI pipeline:
- **Unit tests**: Run on every push/PR
- **Integration tests**: Full workflow testing
- **Docker tests**: Containerized testing
- **Performance tests**: Benchmarking on main branch
- **Coverage reporting**: Code coverage analysis

### Local CI Simulation
```bash
# Run the same tests as CI locally
./run_tests.sh --coverage --valgrind
```

## Coverage Reporting

Generate code coverage reports:
```bash
cd tests
make coverage
```

View coverage report:
```bash
firefox coverage_html/index.html
```

## Memory Testing

### Valgrind Integration
```bash
./run_tests.sh --valgrind
```

This runs memory leak detection and reports:
- Memory leaks
- Invalid memory access
- Uninitialized memory usage

### AddressSanitizer (Alternative)
Compile with AddressSanitizer for faster memory error detection:
```bash
cd tronko-assign
make CFLAGS="-fsanitize=address -g" clean all
```

## Debugging Tests

### Running Individual Tests
```bash
cd tests
make test_readreference
./test_readreference
```

### Debugging with GDB
```bash
cd tests
gdb ./test_readreference
(gdb) run
(gdb) bt  # backtrace on failure
```

### Verbose Output
```bash
./run_tests.sh --verbose
```

## Best Practices

### Test Naming
- Use descriptive test names: `test_function_should_behavior_when_condition`
- Group related tests in the same file
- Use setUp/tearDown for common initialization

### Test Data
- Keep test data small and focused
- Use fixtures for reusable test data
- Clean up temporary files in tearDown

### Assertions
- Use specific assertions (TEST_ASSERT_EQUAL_INT vs TEST_ASSERT_TRUE)
- Include meaningful failure messages
- Test both positive and negative cases

### Coverage Goals
- Aim for >80% line coverage on critical modules
- Focus on testing error conditions
- Test boundary conditions and edge cases

## Adding New Tests

1. **Create test file** in appropriate directory (`unit/` or `integration/`)
2. **Include headers** for Unity and module under test
3. **Implement setUp/tearDown** for test initialization
4. **Write test functions** following naming conventions
5. **Add to Makefile** in TEST_TARGETS
6. **Update CI** if needed for new test dependencies

## Troubleshooting

### Common Issues

**Tests fail to compile**:
- Check include paths in Makefile
- Verify all dependencies are built
- Check for missing headers

**Tests fail unexpectedly**:
- Run individual tests to isolate failures
- Use GDB for debugging
- Check test setup/teardown

**Memory tests fail**:
- Review Valgrind output for specific issues
- Check for unfreed memory in tearDown
- Verify proper cleanup in error paths

### Getting Help

1. Check this README for common solutions
2. Review Unity documentation: https://github.com/ThrowTheSwitch/Unity
3. Look at existing tests for examples
4. Use `./run_tests.sh --help` for options