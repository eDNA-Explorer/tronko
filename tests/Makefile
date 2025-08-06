# Makefile for Tronko Tests
# Simple Makefile-based testing (alternative to CMake)

CC = gcc
CFLAGS = -std=gnu99 -Wall -Wextra -g -O0 -fprofile-arcs -ftest-coverage
INCLUDES = -I../tronko-assign -I../tronko-build -Iunity/src
LIBS = -lm -pthread -lz -lrt -lgcov

# Unity framework
UNITY_SRC = unity/src/unity.c
UNITY_OBJ = unity.o

# Common Tronko source files (without main)
TRONKO_SOURCES = ../tronko-assign/readreference.c \
                 ../tronko-assign/assignment.c \
                 ../tronko-assign/placement.c \
                 ../tronko-assign/logger.c \
                 ../tronko-assign/crash_debug.c \
                 ../tronko-assign/symbol_resolver.c \
                 ../tronko-assign/options.c \
                 ../tronko-assign/allocateMemoryForResults.c \
                 ../tronko-assign/allocatetreememory.c \
                 ../tronko-assign/getSequenceinRoot.c \
                 ../tronko-assign/printAlignments.c \
                 ../tronko-assign/resource_monitor.c

TRONKO_OBJECTS = $(TRONKO_SOURCES:.c=.o)

# Test executables
TEST_TARGETS = test_readreference test_crash_debug test_workflows

.PHONY: all clean test coverage

all: $(TEST_TARGETS)

# Unity framework
$(UNITY_OBJ): $(UNITY_SRC)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Individual test executables
test_readreference: unit/test_readreference.c $(UNITY_OBJ) $(TRONKO_OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)

test_crash_debug: unit/test_crash_debug.c $(UNITY_OBJ) $(TRONKO_OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)

test_workflows: integration/test_workflows.c $(UNITY_OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)

# Compile Tronko source files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Run all tests
test: $(TEST_TARGETS)
	@echo "Running unit tests..."
	@for test in $(TEST_TARGETS); do \
		echo "Running $$test..."; \
		./$$test || exit 1; \
	done
	@echo "All tests passed!"

# Generate coverage report
coverage: test
	@echo "Generating coverage report..."
	gcov $(TRONKO_SOURCES)
	lcov --capture --directory . --output-file coverage.info
	genhtml coverage.info --output-directory coverage_html
	@echo "Coverage report generated in coverage_html/"

# Clean build artifacts
clean:
	rm -f $(TEST_TARGETS) $(UNITY_OBJ) $(TRONKO_OBJECTS)
	rm -f *.gcda *.gcno *.gcov coverage.info
	rm -rf coverage_html

# Install test dependencies (Unity framework)
setup:
	@if [ ! -d "unity" ]; then \
		git clone https://github.com/ThrowTheSwitch/Unity.git unity; \
	fi

# Quick smoke test (basic functionality)
smoke: test_crash_debug
	@echo "Running smoke test..."
	./test_crash_debug
	@echo "Smoke test passed!"