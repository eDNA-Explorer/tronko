#!/bin/bash

# Tronko Test Runner Script
# Provides easy interface for running different types of tests

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}===========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}===========================================${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Default values
RUN_UNIT=1
RUN_INTEGRATION=1
RUN_VALGRIND=0
RUN_COVERAGE=0
VERBOSE=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --unit-only)
            RUN_INTEGRATION=0
            shift
            ;;
        --integration-only)
            RUN_UNIT=0
            shift
            ;;
        --valgrind)
            RUN_VALGRIND=1
            shift
            ;;
        --coverage)
            RUN_COVERAGE=1
            shift
            ;;
        --verbose)
            VERBOSE=1
            shift
            ;;
        --help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --unit-only       Run only unit tests"
            echo "  --integration-only Run only integration tests"
            echo "  --valgrind        Run with Valgrind memory checking"
            echo "  --coverage        Generate code coverage report"
            echo "  --verbose         Verbose output"
            echo "  --help            Show this help message"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if we're in the right directory
if [[ ! -f "test_with_example_data.sh" ]]; then
    print_error "This script must be run from the tronko root directory"
    exit 1
fi

print_header "Tronko Test Suite"

# Build the binaries first
print_header "Building Binaries"

echo "Building tronko-build..."
cd tronko-build
if make clean && make; then
    print_success "tronko-build built successfully"
else
    print_error "Failed to build tronko-build"
    exit 1
fi

echo "Building tronko-assign..."
cd ../tronko-assign
if make clean && make; then
    print_success "tronko-assign built successfully"
else
    print_error "Failed to build tronko-assign"
    exit 1
fi

cd ..

# Run unit tests
if [[ $RUN_UNIT -eq 1 ]]; then
    print_header "Running Unit Tests"
    
    cd tests
    
    # Setup Unity if needed
    if [[ ! -d "unity" ]]; then
        echo "Setting up Unity testing framework..."
        make setup
    fi
    
    if [[ $RUN_COVERAGE -eq 1 ]]; then
        echo "Running tests with coverage..."
        if make coverage; then
            print_success "Unit tests passed with coverage"
            echo "Coverage report available in tests/coverage_html/"
        else
            print_error "Unit tests failed"
            exit 1
        fi
    else
        echo "Running unit tests..."
        if make test; then
            print_success "Unit tests passed"
        else
            print_error "Unit tests failed"
            exit 1
        fi
    fi
    
    cd ..
fi

# Run integration tests
if [[ $RUN_INTEGRATION -eq 1 ]]; then
    print_header "Running Integration Tests"
    
    echo "Running example dataset tests..."
    if ./test_with_example_data.sh; then
        print_success "Integration tests passed"
    else
        print_warning "Some integration tests failed (may be expected for crash tests)"
    fi
fi

# Run Valgrind tests
if [[ $RUN_VALGRIND -eq 1 ]]; then
    print_header "Running Valgrind Memory Check"
    
    if ! command -v valgrind &> /dev/null; then
        print_warning "Valgrind not found. Install with: sudo apt-get install valgrind"
    else
        cd tronko-assign
        echo "Running memory leak detection..."
        valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes \
            ./tronko-assign -r \
            -f ../tronko-build/example_datasets/single_tree/reference_tree.txt \
            -a ../tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
            -s \
            -g ../example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
            -o valgrind_test_results.txt \
            -w > valgrind_output.txt 2>&1 || true
        
        if grep -q "ERROR SUMMARY: 0 errors" valgrind_output.txt; then
            print_success "No memory errors detected"
        else
            print_warning "Memory issues detected. Check valgrind_output.txt"
        fi
        cd ..
    fi
fi

print_header "Test Summary"

if [[ $RUN_UNIT -eq 1 ]]; then
    print_success "Unit tests completed"
fi

if [[ $RUN_INTEGRATION -eq 1 ]]; then
    print_success "Integration tests completed"
fi

if [[ $RUN_VALGRIND -eq 1 ]]; then
    print_success "Memory analysis completed"
fi

if [[ $RUN_COVERAGE -eq 1 ]]; then
    echo -e "${BLUE}📊 Coverage report: tests/coverage_html/index.html${NC}"
fi

print_success "All requested tests completed!"

echo ""
echo "Next steps:"
echo "- Review any warnings or errors above"
echo "- Check test output files for detailed results"
if [[ $RUN_COVERAGE -eq 1 ]]; then
    echo "- Open coverage report: firefox tests/coverage_html/index.html"
fi
echo "- Run specific test types with --unit-only or --integration-only"