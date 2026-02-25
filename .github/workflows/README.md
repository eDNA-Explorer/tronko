# GitHub Actions CI/CD Workflows

This directory contains automated testing workflows for the Tronko project.

## Workflow Overview

### 1. Quick Tests (`quick-test.yml`)
**Triggers**: Every push to any branch, all pull requests
**Purpose**: Fast feedback for developers
**Duration**: ~2-3 minutes

**What it runs:**
- Build verification (tronko-build & tronko-assign)
- Unit tests (crash debugging system)
- Integration tests (binary & data validation)
- Quick functional test (subset of full test suite)

### 2. Comprehensive Tests (`test.yml`)
**Triggers**: Push to main/experimental/develop branches, PRs to main/experimental
**Purpose**: Full validation before merging
**Duration**: ~10-15 minutes

**Jobs:**
- **unit-tests**: Complete unit test suite with coverage
- **integration-tests**: Full integration testing with memory leak detection
- **comprehensive-tests**: End-to-end testing including crash debugging
- **docker-tests**: Containerized testing for environment consistency
- **performance-tests**: Benchmarking (main/experimental branches only)

## Test Categories

### Unit Tests
- **Framework**: Unity (lightweight C testing framework)
- **Coverage**: Crash debugging system, corruption detection, signal handling
- **Files**: `tests/unit/test_*.c`
- **Command**: `make -f Makefile.simple smoke`

### Integration Tests  
- **Coverage**: Binary existence, command-line interfaces, error handling
- **Files**: `tests/integration/test_*.c`
- **Command**: `make -f Makefile.simple test_simple_workflows`

### Functional Tests
- **Coverage**: End-to-end workflows with real data
- **Files**: `test_with_example_data.sh`
- **Includes**: Crash testing with corrupted data

### Performance Tests
- **Coverage**: Memory usage, processing speed, resource monitoring
- **Triggers**: Only on main/experimental branch pushes
- **Output**: Artifacts with performance logs

## Status Badges

Add these to your README.md to show build status:

```markdown
![Quick Tests](https://github.com/USERNAME/REPO/workflows/Quick%20Tests/badge.svg)
![Comprehensive Tests](https://github.com/USERNAME/REPO/workflows/Tests/badge.svg)
```

## Reading Test Results

### ✅ Successful Run
All tests passed, code is ready for review/merge.

### ❌ Failed Run
Click on the failed job to see details:
- **Build failures**: Check compilation errors
- **Unit test failures**: Review test output for specific assertion failures
- **Integration failures**: Check binary/data availability
- **Memory leaks**: Review Valgrind output
- **Performance regressions**: Check benchmark comparisons

### ⚠️ Warnings
- Coverage upload failures (non-blocking)
- Performance tests skipped (not on main/experimental)
- Expected crashes in functional tests

## Local Testing

Run the same tests locally:

```bash
# Quick unit tests
cd tests && make -f Makefile.simple smoke

# Full test suite  
cd tests && make -f Makefile.simple test

# Functional tests
./test_with_example_data.sh

# Using test runner
./run_tests.sh --unit-only
./run_tests.sh --integration-only
./run_tests.sh --coverage
```

## Docker Testing

```bash
# Build and test in Docker (matches CI environment)
docker build -t tronko:test .
docker run --rm -v $PWD:/workspace tronko:test bash -c "
  cd /workspace/tests && make -f Makefile.simple test
"
```

## Workflow Configuration

### Adding New Tests
1. Add test files to `tests/unit/` or `tests/integration/`
2. Update `Makefile.simple` with new test targets
3. Tests will automatically run in CI

### Branch Protection
Recommended branch protection rules for main/experimental:
- Require status checks: "Quick Tests"
- Require up-to-date branches
- Require review from code owners

### Performance Monitoring
Performance tests run on main/experimental and upload artifacts:
- Memory usage benchmarks
- Processing speed metrics  
- Resource utilization data

## Troubleshooting

### Common Issues

**Build Failures:**
- Check for missing dependencies in workflow files
- Verify Makefile targets exist
- Review compiler warnings/errors

**Test Failures:**
- Unit tests: Check test logic and assertions
- Integration tests: Verify test data and binaries
- Functional tests: May include expected crashes (check logs)

**Timeout Issues:**
- Functional tests have 300s timeout
- Quick tests have 120s timeout
- Increase if needed for complex datasets

**Docker Issues:**
- Verify Dockerfile builds locally
- Check volume mounts and permissions
- Review container resource limits

### Getting Help

1. Check workflow logs in GitHub Actions tab
2. Run tests locally to reproduce issues
3. Review this documentation for workflow details
4. Check individual test README files for specific guidance