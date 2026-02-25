# Tronko Documentation

This directory contains comprehensive documentation for advanced Tronko features and development workflows.

## Documentation Overview

### Core Features
- **[Performance Logging](performance-logging.md)**: Comprehensive guide to performance monitoring, resource tracking, and optimization
- **[Crash Debugging](crash-debugging.md)**: Advanced crash detection, root cause analysis, and debugging workflows

### Development
- **[Testing Guide](../tests/README.md)**: Unit testing, integration testing, and CI/CD workflows
- **[GitHub Actions](../.github/workflows/README.md)**: Automated testing and continuous integration setup

## Quick Navigation

### Performance Monitoring
For optimizing Tronko performance and monitoring resource usage:
- [Getting Started with Performance Logging](performance-logging.md#quick-start)
- [Resource Monitoring Options](performance-logging.md#resource-monitoring--r-flag)
- [Timing Analysis](performance-logging.md#timing-information--t-flag)
- [Production Monitoring Best Practices](performance-logging.md#best-practices)

### Debugging Crashes
For investigating and resolving crashes and errors:
- [Crash Detection Overview](crash-debugging.md#crash-detection-capabilities)
- [Understanding Crash Reports](crash-debugging.md#crash-report-generation)
- [Data Corruption Analysis](crash-debugging.md#data-corruption-detection)
- [Common Crash Scenarios](crash-debugging.md#common-crash-scenarios)

### Testing and Development
For developers working on Tronko:
- [Running Tests](../tests/README.md#quick-start)
- [Unit Testing Framework](../tests/README.md#test-structure)
- [CI/CD Workflows](../.github/workflows/README.md#workflow-overview)
- [Docker Development](../tests/README.md#docker-testing)

## Feature Matrix

| Feature | Basic Usage | Advanced Usage | Documentation |
|---------|-------------|----------------|---------------|
| **Performance Logging** | `-V2` | `-V3 -R -T -l performance.log` | [Performance Guide](performance-logging.md) |
| **Crash Debugging** | Always enabled | `-V2` for context | [Crash Debug Guide](crash-debugging.md) |
| **Resource Monitoring** | `-R` | `-R -T` with analysis | [Resource Monitoring](performance-logging.md#resource-monitoring--r-flag) |
| **Unit Testing** | `make -f Makefile.simple smoke` | `./run_tests.sh --coverage` | [Testing Guide](../tests/README.md) |
| **Integration Testing** | `./test_with_example_data.sh` | CI/CD automation | [GitHub Actions](../.github/workflows/README.md) |

## Common Workflows

### Performance Optimization
1. [Enable performance logging](performance-logging.md#basic-performance-monitoring)
2. [Run with resource monitoring](performance-logging.md#resource-monitoring--r-flag)
3. [Analyze timing data](performance-logging.md#performance-analysis-examples)
4. [Optimize based on bottlenecks](performance-logging.md#troubleshooting-performance-issues)

### Debugging Issues
1. [Reproduce with verbose logging](crash-debugging.md#enhanced-crash-reporting)
2. [Check crash reports](crash-debugging.md#crash-report-generation)
3. [Analyze data corruption](crash-debugging.md#data-corruption-detection)
4. [Follow debugging workflows](crash-debugging.md#debugging-workflows)

### Development Testing
1. [Run unit tests](../tests/README.md#unit-tests)
2. [Execute integration tests](../tests/README.md#integration-tests)
3. [Validate with functional tests](../tests/README.md#functional-tests)
4. [Check CI/CD status](../.github/workflows/README.md#status-badges)

## Command Reference

### Logging Commands
```bash
# Basic verbose logging
tronko-assign -V2 [options...]

# Full performance monitoring
tronko-assign -V2 -R -T -l performance.log [options...]

# Debug level logging
tronko-assign -V3 -R -T [options...]
```

### Testing Commands
```bash
# Quick unit tests
cd tests && make -f Makefile.simple smoke

# Full test suite
./run_tests.sh

# Docker testing
docker compose exec tronko-dev bash -c "cd /app/tests && make -f Makefile.simple test"
```

### Debugging Commands
```bash
# Test crash detection
./test_with_example_data.sh

# Analyze crash reports
ls /tmp/tronko_assign_crash_*.crash

# Memory debugging
valgrind --leak-check=full tronko-assign [options...]
```

## Troubleshooting Quick Reference

### Performance Issues
| Problem | Solution | Documentation |
|---------|----------|---------------|
| High memory usage | Monitor with `-R`, check allocations | [Memory Analysis](performance-logging.md#memory-usage-patterns) |
| Slow processing | Use `-T` to identify bottlenecks | [Performance Analysis](performance-logging.md#bottleneck-identification) |
| CPU underutilization | Adjust `-C` thread count | [CPU Optimization](performance-logging.md#cpu-optimization) |

### Crash Issues
| Problem | Solution | Documentation |
|---------|----------|---------------|
| Segmentation fault | Check crash report for context | [Crash Investigation](crash-debugging.md#initial-crash-investigation) |
| Data corruption | Validate input files | [Corruption Detection](crash-debugging.md#data-corruption-detection) |
| Memory errors | Use Valgrind for analysis | [Memory Debugging](crash-debugging.md#memory-debugging) |

### Testing Issues
| Problem | Solution | Documentation |
|---------|----------|---------------|
| Test failures | Check individual test output | [Troubleshooting Tests](../tests/README.md#troubleshooting) |
| CI/CD failures | Review GitHub Actions logs | [CI/CD Troubleshooting](../.github/workflows/README.md#troubleshooting) |
| Docker issues | Verify container setup | [Docker Testing](../tests/README.md#docker-testing) |

## Advanced Topics

### Custom Development
- [Extending the logging system](performance-logging.md#integration-with-external-tools)
- [Adding new crash detection](crash-debugging.md#configuration-and-customization)
- [Creating custom tests](../tests/README.md#adding-new-tests)

### Production Deployment
- [Production logging setup](performance-logging.md#production-logging)
- [Monitoring integration](performance-logging.md#grafana-dashboard)
- [Automated analysis](crash-debugging.md#automated-analysis-script)

### Research and Development
- [Performance benchmarking](performance-logging.md#performance-benchmarking)
- [Algorithm optimization](performance-logging.md#performance-optimization-tips)
- [Debugging complex issues](crash-debugging.md#advanced-features)

## Contributing

When contributing to Tronko:
1. **Add tests** for new features using the [testing framework](../tests/README.md)
2. **Update documentation** for any new logging or debugging features
3. **Follow performance monitoring** best practices during development
4. **Ensure crash safety** by testing edge cases and error conditions

## Getting Help

1. **Documentation**: Start with the relevant guide above
2. **Testing**: Use the test suite to reproduce issues
3. **Logging**: Enable verbose logging to understand program behavior
4. **Community**: Check GitHub issues and discussions
5. **Support**: Include crash reports and performance logs when reporting issues

This documentation provides comprehensive coverage of Tronko's advanced features, enabling effective development, debugging, and production deployment.