# Performance Logging and Monitoring Guide

This guide provides comprehensive documentation for Tronko's advanced performance logging and monitoring capabilities.

## Overview

Tronko includes a sophisticated logging infrastructure designed for:
- **Performance monitoring**: Track memory usage, CPU utilization, and processing speed
- **Debugging assistance**: Detailed milestone tracking and timing information
- **Resource optimization**: Identify bottlenecks and memory leaks
- **Production monitoring**: Comprehensive logging for large-scale deployments

## Quick Start

### Basic Performance Monitoring
```bash
# Enable INFO level logging with resource monitoring
tronko-assign -V2 -R [other options...]

# Enable DEBUG level logging with timing information
tronko-assign -V3 -T [other options...]

# Log to file with full monitoring
tronko-assign -V2 -R -T -l performance.log [other options...]
```

## Logging Levels

### Level 0: ERROR (Default)
Only critical errors are logged.
```bash
tronko-assign [options...]  # Default behavior
```

### Level 1: WARN
Errors and warnings, including potential issues.
```bash
tronko-assign -V1 [options...]
```

### Level 2: INFO (Recommended for Production)
Comprehensive operational information including milestones and progress.
```bash
tronko-assign -V2 [options...]
```

**Example INFO output:**
```
[2025-08-06 16:35:17.655] [INFO] MILESTONE: STARTUP - [0.001s total, 0.001s since last]
[2025-08-06 16:35:17.655] [INFO] MILESTONE: OPTIONS_PARSED - Verbose logging enabled, level=2, cores=1, lines=50000 [0.001s total, 0.000s since last]
[2025-08-06 16:35:17.982] [INFO] MILESTONE: REFERENCE_LOADED - Loaded 1 trees, max_nodename=13, max_taxname=28 [0.328s total, 0.327s since last]
[2025-08-06 16:35:18.040] [INFO] MILESTONE: MEMORY_ALLOCATED - Memory allocated: threads=1, lines_per_batch=50000, total_nodes=2931 [0.386s total, 0.001s since last]
```

### Level 3: DEBUG (Development Only)
Maximum detail including internal algorithm steps and detailed debugging information.
```bash
tronko-assign -V3 [options...]
```

## Resource Monitoring (`-R` flag)

### Memory and CPU Tracking
The `-R` flag enables real-time resource monitoring at each processing milestone:

```bash
tronko-assign -V2 -R [options...]
```

**Example resource output:**
```
[2025-08-06 16:35:19.044] [INFO] MILESTONE: STARTUP - [0.001s total, 0.001s since last] [Memory: 1804 KB, CPU: 103.3%]
[2025-08-06 16:35:19.425] [INFO] MILESTONE: REFERENCE_LOADED - Loaded 1 trees, max_nodename=13, max_taxname=28 [0.382s total, 0.381s since last] [Memory: 53260 KB, CPU: 98.4%]
[2025-08-06 16:35:19.485] [INFO] MILESTONE: MEMORY_ALLOCATED - Memory allocated: threads=1, lines_per_batch=50000, total_nodes=2931 [0.441s total, 0.002s since last] [Memory: 59260 KB, CPU: 97.4%]
```

### Detailed Resource Statistics
Resource monitoring provides additional detailed statistics:

```
[2025-08-06 16:35:19.737] [INFO] Resource Stats [After closing files]: Memory RSS=59.6MB, VM=275.8MB, CPU=97.8%, Wall=0.693s
[2025-08-06 16:35:19.738] [INFO] Resource Stats [Before freeing thread structures]: Memory RSS=59.6MB, VM=275.8MB, CPU=0.0%, Wall=0.000s
[2025-08-06 16:35:19.744] [INFO] Resource Stats [After freeing tree arrays]: Memory RSS=59.6MB, VM=275.7MB, CPU=0.0%, Wall=0.000s
```

### Resource Metrics Explained
- **Memory (KB/MB)**: Current process memory usage
- **RSS**: Resident Set Size (physical memory currently used)
- **VM**: Virtual Memory (total virtual memory used)
- **CPU (%)**: CPU utilization percentage
- **Wall (s)**: Wall clock time elapsed

## Timing Information (`-T` flag)

### Milestone Timing
The `-T` flag enables detailed timing information for performance analysis:

```bash
tronko-assign -V2 -T [options...]
```

**Timing format:**
```
[timestamp] [level] MILESTONE: NAME - Description [total_time since start, time_since_last_milestone]
```

### Key Milestones Tracked

#### Startup and Initialization
- **STARTUP**: Process initialization
- **OPTIONS_PARSED**: Command-line parsing and validation
- **REFERENCE_LOADED**: Database loading and tree construction
- **BWA_INDEX_BUILT**: BWA index creation or loading
- **MEMORY_ALLOCATED**: Memory allocation for processing
- **THREADS_INITIALIZED**: Multi-threading setup

#### Processing Milestones
- **READ_SPECS_DETECTED**: Input file analysis
- **BATCH_START**: Beginning of each processing batch
- **BATCH_LOADED**: Batch data loaded into memory
- **PLACEMENT_COMPLETE**: Phylogenetic placement finished
- **RESULTS_WRITTEN**: Output written to file
- **BATCH_COMPLETE**: Batch processing finished

#### Cleanup Milestones
- **CLEANUP_START**: Beginning cleanup process
- **CLEANUP_COMPLETE**: Cleanup finished
- **PROGRAM_END**: Program termination

### Performance Analysis Examples

#### Memory Usage Patterns
```bash
# Monitor memory allocation and cleanup
tronko-assign -V2 -R -T [options...] 2>&1 | grep -E "(MEMORY_ALLOCATED|Resource Stats|CLEANUP)"
```

#### Processing Speed Analysis
```bash
# Track batch processing performance
tronko-assign -V2 -T [options...] 2>&1 | grep -E "(BATCH_|PLACEMENT_COMPLETE)"
```

#### Bottleneck Identification
```bash
# Find the slowest operations
tronko-assign -V2 -T [options...] 2>&1 | grep "MILESTONE" | sort -k4 -nr
```

## File Logging (`-l` flag)

### Log File Output
```bash
# Log to file (also outputs to stderr)
tronko-assign -V2 -l assignment.log [options...]

# Log with full monitoring
tronko-assign -V3 -R -T -l detailed_performance.log [options...]
```

### Log File Features
- **Dual output**: Logs written to both file and stderr
- **Timestamped entries**: Precise timing for performance analysis
- **Structured format**: Easy parsing for automated analysis
- **Rotation support**: Suitable for long-running processes

## Advanced Use Cases

### Production Monitoring
For production deployments, use INFO level with resource monitoring:
```bash
tronko-assign -V2 -R -l production.log \
  -f database.txt \
  -a reference.fasta \
  -p -1 reads_R1.fastq -2 reads_R2.fastq \
  -o results.txt
```

### Performance Benchmarking
For detailed performance analysis:
```bash
tronko-assign -V2 -R -T -l benchmark.log \
  -f database.txt \
  -a reference.fasta \
  -s -g reads.fasta \
  -o results.txt
```

### Memory Leak Detection
Combined with Valgrind for memory analysis:
```bash
valgrind --leak-check=full --track-origins=yes \
  tronko-assign -V2 -R -l memory_test.log \
  [options...]
```

### Multi-threaded Performance Analysis
Monitor multi-threading efficiency:
```bash
tronko-assign -V2 -R -T -C 8 -l multithread.log \
  [options...]
```

## Log Analysis Tools

### Extracting Performance Metrics
```bash
# Extract milestone timing
grep "MILESTONE" assignment.log | awk '{print $4, $7, $10}'

# Extract memory usage over time
grep "Memory:" assignment.log | awk '{print $2, $8, $10}'

# Extract CPU utilization
grep "CPU:" assignment.log | awk '{print $2, $10}'
```

### Performance Summary Script
```bash
#!/bin/bash
# performance_summary.sh
LOG_FILE=$1

echo "=== Performance Summary ==="
echo "Total runtime:"
grep "PROGRAM_END" $LOG_FILE | awk '{print $7}' | tr -d '[]s'

echo "Peak memory usage:"
grep "Memory:" $LOG_FILE | awk '{print $8}' | sort -n | tail -1

echo "Slowest operations:"
grep "MILESTONE" $LOG_FILE | awk '{print $10, $4}' | sort -nr | head -5
```

### Automated Monitoring
For continuous monitoring, integrate with system monitoring tools:
```bash
# Extract metrics for Prometheus/Grafana
grep "Resource Stats" assignment.log | \
  awk '{print "tronko_memory_rss_mb " $8; print "tronko_cpu_percent " $10}' | \
  tr -d 'MB,%'
```

## Integration with External Tools

### Grafana Dashboard
Monitor Tronko performance in real-time:
1. Configure log shipping to Elasticsearch/Loki
2. Create Grafana queries for memory/CPU metrics
3. Set up alerts for performance degradation

### Log Aggregation
For large-scale deployments:
```bash
# Ship logs to centralized logging
tronko-assign -V2 -R -T [options...] 2>&1 | \
  logger -t tronko-assign
```

### Performance Testing
Automate performance regression testing:
```bash
#!/bin/bash
# performance_test.sh
for cores in 1 2 4 8; do
  echo "Testing with $cores cores..."
  time tronko-assign -V2 -R -T -C $cores \
    -l "perf_${cores}cores.log" \
    [options...]
done
```

## Troubleshooting Performance Issues

### Common Performance Problems

#### High Memory Usage
```bash
# Monitor memory allocation patterns
grep -E "(MEMORY_ALLOCATED|Resource Stats)" assignment.log

# Check for memory leaks
grep "freeing" assignment.log
```

#### Slow Processing
```bash
# Identify bottlenecks
grep "MILESTONE" assignment.log | awk '{print $10, $4}' | sort -nr

# Check BWA index building time
grep "BWA_INDEX_BUILT" assignment.log
```

#### CPU Underutilization
```bash
# Check thread initialization
grep "THREADS_INITIALIZED" assignment.log

# Monitor CPU usage patterns
grep "CPU:" assignment.log | awk '{print $2, $10}'
```

### Performance Optimization Tips

1. **Memory Management**:
   - Use appropriate batch sizes (`-L` flag)
   - Monitor RSS vs VM memory usage
   - Check cleanup milestone timing

2. **CPU Optimization**:
   - Adjust thread count (`-C` flag) based on CPU cores
   - Monitor CPU utilization during processing
   - Consider I/O vs CPU bound operations

3. **I/O Optimization**:
   - Use compressed databases when possible
   - Monitor file loading times
   - Consider SSD storage for large databases

4. **Algorithm Selection**:
   - Compare WFA vs Needleman-Wunsch performance (`-w` flag)
   - Adjust LCA cutoff values (`-c` flag)
   - Monitor placement completion times

## Best Practices

### Production Logging
- Use INFO level (`-V2`) for production systems
- Enable resource monitoring (`-R`) for capacity planning
- Log to files (`-l`) for historical analysis
- Set up log rotation for long-running processes

### Development Logging
- Use DEBUG level (`-V3`) for algorithm development
- Enable timing information (`-T`) for optimization
- Combine with debugging tools (gdb, valgrind)
- Use verbose output for troubleshooting

### Performance Monitoring
- Establish baseline performance metrics
- Monitor trends over time
- Set up alerts for performance degradation
- Regular performance regression testing

This comprehensive logging system provides the visibility needed to optimize Tronko's performance and troubleshoot issues in both development and production environments.