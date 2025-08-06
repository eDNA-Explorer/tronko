# Crash Debugging System Guide

This guide provides comprehensive documentation for Tronko's advanced crash debugging and error analysis capabilities.

## Overview

Tronko includes a sophisticated crash debugging system designed to:
- **Detect crashes**: Automatic detection of segmentation faults and other critical errors
- **Root cause analysis**: Identify the original source of crashes, including data corruption
- **Detailed forensics**: Comprehensive crash reports with stack traces and system state
- **Real-time monitoring**: Track application context and processing state during execution

## Quick Start

### Basic Crash Detection
The crash debugging system is automatically enabled when you run tronko-assign:
```bash
# Crash debugging is always active
tronko-assign [options...]
```

### Enhanced Crash Reporting
For maximum crash analysis detail:
```bash
# Enable verbose logging to see crash context
tronko-assign -V2 -R -T [options...]
```

### Testing Crash Detection
Verify crash debugging functionality with intentionally corrupted data:
```bash
# Run Test 6 from the test suite (includes crash testing)
./test_with_example_data.sh
```

## Crash Detection Capabilities

### Signals Monitored
The crash debugging system automatically detects these critical signals:

- **SIGSEGV**: Segmentation fault (memory access violations)
- **SIGABRT**: Abnormal termination (abort() calls)
- **SIGBUS**: Bus error (alignment or non-existent address access)
- **SIGFPE**: Floating-point exception (division by zero, overflow)
- **SIGILL**: Illegal instruction (corrupted executable or hardware issues)

### Automatic Features
- **Signal handlers**: Installed at program startup
- **Stack trace capture**: Automatic backtrace generation
- **Context preservation**: Application state at crash time
- **Core dump generation**: Traditional debugging support
- **Safe crash handling**: Prevents recursive crashes during reporting

## Crash Report Generation

### Report Locations
Crash reports are automatically written to:
- **Crash files**: `/tmp/tronko_assign_crash_[PID]_[TIMESTAMP].crash`
- **Log output**: Integrated with the logging system (if enabled)
- **Console output**: Summary displayed on stderr

### Report Contents
Each crash report includes:

#### 1. Basic Crash Information
```
=== TRONKO CRASH REPORT ===
Crash Time: Wed Aug  6 15:09:23 2025
Process ID: 1811
Program: /app/tronko-assign/tronko-assign
Working Directory: /app/tronko-assign
Command Line: tronko-assign -r -f reference.txt -a sequences.fasta -p -1 reads1.fasta -2 reads2.fasta -o results.txt

CRASH: SIGSEGV (Segmentation fault) at address (nil) in process 1811
```

#### 2. Application Context
```
Application Context:
  Processing Stage: Reading paired-end input files
  Current File: ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read2.fasta (line 68)
  Current Read: read_001 (batch 1, index 45)
  Current Tree: 0
```

#### 3. Root Cause Analysis (If Data Corruption Detected)
```
  *** DATA CORRUPTION DETECTED ***
  Corrupted File: ../example_datasets/multiple_trees/missingreads_pairedend_150bp_2error_read1.fasta (line 26)
  Corruption Type: Sequence line starts with whitespace
```

#### 4. Stack Trace with Symbol Resolution
```
Stack Trace:
#0  ./tronko-assign(+0x59ba4) [0xaaaab9559ba4]
    Object: ./tronko-assign
#1  linux-vdso.so.1(__kernel_rt_sigreturn+0) [0xffff843e07a0]
    Symbol: __kernel_rt_sigreturn
    Object: linux-vdso.so.1
#2  ./tronko-assign(+0x2ab28) [0xaaaab952ab28]
    Object: ./tronko-assign
```

#### 5. System State Information
```
Register state:
[Platform-specific register dumps]

Memory state:
RSS: 65764 KB
VM Size: 287504 KB
CPU Usage: 97.54%
```

## Data Corruption Detection

### Automatic Validation
The crash debugging system includes real-time data validation that detects:

#### FASTA Format Issues
- **Whitespace-prefixed sequences**: Invalid leading spaces/tabs
- **Empty headers**: Headers without sequence names
- **Malformed structure**: Missing '>' characters or blank lines

#### FASTQ Format Issues  
- **Invalid headers**: Missing '@' characters
- **Whitespace issues**: Leading spaces in sequence/quality lines
- **Short sequences**: Suspiciously short reads (< 10 nucleotides)

#### File Structure Problems
- **Encoding issues**: Non-ASCII characters in critical sections
- **Truncated files**: Unexpected end-of-file conditions
- **Binary corruption**: Invalid characters in text files

### Corruption Types Detected
The system identifies specific corruption patterns:

```
Corruption Type: "Sequence line starts with whitespace"
Corruption Type: "Suspiciously short sequence"  
Corruption Type: "Empty FASTA header"
Corruption Type: "FASTQ line starts with whitespace"
Corruption Type: "Empty FASTQ header"
```

### Root Cause Tracking
The corruption detection system:
1. **Flags corruption** during file reading
2. **Preserves location** (file and line number)
3. **Tracks through processing** even when crash occurs later
4. **Links crashes to original corruption** source

## Context Tracking

### Real-time Application State
The crash debugging system continuously tracks:

#### Processing Context
- **Current stage**: What operation was being performed
- **Current file**: Which file was being processed
- **Current line**: Exact line number in the file
- **Current read**: Which sequence was being processed
- **Batch information**: Processing batch and position

#### Example Context Tracking
```bash
# During normal processing, context is updated:
crash_set_processing_stage("Reading paired-end input files");
crash_set_current_file(opt.read1_file);
crash_set_current_file_line(opt.read1_file, line_number);
crash_set_current_read(read_name, batch_number, read_index);
```

### Thread Safety
All context tracking is thread-safe using pthread mutexes:
- **Atomic updates**: Context changes are atomic
- **No race conditions**: Safe for multi-threaded processing
- **Consistent state**: Context always reflects current operation

## Advanced Features

### Symbol Resolution
Enhanced debugging with function name resolution:

#### addr2line Integration
- **Automatic lookup**: Function names from addresses
- **Source file information**: File and line numbers (with debug symbols)
- **Symbol caching**: Performance optimization for repeated lookups

#### Debug Symbol Support
For enhanced symbol resolution, compile with debug symbols:
```bash
cd tronko-assign
make debug
```

### Core Dump Generation
Traditional debugging support:
- **Automatic core dumps**: Generated alongside crash reports
- **GDB integration**: Use with debuggers
- **System configuration**: Respects system core dump settings

### Crash Report Analysis
Tools for analyzing crash reports:

#### Extract Key Information
```bash
# Get crash summary
grep "CRASH:" /tmp/tronko_assign_crash_*.crash

# Find corruption information
grep -A 5 "DATA CORRUPTION DETECTED" /tmp/tronko_assign_crash_*.crash

# Extract stack trace
sed -n '/Stack Trace:/,/^$/p' /tmp/tronko_assign_crash_*.crash
```

#### Automated Analysis Script
```bash
#!/bin/bash
# analyze_crash.sh
CRASH_FILE=$1

echo "=== Crash Analysis ==="
echo "Crash type:"
grep "CRASH:" $CRASH_FILE | awk -F': ' '{print $2}'

echo "Crash location:"
grep "Current File:" $CRASH_FILE

echo "Root cause:"
if grep -q "DATA CORRUPTION DETECTED" $CRASH_FILE; then
    grep -A 3 "DATA CORRUPTION DETECTED" $CRASH_FILE
else
    echo "No data corruption detected"
fi

echo "Stack trace summary:"
grep -A 5 "Stack Trace:" $CRASH_FILE | head -10
```

## Common Crash Scenarios

### 1. Data Corruption Crashes
**Symptoms**: Crashes during file reading or processing
**Cause**: Malformed input files
**Resolution**: Check crash report for corrupted file/line

**Example**:
```
Current File: reads.fasta (line 123)
*** DATA CORRUPTION DETECTED ***
Corrupted File: reads.fasta (line 89)
Corruption Type: Sequence line starts with whitespace
```

**Fix**: Edit the input file to remove leading whitespace from line 89.

### 2. Memory Exhaustion
**Symptoms**: Crashes during memory allocation
**Cause**: Insufficient system memory or memory leaks
**Resolution**: Monitor memory usage with `-R` flag

**Example**:
```
Processing Stage: Memory allocation
Memory state: RSS: 8GB, VM Size: 12GB
```

**Fix**: Increase system memory or reduce batch size with `-L` flag.

### 3. Threading Issues
**Symptoms**: Crashes in multi-threaded mode
**Cause**: Race conditions or thread safety issues
**Resolution**: Reduce thread count or run single-threaded

**Example**:
```
Processing Stage: Multi-threaded assignment
Current Thread: 4 of 8
```

**Fix**: Use `-C 1` for single-threaded processing to isolate the issue.

### 4. Algorithm Errors
**Symptoms**: Crashes during placement or assignment
**Cause**: Edge cases in phylogenetic algorithms
**Resolution**: Try alternative algorithms

**Example**:
```
Processing Stage: Phylogenetic placement
Current Read: problematic_sequence
```

**Fix**: Use Needleman-Wunsch alignment (`-w` flag) instead of WFA.

## Debugging Workflows

### 1. Initial Crash Investigation
```bash
# Run with verbose logging to capture context
tronko-assign -V2 -R -T [options...] 2>&1 | tee debug.log

# Check for crash reports
ls -la /tmp/tronko_assign_crash_*.crash

# Analyze the most recent crash
./analyze_crash.sh /tmp/tronko_assign_crash_*.crash
```

### 2. Data Validation
```bash
# Test with known good data first
tronko-assign [options...] -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta

# If successful, validate your input data
head -50 your_input_file.fasta  # Check first 50 lines for issues
grep -n "^[[:space:]]" your_input_file.fasta  # Find lines with leading whitespace
```

### 3. Incremental Testing
```bash
# Test with smaller datasets
head -1000 your_input_file.fasta > test_subset.fasta
tronko-assign [options...] -g test_subset.fasta

# Gradually increase dataset size
head -10000 your_input_file.fasta > test_larger.fasta
```

### 4. Memory Debugging
```bash
# Use Valgrind for detailed memory analysis
valgrind --leak-check=full --track-origins=yes \
  tronko-assign -V2 [options...] 2>&1 | tee valgrind.log

# Check for memory leaks
grep "definitely lost" valgrind.log
```

## Integration with Development Tools

### GDB Integration
Debug with GDB using crash information:
```bash
# Compile with debug symbols
cd tronko-assign && make debug

# Run under GDB
gdb ./tronko-assign
(gdb) run [options...]
(gdb) bt  # Get backtrace on crash
(gdb) info registers  # Check register state
```

### Automated Testing
Include crash testing in CI/CD:
```bash
# Test crash detection with corrupted data
./test_with_example_data.sh

# Verify crash reports are generated
ls /tmp/tronko_assign_crash_*.crash || echo "No crashes detected"
```

### Log Integration
Crash information integrates with the logging system:
```bash
# Crashes appear in log files when verbose logging is enabled
tronko-assign -V2 -l assignment.log [options...]
grep "CRASH DETECTED" assignment.log
```

## Configuration and Customization

### Crash Report Directory
Default crash reports go to `/tmp/`. This can be customized by modifying the source code in `crash_debug.c`:

```c
// Default configuration
static crash_config_t g_crash_config = {
    .crash_log_directory = "/tmp",
    .crash_log_prefix = "tronko_assign_crash",
    // ...
};
```

### Signal Handler Customization
The system can be extended to handle additional signals or provide custom crash behavior by modifying `crash_debug.c`.

### Debugging Output Control
Control the level of debugging output:
- **Minimal**: Basic crash detection (always enabled)
- **Standard**: Include application context (default)
- **Verbose**: Add detailed system information (`-V2`)
- **Maximum**: Include all debugging data (`-V3`)

## Best Practices

### For Users
1. **Always check crash reports**: Look in `/tmp/` for `.crash` files
2. **Enable verbose logging**: Use `-V2` to capture context
3. **Validate input data**: Check for common corruption patterns
4. **Test with subsets**: Use smaller datasets to isolate issues
5. **Report crashes**: Include crash reports when reporting bugs

### For Developers
1. **Preserve crash context**: Update context tracking when modifying code
2. **Test crash scenarios**: Include intentional crash tests
3. **Handle edge cases**: Add validation for new input formats
4. **Update documentation**: Keep crash analysis guides current
5. **Monitor production**: Set up alerts for crash detection

### For System Administrators
1. **Monitor crash reports**: Set up log rotation for `/tmp/`
2. **Resource monitoring**: Use `-R` flag for capacity planning
3. **Automated analysis**: Script crash report processing
4. **Alerting systems**: Integrate with monitoring infrastructure
5. **Performance baselines**: Establish normal operation patterns

This comprehensive crash debugging system provides the tools needed to quickly identify, analyze, and resolve issues in Tronko, ensuring reliable operation in both development and production environments.