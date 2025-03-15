# Discarded Reads Logging System Implementation Plan

## Overview

This document outlines a plan to implement a comprehensive logging system for tracking discarded reads in Tronko. The primary goal is to provide detailed information about which reads are being filtered out during the taxonomic assignment process and the specific reasons for their exclusion. This information will be invaluable for troubleshooting, quality control, and refining parameters for optimal performance.

## Problem Statement

Currently, Tronko silently discards reads that fail to meet various quality thresholds at different stages of the pipeline. There is no built-in mechanism to:

1. Identify which reads were discarded
2. Determine at which stage each read was filtered out
3. Understand the specific reason for rejection
4. Track the percentage of reads lost at each filtering stage
5. Correlate discarded reads with specific parameters or thresholds

This lack of information makes it difficult to troubleshoot when a significant portion of reads are unassigned, or to optimize the pipeline for specific datasets.

## Proposed Solution

We propose implementing a comprehensive logging system that captures detailed information about discarded reads at each filtering stage in the Tronko-assign pipeline.

### Key Components

1. **Configurable Logging Level**: Allow users to set the verbosity of discarded read logging
2. **Structured Log Output**: Generate a standardized, tab-delimited log file for easy parsing
3. **Multiple Filtering Stage Tracking**: Capture information at each distinct filtering point
4. **Threshold Reporting**: Record the relevant thresholds and actual values that led to filtering
5. **Summary Statistics**: Provide summary counts and percentages of reads discarded at each stage

## Filtering Stages to Monitor

Based on our analysis of the Tronko codebase, reads are discarded at four main stages:

### 1. BWA Alignment Stage
- **Location**: `tronko-assign.c` 
- **Discard Conditions**: 
  - No hits to reference sequences
  - Mapping quality below threshold
  - Other BWA-specific filtering

### 2. Alignment Quality Stage
- **Location**: `alignment.c` / `alignment_scoring.c`
- **Discard Conditions**:
  - Alignment score below minimum threshold
  - Alignment length below minimum threshold
  - Excessive gaps or mismatches
  - Other alignment quality metrics

### 3. Tree Placement Stage
- **Location**: `placement.c`
- **Discard Conditions**:
  - Ambiguous placement across multiple trees
  - Extremely low likelihood values
  - Inconsistent placements between paired-end reads
  - Other placement-specific criteria

### 4. LCA Determination Stage
- **Location**: `assignment.c`
- **Discard Conditions**:
  - Score below LCA cutoff at all relevant taxonomic levels
  - Assigned to excessively high taxonomic level
  - Other LCA-specific criteria

## Implementation Details

### 1. Core Logging Infrastructure

Create a new module (`discard_logger.c`/`discard_logger.h`) with the following functionality:

```c
/* Initialize the logging system */
void init_discard_logger(const char* log_file_path, int log_level);

/* Log a discarded read with detailed information */
void log_discarded_read(const char* read_name, const char* stage, const char* reason,
                        double actual_value, double threshold_value,
                        const char* additional_info);

/* Generate a summary of all discarded reads */
void summarize_discarded_reads(FILE* output_stream);

/* Clean up and close the logging system */
void close_discard_logger();
```

### 2. Command-Line Interface Extensions

Add new command-line options to `tronko-assign`:

```
-D, Enable discarded reads logging
-L [FILE], Path to write discarded reads log (default: tronko_discarded_reads.log)
-V [INT], Verbosity level for discarded reads log (0-3, default: 1)
```

### 3. Modification of Filtering Code

Insert logging calls at each point where reads are discarded. Example modifications:

#### In BWA Alignment Stage (`tronko-assign.c`):

```c
if (bwa_hit == NULL || bwa_hit->mapq < MIN_MAPQ) {
    // Current: Silently discard read
    // New: Log the discard event
    char reason[100];
    if (bwa_hit == NULL) {
        sprintf(reason, "No BWA hit found");
        log_discarded_read(read_name, "BWA_Alignment", reason, 0.0, 0.0, NULL);
    } else {
        sprintf(reason, "Mapping quality below threshold");
        log_discarded_read(read_name, "BWA_Alignment", reason, 
                         (double)bwa_hit->mapq, (double)MIN_MAPQ, NULL);
    }
    continue; // Proceed with discard
}
```

#### Similar modifications for other filtering stages

### 4. Log File Format

Generate a tab-delimited log file with the following columns:

```
ReadName   Stage   Reason   ActualValue   ThresholdValue   AdditionalInfo   Timestamp
```

Example entries:
```
Read1    BWA_Alignment    No_BWA_hit    0.0    0.0    NULL    2025-03-15T10:32:15
Read2    BWA_Alignment    Low_mapq    12.5    20.0    NULL    2025-03-15T10:32:15
Read3    Alignment_Quality    Low_score    -42.3    -30.0    alignment_length=124    2025-03-15T10:32:16
Read4    LCA_Determination    Below_cutoff    3.2    5.0    taxonomic_level=genus    2025-03-15T10:32:16
```

### 5. Summary Statistics

At the end of the log file, generate summary statistics:

```
# Summary Statistics
Total_Reads_Processed: 10000
Total_Reads_Discarded: 1500 (15.0%)
BWA_Alignment: 800 (8.0%)
Alignment_Quality: 450 (4.5%)
Tree_Placement: 150 (1.5%)
LCA_Determination: 100 (1.0%)
Reads_Successfully_Assigned: 8500 (85.0%)
```

## Code Changes Required

### New Files:
1. `discard_logger.c`
2. `discard_logger.h`

### Modified Files:
1. `tronko-assign.c` - Add command-line options, initialize logger, BWA filtering logging
2. `alignment.c` - Add alignment quality filtering logging
3. `placement.c` - Add tree placement filtering logging
4. `assignment.c` - Add LCA determination filtering logging
5. `options.c` - Add new options for discard logging
6. `options.h` - Add option declarations
7. `Makefile` - Add new source files

## Testing Strategy

1. **Unit Tests**:
   - Create tests to verify logger initialization and functionality
   - Test each logging call with different parameter combinations
   - Verify log file format and content

2. **Integration Tests**:
   - Process test datasets with known filtering characteristics
   - Verify that all discarded reads are correctly logged
   - Validate summary statistics against manual counts

3. **Validation Tests**:
   - Compare discarded read counts with expected outcomes
   - Verify that enabling logging doesn't affect assignment results
   - Test performance impact and optimize if necessary

## Performance Considerations

1. **File I/O Overhead**:
   - Use buffered I/O for logging to minimize performance impact
   - Consider periodic flushing rather than writing each entry immediately

2. **Memory Usage**:
   - Track counts in memory and only write summaries when needed
   - Avoid excessive string operations in logging code

3. **Scalability**:
   - Ensure logging system can handle very large datasets efficiently
   - Consider optional compression for log files

## Implementation Timeline

1. **Phase 1**: Core logging infrastructure (3 days)
   - Implement `discard_logger.c`/`discard_logger.h`
   - Add command-line options
   - Develop basic tests

2. **Phase 2**: BWA and alignment stage integration (2 days)
   - Modify `tronko-assign.c` for BWA filtering
   - Modify `alignment.c` for alignment quality filtering
   - Test with sample datasets

3. **Phase 3**: Tree placement and LCA stage integration (2 days)
   - Modify `placement.c` for tree placement filtering
   - Modify `assignment.c` for LCA determination filtering
   - Complete integration tests

4. **Phase 4**: Summary statistics and optimization (1 day)
   - Implement summary statistics generation
   - Optimize performance
   - Complete validation tests

## Benefits and Impact

This logging system will provide several key benefits:

1. **Troubleshooting**: Quickly identify why specific reads are being discarded
2. **Parameter Optimization**: Fine-tune filtering thresholds based on actual data
3. **Quality Control**: Monitor the distribution of discarded reads for anomalies
4. **Performance Analysis**: Identify bottlenecks in the assignment process
5. **Research Insights**: Gain deeper understanding of read characteristics that lead to poor assignments

## Future Extensions

1. **Visualization Tools**: Develop scripts to visualize discarded read patterns
2. **Parameter Recommendation**: Implement algorithms to suggest optimal parameters based on discard patterns
3. **Real-time Monitoring**: Add capabilities for real-time monitoring of discarded reads during processing
4. **Database Integration**: Allow logging to a database for more complex analysis