---
date: 2026-01-01T12:00:00-08:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Streaming tronko-assign output to Parquet files"
tags: [research, codebase, parquet, output, tronko-assign, performance]
status: complete
last_updated: 2026-01-01
last_updated_by: Claude
---

# Research: Streaming tronko-assign Output to Parquet Files

**Date**: 2026-01-01
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

What would be involved in changing tronko-assign to stream output to a Parquet file instead of the current assignments.txt TSV format?

## Summary

Converting tronko-assign from TSV to Parquet output is **feasible but requires significant changes**:

1. **Library Selection**: Two options exist - Apache Arrow C++ (production-ready, requires C++17) or Carquet (pure C11, simpler integration)
2. **Architecture Change**: The current batch-buffer-then-write pattern is actually well-suited for Parquet's row-group model
3. **Threading**: Concurrent writes to a single Parquet file are NOT safe - would need per-thread files or single-writer pattern
4. **Build Complexity**: Arrow adds significant dependencies; Carquet is lighter but less mature
5. **Schema**: The current 7-column flat schema maps cleanly to Parquet primitive types

**Recommendation**: For minimal disruption, use Carquet (pure C) and write one Parquet file per thread, merging downstream. For production robustness, use Arrow C++ with a single-writer thread collecting results from workers.

## Detailed Findings

### Current Output Mechanism

#### File Operations Flow

The current output flow follows a batch-buffer-then-write pattern:

1. **File Opened** with header at start:
   - Single-end: `tronko-assign.c:1156-1158`
   - Paired-end: `tronko-assign.c:1414-1416`

   ```c
   FILE *results = fopen(opt.results_file,"w");
   fprintf(results,"Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");
   ```

2. **Per-batch Processing**:
   - Memory allocated for `taxonPath` arrays per thread (`tronko-assign.c:1254-1257`)
   - Threads populate their arrays via `runAssignmentOnChunk_WithBWA` (`tronko-assign.c:697, 711, 764`)
   - `pthread_join` barrier ensures all threads complete
   - Main thread writes all results sequentially (`tronko-assign.c:1275-1279`)

   ```c
   for (i=0; i<opt.number_of_cores; i++){
       for (j=0; j<(mstr[i].end-mstr[i].start); j++){
           fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
       }
   }
   ```

3. **File Closed** after all batches complete (`tronko-assign.c:1319` / `1527`)

#### Key Insight
There is **no mutex around file writes** because all writes happen from the main thread after `pthread_join`. This is important for Parquet integration.

### Current Output Schema

| Column | Name | C Type | Parquet Type |
|--------|------|--------|--------------|
| 1 | Readname | `char*` | BYTE_ARRAY (STRING) |
| 2 | Taxonomic_Path | `char*` | BYTE_ARRAY (STRING) |
| 3 | Score | `double` (or `float`) | DOUBLE (or FLOAT) |
| 4 | Forward_Mismatch | `double` | DOUBLE |
| 5 | Reverse_Mismatch | `double` | DOUBLE |
| 6 | Tree_Number | `int` | INT32 |
| 7 | Node_Number | `int` | INT32 |

The `type_of_PP` is `double` by default, `float` if `USE_FLOAT` is defined (`global.h:15-19`).

**Special Values for Taxonomic_Path**:
- `"unassigned"` - No BWA matches found
- `"unassigned Euk or Bac"` - Single match with no taxonomy
- `"unassigned Eukaryote or Bacteria"` - Conflicting multi-match taxonomy

### Parquet Library Options

#### Option 1: Apache Arrow C++ (Recommended for Production)

**Pros**:
- Battle-tested, used across industry
- Full feature support (compression, encryption, nested types)
- Excellent documentation and community support

**Cons**:
- Requires C++17 (would need mixed C/C++ compilation)
- Large dependency tree (Boost, Thrift, etc.)
- Complex build process

**Streaming Write Pattern**:
```cpp
#include "parquet/stream_writer.h"

parquet::StreamWriter os{
    parquet::ParquetFileWriter::Open(outfile, schema, props)};

// Per-row write (from batch buffer)
os << readname << taxon_path << score << fwd_mm << rev_mm
   << tree_num << node_num << parquet::EndRow;
```

**Batch Write Pattern** (better for current architecture):
```cpp
std::unique_ptr<parquet::arrow::FileWriter> writer;
// ... setup ...
for (auto batch : batches) {
    writer->WriteRecordBatch(*batch);  // One row group per batch
}
writer->Close();
```

#### Option 2: Carquet (Pure C Alternative)

**Pros**:
- Pure C11, easy integration with existing codebase
- Minimal dependencies (zstd, zlib - already linked)
- ~15k LOC, MIT licensed
- SIMD optimizations for x86/ARM

**Cons**:
- Newer, less battle-tested
- Flat schemas only (sufficient for this use case)
- No encryption support

**Write Pattern**:
```c
#include "carquet.h"

carquet_writer_t* writer = carquet_writer_create(
    "output.parquet", schema, &opts, &err);

// Write batch of values per column
carquet_writer_write_batch(writer, col_idx, values, count, NULL, NULL);

carquet_writer_close(writer);
```

### Architectural Changes Required

#### Change 1: Schema Definition

Need to define Parquet schema at startup:

```c
// With Carquet (pure C)
carquet_schema_t* create_assignment_schema() {
    carquet_schema_t* schema = carquet_schema_create();
    carquet_schema_add_field(schema, "Readname", CARQUET_TYPE_STRING, false);
    carquet_schema_add_field(schema, "Taxonomic_Path", CARQUET_TYPE_STRING, false);
    carquet_schema_add_field(schema, "Score", CARQUET_TYPE_DOUBLE, false);
    carquet_schema_add_field(schema, "Forward_Mismatch", CARQUET_TYPE_DOUBLE, false);
    carquet_schema_add_field(schema, "Reverse_Mismatch", CARQUET_TYPE_DOUBLE, false);
    carquet_schema_add_field(schema, "Tree_Number", CARQUET_TYPE_INT32, false);
    carquet_schema_add_field(schema, "Node_Number", CARQUET_TYPE_INT32, false);
    return schema;
}
```

#### Change 2: Data Structure Modification

Current: Results stored as pre-formatted strings in `taxonPath[]`
Required: Store structured data for columnar output

```c
// New structure for Parquet output
typedef struct assignmentResult {
    char* readname;
    char* taxonomic_path;
    double score;
    double forward_mismatch;
    double reverse_mismatch;
    int tree_number;
    int node_number;
} assignmentResult;

// In resultsStruct (global.h)
typedef struct resultsStruct {
    // ... existing fields ...
    assignmentResult* results;  // Replace taxonPath
} resultsStruct;
```

**Impact**: Requires modifying `tronko-assign.c:697, 711, 764` where `strcpy(results->taxonPath[iter], resultsPath)` is called.

#### Change 3: Batch Write Logic

Replace the current fprintf loop with columnar batch writes:

```c
// Current (tronko-assign.c:1275-1279)
for (i=0; i<opt.number_of_cores; i++){
    for (j=0; j<(mstr[i].end-mstr[i].start); j++){
        fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
    }
}

// New: Collect columns and write batch
for (i=0; i<opt.number_of_cores; i++){
    int batch_size = mstr[i].end - mstr[i].start;
    assignmentResult* batch = mstr[i].str->results;

    // Write each column
    write_string_column(writer, 0, batch->readnames, batch_size);
    write_string_column(writer, 1, batch->taxon_paths, batch_size);
    write_double_column(writer, 2, batch->scores, batch_size);
    // ... etc
}
```

#### Change 4: Threading Strategy

**Critical**: Concurrent writes to a single Parquet file are NOT safe.

**Options**:

1. **Single Writer (Current Pattern)** - Works well
   - Keep current pattern: threads buffer, main thread writes
   - Each batch becomes a Parquet row group
   - Clean, simple, minimal code changes

2. **Per-Thread Files** - Maximum parallelism
   - Each thread writes to its own `.parquet` file
   - Downstream tools can read multiple files as one dataset
   - Requires output naming convention (e.g., `output_0.parquet`, `output_1.parquet`)

3. **Queue-Based Writer** - Production pattern
   - Worker threads push results to a lock-free queue
   - Dedicated writer thread consumes and writes
   - Most complex but scales best

**Recommendation**: Start with option 1 (current pattern), which requires minimal threading changes.

#### Change 5: Build System Updates

**For Carquet** (simpler):
```makefile
# Add to Makefile
CARQUET_DIR = ../carquet
CFLAGS += -I$(CARQUET_DIR)/include
LDFLAGS += -L$(CARQUET_DIR)/lib -lcarquet -lzstd
```

**For Arrow C++** (more complex):
```makefile
# Need C++ compilation for Arrow wrapper
CXX = g++
CXXFLAGS = -std=c++17 -O3

# Link Arrow libraries
LDFLAGS += -larrow -lparquet -lpthread
```

Would require either:
- Writing a C wrapper around Arrow C++ APIs
- Converting relevant output code to C++
- Using `extern "C"` bridges

### Compression Considerations

Parquet supports multiple compression codecs:

| Codec | Speed | Ratio | Notes |
|-------|-------|-------|-------|
| None | Fastest | 1:1 | Good for debugging |
| Snappy | Fast | ~2:1 | Default for many systems |
| LZ4 | Fast | ~2:1 | Best speed/ratio balance |
| ZSTD | Medium | ~3:1 | Best ratio, good for archival |
| GZIP | Slow | ~3:1 | Wide compatibility |

For metabarcoding output with many repeated taxonomy strings, **ZSTD** or **dictionary encoding** would provide excellent compression (~10x over raw TSV).

### Migration Path

#### Phase 1: Add Carquet as Optional Output
1. Add Carquet as git submodule or vendored dependency
2. Add `-P` flag to enable Parquet output
3. Modify result structures to store typed data
4. Keep TSV as default for backward compatibility

#### Phase 2: Optimize Data Structures
1. Refactor `taxonPath` to structured `assignmentResult`
2. Implement columnar batch writing
3. Add compression configuration

#### Phase 3: Consider Arrow Migration (Optional)
1. If more features needed (nested types, encryption)
2. Create C wrapper around Arrow C++ APIs
3. Update build system for mixed C/C++ compilation

## Code References

- `tronko-assign/tronko-assign.c:1156-1158` - Output file open and header write (single-end)
- `tronko-assign/tronko-assign.c:1414-1416` - Output file open and header write (paired-end)
- `tronko-assign/tronko-assign.c:1275-1279` - Batch result writing (single-end)
- `tronko-assign/tronko-assign.c:1504-1508` - Batch result writing (paired-end)
- `tronko-assign/tronko-assign.c:697, 711, 764` - Result string population in worker
- `tronko-assign/global.h:115-133` - `resultsStruct` definition
- `tronko-assign/global.h:15-19` - `type_of_PP` conditional typedef
- `tronko-assign/Makefile` - Current build configuration

## Architecture Insights

1. **Batch-Buffer Pattern**: The existing architecture buffers results in memory before writing, which aligns well with Parquet's row-group model. Each batch could naturally become one row group.

2. **Thread Isolation**: Each thread has isolated result storage (`mstr[i].str`), avoiding contention. This pattern should be preserved.

3. **No Mid-Stream Writes**: Unlike true streaming, Parquet requires the footer (with all metadata) to be written last. The current pattern of "process all batches, then close file" works well.

4. **Column-Store Fit**: The flat, fixed-schema output maps perfectly to Parquet's columnar format. Dictionary encoding on `Taxonomic_Path` could dramatically reduce file size.

## Open Questions

1. **Backward Compatibility**: Should Parquet be an additional output format (`-P` flag) or replace TSV entirely?

2. **Compression Trade-off**: What compression level is acceptable given the speed requirements of tronko-assign?

3. **Multi-File Output**: Would downstream consumers prefer one large Parquet file or partitioned files (one per batch or per thread)?

4. **Library Choice**: Is C++ acceptable for Arrow integration, or must the codebase remain pure C?

5. **Schema Evolution**: If output schema changes in future versions, how should compatibility be handled?
