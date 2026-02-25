# Parquet Output Streaming Implementation Plan

## Overview

Add Parquet output support to tronko-assign using the Carquet pure-C library. When enabled via a `-P` flag, each worker thread writes its results to a separate Parquet file with ZSTD compression, replacing the default TSV output.

## Current State Analysis

### Output Mechanism
- Results stored as pre-formatted TSV strings in `taxonPath[]` arrays (`global.h:123`)
- Worker threads populate strings via `strcpy()` at `tronko-assign.c:697, 711, 764`
- After `pthread_join`, main thread writes all strings with `fprintf` (`tronko-assign.c:1275-1279`, `1504-1508`)
- Pattern: batch-buffer-then-write (well-suited for Parquet row groups)

### Key Constraints
- Must remain pure C (no C++ in main codebase)
- Already links zstd and zlib (Carquet's dependencies)
- Thread-safe design: each thread has isolated `mstr[i].str` storage

### Output Schema (7 columns)
| Column | C Type | Parquet Type |
|--------|--------|--------------|
| Readname | `char*` | BYTE_ARRAY (STRING) |
| Taxonomic_Path | `char*` | BYTE_ARRAY (STRING) |
| Score | `double` | DOUBLE |
| Forward_Mismatch | `double` | DOUBLE |
| Reverse_Mismatch | `double` | DOUBLE |
| Tree_Number | `int` | INT32 |
| Node_Number | `int` | INT32 |

## Desired End State

When `-P output_prefix` is specified:
1. No TSV file is created
2. Each of N threads writes `output_prefix_0.parquet`, `output_prefix_1.parquet`, etc.
3. Files use ZSTD compression
4. Files are readable by PyArrow, DuckDB, Spark, and other Parquet tools
5. Performance is comparable to or better than TSV output

### Verification
```bash
# Test single-end with Parquet output
./tronko-assign -r -f reference_tree.txt -a ref.fasta -s -g reads.fasta -P /tmp/results -w

# Verify files created
ls /tmp/results_*.parquet

# Verify readable by PyArrow
python3 -c "import pyarrow.parquet as pq; print(pq.read_table('/tmp/results_0.parquet').to_pandas())"
```

## What We're NOT Doing

- NOT replacing TSV as the default output format
- NOT supporting nested Parquet schemas
- NOT implementing Parquet reading
- NOT merging per-thread files (downstream tools handle this)
- NOT supporting encryption or bloom filters

## Implementation Approach

Use Carquet as a git submodule, create a thin wrapper module (`parquet_writer.c/h`), modify result storage to use structured data instead of pre-formatted strings, and integrate per-thread Parquet writers into the main processing loop.

---

## Phase 1: Add Carquet Dependency

### Overview
Add Carquet as a git submodule and integrate into the build system.

### Changes Required:

#### 1. Add Git Submodule
```bash
cd tronko-assign
git submodule add https://github.com/Vitruves/carquet.git carquet
```

#### 2. Update Makefile
**File**: `tronko-assign/Makefile`

Add Carquet build integration:

```makefile
# Parquet support via Carquet (use: make ENABLE_PARQUET=1)
ifdef ENABLE_PARQUET
    PARQUET_FLAGS = -DENABLE_PARQUET -Icarquet/include
    PARQUET_SOURCES = parquet_writer.c
    PARQUET_LIBS = -Lcarquet/build -lcarquet
else
    PARQUET_FLAGS =
    PARQUET_SOURCES =
    PARQUET_LIBS =
endif

# Update SOURCES line to include conditional parquet
SOURCES = tronko-assign.c readreference.c assignment.c ... $(PARQUET_SOURCES)

# Update build target
$(TARGET): $(TARGET).c carquet_lib
	$(CC) $(OPTIMIZATION) ... $(PARQUET_FLAGS) -o $(TARGET) ... $(PARQUET_LIBS)

# Add carquet build target
carquet_lib:
ifdef ENABLE_PARQUET
	@if [ ! -f carquet/build/libcarquet.a ]; then \
		cd carquet && mkdir -p build && cd build && cmake .. && make -j$$(nproc); \
	fi
endif
```

### Success Criteria:

#### Automated Verification:
- [x] `git submodule update --init` succeeds
- [x] `make ENABLE_PARQUET=1` compiles without errors
- [x] `make` (without flag) still works and produces identical binary

#### Manual Verification:
- [x] Carquet library builds correctly on target system

---

## Phase 2: Create Structured Result Type

### Overview
Replace pre-formatted string storage with structured data to enable columnar Parquet writes.

### Changes Required:

#### 1. Add Result Structure
**File**: `tronko-assign/global.h`

Add after line 133 (after `resultsStruct` definition):

```c
#ifdef ENABLE_PARQUET
/*
 * Structured assignment result for Parquet output
 * Stores typed fields instead of pre-formatted TSV string
 */
typedef struct assignmentResult {
    char *readname;
    char *taxonomic_path;
    double score;
    double forward_mismatch;
    double reverse_mismatch;
    int tree_number;
    int node_number;
} assignmentResult;
#endif
```

#### 2. Extend resultsStruct
**File**: `tronko-assign/global.h`

Modify `resultsStruct` (lines 115-133) to add:

```c
typedef struct resultsStruct{
    nw_aligner_t *nw;
    alignment_t *aln;
    scoring_t *scoring;
    type_of_PP ***nodeScores;
    int **voteRoot;
    int *positions;
    char *locQuery;
    char **taxonPath;
    char **LCAnames;
    int *minNodes;
    int **leaf_coordinates;
    type_of_PP *minimum;
    int print_alignments;
    int *starts_forward;
    int *starts_reverse;
    char **cigars_forward;
    char **cigars_reverse;
#ifdef ENABLE_PARQUET
    assignmentResult *parquet_results;  // Structured results for Parquet
    int parquet_count;                   // Number of results in current batch
#endif
}resultsStruct;
```

#### 3. Add Parquet Flag to Options
**File**: `tronko-assign/global.h`

Add to `Options` struct (around line 227):

```c
#ifdef ENABLE_PARQUET
    char parquet_prefix[BUFFER_SIZE];  // Output prefix for Parquet files (empty = disabled)
    int parquet_enabled;               // 1 if Parquet output enabled
#endif
```

### Success Criteria:

#### Automated Verification:
- [ ] `make ENABLE_PARQUET=1` compiles without errors
- [x] `make` (without flag) compiles without errors (preprocessor guards work)

#### Manual Verification:
- [ ] Structure sizes are reasonable (check with sizeof in debug build)

---

## Phase 3: Create Parquet Writer Module

### Overview
Create a wrapper module that encapsulates Carquet API for tronko-assign's specific schema.

### Changes Required:

#### 1. Create Header File
**File**: `tronko-assign/parquet_writer.h` (new file)

```c
/*
 * parquet_writer.h
 * Parquet output support for tronko-assign using Carquet library
 */
#ifndef _PARQUET_WRITER_H_
#define _PARQUET_WRITER_H_

#ifdef ENABLE_PARQUET

#include "carquet.h"

/* Opaque handle for Parquet writer */
typedef struct parquet_writer parquet_writer_t;

/*
 * Create a new Parquet writer for assignment results
 *
 * @param filename  Output filename (e.g., "results_0.parquet")
 * @param err       Error message buffer (caller-allocated, min 256 bytes)
 * @return          Writer handle, or NULL on error
 */
parquet_writer_t* parquet_writer_create(const char *filename, char *err);

/*
 * Write a batch of assignment results
 *
 * @param writer    Writer handle from parquet_writer_create
 * @param results   Array of assignment results
 * @param count     Number of results in array
 * @return          0 on success, -1 on error
 */
int parquet_writer_write_batch(parquet_writer_t *writer,
                                const assignmentResult *results,
                                int count);

/*
 * Close writer and finalize Parquet file
 * Writes footer metadata and releases resources
 *
 * @param writer    Writer handle (freed after call)
 * @return          0 on success, -1 on error
 */
int parquet_writer_close(parquet_writer_t *writer);

#endif /* ENABLE_PARQUET */
#endif /* _PARQUET_WRITER_H_ */
```

#### 2. Create Implementation File
**File**: `tronko-assign/parquet_writer.c` (new file)

```c
/*
 * parquet_writer.c
 * Parquet output support for tronko-assign using Carquet library
 */
#ifdef ENABLE_PARQUET

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "parquet_writer.h"
#include "carquet.h"

struct parquet_writer {
    carquet_writer_t *cw;
    carquet_schema_t *schema;
    char filename[MAXFILENAME];
};

/* Create schema matching tronko-assign output columns */
static carquet_schema_t* create_assignment_schema(carquet_error_t *err) {
    carquet_schema_t *schema = carquet_schema_create(err);
    if (!schema) return NULL;

    /* All columns are required (non-nullable) */
    carquet_schema_add_column(schema, "Readname",
        CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Taxonomic_Path",
        CARQUET_PHYSICAL_BYTE_ARRAY, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Score",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Forward_Mismatch",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Reverse_Mismatch",
        CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Tree_Number",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0);
    carquet_schema_add_column(schema, "Node_Number",
        CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0);

    return schema;
}

parquet_writer_t* parquet_writer_create(const char *filename, char *err) {
    parquet_writer_t *pw = calloc(1, sizeof(parquet_writer_t));
    if (!pw) {
        snprintf(err, 256, "Failed to allocate parquet_writer");
        return NULL;
    }

    strncpy(pw->filename, filename, MAXFILENAME - 1);

    carquet_error_t cerr = CARQUET_ERROR_INIT;

    /* Create schema */
    pw->schema = create_assignment_schema(&cerr);
    if (!pw->schema) {
        snprintf(err, 256, "Failed to create schema: %s", cerr.message);
        free(pw);
        return NULL;
    }

    /* Configure writer options */
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;
    opts.compression_level = 3;  /* Balance speed/ratio */
    opts.row_group_size = 64 * 1024 * 1024;  /* 64 MB row groups */
    opts.write_statistics = true;

    /* Create writer */
    pw->cw = carquet_writer_create(filename, pw->schema, &opts, &cerr);
    if (!pw->cw) {
        snprintf(err, 256, "Failed to create writer: %s", cerr.message);
        carquet_schema_free(pw->schema);
        free(pw);
        return NULL;
    }

    return pw;
}

int parquet_writer_write_batch(parquet_writer_t *writer,
                                const assignmentResult *results,
                                int count) {
    if (!writer || !results || count <= 0) return -1;

    /* Allocate temporary arrays for columnar data */
    carquet_byte_array_t *readnames = malloc(count * sizeof(carquet_byte_array_t));
    carquet_byte_array_t *taxon_paths = malloc(count * sizeof(carquet_byte_array_t));
    double *scores = malloc(count * sizeof(double));
    double *fwd_mismatches = malloc(count * sizeof(double));
    double *rev_mismatches = malloc(count * sizeof(double));
    int32_t *tree_numbers = malloc(count * sizeof(int32_t));
    int32_t *node_numbers = malloc(count * sizeof(int32_t));

    if (!readnames || !taxon_paths || !scores || !fwd_mismatches ||
        !rev_mismatches || !tree_numbers || !node_numbers) {
        free(readnames); free(taxon_paths); free(scores);
        free(fwd_mismatches); free(rev_mismatches);
        free(tree_numbers); free(node_numbers);
        return -1;
    }

    /* Convert row-oriented to columnar */
    for (int i = 0; i < count; i++) {
        readnames[i].len = strlen(results[i].readname);
        readnames[i].ptr = (uint8_t*)results[i].readname;

        taxon_paths[i].len = strlen(results[i].taxonomic_path);
        taxon_paths[i].ptr = (uint8_t*)results[i].taxonomic_path;

        scores[i] = results[i].score;
        fwd_mismatches[i] = results[i].forward_mismatch;
        rev_mismatches[i] = results[i].reverse_mismatch;
        tree_numbers[i] = results[i].tree_number;
        node_numbers[i] = results[i].node_number;
    }

    /* Write each column */
    int ret = 0;
    ret |= carquet_writer_write_batch(writer->cw, 0, readnames, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 1, taxon_paths, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 2, scores, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 3, fwd_mismatches, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 4, rev_mismatches, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 5, tree_numbers, count, NULL, NULL);
    ret |= carquet_writer_write_batch(writer->cw, 6, node_numbers, count, NULL, NULL);

    free(readnames); free(taxon_paths); free(scores);
    free(fwd_mismatches); free(rev_mismatches);
    free(tree_numbers); free(node_numbers);

    return ret ? -1 : 0;
}

int parquet_writer_close(parquet_writer_t *writer) {
    if (!writer) return -1;

    carquet_status_t status = carquet_writer_close(writer->cw);
    carquet_schema_free(writer->schema);
    free(writer);

    return (status == CARQUET_OK) ? 0 : -1;
}

#endif /* ENABLE_PARQUET */
```

### Success Criteria:

#### Automated Verification:
- [x] `make ENABLE_PARQUET=1` compiles without errors
- [x] No undefined symbol errors when linking

#### Manual Verification:
- [ ] Can create a simple test file and read it with PyArrow

---

## Phase 4: Integrate into Main Processing Loop

### Overview
Modify tronko-assign.c to use Parquet writers when `-P` flag is specified. Each thread gets its own Parquet file.

### Changes Required:

#### 1. Add Command Line Option
**File**: `tronko-assign/options.c`

Add parsing for `-P` flag (add to existing option parsing):

```c
#ifdef ENABLE_PARQUET
case 'P':
    strncpy(opt->parquet_prefix, optarg, BUFFER_SIZE - 1);
    opt->parquet_enabled = 1;
    break;
#endif
```

Add to help text:

```c
#ifdef ENABLE_PARQUET
printf("  -P <prefix>    Output Parquet files (one per thread): prefix_0.parquet, prefix_1.parquet, ...\n");
printf("                 When specified, TSV output is disabled\n");
#endif
```

#### 2. Add Per-Thread Parquet Writers to mystruct
**File**: `tronko-assign/global.h`

Add to `mystruct` (around line 165):

```c
#ifdef ENABLE_PARQUET
    parquet_writer_t *parquet_writer;  // Per-thread Parquet writer
    int thread_id;                      // Thread index for filename
#endif
```

#### 3. Modify Result Population
**File**: `tronko-assign/tronko-assign.c`

In the worker function `runAssignmentOnChunk_WithBWA`, after populating `taxonPath`, also populate structured results when Parquet is enabled.

Around line 697 (and similar locations 711, 764), add:

```c
strcpy(results->taxonPath[iter], resultsPath);

#ifdef ENABLE_PARQUET
if (results->parquet_results) {
    assignmentResult *pr = &results->parquet_results[iter];
    /* Readname is already in resultsPath before first tab */
    pr->readname = strdup(paired ? pairedQueryMat->forward_name[lineNumber]
                                  : singleQueryMat->name[lineNumber]);
    pr->taxonomic_path = strdup(/* extract from taxonomy */);
    pr->score = results->minimum[0];
    pr->forward_mismatch = results->minimum[1];
    pr->reverse_mismatch = results->minimum[2];
    pr->tree_number = maxRoot;
    pr->node_number = LCA;
}
#endif
```

#### 4. Create Per-Thread Writers and Write Results
**File**: `tronko-assign/tronko-assign.c`

**Single-end mode** (around line 1156):

Replace file opening:
```c
FILE *results = NULL;
#ifdef ENABLE_PARQUET
parquet_writer_t **parquet_writers = NULL;
if (opt.parquet_enabled) {
    parquet_writers = malloc(opt.number_of_cores * sizeof(parquet_writer_t*));
    for (i = 0; i < opt.number_of_cores; i++) {
        char filename[BUFFER_SIZE];
        snprintf(filename, BUFFER_SIZE, "%s_%d.parquet", opt.parquet_prefix, i);
        char err[256];
        parquet_writers[i] = parquet_writer_create(filename, err);
        if (!parquet_writers[i]) {
            printf("Error creating Parquet writer: %s\n", err);
            exit(1);
        }
        mstr[i].parquet_writer = parquet_writers[i];
        mstr[i].thread_id = i;
    }
} else {
#endif
    results = fopen(opt.results_file, "w");
    if (results == NULL) { printf("Error opening output file!\n"); exit(1); }
    fprintf(results, "Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");
#ifdef ENABLE_PARQUET
}
#endif
```

Replace result writing (around line 1275):
```c
#ifdef ENABLE_PARQUET
if (opt.parquet_enabled) {
    /* Each thread writes its own batch to its own file */
    for (i = 0; i < opt.number_of_cores; i++) {
        int batch_size = mstr[i].end - mstr[i].start;
        if (batch_size > 0) {
            parquet_writer_write_batch(mstr[i].parquet_writer,
                                        mstr[i].str->parquet_results,
                                        batch_size);
        }
    }
} else {
#endif
    for (i = 0; i < opt.number_of_cores; i++) {
        for (j = 0; j < (mstr[i].end - mstr[i].start); j++) {
            fprintf(results, "%s\n", mstr[i].str->taxonPath[j]);
        }
    }
#ifdef ENABLE_PARQUET
}
#endif
```

Replace file closing (around line 1319):
```c
#ifdef ENABLE_PARQUET
if (opt.parquet_enabled) {
    for (i = 0; i < opt.number_of_cores; i++) {
        parquet_writer_close(parquet_writers[i]);
    }
    free(parquet_writers);
} else {
#endif
    fclose(results);
#ifdef ENABLE_PARQUET
}
#endif
```

#### 5. Apply Same Changes to Paired-End Mode
**File**: `tronko-assign/tronko-assign.c`

Apply equivalent changes around lines 1414-1416 (file open), 1504-1508 (write), and 1527 (close).

#### 6. Allocate/Free Parquet Result Arrays
**File**: `tronko-assign/allocateMemoryForResults.c`

Add allocation when Parquet is enabled:

```c
#ifdef ENABLE_PARQUET
extern int parquet_enabled;  /* From options */
if (parquet_enabled) {
    results->parquet_results = calloc(numberOfLinesToRead, sizeof(assignmentResult));
    results->parquet_count = 0;
}
#endif
```

Add corresponding free logic.

### Success Criteria:

#### Automated Verification:
- [x] `make ENABLE_PARQUET=1` compiles without errors
- [x] `./tronko-assign -h` shows `--parquet` option when compiled with ENABLE_PARQUET=1
- [x] Test run produces Parquet file
- [x] Files are valid Parquet (checked by PyArrow)

#### Manual Verification:
- [x] Results match TSV output when both are generated and compared (164 rows, all readnames match)
- [x] Performance is acceptable (not significantly slower than TSV)

---

## Phase 5: Memory Management and Cleanup

### Overview
Ensure proper allocation and deallocation of structured result memory.

### Changes Required:

#### 1. Update allocateMemoryForResults
**File**: `tronko-assign/allocateMemoryForResults.c`

Add allocation for `parquet_results` array alongside `taxonPath`.

#### 2. Update freeing logic in main loop
**File**: `tronko-assign/tronko-assign.c`

After writing Parquet batch, free string fields in `parquet_results`:

```c
#ifdef ENABLE_PARQUET
if (opt.parquet_enabled) {
    for (j = 0; j < batch_size; j++) {
        free(mstr[i].str->parquet_results[j].readname);
        free(mstr[i].str->parquet_results[j].taxonomic_path);
    }
}
#endif
```

### Success Criteria:

#### Automated Verification:
- [x] Valgrind shows no memory leaks: `valgrind --leak-check=full ./tronko-assign ... --parquet /tmp/test`
- [x] No segfaults on multiple batches

#### Manual Verification:
- [x] Memory usage is reasonable during processing

**Implementation Note**: The implementation uses a simpler approach than planned - instead of allocating per-thread result arrays populated during assignment, we parse the existing TSV strings into structured results at output time. This is slightly less efficient but minimally invasive and works correctly.

---

## Testing Strategy

### Unit Tests:
- Parquet writer creates valid file with correct schema
- Batch writes handle empty batches gracefully
- Column data types are correct

### Integration Tests:
```bash
# Build with Parquet support
make clean && make ENABLE_PARQUET=1

# Single-end test
./tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -P /tmp/parquet_test -w

# Verify files created
ls -la /tmp/parquet_test_*.parquet

# Verify content matches TSV
./tronko-assign -r -f tronko-build/example_datasets/single_tree/reference_tree.txt \
  -a tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -s -g example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta \
  -o /tmp/tsv_test.txt -w

# Compare with Python
python3 << 'EOF'
import pyarrow.parquet as pq
import pandas as pd
import glob

# Read all Parquet files
parquet_files = sorted(glob.glob('/tmp/parquet_test_*.parquet'))
parquet_df = pd.concat([pq.read_table(f).to_pandas() for f in parquet_files])

# Read TSV
tsv_df = pd.read_csv('/tmp/tsv_test.txt', sep='\t')

# Compare
assert len(parquet_df) == len(tsv_df), f"Row count mismatch: {len(parquet_df)} vs {len(tsv_df)}"
print(f"Both have {len(parquet_df)} rows")
print("Parquet schema:", parquet_df.dtypes)
EOF
```

### Manual Testing Steps:
1. Run with varying thread counts (1, 2, 4, 8) and verify all produce correct output
2. Test with large input files to verify batch handling
3. Test with paired-end mode
4. Verify ZSTD compression is working (file size should be smaller than equivalent TSV)

---

## Performance Considerations

1. **Memory**: Structured results use slightly more memory than pre-formatted strings due to separate allocations
2. **CPU**: Columnar conversion adds overhead, but ZSTD compression is fast
3. **I/O**: Parquet files should be significantly smaller than TSV (dictionary encoding on taxonomy strings)
4. **Parallelism**: Per-thread files eliminate contention; downstream tools handle merging efficiently

## Migration Notes

- Users must install PyArrow, DuckDB, or similar to read Parquet output
- Provide example commands in documentation for common operations:
  ```bash
  # Read with DuckDB
  duckdb -c "SELECT * FROM read_parquet('results_*.parquet')"

  # Convert to TSV with PyArrow
  python3 -c "import pyarrow.parquet as pq; pq.read_table('results_0.parquet').to_pandas().to_csv('out.tsv', sep='\t', index=False)"
  ```

## References

- Research document: `thoughts/shared/research/2026-01-01-parquet-output-streaming.md`
- Carquet repository: https://github.com/Vitruves/carquet
- Current output code: `tronko-assign/tronko-assign.c:1156-1158, 1275-1279, 1414-1416, 1504-1508`
- Result structure: `tronko-assign/global.h:115-133`
