---
date: 2025-12-29T12:00:00-08:00
researcher: Claude
git_commit: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
branch: main
repository: tronko
topic: "How tronko-assign processes FASTA files and streams input/output"
tags: [research, codebase, tronko-assign, streaming, io, threading]
status: complete
last_updated: 2025-12-29
last_updated_by: Claude
---

# Research: How tronko-assign Processes FASTA Files and Handles I/O

**Date**: 2025-12-29T12:00:00-08:00
**Researcher**: Claude
**Git Commit**: 19a83663c009bdbd2e5207c5e5234d62c80da8f9
**Branch**: main
**Repository**: tronko

## Research Question

Once tronko has loaded the reference library -- how does it actually perform assignments from the FASTA files? Does it stream the FASTAs? Does it stream its output as well?

## Summary

**tronko-assign uses a "chunked streaming" architecture for both input and output:**

1. **Reference Library**: Loaded entirely into memory at startup (NOT streamed)
2. **FASTA Input**: Processed in configurable batches (default 50,000 lines) - a hybrid between pure streaming and full batch loading
3. **Output**: Written after each batch completes (batched streaming) - results accumulated in memory during parallel processing, then flushed sequentially

The key insight is that tronko-assign is designed for memory efficiency with large datasets while still leveraging parallelism. It doesn't load all query sequences at once, nor does it process them one-by-one. Instead, it uses a batch-parallel model.

## Detailed Findings

### Reference Library Loading (Fully In-Memory)

The reference database (`reference_tree.txt`) is loaded completely into memory before any query processing begins:

- **Entry point**: `tronko-assign.c:780-788` opens the file and calls `readReferenceTree()`
- **Implementation**: `readreference.c:311-433` reads the entire file using `gzgets()` in a loop
- **Data structures populated**:
  - `treeArr` - Phylogenetic tree nodes with posterior probabilities
  - `taxonomyArr` - 4D array `[tree][species][level][char]` for taxonomy
  - `numbaseArr`, `rootArr`, `numspecArr` - Per-tree metadata

After loading, posterior probabilities are log-transformed at `tronko-assign.c:825` via `store_PPs_Arr()` for efficient scoring.

### FASTA Input Processing (Chunked Streaming)

Query FASTA/FASTQ files are processed in **configurable batches**, not fully loaded or purely streamed:

#### Two-Pass Reading Strategy
1. **First pass** (`readreference.c:454-477`): `find_specs_for_reads()` scans entire file to find maximum sequence and name lengths
2. **Memory allocation**: Arrays pre-allocated based on max lengths
3. **Second pass**: File re-opened and read in batches

#### Batch Reading Implementation
- **Batch size**: Configurable via `-L` flag, default 50,000 lines (`tronko-assign.c:763`)
- **Reading function**: `readInXNumberOfLines()` at `readreference.c:144-300`
- **Storage**: Global `singleQueryMat` or `pairedQueryMat` structures (`global.h:58-68`)

```c
// Main loop structure (tronko-assign.c:935-981)
while (1) {
    returnLineNumber = readInXNumberOfLines(numberOfLinesToRead/2, seqinfile, ...);
    if (returnLineNumber == 0) break;  // EOF

    // Divide work among threads
    // Spawn threads
    // Wait for threads
    // Write results
    // Free batch memory
}
```

#### Compression Support
Uses `gzopen()`/`gzgets()` for transparent gzip support - same code handles both compressed and uncompressed files.

### Assignment Algorithm Flow

For each batch of sequences:

1. **BWA-MEM Alignment** (`tronko-assign.c:271`, `bwa_source_files/fastmap.c:572`)
   - Aligns query sequences against leaf node sequences
   - Finds candidate phylogenetic placements
   - Results stored as `(tree, node)` coordinates

2. **Fine-Grained Alignment** (`placement.c:60-930`)
   - Extracts leaf sequence from posterior probabilities (`getSequenceinRoot.c:138-167`)
   - Performs WFA or Needleman-Wunsch alignment (`placement.c:115-125`)
   - Maps query positions to reference MSA positions

3. **Score Calculation** (`assignment.c:24-65`, `assignment.c:143-210`)
   - Traverses tree recursively from matched leaves
   - Calculates likelihood scores using log-transformed posterior probabilities
   - Each node gets a score based on sequence alignment

4. **LCA Computation** (`placement.c:862-926`, `tronko-assign.c:508-554`)
   - Finds all nodes within confidence interval of max score
   - Computes Lowest Common Ancestor of high-scoring nodes
   - Cross-tree reconciliation for multi-tree databases

5. **Taxonomy Assignment** (`tronko-assign.c:560-684`)
   - Extracts taxonomic path from LCA node's `taxIndex`
   - Builds output string with scores

### Output Handling (Batched Streaming)

Output uses a **batched streaming pattern**:

1. **File opened once** at startup (`tronko-assign.c:899`)
2. **Header written immediately** (`tronko-assign.c:901`)
3. **Per-batch accumulation**: Results stored in thread-local `taxonPath` arrays
4. **Sequential write after batch** (`tronko-assign.c:970-974`):
   ```c
   for (i=0; i<opt.number_of_cores; i++) {
       for (j=0; j<(mstr[i].end-mstr[i].start); j++) {
           fprintf(results, "%s\n", mstr[i].str->taxonPath[j]);
       }
   }
   ```
5. **No explicit flushing** - relies on C library buffering
6. **File closed at end** (`tronko-assign.c:982`)

### Threading Model

Work is parallelized within each batch:

- **Thread creation**: `pthread_create()` at `tronko-assign.c:964-966`
- **Work division**: Sequences divided evenly among threads (`tronko-assign.c:944-958`)
- **Thread-local storage**: Each thread writes to its own `mstr[i].str->taxonPath` array
- **Synchronization**: `pthread_join()` waits for all threads before output
- **No locks needed**: Data race prevention via non-overlapping access patterns

## Code References

| Component | File | Lines | Description |
|-----------|------|-------|-------------|
| Reference loading | `tronko-assign.c` | 780-788 | Opens and reads reference file |
| Reference parsing | `readreference.c` | 311-433 | `readReferenceTree()` implementation |
| Query batch reading | `readreference.c` | 144-300 | `readInXNumberOfLines()` for FASTA |
| FASTQ reading | `readreference.c` | 2-141 | `readInXNumberOfLines_fastq()` |
| Main processing loop | `tronko-assign.c` | 935-981 | Single-end batch loop |
| Thread worker | `tronko-assign.c` | 174-742 | `runAssignmentOnChunk_WithBWA()` |
| BWA alignment | `bwa_source_files/fastmap.c` | 572+ | `main_mem()` function |
| Phylogenetic placement | `placement.c` | 60-930 | `place_paired()` function |
| Score calculation | `assignment.c` | 143-210 | `getscore_Arr()` function |
| Output writing | `tronko-assign.c` | 970-974 | Post-batch fprintf loop |

## Architecture Insights

### Why Chunked Streaming?

1. **Memory efficiency**: Doesn't require loading all query sequences into memory
2. **Parallelism**: Can still use multiple threads within each batch
3. **Deterministic output**: Results written in consistent order (by thread, then by sequence)
4. **Progress visibility**: Could theoretically track progress per-batch

### Key Design Decisions

1. **Two-pass file reading**: Trades I/O time for optimal memory allocation
2. **Global query storage**: Avoids per-thread memory allocation overhead
3. **Thread-local results**: Eliminates need for synchronization during processing
4. **Sequential output after parallel processing**: Ensures deterministic ordering

### Configuration Options

| Flag | Default | Purpose |
|------|---------|---------|
| `-L` | 50000 | Batch size (lines to read per iteration) |
| `-C` | 1 | Number of processing threads |
| `-q` | off | FASTQ mode (adjusts line counting: /4 vs /2) |

## Open Questions

1. **Memory scaling**: What's the relationship between batch size and peak memory usage?
2. **Optimal batch size**: Has the default 50,000 been tuned for specific hardware?
3. **Output buffering**: Would explicit `fflush()` after batches improve data durability?
4. **Two-pass trade-off**: Could max lengths be estimated or user-provided to avoid the first pass?
