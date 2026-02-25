# tronko-assign Examples

This directory contains example scripts demonstrating various tronko-assign features.

## reproducibility_examples.sh

Demonstrates three approaches to reproducibility:

1. **Single-threaded mode** - 100% deterministic but slow
2. **Multi-threaded mode** - Fast but ~3% variance
3. **Variance testing** - Measure variance on your dataset

### Usage

1. Edit the script to set file paths:
   ```bash
   vim reproducibility_examples.sh
   # Update paths at the top:
   REF_TRKB="path/to/reference_tree.trkb"
   REF_FASTA="path/to/reference.fasta"
   FORWARD_READS="path/to/forward.fasta"
   REVERSE_READS="path/to/reverse.fasta"
   ```

2. Run the script:
   ```bash
   ./reproducibility_examples.sh
   ```

3. The script will run each example interactively, showing commands and waiting for confirmation.

### Requirements

- tronko-assign binary in parent directory
- Reference database (.trkb file from tronko-build)
- Reference sequences (.fasta file)
- Query reads (paired-end FASTA files)

### Expected Output

- `deterministic_results.txt` - Results from single-threaded mode
- `fast_results.txt` - Results from multi-threaded mode
- `/tmp/determinism_test/` - Variance test outputs

See the main README.md "Reproducibility and Variance" section for more details.
