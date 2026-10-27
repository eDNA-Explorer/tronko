# AncestralClust Golden Test Data

## How golden output was generated

```bash
# From the tronko-fork root directory, using single-threaded mode for determinism (srand(42)):
./ancestralclust/ancestralclust \
  -f -u \
  -i tronko-build/example_datasets/single_tree/Charadriiformes.fasta \
  -b 50 -r 100 -p 10 \
  -c 1 \
  -d ancestralclust/test_data/golden_output
```

Input: Charadriiformes.fasta (1466 sequences, COI barcode region)
Parameters: 50 initial bins, 100 seed sequences, 10 descendant threshold, NW alignment (-u), single-threaded (-c 1)
Output: 13 clusters (0.fasta through 12.fasta) + output.clstr

## Verification

Run `test_consistency.sh` from the tronko-fork root to verify optimized builds produce identical output.

## Regenerating golden output

If the clustering algorithm itself changes (not just performance optimizations), regenerate:
```bash
rm -rf ancestralclust/test_data/golden_output/*
# Run the command above, then commit the new golden output
```
