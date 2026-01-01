# Tronko Experiments Log

This document records benchmark results and experiments conducted to evaluate performance characteristics and optimization tradeoffs in tronko-assign.

## Table of Contents

- [2026-01-01: Memory Optimization (Float vs Double Precision)](#2026-01-01-memory-optimization-float-vs-double-precision)

---

## 2026-01-01: Memory Optimization (Float vs Double Precision)

### Objective

Evaluate the tradeoff between memory usage and taxonomic accuracy when using 32-bit floats (`OPTIMIZE_MEMORY`) versus 64-bit doubles for posterior storage.

### Test Configuration

| Parameter | Value |
|-----------|-------|
| Reference Database | 16S_Bacteria (17,868 trees, 6.2 GB binary) |
| Query Data | 1,322,726 paired-end reads |
| Hardware | 4 cores |
| Database Format | Binary (.trkb) |

### Build Commands

```bash
# Double precision (default)
make clean && make

# Float precision (memory optimized)
make clean && make OPTIMIZE_MEMORY=1
```

### Results

#### Performance Comparison

| Metric | Double Precision | Float (OPTIMIZE_MEMORY) | Difference |
|--------|------------------|-------------------------|------------|
| Wall Time | 81m 38s | 82m 15s | +0.8% |
| CPU Time | 317m 26s | 322m 35s | +1.6% |
| Peak Memory (RSS) | 14,655 MB | 8,446 MB | **-42%** |
| Reference Load Time | 2.99s | 1.94s | -35% |

#### Accuracy Comparison

| Metric | Value |
|--------|-------|
| Total Reads | 1,322,726 |
| Identical Assignments | 1,257,696 (95.1%) |
| Different Assignments | 65,030 (4.9%) |

#### Nature of Differences

When assignments differed between double and float precision:

| Type | Count | % of Differences |
|------|-------|------------------|
| Float **less specific** (higher LCA) | 57,630 | 88.6% |
| Float **more specific** (lower LCA) | 5,321 | 8.2% |
| Same depth, different taxon | 2,079 | 3.2% |

Average resolution loss when less specific: **2.5 taxonomic levels**

#### Example Differences

```
Double: Bacteria;Actinobacteria;Thermoleophilia;Solirubrobacterales;Conexibacteraceae;Conexibacter;Conexibacter sp. BS10
Float:  Bacteria;Actinobacteria;Thermoleophilia;Solirubrobacterales;Conexibacteraceae;Conexibacter

Double: Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae;NA;Sphingomonadaceae bacterium 417
Float:  Bacteria;Proteobacteria;Alphaproteobacteria;Sphingomonadales;Sphingomonadaceae
```

### Conclusion

**Recommendation: Use double precision (default) for production workflows.**

The 42% memory reduction from float precision comes at the cost of ~5% of reads receiving less specific taxonomic assignments (averaging 2.5 levels less resolution). For a database of this size (17,868 trees), the 14.6 GB memory requirement is manageable on modern systems.

**Use float precision only when:**
- Memory is severely constrained (< 10 GB available)
- Running exploratory/preliminary analysis
- Filtering candidates for re-analysis with double precision

### Reproduction

```bash
# Run double precision benchmark
make clean && make
time tronko-assign -r \
  -f reference_tree.trkb \
  -a 16S_Bacteria.fasta \
  -p -1 paired_F.fasta -2 paired_R.fasta \
  -o results_double.txt \
  -V 2 -R -C 4

# Run float precision benchmark
make clean && make OPTIMIZE_MEMORY=1
time tronko-assign -r \
  -f reference_tree.trkb \
  -a 16S_Bacteria.fasta \
  -p -1 paired_F.fasta -2 paired_R.fasta \
  -o results_float.txt \
  -V 2 -R -C 4

# Compare results
diff <(cut -f2 results_double.txt) <(cut -f2 results_float.txt) | wc -l
```

---

*Add new experiments above this line, following the same format.*
