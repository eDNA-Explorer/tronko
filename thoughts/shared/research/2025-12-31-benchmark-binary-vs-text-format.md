# Benchmarking Binary vs Text Reference Format

## Objective

Compare performance characteristics of gzipped binary (`.trkb.gz`) vs gzipped text (`.txt.gz`) reference formats in tronko-assign, without the `OPTIMIZE_MEMORY` flag (using double precision).

## Test Environment

- **Database**: 16S Bacteria reference
- **Build**: Default (no `OPTIMIZE_MEMORY`, uses double precision)
- **Hardware**: Record CPU model, RAM, storage type (SSD/HDD/NVMe)

## File Sizes

| Format | Path | Size |
|--------|------|------|
| Text gzipped | `example_datasets/16S_Bacteria/reference_tree.txt.gz` | 758 MB |
| Binary gzipped | `example_datasets/16S_Bacteria/reference_tree.trkb.gz` | 352 MB |
| Binary uncompressed | `example_datasets/16S_Bacteria/reference_tree.trkb` | 5.8 GB |

## Metrics to Capture

1. **Load time**: Wall clock time to load reference database
2. **Peak memory**: Maximum RSS during loading
3. **Final memory**: RSS after reference is fully loaded
4. **CPU time**: User + system time during loading
5. **Assignment throughput**: Reads processed per second (end-to-end)

## Benchmark Method

### 1. Prepare Test Data

```bash
cd /home/jimjeffers/Work/tronko

# Ensure clean build without OPTIMIZE_MEMORY
cd tronko-assign && make clean && make && cd ..

# Verify files exist
ls -lh example_datasets/16S_Bacteria/reference_tree.txt.gz
ls -lh example_datasets/16S_Bacteria/reference_tree.trkb.gz

# Use existing test reads or create synthetic ones
# For meaningful comparison, use a larger read set (e.g., 10k-100k reads)
TEST_READS="example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta"
```

### 2. Clear System Caches (requires root)

```bash
# Drop filesystem caches between runs for fair comparison
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches
```

### 3. Run Text Format Benchmark

```bash
# Clear caches
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches

# Run with TSV logging
./tronko-assign/tronko-assign -r \
  -f example_datasets/16S_Bacteria/reference_tree.txt.gz \
  -a example_datasets/16S_Bacteria/16S_Bacteria.fasta \
  -s -g "$TEST_READS" \
  -o /tmp/results_text.txt -w \
  --tsv-log /tmp/benchmark_text.tsv \
  2>&1 | tee /tmp/benchmark_text.log

# Extract key metrics
grep "REFERENCE_LOADED\|FINAL" /tmp/benchmark_text.tsv
```

### 4. Run Binary Gzipped Format Benchmark

```bash
# Clear caches
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches

# Run with TSV logging
./tronko-assign/tronko-assign -r \
  -f example_datasets/16S_Bacteria/reference_tree.trkb.gz \
  -a example_datasets/16S_Bacteria/16S_Bacteria.fasta \
  -s -g "$TEST_READS" \
  -o /tmp/results_binary.txt -w \
  --tsv-log /tmp/benchmark_binary.tsv \
  2>&1 | tee /tmp/benchmark_binary.log

# Extract key metrics
grep "REFERENCE_LOADED\|FINAL" /tmp/benchmark_binary.tsv
```

### 5. Run Binary Uncompressed Format Benchmark (Baseline)

```bash
# Clear caches
sync && echo 3 | sudo tee /proc/sys/vm/drop_caches

# Run with TSV logging
./tronko-assign/tronko-assign -r \
  -f example_datasets/16S_Bacteria/reference_tree.trkb \
  -a example_datasets/16S_Bacteria/16S_Bacteria.fasta \
  -s -g "$TEST_READS" \
  -o /tmp/results_binary_raw.txt -w \
  --tsv-log /tmp/benchmark_binary_raw.tsv \
  2>&1 | tee /tmp/benchmark_binary_raw.log

# Extract key metrics
grep "REFERENCE_LOADED\|FINAL" /tmp/benchmark_binary_raw.tsv
```

### 6. Repeat for Statistical Significance

Run each benchmark 3-5 times and compute mean/stddev:

```bash
for i in 1 2 3; do
  echo "=== Run $i ==="
  sync && echo 3 | sudo tee /proc/sys/vm/drop_caches

  ./tronko-assign/tronko-assign -r \
    -f example_datasets/16S_Bacteria/reference_tree.trkb.gz \
    -a example_datasets/16S_Bacteria/16S_Bacteria.fasta \
    -s -g "$TEST_READS" \
    -o /tmp/results_$i.txt -w \
    --tsv-log /tmp/benchmark_run_$i.tsv 2>&1

  grep "REFERENCE_LOADED" /tmp/benchmark_run_$i.tsv
done
```

## Analysis Script

```bash
#!/bin/bash
# analyze_benchmarks.sh

echo "Format,LoadTime(s),PeakRSS(MB),FinalRSS(MB),CPUUser(s),CPUSys(s)"

for format in text binary binary_raw; do
  file="/tmp/benchmark_${format}.tsv"
  if [[ -f "$file" ]]; then
    # Extract REFERENCE_LOADED line
    line=$(grep "REFERENCE_LOADED" "$file")
    load_time=$(echo "$line" | cut -f1)
    rss=$(echo "$line" | cut -f3)
    peak=$(echo "$line" | cut -f5)
    cpu_user=$(echo "$line" | cut -f6)
    cpu_sys=$(echo "$line" | cut -f7)
    echo "${format},${load_time},${peak},${rss},${cpu_user},${cpu_sys}"
  fi
done
```

## Expected Results Template

| Metric | Text (.txt.gz) | Binary (.trkb.gz) | Binary (.trkb) |
|--------|----------------|-------------------|----------------|
| File size | 758 MB | 352 MB | 5.8 GB |
| Load time | ? sec | ? sec | ? sec |
| Peak RSS | ? GB | ? GB | ? GB |
| CPU user | ? sec | ? sec | ? sec |
| CPU sys | ? sec | ? sec | ? sec |

## Verification

Ensure all formats produce identical assignment results:

```bash
# Compare taxonomic assignments (columns 1-2)
cut -f1,2 /tmp/results_text.txt | sort > /tmp/taxa_text.txt
cut -f1,2 /tmp/results_binary.txt | sort > /tmp/taxa_binary.txt
cut -f1,2 /tmp/results_binary_raw.txt | sort > /tmp/taxa_raw.txt

diff /tmp/taxa_text.txt /tmp/taxa_binary.txt && echo "Text vs Binary: MATCH"
diff /tmp/taxa_binary.txt /tmp/taxa_raw.txt && echo "Binary gz vs raw: MATCH"
```

## Notes

- **I/O vs CPU tradeoff**: Gzipped formats trade CPU (decompression) for reduced I/O
- **Cold vs warm cache**: First run after cache clear shows true I/O performance
- **Memory baseline**: Both gzipped formats should use identical memory (same decompressed data)
- **Text parsing overhead**: Text format requires string parsing; binary is direct memory copy

## Hypothesis

Based on format analysis:
1. Binary gzipped should load **faster** than text gzipped (no string parsing)
2. Binary uncompressed should be **fastest** on NVMe/SSD (no decompression overhead)
3. Binary gzipped may be **faster** than uncompressed on slow storage (I/O bound)
4. All formats should use **identical memory** after loading (~12.3 GB for 16S Bacteria)
