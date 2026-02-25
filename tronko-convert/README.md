# tronko-convert

A standalone utility for converting tronko reference databases between text and binary formats.

## Overview

`tronko-convert` converts tronko-build's text-based `reference_tree.txt` format to an optimized binary format (`.trkb`) for faster loading in tronko-assign. The binary format provides:

- **~2.7x smaller file size** (40MB text → 15MB binary for typical databases)
- **50-70x faster loading** (projected, via bulk reads instead of line-by-line parsing)
- **Data integrity validation** via CRC-32 checksums

## Building

```bash
cd tronko-convert
make          # Release build with -O3 optimization
make debug    # Debug build with symbols
make clean    # Remove build artifacts
```

**Dependencies**: gcc, zlib (`-lz`), math library (`-lm`)

## Usage

### Convert text to binary (default)

```bash
./tronko-convert -i reference_tree.txt -o reference_tree.trkb
```

### Convert binary back to text

```bash
./tronko-convert -i reference_tree.trkb -o reference_tree.txt -t
```

### Options

| Flag | Description |
|------|-------------|
| `-i <file>` | Input file (required) |
| `-o <file>` | Output file (required) |
| `-t` | Output as text format (default: binary) |
| `-v` | Verbose output with progress information |
| `-h` | Show help message |

### Examples

```bash
# Convert with verbose output
./tronko-convert -i reference_tree.txt -o reference_tree.trkb -v

# Round-trip conversion for validation
./tronko-convert -i reference_tree.txt -o test.trkb -v
./tronko-convert -i test.trkb -o roundtrip.txt -t -v
diff reference_tree.txt roundtrip.txt  # Will show float precision differences
```

## Testing

```bash
make test  # Runs full conversion test with example dataset
```

Or use the standalone test script:

```bash
./test_roundtrip.sh
```

## Input Format Detection

The tool automatically detects input format based on file contents:

| Magic Bytes | Format |
|-------------|--------|
| `0x89 'T' 'R' 'K'` | Binary (.trkb) |
| `0x1f 0x8b` | Gzipped text |
| ASCII digits | Plain text |

Both gzipped and uncompressed text files are supported as input.

## Binary Format Specification (v1.0)

### File Structure

```
+---------------------------+
| File Header (64 bytes)    |
+---------------------------+
| Global Metadata (16 bytes)|
+---------------------------+
| Tree Metadata Section     |
+---------------------------+
| Taxonomy String Table     |
+---------------------------+
| Node Structure Section    |
+---------------------------+
| Posterior Data Section    |
+---------------------------+
| Footer (8 bytes)          |
+---------------------------+
```

### File Header (64 bytes)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint8[4] | Magic: `0x89 'T' 'R' 'K'` |
| 4 | 1 | uint8 | Version major (1) |
| 5 | 1 | uint8 | Version minor (0) |
| 6 | 1 | uint8 | Endianness: 0x01=little |
| 7 | 1 | uint8 | Precision: 0x01=float |
| 8 | 4 | uint32 | Header CRC-32 (of bytes 0-7) |
| 12 | 4 | uint32 | Reserved |
| 16 | 8 | uint64 | Taxonomy section offset |
| 24 | 8 | uint64 | Node section offset |
| 32 | 8 | uint64 | Posterior section offset |
| 40 | 8 | uint64 | Total file size |
| 48 | 16 | - | Reserved |

### Global Metadata (16 bytes)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | Number of trees |
| 4 | 4 | int32 | Max node name length |
| 8 | 4 | int32 | Max taxonomy name length |
| 12 | 4 | int32 | Max taxonomy line length |

### Tree Metadata (12 bytes per tree)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | numbase (MSA alignment length) |
| 4 | 4 | int32 | root (root node index) |
| 8 | 4 | int32 | numspec (number of species) |

### Taxonomy Section

Per tree:
- 4 bytes: taxonomy data size for this tree
- 4 bytes: reserved

Per species (numspec per tree), per taxonomy level (7 levels):
- 2 bytes: string length (including null terminator)
- N bytes: null-terminated string

Taxonomy levels: domain, phylum, class, order, family, genus, species

### Node Section

Per tree:
- 4 bytes: node count (should equal `2 * numspec - 1`)

Per node (32 bytes fixed):

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | int32 | up[0] (first child, -1 for leaf) |
| 4 | 4 | int32 | up[1] (second child, -1 for leaf) |
| 8 | 4 | int32 | down (parent index) |
| 12 | 4 | int32 | depth |
| 16 | 4 | int32 | taxIndex[0] (species index) |
| 20 | 4 | int32 | taxIndex[1] (taxonomy level) |
| 24 | 4 | uint32 | Name offset (0 if internal node) |
| 28 | 4 | uint32 | Reserved |

Followed by name table (length-prefixed strings for leaf nodes only).

### Posterior Section

Per tree, per node, per position (numbase):
- 4 floats (16 bytes): posterior probabilities for A, C, G, T

Total size per tree: `(2 * numspec - 1) * numbase * 4 * sizeof(float)` bytes

### Footer (8 bytes)

| Offset | Size | Type | Description |
|--------|------|------|-------------|
| 0 | 4 | uint32 | CRC-32 of all preceding data |
| 4 | 4 | uint32 | Magic footer: `0x454E4421` ("END!") |

## Data Type Considerations

### Float vs Double Precision

The binary format uses 32-bit floats for posterior probabilities instead of 64-bit doubles used in the text format. This:

- Reduces file size by ~50% for posterior data
- Maintains sufficient precision for phylogenetic placement
- Matches the `USE_FLOAT_PP` compile option already available in tronko-assign

When converting binary back to text, posterior values will have reduced precision compared to the original. This is expected and acceptable for the assignment algorithm.

### Endianness

The format uses little-endian byte order throughout, which is standard for modern x86/AMD64 and ARM systems. Big-endian systems are not supported.

## File Size Calculation

For a database with T trees, S species, and P positions:

```
Header:     64 bytes
Global:     16 bytes
Tree meta:  T * 12 bytes
Taxonomy:   ~(S * 7 * avg_name_len) bytes
Nodes:      T * ((2S-1) * 32 + leaf_names) bytes
Posteriors: T * (2S-1) * P * 16 bytes
Footer:     8 bytes
```

The posterior section dominates for typical databases.

## Error Handling

The tool uses error-checked I/O wrappers that will:

- Print descriptive error messages to stderr
- Exit with non-zero status on any I/O failure
- Validate magic numbers and format versions on read

## Source Files

| File | Description |
|------|-------------|
| `tronko-convert.c` | Main entry point, CLI parsing |
| `format_common.h/c` | Data structures, format detection, memory management |
| `format_text.h/c` | Text format reader (gzip-aware) and writer |
| `format_binary.h/c` | Binary format reader and writer |
| `utils.h/c` | Error-checked I/O, little-endian helpers |
| `crc32.h/c` | CRC-32 checksum implementation |

## Related Documentation

- Original research: `thoughts/shared/research/2025-12-30-binary-format-conversion-tool.md`
- Implementation plan: `thoughts/shared/plans/2025-12-31-binary-format-converter-phase1.md`
- Text format output: `tronko-build/printtree.c`
- Text format parsing: `tronko-assign/readreference.c`

## Future Work (Phase 2+)

- Integration with tronko-assign for auto-detection and binary loading
- Optional native binary output in tronko-build
- Memory-mapped I/O for even faster loading
