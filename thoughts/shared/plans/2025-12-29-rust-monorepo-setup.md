# Rust Monorepo Setup Implementation Plan

## Overview

Set up a Rust monorepo in the `rust/` subdirectory as the first step toward porting tronko to Rust. This plan establishes the workspace foundation, implements FASTA/FASTQ parsing as the first ported component, and creates a CLI interface that demonstrates the Rust implementation works correctly.

## Current State Analysis

The tronko project consists of two C modules:
- **tronko-build**: Database builder (11 source files, ~4,500 lines)
- **tronko-assign**: Sequence assigner (63 source files including BWA/WFA2, ~15,000+ lines)

Both use gcc with `-O3` optimization and link against `-lm -pthread -lz` (and `-lrt` for tronko-assign).

### Key Discoveries:
- No existing Rust code in the repository
- CI runs on Ubuntu 24.04 (`/.github/workflows/build.yml`)
- FASTA parsing in `tronko-build/readfasta.c` and `tronko-assign/readreference.c` uses zlib for gzip support
- Example datasets exist for testing: `example_datasets/` and `tronko-build/example_datasets/`

## Desired End State

After completing this plan:
1. A `rust/` directory with a Cargo workspace containing multiple crates
2. A working FASTA/FASTQ parser using `needletail` that can read gzipped files
3. A CLI tool `tronko-rs` that can parse and count sequences from FASTA/FASTQ files
4. CI/CD that builds and tests the Rust code alongside the C code
5. Validation that Rust output matches C implementation behavior

### Verification:
```bash
# Build and test Rust workspace
cd rust && cargo build --release && cargo test

# Verify FASTA parsing matches C behavior
./rust/target/release/tronko-rs parse --input tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta
# Should output sequence count, total length, etc.
```

## What We're NOT Doing

- **Not porting tree structures yet** - Arena-based trees are the hardest part; deferred to future phases
- **Not porting alignment algorithms** - BWA/WFA2 will remain in C via FFI initially
- **Not porting likelihood calculations** - Mathematical complexity deferred
- **Not replacing tronko-build or tronko-assign** - This is additive, not a replacement
- **Not creating FFI bindings** - Pure Rust implementation for the components we port

## Implementation Approach

Following the research recommendation of "incremental pure Rust" approach for low-risk components:
1. Start with I/O (FASTA parsing) - lowest risk, immediate validation
2. Build a CLI that uses the Rust parser - proves the library works
3. Set up CI to catch regressions early

---

## Phase 1: Rust Workspace Foundation

### Overview
Create the Cargo workspace structure with initial crates for a shared library and a CLI binary.

### Changes Required:

#### 1. Create Workspace Root
**File**: `rust/Cargo.toml`

```toml
[workspace]
resolver = "2"
members = [
    "crates/tronko-core",
    "crates/tronko-rs",
]

[workspace.package]
version = "0.1.0"
edition = "2024"
authors = ["CALeDNA <caledna@ucsc.edu>"]
license = "MIT"
repository = "https://github.com/CALeDNA/tronko"

[workspace.dependencies]
# Bioinformatics
needletail = "0.6"

# Error handling
thiserror = "2"
anyhow = "1"

# CLI
clap = { version = "4", features = ["derive"] }

# Serialization (for future phases)
serde = { version = "1", features = ["derive"] }
serde_json = "1"

# Testing
pretty_assertions = "1"
```

#### 2. Create Core Library Crate
**File**: `rust/crates/tronko-core/Cargo.toml`

```toml
[package]
name = "tronko-core"
version.workspace = true
edition.workspace = true
authors.workspace = true
license.workspace = true
description = "Core library for tronko phylogenetic assignment"

[dependencies]
needletail.workspace = true
thiserror.workspace = true

[dev-dependencies]
pretty_assertions.workspace = true
```

**File**: `rust/crates/tronko-core/src/lib.rs`

```rust
//! tronko-core: Core library for phylogenetic sequence analysis
//!
//! This crate provides the foundational data structures and algorithms
//! for the tronko phylogenetic assignment system.

pub mod fasta;

pub use fasta::{FastaReader, Sequence};
```

#### 3. Create CLI Binary Crate
**File**: `rust/crates/tronko-rs/Cargo.toml`

```toml
[package]
name = "tronko-rs"
version.workspace = true
edition.workspace = true
authors.workspace = true
license.workspace = true
description = "CLI for tronko phylogenetic assignment (Rust implementation)"

[[bin]]
name = "tronko-rs"
path = "src/main.rs"

[dependencies]
tronko-core = { path = "../tronko-core" }
clap.workspace = true
anyhow.workspace = true
```

**File**: `rust/crates/tronko-rs/src/main.rs`

```rust
use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "tronko-rs")]
#[command(about = "Rust implementation of tronko phylogenetic assignment")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Parse and analyze FASTA/FASTQ files
    Parse {
        /// Input file (FASTA or FASTQ, optionally gzipped)
        #[arg(short, long)]
        input: std::path::PathBuf,

        /// Output format
        #[arg(short, long, default_value = "summary")]
        format: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Parse { input, format } => {
            cmd_parse(&input, &format)
        }
    }
}

fn cmd_parse(input: &std::path::Path, format: &str) -> Result<()> {
    use tronko_core::FastaReader;

    let reader = FastaReader::from_path(input)?;
    let mut count = 0;
    let mut total_len = 0;

    for result in reader {
        let seq = result?;
        count += 1;
        total_len += seq.sequence.len();
    }

    match format {
        "summary" => {
            println!("Sequences: {}", count);
            println!("Total bases: {}", total_len);
            println!("Average length: {:.1}", total_len as f64 / count as f64);
        }
        "json" => {
            println!(r#"{{"sequences": {}, "total_bases": {}, "avg_length": {:.1}}}"#,
                     count, total_len, total_len as f64 / count as f64);
        }
        _ => {
            anyhow::bail!("Unknown format: {}", format);
        }
    }

    Ok(())
}
```

#### 4. Create Directory Structure

```bash
mkdir -p rust/crates/tronko-core/src
mkdir -p rust/crates/tronko-rs/src
```

#### 5. Add .gitignore for Rust
**File**: `rust/.gitignore`

```
/target/
Cargo.lock
*.swp
*.swo
.DS_Store
```

### Success Criteria:

#### Automated Verification:
- [x] Workspace compiles: `cd rust && cargo build`
- [x] Workspace structure valid: `cd rust && cargo metadata --format-version=1 | head -1`

#### Manual Verification:
- [x] Directory structure looks correct in file explorer
- [x] `cargo check` shows no warnings

---

## Phase 2: FASTA/FASTQ Parsing with needletail

### Overview
Implement FASTA/FASTQ parsing using the `needletail` crate, which provides high-performance parsing comparable to C implementations with native gzip support.

### Changes Required:

#### 1. FASTA Reader Module
**File**: `rust/crates/tronko-core/src/fasta.rs`

```rust
//! FASTA/FASTQ file parsing using needletail
//!
//! This module provides a unified interface for reading both FASTA and FASTQ
//! files, with transparent support for gzip compression.

use std::path::Path;
use needletail::{parse_fastx_file, FastxReader};
use thiserror::Error;

/// Errors that can occur during FASTA/FASTQ parsing
#[derive(Error, Debug)]
pub enum FastaError {
    #[error("Failed to open file: {0}")]
    FileOpen(#[from] std::io::Error),

    #[error("Failed to parse sequence: {0}")]
    Parse(String),

    #[error("Invalid sequence data")]
    InvalidSequence,
}

/// A parsed sequence record
#[derive(Debug, Clone)]
pub struct Sequence {
    /// Sequence identifier (header without '>' or '@')
    pub id: String,
    /// The nucleotide sequence
    pub sequence: Vec<u8>,
    /// Quality scores (only for FASTQ)
    pub quality: Option<Vec<u8>>,
}

impl Sequence {
    /// Returns the sequence as a string (ASCII)
    pub fn sequence_str(&self) -> &str {
        std::str::from_utf8(&self.sequence).unwrap_or("")
    }

    /// Returns the sequence length
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns true if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

/// Reader for FASTA/FASTQ files
///
/// Automatically handles gzip-compressed files based on file extension.
///
/// # Example
///
/// ```no_run
/// use tronko_core::FastaReader;
///
/// let reader = FastaReader::from_path("sequences.fasta.gz").unwrap();
/// for result in reader {
///     let seq = result.unwrap();
///     println!("{}: {} bp", seq.id, seq.len());
/// }
/// ```
pub struct FastaReader {
    inner: Box<dyn FastxReader>,
}

impl FastaReader {
    /// Create a new reader from a file path
    ///
    /// Supports both FASTA and FASTQ formats, optionally gzip-compressed.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastaError> {
        let reader = parse_fastx_file(path.as_ref())
            .map_err(|e| FastaError::Parse(e.to_string()))?;

        Ok(Self { inner: reader })
    }
}

impl Iterator for FastaReader {
    type Item = Result<Sequence, FastaError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next() {
            Some(Ok(record)) => {
                let id = String::from_utf8_lossy(record.id()).to_string();
                let sequence = record.seq().to_vec();
                let quality = record.qual().map(|q| q.to_vec());

                Some(Ok(Sequence { id, sequence, quality }))
            }
            Some(Err(e)) => Some(Err(FastaError::Parse(e.to_string()))),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_len() {
        let seq = Sequence {
            id: "test".to_string(),
            sequence: b"ACGT".to_vec(),
            quality: None,
        };
        assert_eq!(seq.len(), 4);
        assert!(!seq.is_empty());
    }

    #[test]
    fn test_sequence_str() {
        let seq = Sequence {
            id: "test".to_string(),
            sequence: b"ACGTNN".to_vec(),
            quality: None,
        };
        assert_eq!(seq.sequence_str(), "ACGTNN");
    }
}
```

#### 2. Integration Tests
**File**: `rust/crates/tronko-core/tests/fasta_integration.rs`

```rust
//! Integration tests for FASTA parsing
//!
//! These tests use the example datasets from the main tronko repository.

use tronko_core::FastaReader;
use std::path::PathBuf;

fn example_data_path(name: &str) -> PathBuf {
    // Navigate from rust/crates/tronko-core to repository root
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop(); // crates/
    path.pop(); // rust/
    path.push(name);
    path
}

#[test]
fn test_parse_charadriiformes_msa() {
    let path = example_data_path("tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta");

    if !path.exists() {
        eprintln!("Skipping test: example file not found at {:?}", path);
        return;
    }

    let reader = FastaReader::from_path(&path).expect("Failed to open FASTA file");

    let sequences: Vec<_> = reader.collect::<Result<Vec<_>, _>>().expect("Failed to parse");

    assert!(!sequences.is_empty(), "Expected at least one sequence");

    // All sequences should have non-empty IDs
    for seq in &sequences {
        assert!(!seq.id.is_empty(), "Sequence ID should not be empty");
        assert!(!seq.is_empty(), "Sequence should not be empty");
    }
}

#[test]
fn test_parse_single_end_reads() {
    let path = example_data_path("example_datasets/single_tree/missingreads_singleend_150bp_2error.fasta");

    if !path.exists() {
        eprintln!("Skipping test: example file not found at {:?}", path);
        return;
    }

    let reader = FastaReader::from_path(&path).expect("Failed to open FASTA file");

    let sequences: Vec<_> = reader.collect::<Result<Vec<_>, _>>().expect("Failed to parse");

    assert!(!sequences.is_empty(), "Expected at least one sequence");

    // Single-end reads should be ~150bp (with some variation)
    for seq in &sequences {
        assert!(seq.len() >= 100 && seq.len() <= 200,
                "Expected ~150bp reads, got {} bp", seq.len());
    }
}
```

### Success Criteria:

#### Automated Verification:
- [x] Library compiles: `cd rust && cargo build -p tronko-core`
- [x] Unit tests pass: `cd rust && cargo test -p tronko-core`
- [x] Integration tests pass with example data: `cd rust && cargo test -p tronko-core --test fasta_integration`
- [x] No clippy warnings: `cd rust && cargo clippy -p tronko-core -- -D warnings`

#### Manual Verification:
- [x] Parse Charadriiformes_MSA.fasta and verify sequence count matches `grep -c "^>" file`
- [ ] Parse a gzipped FASTA to verify gzip support works

---

## Phase 3: CLI Interface

### Overview
Complete the CLI implementation and verify it produces correct output for the example datasets.

### Changes Required:

#### 1. Enhanced CLI with Verbose Output
**File**: `rust/crates/tronko-rs/src/main.rs` (update)

Replace the existing `cmd_parse` function with enhanced functionality:

```rust
use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "tronko-rs")]
#[command(version)]
#[command(about = "Rust implementation of tronko phylogenetic assignment")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Parse and analyze FASTA/FASTQ files
    Parse {
        /// Input file (FASTA or FASTQ, optionally gzipped)
        #[arg(short, long)]
        input: PathBuf,

        /// Output format: summary, json, or ids
        #[arg(short, long, default_value = "summary")]
        format: String,

        /// Show first N sequences (0 = all)
        #[arg(short = 'n', long, default_value = "0")]
        limit: usize,
    },

    /// Validate a FASTA/FASTQ file
    Validate {
        /// Input file to validate
        #[arg(short, long)]
        input: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Parse { input, format, limit } => {
            cmd_parse(&input, &format, limit)
        }
        Commands::Validate { input } => {
            cmd_validate(&input)
        }
    }
}

fn cmd_parse(input: &std::path::Path, format: &str, limit: usize) -> Result<()> {
    use tronko_core::FastaReader;

    let reader = FastaReader::from_path(input)?;
    let mut count = 0;
    let mut total_len = 0;
    let mut min_len = usize::MAX;
    let mut max_len = 0;
    let mut ids = Vec::new();

    for result in reader {
        let seq = result?;
        count += 1;
        let len = seq.len();
        total_len += len;
        min_len = min_len.min(len);
        max_len = max_len.max(len);

        if format == "ids" && (limit == 0 || ids.len() < limit) {
            ids.push(seq.id.clone());
        }
    }

    if count == 0 {
        anyhow::bail!("No sequences found in file");
    }

    match format {
        "summary" => {
            println!("File: {}", input.display());
            println!("Sequences: {}", count);
            println!("Total bases: {}", total_len);
            println!("Average length: {:.1}", total_len as f64 / count as f64);
            println!("Min length: {}", min_len);
            println!("Max length: {}", max_len);
        }
        "json" => {
            println!(
                r#"{{"file": "{}", "sequences": {}, "total_bases": {}, "avg_length": {:.1}, "min_length": {}, "max_length": {}}}"#,
                input.display(), count, total_len,
                total_len as f64 / count as f64, min_len, max_len
            );
        }
        "ids" => {
            for id in &ids {
                println!("{}", id);
            }
            if limit > 0 && count > limit {
                println!("... and {} more", count - limit);
            }
        }
        _ => {
            anyhow::bail!("Unknown format: {}. Use: summary, json, or ids", format);
        }
    }

    Ok(())
}

fn cmd_validate(input: &std::path::Path) -> Result<()> {
    use tronko_core::FastaReader;

    let reader = FastaReader::from_path(input)?;
    let mut count = 0;
    let mut errors = Vec::new();

    for (i, result) in reader.enumerate() {
        match result {
            Ok(seq) => {
                count += 1;
                // Check for empty sequences
                if seq.is_empty() {
                    errors.push(format!("Sequence {} ({}) is empty", i + 1, seq.id));
                }
                // Check for non-ACGTN characters
                for (j, &b) in seq.sequence.iter().enumerate() {
                    if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n' | b'-') {
                        errors.push(format!(
                            "Sequence {} ({}): invalid character '{}' at position {}",
                            i + 1, seq.id, b as char, j + 1
                        ));
                        break;
                    }
                }
            }
            Err(e) => {
                errors.push(format!("Failed to parse sequence {}: {}", i + 1, e));
            }
        }
    }

    if errors.is_empty() {
        println!("Valid: {} sequences parsed successfully", count);
        Ok(())
    } else {
        println!("Invalid: found {} errors", errors.len());
        for error in &errors {
            println!("  - {}", error);
        }
        anyhow::bail!("Validation failed with {} errors", errors.len());
    }
}
```

#### 2. Update CI to Build Rust
**File**: `.github/workflows/build.yml` (add after line 48)

Add the following steps to the existing workflow:

```yaml
      # Rust build steps
      - name: Install Rust toolchain
        uses: dtolnay/rust-toolchain@stable
        with:
          components: clippy

      - name: Cache Rust dependencies
        uses: Swatinem/rust-cache@v2
        with:
          workspaces: rust

      - name: Build Rust workspace
        run: cargo build --release
        working-directory: ./rust

      - name: Run Rust tests
        run: cargo test
        working-directory: ./rust

      - name: Run Clippy
        run: cargo clippy -- -D warnings
        working-directory: ./rust

      - name: Copy Rust binaries
        run: cp rust/target/release/tronko-rs bin/

      - name: Test tronko-rs parse
        run: |
          ./bin/tronko-rs parse --input tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta
```

### Success Criteria:

#### Automated Verification:
- [x] CLI builds: `cd rust && cargo build -p tronko-rs --release`
- [x] CLI help works: `./rust/target/release/tronko-rs --help`
- [x] Parse command works: `./rust/target/release/tronko-rs parse --input tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`
- [x] Validate command works: `./rust/target/release/tronko-rs validate --input tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`
- [x] JSON output is valid: `./rust/target/release/tronko-rs parse --input ... --format json | python3 -m json.tool`
- [ ] CI workflow passes: GitHub Actions shows green checkmark

#### Manual Verification:
- [x] Sequence count from `tronko-rs parse` matches `grep -c "^>" file.fasta`
- [x] Binary size is reasonable (< 10MB for release build)
- [ ] Parse performance is comparable to C (subjective, not a blocker)

---

## Testing Strategy

### Unit Tests:
- `Sequence` struct methods (len, is_empty, sequence_str)
- Error handling for invalid files

### Integration Tests:
- Parse example FASTA files from repository
- Verify sequence counts match expected values
- Test gzip support (if gzipped example files exist)

### Manual Testing Steps:
1. Build release binary: `cd rust && cargo build --release`
2. Parse MSA file: `./target/release/tronko-rs parse -i ../../tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`
3. Compare count with grep: `grep -c "^>" ../../tronko-build/example_datasets/single_tree/Charadriiformes_MSA.fasta`
4. Test JSON output: `./target/release/tronko-rs parse -i ... -f json`
5. Test validation: `./target/release/tronko-rs validate -i ...`

---

## Performance Considerations

- `needletail` uses SIMD acceleration where available - performance should be comparable to or better than C
- Using `Vec<u8>` for sequences avoids UTF-8 validation overhead
- Release builds with `-O3` equivalent optimizations via `--release` flag

---

## Migration Notes

This is additive - no existing functionality is changed. The Rust code exists alongside the C code.

Future phases will:
1. Add FFI bindings to call C code from Rust
2. Incrementally port more components
3. Eventually provide Rust alternatives to the C executables

---

## References

- Feasibility research: `thoughts/shared/research/2025-12-29-rust-port-feasibility.md`
- Data flow analysis: `thoughts/shared/research/2025-12-29-tronko-assign-data-flow.md`
- needletail crate: https://crates.io/crates/needletail
- clap crate: https://crates.io/crates/clap
