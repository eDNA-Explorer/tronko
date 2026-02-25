# tronko-bench Benchmarking Workflow - Unified Implementation Plan

## Overview

Extend the tronko-bench TUI to support a complete benchmarking workflow: creating benchmarks with configurable parameters, storing results persistently, viewing historical runs, and comparing up to two benchmarks side-by-side.

This plan unifies two earlier drafts, taking the best aspects of each:
- **Detailed data model and storage** from the original plan
- **Granular wizard steps for project/marker discovery** from the TUI-focused plan
- **Explicit binary selection** for OPTIMIZE_MEMORY builds
- **Both quick hash comparison and detailed score analysis**

## Current State

The TUI currently has three tabs (Databases, Queries, MSA Files) with:
- Basic list navigation and fuzzy search
- Database conversion to binary format
- Config file generation

**Not yet supported:**
- Running actual benchmarks
- Storing/retrieving benchmark results
- Comparing benchmark runs

## Goals

1. **Build Benchmark Workflow**: Multi-step wizard for selecting database, project, marker, query mode, binary, and execution parameters
2. **Result Storage**: Persist benchmark metadata, performance metrics, and outputs to disk
3. **History View**: Display past benchmarks with full details and performance summaries
4. **Comparison Mode**: Select up to two benchmarks to compare speed, memory, and result differences

## Non-Goals

- Real-time progress streaming (v1 uses blocking execution with status updates)
- Automated benchmark scheduling/batching
- Cloud storage of results
- Graphical charts (text-based summaries only)
- Automatic rebuilding with OPTIMIZE_MEMORY (user must pre-build both binaries)

---

## Data Model

### Benchmark Configuration

```rust
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use chrono::{DateTime, Utc};

/// A complete benchmark configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    /// Unique identifier (UUID v4)
    pub id: String,

    /// Human-readable name (auto-generated from selections)
    pub name: String,

    /// When the benchmark was created
    pub created_at: DateTime<Utc>,

    /// Reference database configuration
    pub database: DatabaseConfig,

    /// Query files configuration
    pub queries: QueryConfig,

    /// Execution parameters
    pub params: ExecutionParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseConfig {
    /// Path to reference database file
    pub reference_path: PathBuf,

    /// Path to associated FASTA file (for BWA index)
    pub fasta_path: PathBuf,

    /// Database format: "text", "text (gzipped)", "binary", "binary (gzipped)"
    pub format: String,

    /// Number of trees in database
    pub num_trees: u32,

    /// Size on disk in bytes
    pub size_bytes: u64,

    /// Parent directory name (for display context, e.g., "16S_Bacteria")
    pub parent_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryConfig {
    /// Project name (directory name from example_project)
    pub project: String,

    /// Marker name (e.g., "16S", "ITS", "COI")
    pub marker: String,

    /// Read mode
    pub read_mode: ReadMode,

    /// Paths to query files
    pub files: QueryFiles,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum ReadMode {
    Paired,
    UnpairedForward,
    UnpairedReverse,
}

impl ReadMode {
    pub fn label(&self) -> &'static str {
        match self {
            ReadMode::Paired => "Paired-end",
            ReadMode::UnpairedForward => "Unpaired (Forward)",
            ReadMode::UnpairedReverse => "Unpaired (Reverse)",
        }
    }

    pub fn tronko_flags(&self) -> &'static str {
        match self {
            ReadMode::Paired => "-p",
            ReadMode::UnpairedForward | ReadMode::UnpairedReverse => "-s",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum QueryFiles {
    Single { path: PathBuf },
    Paired { forward: PathBuf, reverse: PathBuf },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExecutionParams {
    /// Path to tronko-assign binary
    pub tronko_assign_path: PathBuf,

    /// Whether this is an OPTIMIZE_MEMORY build
    pub optimize_memory: bool,

    /// Number of CPU cores to use (default: 1)
    pub cores: u32,

    /// Use Needleman-Wunsch instead of WFA
    pub use_needleman_wunsch: bool,

    /// Batch size (lines per batch, default: 50000)
    pub batch_size: u32,
}

impl Default for ExecutionParams {
    fn default() -> Self {
        Self {
            tronko_assign_path: PathBuf::new(),
            optimize_memory: false,
            cores: 1,
            use_needleman_wunsch: false,
            batch_size: 50000,
        }
    }
}
```

### Benchmark Results

```rust
/// Results from a completed benchmark run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    /// Reference to the config ID
    pub config_id: String,

    /// When the benchmark started
    pub started_at: DateTime<Utc>,

    /// When the benchmark completed
    pub completed_at: DateTime<Utc>,

    /// Exit status
    pub status: BenchmarkStatus,

    /// Performance metrics (if successful)
    pub metrics: Option<PerformanceMetrics>,

    /// Assignment statistics (if successful)
    pub assignment_stats: Option<AssignmentStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BenchmarkStatus {
    Success,
    Failed { error: String },
    Cancelled,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    /// Total wall clock time in seconds
    pub wall_time_secs: f64,

    /// Peak RSS memory in bytes
    pub peak_rss_bytes: u64,

    /// Final RSS memory in bytes
    pub final_rss_bytes: u64,

    /// Reads processed per second
    pub reads_per_second: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssignmentStats {
    /// Total reads processed
    pub total_reads: u64,

    /// Reads successfully assigned
    pub assigned_reads: u64,

    /// Unique taxonomic paths found
    pub unique_taxa: u64,

    /// MD5 hash of sorted assignments (for quick comparison)
    pub assignments_hash: String,
}
```

### Storage Structure

Benchmarks stored in `~/.tronko-bench/` (cross-platform via `directories` crate):

```
~/.tronko-bench/
├── config.toml              # Global settings (future use)
├── index.json               # Quick-access index of all benchmarks
└── benchmarks/
    └── {uuid}/
        ├── config.json      # BenchmarkConfig
        ├── result.json      # BenchmarkResult (after completion)
        ├── output.tsv       # tronko-assign output
        ├── metrics.tsv      # --tsv-log output
        └── stderr.log       # Captured stderr
```

---

## UI Mockups

This section provides comprehensive ASCII mockups of all new and modified TUI interfaces.

### Overall Layout

The TUI maintains a consistent layout across all modes:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ Tab Bar (with item counts)                                                   │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  Left Panel (60%)              │  Right Panel (40%)                          │
│  - List of items               │  - Details for selected item                │
│                                │                                             │
│                                │                                             │
├──────────────────────────────────────────────────────────────────────────────┤
│ Status Bar (context-sensitive hints)                                         │
└──────────────────────────────────────────────────────────────────────────────┘
```

Overlays (wizard, comparison, help) render on top of this base layout.

---

### Enhanced Databases Tab

Shows parent directory context for each database:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ Databases (4)          Queries (12)         MSA Files (3)       Benchmarks   │
├────────────────────────────────────────────────┬─────────────────────────────┤
│  Reference Databases                           │  Details                    │
│ ───────────────────────────────────────────────│─────────────────────────────│
│                                                │                             │
│ > 16S_Bacteria/reference_tree.trkb       156   │  Name: reference_tree.trkb  │
│     5.78 GiB  [BIN]                            │                             │
│   16S_Bacteria/reference_tree.trkb.gz    156   │  Format: binary (.trkb)     │
│     352 MiB   [BIN.gz]                         │                             │
│   16S_Bacteria/reference_tree.txt.gz     156   │  Compressed: No             │
│     757 MiB   [TXT.gz]                         │                             │
│   single_tree/reference_tree.txt         1     │  Trees: 156                 │
│     38.1 MiB  [TXT]                            │                             │
│                                                │  Size: 5.78 GiB             │
│                                                │                             │
│                                                │  Path:                      │
│                                                │  /home/user/tronko/         │
│                                                │  example_datasets/          │
│                                                │  16S_Bacteria/              │
│                                                │  reference_tree.trkb        │
│                                                │                             │
│                                                │  Press 'c' to convert       │
│                                                │  to binary format           │
│                                                │                             │
├────────────────────────────────────────────────┴─────────────────────────────┤
│ /home/user/tronko | Tab:switch j/k:nav /:search b:benchmark c:convert ?:help │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmarks Tab - Empty State

When no benchmarks have been run yet:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ Databases (4)          Queries (12)         MSA Files (3)       Benchmarks   │
├────────────────────────────────────────────────┬─────────────────────────────┤
│  Benchmarks                                    │  Details                    │
│ ───────────────────────────────────────────────│─────────────────────────────│
│                                                │                             │
│                                                │                             │
│                                                │                             │
│          No benchmarks yet.                    │     No benchmark selected   │
│                                                │                             │
│          Press 'b' to create your first        │                             │
│          benchmark.                            │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
│                                                │                             │
├────────────────────────────────────────────────┴─────────────────────────────┤
│ b:new benchmark  ?:help  q:quit                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmarks Tab - With Items

Showing benchmark list with selection checkboxes:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│ Databases (4)          Queries (12)         MSA Files (3)       Benchmarks(5)│
├────────────────────────────────────────────────┬─────────────────────────────┤
│  Benchmarks                                    │  Details                    │
│ ───────────────────────────────────────────────│─────────────────────────────│
│                                                │                             │
│ > [ ] 16S paired (binary)         45.2s       │  Name:                      │
│       64.8 MB   2 min ago         [OK]        │  16S paired (binary)        │
│   [*] 16S paired (text)           89.7s       │                             │
│       128.4 MB  5 min ago         [OK]        │  Status: Success            │
│   [*] 16S paired (binary+optmem)  38.1s       │                             │
│       42.3 MB   1 hour ago        [OK]        │  Duration: 45.2 seconds     │
│   [ ] ITS single (binary)         12.4s       │  Peak Memory: 64.8 MB       │
│       28.1 MB   2 hours ago       [OK]        │                             │
│   [ ] 16S paired (text)           N/A         │  ─────────────────────────  │
│       N/A       3 hours ago       [FAIL]      │  Database:                  │
│                                                │    16S_Bacteria/            │
│                                                │    reference_tree.trkb      │
│                                                │    Format: binary           │
│                                                │    Trees: 156               │
│                                                │                             │
│                                                │  Queries:                   │
│                                                │    Project: 16S             │
│                                                │    Marker: 16S              │
│                                                │    Mode: paired-end         │
│                                                │                             │
│                                                │  Binary: standard           │
│                                                │                             │
│                                                │  Results:                   │
│                                                │    Assigned: 9,847 (98.5%)  │
│                                                │    Unique Taxa: 142         │
│                                                │                             │
│                                                │  Timestamp:                 │
│                                                │    2026-01-01 09:15:23      │
│                                                │                             │
├────────────────────────────────────────────────┴─────────────────────────────┤
│ Space:select(2) b:new d:delete | 2 selected: press 'c' to compare            │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 1: Select Database

Full-screen overlay with fuzzy search:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 1 of 6: Select Reference Database                                 │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  [1.DB] → 2.Project → 3.Marker → 4.Query → 5.Binary → 6.Confirm         │ │
│  │                                                                         │ │
│  │  Filter: 16S_______________                                             │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  > 16S_Bacteria/reference_tree.trkb                                     │ │
│  │      156 trees    5.78 GiB    [BIN]                                     │ │
│  │                                                                         │ │
│  │    16S_Bacteria/reference_tree.trkb.gz                                  │ │
│  │      156 trees    352 MiB     [BIN.gz]                                  │ │
│  │                                                                         │ │
│  │    16S_Bacteria/reference_tree.txt.gz                                   │ │
│  │      156 trees    757 MiB     [TXT.gz]                                  │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓ navigate   Enter select   / search   Esc cancel                     │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 2: Select Project

Shows discovered projects with their available markers:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 2 of 6: Select Project                                            │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → [2.Project] → 3.Marker → 4.Query → 5.Binary → 6.Confirm       │ │
│  │                                                                         │ │
│  │  Selected Database: 16S_Bacteria/reference_tree.trkb                    │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Available Projects (from example_project/):                            │ │
│  │                                                                         │ │
│  │  > 16S                                                                  │ │
│  │      Markers: 16S                                                       │ │
│  │      Modes: paired, unpaired_F, unpaired_R                              │ │
│  │                                                                         │ │
│  │    43f82f91-3857-443d-8c15-e3693be268d6                                 │ │
│  │      Markers: 16S, ITS                                                  │ │
│  │      Modes: paired                                                      │ │
│  │                                                                         │ │
│  │    c8939113-7f18-4103-83d8-50a72f8b4297                                 │ │
│  │      Markers: 16S                                                       │ │
│  │      Modes: paired, unpaired_F                                          │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓ navigate   Enter select   ← back   Esc cancel                       │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 3: Select Marker

Shows markers available in the selected project:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 3 of 6: Select Marker                                             │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → 2.Project ✓ → [3.Marker] → 4.Query → 5.Binary → 6.Confirm     │ │
│  │                                                                         │ │
│  │  Database: 16S_Bacteria/reference_tree.trkb                             │ │
│  │  Project:  16S                                                          │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Available Markers:                                                     │ │
│  │                                                                         │ │
│  │  > 16S                                                                  │ │
│  │      Available modes:                                                   │ │
│  │        • paired     (paired_F.fasta.zst + paired_R.fasta.zst)           │ │
│  │        • unpaired_F (unpaired_F.fasta.zst)                              │ │
│  │        • unpaired_R (unpaired_R.fasta.zst)                              │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓ navigate   Enter select   ← back   Esc cancel                       │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 4: Select Query Mode

Choose between paired and unpaired reads:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 4 of 6: Select Query Type                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → 2.Project ✓ → 3.Marker ✓ → [4.Query] → 5.Binary → 6.Confirm   │ │
│  │                                                                         │ │
│  │  Database: 16S_Bacteria/reference_tree.trkb                             │ │
│  │  Project:  16S                                                          │ │
│  │  Marker:   16S                                                          │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Select read type:                                                      │ │
│  │                                                                         │ │
│  │  > Paired-end reads                                                     │ │
│  │      Uses: paired_F.fasta.zst + paired_R.fasta.zst                      │ │
│  │      tronko-assign flags: -p -1 <forward> -2 <reverse>                  │ │
│  │                                                                         │ │
│  │    Unpaired (Forward) reads                                             │ │
│  │      Uses: unpaired_F.fasta.zst                                         │ │
│  │      tronko-assign flags: -s -g <file>                                  │ │
│  │                                                                         │ │
│  │    Unpaired (Reverse) reads                                             │ │
│  │      Uses: unpaired_R.fasta.zst                                         │ │
│  │      tronko-assign flags: -s -g <file>                                  │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓ navigate   Enter select   ← back   Esc cancel                       │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 5: Select Binary

Choose between standard and OPTIMIZE_MEMORY builds:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 5 of 6: Select tronko-assign Binary                               │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → 2.Project ✓ → 3.Marker ✓ → 4.Query ✓ → [5.Binary] → 6.Confirm │ │
│  │                                                                         │ │
│  │  Database: 16S_Bacteria/reference_tree.trkb                             │ │
│  │  Project:  16S / Marker: 16S / Mode: paired-end                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Select tronko-assign binary:                                           │ │
│  │                                                                         │ │
│  │  > tronko-assign (standard)                                             │ │
│  │      Path: ./tronko-assign/tronko-assign                                │ │
│  │      Precision: double (64-bit floats)                                  │ │
│  │      Memory: Higher memory usage                                        │ │
│  │                                                                         │ │
│  │    tronko-assign-optimized (OPTIMIZE_MEMORY)                            │ │
│  │      Path: ./tronko-assign/tronko-assign-optimized                      │ │
│  │      Precision: float (32-bit floats)                                   │ │
│  │      Memory: ~50% less memory usage                                     │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Tip: Build optimized version with: cd tronko-assign && make clean &&   │ │
│  │       make OPTIMIZE_MEMORY=1 && mv tronko-assign tronko-assign-optimized│ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓ navigate   Enter select   ← back   Esc cancel                       │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 5 (No Binaries Found)

Error state when no tronko-assign binaries are discovered:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 5 of 6: Select tronko-assign Binary                               │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → 2.Project ✓ → 3.Marker ✓ → 4.Query ✓ → [5.Binary] → 6.Confirm │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │   ⚠  No tronko-assign binaries found!                                   │ │
│  │                                                                         │ │
│  │   Searched locations:                                                   │ │
│  │     • ./tronko-assign/tronko-assign                                     │ │
│  │     • ./tronko-assign/tronko-assign-optimized                           │ │
│  │     • System PATH                                                       │ │
│  │                                                                         │ │
│  │   To build tronko-assign:                                               │ │
│  │                                                                         │ │
│  │     cd tronko-assign && make                                            │ │
│  │                                                                         │ │
│  │   To build with memory optimization:                                    │ │
│  │                                                                         │ │
│  │     cd tronko-assign && make clean && make OPTIMIZE_MEMORY=1            │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ← back   Esc cancel                                                    │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Wizard - Step 6: Confirm & Run

Final confirmation with full command preview:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 6 of 6: Confirm & Run                                             │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  1.DB ✓ → 2.Project ✓ → 3.Marker ✓ → 4.Query ✓ → 5.Binary ✓ → [6.Run]   │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  BENCHMARK CONFIGURATION                                                │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Database:                                                              │ │
│  │    Path:   16S_Bacteria/reference_tree.trkb                             │ │
│  │    Format: binary (.trkb)                                               │ │
│  │    Trees:  156                                                          │ │
│  │    Size:   5.78 GiB                                                     │ │
│  │                                                                         │ │
│  │  Query Files:                                                           │ │
│  │    Project: 16S                                                         │ │
│  │    Marker:  16S                                                         │ │
│  │    Mode:    paired-end                                                  │ │
│  │    Files:   paired_F.fasta.zst, paired_R.fasta.zst                      │ │
│  │                                                                         │ │
│  │  Binary:                                                                │ │
│  │    Path:       ./tronko-assign/tronko-assign                            │ │
│  │    Optimized:  No (double precision)                                    │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Command to execute:                                                    │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  tronko-assign -r \                                                     │ │
│  │    -f example_datasets/16S_Bacteria/reference_tree.trkb \               │ │
│  │    -a example_datasets/16S_Bacteria/16S_Bacteria.fasta \                │ │
│  │    -p -1 example_project/16S/16S/paired/paired_F.fasta.zst \            │ │
│  │       -2 example_project/16S/16S/paired/paired_R.fasta.zst \            │ │
│  │    -o ~/.tronko-bench/benchmarks/<uuid>/output.tsv \                    │ │
│  │    -R -T --tsv-log ~/.tronko-bench/benchmarks/<uuid>/metrics.tsv        │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │              [ Run Benchmark ]              [ Cancel ]                  │ │
│  │                    ▲                                                    │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Enter run   Tab switch button   ← back   Esc cancel                    │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Running State

Shown while tronko-assign is executing:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Running Benchmark ─────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  16S paired (binary)                                                    │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │     ████████████████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  Running...  │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  Elapsed:  00:01:23                                                     │ │
│  │                                                                         │ │
│  │  Status:   Processing reads...                                          │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  This may take several minutes for large datasets.                      │ │
│  │  The TUI will update when the benchmark completes.                      │ │
│  │                                                                         │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Ctrl+C cancel benchmark                                                │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Complete - Success

Shown immediately after benchmark completes:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Benchmark Complete ────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  ✓ 16S paired (binary)                                                  │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Status:       SUCCESS                                                  │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Performance                                                            │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Duration:     45.2 seconds                                             │ │
│  │  Peak Memory:  64.8 MB                                                  │ │
│  │  Reads/sec:    221.2                                                    │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Results                                                                │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Total Reads:  10,000                                                   │ │
│  │  Assigned:     9,847 (98.5%)                                            │ │
│  │  Unique Taxa:  142                                                      │ │
│  │                                                                         │ │
│  │  Output saved to:                                                       │ │
│  │  ~/.tronko-bench/benchmarks/abc123.../output.tsv                        │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Enter close and view in Benchmarks tab                                 │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Benchmark Complete - Failure

Shown when benchmark fails:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Benchmark Failed ──────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  ✗ 16S paired (text)                                                    │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  Status:       FAILED                                                   │ │
│  │  Duration:     12.3 seconds                                             │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Error                                                                  │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │                                                                         │ │
│  │  tronko-assign exited with code 1                                       │ │
│  │                                                                         │ │
│  │  Stderr output:                                                         │ │
│  │  ┌─────────────────────────────────────────────────────────────────┐    │ │
│  │  │ Error: Could not open reference FASTA file                      │    │ │
│  │  │ Path: example_datasets/16S_Bacteria/16S_Bacteria.fasta          │    │ │
│  │  │ Reason: No such file or directory                               │    │ │
│  │  └─────────────────────────────────────────────────────────────────┘    │ │
│  │                                                                         │ │
│  │  Full stderr saved to:                                                  │ │
│  │  ~/.tronko-bench/benchmarks/abc123.../stderr.log                        │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Enter close                                                            │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Comparison View - Identical Results

When two benchmarks have identical assignments:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Benchmark Comparison ──────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Comparing 2 benchmarks                                                 │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  BENCHMARK A                      │  BENCHMARK B                        │ │
│  │  16S paired (binary)              │  16S paired (text)                  │ │
│  │  ─────────────────────────────────┼──────────────────────────────────   │ │
│  │  Format: binary                   │  Format: text (gzipped)             │ │
│  │  OPTIMIZE_MEMORY: No              │  OPTIMIZE_MEMORY: No                │ │
│  │  Query: paired-end                │  Query: paired-end                  │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  PERFORMANCE COMPARISON                                                 │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Metric           A              B              Diff        Winner      │ │
│  │  ────────────────────────────────────────────────────────────────────   │ │
│  │  Wall Time        45.2s          89.7s          +99%        ← A         │ │
│  │  Peak Memory      64.8 MB        128.4 MB       +98%        ← A         │ │
│  │  Reads/sec        221.2          111.5          -50%        ← A         │ │
│  │                                                                         │ │
│  │  Performance Winner: BENCHMARK A                                        │ │
│  │    • 2.0x faster                                                        │ │
│  │    • Uses 50% less memory                                               │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  RESULTS COMPARISON                                                     │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Metric           A              B                                      │ │
│  │  ────────────────────────────────────────────────────────────────────   │ │
│  │  Assigned         9,847 (98.5%)  9,847 (98.5%)                          │ │
│  │  Unique Taxa      142            142                                    │ │
│  │                                                                         │ │
│  │  ✓ ASSIGNMENTS ARE IDENTICAL                                            │ │
│  │    Hash: d41d8cd98f00b204e9800998ecf8427e                               │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  VERDICT                                                                │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Results are IDENTICAL. Binary format (A) is 2x faster with 50%         │ │
│  │  less memory. Recommend using binary format for this workload.          │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Esc close                                                              │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Comparison View - Different Results

When two benchmarks have different assignments (warning state):

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Benchmark Comparison ──────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Comparing 2 benchmarks                                                 │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  BENCHMARK A                      │  BENCHMARK B                        │ │
│  │  16S paired (binary)              │  16S paired (binary+optmem)         │ │
│  │  ─────────────────────────────────┼──────────────────────────────────   │ │
│  │  Format: binary                   │  Format: binary                     │ │
│  │  OPTIMIZE_MEMORY: No              │  OPTIMIZE_MEMORY: Yes               │ │
│  │  Query: paired-end                │  Query: paired-end                  │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  PERFORMANCE COMPARISON                                                 │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Metric           A              B              Diff        Winner      │ │
│  │  ────────────────────────────────────────────────────────────────────   │ │
│  │  Wall Time        45.2s          38.1s          -16%        B →         │ │
│  │  Peak Memory      64.8 MB        42.3 MB        -35%        B →         │ │
│  │  Reads/sec        221.2          262.5          +19%        B →         │ │
│  │                                                                         │ │
│  │  Performance Winner: BENCHMARK B                                        │ │
│  │    • 1.2x faster                                                        │ │
│  │    • Uses 35% less memory                                               │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  RESULTS COMPARISON                                                     │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  Metric           A              B                                      │ │
│  │  ────────────────────────────────────────────────────────────────────   │ │
│  │  Assigned         9,847 (98.5%)  9,842 (98.4%)                          │ │
│  │  Unique Taxa      142            141                                    │ │
│  │                                                                         │ │
│  │  ⚠ ASSIGNMENTS DIFFER                                                   │ │
│  │    Hash A: d41d8cd98f00b204e9800998ecf8427e                             │ │
│  │    Hash B: 7d793037a0760186574b0282f2f435e7                             │ │
│  │                                                                         │ │
│  │    Differences detected (5 reads with different assignments)            │ │
│  │    This may be due to floating-point precision differences.             │ │
│  │                                                                         │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │  VERDICT                                                                │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  ⚠ Results DIFFER slightly. OPTIMIZE_MEMORY build (B) is faster and     │ │
│  │  uses less memory, but produces slightly different assignments due      │ │
│  │  to reduced floating-point precision. Review differences carefully.     │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Esc close                                                              │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

### Help Overlay

Updated help showing all new keybindings:

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ Help ──────────────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  tronko-bench TUI                                                       │ │
│  │  ═══════════════════════════════════════════════════════════════════    │ │
│  │                                                                         │ │
│  │  NAVIGATION                                                             │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Tab / Shift+Tab     Switch between tabs                                │ │
│  │  j / k / ↑ / ↓       Navigate list items                                │ │
│  │  /                   Fuzzy search in current list                       │ │
│  │  Esc                 Close overlay / Cancel                             │ │
│  │                                                                         │ │
│  │  DATABASES TAB                                                          │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  c                   Convert selected database to binary                │ │
│  │  b                   Start new benchmark wizard                         │ │
│  │  g                   Generate benchmark.json config                     │ │
│  │  r                   Rescan directory                                   │ │
│  │                                                                         │ │
│  │  BENCHMARKS TAB                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  b                   Start new benchmark wizard                         │ │
│  │  Space               Select/deselect benchmark (max 2)                  │ │
│  │  c                   Compare selected benchmarks (when 2 selected)      │ │
│  │  d                   Delete selected benchmark                          │ │
│  │                                                                         │ │
│  │  WIZARD                                                                 │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Enter               Select item / Confirm                              │ │
│  │  ← / Backspace       Go back to previous step                           │ │
│  │  Esc                 Cancel wizard                                      │ │
│  │                                                                         │ │
│  │  OTHER                                                                  │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ?                   Show this help                                     │ │
│  │  q                   Quit                                               │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Press any key to close                                                 │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

---

## UI Design

### Tab Structure

**Keep existing tabs** and add a Benchmarks tab:

```
┌──────────────────────────────────────────────────────────────────────┐
│ [Databases (5)]  [Queries (12)]  [MSA Files (3)]  [Benchmarks (7)]   │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  (Content varies by tab)                                             │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ Status bar with context-sensitive hints                              │
└──────────────────────────────────────────────────────────────────────┘
```

### Benchmarks Tab - List View

```
┌─ Benchmarks ─────────────────────────────┬─ Details ─────────────────┐
│ [ ] 16S_Bacteria paired (binary) 2m ago  │ Name: 16S_Bacteria paired │
│ [*] 16S_Bacteria paired (text)   5m ago  │                           │
│ [*] ITS_Fungi single (binary)    1h ago  │ Database:                 │
│ [ ] 16S_Bacteria paired (optmem) 2h ago  │   16S_Bacteria/           │
│                                          │   reference_tree.trkb     │
│                                          │   Format: binary          │
│                                          │   Trees: 156              │
│                                          │                           │
│                                          │ Queries:                  │
│                                          │   Project: example_proj   │
│                                          │   Marker: 16S             │
│                                          │   Mode: paired-end        │
│                                          │                           │
│                                          │ Performance:              │
│                                          │   Time: 45.2s             │
│                                          │   Peak Memory: 64.8 MB    │
│                                          │                           │
│                                          │ Results:                  │
│                                          │   Assigned: 9,847 (98.5%) │
├──────────────────────────────────────────┴───────────────────────────┤
│ Space:select  b:new  Enter:details  d:delete  2 selected → c:compare │
└──────────────────────────────────────────────────────────────────────┘
```

### Benchmark Wizard (6 Steps)

**Step 1: Select Reference Database**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 1/6: Select Reference Database                                  │
│ ────────────────────────────────────────────────────────────────────│
│ [1.DB] → 2.Project → 3.Marker → 4.Query → 5.Binary → 6.Confirm       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ > 16S_Bacteria/reference_tree.trkb      156 trees  5.78 GiB  [BIN]   │
│   16S_Bacteria/reference_tree.trkb.gz   156 trees  352 MiB   [BIN.gz]│
│   16S_Bacteria/reference_tree.txt.gz    156 trees  757 MiB   [TXT.gz]│
│   single_tree/reference_tree.txt        1 tree     38 MiB    [TXT]   │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ ↑↓:navigate  Enter:select  /:search  Esc:cancel                      │
└──────────────────────────────────────────────────────────────────────┘
```

**Step 2: Select Project**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 2/6: Select Project                                             │
│ ────────────────────────────────────────────────────────────────────│
│ 1.DB → [2.Project] → 3.Marker → 4.Query → 5.Binary → 6.Confirm       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ > 16S (markers: 16S)                                                 │
│   43f82f91-3857-443d-8c15-e3693be268d6 (markers: 16S, ITS)           │
│   c8939113-7f18-4103-83d8-50a72f8b4297 (markers: 16S)                │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ ↑↓:navigate  Enter:select  ←:back  Esc:cancel                        │
└──────────────────────────────────────────────────────────────────────┘
```

**Step 3: Select Marker**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 3/6: Select Marker                                              │
│ ────────────────────────────────────────────────────────────────────│
│ 1.DB → 2.Project → [3.Marker] → 4.Query → 5.Binary → 6.Confirm       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ > 16S (paired, unpaired_F, unpaired_R)                               │
│   ITS (paired)                                                       │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ ↑↓:navigate  Enter:select  ←:back  Esc:cancel                        │
└──────────────────────────────────────────────────────────────────────┘
```

**Step 4: Select Query Mode**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 4/6: Select Query Type                                          │
│ ────────────────────────────────────────────────────────────────────│
│ 1.DB → 2.Project → 3.Marker → [4.Query] → 5.Binary → 6.Confirm       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ > Paired-end reads                                                   │
│   Unpaired (Forward) reads                                           │
│   Unpaired (Reverse) reads                                           │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ ↑↓:navigate  Enter:select  ←:back  Esc:cancel                        │
└──────────────────────────────────────────────────────────────────────┘
```

**Step 5: Select Binary**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 5/6: Select tronko-assign Binary                                │
│ ────────────────────────────────────────────────────────────────────│
│ 1.DB → 2.Project → 3.Marker → 4.Query → [5.Binary] → 6.Confirm       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ > tronko-assign (standard, double precision)                         │
│   tronko-assign-optimized (OPTIMIZE_MEMORY, float precision)         │
│                                                                      │
│ Note: Build optimized version with: make OPTIMIZE_MEMORY=1           │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ ↑↓:navigate  Enter:select  ←:back  Esc:cancel                        │
└──────────────────────────────────────────────────────────────────────┘
```

**Step 6: Confirm and Run**
```
┌─ New Benchmark ──────────────────────────────────────────────────────┐
│ Step 6/6: Confirm & Run                                              │
│ ────────────────────────────────────────────────────────────────────│
│ 1.DB → 2.Project → 3.Marker → 4.Query → 5.Binary → [6.Confirm]       │
├──────────────────────────────────────────────────────────────────────┤
│                                                                      │
│ Database:    16S_Bacteria/reference_tree.trkb                        │
│              Format: binary, Trees: 156, Size: 5.78 GiB              │
│                                                                      │
│ Queries:     16S/paired/                                             │
│              Project: 16S, Mode: paired-end                          │
│                                                                      │
│ Binary:      tronko-assign (standard)                                │
│                                                                      │
│ Command:                                                             │
│   tronko-assign -r -f reference_tree.trkb -a 16S.fasta \             │
│     -p -1 paired_F.fasta.zst -2 paired_R.fasta.zst -o output.tsv     │
│                                                                      │
│                    [ Run Benchmark ]    [ Cancel ]                   │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ Enter:run  ←:back  Esc:cancel                                        │
└──────────────────────────────────────────────────────────────────────┘
```

### Comparison View

When two benchmarks are selected via spacebar:

```
┌─ Benchmark Comparison ───────────────────────────────────────────────┐
│                                                                      │
│  Benchmark A                         vs   Benchmark B                │
│  16S_Bacteria paired (binary)             16S_Bacteria paired (text) │
│  ────────────────────────────────────────────────────────────────── │
│                                                                      │
│  PERFORMANCE                                                         │
│  ───────────                                                         │
│  Metric          Benchmark A         Benchmark B         Winner      │
│  ─────────────────────────────────────────────────────────────────  │
│  Wall Time       45.2s               89.7s (+99%)        ← A         │
│  Peak Memory     64.8 MB             128.4 MB (+98%)     ← A         │
│  Reads/sec       221.2               111.5 (-50%)        ← A         │
│                                                                      │
│  Winner: Benchmark A (faster, uses less memory)                      │
│                                                                      │
│  RESULTS                                                             │
│  ───────                                                             │
│  Assigned        9,847 (98.5%)       9,847 (98.5%)                   │
│  Unique Taxa     142                 142                             │
│                                                                      │
│  Assignment Comparison:                                              │
│    ✓ Assignments are IDENTICAL (hash match)                          │
│                                                                      │
│  VERDICT: Results identical. Binary format is 2x faster.             │
│                                                                      │
├──────────────────────────────────────────────────────────────────────┤
│ Esc:close                                                            │
└──────────────────────────────────────────────────────────────────────┘
```

---

## Implementation Phases

### Phase 1: Data Model & Storage Layer

**New files:**
- `src/benchmark/mod.rs` - Module exports
- `src/benchmark/config.rs` - BenchmarkConfig, DatabaseConfig, QueryConfig, ExecutionParams
- `src/benchmark/result.rs` - BenchmarkResult, PerformanceMetrics, AssignmentStats
- `src/benchmark/storage.rs` - BenchmarkStorage for persistence

**Dependencies to add:**
```toml
uuid = { version = "1", features = ["v4", "serde"] }
chrono = { version = "0.4", features = ["serde"] }
directories = "5"
md5 = "0.7"
```

**Key implementation:**

```rust
// src/benchmark/storage.rs
use directories::ProjectDirs;

pub struct BenchmarkStorage {
    base_dir: PathBuf,
}

impl BenchmarkStorage {
    pub fn new() -> Result<Self> {
        let dirs = ProjectDirs::from("", "", "tronko-bench")
            .ok_or_else(|| anyhow::anyhow!("Could not determine config directory"))?;
        let base_dir = dirs.data_dir().to_path_buf();
        std::fs::create_dir_all(&base_dir)?;
        std::fs::create_dir_all(base_dir.join("benchmarks"))?;
        Ok(Self { base_dir })
    }

    pub fn list_benchmarks(&self) -> Result<Vec<BenchmarkSummary>>;
    pub fn get_benchmark(&self, id: &str) -> Result<(BenchmarkConfig, Option<BenchmarkResult>)>;
    pub fn save_config(&self, config: &BenchmarkConfig) -> Result<PathBuf>;
    pub fn save_result(&self, config_id: &str, result: &BenchmarkResult) -> Result<()>;
    pub fn delete_benchmark(&self, id: &str) -> Result<()>;
    pub fn get_output_paths(&self, id: &str) -> (PathBuf, PathBuf, PathBuf);  // output, metrics, stderr
}
```

**Success Criteria:**
- [ ] Can create and serialize BenchmarkConfig to JSON
- [ ] Storage directory created at `~/.local/share/tronko-bench/` (Linux) or equivalent
- [ ] Index file maintained for quick listing
- [ ] Code compiles: `cargo build -p tronko-bench`

---

### Phase 2: Query Directory Discovery

**New file:**
- `src/benchmark/discovery.rs` - Scan for projects, markers, and query files

**Implementation:**

```rust
/// A discovered project with query files
#[derive(Debug, Clone)]
pub struct DiscoveredProject {
    pub name: String,
    pub path: PathBuf,
    pub markers: Vec<DiscoveredMarker>,
}

#[derive(Debug, Clone)]
pub struct DiscoveredMarker {
    pub name: String,
    pub path: PathBuf,
    pub has_paired: bool,
    pub has_unpaired_f: bool,
    pub has_unpaired_r: bool,
    pub paired_f: Option<PathBuf>,
    pub paired_r: Option<PathBuf>,
    pub unpaired_f: Option<PathBuf>,
    pub unpaired_r: Option<PathBuf>,
}

/// Scan for query directories in example_project/
pub fn discover_projects(base_dir: &Path) -> Vec<DiscoveredProject> {
    // Scan base_dir for project directories
    // Within each project, find marker directories (16S, ITS, etc.)
    // Within each marker, look for paired/, unpaired_F/, unpaired_R/
    // Detect files: *.fasta, *.fasta.zst, *.fasta.gz, *.fastq, etc.
}

/// Discover available tronko-assign binaries
pub fn discover_binaries(scan_dir: &Path) -> Vec<TronkoBinary> {
    // Look for:
    //   - ./tronko-assign/tronko-assign (standard)
    //   - ./tronko-assign/tronko-assign-optimized (OPTIMIZE_MEMORY)
    //   - System PATH
}

#[derive(Debug, Clone)]
pub struct TronkoBinary {
    pub path: PathBuf,
    pub optimize_memory: bool,
    pub label: String,
}
```

**Success Criteria:**
- [ ] Correctly discovers projects in example_project/
- [ ] Detects paired vs unpaired directories
- [ ] Finds .fasta.zst, .fasta.gz, .fasta files
- [ ] Discovers both standard and optimized tronko-assign binaries

---

### Phase 3: Enhanced Database Display

**Modify:**
- `src/main.rs` - Add `context` and `size_bytes` fields to DatabaseInfo
- `src/tui/ui/list.rs` - Update rendering to show context

**DatabaseInfo changes:**
```rust
pub struct DatabaseInfo {
    pub name: String,
    pub context: String,      // Parent directory name
    pub format: String,       // "text", "binary"
    pub compressed: bool,     // Whether .gz
    pub num_trees: String,
    pub size: String,         // Human-readable
    pub size_bytes: u64,      // Raw bytes
    pub path: PathBuf,
}
```

**Display format:**
```
16S_Bacteria/reference_tree.trkb      156 trees   5.78 GiB  [BIN]
16S_Bacteria/reference_tree.trkb.gz   156 trees   352 MiB   [BIN.gz]
single_tree/reference_tree.txt        1 tree      38 MiB    [TXT]
```

**Success Criteria:**
- [ ] Database list shows parent directory for context
- [ ] Format shows binary/text and compression status
- [ ] Code compiles

---

### Phase 4: App State Extension

**Modify:**
- `src/tui/app.rs` - Add Benchmarks tab, modes, wizard state
- `src/tui/message.rs` - Add new message variants

**New App state:**

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Databases,
    Queries,
    MsaFiles,
    Benchmarks,  // NEW
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AppMode {
    Normal,
    Search,
    Help,
    BenchmarkWizard,
    Comparison,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WizardStep {
    SelectDatabase,
    SelectProject,
    SelectMarker,
    SelectQueryMode,
    SelectBinary,
    Confirm,
}

pub struct WizardState {
    pub step: WizardStep,
    pub list_index: usize,
    pub selected_db_index: Option<usize>,
    pub discovered_projects: Vec<DiscoveredProject>,
    pub selected_project_index: Option<usize>,
    pub selected_marker_index: Option<usize>,
    pub selected_query_mode: Option<ReadMode>,
    pub tronko_binaries: Vec<TronkoBinary>,
    pub selected_binary_index: Option<usize>,
}

pub struct App {
    // ... existing fields ...

    // NEW
    pub mode: AppMode,
    pub storage: BenchmarkStorage,
    pub benchmarks: Vec<BenchmarkSummary>,
    pub wizard: WizardState,
    pub selected_benchmarks: Vec<String>,  // IDs, max 2
}
```

**New messages:**
```rust
pub enum Message {
    // ... existing ...

    // Wizard
    StartBenchmarkWizard,
    WizardSelectItem,
    WizardNext,
    WizardBack,
    WizardCancel,

    // Benchmarks
    ToggleBenchmarkSelection,
    DeleteBenchmark,
    EnterComparison,
    ExitComparison,
    RefreshBenchmarks,
}
```

**Success Criteria:**
- [ ] Benchmarks tab appears
- [ ] Mode transitions work
- [ ] Code compiles

---

### Phase 5: Benchmark Tab UI

**New/modify files:**
- `src/tui/ui/list.rs` - Add benchmark list rendering
- `src/tui/ui/detail.rs` - Add benchmark detail rendering

**Implementation:**

```rust
// In list.rs
fn render_benchmark_items(app: &App) -> Vec<ListItem<'static>> {
    app.benchmarks.iter().enumerate().map(|(idx, bench)| {
        let selected = app.selected_benchmarks.contains(&bench.id);
        let indicator = if selected { "[*]" } else { "[ ]" };

        // Format: "[ ] name (format) relative_time"
        ListItem::new(Line::from(vec![
            Span::styled(indicator, ...),
            Span::raw(format!(" {} ", bench.name)),
            Span::styled(format!("{:.1}s", bench.wall_time_secs), ...),
            Span::styled(format!("{} ", format_size(bench.peak_memory, BINARY)), ...),
            Span::styled(if bench.success { "[OK]" } else { "[FAIL]" }, ...),
        ]))
    }).collect()
}
```

**Success Criteria:**
- [ ] Benchmark list displays with selection indicators
- [ ] Detail panel shows full benchmark info
- [ ] Can select/deselect with spacebar (max 2)

---

### Phase 6: Benchmark Wizard UI

**New file:**
- `src/tui/ui/wizard.rs` - Multi-step wizard overlay

**Implementation:**
- Progress indicator showing current step
- Appropriate list for each step
- Navigation: Enter to select, ← to go back, Esc to cancel

**Success Criteria:**
- [ ] Wizard opens with 'b' key
- [ ] Can navigate through all 6 steps
- [ ] Selections preserved when going back
- [ ] Cancel returns to normal mode

---

### Phase 7: Comparison UI

**New file:**
- `src/tui/ui/comparison.rs` - Side-by-side comparison overlay

**Implementation:**
- Two-column layout with benchmark details
- Performance comparison with percentage differences
- Winner indicators
- Result comparison (hash match or diff)
- Summary verdict

**Success Criteria:**
- [ ] Comparison view opens when 2 benchmarks selected
- [ ] Shows performance differences with percentages
- [ ] Indicates which benchmark wins each metric
- [ ] Shows whether assignments are identical

---

### Phase 8: Event Handling

**Modify:**
- `src/tui/event.rs` - Mode-aware key handling
- `src/tui/mod.rs` - Handle new messages

**Key bindings by mode:**

| Mode | Key | Action |
|------|-----|--------|
| Normal (Benchmarks tab) | `b` | Start wizard |
| Normal (Benchmarks tab) | `Space` | Toggle selection |
| Normal (Benchmarks tab) | `d` | Delete benchmark |
| Normal (2 selected) | `c` or auto | Enter comparison |
| Wizard | `Enter` | Select/confirm |
| Wizard | `←` / `Backspace` | Go back |
| Wizard | `Esc` | Cancel |
| Comparison | `Esc` | Exit comparison |
| All | `?` | Help |
| All | `q` | Quit |

**Success Criteria:**
- [ ] All keybindings work in appropriate modes
- [ ] Mode transitions correct

---

### Phase 9: Benchmark Execution

**New file:**
- `src/benchmark/runner.rs` - Execute tronko-assign and capture results

**Implementation:**

```rust
pub struct BenchmarkRunner;

impl BenchmarkRunner {
    pub fn run(
        config: &BenchmarkConfig,
        output_path: &Path,
        metrics_path: &Path,
        stderr_path: &Path,
    ) -> Result<BenchmarkResult> {
        let mut cmd = Command::new(&config.params.tronko_assign_path);

        cmd.arg("-r")
           .arg("-f").arg(&config.database.reference_path)
           .arg("-a").arg(&config.database.fasta_path)
           .arg("-o").arg(output_path)
           .arg("-R")  // Resource monitoring
           .arg("-T")  // Timing
           .arg("--tsv-log").arg(metrics_path)
           .arg("-C").arg(config.params.cores.to_string())
           .arg("-L").arg(config.params.batch_size.to_string());

        if config.params.use_needleman_wunsch {
            cmd.arg("-w");
        }

        match &config.queries.files {
            QueryFiles::Single { path } => {
                cmd.arg("-s").arg("-g").arg(path);
            }
            QueryFiles::Paired { forward, reverse } => {
                cmd.arg("-p")
                   .arg("-1").arg(forward)
                   .arg("-2").arg(reverse);
            }
        }

        let start = Instant::now();
        let output = cmd.output()?;
        let duration = start.elapsed();

        // Parse results...
    }
}
```

**New file:**
- `src/benchmark/parser.rs` - Parse TSV metrics and assignment output

```rust
pub fn parse_metrics_tsv(path: &Path) -> Result<PerformanceMetrics> {
    // Parse --tsv-log output to extract peak RSS, etc.
}

pub fn parse_assignments(path: &Path) -> Result<AssignmentStats> {
    // Count total, assigned, unique taxa
    // Compute MD5 hash of sorted assignments
}
```

**Success Criteria:**
- [ ] tronko-assign command built correctly
- [ ] Output, metrics, stderr captured to files
- [ ] Metrics parsed from TSV log
- [ ] Assignment stats computed with hash

---

### Phase 10: Help & Status Updates

**Modify:**
- `src/tui/ui/help.rs` - Add new keybindings
- `src/tui/ui/status.rs` - Context-aware hints

**Success Criteria:**
- [ ] Help shows all keyboard shortcuts
- [ ] Status bar shows relevant hints per tab/mode

---

## File Summary

**New files:**
```
src/benchmark/
├── mod.rs
├── config.rs
├── result.rs
├── storage.rs
├── discovery.rs
├── runner.rs
├── parser.rs
└── comparison.rs

src/tui/ui/
├── wizard.rs
└── comparison.rs
```

**Modified files:**
```
Cargo.toml           (add uuid, chrono, directories, md5)
src/main.rs          (add benchmark module, update DatabaseInfo)
src/tui/
├── app.rs           (add modes, wizard state, benchmark list)
├── message.rs       (add new messages)
├── event.rs         (mode-aware key handling)
├── mod.rs           (handle new messages)
└── ui/
    ├── mod.rs       (mode-based rendering)
    ├── list.rs      (benchmark list + multi-select)
    ├── detail.rs    (benchmark details)
    ├── status.rs    (context hints)
    └── help.rs      (new keybindings)
```

---

## Testing Strategy

### Unit Tests
- `storage.rs`: Test save/load/delete operations
- `parser.rs`: Test TSV parsing with sample files
- `discovery.rs`: Test project/marker detection

### Integration Tests
- Full wizard flow → execution → result storage
- Comparison of two benchmarks with known outputs

### Manual Testing Checklist
- [ ] Create benchmark with binary database
- [ ] Create benchmark with gzipped text database
- [ ] Create benchmark with paired-end reads
- [ ] Create benchmark with unpaired reads
- [ ] Use OPTIMIZE_MEMORY build
- [ ] Cancel wizard at each step
- [ ] Compare two successful benchmarks
- [ ] Compare benchmarks with identical results
- [ ] Delete a benchmark
- [ ] Restart app, verify benchmarks persist

---

## Comparison: Original Plans

| Aspect | Original Plan | TUI Plan | Unified |
|--------|---------------|----------|---------|
| Storage location | `~/.tronko-bench/` | Scan directory | `~/.tronko-bench/` |
| Tab structure | Replace old tabs | Add Benchmarks tab | Add Benchmarks tab |
| Wizard steps | 4 steps | 6 steps (granular) | 6 steps |
| Binary selection | Checkbox in config | Dedicated step | Dedicated step |
| Result comparison | Detailed score stats | MD5 hash only | Hash + verdict |
| Progress display | Live output | Blocking | Blocking (v1) |
| File structure | Modular | Single file | Modular |

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| tronko-assign not found | Scan common locations, show clear error |
| Benchmark takes too long | Show status message, allow Ctrl+C to cancel |
| Large output files | Stream parse, don't load entire file |
| Concurrent modification | Single-user assumed; file locking in future |
| OPTIMIZE_MEMORY binary missing | Clear note in wizard with build instructions |

---

## Success Metrics

1. **Usability**: User can create and run a benchmark in < 1 minute via wizard
2. **Accuracy**: Comparison correctly identifies identical vs different assignments
3. **Persistence**: Benchmarks survive app restarts
4. **Performance**: TUI remains responsive during benchmark (uses status updates)
