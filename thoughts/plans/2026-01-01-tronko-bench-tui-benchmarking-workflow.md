# Tronko-Bench TUI Benchmarking Workflow Implementation Plan

## Overview

Extend the tronko-bench TUI to support a complete benchmarking workflow: running new benchmarks with user-selected configurations, persisting benchmark results, viewing historical benchmarks, and comparing up to two benchmarks side-by-side.

## Current State Analysis

### Existing TUI Structure
- **Tabs**: Databases, Queries, MSA Files (`tui/app.rs:1-115`)
- **UI Components**: tabs, list, detail, status, help (`tui/ui/mod.rs:1-84`)
- **Message pattern**: Event-driven with `Message` enum (`tui/message.rs:1-27`)
- **Scanning**: Uses `scan_databases()` and `scan_query_files()` from `main.rs`

### Directory Structure Patterns
- **Reference databases**: `example_datasets/<marker>/reference_tree.*` (formats: `.txt`, `.txt.gz`, `.trkb`, `.trkb.gz`)
- **Query files**: `example_project/<project_uuid>/<marker>/paired/` or `unpaired_F/` or `unpaired_R/`
  - Paired: `paired_F.fasta.zst`, `paired_R.fasta.zst`
  - Unpaired: `unpaired_F.fasta.zst`, `unpaired_R.fasta.zst`

### tronko-assign Options
- `-s` single-end, `-p` paired-end
- `-g` single-end file path
- `-1`/`-2` paired-end file paths
- `-R` resource monitoring, `-T` timing, `--tsv-log` memory export
- `OPTIMIZE_MEMORY` is compile-time (needs two binaries or rebuild)

## Desired End State

A TUI with the following capabilities:

1. **Benchmark Creation Wizard**: Multi-step wizard to configure and run a benchmark
2. **Benchmark History**: Persistent list of all past benchmarks with full details
3. **Comparison Mode**: Select 2 benchmarks via spacebar, view side-by-side comparison

### Verification Criteria
- User can run benchmarks via TUI and see results persist
- User can navigate historical benchmarks and see performance details
- User can select 2 benchmarks and compare speed, memory, and result differences

## What We're NOT Doing

- Automatic rebuilding of tronko-assign with/without OPTIMIZE_MEMORY (user must pre-build both)
- Result assignment diffing at the individual read level (only summary differences)
- Real-time progress updates during benchmark execution (blocking operation)
- Integration with external benchmarking tools (hyperfine, etc.)

## Implementation Approach

We'll extend the existing TUI with new tabs and modes while preserving the current architecture. Key additions:
- **Benchmarks tab**: List historical benchmarks
- **Benchmark wizard**: Multi-step overlay for creating benchmarks
- **Comparison view**: Side-by-side when 2 benchmarks selected
- **Persistence**: JSON file storing benchmark results

## Phase 1: Data Structures and Persistence

### Overview
Define the data structures for benchmarks and implement JSON persistence.

### Changes Required:

#### 1. Create benchmark types module
**File**: `rust/crates/tronko-bench/src/benchmark.rs` (new file)
**Changes**: Define benchmark configuration and result types

```rust
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::time::Duration;

/// Configuration for a benchmark run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    /// Unique identifier
    pub id: String,
    /// Human-readable name
    pub name: String,
    /// Path to reference database
    pub reference_db_path: PathBuf,
    /// Database format (text, text (gzipped), binary (.trkb), binary (gzipped))
    pub db_format: String,
    /// Number of trees in database
    pub num_trees: u32,
    /// Database file size
    pub db_size_bytes: u64,
    /// Path to FASTA file for BWA
    pub fasta_path: PathBuf,
    /// Query file paths (single for single-end, two for paired-end)
    pub query_paths: Vec<PathBuf>,
    /// Is this paired-end?
    pub paired_end: bool,
    /// Marker name (e.g., "16S_Bacteria")
    pub marker: String,
    /// Project name (directory name from example_project)
    pub project: String,
    /// Whether OPTIMIZE_MEMORY build was used
    pub optimize_memory: bool,
    /// Path to tronko-assign binary used
    pub tronko_assign_path: PathBuf,
}

/// Results from a benchmark run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    /// Configuration used
    pub config: BenchmarkConfig,
    /// When the benchmark was run
    pub timestamp: String,
    /// Total wall-clock time
    pub duration_secs: f64,
    /// Peak memory usage in bytes (if available)
    pub peak_memory_bytes: Option<u64>,
    /// Number of reads processed
    pub reads_processed: Option<u64>,
    /// Path to output file
    pub output_path: PathBuf,
    /// Path to TSV memory log (if --tsv-log was used)
    pub tsv_log_path: Option<PathBuf>,
    /// Whether the run succeeded
    pub success: bool,
    /// Error message if failed
    pub error_message: Option<String>,
    /// Output file size
    pub output_size_bytes: Option<u64>,
    /// Summary statistics from output
    pub output_summary: Option<OutputSummary>,
}

/// Summary of assignment output for comparison
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputSummary {
    /// Total lines in output
    pub total_lines: u64,
    /// Lines with assignments (non-empty taxonomy)
    pub assigned_lines: u64,
    /// Unique taxonomic assignments
    pub unique_taxa: u64,
    /// MD5 hash of sorted assignments for quick comparison
    pub assignments_hash: String,
}

/// Collection of all benchmarks
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct BenchmarkStore {
    pub benchmarks: Vec<BenchmarkResult>,
}

impl BenchmarkStore {
    pub fn load(path: &PathBuf) -> anyhow::Result<Self> {
        if path.exists() {
            let content = std::fs::read_to_string(path)?;
            Ok(serde_json::from_str(&content)?)
        } else {
            Ok(Self::default())
        }
    }

    pub fn save(&self, path: &PathBuf) -> anyhow::Result<()> {
        let content = serde_json::to_string_pretty(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }

    pub fn add(&mut self, result: BenchmarkResult) {
        self.benchmarks.push(result);
    }
}
```

#### 2. Update main.rs to expose the module
**File**: `rust/crates/tronko-bench/src/main.rs`
**Changes**: Add module declaration at the top

```rust
mod benchmark;
mod tui;
// ... rest of existing code
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles without errors: `cd rust && cargo build -p tronko-bench`
- [ ] No clippy warnings: `cd rust && cargo clippy -p tronko-bench`

#### Manual Verification:
- [ ] N/A (data structures only)

---

## Phase 2: Query File Discovery

### Overview
Add functionality to discover and categorize query files from example_project directory structure.

### Changes Required:

#### 1. Add query discovery types and functions
**File**: `rust/crates/tronko-bench/src/benchmark.rs`
**Changes**: Add query discovery types after the existing types

```rust
/// A discovered project with query files
#[derive(Debug, Clone)]
pub struct DiscoveredProject {
    /// Project directory name (UUID or identifier)
    pub name: String,
    /// Full path to project directory
    pub path: PathBuf,
    /// Markers found in this project
    pub markers: Vec<DiscoveredMarker>,
}

/// A marker within a project (e.g., 16S, COI)
#[derive(Debug, Clone)]
pub struct DiscoveredMarker {
    /// Marker name (directory name)
    pub name: String,
    /// Full path to marker directory
    pub path: PathBuf,
    /// Whether paired-end files are available
    pub has_paired: bool,
    /// Whether unpaired forward files are available
    pub has_unpaired_f: bool,
    /// Whether unpaired reverse files are available
    pub has_unpaired_r: bool,
    /// Paired forward file path (if exists)
    pub paired_f: Option<PathBuf>,
    /// Paired reverse file path (if exists)
    pub paired_r: Option<PathBuf>,
    /// Unpaired forward file path (if exists)
    pub unpaired_f: Option<PathBuf>,
    /// Unpaired reverse file path (if exists)
    pub unpaired_r: Option<PathBuf>,
}

/// Query mode selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryMode {
    Paired,
    UnpairedF,
    UnpairedR,
}

impl QueryMode {
    pub fn label(&self) -> &'static str {
        match self {
            QueryMode::Paired => "Paired-end",
            QueryMode::UnpairedF => "Unpaired (Forward)",
            QueryMode::UnpairedR => "Unpaired (Reverse)",
        }
    }
}

/// Scan example_project directory for available query files
pub fn discover_projects(base_dir: &std::path::Path) -> Vec<DiscoveredProject> {
    let mut projects = Vec::new();

    let Ok(entries) = std::fs::read_dir(base_dir) else {
        return projects;
    };

    for entry in entries.filter_map(|e| e.ok()) {
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let project_name = path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        let markers = discover_markers(&path);
        if !markers.is_empty() {
            projects.push(DiscoveredProject {
                name: project_name,
                path,
                markers,
            });
        }
    }

    projects.sort_by(|a, b| a.name.cmp(&b.name));
    projects
}

fn discover_markers(project_dir: &std::path::Path) -> Vec<DiscoveredMarker> {
    let mut markers = Vec::new();

    let Ok(entries) = std::fs::read_dir(project_dir) else {
        return markers;
    };

    for entry in entries.filter_map(|e| e.ok()) {
        let path = entry.path();
        if !path.is_dir() {
            continue;
        }

        let marker_name = path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Look for paired/unpaired subdirectories
        let paired_dir = path.join("paired");
        let unpaired_f_dir = path.join("unpaired_F");
        let unpaired_r_dir = path.join("unpaired_R");

        // Find fasta files (support .fasta, .fasta.zst, .fasta.gz)
        let paired_f = find_query_file(&paired_dir, "paired_F", "paired_f");
        let paired_r = find_query_file(&paired_dir, "paired_R", "paired_r");
        let unpaired_f = find_query_file(&unpaired_f_dir, "unpaired_F", "unpaired_f");
        let unpaired_r = find_query_file(&unpaired_r_dir, "unpaired_R", "unpaired_r");

        let has_paired = paired_f.is_some() && paired_r.is_some();
        let has_unpaired_f = unpaired_f.is_some();
        let has_unpaired_r = unpaired_r.is_some();

        if has_paired || has_unpaired_f || has_unpaired_r {
            markers.push(DiscoveredMarker {
                name: marker_name,
                path,
                has_paired,
                has_unpaired_f,
                has_unpaired_r,
                paired_f,
                paired_r,
                unpaired_f,
                unpaired_r,
            });
        }
    }

    markers.sort_by(|a, b| a.name.cmp(&b.name));
    markers
}

fn find_query_file(dir: &std::path::Path, name1: &str, name2: &str) -> Option<PathBuf> {
    if !dir.exists() {
        return None;
    }

    let extensions = [".fasta.zst", ".fasta.gz", ".fasta", ".fa.zst", ".fa.gz", ".fa"];

    for ext in extensions {
        for name in [name1, name2] {
            let path = dir.join(format!("{}{}", name, ext));
            if path.exists() {
                return Some(path);
            }
        }
    }

    None
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`
- [ ] Unit tests pass: `cd rust && cargo test -p tronko-bench`

#### Manual Verification:
- [ ] N/A (discovery functions only)

---

## Phase 3: Enhanced Reference Database Display

### Overview
Update the database list to show contextual information (parent directory, format, compression, tree count, size).

### Changes Required:

#### 1. Update DatabaseInfo to include parent directory context
**File**: `rust/crates/tronko-bench/src/main.rs`
**Changes**: Add `context_path` field to `DatabaseInfo` struct and update `scan_databases()`

```rust
/// Information about a reference database
#[derive(Debug, Clone, Serialize, Deserialize, Tabled)]
pub struct DatabaseInfo {
    #[tabled(rename = "Name")]
    pub name: String,
    #[tabled(rename = "Context")]
    pub context: String,  // NEW: parent directory for context
    #[tabled(rename = "Format")]
    pub format: String,
    #[tabled(rename = "Compressed")]
    pub compressed: bool,  // NEW: whether gzipped
    #[tabled(rename = "Trees")]
    pub num_trees: String,
    #[tabled(rename = "Size")]
    pub size: String,
    #[tabled(skip)]
    pub size_bytes: u64,  // NEW: raw size for calculations
    #[tabled(rename = "Path")]
    #[tabled(display_with = "truncate_path")]
    pub path: PathBuf,
}
```

Update `scan_databases()` to populate the new fields:

```rust
/// Scan a directory for databases
pub fn scan_databases(dir: &Path, max_depth: usize) -> Vec<DatabaseInfo> {
    let mut databases = Vec::new();

    for entry in WalkDir::new(dir)
        .max_depth(max_depth)
        .into_iter()
        .filter_map(|e| e.ok())
    {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }

        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Get parent directory for context
        let context = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();

        // Check for binary format first
        if let Some(num_trees) = detect_binary_database(path) {
            let size = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
            let compressed = name.ends_with(".gz");
            databases.push(DatabaseInfo {
                name,
                context,
                format: "binary".to_string(),
                compressed,
                num_trees: num_trees.to_string(),
                size: format_size(size, BINARY),
                size_bytes: size,
                path: path.to_path_buf(),
            });
            continue;
        }

        // Check for text format
        let is_text_db = name.contains("reference_tree")
            || name.ends_with("_tree.txt")
            || name.ends_with("_tree.txt.gz");

        if is_text_db {
            if let Some(num_trees) = detect_text_database(path) {
                let size = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
                let compressed = name.ends_with(".gz");
                databases.push(DatabaseInfo {
                    name,
                    context,
                    format: "text".to_string(),
                    compressed,
                    num_trees: num_trees.to_string(),
                    size: format_size(size, BINARY),
                    size_bytes: size,
                    path: path.to_path_buf(),
                });
            }
        }
    }

    databases.sort_by(|a, b| a.path.cmp(&b.path));
    databases
}
```

#### 2. Update list rendering to show enhanced info
**File**: `rust/crates/tronko-bench/src/tui/ui/list.rs`
**Changes**: Update `render_database_items()` to show context and compression status

```rust
fn render_database_items(app: &App) -> Vec<ListItem<'static>> {
    let indices: Vec<usize> = if app.filtered_indices.is_empty() {
        (0..app.scan_results.databases.len()).collect()
    } else {
        app.filtered_indices.clone()
    };

    indices
        .into_iter()
        .filter_map(|i| app.scan_results.databases.get(i))
        .map(|db| {
            let format_color = if db.format == "binary" {
                Color::Green
            } else {
                Color::Yellow
            };

            let compression_indicator = if db.compressed { "gz" } else { "" };

            // Show context/name together for clarity (e.g., "16S_Bacteria/reference_tree.trkb")
            let display_name = if db.context.is_empty() {
                db.name.clone()
            } else {
                format!("{}/{}", db.context, db.name)
            };

            ListItem::new(Line::from(vec![
                Span::styled(
                    format!("{:<40}", truncate(&display_name, 40)),
                    Style::default(),
                ),
                Span::styled(
                    format!(" {:>4} ", db.num_trees),
                    Style::default().fg(Color::Cyan),
                ),
                Span::styled(
                    format!("{:>10} ", db.size),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::styled(
                    format!("[{}{}]", if db.format == "binary" { "BIN" } else { "TXT" },
                            if db.compressed { ".gz" } else { "" }),
                    Style::default().fg(format_color),
                ),
            ]))
        })
        .collect()
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Database list shows parent directory context (e.g., "16S_Bacteria/reference_tree.trkb")
- [ ] Compression status shown in format indicator ([BIN.gz] vs [BIN])

---

## Phase 4: App State Extension

### Overview
Extend the App state to support benchmark workflow, history, and comparison mode.

### Changes Required:

#### 1. Create new app mode enum and state
**File**: `rust/crates/tronko-bench/src/tui/app.rs`
**Changes**: Add new modes, tabs, and state fields

```rust
use std::path::PathBuf;
use crate::ScanResults;
use crate::benchmark::{BenchmarkStore, BenchmarkResult, DiscoveredProject, QueryMode};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Databases,
    Queries,
    MsaFiles,
    Benchmarks,  // NEW: Historical benchmarks
}

impl Tab {
    pub fn all() -> &'static [Tab] {
        &[Tab::Databases, Tab::Queries, Tab::MsaFiles, Tab::Benchmarks]
    }

    pub fn title(&self) -> &'static str {
        match self {
            Tab::Databases => "Databases",
            Tab::Queries => "Queries",
            Tab::MsaFiles => "MSA Files",
            Tab::Benchmarks => "Benchmarks",
        }
    }
}

/// Mode the app is currently in
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AppMode {
    Normal,
    Search,
    Help,
    BenchmarkWizard,  // NEW: Multi-step benchmark creation
    Comparison,       // NEW: Comparing two benchmarks
}

/// Step in the benchmark wizard
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WizardStep {
    SelectDatabase,
    SelectProject,
    SelectMarker,
    SelectQueryMode,
    SelectBinary,
    Confirm,
}

impl WizardStep {
    pub fn title(&self) -> &'static str {
        match self {
            WizardStep::SelectDatabase => "Select Reference Database",
            WizardStep::SelectProject => "Select Project",
            WizardStep::SelectMarker => "Select Marker",
            WizardStep::SelectQueryMode => "Select Query Type",
            WizardStep::SelectBinary => "Select tronko-assign Binary",
            WizardStep::Confirm => "Confirm & Run",
        }
    }

    pub fn next(&self) -> Option<WizardStep> {
        match self {
            WizardStep::SelectDatabase => Some(WizardStep::SelectProject),
            WizardStep::SelectProject => Some(WizardStep::SelectMarker),
            WizardStep::SelectMarker => Some(WizardStep::SelectQueryMode),
            WizardStep::SelectQueryMode => Some(WizardStep::SelectBinary),
            WizardStep::SelectBinary => Some(WizardStep::Confirm),
            WizardStep::Confirm => None,
        }
    }

    pub fn prev(&self) -> Option<WizardStep> {
        match self {
            WizardStep::SelectDatabase => None,
            WizardStep::SelectProject => Some(WizardStep::SelectDatabase),
            WizardStep::SelectMarker => Some(WizardStep::SelectProject),
            WizardStep::SelectQueryMode => Some(WizardStep::SelectMarker),
            WizardStep::SelectBinary => Some(WizardStep::SelectQueryMode),
            WizardStep::Confirm => Some(WizardStep::SelectBinary),
        }
    }
}

/// State for the benchmark wizard
#[derive(Debug, Clone)]
pub struct WizardState {
    pub step: WizardStep,
    pub selected_db_index: Option<usize>,
    pub discovered_projects: Vec<DiscoveredProject>,
    pub selected_project_index: Option<usize>,
    pub selected_marker_index: Option<usize>,
    pub selected_query_mode: Option<QueryMode>,
    pub tronko_binaries: Vec<TronkoBinary>,
    pub selected_binary_index: Option<usize>,
    pub list_index: usize,  // Current selection in step's list
}

/// Discovered tronko-assign binary
#[derive(Debug, Clone)]
pub struct TronkoBinary {
    pub path: PathBuf,
    pub optimize_memory: bool,
    pub label: String,
}

impl Default for WizardState {
    fn default() -> Self {
        Self {
            step: WizardStep::SelectDatabase,
            selected_db_index: None,
            discovered_projects: Vec::new(),
            selected_project_index: None,
            selected_marker_index: None,
            selected_query_mode: None,
            tronko_binaries: Vec::new(),
            selected_binary_index: None,
            list_index: 0,
        }
    }
}

pub struct App {
    // Data
    pub scan_results: ScanResults,
    pub scan_dir: PathBuf,
    pub max_depth: usize,
    pub benchmark_store: BenchmarkStore,
    pub benchmarks_path: PathBuf,

    // Navigation
    pub current_tab: Tab,
    pub selected_index: usize,
    pub mode: AppMode,  // NEW: replaces search_mode, show_help

    // Search (kept for backwards compat)
    pub search_mode: bool,
    pub search_query: String,
    pub filtered_indices: Vec<usize>,

    // Status
    pub loading: bool,
    pub status_message: Option<String>,
    pub should_quit: bool,

    // Help (kept for backwards compat)
    pub show_help: bool,

    // NEW: Wizard state
    pub wizard: WizardState,

    // NEW: Comparison state
    pub selected_benchmarks: Vec<usize>,  // indices of selected benchmarks (max 2)
}

impl App {
    pub fn new(scan_dir: PathBuf, max_depth: usize) -> Self {
        let benchmarks_path = scan_dir.join(".tronko-bench-results.json");
        let benchmark_store = BenchmarkStore::load(&benchmarks_path).unwrap_or_default();

        Self {
            scan_results: ScanResults {
                databases: Vec::new(),
                query_files: Vec::new(),
                msa_files: Vec::new(),
            },
            scan_dir,
            max_depth,
            benchmark_store,
            benchmarks_path,
            current_tab: Tab::Databases,
            selected_index: 0,
            mode: AppMode::Normal,
            search_mode: false,
            search_query: String::new(),
            filtered_indices: Vec::new(),
            loading: true,
            status_message: Some("Scanning...".to_string()),
            should_quit: false,
            show_help: false,
            wizard: WizardState::default(),
            selected_benchmarks: Vec::new(),
        }
    }

    // ... keep existing methods, add new ones below

    pub fn current_list_len(&self) -> usize {
        if !self.filtered_indices.is_empty() {
            return self.filtered_indices.len();
        }
        match self.current_tab {
            Tab::Databases => self.scan_results.databases.len(),
            Tab::Queries => self.scan_results.query_files.len(),
            Tab::MsaFiles => self.scan_results.msa_files.len(),
            Tab::Benchmarks => self.benchmark_store.benchmarks.len(),
        }
    }

    pub fn toggle_benchmark_selection(&mut self) {
        if self.current_tab != Tab::Benchmarks {
            return;
        }

        let idx = self.selected_index;
        if self.selected_benchmarks.contains(&idx) {
            self.selected_benchmarks.retain(|&i| i != idx);
        } else if self.selected_benchmarks.len() < 2 {
            self.selected_benchmarks.push(idx);
        }

        // Enter comparison mode if 2 selected
        if self.selected_benchmarks.len() == 2 {
            self.mode = AppMode::Comparison;
        }
    }

    pub fn exit_comparison(&mut self) {
        self.mode = AppMode::Normal;
        self.selected_benchmarks.clear();
    }

    pub fn start_benchmark_wizard(&mut self) {
        self.wizard = WizardState::default();
        self.mode = AppMode::BenchmarkWizard;
    }

    pub fn cancel_wizard(&mut self) {
        self.mode = AppMode::Normal;
    }

    pub fn save_benchmarks(&self) -> anyhow::Result<()> {
        self.benchmark_store.save(&self.benchmarks_path)
    }

    // Keep existing methods...
    pub fn select_next(&mut self) {
        let max = self.current_list_len();
        if max > 0 {
            self.selected_index = (self.selected_index + 1).min(max - 1);
        }
    }

    pub fn select_previous(&mut self) {
        self.selected_index = self.selected_index.saturating_sub(1);
    }

    pub fn next_tab(&mut self) {
        let tabs = Tab::all();
        let current = tabs.iter().position(|t| *t == self.current_tab).unwrap_or(0);
        self.current_tab = tabs[(current + 1) % tabs.len()];
        self.selected_index = 0;
        self.clear_search();
    }

    pub fn prev_tab(&mut self) {
        let tabs = Tab::all();
        let current = tabs.iter().position(|t| *t == self.current_tab).unwrap_or(0);
        self.current_tab = tabs[(current + tabs.len() - 1) % tabs.len()];
        self.selected_index = 0;
        self.clear_search();
    }

    pub fn clear_search(&mut self) {
        self.search_mode = false;
        self.search_query.clear();
        self.filtered_indices.clear();
    }
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] N/A (state changes only)

---

## Phase 5: Benchmark Tab UI

### Overview
Implement the UI for displaying historical benchmarks in a new tab.

### Changes Required:

#### 1. Update tabs rendering to include Benchmarks
**File**: `rust/crates/tronko-bench/src/tui/ui/tabs.rs`
**Changes**: Already handled by Tab::all() change in Phase 4

#### 2. Update list rendering to handle Benchmarks tab
**File**: `rust/crates/tronko-bench/src/tui/ui/list.rs`
**Changes**: Add case for Benchmarks tab in render_list

```rust
pub fn render_list(f: &mut Frame, area: Rect, app: &App) {
    let (items, empty_msg) = match app.current_tab {
        Tab::Databases => (render_database_items(app), "No databases found."),
        Tab::Queries => (render_query_items(app, false), "No query files found."),
        Tab::MsaFiles => (render_query_items(app, true), "No MSA files found."),
        Tab::Benchmarks => (render_benchmark_items(app), "No benchmarks yet. Press 'b' to create one."),
    };
    // ... rest unchanged
}

fn render_benchmark_items(app: &App) -> Vec<ListItem<'static>> {
    app.benchmark_store.benchmarks
        .iter()
        .enumerate()
        .map(|(idx, bench)| {
            let selected = app.selected_benchmarks.contains(&idx);
            let select_indicator = if selected { "[*]" } else { "[ ]" };

            let status_color = if bench.success {
                Color::Green
            } else {
                Color::Red
            };

            let duration = format!("{:.1}s", bench.duration_secs);
            let memory = bench.peak_memory_bytes
                .map(|b| humansize::format_size(b, humansize::BINARY))
                .unwrap_or_else(|| "N/A".to_string());

            ListItem::new(Line::from(vec![
                Span::styled(
                    format!("{} ", select_indicator),
                    Style::default().fg(if selected { Color::Cyan } else { Color::DarkGray }),
                ),
                Span::styled(
                    format!("{:<25}", truncate(&bench.config.name, 25)),
                    Style::default(),
                ),
                Span::styled(
                    format!(" {:>8} ", duration),
                    Style::default().fg(Color::Yellow),
                ),
                Span::styled(
                    format!("{:>10} ", memory),
                    Style::default().fg(Color::Magenta),
                ),
                Span::styled(
                    if bench.success { "[OK]" } else { "[FAIL]" },
                    Style::default().fg(status_color),
                ),
            ]))
        })
        .collect()
}
```

#### 3. Update detail rendering for Benchmarks tab
**File**: `rust/crates/tronko-bench/src/tui/ui/detail.rs`
**Changes**: Add benchmark detail rendering

```rust
pub fn render_detail(f: &mut Frame, area: Rect, app: &App) {
    let content = match app.current_tab {
        Tab::Databases => render_database_detail(app),
        Tab::Queries | Tab::MsaFiles => render_query_detail(app),
        Tab::Benchmarks => render_benchmark_detail(app),
    };
    // ... rest unchanged
}

fn render_benchmark_detail(app: &App) -> Vec<Line<'static>> {
    let bench = match app.benchmark_store.benchmarks.get(app.selected_index) {
        Some(b) => b,
        None => return vec![Line::from("No benchmark selected")],
    };

    let mut lines = vec![
        Line::from(vec![
            Span::styled("Name: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(bench.config.name.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Status: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::styled(
                if bench.success { "Success" } else { "Failed" },
                Style::default().fg(if bench.success { Color::Green } else { Color::Red }),
            ),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Duration: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(format!("{:.2} seconds", bench.duration_secs)),
        ]),
    ];

    if let Some(mem) = bench.peak_memory_bytes {
        lines.push(Line::from(""));
        lines.push(Line::from(vec![
            Span::styled("Peak Memory: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(humansize::format_size(mem, humansize::BINARY)),
        ]));
    }

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled("Database: ", Style::default().add_modifier(Modifier::BOLD)),
        Span::raw(format!("{} ({})",
            bench.config.reference_db_path.file_name().unwrap_or_default().to_string_lossy(),
            bench.config.db_format)),
    ]));

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled("Query Type: ", Style::default().add_modifier(Modifier::BOLD)),
        Span::raw(if bench.config.paired_end { "Paired-end" } else { "Single-end" }),
    ]));

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled("Marker: ", Style::default().add_modifier(Modifier::BOLD)),
        Span::raw(bench.config.marker.clone()),
    ]));

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled("Timestamp: ", Style::default().add_modifier(Modifier::BOLD)),
        Span::raw(bench.timestamp.clone()),
    ]));

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled(
            "Space to select, compare 2 benchmarks",
            Style::default().fg(Color::DarkGray),
        ),
    ]));

    lines
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Benchmarks tab appears in TUI
- [ ] Benchmark list shows name, duration, memory, status
- [ ] Detail panel shows full benchmark info when selected

---

## Phase 6: Benchmark Wizard UI

### Overview
Implement the multi-step wizard overlay for creating new benchmarks.

### Changes Required:

#### 1. Create wizard UI module
**File**: `rust/crates/tronko-bench/src/tui/ui/wizard.rs` (new file)
**Changes**: Implement wizard overlay rendering

```rust
use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph},
    Frame,
};
use crate::tui::app::{App, WizardStep};

pub fn render_wizard(f: &mut Frame, app: &App) {
    let area = centered_rect(80, 80, f.area());
    f.render_widget(Clear, area);

    // Split into header and content
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),  // Title
            Constraint::Length(2),  // Progress
            Constraint::Min(10),    // Content
            Constraint::Length(2),  // Instructions
        ])
        .split(area);

    // Title
    let title = Paragraph::new(app.wizard.step.title())
        .style(Style::default().add_modifier(Modifier::BOLD))
        .block(Block::default().borders(Borders::ALL).title(" New Benchmark "));
    f.render_widget(title, chunks[0]);

    // Progress indicator
    let progress = render_progress(&app.wizard.step);
    f.render_widget(progress, chunks[1]);

    // Content based on step
    match app.wizard.step {
        WizardStep::SelectDatabase => render_database_selection(f, chunks[2], app),
        WizardStep::SelectProject => render_project_selection(f, chunks[2], app),
        WizardStep::SelectMarker => render_marker_selection(f, chunks[2], app),
        WizardStep::SelectQueryMode => render_query_mode_selection(f, chunks[2], app),
        WizardStep::SelectBinary => render_binary_selection(f, chunks[2], app),
        WizardStep::Confirm => render_confirmation(f, chunks[2], app),
    }

    // Instructions
    let instructions = Paragraph::new("↑↓: Navigate  Enter: Select  Esc: Cancel  ←: Back")
        .style(Style::default().fg(Color::DarkGray));
    f.render_widget(instructions, chunks[3]);
}

fn render_progress(step: &WizardStep) -> Paragraph<'static> {
    let steps = [
        ("1.DB", WizardStep::SelectDatabase),
        ("2.Project", WizardStep::SelectProject),
        ("3.Marker", WizardStep::SelectMarker),
        ("4.Query", WizardStep::SelectQueryMode),
        ("5.Binary", WizardStep::SelectBinary),
        ("6.Run", WizardStep::Confirm),
    ];

    let spans: Vec<Span> = steps
        .iter()
        .map(|(label, s)| {
            let style = if s == step {
                Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD)
            } else {
                Style::default().fg(Color::DarkGray)
            };
            Span::styled(format!(" {} ", label), style)
        })
        .collect();

    Paragraph::new(Line::from(spans))
}

fn render_database_selection(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = app.scan_results.databases
        .iter()
        .map(|db| {
            let display = format!(
                "{}/{} - {} trees, {} [{}{}]",
                db.context,
                db.name,
                db.num_trees,
                db.size,
                db.format,
                if db.compressed { ".gz" } else { "" }
            );
            ListItem::new(display)
        })
        .collect();

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray).add_modifier(Modifier::BOLD))
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.list_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_project_selection(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = app.wizard.discovered_projects
        .iter()
        .map(|proj| {
            let markers: Vec<_> = proj.markers.iter().map(|m| m.name.as_str()).collect();
            ListItem::new(format!("{} (markers: {})", proj.name, markers.join(", ")))
        })
        .collect();

    if items.is_empty() {
        let msg = Paragraph::new("No projects found in scan directory")
            .style(Style::default().fg(Color::Yellow))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(msg, area);
        return;
    }

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray).add_modifier(Modifier::BOLD))
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.list_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_marker_selection(f: &mut Frame, area: Rect, app: &App) {
    let project = match app.wizard.selected_project_index
        .and_then(|i| app.wizard.discovered_projects.get(i))
    {
        Some(p) => p,
        None => {
            f.render_widget(Paragraph::new("No project selected"), area);
            return;
        }
    };

    let items: Vec<ListItem> = project.markers
        .iter()
        .map(|marker| {
            let mut modes = Vec::new();
            if marker.has_paired { modes.push("paired"); }
            if marker.has_unpaired_f { modes.push("unpaired_F"); }
            if marker.has_unpaired_r { modes.push("unpaired_R"); }
            ListItem::new(format!("{} ({})", marker.name, modes.join(", ")))
        })
        .collect();

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray).add_modifier(Modifier::BOLD))
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.list_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_query_mode_selection(f: &mut Frame, area: Rect, app: &App) {
    let marker = match app.wizard.selected_project_index
        .and_then(|pi| app.wizard.discovered_projects.get(pi))
        .and_then(|p| app.wizard.selected_marker_index.and_then(|mi| p.markers.get(mi)))
    {
        Some(m) => m,
        None => {
            f.render_widget(Paragraph::new("No marker selected"), area);
            return;
        }
    };

    let mut items = Vec::new();
    if marker.has_paired {
        items.push(ListItem::new("Paired-end reads"));
    }
    if marker.has_unpaired_f {
        items.push(ListItem::new("Unpaired (Forward) reads"));
    }
    if marker.has_unpaired_r {
        items.push(ListItem::new("Unpaired (Reverse) reads"));
    }

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray).add_modifier(Modifier::BOLD))
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.list_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_binary_selection(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = app.wizard.tronko_binaries
        .iter()
        .map(|bin| ListItem::new(bin.label.clone()))
        .collect();

    if items.is_empty() {
        let msg = Paragraph::new("No tronko-assign binaries found.\nBuild with: cd tronko-assign && make\nFor memory optimized: make OPTIMIZE_MEMORY=1")
            .style(Style::default().fg(Color::Yellow))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(msg, area);
        return;
    }

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL))
        .highlight_style(Style::default().bg(Color::DarkGray).add_modifier(Modifier::BOLD))
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.list_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_confirmation(f: &mut Frame, area: Rect, app: &App) {
    // Build summary of selected options
    let mut lines = vec![
        Line::from(Span::styled("Benchmark Configuration:", Style::default().add_modifier(Modifier::BOLD))),
        Line::from(""),
    ];

    if let Some(idx) = app.wizard.selected_db_index {
        if let Some(db) = app.scan_results.databases.get(idx) {
            lines.push(Line::from(format!("Database: {}/{}", db.context, db.name)));
            lines.push(Line::from(format!("  Format: {}{}", db.format, if db.compressed { " (gzipped)" } else { "" })));
            lines.push(Line::from(format!("  Trees: {}, Size: {}", db.num_trees, db.size)));
        }
    }

    lines.push(Line::from(""));

    if let Some(proj_idx) = app.wizard.selected_project_index {
        if let Some(proj) = app.wizard.discovered_projects.get(proj_idx) {
            lines.push(Line::from(format!("Project: {}", proj.name)));
        }
    }

    if let Some(query_mode) = &app.wizard.selected_query_mode {
        lines.push(Line::from(format!("Query Mode: {}", query_mode.label())));
    }

    lines.push(Line::from(""));

    if let Some(idx) = app.wizard.selected_binary_index {
        if let Some(bin) = app.wizard.tronko_binaries.get(idx) {
            lines.push(Line::from(format!("Binary: {}", bin.label)));
        }
    }

    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled(
        "Press Enter to run benchmark, Esc to cancel",
        Style::default().fg(Color::Cyan),
    )));

    let paragraph = Paragraph::new(lines)
        .block(Block::default().borders(Borders::ALL));
    f.render_widget(paragraph, area);
}

fn centered_rect(percent_x: u16, percent_y: u16, r: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}
```

#### 2. Update UI mod.rs to include wizard
**File**: `rust/crates/tronko-bench/src/tui/ui/mod.rs`
**Changes**: Add wizard module and render call

```rust
pub mod tabs;
pub mod list;
pub mod detail;
pub mod status;
pub mod help;
pub mod wizard;      // NEW
pub mod comparison;  // NEW (for Phase 7)

use ratatui::{/* ... existing imports ... */};
use crate::tui::app::{App, AppMode};

pub fn render(f: &mut Frame, app: &App) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(2),
            Constraint::Min(3),
            Constraint::Length(1),
        ])
        .split(f.area());

    tabs::render_tabs(f, chunks[0], app);

    // Split main area for list + detail
    let main_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(60),
            Constraint::Percentage(40),
        ])
        .split(chunks[1]);

    list::render_list(f, main_chunks[0], app);
    detail::render_detail(f, main_chunks[1], app);

    status::render_status(f, chunks[2], app);

    // Render overlays based on mode
    match app.mode {
        AppMode::Help => help::render_help(f),
        AppMode::Search => render_search_overlay(f, &app.search_query),
        AppMode::BenchmarkWizard => wizard::render_wizard(f, app),
        AppMode::Comparison => comparison::render_comparison(f, app),
        AppMode::Normal => {}
    }
}

// ... rest of file unchanged
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Pressing 'b' opens the benchmark wizard
- [ ] Can navigate through wizard steps with arrow keys and Enter
- [ ] Wizard shows appropriate options at each step
- [ ] Can cancel wizard with Esc

---

## Phase 7: Comparison View UI

### Overview
Implement the side-by-side comparison view for two benchmarks.

### Changes Required:

#### 1. Create comparison UI module
**File**: `rust/crates/tronko-bench/src/tui/ui/comparison.rs` (new file)
**Changes**: Implement comparison overlay

```rust
use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph, Row, Table},
    Frame,
};
use crate::tui::app::App;
use crate::benchmark::BenchmarkResult;

pub fn render_comparison(f: &mut Frame, app: &App) {
    let area = centered_rect(90, 85, f.area());
    f.render_widget(Clear, area);

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),   // Title
            Constraint::Min(10),     // Main comparison
            Constraint::Length(2),   // Instructions
        ])
        .split(area);

    // Title
    let title = Paragraph::new("Benchmark Comparison")
        .style(Style::default().add_modifier(Modifier::BOLD))
        .block(Block::default().borders(Borders::ALL).title(" Compare "));
    f.render_widget(title, chunks[0]);

    // Get the two benchmarks
    let (bench1, bench2) = match (
        app.selected_benchmarks.get(0).and_then(|i| app.benchmark_store.benchmarks.get(*i)),
        app.selected_benchmarks.get(1).and_then(|i| app.benchmark_store.benchmarks.get(*i)),
    ) {
        (Some(b1), Some(b2)) => (b1, b2),
        _ => {
            let msg = Paragraph::new("Select exactly 2 benchmarks to compare")
                .block(Block::default().borders(Borders::ALL));
            f.render_widget(msg, chunks[1]);
            return;
        }
    };

    render_comparison_table(f, chunks[1], bench1, bench2);

    // Instructions
    let instructions = Paragraph::new("Esc: Close comparison")
        .style(Style::default().fg(Color::DarkGray));
    f.render_widget(instructions, chunks[2]);
}

fn render_comparison_table(f: &mut Frame, area: Rect, b1: &BenchmarkResult, b2: &BenchmarkResult) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(12),  // Basic stats
            Constraint::Min(5),      // Summary/winner
        ])
        .split(area);

    // Create comparison rows
    let rows = vec![
        create_row("Name", &b1.config.name, &b2.config.name, None),
        create_row("Status",
            if b1.success { "Success" } else { "Failed" },
            if b2.success { "Success" } else { "Failed" },
            None),
        create_duration_row("Duration", b1.duration_secs, b2.duration_secs),
        create_memory_row("Peak Memory", b1.peak_memory_bytes, b2.peak_memory_bytes),
        create_row("Database Format", &b1.config.db_format, &b2.config.db_format, None),
        create_row("Compressed",
            if b1.config.reference_db_path.to_string_lossy().contains(".gz") { "Yes" } else { "No" },
            if b2.config.reference_db_path.to_string_lossy().contains(".gz") { "Yes" } else { "No" },
            None),
        create_row("Query Mode",
            if b1.config.paired_end { "Paired" } else { "Single" },
            if b2.config.paired_end { "Paired" } else { "Single" },
            None),
        create_row("OPTIMIZE_MEMORY",
            if b1.config.optimize_memory { "Yes" } else { "No" },
            if b2.config.optimize_memory { "Yes" } else { "No" },
            None),
    ];

    let table = Table::new(
        rows,
        [
            Constraint::Length(20),  // Label
            Constraint::Percentage(35), // Benchmark 1
            Constraint::Percentage(35), // Benchmark 2
            Constraint::Length(10),  // Winner indicator
        ],
    )
    .header(Row::new(vec!["Metric", "Benchmark 1", "Benchmark 2", "Winner"])
        .style(Style::default().add_modifier(Modifier::BOLD)))
    .block(Block::default().borders(Borders::ALL).title(" Comparison "));

    f.render_widget(table, chunks[0]);

    // Summary
    let summary = create_summary(b1, b2);
    let summary_widget = Paragraph::new(summary)
        .block(Block::default().borders(Borders::ALL).title(" Summary "));
    f.render_widget(summary_widget, chunks[1]);
}

fn create_row<'a>(label: &'a str, v1: &str, v2: &str, winner: Option<u8>) -> Row<'a> {
    let winner_str = match winner {
        Some(1) => "← #1",
        Some(2) => "#2 →",
        _ => "",
    };
    Row::new(vec![
        label.to_string(),
        v1.to_string(),
        v2.to_string(),
        winner_str.to_string(),
    ])
}

fn create_duration_row<'a>(label: &'a str, d1: f64, d2: f64) -> Row<'a> {
    let winner = if d1 < d2 { Some(1) } else if d2 < d1 { Some(2) } else { None };
    let diff = ((d1 - d2).abs() / d1.max(d2) * 100.0) as i32;

    let v1 = format!("{:.2}s", d1);
    let v2 = format!("{:.2}s ({}{}%)", d2, if d2 < d1 { "-" } else { "+" }, diff);

    let winner_str = match winner {
        Some(1) => format!("← {}% faster", diff),
        Some(2) => format!("{}% faster →", diff),
        _ => "Tie".to_string(),
    };

    Row::new(vec![label.to_string(), v1, v2, winner_str])
        .style(Style::default().fg(Color::Yellow))
}

fn create_memory_row<'a>(label: &'a str, m1: Option<u64>, m2: Option<u64>) -> Row<'a> {
    match (m1, m2) {
        (Some(mem1), Some(mem2)) => {
            let winner = if mem1 < mem2 { Some(1) } else if mem2 < mem1 { Some(2) } else { None };
            let diff = ((mem1 as f64 - mem2 as f64).abs() / mem1.max(mem2) as f64 * 100.0) as i32;

            let v1 = humansize::format_size(mem1, humansize::BINARY);
            let v2 = format!("{} ({}{}%)",
                humansize::format_size(mem2, humansize::BINARY),
                if mem2 < mem1 { "-" } else { "+" },
                diff);

            let winner_str = match winner {
                Some(1) => format!("← {}% less", diff),
                Some(2) => format!("{}% less →", diff),
                _ => "Tie".to_string(),
            };

            Row::new(vec![label.to_string(), v1, v2, winner_str])
                .style(Style::default().fg(Color::Magenta))
        }
        _ => Row::new(vec![label.to_string(), "N/A".to_string(), "N/A".to_string(), "".to_string()]),
    }
}

fn create_summary(b1: &BenchmarkResult, b2: &BenchmarkResult) -> Vec<Line<'static>> {
    let mut lines = Vec::new();

    // Speed comparison
    if b1.duration_secs < b2.duration_secs {
        let speedup = ((b2.duration_secs / b1.duration_secs - 1.0) * 100.0) as i32;
        lines.push(Line::from(Span::styled(
            format!("Speed: Benchmark 1 is {}% faster", speedup),
            Style::default().fg(Color::Green),
        )));
    } else if b2.duration_secs < b1.duration_secs {
        let speedup = ((b1.duration_secs / b2.duration_secs - 1.0) * 100.0) as i32;
        lines.push(Line::from(Span::styled(
            format!("Speed: Benchmark 2 is {}% faster", speedup),
            Style::default().fg(Color::Green),
        )));
    } else {
        lines.push(Line::from("Speed: Equal"));
    }

    // Memory comparison
    match (b1.peak_memory_bytes, b2.peak_memory_bytes) {
        (Some(m1), Some(m2)) if m1 < m2 => {
            let savings = (((m2 - m1) as f64 / m2 as f64) * 100.0) as i32;
            lines.push(Line::from(Span::styled(
                format!("Memory: Benchmark 1 uses {}% less memory", savings),
                Style::default().fg(Color::Green),
            )));
        }
        (Some(m1), Some(m2)) if m2 < m1 => {
            let savings = (((m1 - m2) as f64 / m1 as f64) * 100.0) as i32;
            lines.push(Line::from(Span::styled(
                format!("Memory: Benchmark 2 uses {}% less memory", savings),
                Style::default().fg(Color::Green),
            )));
        }
        _ => lines.push(Line::from("Memory: No comparison available")),
    }

    // Results comparison (if output summaries exist)
    match (&b1.output_summary, &b2.output_summary) {
        (Some(s1), Some(s2)) => {
            if s1.assignments_hash == s2.assignments_hash {
                lines.push(Line::from(Span::styled(
                    "Results: Identical assignments!",
                    Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
                )));
            } else {
                lines.push(Line::from(Span::styled(
                    "Results: Different assignments detected",
                    Style::default().fg(Color::Yellow),
                )));
                lines.push(Line::from(format!(
                    "  B1: {} assigned/{} total, {} unique taxa",
                    s1.assigned_lines, s1.total_lines, s1.unique_taxa
                )));
                lines.push(Line::from(format!(
                    "  B2: {} assigned/{} total, {} unique taxa",
                    s2.assigned_lines, s2.total_lines, s2.unique_taxa
                )));
            }
        }
        _ => lines.push(Line::from("Results: No output summary available")),
    }

    lines
}

fn centered_rect(percent_x: u16, percent_y: u16, r: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(popup_layout[1])[1]
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Selecting 2 benchmarks with spacebar opens comparison view
- [ ] Comparison shows side-by-side metrics
- [ ] Winner indicators show which benchmark is better for each metric
- [ ] Summary section answers: which is faster, which uses less memory, are results identical

---

## Phase 8: Event Handling Updates

### Overview
Update event handling to support new messages for wizard and comparison modes.

### Changes Required:

#### 1. Update Message enum
**File**: `rust/crates/tronko-bench/src/tui/message.rs`
**Changes**: Add new message variants

```rust
#[derive(Debug)]
pub enum Message {
    // Navigation
    NextTab,
    PrevTab,
    SelectNext,
    SelectPrev,

    // Search
    EnterSearch,
    ExitSearch,
    SearchInput(char),
    SearchBackspace,

    // Actions
    Rescan,
    GenerateConfig,
    ConvertSelected,

    // Help
    ShowHelp,
    HideHelp,

    // NEW: Benchmark wizard
    StartBenchmarkWizard,
    WizardNext,
    WizardPrev,
    WizardCancel,
    WizardConfirm,
    WizardSelectItem,

    // NEW: Comparison
    ToggleBenchmarkSelection,
    ExitComparison,

    // App
    Quit,
}
```

#### 2. Update event handler
**File**: `rust/crates/tronko-bench/src/tui/event.rs`
**Changes**: Handle new key events based on mode

```rust
use crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use std::time::Duration;
use crate::tui::app::{App, AppMode, Tab};
use crate::tui::message::Message;

pub fn poll_event(timeout: Duration) -> std::io::Result<Option<Event>> {
    if event::poll(timeout)? {
        Ok(Some(event::read()?))
    } else {
        Ok(None)
    }
}

pub fn handle_key_event(app: &App, key: event::KeyEvent) -> Option<Message> {
    if key.kind != KeyEventKind::Press {
        return None;
    }

    // Handle based on current mode
    match app.mode {
        AppMode::Help => Some(Message::HideHelp),
        AppMode::Search => handle_search_mode(key),
        AppMode::BenchmarkWizard => handle_wizard_mode(key),
        AppMode::Comparison => handle_comparison_mode(key),
        AppMode::Normal => handle_normal_mode(app, key),
    }
}

fn handle_search_mode(key: event::KeyEvent) -> Option<Message> {
    match key.code {
        KeyCode::Esc => Some(Message::ExitSearch),
        KeyCode::Enter => Some(Message::ExitSearch),
        KeyCode::Backspace => Some(Message::SearchBackspace),
        KeyCode::Char(c) => Some(Message::SearchInput(c)),
        _ => None,
    }
}

fn handle_wizard_mode(key: event::KeyEvent) -> Option<Message> {
    match key.code {
        KeyCode::Esc => Some(Message::WizardCancel),
        KeyCode::Enter => Some(Message::WizardSelectItem),
        KeyCode::Up | KeyCode::Char('k') => Some(Message::SelectPrev),
        KeyCode::Down | KeyCode::Char('j') => Some(Message::SelectNext),
        KeyCode::Left | KeyCode::Char('h') => Some(Message::WizardPrev),
        KeyCode::Right | KeyCode::Char('l') => Some(Message::WizardNext),
        _ => None,
    }
}

fn handle_comparison_mode(key: event::KeyEvent) -> Option<Message> {
    match key.code {
        KeyCode::Esc | KeyCode::Char('q') => Some(Message::ExitComparison),
        _ => None,
    }
}

fn handle_normal_mode(app: &App, key: event::KeyEvent) -> Option<Message> {
    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => Some(Message::Quit),
        KeyCode::Tab => Some(Message::NextTab),
        KeyCode::BackTab => Some(Message::PrevTab),
        KeyCode::Up | KeyCode::Char('k') => Some(Message::SelectPrev),
        KeyCode::Down | KeyCode::Char('j') => Some(Message::SelectNext),
        KeyCode::Char('/') => Some(Message::EnterSearch),
        KeyCode::Char('r') => Some(Message::Rescan),
        KeyCode::Char('g') => Some(Message::GenerateConfig),
        KeyCode::Char('c') => Some(Message::ConvertSelected),
        KeyCode::Char('?') => Some(Message::ShowHelp),
        KeyCode::Char('b') => Some(Message::StartBenchmarkWizard),
        KeyCode::Char(' ') if app.current_tab == Tab::Benchmarks => {
            Some(Message::ToggleBenchmarkSelection)
        }
        KeyCode::Char('l') | KeyCode::Right => {
            if key.modifiers.contains(KeyModifiers::SHIFT) {
                Some(Message::NextTab)
            } else {
                None
            }
        }
        KeyCode::Char('h') | KeyCode::Left => {
            if key.modifiers.contains(KeyModifiers::SHIFT) {
                Some(Message::PrevTab)
            } else {
                None
            }
        }
        _ => None,
    }
}
```

#### 3. Update message handler in mod.rs
**File**: `rust/crates/tronko-bench/src/tui/mod.rs`
**Changes**: Handle new messages in `handle_message()`

```rust
fn handle_message(app: &mut App, msg: Message) -> Result<()> {
    match msg {
        Message::Quit => app.should_quit = true,
        Message::NextTab => app.next_tab(),
        Message::PrevTab => app.prev_tab(),
        Message::SelectNext => {
            if app.mode == AppMode::BenchmarkWizard {
                wizard_select_next(app);
            } else {
                app.select_next();
            }
        }
        Message::SelectPrev => {
            if app.mode == AppMode::BenchmarkWizard {
                wizard_select_prev(app);
            } else {
                app.select_previous();
            }
        }
        Message::EnterSearch => {
            app.mode = AppMode::Search;
            app.search_mode = true;
            app.search_query.clear();
        }
        Message::ExitSearch => {
            app.mode = AppMode::Normal;
            app.search_mode = false;
        }
        Message::SearchInput(c) => {
            app.search_query.push(c);
            update_search_filter(app);
        }
        Message::SearchBackspace => {
            app.search_query.pop();
            update_search_filter(app);
        }
        Message::Rescan => {
            app.loading = true;
            app.status_message = Some("Rescanning...".to_string());
            let results = perform_scan(&app.scan_dir, app.max_depth);
            app.scan_results = results;
            app.loading = false;
            app.status_message = Some("Scan complete".to_string());
            app.clear_search();
            app.selected_index = 0;
        }
        Message::ShowHelp => {
            app.mode = AppMode::Help;
            app.show_help = true;
        }
        Message::HideHelp => {
            app.mode = AppMode::Normal;
            app.show_help = false;
        }
        Message::GenerateConfig => {
            // ... existing implementation
        }
        Message::ConvertSelected => {
            // ... existing implementation
        }

        // NEW: Wizard messages
        Message::StartBenchmarkWizard => {
            app.start_benchmark_wizard();
            // Discover projects
            let projects_dir = app.scan_dir.join("example_project");
            app.wizard.discovered_projects = crate::benchmark::discover_projects(&projects_dir);
            // Discover tronko-assign binaries
            app.wizard.tronko_binaries = discover_tronko_binaries(&app.scan_dir);
        }
        Message::WizardNext => wizard_advance(app),
        Message::WizardPrev => wizard_go_back(app),
        Message::WizardCancel => app.cancel_wizard(),
        Message::WizardSelectItem => wizard_select_item(app),
        Message::WizardConfirm => {
            if let Some(result) = run_benchmark(app) {
                app.benchmark_store.add(result);
                let _ = app.save_benchmarks();
                app.status_message = Some("Benchmark completed and saved".to_string());
            }
            app.cancel_wizard();
        }

        // NEW: Comparison messages
        Message::ToggleBenchmarkSelection => app.toggle_benchmark_selection(),
        Message::ExitComparison => app.exit_comparison(),
    }
    Ok(())
}

fn wizard_select_next(app: &mut App) {
    let max = wizard_current_list_len(app);
    if max > 0 {
        app.wizard.list_index = (app.wizard.list_index + 1).min(max - 1);
    }
}

fn wizard_select_prev(app: &mut App) {
    app.wizard.list_index = app.wizard.list_index.saturating_sub(1);
}

fn wizard_current_list_len(app: &App) -> usize {
    match app.wizard.step {
        WizardStep::SelectDatabase => app.scan_results.databases.len(),
        WizardStep::SelectProject => app.wizard.discovered_projects.len(),
        WizardStep::SelectMarker => {
            app.wizard.selected_project_index
                .and_then(|i| app.wizard.discovered_projects.get(i))
                .map(|p| p.markers.len())
                .unwrap_or(0)
        }
        WizardStep::SelectQueryMode => {
            app.wizard.selected_project_index
                .and_then(|pi| app.wizard.discovered_projects.get(pi))
                .and_then(|p| app.wizard.selected_marker_index.and_then(|mi| p.markers.get(mi)))
                .map(|m| {
                    let mut count = 0;
                    if m.has_paired { count += 1; }
                    if m.has_unpaired_f { count += 1; }
                    if m.has_unpaired_r { count += 1; }
                    count
                })
                .unwrap_or(0)
        }
        WizardStep::SelectBinary => app.wizard.tronko_binaries.len(),
        WizardStep::Confirm => 0,
    }
}

fn wizard_select_item(app: &mut App) {
    use crate::tui::app::WizardStep;
    use crate::benchmark::QueryMode;

    match app.wizard.step {
        WizardStep::SelectDatabase => {
            app.wizard.selected_db_index = Some(app.wizard.list_index);
            wizard_advance(app);
        }
        WizardStep::SelectProject => {
            app.wizard.selected_project_index = Some(app.wizard.list_index);
            wizard_advance(app);
        }
        WizardStep::SelectMarker => {
            app.wizard.selected_marker_index = Some(app.wizard.list_index);
            wizard_advance(app);
        }
        WizardStep::SelectQueryMode => {
            // Map list index to QueryMode based on what's available
            let marker = app.wizard.selected_project_index
                .and_then(|pi| app.wizard.discovered_projects.get(pi))
                .and_then(|p| app.wizard.selected_marker_index.and_then(|mi| p.markers.get(mi)));

            if let Some(m) = marker {
                let mut modes = Vec::new();
                if m.has_paired { modes.push(QueryMode::Paired); }
                if m.has_unpaired_f { modes.push(QueryMode::UnpairedF); }
                if m.has_unpaired_r { modes.push(QueryMode::UnpairedR); }

                if let Some(mode) = modes.get(app.wizard.list_index) {
                    app.wizard.selected_query_mode = Some(*mode);
                    wizard_advance(app);
                }
            }
        }
        WizardStep::SelectBinary => {
            app.wizard.selected_binary_index = Some(app.wizard.list_index);
            wizard_advance(app);
        }
        WizardStep::Confirm => {
            // Run the benchmark
            if let Some(result) = run_benchmark(app) {
                app.benchmark_store.add(result);
                let _ = app.save_benchmarks();
                app.status_message = Some("Benchmark completed and saved".to_string());
            }
            app.cancel_wizard();
        }
    }
}

fn wizard_advance(app: &mut App) {
    if let Some(next) = app.wizard.step.next() {
        app.wizard.step = next;
        app.wizard.list_index = 0;
    }
}

fn wizard_go_back(app: &mut App) {
    if let Some(prev) = app.wizard.step.prev() {
        app.wizard.step = prev;
        app.wizard.list_index = 0;
    }
}

fn discover_tronko_binaries(scan_dir: &Path) -> Vec<TronkoBinary> {
    let mut binaries = Vec::new();

    // Look in common locations
    let candidates = [
        scan_dir.join("tronko-assign/tronko-assign"),
        scan_dir.join("tronko-assign/tronko-assign-optimized"),
        PathBuf::from("./tronko-assign/tronko-assign"),
    ];

    for path in candidates {
        if path.exists() && path.is_file() {
            // Check if it's the optimized version by looking at the name or running --help
            let optimize_memory = path.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.contains("optimized"))
                .unwrap_or(false);

            let label = if optimize_memory {
                format!("{} (OPTIMIZE_MEMORY)", path.file_name().unwrap_or_default().to_string_lossy())
            } else {
                format!("{} (standard)", path.file_name().unwrap_or_default().to_string_lossy())
            };

            binaries.push(TronkoBinary {
                path: path.clone(),
                optimize_memory,
                label,
            });
        }
    }

    binaries
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] All keyboard shortcuts work in appropriate modes
- [ ] Mode transitions work correctly (normal → wizard → back to normal)
- [ ] Wizard navigation with arrow keys and Enter works

---

## Phase 9: Benchmark Execution

### Overview
Implement the actual benchmark execution logic that runs tronko-assign and captures results.

### Changes Required:

#### 1. Add benchmark execution function
**File**: `rust/crates/tronko-bench/src/benchmark.rs`
**Changes**: Add execution and result parsing functions

```rust
use std::process::Command;
use std::time::Instant;
use std::io::{BufRead, BufReader};
use std::collections::HashSet;

/// Run a benchmark with the given configuration
pub fn execute_benchmark(config: BenchmarkConfig) -> BenchmarkResult {
    let timestamp = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    // Create output paths
    let output_dir = std::env::temp_dir().join("tronko-bench");
    std::fs::create_dir_all(&output_dir).ok();

    let output_path = output_dir.join(format!("output_{}.txt", config.id));
    let tsv_log_path = output_dir.join(format!("memlog_{}.tsv", config.id));

    // Build command
    let mut cmd = Command::new(&config.tronko_assign_path);
    cmd.arg("-r")
       .arg("-f").arg(&config.reference_db_path)
       .arg("-a").arg(&config.fasta_path)
       .arg("-o").arg(&output_path)
       .arg("-R")  // Resource monitoring
       .arg("-T")  // Timing
       .arg("--tsv-log").arg(&tsv_log_path);

    if config.paired_end {
        cmd.arg("-p")
           .arg("-1").arg(&config.query_paths[0])
           .arg("-2").arg(&config.query_paths[1]);
    } else {
        cmd.arg("-s")
           .arg("-g").arg(&config.query_paths[0]);
    }

    // Run and time
    let start = Instant::now();
    let output = cmd.output();
    let duration = start.elapsed();

    match output {
        Ok(out) => {
            let success = out.status.success();
            let error_message = if !success {
                Some(String::from_utf8_lossy(&out.stderr).to_string())
            } else {
                None
            };

            // Parse peak memory from tsv log
            let peak_memory_bytes = parse_peak_memory(&tsv_log_path);

            // Get output file size and summary
            let output_size_bytes = std::fs::metadata(&output_path).ok().map(|m| m.len());
            let output_summary = if success {
                parse_output_summary(&output_path)
            } else {
                None
            };

            BenchmarkResult {
                config,
                timestamp,
                duration_secs: duration.as_secs_f64(),
                peak_memory_bytes,
                reads_processed: None,
                output_path,
                tsv_log_path: Some(tsv_log_path),
                success,
                error_message,
                output_size_bytes,
                output_summary,
            }
        }
        Err(e) => {
            BenchmarkResult {
                config,
                timestamp,
                duration_secs: duration.as_secs_f64(),
                peak_memory_bytes: None,
                reads_processed: None,
                output_path,
                tsv_log_path: None,
                success: false,
                error_message: Some(format!("Failed to execute: {}", e)),
                output_size_bytes: None,
                output_summary: None,
            }
        }
    }
}

fn parse_peak_memory(tsv_path: &PathBuf) -> Option<u64> {
    let file = std::fs::File::open(tsv_path).ok()?;
    let reader = BufReader::new(file);

    let mut peak = 0u64;
    for line in reader.lines().skip(1) {  // Skip header
        if let Ok(line) = line {
            let parts: Vec<&str> = line.split('\t').collect();
            // Assuming format: timestamp, rss_kb, ...
            if let Some(rss_str) = parts.get(1) {
                if let Ok(rss_kb) = rss_str.parse::<u64>() {
                    peak = peak.max(rss_kb * 1024);  // Convert KB to bytes
                }
            }
        }
    }

    if peak > 0 { Some(peak) } else { None }
}

fn parse_output_summary(output_path: &PathBuf) -> Option<OutputSummary> {
    let file = std::fs::File::open(output_path).ok()?;
    let reader = BufReader::new(file);

    let mut total_lines = 0u64;
    let mut assigned_lines = 0u64;
    let mut taxa = HashSet::new();
    let mut assignments = Vec::new();

    for line in reader.lines() {
        if let Ok(line) = line {
            total_lines += 1;

            // Parse TSV: read_id, taxonomy, score, etc.
            let parts: Vec<&str> = line.split('\t').collect();
            if let Some(taxonomy) = parts.get(1) {
                if !taxonomy.is_empty() && *taxonomy != "NA" && *taxonomy != "unassigned" {
                    assigned_lines += 1;
                    taxa.insert(taxonomy.to_string());
                    assignments.push(taxonomy.to_string());
                }
            }
        }
    }

    // Create a hash of sorted assignments for comparison
    assignments.sort();
    let hash_input = assignments.join("\n");
    let assignments_hash = format!("{:x}", md5::compute(hash_input.as_bytes()));

    Some(OutputSummary {
        total_lines,
        assigned_lines,
        unique_taxa: taxa.len() as u64,
        assignments_hash,
    })
}
```

#### 2. Add md5 dependency
**File**: `rust/crates/tronko-bench/Cargo.toml`
**Changes**: Add md5 and chrono dependencies

```toml
[dependencies]
# ... existing deps
md5 = "0.7"
chrono = "0.4"
```

#### 3. Wire up execution in mod.rs
**File**: `rust/crates/tronko-bench/src/tui/mod.rs`
**Changes**: Add run_benchmark function

```rust
use crate::benchmark::{BenchmarkConfig, BenchmarkResult, execute_benchmark, QueryMode};
use uuid::Uuid;

fn run_benchmark(app: &App) -> Option<BenchmarkResult> {
    // Build config from wizard state
    let db_idx = app.wizard.selected_db_index?;
    let db = app.scan_results.databases.get(db_idx)?;

    let proj_idx = app.wizard.selected_project_index?;
    let proj = app.wizard.discovered_projects.get(proj_idx)?;

    let marker_idx = app.wizard.selected_marker_index?;
    let marker = proj.markers.get(marker_idx)?;

    let query_mode = app.wizard.selected_query_mode?;
    let bin_idx = app.wizard.selected_binary_index?;
    let binary = app.wizard.tronko_binaries.get(bin_idx)?;

    // Get query paths based on mode
    let (query_paths, paired_end) = match query_mode {
        QueryMode::Paired => {
            let f = marker.paired_f.clone()?;
            let r = marker.paired_r.clone()?;
            (vec![f, r], true)
        }
        QueryMode::UnpairedF => {
            let f = marker.unpaired_f.clone()?;
            (vec![f], false)
        }
        QueryMode::UnpairedR => {
            let r = marker.unpaired_r.clone()?;
            (vec![r], false)
        }
    };

    // Find FASTA file for BWA
    let fasta_path = crate::find_associated_fasta(&db.path)?;

    let config = BenchmarkConfig {
        id: Uuid::new_v4().to_string(),
        name: format!("{}/{} - {} {}",
            db.context,
            marker.name,
            if paired_end { "paired" } else { "single" },
            if binary.optimize_memory { "(opt)" } else { "" }
        ),
        reference_db_path: db.path.clone(),
        db_format: db.format.clone(),
        num_trees: db.num_trees.parse().unwrap_or(0),
        db_size_bytes: db.size_bytes,
        fasta_path,
        query_paths,
        paired_end,
        marker: marker.name.clone(),
        project: proj.name.clone(),
        optimize_memory: binary.optimize_memory,
        tronko_assign_path: binary.path.clone(),
    };

    Some(execute_benchmark(config))
}
```

#### 4. Add uuid dependency
**File**: `rust/crates/tronko-bench/Cargo.toml`
**Changes**: Add uuid dependency

```toml
[dependencies]
# ... existing deps
uuid = { version = "1", features = ["v4"] }
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Completing wizard runs tronko-assign
- [ ] Results are captured (duration, memory, output summary)
- [ ] Results are persisted to JSON file
- [ ] New benchmark appears in Benchmarks tab

---

## Phase 10: Help and Status Updates

### Overview
Update help text and status bar to reflect new functionality.

### Changes Required:

#### 1. Update help overlay
**File**: `rust/crates/tronko-bench/src/tui/ui/help.rs`
**Changes**: Add new keyboard shortcuts

```rust
let help_text = vec![
    Line::from(vec![
        Span::styled("tronko-bench TUI", Style::default().add_modifier(Modifier::BOLD)),
    ]),
    Line::from(""),
    Line::from(Span::styled("Navigation", Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan))),
    Line::from(vec![
        Span::styled("Tab/Shift+Tab ", Style::default().fg(Color::Cyan)),
        Span::raw("Switch tabs"),
    ]),
    Line::from(vec![
        Span::styled("j/k or ↑↓    ", Style::default().fg(Color::Cyan)),
        Span::raw("Navigate list"),
    ]),
    Line::from(vec![
        Span::styled("/            ", Style::default().fg(Color::Cyan)),
        Span::raw("Fuzzy search"),
    ]),
    Line::from(""),
    Line::from(Span::styled("Actions", Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan))),
    Line::from(vec![
        Span::styled("b            ", Style::default().fg(Color::Cyan)),
        Span::raw("New benchmark wizard"),
    ]),
    Line::from(vec![
        Span::styled("Space        ", Style::default().fg(Color::Cyan)),
        Span::raw("Select benchmark for comparison (max 2)"),
    ]),
    Line::from(vec![
        Span::styled("c            ", Style::default().fg(Color::Cyan)),
        Span::raw("Convert database to binary"),
    ]),
    Line::from(vec![
        Span::styled("g            ", Style::default().fg(Color::Cyan)),
        Span::raw("Generate benchmark.json"),
    ]),
    Line::from(vec![
        Span::styled("r            ", Style::default().fg(Color::Cyan)),
        Span::raw("Rescan directory"),
    ]),
    Line::from(""),
    Line::from(Span::styled("Other", Style::default().add_modifier(Modifier::BOLD).fg(Color::Cyan))),
    Line::from(vec![
        Span::styled("?            ", Style::default().fg(Color::Cyan)),
        Span::raw("Show this help"),
    ]),
    Line::from(vec![
        Span::styled("q/Esc        ", Style::default().fg(Color::Cyan)),
        Span::raw("Quit / Close overlay"),
    ]),
    Line::from(""),
    Line::from(Span::styled(
        "Press any key to close",
        Style::default().fg(Color::DarkGray),
    )),
];
```

#### 2. Update status bar
**File**: `rust/crates/tronko-bench/src/tui/ui/status.rs`
**Changes**: Update hints based on current tab

```rust
pub fn render_status(f: &mut Frame, area: Rect, app: &App) {
    let status = if app.loading {
        Line::from(vec![
            Span::styled(" Scanning...", Style::default().fg(Color::Yellow)),
        ])
    } else if let Some(ref msg) = app.status_message {
        Line::from(vec![
            Span::styled(format!(" {}", msg), Style::default().fg(Color::Cyan)),
        ])
    } else {
        let hints = match app.current_tab {
            Tab::Databases => "Tab:switch  j/k:nav  /:search  b:benchmark  c:convert  r:rescan  ?:help",
            Tab::Queries | Tab::MsaFiles => "Tab:switch  j/k:nav  /:search  r:rescan  ?:help",
            Tab::Benchmarks => "Tab:switch  j/k:nav  Space:select  b:new  ?:help",
        };

        Line::from(vec![
            Span::raw(format!(" {} ", app.scan_dir.display())),
            Span::styled(" | ", Style::default().fg(Color::DarkGray)),
            Span::raw(hints),
        ])
    };

    let paragraph = Paragraph::new(status)
        .style(Style::default().bg(Color::DarkGray));

    f.render_widget(paragraph, area);
}
```

### Success Criteria:

#### Automated Verification:
- [ ] Code compiles: `cd rust && cargo build -p tronko-bench`

#### Manual Verification:
- [ ] Help screen shows all new keyboard shortcuts
- [ ] Status bar shows context-appropriate hints for each tab

---

## Testing Strategy

### Unit Tests:
- Test `discover_projects()` with mock directory structure
- Test `parse_output_summary()` with sample output files
- Test `BenchmarkStore` save/load roundtrip

### Integration Tests:
- Run full benchmark wizard workflow
- Verify benchmark results persist across app restarts
- Test comparison mode with two benchmarks

### Manual Testing Steps:
1. Launch TUI with `tronko-bench tui -d /path/to/tronko`
2. Navigate to Databases tab, verify context shown
3. Press 'b' to open benchmark wizard
4. Select a database, project, marker, query mode, and binary
5. Confirm and verify benchmark runs
6. Navigate to Benchmarks tab, verify result appears
7. Run a second benchmark with different settings
8. Select both benchmarks with spacebar
9. Verify comparison view shows correct metrics
10. Restart app, verify benchmarks persist

## Performance Considerations

- Benchmark execution is blocking; consider async for large datasets
- Output parsing should handle multi-GB files efficiently (line-by-line streaming)
- Benchmark store JSON may grow large; consider pagination or separate files per benchmark

## References

- Existing TUI implementation: `rust/crates/tronko-bench/src/tui/`
- tronko-assign options: `tronko-assign/options.c`
- OPTIMIZE_MEMORY flag: `tronko-assign/Makefile:16-17`
