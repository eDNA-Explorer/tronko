# tronko-bench Benchmarking Workflow Plan

## Overview

Extend the tronko-bench TUI to support a complete benchmarking workflow: building benchmarks with configurable parameters, storing results, viewing historical runs, and comparing up to two benchmarks side-by-side.

## Current State

The TUI currently has three tabs (Databases, Queries, MSA Files) with basic navigation, fuzzy search, database conversion, and config generation. It does not yet support:
- Running actual benchmarks
- Storing/retrieving benchmark results
- Comparing benchmark runs

## Goals

1. **Build Benchmark Workflow**: Guide users through selecting a reference database, query files (with paired/unpaired detection), and optimization flags
2. **Result Storage**: Persist benchmark metadata and performance metrics to disk
3. **History View**: Display past benchmarks with full details and performance summaries
4. **Comparison Mode**: Select up to two benchmarks to compare speed, memory, and result differences

## Non-Goals

- Real-time progress streaming during benchmark execution (blocking is acceptable for v1)
- Automated benchmark scheduling/batching
- Cloud storage of results
- Graphical charts (text-based summaries only)

---

## Data Model

### Benchmark Configuration

```rust
/// A complete benchmark configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    /// Unique identifier (UUID v4)
    pub id: String,

    /// Human-readable name (auto-generated or user-provided)
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

    /// Database format: "text", "text (gzipped)", "binary (.trkb)"
    pub format: String,

    /// Number of trees in database
    pub num_trees: u32,

    /// Size on disk in bytes
    pub size_bytes: u64,

    /// Parent directory name (for display context)
    pub parent_dir: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryConfig {
    /// Base directory containing query files
    pub base_dir: PathBuf,

    /// Project name (first directory component)
    pub project: String,

    /// Marker name (e.g., "16S", "ITS")
    pub marker: String,

    /// Read mode
    pub read_mode: ReadMode,

    /// Paths to query files
    pub files: QueryFiles,

    /// Total sequences across all files
    pub total_sequences: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ReadMode {
    SingleEnd,
    PairedEnd,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum QueryFiles {
    Single { path: PathBuf },
    Paired { forward: PathBuf, reverse: PathBuf },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExecutionParams {
    /// Use OPTIMIZE_MEMORY build (float vs double precision)
    pub optimize_memory: bool,

    /// Number of CPU cores to use
    pub cores: u32,

    /// Use Needleman-Wunsch instead of WFA
    pub use_needleman_wunsch: bool,

    /// Batch size (lines per batch)
    pub batch_size: u32,
}
```

### Benchmark Results

```rust
/// Results from a completed benchmark run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    /// Reference to the config that produced this result
    pub config_id: String,

    /// When the benchmark started
    pub started_at: DateTime<Utc>,

    /// When the benchmark completed
    pub completed_at: DateTime<Utc>,

    /// Exit status
    pub status: BenchmarkStatus,

    /// Performance metrics
    pub metrics: Option<PerformanceMetrics>,

    /// Path to the output assignments file
    pub output_path: PathBuf,

    /// Path to the TSV memory log
    pub metrics_path: PathBuf,

    /// Assignment summary statistics
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

    /// Peak RSS memory in MB
    pub peak_rss_mb: f64,

    /// Final RSS memory in MB
    pub final_rss_mb: f64,

    /// CPU user time in seconds
    pub cpu_user_secs: f64,

    /// CPU system time in seconds
    pub cpu_sys_secs: f64,

    /// Memory at tree allocation phase
    pub tree_alloc_rss_mb: f64,

    /// Memory after reference loaded
    pub reference_loaded_rss_mb: f64,

    /// Number of batches processed
    pub batches_processed: u32,

    /// Reads processed per second
    pub reads_per_second: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssignmentStats {
    /// Total reads processed
    pub total_reads: u64,

    /// Reads successfully assigned
    pub assigned_reads: u64,

    /// Reads marked as unassigned
    pub unassigned_reads: u64,

    /// Unique taxonomic paths found
    pub unique_taxa: u64,

    /// Distribution of assignments by tree
    pub tree_distribution: HashMap<u32, u64>,
}
```

### Storage Format

Benchmarks will be stored in `~/.tronko-bench/` or a configurable directory:

```
~/.tronko-bench/
├── config.toml              # Global settings
├── benchmarks/
│   ├── index.json           # Quick-access index of all benchmarks
│   ├── {uuid}/
│   │   ├── config.json      # BenchmarkConfig
│   │   ├── result.json      # BenchmarkResult (after completion)
│   │   ├── output.tsv       # tronko-assign output
│   │   ├── metrics.tsv      # --tsv-log output
│   │   └── stderr.log       # Captured stderr
```

---

## UI Design

### New Tab Structure

Replace the current 3-tab layout with a mode-based navigation:

```
┌─────────────────────────────────────────────────────────────────┐
│ [Benchmarks]  [New Benchmark]  [Compare]                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  (Content varies by mode)                                       │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│ Status bar with context-sensitive hints                         │
└─────────────────────────────────────────────────────────────────┘
```

### Mode 1: Benchmarks List (Default)

Shows all past benchmark runs with key metrics:

```
┌─ Benchmarks ────────────────────────────┬─ Details ─────────────────┐
│ > 16S_Bacteria paired (binary) 2m ago   │ Name: 16S_Bacteria paired │
│   16S_Bacteria paired (text)   5m ago   │ Database: 16S_Bacteria/   │
│   ITS_Fungi single (binary)    1h ago   │   reference_tree.trkb     │
│   16S_Bacteria paired (optmem) 2h ago   │   Format: binary          │
│                                         │   Trees: 156              │
│                                         │   Size: 38.08 MiB         │
│                                         │                           │
│                                         │ Queries: example_proj/16S │
│                                         │   Mode: paired-end        │
│                                         │   Sequences: 10,000       │
│                                         │                           │
│                                         │ Performance:              │
│                                         │   Time: 45.2s             │
│                                         │   Peak Memory: 64.8 MB    │
│                                         │   Reads/sec: 221.2        │
│                                         │                           │
│                                         │ Results:                  │
│                                         │   Assigned: 9,847 (98.5%) │
│                                         │   Unassigned: 153 (1.5%)  │
├─────────────────────────────────────────┴───────────────────────────┤
│ j/k:nav  Space:select  Enter:view  n:new  d:delete  ?:help  q:quit  │
└─────────────────────────────────────────────────────────────────────┘
```

List item format: `{marker}_{db_name} {read_mode} ({format}{+optmem?}) {relative_time}`

### Mode 2: New Benchmark Wizard

A multi-step wizard flow:

**Step 1: Select Reference Database**
```
┌─ New Benchmark: Select Database ────────┬─ Database Info ───────────┐
│ Filter: ____________________________    │                           │
│                                         │ Path: /data/16S_Bacteria/ │
│ > 16S_Bacteria/reference_tree.trkb     │   reference_tree.trkb     │
│   16S_Bacteria/reference_tree.txt.gz   │                           │
│   ITS_Fungi/reference_tree.trkb        │ Format: binary (.trkb)    │
│   18S_Eukaryotes/reference_tree.txt    │ Trees: 156                │
│   COI_Animals/reference_tree.trkb      │ Size: 38.08 MiB           │
│                                         │                           │
│                                         │ Associated FASTA:         │
│                                         │   16S_Bacteria.fasta      │
│                                         │   (auto-detected)         │
├─────────────────────────────────────────┴───────────────────────────┤
│ Step 1/4  j/k:nav  Enter:select  Esc:cancel  /:search               │
└─────────────────────────────────────────────────────────────────────┘
```

**Step 2: Select Query Directory**
```
┌─ New Benchmark: Select Queries ─────────┬─ Query Info ──────────────┐
│ Filter: ____________________________    │                           │
│                                         │ Project: example_project  │
│ > example_project/16S/paired/          │ Marker: 16S               │
│   example_project/16S/unpaired/        │ Type: paired-end          │
│   example_project/ITS/paired/          │                           │
│   test_data/16S/single/                │ Files:                    │
│   validation/16S/paired/               │   R1.fastq.gz (5,000 seq) │
│                                         │   R2.fastq.gz (5,000 seq) │
│                                         │                           │
│                                         │ Total: 10,000 sequences   │
├─────────────────────────────────────────┴───────────────────────────┤
│ Step 2/4  j/k:nav  Enter:select  Esc:back  /:search                 │
└─────────────────────────────────────────────────────────────────────┘
```

**Step 3: Configure Parameters**
```
┌─ New Benchmark: Configuration ──────────────────────────────────────┐
│                                                                     │
│   Memory Optimization:  [ ] Enable (uses float instead of double)   │
│                                                                     │
│   CPU Cores:            [4____]                                     │
│                                                                     │
│   Alignment Algorithm:  (•) WFA  ( ) Needleman-Wunsch               │
│                                                                     │
│   Batch Size:           [50000_]                                    │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│ Step 3/4  Tab:next field  Space:toggle  Enter:continue  Esc:back    │
└─────────────────────────────────────────────────────────────────────┘
```

**Step 4: Confirm and Run**
```
┌─ New Benchmark: Confirm ────────────────────────────────────────────┐
│                                                                     │
│   Database:    16S_Bacteria/reference_tree.trkb (binary, 156 trees) │
│   FASTA:       16S_Bacteria/16S_Bacteria.fasta                      │
│   Queries:     example_project/16S/paired/ (10,000 paired reads)    │
│   Memory Opt:  No                                                   │
│   Cores:       4                                                    │
│   Algorithm:   WFA                                                  │
│                                                                     │
│   Command:                                                          │
│   tronko-assign -r -f reference_tree.trkb -a 16S.fasta \            │
│     -p -1 R1.fastq.gz -2 R2.fastq.gz -C 4 -o output.tsv             │
│                                                                     │
│                    [ Run Benchmark ]    [ Cancel ]                  │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│ Step 4/4  Enter:run  Esc:back  Tab:switch button                    │
└─────────────────────────────────────────────────────────────────────┘
```

**Running State**
```
┌─ Running Benchmark ─────────────────────────────────────────────────┐
│                                                                     │
│   ████████████████████░░░░░░░░░░░░░░░░░░░░  42%                     │
│                                                                     │
│   Status: Processing batch 21/50                                    │
│   Elapsed: 00:01:23                                                 │
│   Memory: 64.2 MB                                                   │
│                                                                     │
│   Live output:                                                      │
│   [INFO] Loaded tree 0 with 2931 nodes                              │
│   [INFO] BWA index built successfully                               │
│   [INFO] Processing batch 21...                                     │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│ Esc:cancel (will stop tronko-assign)                                │
└─────────────────────────────────────────────────────────────────────┘
```

### Mode 3: Comparison View

When two benchmarks are selected (via Space in list view):

```
┌─ Comparison ────────────────────────────────────────────────────────┐
│                                                                     │
│  Benchmark A                      vs    Benchmark B                 │
│  16S_Bacteria paired (binary)           16S_Bacteria paired (text)  │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  PERFORMANCE                                                        │
│  ───────────                                                        │
│  Wall Time:     45.2s                   89.7s (+98.9%)              │
│  Peak Memory:   64.8 MB                 128.4 MB (+98.1%)           │
│  Reads/sec:     221.2                   111.5 (-49.6%)              │
│                                                                     │
│  Winner: Benchmark A (faster, less memory)                          │
│                                                                     │
│  RESULTS                                                            │
│  ───────                                                            │
│  Assigned:      9,847 (98.5%)           9,847 (98.5%)               │
│  Unassigned:    153 (1.5%)              153 (1.5%)                  │
│                                                                     │
│  Assignment Comparison:                                             │
│    Identical assignments: 9,847 (100%)                              │
│    Different assignments: 0 (0%)                                    │
│                                                                     │
│  Score Comparison:                                                  │
│    Avg score diff: 0.0023                                           │
│    Max score diff: 0.0089                                           │
│    Scores within 0.01: 9,721 (98.7%)                                │
│                                                                     │
│  VERDICT: Results are identical. Binary format is 2x faster.        │
│                                                                     │
├─────────────────────────────────────────────────────────────────────┤
│ Esc:back to list  d:detailed diff  e:export comparison              │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Implementation Phases

### Phase 1: Data Model & Storage Layer

**Files to create:**
- `rust/crates/tronko-bench/src/benchmark/mod.rs` - Module exports
- `rust/crates/tronko-bench/src/benchmark/config.rs` - BenchmarkConfig, DatabaseConfig, QueryConfig
- `rust/crates/tronko-bench/src/benchmark/result.rs` - BenchmarkResult, PerformanceMetrics
- `rust/crates/tronko-bench/src/benchmark/storage.rs` - Load/save benchmarks to disk
- `rust/crates/tronko-bench/src/benchmark/parser.rs` - Parse TSV metrics and assignment output

**Dependencies to add:**
```toml
uuid = { version = "1", features = ["v4", "serde"] }
chrono = { version = "0.4", features = ["serde"] }
directories = "5"  # For ~/.tronko-bench location
```

**Key functions:**
```rust
impl BenchmarkStorage {
    pub fn new() -> Result<Self>;
    pub fn list_benchmarks(&self) -> Result<Vec<BenchmarkSummary>>;
    pub fn get_benchmark(&self, id: &str) -> Result<(BenchmarkConfig, Option<BenchmarkResult>)>;
    pub fn save_config(&self, config: &BenchmarkConfig) -> Result<()>;
    pub fn save_result(&self, result: &BenchmarkResult) -> Result<()>;
    pub fn delete_benchmark(&self, id: &str) -> Result<()>;
}

impl MetricsParser {
    pub fn parse_tsv_log(path: &Path) -> Result<PerformanceMetrics>;
    pub fn parse_assignments(path: &Path) -> Result<AssignmentStats>;
}
```

**Success Criteria:**
- [ ] Can create and serialize BenchmarkConfig to JSON
- [ ] Can parse tronko-assign TSV metrics log
- [ ] Can parse assignment output and compute stats
- [ ] Storage directory created with proper structure
- [ ] Index file maintained for quick listing

### Phase 2: Query Directory Detection

**Files to modify:**
- `rust/crates/tronko-bench/src/main.rs` - Add query directory scanning
- `rust/crates/tronko-bench/src/benchmark/discovery.rs` - New file for discovery logic

**Query directory detection logic:**

```rust
/// Scan for query directories with project/marker/mode structure
pub fn scan_query_directories(root: &Path, max_depth: usize) -> Vec<QueryDirectory> {
    // Look for directories containing FASTA/FASTQ files
    // Detect paired vs single by:
    //   - "paired" or "unpaired" in path
    //   - Presence of R1/R2, _1/_2, forward/reverse file pairs
    //   - Single file = single-end
}

#[derive(Debug, Clone)]
pub struct QueryDirectory {
    pub path: PathBuf,
    pub project: String,       // First directory component
    pub marker: String,        // e.g., "16S", "ITS"
    pub read_mode: ReadMode,
    pub files: Vec<QueryFile>,
    pub total_sequences: u64,
}

#[derive(Debug, Clone)]
pub struct QueryFile {
    pub path: PathBuf,
    pub name: String,
    pub file_type: String,     // "FASTA" or "FASTQ"
    pub compressed: bool,      // .gz
    pub sequences: u64,
    pub pair_id: Option<PairId>, // R1, R2, or None
}

#[derive(Debug, Clone)]
pub enum PairId { Forward, Reverse }
```

**Paired-end detection heuristics:**
1. Directory name contains "paired" → PairedEnd
2. Directory name contains "unpaired" or "single" → SingleEnd
3. Files matching patterns: `*_R1*`/`*_R2*`, `*_1.*`/`*_2.*`, `*forward*`/`*reverse*` → PairedEnd
4. Single FASTA/FASTQ file → SingleEnd

**Success Criteria:**
- [ ] Detects paired vs single-end directories correctly
- [ ] Extracts project and marker from path
- [ ] Counts sequences in each file (including gzipped)
- [ ] Groups paired files correctly

### Phase 3: Enhanced Database Display

**Files to modify:**
- `rust/crates/tronko-bench/src/tui/ui/list.rs` - Update database list rendering

**Display format for databases:**
```
{parent_dir}/{filename}  {num_trees} trees  {size}  [{format}]
```

Examples:
```
16S_Bacteria/reference_tree.trkb     156 trees   38.08 MiB  [BIN]
16S_Bacteria/reference_tree.txt.gz   156 trees  102.34 MiB  [TXT.GZ]
ITS_Fungi/reference_tree.txt         89 trees    24.12 MiB  [TXT]
```

**Changes to DatabaseInfo:**
```rust
pub struct DatabaseInfo {
    pub name: String,
    pub parent_dir: String,      // NEW: Parent directory name
    pub format: String,
    pub num_trees: String,
    pub size: String,
    pub size_bytes: u64,         // NEW: Raw size for sorting
    pub path: PathBuf,
}
```

**Success Criteria:**
- [ ] Database list shows parent directory for context
- [ ] Format clearly indicates binary/text/gzipped
- [ ] Can sort by size, trees, or name

### Phase 4: TUI Mode System Refactor

**Files to modify:**
- `rust/crates/tronko-bench/src/tui/app.rs` - New mode-based state
- `rust/crates/tronko-bench/src/tui/message.rs` - New messages
- `rust/crates/tronko-bench/src/tui/event.rs` - Mode-aware key handling
- `rust/crates/tronko-bench/src/tui/ui/mod.rs` - Mode-based rendering

**New App state:**

```rust
pub enum AppMode {
    /// Viewing list of past benchmarks
    BenchmarkList {
        benchmarks: Vec<BenchmarkSummary>,
        selected_index: usize,
        selected_for_compare: Vec<String>,  // Up to 2 benchmark IDs
    },

    /// Creating a new benchmark
    NewBenchmark(WizardState),

    /// Running a benchmark
    Running {
        config: BenchmarkConfig,
        start_time: Instant,
        output_lines: Vec<String>,
        child_process: Option<Child>,
    },

    /// Comparing two benchmarks
    Comparison {
        benchmark_a: (BenchmarkConfig, BenchmarkResult),
        benchmark_b: (BenchmarkConfig, BenchmarkResult),
        diff: Option<AssignmentDiff>,
    },
}

pub enum WizardState {
    SelectDatabase {
        databases: Vec<DatabaseInfo>,
        selected_index: usize,
        filter: String,
    },
    SelectQueries {
        selected_db: DatabaseConfig,
        query_dirs: Vec<QueryDirectory>,
        selected_index: usize,
        filter: String,
    },
    ConfigureParams {
        selected_db: DatabaseConfig,
        selected_queries: QueryConfig,
        params: ExecutionParams,
        focused_field: usize,
    },
    Confirm {
        config: BenchmarkConfig,
        focused_button: usize,  // 0 = Run, 1 = Cancel
    },
}
```

**New messages:**
```rust
pub enum Message {
    // ... existing messages ...

    // Mode transitions
    EnterNewBenchmark,
    ExitWizard,

    // Wizard navigation
    WizardNext,
    WizardBack,
    WizardSelectItem,
    WizardToggleParam,
    WizardUpdateField(String),

    // Benchmark actions
    StartBenchmark,
    CancelBenchmark,
    BenchmarkComplete(BenchmarkResult),
    BenchmarkFailed(String),

    // Comparison
    ToggleSelectForCompare,
    EnterComparison,
    ExitComparison,

    // List actions
    DeleteBenchmark,
    ViewBenchmarkDetails,
}
```

**Success Criteria:**
- [ ] Mode transitions work correctly
- [ ] Wizard flow navigates forward/backward
- [ ] Esc cancels wizard and returns to list
- [ ] State is preserved when navigating back in wizard

### Phase 5: Benchmark Execution

**Files to create:**
- `rust/crates/tronko-bench/src/benchmark/runner.rs` - Execute tronko-assign

**tronko-assign invocation:**

```rust
pub struct BenchmarkRunner {
    tronko_assign_path: PathBuf,
    tronko_assign_optimized_path: Option<PathBuf>,
}

impl BenchmarkRunner {
    /// Find tronko-assign binaries
    pub fn new() -> Result<Self>;

    /// Build command line arguments
    pub fn build_command(&self, config: &BenchmarkConfig) -> Command;

    /// Run benchmark and capture output
    pub fn run(&self, config: &BenchmarkConfig) -> Result<BenchmarkResult>;
}
```

**Command building:**
```rust
fn build_command(&self, config: &BenchmarkConfig) -> Command {
    let binary = if config.params.optimize_memory {
        &self.tronko_assign_optimized_path
    } else {
        &self.tronko_assign_path
    };

    let mut cmd = Command::new(binary);
    cmd.arg("-r")  // Reference mode
       .arg("-f").arg(&config.database.reference_path)
       .arg("-a").arg(&config.database.fasta_path)
       .arg("-o").arg(&output_path)
       .arg("--tsv-log").arg(&metrics_path)
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

    // Detect FASTQ vs FASTA
    if is_fastq(&config.queries.files) {
        cmd.arg("-q");
    }

    cmd
}
```

**Success Criteria:**
- [ ] Correctly builds tronko-assign command
- [ ] Captures stdout/stderr during execution
- [ ] Parses metrics TSV after completion
- [ ] Stores results to benchmark directory
- [ ] Handles cancellation gracefully (SIGTERM)

### Phase 6: Comparison Engine

**Files to create:**
- `rust/crates/tronko-bench/src/benchmark/comparison.rs` - Compare two benchmark results

**Comparison logic:**

```rust
#[derive(Debug, Clone, Serialize)]
pub struct BenchmarkComparison {
    pub benchmark_a_id: String,
    pub benchmark_b_id: String,

    pub performance: PerformanceComparison,
    pub assignments: AssignmentComparison,
    pub verdict: ComparisonVerdict,
}

#[derive(Debug, Clone, Serialize)]
pub struct PerformanceComparison {
    pub wall_time_diff_pct: f64,    // Positive = A faster
    pub peak_memory_diff_pct: f64,  // Positive = A uses less
    pub throughput_diff_pct: f64,   // Positive = A faster
    pub winner: Option<char>,       // 'A', 'B', or None (tie)
}

#[derive(Debug, Clone, Serialize)]
pub struct AssignmentComparison {
    /// Reads with identical taxonomic assignments
    pub identical_assignments: u64,

    /// Reads with different taxonomic assignments
    pub different_assignments: u64,

    /// Assignment match percentage
    pub match_percentage: f64,

    /// Score statistics for identically-assigned reads
    pub score_diff_avg: f64,
    pub score_diff_max: f64,
    pub score_diff_std: f64,

    /// Percentage of scores within various thresholds
    pub scores_within_001: f64,
    pub scores_within_01: f64,
}

#[derive(Debug, Clone, Serialize)]
pub enum ComparisonVerdict {
    /// Assignments identical, recommend faster option
    IdenticalRecommendFaster { winner: char },

    /// Assignments differ, need investigation
    AssignmentsDiffer { diff_count: u64, diff_pct: f64 },

    /// Same performance, same results
    Equivalent,
}

impl BenchmarkComparison {
    pub fn compare(
        config_a: &BenchmarkConfig,
        result_a: &BenchmarkResult,
        config_b: &BenchmarkConfig,
        result_b: &BenchmarkResult,
    ) -> Result<Self>;
}
```

**Assignment comparison algorithm:**
1. Parse both output TSV files into `HashMap<read_name, Assignment>`
2. For each read in A, check if B has same taxonomic path
3. For matching paths, compute score difference
4. Generate statistics

**Success Criteria:**
- [ ] Correctly identifies identical vs different assignments
- [ ] Computes score difference statistics
- [ ] Determines performance winner
- [ ] Generates actionable verdict

### Phase 7: UI Polish & Help

**Updates:**
- Add help overlay for each mode
- Add keyboard shortcuts display in status bar
- Add loading indicators
- Add error display

**Key bindings by mode:**

| Mode | Key | Action |
|------|-----|--------|
| List | `n` | New benchmark |
| List | `Enter` | View details |
| List | `Space` | Toggle select for compare |
| List | `c` | Compare selected (when 2 selected) |
| List | `d` | Delete benchmark |
| List | `/` | Filter list |
| Wizard | `Enter` | Select/Next |
| Wizard | `Esc` | Back/Cancel |
| Wizard | `Tab` | Next field |
| Wizard | `Space` | Toggle checkbox |
| Running | `Esc` | Cancel benchmark |
| Compare | `Esc` | Back to list |
| Compare | `e` | Export comparison |
| All | `?` | Help |
| All | `q` | Quit (with confirmation if running) |

---

## Testing Strategy

### Unit Tests
- `benchmark/parser.rs`: Test TSV parsing with sample files
- `benchmark/storage.rs`: Test CRUD operations on temp directory
- `benchmark/comparison.rs`: Test comparison logic with known outputs

### Integration Tests
- Run mini benchmark with test data
- Verify storage persistence across restarts
- Test comparison with known-different outputs

### Manual Testing Checklist
- [ ] Create benchmark with binary database
- [ ] Create benchmark with text database
- [ ] Create benchmark with gzipped database
- [ ] Create benchmark with paired-end reads
- [ ] Create benchmark with single-end reads
- [ ] Toggle memory optimization
- [ ] Cancel running benchmark
- [ ] Compare two successful benchmarks
- [ ] Compare benchmarks with different assignments
- [ ] Delete benchmark
- [ ] Filter benchmark list

---

## File Summary

**New files:**
```
rust/crates/tronko-bench/src/
├── benchmark/
│   ├── mod.rs
│   ├── config.rs
│   ├── result.rs
│   ├── storage.rs
│   ├── parser.rs
│   ├── discovery.rs
│   ├── runner.rs
│   └── comparison.rs
└── tui/
    └── ui/
        ├── wizard.rs      (new)
        ├── running.rs     (new)
        └── comparison.rs  (new)
```

**Modified files:**
```
rust/crates/tronko-bench/
├── Cargo.toml           (add uuid, chrono, directories)
├── src/main.rs          (add benchmark module, update CLI)
└── src/tui/
    ├── app.rs           (new mode-based state)
    ├── message.rs       (new messages)
    ├── event.rs         (mode-aware handling)
    └── ui/
        ├── mod.rs       (mode-based rendering)
        ├── list.rs      (benchmark list + multi-select)
        └── detail.rs    (benchmark details)
```

---

## Dependencies

```toml
[dependencies]
# Existing
ratatui.workspace = true
crossterm.workspace = true
serde.workspace = true
serde_json.workspace = true
anyhow.workspace = true
chrono.workspace = true

# New
uuid = { version = "1", features = ["v4", "serde"] }
directories = "5"
```

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| tronko-assign not found | Scan common locations, prompt user to configure path |
| Benchmark takes too long | Show progress, allow cancellation |
| Large output files | Stream parse, don't load entire file to memory |
| Comparison memory usage | Compare in chunks if files are huge |
| Concurrent modification | Use file locking for index.json |

---

## Success Metrics

1. **Usability**: User can create and run a benchmark in < 1 minute
2. **Accuracy**: Comparison correctly identifies identical vs different assignments
3. **Performance**: TUI remains responsive during benchmark execution
4. **Reliability**: No data loss on crash/cancel
