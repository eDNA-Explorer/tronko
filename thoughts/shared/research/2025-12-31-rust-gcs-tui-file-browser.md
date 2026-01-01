---
date: 2025-12-31T12:00:00-07:00
researcher: Claude
git_commit: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
branch: experimental
repository: tronko
topic: "Building a GCS File Browser TUI in Rust"
tags: [research, rust, gcs, tui, ratatui, google-cloud-storage, fuzzy-search]
status: complete
last_updated: 2025-12-31
last_updated_by: Claude
---

# Research: Building a GCS File Browser TUI in Rust

**Date**: 2025-12-31T12:00:00-07:00
**Researcher**: Claude
**Git Commit**: 77ade9ecf15855e6902a77d6fb11bfabfd7428e0
**Branch**: experimental
**Repository**: tronko

## Research Question

How can we create a simple file browser for Google Cloud Storage (GCS) in Rust with the following features:
- List all buckets
- Select a bucket and drill down into directories
- List directory contents and total size
- Fuzzy search utility
- Select and download files to a customizable local directory (defaults to launch directory)

---

## Summary

Building a GCS file browser TUI in Rust is highly feasible using:
- **Ratatui** - The dominant TUI framework with excellent ecosystem
- **cloud-storage** or **google-cloud-storage** crate - For GCS API interaction
- **nucleo** - High-performance fuzzy matching from the Helix editor team
- **tokio** - Async runtime for non-blocking GCS operations

No existing GCS TUI browser exists in Rust, making this a novel project. The closest analog is **STU** (S3 browser), which provides excellent architectural reference.

---

## Recommended Technology Stack

| Component | Recommendation | Rationale |
|-----------|---------------|-----------|
| **TUI Framework** | `ratatui` | 17k+ stars, excellent docs, active community |
| **GCS Client** | `cloud-storage` | Simpler API, good async support |
| **Fuzzy Search** | `nucleo-matcher` | 6x faster than alternatives, excellent Unicode |
| **Async Runtime** | `tokio` | Standard for async Rust |
| **Terminal Backend** | `crossterm` | Cross-platform, default for ratatui |
| **Tree Widget** | `tui-tree-widget` | Hierarchical navigation |

---

## Detailed Findings

### 1. GCS API Operations

The GCS JSON API v1 provides all necessary endpoints:

| Operation | REST Endpoint | Key Parameters |
|-----------|---------------|----------------|
| List buckets | `GET /storage/v1/b?project=PROJECT` | `prefix`, `maxResults`, `pageToken` |
| List objects | `GET /storage/v1/b/{bucket}/o` | `prefix`, `delimiter`, `maxResults` |
| Get metadata | `GET /storage/v1/b/{bucket}/o/{object}` | `projection=full` |
| Download | `GET /storage/v1/b/{bucket}/o/{object}?alt=media` | Supports Range headers |

**Directory-like browsing** uses `prefix` and `delimiter=/` parameters. The response includes:
- `items[]` - Objects at current level
- `prefixes[]` - "Subdirectories" to drill into

**Authentication options**:
1. Application Default Credentials (ADC) via `gcloud auth application-default login`
2. Service Account JSON file via `GOOGLE_APPLICATION_CREDENTIALS` env var
3. GCE metadata server (when running on Google Cloud)

---

### 2. Rust GCS Client Libraries

#### Recommended: `cloud-storage` crate

```toml
[dependencies]
cloud-storage = "0.11"
```

```rust
use cloud_storage::{Client, Object, ListRequest};

// List buckets
let client = Client::default();
let buckets = client.bucket().list().await?;

// List objects with prefix (directory browsing)
let request = ListRequest {
    prefix: Some("folder/".to_string()),
    delimiter: Some("/".to_string()),
    max_results: Some(100),
    ..Default::default()
};
let objects = Object::list("my-bucket", request).await?;

// Download file
let bytes = Object::download("my-bucket", "path/to/file.txt").await?;

// Stream download for large files
let mut stream = Object::download_streamed("my-bucket", "large-file.bin").await?;
while let Some(byte) = stream.next().await {
    // Write to local file
}
```

#### Alternative: Official `google-cloud-storage` (requires Rust 1.85+)

```toml
[dependencies]
google-cloud-storage = "1.5"
```

```rust
use google_cloud_storage::client::{Storage, StorageControl};

let storage = Storage::builder().build().await?;
let control = StorageControl::builder().build().await?;
```

---

### 3. Ratatui TUI Framework

Ratatui uses **immediate mode rendering** - redraw entire UI each frame. Key concepts:

#### Layouts
```rust
use ratatui::layout::{Layout, Direction, Constraint};

// Three-pane layout: sidebar | main | preview
let layout = Layout::default()
    .direction(Direction::Horizontal)
    .constraints([
        Constraint::Percentage(20),  // Bucket list
        Constraint::Percentage(50),  // File browser
        Constraint::Percentage(30),  // Preview/details
    ])
    .split(frame.area());
```

#### Built-in Widgets
- **List** - Scrollable, selectable items (for bucket/file lists)
- **Table** - Rows with columns (for file details with size)
- **Paragraph** - Text display (for status, help)
- **Block** - Borders and titles

#### Event Handling (with crossterm)
```rust
use crossterm::event::{self, Event, KeyCode, KeyEventKind};

fn handle_events(app: &mut App) -> io::Result<bool> {
    if event::poll(Duration::from_millis(100))? {
        if let Event::Key(key) = event::read()? {
            if key.kind == KeyEventKind::Press {
                match key.code {
                    KeyCode::Char('q') => return Ok(true),
                    KeyCode::Up => app.previous(),
                    KeyCode::Down => app.next(),
                    KeyCode::Enter => app.select(),
                    KeyCode::Backspace => app.go_up(),
                    KeyCode::Char('/') => app.start_search(),
                    _ => {}
                }
            }
        }
    }
    Ok(false)
}
```

#### Useful Extensions
- **ratatui-explorer** - Ready-made file browser widget
- **tui-tree-widget** - Hierarchical tree navigation
- **tui-input** - Text input widget (for search)

---

### 4. Fuzzy Search with Nucleo

Nucleo is 6x faster than fuzzy-matcher/skim and has excellent Unicode support:

```rust
use nucleo_matcher::{Matcher, Config};
use nucleo_matcher::pattern::{Pattern, CaseMatching, Normalization};

fn fuzzy_filter(items: &[String], query: &str) -> Vec<(usize, i64)> {
    let mut matcher = Matcher::new(Config::DEFAULT);
    let pattern = Pattern::parse(query, CaseMatching::Ignore, Normalization::Smart);

    let mut results: Vec<_> = items.iter()
        .enumerate()
        .filter_map(|(idx, item)| {
            pattern.score(item.chars(), &mut matcher)
                .map(|score| (idx, score as i64))
        })
        .collect();

    results.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by score descending
    results
}
```

For highlighting matched characters:
```rust
let indices = pattern.indices(item.chars(), &mut matcher);
// indices contains positions of matched characters for styling
```

---

### 5. Recommended Application Architecture

Based on STU (S3 browser) and The Elm Architecture (TEA):

```
gcs-browser/
├── Cargo.toml
├── src/
│   ├── main.rs           # Entry point, tokio runtime
│   ├── app.rs            # Application state (Model)
│   ├── event.rs          # Event loop, keyboard handling
│   ├── message.rs        # Message enum for state transitions
│   ├── ui/
│   │   ├── mod.rs
│   │   ├── browser.rs    # File/bucket list widget
│   │   ├── search.rs     # Fuzzy search overlay
│   │   ├── status.rs     # Status bar with path, size
│   │   └── help.rs       # Help/keybindings overlay
│   ├── gcs/
│   │   ├── mod.rs
│   │   ├── client.rs     # Async GCS wrapper
│   │   └── types.rs      # Bucket, Object, Folder types
│   └── download.rs       # Download manager with progress
```

#### Core State (app.rs)
```rust
pub struct App {
    // Navigation
    pub current_view: View,
    pub buckets: Vec<Bucket>,
    pub current_bucket: Option<String>,
    pub current_prefix: String,
    pub entries: Vec<Entry>,
    pub selected_index: usize,

    // Search
    pub search_mode: bool,
    pub search_query: String,
    pub filtered_indices: Vec<usize>,

    // Download
    pub download_dir: PathBuf,
    pub downloads: Vec<Download>,

    // Status
    pub loading: bool,
    pub error: Option<String>,
    pub total_size: u64,
}

pub enum View {
    BucketList,
    ObjectBrowser,
    Help,
}

pub enum Entry {
    Folder { name: String, prefix: String },
    Object { name: String, size: u64, modified: DateTime<Utc> },
}
```

#### Message Enum (message.rs)
```rust
pub enum Message {
    // Navigation
    SelectBucket(String),
    NavigateToPrefix(String),
    GoUp,

    // Data loading
    BucketsLoaded(Vec<Bucket>),
    ObjectsLoaded { entries: Vec<Entry>, total_size: u64 },
    LoadError(String),

    // Search
    StartSearch,
    UpdateSearch(String),
    EndSearch,

    // Download
    DownloadSelected,
    DownloadProgress { id: usize, progress: f64 },
    DownloadComplete(usize),

    // App
    Quit,
}
```

#### Main Loop with Async (main.rs)
```rust
#[tokio::main]
async fn main() -> Result<()> {
    let mut terminal = ratatui::init();
    let mut app = App::new(std::env::current_dir()?);

    // Load buckets on startup
    let buckets = gcs::list_buckets().await?;
    app.handle(Message::BucketsLoaded(buckets));

    loop {
        terminal.draw(|f| ui::render(&app, f))?;

        if crossterm::event::poll(Duration::from_millis(100))? {
            if let Event::Key(key) = crossterm::event::read()? {
                if let Some(msg) = handle_key(&app, key) {
                    if matches!(msg, Message::Quit) {
                        break;
                    }
                    // Handle async messages
                    if let Some(future) = app.handle(msg) {
                        tokio::spawn(async move {
                            // Send result back via channel
                        });
                    }
                }
            }
        }
    }

    ratatui::restore();
    Ok(())
}
```

---

### 6. Key Implementation Details

#### Calculating Directory Sizes
```rust
async fn calculate_prefix_size(bucket: &str, prefix: &str) -> u64 {
    let request = ListRequest {
        prefix: Some(prefix.to_string()),
        ..Default::default()
    };

    let mut total = 0u64;
    let objects = Object::list(bucket, request).await?;
    for page in objects {
        for obj in page.items {
            total += obj.size;
        }
    }
    total
}
```

#### Download with Progress
```rust
use tokio::fs::File;
use tokio::io::AsyncWriteExt;

async fn download_file(
    bucket: &str,
    object: &str,
    dest: &Path,
    progress_tx: mpsc::Sender<f64>,
) -> Result<()> {
    let metadata = Object::read(bucket, object).await?;
    let total_size = metadata.size;

    let mut stream = Object::download_streamed(bucket, object).await?;
    let mut file = File::create(dest).await?;
    let mut downloaded = 0u64;

    while let Some(chunk) = stream.next().await {
        let chunk = chunk?;
        downloaded += chunk.len() as u64;
        file.write_all(&chunk).await?;
        progress_tx.send(downloaded as f64 / total_size as f64).await?;
    }

    Ok(())
}
```

#### Fuzzy Search Integration
```rust
impl App {
    pub fn update_search(&mut self, query: String) {
        self.search_query = query;

        if self.search_query.is_empty() {
            self.filtered_indices = (0..self.entries.len()).collect();
            return;
        }

        let names: Vec<_> = self.entries.iter()
            .map(|e| e.name())
            .collect();

        self.filtered_indices = fuzzy_filter(&names, &self.search_query)
            .into_iter()
            .map(|(idx, _)| idx)
            .collect();

        self.selected_index = 0;
    }
}
```

---

### 7. Cargo.toml

```toml
[package]
name = "gcs-browser"
version = "0.1.0"
edition = "2021"

[dependencies]
# TUI
ratatui = "0.28"
crossterm = { version = "0.28", features = ["event-stream"] }
tui-input = "0.10"

# Async
tokio = { version = "1", features = ["full"] }
futures = "0.3"

# GCS
cloud-storage = "0.11"

# Fuzzy search
nucleo-matcher = "0.3"

# Utilities
chrono = { version = "0.4", features = ["serde"] }
humansize = "2"
anyhow = "1"
directories = "5"  # For default download location
```

---

## Reference Projects

| Project | GitHub | Relevance |
|---------|--------|-----------|
| **STU** | [lusingander/stu](https://github.com/lusingander/stu) | S3 browser TUI - closest analog |
| **Yazi** | [sxyazi/yazi](https://github.com/sxyazi/yazi) | Sophisticated async file manager |
| **Joshuto** | [kamiyaa/joshuto](https://github.com/kamiyaa/joshuto) | Ranger-like architecture |
| **Television** | [alexpasmantier/television](https://github.com/alexpasmantier/television) | Nucleo + ratatui fuzzy finder |
| **ratatui-explorer** | [tatounee/ratatui-explorer](https://github.com/tatounee/ratatui-explorer) | Reusable file browser widget |

---

## Proposed Keybindings

| Key | Action |
|-----|--------|
| `j` / `↓` | Move down |
| `k` / `↑` | Move up |
| `Enter` | Select bucket / Enter directory |
| `Backspace` / `h` | Go up one level |
| `/` | Start fuzzy search |
| `Esc` | Cancel search / Go back |
| `d` | Download selected file |
| `D` | Download selected directory (recursive) |
| `Space` | Toggle selection (for multi-download) |
| `?` | Show help |
| `q` | Quit |

---

## Open Questions

1. **Multi-project support**: Should the app support switching between GCP projects, or assume current project from ADC?

2. **File preview**: Should we support previewing text files inline? This adds complexity but is useful.

3. **Upload support**: Is upload functionality needed, or just browsing/downloading?

4. **Caching**: Should we cache bucket/object listings for faster navigation, or always fetch fresh?

5. **Large directory handling**: For directories with 10k+ objects, should we implement virtual scrolling or pagination?

---

## Implementation Phases

### Phase 1: Core Browser
- Bucket listing and selection
- Object listing with prefix navigation
- Basic keyboard navigation
- Status bar with current path

### Phase 2: Search & Details
- Fuzzy search overlay
- File size display
- Directory size calculation
- Human-readable size formatting

### Phase 3: Downloads
- Single file download
- Download progress indicator
- Customizable download directory
- Multi-file selection and batch download

### Phase 4: Polish
- Help screen
- Error handling and retry
- Configuration file support
- Color theming

---

## Sources

- [Google Cloud Storage JSON API v1](https://cloud.google.com/storage/docs/json_api/v1)
- [gcloud storage ls Reference](https://cloud.google.com/sdk/gcloud/reference/storage/ls)
- [Ratatui Official Website](https://ratatui.rs/)
- [cloud-storage Crate](https://docs.rs/cloud-storage)
- [nucleo-matcher Crate](https://docs.rs/nucleo-matcher)
- [STU - S3 TUI Explorer](https://github.com/lusingander/stu)
- [Yazi File Manager](https://github.com/sxyazi/yazi)
- [Television Fuzzy Finder](https://github.com/alexpasmantier/television)
- [Ratatui Templates](https://ratatui.rs/templates/)
- [The Elm Architecture for Ratatui](https://ratatui.rs/concepts/application-patterns/the-elm-architecture/)
