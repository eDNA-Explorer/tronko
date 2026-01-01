# GCS Browser TUI Implementation Plan

## Overview

Build a terminal-based file browser for Google Cloud Storage (GCS) in Rust. The TUI will allow users to:
- List all buckets in their current GCP project
- Navigate directories within buckets
- View file sizes and directory totals
- Search files with fuzzy matching
- Preview text files inline
- Download files to a local directory

This will be a new crate (`gcs-browser`) in the existing Rust workspace at `rust/`.

## Current State Analysis

### Existing Infrastructure
- Rust workspace at `rust/` with three crates: `tronko-core`, `tronko-rs`, `tronko-bench`
- Workspace uses edition 2024, requires Rust 1.75+
- Standard workspace dependency management pattern established

### What We're Building
A standalone TUI application following The Elm Architecture (TEA) pattern with:
- **Model**: Application state (`App` struct)
- **Update**: Message-based state transitions
- **View**: Immediate-mode rendering with ratatui

### Key Design Decisions (Resolved)
| Decision | Choice | Rationale |
|----------|--------|-----------|
| GCP Project | Single (ADC) | Simpler implementation, covers primary use case |
| File Preview | Basic text preview | Show first ~100 lines of text files |
| Upload Support | Download only | Safer, simpler, covers primary use case |
| Crate Name | `gcs-browser` | Generic name, not tronko-specific |

## Desired End State

A working TUI that:
1. Launches and authenticates via Application Default Credentials
2. Lists all buckets in the user's GCP project
3. Allows navigating into buckets and directories
4. Displays files with human-readable sizes
5. Shows directory total sizes
6. Supports fuzzy search with `/` key
7. Shows text file previews in a side panel
8. Downloads files with progress indication
9. Handles errors gracefully with user-friendly messages

### Verification
- `cargo build --release -p gcs-browser` completes successfully
- `./target/release/gcs-browser` launches without errors
- Can browse a real GCS bucket and download files

## What We're NOT Doing

- Multi-project switching
- File upload functionality
- Image/binary file preview
- Configuration file support (hardcoded defaults)
- Color theming customization
- Caching of bucket/object listings
- Virtual scrolling for 10k+ object directories
- Recursive directory download (single files only initially)

---

## Implementation Approach

We'll follow the phased approach from the research document, building incrementally with testable milestones. Each phase delivers working functionality.

**Architecture Pattern**: The Elm Architecture (TEA)
- Unidirectional data flow
- Messages trigger state updates
- Pure rendering from state

**Async Strategy**:
- Tokio runtime for async GCS operations
- Channel-based communication between async tasks and UI
- Non-blocking UI updates during downloads

---

## Phase 1: Project Setup & Core Navigation

### Overview
Set up the crate structure, basic TUI framework, and bucket listing with keyboard navigation.

### Changes Required:

#### 1. Workspace Configuration
**File**: `rust/Cargo.toml`
**Changes**: Add gcs-browser to workspace members and new dependencies

```toml
[workspace]
resolver = "2"
members = [
    "crates/tronko-core",
    "crates/tronko-rs",
    "crates/tronko-bench",
    "crates/gcs-browser",  # Add this
]

[workspace.dependencies]
# ... existing deps ...

# TUI (add these)
ratatui = "0.29"
crossterm = { version = "0.28", features = ["event-stream"] }
tui-input = "0.11"

# GCS
cloud-storage = "0.11"

# Fuzzy search
nucleo-matcher = "0.3"

# Async
tokio = { version = "1", features = ["full"] }
futures = "0.3"

# Utilities
chrono = { version = "0.4", features = ["serde"] }
humansize = "2"
directories = "5"
```

#### 2. Create Crate Structure
**Directory**: `rust/crates/gcs-browser/`

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
│   │   └── status.rs     # Status bar with path, size
│   └── gcs/
│       ├── mod.rs
│       ├── client.rs     # Async GCS wrapper
│       └── types.rs      # Bucket, Object, Folder types
```

#### 3. Crate Cargo.toml
**File**: `rust/crates/gcs-browser/Cargo.toml`

```toml
[package]
name = "gcs-browser"
version.workspace = true
edition.workspace = true
authors.workspace = true
license.workspace = true
description = "Terminal file browser for Google Cloud Storage"

[[bin]]
name = "gcs-browser"
path = "src/main.rs"

[dependencies]
# TUI
ratatui.workspace = true
crossterm.workspace = true

# Async
tokio.workspace = true
futures.workspace = true

# GCS
cloud-storage.workspace = true

# Utilities
anyhow.workspace = true
chrono.workspace = true
humansize.workspace = true
directories.workspace = true
```

#### 4. Core Types
**File**: `rust/crates/gcs-browser/src/gcs/types.rs`

```rust
use chrono::{DateTime, Utc};

#[derive(Debug, Clone)]
pub struct Bucket {
    pub name: String,
    pub location: String,
    pub created: DateTime<Utc>,
}

#[derive(Debug, Clone)]
pub enum Entry {
    Folder {
        name: String,
        prefix: String,
    },
    Object {
        name: String,
        full_path: String,
        size: u64,
        modified: DateTime<Utc>,
    },
}

impl Entry {
    pub fn name(&self) -> &str {
        match self {
            Entry::Folder { name, .. } => name,
            Entry::Object { name, .. } => name,
        }
    }

    pub fn is_folder(&self) -> bool {
        matches!(self, Entry::Folder { .. })
    }
}
```

#### 5. Application State
**File**: `rust/crates/gcs-browser/src/app.rs`

```rust
use std::path::PathBuf;
use crate::gcs::types::{Bucket, Entry};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum View {
    BucketList,
    ObjectBrowser,
}

pub struct App {
    // Navigation
    pub current_view: View,
    pub buckets: Vec<Bucket>,
    pub current_bucket: Option<String>,
    pub current_prefix: String,
    pub entries: Vec<Entry>,
    pub selected_index: usize,

    // Breadcrumb path for display
    pub path_segments: Vec<String>,

    // Download
    pub download_dir: PathBuf,

    // Status
    pub loading: bool,
    pub error: Option<String>,
    pub total_size: u64,

    // For quitting
    pub should_quit: bool,
}

impl App {
    pub fn new(download_dir: PathBuf) -> Self {
        Self {
            current_view: View::BucketList,
            buckets: Vec::new(),
            current_bucket: None,
            current_prefix: String::new(),
            entries: Vec::new(),
            selected_index: 0,
            path_segments: Vec::new(),
            download_dir,
            loading: true,
            error: None,
            total_size: 0,
            should_quit: false,
        }
    }

    pub fn select_next(&mut self) {
        let max = match self.current_view {
            View::BucketList => self.buckets.len(),
            View::ObjectBrowser => self.entries.len(),
        };
        if max > 0 {
            self.selected_index = (self.selected_index + 1).min(max - 1);
        }
    }

    pub fn select_previous(&mut self) {
        self.selected_index = self.selected_index.saturating_sub(1);
    }

    pub fn current_path(&self) -> String {
        match &self.current_bucket {
            Some(bucket) => {
                if self.current_prefix.is_empty() {
                    format!("gs://{}/", bucket)
                } else {
                    format!("gs://{}/{}", bucket, self.current_prefix)
                }
            }
            None => "Buckets".to_string(),
        }
    }
}
```

#### 6. Message Enum
**File**: `rust/crates/gcs-browser/src/message.rs`

```rust
use crate::gcs::types::{Bucket, Entry};

#[derive(Debug)]
pub enum Message {
    // Navigation
    SelectBucket(String),
    NavigateToPrefix(String),
    GoUp,
    Select,

    // Data loading
    BucketsLoaded(Vec<Bucket>),
    ObjectsLoaded { entries: Vec<Entry>, total_size: u64 },
    LoadError(String),

    // App control
    Quit,
}
```

#### 7. GCS Client Wrapper
**File**: `rust/crates/gcs-browser/src/gcs/client.rs`

```rust
use anyhow::Result;
use cloud_storage::{Client, ListRequest, Object};
use crate::gcs::types::{Bucket, Entry};

pub struct GcsClient {
    client: Client,
}

impl GcsClient {
    pub fn new() -> Result<Self> {
        Ok(Self {
            client: Client::default(),
        })
    }

    pub async fn list_buckets(&self) -> Result<Vec<Bucket>> {
        let buckets = self.client.bucket().list().await?;
        Ok(buckets
            .into_iter()
            .map(|b| Bucket {
                name: b.name,
                location: b.location.unwrap_or_default(),
                created: b.time_created,
            })
            .collect())
    }

    pub async fn list_objects(&self, bucket: &str, prefix: &str) -> Result<(Vec<Entry>, u64)> {
        let request = ListRequest {
            prefix: if prefix.is_empty() { None } else { Some(prefix.to_string()) },
            delimiter: Some("/".to_string()),
            max_results: Some(1000),
            ..Default::default()
        };

        let response = Object::list(bucket, request).await?;
        let mut entries = Vec::new();
        let mut total_size = 0u64;

        for page in response {
            // Add folders (prefixes)
            for prefix in page.prefixes {
                let name = prefix
                    .trim_end_matches('/')
                    .rsplit('/')
                    .next()
                    .unwrap_or(&prefix)
                    .to_string();
                entries.push(Entry::Folder {
                    name,
                    prefix: prefix.clone(),
                });
            }

            // Add objects
            for obj in page.items {
                // Skip "directory marker" objects
                if obj.name.ends_with('/') {
                    continue;
                }
                let name = obj.name
                    .rsplit('/')
                    .next()
                    .unwrap_or(&obj.name)
                    .to_string();
                total_size += obj.size;
                entries.push(Entry::Object {
                    name,
                    full_path: obj.name,
                    size: obj.size,
                    modified: obj.updated.unwrap_or(obj.time_created),
                });
            }
        }

        // Sort: folders first, then alphabetically
        entries.sort_by(|a, b| {
            match (a.is_folder(), b.is_folder()) {
                (true, false) => std::cmp::Ordering::Less,
                (false, true) => std::cmp::Ordering::Greater,
                _ => a.name().cmp(b.name()),
            }
        });

        Ok((entries, total_size))
    }
}
```

#### 8. GCS Module
**File**: `rust/crates/gcs-browser/src/gcs/mod.rs`

```rust
pub mod client;
pub mod types;

pub use client::GcsClient;
pub use types::{Bucket, Entry};
```

#### 9. Event Handling
**File**: `rust/crates/gcs-browser/src/event.rs`

```rust
use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use std::time::Duration;
use crate::app::{App, View};
use crate::message::Message;

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

    match key.code {
        KeyCode::Char('q') | KeyCode::Esc => Some(Message::Quit),
        KeyCode::Up | KeyCode::Char('k') => None, // Handled directly
        KeyCode::Down | KeyCode::Char('j') => None, // Handled directly
        KeyCode::Enter | KeyCode::Char('l') => Some(Message::Select),
        KeyCode::Backspace | KeyCode::Char('h') => {
            if app.current_view == View::ObjectBrowser {
                Some(Message::GoUp)
            } else {
                None
            }
        }
        _ => None,
    }
}
```

#### 10. UI Browser Widget
**File**: `rust/crates/gcs-browser/src/ui/browser.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, List, ListItem, ListState},
    Frame,
};
use humansize::{format_size, BINARY};
use crate::app::{App, View};

pub fn render_browser(f: &mut Frame, area: Rect, app: &App) {
    match app.current_view {
        View::BucketList => render_bucket_list(f, area, app),
        View::ObjectBrowser => render_object_list(f, area, app),
    }
}

fn render_bucket_list(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = app
        .buckets
        .iter()
        .map(|b| {
            ListItem::new(Line::from(vec![
                Span::styled("  ", Style::default().fg(Color::Yellow)),
                Span::raw(&b.name),
                Span::styled(
                    format!("  ({})", b.location),
                    Style::default().fg(Color::DarkGray),
                ),
            ]))
        })
        .collect();

    let list = List::new(items)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Buckets "),
        )
        .highlight_style(
            Style::default()
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        )
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.selected_index));

    f.render_stateful_widget(list, area, &mut state);
}

fn render_object_list(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = app
        .entries
        .iter()
        .map(|entry| {
            match entry {
                crate::gcs::Entry::Folder { name, .. } => {
                    ListItem::new(Line::from(vec![
                        Span::styled("  ", Style::default().fg(Color::Blue)),
                        Span::styled(name, Style::default().fg(Color::Blue)),
                        Span::raw("/"),
                    ]))
                }
                crate::gcs::Entry::Object { name, size, .. } => {
                    ListItem::new(Line::from(vec![
                        Span::styled("  ", Style::default().fg(Color::White)),
                        Span::raw(name),
                        Span::styled(
                            format!("  {}", format_size(*size, BINARY)),
                            Style::default().fg(Color::DarkGray),
                        ),
                    ]))
                }
            }
        })
        .collect();

    let title = format!(" {} ", app.current_path());
    let list = List::new(items)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(title),
        )
        .highlight_style(
            Style::default()
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        )
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.selected_index));

    f.render_stateful_widget(list, area, &mut state);
}
```

#### 11. UI Status Bar
**File**: `rust/crates/gcs-browser/src/ui/status.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::Paragraph,
    Frame,
};
use humansize::{format_size, BINARY};
use crate::app::App;

pub fn render_status(f: &mut Frame, area: Rect, app: &App) {
    let status = if app.loading {
        Line::from(vec![
            Span::styled(" Loading...", Style::default().fg(Color::Yellow)),
        ])
    } else if let Some(ref error) = app.error {
        Line::from(vec![
            Span::styled(format!(" Error: {}", error), Style::default().fg(Color::Red)),
        ])
    } else {
        let item_count = match app.current_view {
            crate::app::View::BucketList => app.buckets.len(),
            crate::app::View::ObjectBrowser => app.entries.len(),
        };
        Line::from(vec![
            Span::raw(format!(" {} items", item_count)),
            Span::styled(" | ", Style::default().fg(Color::DarkGray)),
            Span::raw(format!("Total: {}", format_size(app.total_size, BINARY))),
            Span::styled(" | ", Style::default().fg(Color::DarkGray)),
            Span::raw("q:quit  Enter:select  Backspace:back"),
        ])
    };

    let paragraph = Paragraph::new(status)
        .style(Style::default().bg(Color::DarkGray));

    f.render_widget(paragraph, area);
}
```

#### 12. UI Module
**File**: `rust/crates/gcs-browser/src/ui/mod.rs`

```rust
pub mod browser;
pub mod status;

use ratatui::{
    layout::{Constraint, Direction, Layout},
    Frame,
};
use crate::app::App;

pub fn render(f: &mut Frame, app: &App) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(3),      // Browser
            Constraint::Length(1),   // Status bar
        ])
        .split(f.area());

    browser::render_browser(f, chunks[0], app);
    status::render_status(f, chunks[1], app);
}
```

#### 13. Main Entry Point
**File**: `rust/crates/gcs-browser/src/main.rs`

```rust
mod app;
mod event;
mod gcs;
mod message;
mod ui;

use std::time::Duration;
use anyhow::Result;
use crossterm::event::{Event, KeyCode, KeyEventKind};
use ratatui::DefaultTerminal;

use app::{App, View};
use gcs::GcsClient;
use message::Message;

#[tokio::main]
async fn main() -> Result<()> {
    // Initialize terminal
    let mut terminal = ratatui::init();
    terminal.clear()?;

    // Initialize app with current directory as download location
    let download_dir = std::env::current_dir()?;
    let mut app = App::new(download_dir);

    // Load buckets on startup
    let result = run_app(&mut terminal, &mut app).await;

    // Restore terminal
    ratatui::restore();

    result
}

async fn run_app(terminal: &mut DefaultTerminal, app: &mut App) -> Result<()> {
    let client = GcsClient::new()?;

    // Initial bucket load
    match client.list_buckets().await {
        Ok(buckets) => {
            app.buckets = buckets;
            app.loading = false;
        }
        Err(e) => {
            app.error = Some(e.to_string());
            app.loading = false;
        }
    }

    loop {
        terminal.draw(|f| ui::render(f, app))?;

        if let Some(event) = event::poll_event(Duration::from_millis(100))? {
            if let Event::Key(key) = event {
                if key.kind == KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char('q') | KeyCode::Esc => break,
                        KeyCode::Up | KeyCode::Char('k') => app.select_previous(),
                        KeyCode::Down | KeyCode::Char('j') => app.select_next(),
                        KeyCode::Enter | KeyCode::Char('l') => {
                            handle_select(app, &client).await?;
                        }
                        KeyCode::Backspace | KeyCode::Char('h') => {
                            handle_go_up(app, &client).await?;
                        }
                        _ => {}
                    }
                }
            }
        }

        if app.should_quit {
            break;
        }
    }

    Ok(())
}

async fn handle_select(app: &mut App, client: &GcsClient) -> Result<()> {
    match app.current_view {
        View::BucketList => {
            if let Some(bucket) = app.buckets.get(app.selected_index) {
                let bucket_name = bucket.name.clone();
                app.current_bucket = Some(bucket_name.clone());
                app.current_prefix = String::new();
                app.path_segments = vec![bucket_name.clone()];
                app.current_view = View::ObjectBrowser;
                app.loading = true;
                app.selected_index = 0;

                match client.list_objects(&bucket_name, "").await {
                    Ok((entries, total_size)) => {
                        app.entries = entries;
                        app.total_size = total_size;
                        app.loading = false;
                    }
                    Err(e) => {
                        app.error = Some(e.to_string());
                        app.loading = false;
                    }
                }
            }
        }
        View::ObjectBrowser => {
            if let Some(entry) = app.entries.get(app.selected_index) {
                if let gcs::Entry::Folder { prefix, name } = entry {
                    let bucket = app.current_bucket.as_ref().unwrap();
                    app.current_prefix = prefix.clone();
                    app.path_segments.push(name.clone());
                    app.loading = true;
                    app.selected_index = 0;

                    match client.list_objects(bucket, prefix).await {
                        Ok((entries, total_size)) => {
                            app.entries = entries;
                            app.total_size = total_size;
                            app.loading = false;
                        }
                        Err(e) => {
                            app.error = Some(e.to_string());
                            app.loading = false;
                        }
                    }
                }
                // TODO: Handle file selection (preview or download)
            }
        }
    }
    Ok(())
}

async fn handle_go_up(app: &mut App, client: &GcsClient) -> Result<()> {
    if app.current_view == View::ObjectBrowser {
        if app.current_prefix.is_empty() {
            // Go back to bucket list
            app.current_view = View::BucketList;
            app.current_bucket = None;
            app.path_segments.clear();
            app.entries.clear();
            app.selected_index = 0;
            app.total_size = 0;
        } else {
            // Go up one directory level
            let parts: Vec<&str> = app.current_prefix.trim_end_matches('/').split('/').collect();
            if parts.len() > 1 {
                app.current_prefix = parts[..parts.len() - 1].join("/") + "/";
            } else {
                app.current_prefix = String::new();
            }
            app.path_segments.pop();
            app.selected_index = 0;

            let bucket = app.current_bucket.as_ref().unwrap();
            app.loading = true;

            match client.list_objects(bucket, &app.current_prefix).await {
                Ok((entries, total_size)) => {
                    app.entries = entries;
                    app.total_size = total_size;
                    app.loading = false;
                }
                Err(e) => {
                    app.error = Some(e.to_string());
                    app.loading = false;
                }
            }
        }
    }
    Ok(())
}
```

### Success Criteria:

#### Automated Verification:
- [x] `cd rust && cargo build -p gcs-browser` compiles without errors
- [x] `cd rust && cargo clippy -p gcs-browser` has no warnings (dead code warnings expected for future phase fields)
- [x] `cd rust && cargo test -p gcs-browser` passes (no tests yet, but should not error)

#### Manual Verification:
- [ ] `./target/release/gcs-browser` launches and shows bucket list
- [ ] Can navigate up/down with j/k or arrow keys
- [ ] Can press Enter to enter a bucket
- [ ] Can navigate directories within the bucket
- [ ] Backspace returns to parent directory/bucket list
- [ ] 'q' quits the application cleanly
- [ ] Status bar shows item count and keybinding hints

---

## Phase 2: Fuzzy Search

### Overview
Add fuzzy search functionality to filter files/folders using the `/` key and nucleo-matcher.

### Changes Required:

#### 1. Add Dependencies
**File**: `rust/crates/gcs-browser/Cargo.toml`
**Changes**: Add nucleo-matcher and tui-input

```toml
[dependencies]
# ... existing deps ...
nucleo-matcher.workspace = true
tui-input.workspace = true
```

#### 2. Extend App State
**File**: `rust/crates/gcs-browser/src/app.rs`
**Changes**: Add search-related fields

```rust
// Add to App struct:
pub search_mode: bool,
pub search_query: String,
pub filtered_indices: Vec<usize>,

// Add to App::new():
search_mode: false,
search_query: String::new(),
filtered_indices: Vec::new(),
```

#### 3. Create Search Module
**File**: `rust/crates/gcs-browser/src/search.rs`

```rust
use nucleo_matcher::{Config, Matcher};
use nucleo_matcher::pattern::{CaseMatching, Normalization, Pattern};

pub fn fuzzy_filter(items: &[String], query: &str) -> Vec<usize> {
    if query.is_empty() {
        return (0..items.len()).collect();
    }

    let mut matcher = Matcher::new(Config::DEFAULT);
    let pattern = Pattern::parse(query, CaseMatching::Ignore, Normalization::Smart);

    let mut results: Vec<(usize, u32)> = items
        .iter()
        .enumerate()
        .filter_map(|(idx, item)| {
            let mut buf = Vec::new();
            pattern
                .score(nucleo_matcher::Utf32Str::new(item, &mut buf), &mut matcher)
                .map(|score| (idx, score))
        })
        .collect();

    results.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by score descending
    results.into_iter().map(|(idx, _)| idx).collect()
}
```

#### 4. Add Search UI Component
**File**: `rust/crates/gcs-browser/src/ui/search.rs`

```rust
use ratatui::{
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Style},
    widgets::{Block, Borders, Clear, Paragraph},
    Frame,
};
use tui_input::Input;

pub fn render_search_overlay(f: &mut Frame, query: &str) {
    let area = centered_rect(60, 3, f.area());

    f.render_widget(Clear, area);

    let input = Paragraph::new(format!("/{}", query))
        .style(Style::default().fg(Color::White))
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Search ")
                .border_style(Style::default().fg(Color::Cyan)),
        );

    f.render_widget(input, area);
}

fn centered_rect(percent_x: u16, height: u16, r: Rect) -> Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length((r.height - height) / 2),
            Constraint::Length(height),
            Constraint::Min(0),
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

#### 5. Update Main Loop
**File**: `rust/crates/gcs-browser/src/main.rs`
**Changes**: Handle '/' to enter search mode, Esc to exit, character input

```rust
// In the key event handling:
KeyCode::Char('/') if !app.search_mode => {
    app.search_mode = true;
    app.search_query.clear();
}
KeyCode::Esc if app.search_mode => {
    app.search_mode = false;
    app.search_query.clear();
    app.filtered_indices = (0..app.entries.len()).collect();
}
KeyCode::Char(c) if app.search_mode => {
    app.search_query.push(c);
    update_search_filter(app);
}
KeyCode::Backspace if app.search_mode && !app.search_query.is_empty() => {
    app.search_query.pop();
    update_search_filter(app);
}
```

#### 6. Update UI Rendering
**File**: `rust/crates/gcs-browser/src/ui/mod.rs`
**Changes**: Render search overlay when in search mode

```rust
pub fn render(f: &mut Frame, app: &App) {
    // ... existing rendering ...

    if app.search_mode {
        search::render_search_overlay(f, &app.search_query);
    }
}
```

#### 7. Update Browser to Use Filtered List
**File**: `rust/crates/gcs-browser/src/ui/browser.rs`
**Changes**: When filtered_indices is non-empty, only show those items

### Success Criteria:

#### Automated Verification:
- [x] `cargo build -p gcs-browser` compiles
- [x] `cargo test -p gcs-browser` passes
- [x] Search module unit tests pass (6 tests)

#### Manual Verification:
- [ ] Press `/` opens search input overlay
- [ ] Typing filters the list in real-time
- [ ] Fuzzy matching works (e.g., "abc" matches "a-big-cat.txt")
- [ ] Esc closes search and shows all items
- [ ] Enter with search active selects the highlighted filtered item

---

## Phase 3: File Preview

### Overview
Add a preview panel on the right side that shows the first ~100 lines of text files when selected.

### Changes Required:

#### 1. Extend App State
**File**: `rust/crates/gcs-browser/src/app.rs`

```rust
// Add to App struct:
pub preview_content: Option<String>,
pub preview_loading: bool,
pub preview_error: Option<String>,
```

#### 2. Add Preview Fetching
**File**: `rust/crates/gcs-browser/src/gcs/client.rs`

```rust
impl GcsClient {
    pub async fn preview_object(&self, bucket: &str, object: &str, max_bytes: usize) -> Result<String> {
        // Fetch first max_bytes of the object
        let bytes = Object::download_range(bucket, object, 0..max_bytes as u64).await?;

        // Try to convert to UTF-8, replacing invalid sequences
        let content = String::from_utf8_lossy(&bytes).to_string();

        // Take first ~100 lines
        let lines: Vec<&str> = content.lines().take(100).collect();
        Ok(lines.join("\n"))
    }
}
```

#### 3. Create Preview UI Component
**File**: `rust/crates/gcs-browser/src/ui/preview.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Style},
    widgets::{Block, Borders, Paragraph, Wrap},
    Frame,
};
use crate::app::App;

pub fn render_preview(f: &mut Frame, area: Rect, app: &App) {
    let (content, style) = if app.preview_loading {
        ("Loading...".to_string(), Style::default().fg(Color::Yellow))
    } else if let Some(ref error) = app.preview_error {
        (format!("Error: {}", error), Style::default().fg(Color::Red))
    } else if let Some(ref content) = app.preview_content {
        (content.clone(), Style::default().fg(Color::White))
    } else {
        ("Select a file to preview".to_string(), Style::default().fg(Color::DarkGray))
    };

    let paragraph = Paragraph::new(content)
        .style(style)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Preview "),
        )
        .wrap(Wrap { trim: false });

    f.render_widget(paragraph, area);
}
```

#### 4. Update Layout
**File**: `rust/crates/gcs-browser/src/ui/mod.rs`
**Changes**: Split into two columns when in ObjectBrowser view

```rust
pub fn render(f: &mut Frame, app: &App) {
    let main_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(3),
            Constraint::Length(1),
        ])
        .split(f.area());

    if app.current_view == View::ObjectBrowser {
        // Two-column layout for object browser
        let columns = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([
                Constraint::Percentage(50),
                Constraint::Percentage(50),
            ])
            .split(main_chunks[0]);

        browser::render_browser(f, columns[0], app);
        preview::render_preview(f, columns[1], app);
    } else {
        browser::render_browser(f, main_chunks[0], app);
    }

    status::render_status(f, main_chunks[1], app);
}
```

#### 5. Trigger Preview on Selection Change
**File**: `rust/crates/gcs-browser/src/main.rs`
**Changes**: When selection changes to a file, fetch preview asynchronously

### Success Criteria:

#### Automated Verification:
- [x] `cargo build -p gcs-browser` compiles
- [x] `cargo test -p gcs-browser` passes

#### Manual Verification:
- [ ] Preview pane appears when browsing objects (not buckets)
- [ ] Selecting a text file shows its contents in preview
- [ ] Preview shows "Loading..." while fetching
- [ ] Large files only show first ~100 lines
- [ ] Binary files show error or placeholder message
- [ ] Preview updates as you navigate with j/k

---

## Phase 4: Downloads with Progress

### Overview
Implement file downloading with a progress bar when pressing 'd' on a file.

### Changes Required:

#### 1. Extend App State
**File**: `rust/crates/gcs-browser/src/app.rs`

```rust
// Add to App struct:
pub downloads: Vec<Download>,
pub active_download: Option<ActiveDownload>,

#[derive(Debug, Clone)]
pub struct Download {
    pub filename: String,
    pub status: DownloadStatus,
}

#[derive(Debug, Clone)]
pub enum DownloadStatus {
    InProgress { progress: f64 },
    Complete,
    Failed(String),
}

pub struct ActiveDownload {
    pub filename: String,
    pub progress: f64,
}
```

#### 2. Add Download Function
**File**: `rust/crates/gcs-browser/src/gcs/client.rs`

```rust
use tokio::sync::mpsc;
use tokio::fs::File;
use tokio::io::AsyncWriteExt;

impl GcsClient {
    pub async fn download_object(
        &self,
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
            let _ = progress_tx.send(downloaded as f64 / total_size as f64).await;
        }

        Ok(())
    }
}
```

#### 3. Create Download Progress UI
**File**: `rust/crates/gcs-browser/src/ui/download.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Style},
    widgets::{Block, Borders, Gauge},
    Frame,
};

pub fn render_download_progress(f: &mut Frame, area: Rect, filename: &str, progress: f64) {
    let gauge = Gauge::default()
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(format!(" Downloading: {} ", filename)),
        )
        .gauge_style(Style::default().fg(Color::Cyan))
        .percent((progress * 100.0) as u16)
        .label(format!("{:.1}%", progress * 100.0));

    f.render_widget(gauge, area);
}
```

#### 4. Handle 'd' Key for Download
**File**: `rust/crates/gcs-browser/src/main.rs`
**Changes**: Add key handler and async download logic

```rust
KeyCode::Char('d') => {
    if let Some(entry) = app.entries.get(app.selected_index) {
        if let gcs::Entry::Object { full_path, name, .. } = entry {
            let dest = app.download_dir.join(name);
            // Start async download with progress channel
            // ...
        }
    }
}
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p gcs-browser` compiles
- [ ] `cargo test -p gcs-browser` passes

#### Manual Verification:
- [ ] Press 'd' on a file starts download
- [ ] Progress bar appears showing download progress
- [ ] File is saved to current working directory
- [ ] Large files show incremental progress
- [ ] Download completion shows success message
- [ ] Can continue browsing during download (non-blocking)

---

## Phase 5: Help Screen & Polish

### Overview
Add help screen ('?'), improve error handling, and polish the UI.

### Changes Required:

#### 1. Add Help View
**File**: `rust/crates/gcs-browser/src/ui/help.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
    Frame,
};

pub fn render_help(f: &mut Frame) {
    let area = centered_rect(60, 20, f.area());
    f.render_widget(Clear, area);

    let help_text = vec![
        Line::from(vec![
            Span::styled("GCS Browser", Style::default().add_modifier(Modifier::BOLD)),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("j/↓    ", Style::default().fg(Color::Cyan)),
            Span::raw("Move down"),
        ]),
        Line::from(vec![
            Span::styled("k/↑    ", Style::default().fg(Color::Cyan)),
            Span::raw("Move up"),
        ]),
        Line::from(vec![
            Span::styled("Enter  ", Style::default().fg(Color::Cyan)),
            Span::raw("Select bucket / Enter directory"),
        ]),
        Line::from(vec![
            Span::styled("h/←    ", Style::default().fg(Color::Cyan)),
            Span::raw("Go back / Up one level"),
        ]),
        Line::from(vec![
            Span::styled("/      ", Style::default().fg(Color::Cyan)),
            Span::raw("Fuzzy search"),
        ]),
        Line::from(vec![
            Span::styled("d      ", Style::default().fg(Color::Cyan)),
            Span::raw("Download selected file"),
        ]),
        Line::from(vec![
            Span::styled("?      ", Style::default().fg(Color::Cyan)),
            Span::raw("Show this help"),
        ]),
        Line::from(vec![
            Span::styled("q/Esc  ", Style::default().fg(Color::Cyan)),
            Span::raw("Quit"),
        ]),
        Line::from(""),
        Line::from(Span::styled(
            "Press any key to close",
            Style::default().fg(Color::DarkGray),
        )),
    ];

    let paragraph = Paragraph::new(help_text)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Help ")
                .border_style(Style::default().fg(Color::Cyan)),
        );

    f.render_widget(paragraph, area);
}
```

#### 2. Add Error Recovery
- Retry failed GCS operations with exponential backoff
- Clear error state on user action
- Show actionable error messages

#### 3. Polish Status Bar
- Show download directory location
- Better keybinding hints based on context
- Show "Authenticated as: PROJECT_ID" on startup

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build --release -p gcs-browser` compiles
- [ ] `cargo clippy -p gcs-browser` has no warnings
- [ ] `cargo test -p gcs-browser` passes

#### Manual Verification:
- [ ] Press '?' shows help overlay
- [ ] Any key closes help overlay
- [ ] Errors are displayed clearly and don't crash the app
- [ ] Can recover from network errors by retrying
- [ ] Application feels responsive during all operations

---

## Testing Strategy

### Unit Tests:
- `search.rs`: Test fuzzy matching edge cases
- `gcs/types.rs`: Test Entry sorting and name extraction
- `app.rs`: Test navigation state transitions

### Integration Tests:
- Mock GCS client for offline testing
- Test full navigation flow without real GCS

### Manual Testing Steps:
1. Run `gcloud auth application-default login` to authenticate
2. Launch `gcs-browser` and verify buckets load
3. Navigate into a bucket with various file types
4. Test fuzzy search with partial matches
5. Download a file and verify contents
6. Test error handling by disconnecting network mid-browse

---

## Performance Considerations

- **Lazy loading**: Only fetch directory contents when navigating into them
- **Streaming downloads**: Use chunked streaming to avoid loading entire files into memory
- **Debounced search**: Consider debouncing search input if performance issues arise
- **Pagination**: For very large directories (10k+ objects), implement pagination rather than loading all at once

---

## References

- Research document: `thoughts/shared/research/2025-12-31-rust-gcs-tui-file-browser.md`
- Existing workspace: `rust/Cargo.toml`
- STU (S3 browser): https://github.com/lusingander/stu
- Ratatui documentation: https://ratatui.rs/
- cloud-storage crate: https://docs.rs/cloud-storage
