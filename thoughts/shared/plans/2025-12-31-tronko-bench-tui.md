# tronko-bench TUI Implementation Plan

## Overview

Transform the existing `tronko-bench` CLI into an interactive Terminal User Interface (TUI). The TUI will allow users to:
- Browse discovered reference databases and query files interactively
- Navigate between different views (databases, queries, MSA files)
- Filter entries with fuzzy search
- View file details in a preview panel
- Trigger database conversions with real-time progress display
- Generate benchmark configurations interactively

This will be a refactor of the existing `tronko-bench` crate, adding a new `--tui` flag or `tui` subcommand to launch the interactive mode.

## Current State Analysis

### Existing Infrastructure
- Rust workspace at `rust/` already has TUI dependencies configured:
  - `ratatui = "0.29"` - TUI framework
  - `crossterm = "0.28"` with `event-stream` feature - terminal backend
  - `tui-input = "0.11"` - text input widget
  - `nucleo-matcher = "0.3"` - fuzzy search
- `tronko-bench` crate already has:
  - Well-defined data structures: `DatabaseInfo`, `QueryFileInfo`, `ScanResults`
  - Scanning logic in `scan_databases()` and `scan_query_files()`
  - Conversion logic with progress tracking in `convert_database()`
  - File format detection for `.txt` and `.trkb` databases
  - Human-readable size formatting via `humansize`

### What We're Building
An interactive TUI application following The Elm Architecture (TEA) pattern with:
- **Model**: Application state (`App` struct)
- **Update**: Message-based state transitions
- **View**: Immediate-mode rendering with ratatui

### Key Design Decisions (Resolved)
| Decision | Choice | Rationale |
|----------|--------|-----------|
| Entry point | `--tui` flag on existing CLI | Preserves CLI compatibility, single binary |
| Initial scan | Auto-scan on startup | Matches current `scan` command behavior |
| Async operations | Sync with progress callback | Conversion is blocking, no network I/O |
| Navigation style | Vim-like (hjkl) + arrows | Consistent with gcs-browser |
| Preview panel | File metadata only | No file content preview needed |

## Desired End State

A working TUI that:
1. Launches with `tronko-bench tui` or `tronko-bench --tui`
2. Auto-scans the current directory (or specified `-d` directory)
3. Displays three tabs: Databases, Queries, MSA Files
4. Allows navigating up/down with j/k or arrow keys
5. Supports fuzzy search with `/` key
6. Shows file details in a right-side panel
7. Pressing `c` on an unconverted database starts conversion
8. Pressing `g` generates a benchmark.json config file
9. Pressing `r` rescans the directory
10. Handles errors gracefully with status bar messages

### Verification
- `cargo build --release -p tronko-bench` completes successfully
- `./target/release/tronko-bench tui` launches the TUI
- All existing CLI commands continue to work unchanged
- Can browse files, filter with search, and convert databases

## What We're NOT Doing

- Real-time file watching (manual rescan only)
- Multi-directory scanning (single root at a time)
- Editing benchmark.json configs (read-only generation)
- File content preview (just metadata)
- Mouse support (keyboard only)
- Configuration persistence (runtime only)
- Custom themes or color schemes

---

## Implementation Approach

We'll refactor incrementally, keeping the existing CLI functional while adding TUI capabilities. Each phase delivers working functionality.

**Architecture Pattern**: The Elm Architecture (TEA)
- Unidirectional data flow
- Messages trigger state updates
- Pure rendering from state

**Execution Strategy**:
- Synchronous main loop (no tokio needed - no network I/O)
- Conversion runs in foreground with progress updates
- Non-blocking terminal event polling

---

## Phase 1: Basic TUI Shell & Tab Navigation

### Overview
Add TUI entry point, basic layout with tabs, and keyboard navigation framework.

### Changes Required:

#### 1. Add TUI Dependencies to Cargo.toml
**File**: `rust/crates/tronko-bench/Cargo.toml`
**Changes**: Add TUI-related dependencies from workspace

```toml
[dependencies]
# ... existing deps ...

# TUI
ratatui.workspace = true
crossterm.workspace = true
tui-input.workspace = true
nucleo-matcher.workspace = true
```

#### 2. Create TUI Module Structure
**Directory**: `rust/crates/tronko-bench/src/tui/`

```
tui/
├── mod.rs           # Module exports and TUI entry point
├── app.rs           # Application state (Model)
├── event.rs         # Event loop, keyboard handling
├── message.rs       # Message enum for state transitions
└── ui/
    ├── mod.rs       # Main render function
    ├── tabs.rs      # Tab bar widget
    ├── list.rs      # Database/query list widget
    ├── detail.rs    # Detail panel widget
    └── status.rs    # Status bar widget
```

#### 3. Add TUI Subcommand
**File**: `rust/crates/tronko-bench/src/main.rs`
**Changes**: Add `Tui` variant to `Commands` enum

```rust
#[derive(Subcommand)]
enum Commands {
    // ... existing commands ...

    /// Launch interactive TUI mode
    Tui {
        /// Directory to scan (defaults to current directory)
        #[arg(short, long, default_value = ".")]
        dir: PathBuf,

        /// Maximum depth to scan
        #[arg(short = 'D', long, default_value = "5")]
        max_depth: usize,
    },
}
```

#### 4. Application State
**File**: `rust/crates/tronko-bench/src/tui/app.rs`

```rust
use std::path::PathBuf;
use crate::{DatabaseInfo, QueryFileInfo, ScanResults};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tab {
    Databases,
    Queries,
    MsaFiles,
}

impl Tab {
    pub fn all() -> &'static [Tab] {
        &[Tab::Databases, Tab::Queries, Tab::MsaFiles]
    }

    pub fn title(&self) -> &'static str {
        match self {
            Tab::Databases => "Databases",
            Tab::Queries => "Queries",
            Tab::MsaFiles => "MSA Files",
        }
    }
}

pub struct App {
    // Data
    pub scan_results: ScanResults,
    pub scan_dir: PathBuf,
    pub max_depth: usize,

    // Navigation
    pub current_tab: Tab,
    pub selected_index: usize,

    // Search
    pub search_mode: bool,
    pub search_query: String,
    pub filtered_indices: Vec<usize>,

    // Status
    pub loading: bool,
    pub status_message: Option<String>,
    pub should_quit: bool,
}

impl App {
    pub fn new(scan_dir: PathBuf, max_depth: usize) -> Self {
        Self {
            scan_results: ScanResults {
                databases: Vec::new(),
                query_files: Vec::new(),
                msa_files: Vec::new(),
            },
            scan_dir,
            max_depth,
            current_tab: Tab::Databases,
            selected_index: 0,
            search_mode: false,
            search_query: String::new(),
            filtered_indices: Vec::new(),
            loading: true,
            status_message: Some("Scanning...".to_string()),
            should_quit: false,
        }
    }

    pub fn current_list_len(&self) -> usize {
        if !self.filtered_indices.is_empty() {
            return self.filtered_indices.len();
        }
        match self.current_tab {
            Tab::Databases => self.scan_results.databases.len(),
            Tab::Queries => self.scan_results.query_files.len(),
            Tab::MsaFiles => self.scan_results.msa_files.len(),
        }
    }

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

#### 5. Message Enum
**File**: `rust/crates/tronko-bench/src/tui/message.rs`

```rust
use crate::ScanResults;

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

    // Data
    ScanComplete(ScanResults),
    ScanError(String),

    // Status
    SetStatus(String),
    ClearStatus,

    // App
    Quit,
}
```

#### 6. Event Handling
**File**: `rust/crates/tronko-bench/src/tui/event.rs`

```rust
use crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use std::time::Duration;
use crate::tui::app::App;
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

    // Search mode handling
    if app.search_mode {
        return match key.code {
            KeyCode::Esc => Some(Message::ExitSearch),
            KeyCode::Enter => Some(Message::ExitSearch),
            KeyCode::Backspace => Some(Message::SearchBackspace),
            KeyCode::Char(c) => Some(Message::SearchInput(c)),
            _ => None,
        };
    }

    // Normal mode
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

#### 7. Tab Bar Widget
**File**: `rust/crates/tronko-bench/src/tui/ui/tabs.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Tabs as RataTabs},
    Frame,
};
use crate::tui::app::{App, Tab};

pub fn render_tabs(f: &mut Frame, area: Rect, app: &App) {
    let titles: Vec<Line> = Tab::all()
        .iter()
        .map(|t| {
            let count = match t {
                Tab::Databases => app.scan_results.databases.len(),
                Tab::Queries => app.scan_results.query_files.len(),
                Tab::MsaFiles => app.scan_results.msa_files.len(),
            };
            Line::from(format!("{} ({})", t.title(), count))
        })
        .collect();

    let selected = Tab::all()
        .iter()
        .position(|t| *t == app.current_tab)
        .unwrap_or(0);

    let tabs = RataTabs::new(titles)
        .block(Block::default().borders(Borders::BOTTOM))
        .select(selected)
        .style(Style::default().fg(Color::DarkGray))
        .highlight_style(
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        );

    f.render_widget(tabs, area);
}
```

#### 8. List Widget
**File**: `rust/crates/tronko-bench/src/tui/ui/list.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, List, ListItem, ListState},
    Frame,
};
use crate::tui::app::{App, Tab};

pub fn render_list(f: &mut Frame, area: Rect, app: &App) {
    let items: Vec<ListItem> = match app.current_tab {
        Tab::Databases => render_database_items(app),
        Tab::Queries => render_query_items(app, false),
        Tab::MsaFiles => render_query_items(app, true),
    };

    let title = match app.current_tab {
        Tab::Databases => " Reference Databases ",
        Tab::Queries => " Query Files ",
        Tab::MsaFiles => " MSA Files ",
    };

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
            let format_color = if db.format.contains("binary") {
                Color::Green
            } else {
                Color::Yellow
            };

            ListItem::new(Line::from(vec![
                Span::styled(
                    format!("{:<30}", truncate(&db.name, 30)),
                    Style::default(),
                ),
                Span::styled(
                    format!(" {:>6} trees ", db.num_trees),
                    Style::default().fg(Color::Cyan),
                ),
                Span::styled(
                    format!(" {:>10} ", db.size),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::styled(
                    if db.format.contains("binary") { "[BIN]" } else { "[TXT]" },
                    Style::default().fg(format_color),
                ),
            ]))
        })
        .collect()
}

fn render_query_items(app: &App, msa: bool) -> Vec<ListItem<'static>> {
    let files = if msa {
        &app.scan_results.msa_files
    } else {
        &app.scan_results.query_files
    };

    let indices: Vec<usize> = if app.filtered_indices.is_empty() {
        (0..files.len()).collect()
    } else {
        app.filtered_indices.clone()
    };

    indices
        .into_iter()
        .filter_map(|i| files.get(i))
        .map(|qf| {
            let type_color = if qf.file_type == "FASTA" {
                Color::Blue
            } else {
                Color::Magenta
            };

            ListItem::new(Line::from(vec![
                Span::styled(
                    format!("{:<35}", truncate(&qf.name, 35)),
                    Style::default(),
                ),
                Span::styled(
                    format!(" {:>8} seqs ", qf.num_sequences),
                    Style::default().fg(Color::Cyan),
                ),
                Span::styled(
                    format!(" {:>10} ", qf.size),
                    Style::default().fg(Color::DarkGray),
                ),
                Span::styled(
                    format!("[{}]", qf.file_type),
                    Style::default().fg(type_color),
                ),
            ]))
        })
        .collect()
}

fn truncate(s: &str, max: usize) -> String {
    if s.len() > max {
        format!("{}...", &s[..max - 3])
    } else {
        s.to_string()
    }
}
```

#### 9. Status Bar Widget
**File**: `rust/crates/tronko-bench/src/tui/ui/status.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Style},
    text::{Line, Span},
    widgets::Paragraph,
    Frame,
};
use crate::tui::app::App;

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
        Line::from(vec![
            Span::raw(format!(" {} ", app.scan_dir.display())),
            Span::styled(" | ", Style::default().fg(Color::DarkGray)),
            Span::raw("Tab:switch  j/k:nav  /:search  r:rescan  c:convert  g:config  q:quit"),
        ])
    };

    let paragraph = Paragraph::new(status)
        .style(Style::default().bg(Color::DarkGray));

    f.render_widget(paragraph, area);
}
```

#### 10. Main UI Render Function
**File**: `rust/crates/tronko-bench/src/tui/ui/mod.rs`

```rust
pub mod tabs;
pub mod list;
pub mod detail;
pub mod status;

use ratatui::{
    layout::{Constraint, Direction, Layout},
    Frame,
};
use crate::tui::app::App;

pub fn render(f: &mut Frame, app: &App) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(2),    // Tabs
            Constraint::Min(3),       // Main content
            Constraint::Length(1),    // Status bar
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

    // Render search overlay if active
    if app.search_mode {
        render_search_overlay(f, &app.search_query);
    }
}

fn render_search_overlay(f: &mut Frame, query: &str) {
    use ratatui::{
        layout::{Alignment, Constraint, Direction, Layout, Rect},
        style::{Color, Style},
        widgets::{Block, Borders, Clear, Paragraph},
    };

    let area = centered_rect(50, 3, f.area());
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

fn centered_rect(percent_x: u16, height: u16, r: ratatui::layout::Rect) -> ratatui::layout::Rect {
    let popup_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length((r.height.saturating_sub(height)) / 2),
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

#### 11. Detail Panel Widget
**File**: `rust/crates/tronko-bench/src/tui/ui/detail.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph, Wrap},
    Frame,
};
use crate::tui::app::{App, Tab};

pub fn render_detail(f: &mut Frame, area: Rect, app: &App) {
    let content = match app.current_tab {
        Tab::Databases => render_database_detail(app),
        Tab::Queries | Tab::MsaFiles => render_query_detail(app),
    };

    let paragraph = Paragraph::new(content)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Details "),
        )
        .wrap(Wrap { trim: false });

    f.render_widget(paragraph, area);
}

fn render_database_detail(app: &App) -> Vec<Line<'static>> {
    let db = match get_selected_database(app) {
        Some(db) => db,
        None => return vec![Line::from("No database selected")],
    };

    vec![
        Line::from(vec![
            Span::styled("Name: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(db.name.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Format: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(db.format.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Trees: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(db.num_trees.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Size: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(db.size.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Path: ", Style::default().add_modifier(Modifier::BOLD)),
        ]),
        Line::from(db.path.display().to_string()),
        Line::from(""),
        if db.format == "text" {
            Line::from(vec![
                Span::styled(
                    "Press 'c' to convert to binary",
                    Style::default().fg(Color::Yellow),
                ),
            ])
        } else {
            Line::from("")
        },
    ]
}

fn render_query_detail(app: &App) -> Vec<Line<'static>> {
    let qf = match get_selected_query(app) {
        Some(qf) => qf,
        None => return vec![Line::from("No file selected")],
    };

    vec![
        Line::from(vec![
            Span::styled("Name: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(qf.name.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Type: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(qf.file_type.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Sequences: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(qf.num_sequences.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Size: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(qf.size.clone()),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Path: ", Style::default().add_modifier(Modifier::BOLD)),
        ]),
        Line::from(qf.path.display().to_string()),
    ]
}

fn get_selected_database(app: &App) -> Option<&crate::DatabaseInfo> {
    let idx = if app.filtered_indices.is_empty() {
        app.selected_index
    } else {
        *app.filtered_indices.get(app.selected_index)?
    };
    app.scan_results.databases.get(idx)
}

fn get_selected_query(app: &App) -> Option<&crate::QueryFileInfo> {
    let files = match app.current_tab {
        Tab::Queries => &app.scan_results.query_files,
        Tab::MsaFiles => &app.scan_results.msa_files,
        _ => return None,
    };

    let idx = if app.filtered_indices.is_empty() {
        app.selected_index
    } else {
        *app.filtered_indices.get(app.selected_index)?
    };
    files.get(idx)
}
```

#### 12. TUI Module Entry Point
**File**: `rust/crates/tronko-bench/src/tui/mod.rs`

```rust
pub mod app;
pub mod event;
pub mod message;
pub mod ui;

use std::path::PathBuf;
use std::time::Duration;
use anyhow::Result;
use crossterm::event::Event;
use ratatui::DefaultTerminal;

use app::App;
use message::Message;
use crate::{scan_databases, scan_query_files, ScanResults};

pub fn run_tui(dir: PathBuf, max_depth: usize) -> Result<()> {
    let mut terminal = ratatui::init();
    terminal.clear()?;

    let mut app = App::new(dir, max_depth);

    // Initial scan
    let results = perform_scan(&app.scan_dir, app.max_depth);
    app.scan_results = results;
    app.loading = false;
    app.status_message = None;

    let result = run_app(&mut terminal, &mut app);

    ratatui::restore();
    result
}

fn run_app(terminal: &mut DefaultTerminal, app: &mut App) -> Result<()> {
    loop {
        terminal.draw(|f| ui::render(f, app))?;

        if let Some(event) = event::poll_event(Duration::from_millis(100))? {
            if let Event::Key(key) = event {
                if let Some(msg) = event::handle_key_event(app, key) {
                    handle_message(app, msg)?;
                }
            }
        }

        if app.should_quit {
            break;
        }
    }

    Ok(())
}

fn handle_message(app: &mut App, msg: Message) -> Result<()> {
    match msg {
        Message::Quit => app.should_quit = true,
        Message::NextTab => app.next_tab(),
        Message::PrevTab => app.prev_tab(),
        Message::SelectNext => app.select_next(),
        Message::SelectPrev => app.select_previous(),
        Message::EnterSearch => {
            app.search_mode = true;
            app.search_query.clear();
        }
        Message::ExitSearch => {
            app.search_mode = false;
            // Keep filtered results visible
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
        Message::GenerateConfig => {
            // TODO: Implement in Phase 3
            app.status_message = Some("Config generation not yet implemented".to_string());
        }
        Message::ConvertSelected => {
            // TODO: Implement in Phase 2
            app.status_message = Some("Conversion not yet implemented".to_string());
        }
        _ => {}
    }
    Ok(())
}

fn perform_scan(dir: &PathBuf, max_depth: usize) -> ScanResults {
    let databases = scan_databases(dir, max_depth);
    let (query_files, msa_files) = scan_query_files(dir, max_depth);
    ScanResults {
        databases,
        query_files,
        msa_files,
    }
}

fn update_search_filter(app: &mut App) {
    use nucleo_matcher::{Config, Matcher};
    use nucleo_matcher::pattern::{CaseMatching, Normalization, Pattern};

    if app.search_query.is_empty() {
        app.filtered_indices.clear();
        return;
    }

    let items: Vec<String> = match app.current_tab {
        app::Tab::Databases => app.scan_results.databases.iter().map(|d| d.name.clone()).collect(),
        app::Tab::Queries => app.scan_results.query_files.iter().map(|q| q.name.clone()).collect(),
        app::Tab::MsaFiles => app.scan_results.msa_files.iter().map(|m| m.name.clone()).collect(),
    };

    let mut matcher = Matcher::new(Config::DEFAULT);
    let pattern = Pattern::parse(&app.search_query, CaseMatching::Ignore, Normalization::Smart);

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

    results.sort_by(|a, b| b.1.cmp(&a.1));
    app.filtered_indices = results.into_iter().map(|(idx, _)| idx).collect();
    app.selected_index = 0;
}
```

#### 13. Wire Up TUI Command in Main
**File**: `rust/crates/tronko-bench/src/main.rs`
**Changes**: Add module declaration and command handler

```rust
// Add at top of file after other mod declarations
mod tui;

// Add in main() match block
Commands::Tui { dir, max_depth } => {
    let dir = dir.canonicalize().context("Failed to resolve directory")?;
    tui::run_tui(dir, max_depth)?;
}
```

### Success Criteria:

#### Automated Verification:
- [x] `cd rust && cargo build -p tronko-bench` compiles without errors
- [x] `cd rust && cargo clippy -p tronko-bench` has no warnings (only pre-existing warnings in main.rs)
- [x] `cd rust && cargo test -p tronko-bench` passes

#### Manual Verification:
- [x] `tronko-bench tui` launches and shows database list
- [x] Tab key switches between Databases, Queries, and MSA Files tabs
- [x] j/k keys navigate up/down in the list
- [x] Selected item details show in right panel
- [x] 'q' quits the application cleanly
- [x] Status bar shows key bindings

---

## Phase 2: Database Conversion with Progress

### Overview
Add the ability to convert text databases to binary format with real-time progress display.

### Changes Required:

#### 1. Add Conversion State to App
**File**: `rust/crates/tronko-bench/src/tui/app.rs`

```rust
// Add to App struct:
pub conversion_in_progress: bool,
pub conversion_progress: f64,
pub conversion_target: Option<String>,
```

#### 2. Implement ConvertSelected Handler
**File**: `rust/crates/tronko-bench/src/tui/mod.rs`

```rust
Message::ConvertSelected => {
    if app.current_tab != app::Tab::Databases {
        app.status_message = Some("Conversion only available for databases".to_string());
        return Ok(());
    }

    // Get selected database
    let idx = if app.filtered_indices.is_empty() {
        app.selected_index
    } else {
        match app.filtered_indices.get(app.selected_index) {
            Some(&i) => i,
            None => return Ok(()),
        }
    };

    let db = match app.scan_results.databases.get(idx) {
        Some(db) => db.clone(),
        None => return Ok(()),
    };

    // Only convert text databases
    if db.format != "text" {
        app.status_message = Some("Database is already in binary format".to_string());
        return Ok(());
    }

    // Find converter
    let converter = match crate::find_tronko_convert(None) {
        Some(c) => c,
        None => {
            app.status_message = Some("tronko-convert not found!".to_string());
            return Ok(());
        }
    };

    let target = crate::get_binary_path(&db.path, None);
    app.conversion_target = Some(db.name.clone());
    app.status_message = Some(format!("Converting {}...", db.name));

    // Run conversion (blocking)
    let result = crate::convert_database(&converter, &db.path, &target, false);

    if result.success {
        app.status_message = Some(format!(
            "Converted {} -> {} (saved {})",
            db.name,
            target.file_name().unwrap_or_default().to_string_lossy(),
            humansize::format_size(
                result.text_size.saturating_sub(result.binary_size.unwrap_or(0)),
                humansize::BINARY
            )
        ));
        // Rescan to update list
        app.scan_results = perform_scan(&app.scan_dir, app.max_depth);
    } else {
        app.status_message = Some(format!("Conversion failed: {}", result.message));
    }

    app.conversion_target = None;
}
```

### Success Criteria:

#### Manual Verification:
- [x] Selecting a text database and pressing 'c' starts conversion
- [x] Status bar shows conversion progress message
- [x] On completion, database list refreshes to show binary format
- [x] Error messages display for failed conversions

---

## Phase 3: Config Generation & Help Screen

### Overview
Add benchmark.json generation and a help overlay accessible via '?'.

### Changes Required:

#### 1. Add Help State to App
```rust
// Add to App struct:
pub show_help: bool,
```

#### 2. Implement Help Overlay
**File**: `rust/crates/tronko-bench/src/tui/ui/help.rs`

```rust
use ratatui::{
    layout::Rect,
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, Paragraph},
    Frame,
};

pub fn render_help(f: &mut Frame) {
    let area = centered_rect(60, 18, f.area());
    f.render_widget(Clear, area);

    let help_text = vec![
        Line::from(vec![
            Span::styled("tronko-bench TUI", Style::default().add_modifier(Modifier::BOLD)),
        ]),
        Line::from(""),
        Line::from(vec![
            Span::styled("Tab      ", Style::default().fg(Color::Cyan)),
            Span::raw("Next tab"),
        ]),
        Line::from(vec![
            Span::styled("Shift+Tab", Style::default().fg(Color::Cyan)),
            Span::raw("Previous tab"),
        ]),
        Line::from(vec![
            Span::styled("j/k or ↑/↓", Style::default().fg(Color::Cyan)),
            Span::raw("Navigate list"),
        ]),
        Line::from(vec![
            Span::styled("/        ", Style::default().fg(Color::Cyan)),
            Span::raw("Fuzzy search"),
        ]),
        Line::from(vec![
            Span::styled("r        ", Style::default().fg(Color::Cyan)),
            Span::raw("Rescan directory"),
        ]),
        Line::from(vec![
            Span::styled("c        ", Style::default().fg(Color::Cyan)),
            Span::raw("Convert selected database to binary"),
        ]),
        Line::from(vec![
            Span::styled("g        ", Style::default().fg(Color::Cyan)),
            Span::raw("Generate benchmark.json"),
        ]),
        Line::from(vec![
            Span::styled("?        ", Style::default().fg(Color::Cyan)),
            Span::raw("Show this help"),
        ]),
        Line::from(vec![
            Span::styled("q/Esc    ", Style::default().fg(Color::Cyan)),
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

#### 3. Implement Config Generation
```rust
Message::GenerateConfig => {
    let config_path = app.scan_dir.join("benchmark.json");
    let config = crate::BenchmarkConfig {
        databases: app.scan_results.databases
            .iter()
            .map(|db| crate::DatabaseEntry {
                name: db.name.clone(),
                reference_path: db.path.clone(),
                fasta_path: crate::find_associated_fasta(&db.path),
            })
            .collect(),
        query_sets: app.scan_results.query_files
            .iter()
            .map(|q| crate::QuerySetEntry {
                name: q.name.clone(),
                path: q.path.clone(),
                file_type: q.file_type.clone(),
            })
            .collect(),
    };

    match serde_json::to_string_pretty(&config) {
        Ok(json) => {
            match std::fs::write(&config_path, &json) {
                Ok(_) => {
                    app.status_message = Some(format!(
                        "Wrote {} ({} databases, {} queries)",
                        config_path.display(),
                        config.databases.len(),
                        config.query_sets.len()
                    ));
                }
                Err(e) => {
                    app.status_message = Some(format!("Failed to write config: {}", e));
                }
            }
        }
        Err(e) => {
            app.status_message = Some(format!("Failed to serialize config: {}", e));
        }
    }
}
```

### Success Criteria:

#### Manual Verification:
- [x] Press '?' shows help overlay
- [x] Any key closes help overlay
- [x] Press 'g' generates benchmark.json in scan directory
- [x] Status message shows success/failure of config generation

---

## Phase 4: Polish & Error Handling

### Overview
Improve error display, add empty state handling, and polish the UI.

### Changes Required:

1. **Empty state messages** - Show helpful text when no files found
2. **Error recovery** - Clear errors after timeout or on next action
3. **Path display** - Better truncation in list view
4. **Scroll indicators** - Show when list has more items above/below
5. **Selection persistence** - Remember selected index when switching tabs

### Success Criteria:

#### Automated Verification:
- [x] `cargo build --release -p tronko-bench` compiles
- [x] `cargo clippy -p tronko-bench` has no warnings (only pre-existing warnings in main.rs)
- [x] `cargo test -p tronko-bench` passes

#### Manual Verification:
- [x] Empty directories show "No files found" message
- [x] Status bar shows ?:help hint
- [x] Long paths are truncated intelligently
- [x] Application feels responsive

---

## Testing Strategy

### Unit Tests:
- `tui/app.rs`: Test tab navigation, selection bounds
- Search filtering: Test fuzzy match logic

### Integration Tests:
- Mock scan results for UI testing
- Test message handling without terminal

### Manual Testing Steps:
1. Run `tronko-bench tui -d /path/to/datasets`
2. Verify all three tabs show correct counts
3. Navigate with keyboard, verify selection updates
4. Test fuzzy search filters correctly
5. Convert a text database, verify it appears as binary
6. Generate config, verify JSON is valid
7. Test on empty directory
8. Test error handling (invalid permissions, etc.)

---

## Performance Considerations

- **Lazy scanning**: Only scan on startup and manual rescan
- **Efficient filtering**: nucleo-matcher is optimized for real-time fuzzy search
- **Minimal redraws**: Only redraw on state changes
- **No file content loading**: Only metadata is displayed

---

## Migration Path

The TUI is additive - existing CLI commands remain unchanged:
- `tronko-bench scan` - Works as before
- `tronko-bench databases` - Works as before
- `tronko-bench queries` - Works as before
- `tronko-bench config` - Works as before
- `tronko-bench convert` - Works as before
- `tronko-bench tui` - NEW: Launches interactive mode

Users can choose their preferred workflow without breaking existing scripts.

---

## References

- Reference TUI plan: `thoughts/shared/plans/2025-12-31-gcs-browser-tui.md`
- Existing crate: `rust/crates/tronko-bench/`
- Workspace config: `rust/Cargo.toml`
- Ratatui documentation: https://ratatui.rs/
- nucleo-matcher: https://docs.rs/nucleo-matcher
