# Benchmark Wizard Optimization Toggles Implementation Plan

## Overview

Update the tronko-bench TUI wizard to support the new toggleable optimization features from the Tier 1 optimizations implementation. This enables users to configure runtime optimization flags during benchmark creation for A/B comparison testing.

## Current State Analysis

### Tier 1 Optimizations Available

**Compile-time flags** (via Makefile - binary selection already handles this):
- `NATIVE_ARCH=1` - Architecture-specific optimization
- `FAST_MATH=1` - Fast math optimization
- `ENABLE_OPENMP=1` - OpenMP support
- `OPTIMIZE_MEMORY=1` - Already supported via binary selection

**Runtime flags** (new command-line options to support):
- `--early-termination` / `--no-early-termination`
- `--strike-box [FLOAT]` (default: 1.0)
- `--max-strikes [INT]` (default: 6)
- `--enable-pruning` / `--disable-pruning`
- `--pruning-factor [FLOAT]` (default: 2.0)

### Current Wizard Implementation

**6-step wizard flow:**
1. Select Database
2. Select Project
3. Select Marker
4. Select Query Mode
5. Select Binary (standard vs OPTIMIZE_MEMORY)
6. Confirm & Run

**Current ExecutionParams** (`rust/crates/tronko-bench/src/benchmark/config.rs:160-188`):
```rust
pub struct ExecutionParams {
    pub tronko_assign_path: PathBuf,
    pub optimize_memory: bool,
    pub cores: u32,
    pub use_needleman_wunsch: bool,
    pub batch_size: u32,
}
```

**Current build_command()** (`rust/crates/tronko-bench/src/benchmark/execution.rs:78-122`):
- Does not pass any optimization flags

## Desired End State

After implementation:
1. Users can toggle early termination and pruning optimizations in the wizard
2. Benchmark names reflect optimization settings for easy identification
3. Comparison view highlights optimization differences between benchmarks
4. Command preview in Confirm step shows the optimization flags
5. All changes are backward compatible - existing benchmarks continue to work

### Verification
```
# Start benchmark wizard
tronko-bench

# In wizard, after selecting binary:
# - See new "Optimizations" step or section
# - Toggle early termination ON
# - Toggle pruning ON
# - See command preview with --early-termination --enable-pruning

# After running:
# - Benchmark name shows "+opt" or similar suffix
# - Comparison view shows optimization settings
```

## What We're NOT Doing

- Compile-time flag configuration (user must pre-build binaries)
- Automatic optimization profiling
- Changes to tronko-assign itself

---

## Design Decision: Wizard Step Structure

### Option A: Add New "Optimizations" Step (Recommended)

Add a 7th step between Binary and Confirm:
```
DB → Project → Marker → Mode → Binary → Optimizations → Confirm
```

**Pros:**
- Clear separation of concerns
- Easy to skip with "Use Defaults" shortcut
- Room for future optimization options

**Cons:**
- One more step in wizard

### Option B: Integrate into Binary Step

Combine binary selection with optimization toggles in same step.

**Pros:**
- Fewer steps
- Logical grouping (both affect execution)

**Cons:**
- Cluttered UI
- Binary selection is different from runtime flags

### Option C: Integrate into Confirm Step

Add toggles to the confirmation screen with editable options.

**Pros:**
- Single place to review everything
- Can adjust before running

**Cons:**
- Makes Confirm step more complex
- Mixed purpose (review vs configure)

**Decision: Option A** - Add a dedicated Optimizations step for clarity and extensibility.

---

## Implementation Approach

### Phase 1: Data Model Updates

Update `ExecutionParams` to include optimization flags, using a nested struct for organization.

### Phase 2: Wizard Step Addition

Add new `SelectOptimizations` wizard step with toggle UI.

### Phase 3: Command Builder Update

Update `build_command()` to pass optimization flags.

### Phase 4: Display Updates

Update benchmark naming, confirmation display, and comparison view.

---

## Phase 1: Data Model Updates

### Overview
Add optimization configuration to the ExecutionParams struct, including presets and advanced parameter support.

### Changes Required:

#### 1. Add OptimizationPreset Enum
**File**: `rust/crates/tronko-bench/src/benchmark/config.rs`

Add after line 188 (after `ExecutionParams` impl):
```rust
/// Preset optimization profiles for common use cases
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum OptimizationPreset {
    /// No optimizations - baseline for comparison
    Baseline,
    /// Balanced - conservative optimizations that maintain accuracy
    Balanced,
    /// Aggressive - maximum speed, may have minor accuracy differences
    Aggressive,
    /// Custom - user-defined settings
    Custom,
}

impl OptimizationPreset {
    pub fn all() -> &'static [OptimizationPreset] {
        &[
            OptimizationPreset::Baseline,
            OptimizationPreset::Balanced,
            OptimizationPreset::Aggressive,
            OptimizationPreset::Custom,
        ]
    }

    pub fn label(&self) -> &'static str {
        match self {
            OptimizationPreset::Baseline => "Baseline",
            OptimizationPreset::Balanced => "Balanced",
            OptimizationPreset::Aggressive => "Aggressive",
            OptimizationPreset::Custom => "Custom",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            OptimizationPreset::Baseline => "No optimizations - use for baseline comparison",
            OptimizationPreset::Balanced => "Conservative optimizations, maintains accuracy",
            OptimizationPreset::Aggressive => "Maximum speed, may have minor accuracy differences",
            OptimizationPreset::Custom => "Configure individual optimization settings",
        }
    }

    /// Get the OptimizationConfig for this preset
    pub fn to_config(&self) -> OptimizationConfig {
        match self {
            OptimizationPreset::Baseline => OptimizationConfig::baseline(),
            OptimizationPreset::Balanced => OptimizationConfig::balanced(),
            OptimizationPreset::Aggressive => OptimizationConfig::aggressive(),
            OptimizationPreset::Custom => OptimizationConfig::default(),
        }
    }
}

impl Default for OptimizationPreset {
    fn default() -> Self {
        OptimizationPreset::Baseline
    }
}

/// Runtime optimization settings for tronko-assign
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationConfig {
    /// The preset this config was derived from (for display purposes)
    #[serde(default)]
    pub preset: OptimizationPreset,

    /// Enable early termination during tree traversal
    pub early_termination: bool,

    /// Strike box size as multiplier of Cinterval (only used if early_termination is true)
    pub strike_box: f64,

    /// Maximum strikes before termination (only used if early_termination is true)
    pub max_strikes: u32,

    /// Enable subtree pruning
    pub enable_pruning: bool,

    /// Pruning threshold = factor * Cinterval (only used if enable_pruning is true)
    pub pruning_factor: f64,
}

impl Default for OptimizationConfig {
    fn default() -> Self {
        Self::baseline()
    }
}

impl OptimizationConfig {
    /// Baseline preset - no optimizations
    pub fn baseline() -> Self {
        Self {
            preset: OptimizationPreset::Baseline,
            early_termination: false,
            strike_box: 1.0,
            max_strikes: 6,
            enable_pruning: false,
            pruning_factor: 2.0,
        }
    }

    /// Balanced preset - conservative optimizations
    pub fn balanced() -> Self {
        Self {
            preset: OptimizationPreset::Balanced,
            early_termination: true,
            strike_box: 1.0,       // Same as Cinterval (conservative)
            max_strikes: 6,        // Default value
            enable_pruning: false, // Pruning disabled for accuracy
            pruning_factor: 2.0,
        }
    }

    /// Aggressive preset - maximum speed
    pub fn aggressive() -> Self {
        Self {
            preset: OptimizationPreset::Aggressive,
            early_termination: true,
            strike_box: 0.5,       // More aggressive early termination
            max_strikes: 4,        // Fewer strikes before giving up
            enable_pruning: true,  // Enable pruning
            pruning_factor: 1.5,   // More aggressive pruning threshold
        }
    }

    /// Check if any optimizations are enabled
    pub fn any_enabled(&self) -> bool {
        self.early_termination || self.enable_pruning
    }

    /// Get a short label for display (e.g., "+ET+PR" for early-term + pruning)
    pub fn short_label(&self) -> String {
        match self.preset {
            OptimizationPreset::Baseline => String::new(),
            OptimizationPreset::Balanced => "+bal".to_string(),
            OptimizationPreset::Aggressive => "+agg".to_string(),
            OptimizationPreset::Custom => {
                let mut parts = Vec::new();
                if self.early_termination {
                    parts.push("ET");
                }
                if self.enable_pruning {
                    parts.push("PR");
                }
                if parts.is_empty() {
                    String::new()
                } else {
                    format!("+{}", parts.join("+"))
                }
            }
        }
    }

    /// Check if using non-default advanced parameters
    pub fn has_custom_params(&self) -> bool {
        (self.strike_box - 1.0).abs() > 0.001
            || self.max_strikes != 6
            || (self.pruning_factor - 2.0).abs() > 0.001
    }

    /// Get command-line arguments for these settings
    pub fn to_args(&self) -> Vec<String> {
        let mut args = Vec::new();

        if self.early_termination {
            args.push("--early-termination".to_string());
            // Only include non-default values
            if (self.strike_box - 1.0).abs() > 0.001 {
                args.push(format!("--strike-box={}", self.strike_box));
            }
            if self.max_strikes != 6 {
                args.push(format!("--max-strikes={}", self.max_strikes));
            }
        }

        if self.enable_pruning {
            args.push("--enable-pruning".to_string());
            if (self.pruning_factor - 2.0).abs() > 0.001 {
                args.push(format!("--pruning-factor={}", self.pruning_factor));
            }
        }

        args
    }
}
```

#### 2. Update ExecutionParams Struct
**File**: `rust/crates/tronko-bench/src/benchmark/config.rs`

Add new field to `ExecutionParams` (after line 175, before closing brace):
```rust
    /// Runtime optimization settings
    pub optimizations: OptimizationConfig,
```

Update the `Default` impl (line 178-187):
```rust
impl Default for ExecutionParams {
    fn default() -> Self {
        Self {
            tronko_assign_path: PathBuf::new(),
            optimize_memory: false,
            cores: 1,
            use_needleman_wunsch: false,
            batch_size: 50000,
            optimizations: OptimizationConfig::default(),
        }
    }
}
```

#### 3. Update Benchmark Name Generation
**File**: `rust/crates/tronko-bench/src/benchmark/config.rs`

Update `generate_name()` method (around line 34-57) to include optimization suffix:
```rust
pub fn generate_name(&self) -> String {
    let mode = match self.queries.read_mode {
        ReadMode::Paired => "paired",
        ReadMode::UnpairedForward => "single-F",
        ReadMode::UnpairedReverse => "single-R",
    };

    let format = if self.database.format.contains("binary") {
        "binary"
    } else {
        "text"
    };

    let optmem = if self.params.optimize_memory {
        "+optmem"
    } else {
        ""
    };

    // Add optimization label
    let opt_label = self.params.optimizations.short_label();

    format!(
        "{} {} ({}{}{})",
        self.queries.marker, mode, format, optmem, opt_label
    )
}
```

#### 4. Update mod.rs Exports
**File**: `rust/crates/tronko-bench/src/benchmark/mod.rs`

Add `OptimizationConfig` to the exports if not already using a wildcard.

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors
- [ ] `cargo test -p tronko-bench` passes
- [ ] Existing benchmark JSON files deserialize correctly (backward compat via Default)

#### Manual Verification:
- [ ] `OptimizationConfig::default()` returns disabled optimizations
- [ ] `OptimizationConfig::short_label()` returns correct labels
- [ ] `OptimizationConfig::to_args()` returns correct command-line args

---

## Phase 2: Wizard Step Addition

### Overview
Add a new wizard step for configuring runtime optimizations.

### Changes Required:

#### 1. Add WizardStep Variant
**File**: `rust/crates/tronko-bench/src/tui/app.rs`

Update `WizardStep` enum (around line 44-52):
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WizardStep {
    SelectDatabase,
    SelectProject,
    SelectMarker,
    SelectQueryMode,
    SelectBinary,
    SelectOptimizations,  // NEW
    Confirm,
}
```

Update `WizardStep::number()` (around line 55-64):
```rust
pub fn number(&self) -> usize {
    match self {
        WizardStep::SelectDatabase => 1,
        WizardStep::SelectProject => 2,
        WizardStep::SelectMarker => 3,
        WizardStep::SelectQueryMode => 4,
        WizardStep::SelectBinary => 5,
        WizardStep::SelectOptimizations => 6,
        WizardStep::Confirm => 7,
    }
}
```

Update `WizardStep::title()` (around line 66-75):
```rust
pub fn title(&self) -> &'static str {
    match self {
        WizardStep::SelectDatabase => "Select Reference Database",
        WizardStep::SelectProject => "Select Project",
        WizardStep::SelectMarker => "Select Marker",
        WizardStep::SelectQueryMode => "Select Query Type",
        WizardStep::SelectBinary => "Select tronko-assign Binary",
        WizardStep::SelectOptimizations => "Configure Optimizations",
        WizardStep::Confirm => "Confirm & Run",
    }
}
```

Update `WizardStep::next()` (around line 77-86):
```rust
pub fn next(&self) -> Option<WizardStep> {
    match self {
        WizardStep::SelectDatabase => Some(WizardStep::SelectProject),
        WizardStep::SelectProject => Some(WizardStep::SelectMarker),
        WizardStep::SelectMarker => Some(WizardStep::SelectQueryMode),
        WizardStep::SelectQueryMode => Some(WizardStep::SelectBinary),
        WizardStep::SelectBinary => Some(WizardStep::SelectOptimizations),
        WizardStep::SelectOptimizations => Some(WizardStep::Confirm),
        WizardStep::Confirm => None,
    }
}
```

Update `WizardStep::prev()` (around line 88-97):
```rust
pub fn prev(&self) -> Option<WizardStep> {
    match self {
        WizardStep::SelectDatabase => None,
        WizardStep::SelectProject => Some(WizardStep::SelectDatabase),
        WizardStep::SelectMarker => Some(WizardStep::SelectProject),
        WizardStep::SelectQueryMode => Some(WizardStep::SelectMarker),
        WizardStep::SelectBinary => Some(WizardStep::SelectQueryMode),
        WizardStep::SelectOptimizations => Some(WizardStep::SelectBinary),
        WizardStep::Confirm => Some(WizardStep::SelectOptimizations),
    }
}
```

#### 2. Add Optimization State to WizardState
**File**: `rust/crates/tronko-bench/src/tui/app.rs`

Add to `WizardState` struct (around line 100-113):
```rust
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
    pub confirm_button_selected: bool,

    // NEW: Optimization settings
    pub optimization_config: OptimizationConfig,
    pub selected_preset_index: usize,     // Which preset is highlighted (0-3)
    pub advanced_mode: bool,              // Whether advanced editing is active
    pub advanced_field_index: usize,      // Which advanced field is selected (0-4)
    pub editing_field: bool,              // Whether currently editing a field value
    pub edit_buffer: String,              // Buffer for field editing
}
```

Update `Default` impl:
```rust
impl Default for WizardState {
    fn default() -> Self {
        Self {
            step: WizardStep::SelectDatabase,
            list_index: 0,
            selected_db_index: None,
            discovered_projects: Vec::new(),
            selected_project_index: None,
            selected_marker_index: None,
            selected_query_mode: None,
            tronko_binaries: Vec::new(),
            selected_binary_index: None,
            confirm_button_selected: true,
            optimization_config: OptimizationConfig::default(),
            selected_preset_index: 0,
            advanced_mode: false,
            advanced_field_index: 0,
            editing_field: false,
            edit_buffer: String::new(),
        }
    }
}
```

#### 3. Update Wizard Progress Display
**File**: `rust/crates/tronko-bench/src/tui/ui/wizard.rs`

Update the progress indicator (around line 64):
```rust
fn render_progress(f: &mut Frame, area: Rect, _app: &App) {
    let steps = ["DB", "Project", "Marker", "Mode", "Binary", "Opts", "Confirm"];
    // ... rest unchanged
}
```

#### 4. Add Optimization Step Rendering
**File**: `rust/crates/tronko-bench/src/tui/ui/wizard.rs`

Add to `render_step_content()` match (around line 88-96):
```rust
fn render_step_content(f: &mut Frame, area: Rect, app: &App) {
    match app.wizard.step {
        WizardStep::SelectDatabase => render_database_selection(f, area, app),
        WizardStep::SelectProject => render_project_selection(f, area, app),
        WizardStep::SelectMarker => render_marker_selection(f, area, app),
        WizardStep::SelectQueryMode => render_query_mode_selection(f, area, app),
        WizardStep::SelectBinary => render_binary_selection(f, area, app),
        WizardStep::SelectOptimizations => render_optimization_selection(f, area, app),
        WizardStep::Confirm => render_confirmation(f, area, app),
    }
}
```

Add new rendering function:
```rust
fn render_optimization_selection(f: &mut Frame, area: Rect, app: &App) {
    let config = &app.wizard.optimization_config;

    // Split area: presets on left, details/advanced on right
    let main_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage(45),  // Preset selection
            Constraint::Percentage(55),  // Details / Advanced
        ])
        .split(area);

    // Left panel: Preset selection
    render_preset_list(f, main_chunks[0], app);

    // Right panel: Details or Advanced editing
    if app.wizard.advanced_mode {
        render_advanced_options(f, main_chunks[1], app);
    } else {
        render_preset_details(f, main_chunks[1], app);
    }
}

fn render_preset_list(f: &mut Frame, area: Rect, app: &App) {
    let presets = OptimizationPreset::all();

    let items: Vec<ListItem> = presets
        .iter()
        .enumerate()
        .map(|(i, preset)| {
            let selected = i == app.wizard.selected_preset_index;
            let indicator = if selected && app.wizard.optimization_config.preset == *preset {
                "(●)"  // Currently active
            } else if app.wizard.optimization_config.preset == *preset {
                "(○)"  // This preset is applied
            } else {
                "( )"
            };

            ListItem::new(vec![
                Line::from(vec![
                    Span::styled(
                        format!("{} ", indicator),
                        if app.wizard.optimization_config.preset == *preset {
                            Style::default().fg(Color::Green)
                        } else {
                            Style::default().fg(Color::DarkGray)
                        }
                    ),
                    Span::styled(
                        preset.label(),
                        Style::default().add_modifier(Modifier::BOLD)
                    ),
                ]),
                Line::from(vec![
                    Span::raw("    "),
                    Span::styled(preset.description(), Style::default().fg(Color::DarkGray)),
                ]),
            ])
        })
        .collect();

    let title = if app.wizard.advanced_mode {
        " Presets (Tab to exit advanced) "
    } else {
        " Select Optimization Preset "
    };

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL).title(title))
        .highlight_style(
            Style::default()
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        )
        .highlight_symbol("> ");

    let mut state = ListState::default();
    if !app.wizard.advanced_mode {
        state.select(Some(app.wizard.selected_preset_index));
    }

    f.render_stateful_widget(list, area, &mut state);
}

fn render_preset_details(f: &mut Frame, area: Rect, app: &App) {
    let config = &app.wizard.optimization_config;

    let mut lines = vec![
        Line::from(Span::styled(
            format!("Preset: {}", config.preset.label()),
            Style::default().add_modifier(Modifier::BOLD)
        )),
        Line::from(""),
    ];

    // Show current settings
    lines.push(Line::from(vec![
        Span::styled("Early Termination: ", Style::default().fg(Color::Cyan)),
        Span::styled(
            if config.early_termination { "ON" } else { "OFF" },
            if config.early_termination { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }
        ),
    ]));

    if config.early_termination {
        lines.push(Line::from(vec![
            Span::raw("  Strike Box: "),
            Span::styled(format!("{:.1}", config.strike_box), Style::default().fg(Color::Yellow)),
        ]));
        lines.push(Line::from(vec![
            Span::raw("  Max Strikes: "),
            Span::styled(format!("{}", config.max_strikes), Style::default().fg(Color::Yellow)),
        ]));
    }

    lines.push(Line::from(""));
    lines.push(Line::from(vec![
        Span::styled("Subtree Pruning: ", Style::default().fg(Color::Cyan)),
        Span::styled(
            if config.enable_pruning { "ON" } else { "OFF" },
            if config.enable_pruning { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }
        ),
    ]));

    if config.enable_pruning {
        lines.push(Line::from(vec![
            Span::raw("  Pruning Factor: "),
            Span::styled(format!("{:.1}", config.pruning_factor), Style::default().fg(Color::Yellow)),
        ]));
    }

    // Command preview
    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled("─".repeat(30), Style::default().fg(Color::DarkGray))));
    lines.push(Line::from(Span::styled("Command flags:", Style::default().fg(Color::DarkGray))));

    let args = config.to_args();
    if args.is_empty() {
        lines.push(Line::from(Span::styled("  (none)", Style::default().fg(Color::DarkGray))));
    } else {
        for arg in args {
            lines.push(Line::from(Span::styled(format!("  {}", arg), Style::default().fg(Color::Cyan))));
        }
    }

    // Hint for advanced mode
    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled(
        "Press 'a' for advanced options",
        Style::default().fg(Color::DarkGray)
    )));

    let paragraph = Paragraph::new(lines)
        .block(Block::default().borders(Borders::ALL).title(" Settings "));
    f.render_widget(paragraph, area);
}

fn render_advanced_options(f: &mut Frame, area: Rect, app: &App) {
    let config = &app.wizard.optimization_config;

    // Advanced parameter fields
    let fields = vec![
        ("Early Termination", format!("{}", if config.early_termination { "ON" } else { "OFF" }), true),
        ("  Strike Box", format!("{:.2}", config.strike_box), config.early_termination),
        ("  Max Strikes", format!("{}", config.max_strikes), config.early_termination),
        ("Subtree Pruning", format!("{}", if config.enable_pruning { "ON" } else { "OFF" }), true),
        ("  Pruning Factor", format!("{:.2}", config.pruning_factor), config.enable_pruning),
    ];

    let items: Vec<ListItem> = fields
        .iter()
        .enumerate()
        .map(|(i, (label, value, enabled))| {
            let is_selected = i == app.wizard.advanced_field_index;
            let is_editing = is_selected && app.wizard.editing_field;

            let display_value = if is_editing {
                format!("{}▌", app.wizard.edit_buffer)
            } else {
                value.clone()
            };

            let value_style = if !enabled {
                Style::default().fg(Color::DarkGray)
            } else if is_editing {
                Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)
            } else if is_selected {
                Style::default().fg(Color::Cyan)
            } else {
                Style::default().fg(Color::White)
            };

            ListItem::new(Line::from(vec![
                Span::styled(
                    format!("{:<18}", label),
                    if *enabled { Style::default() } else { Style::default().fg(Color::DarkGray) }
                ),
                Span::styled(display_value, value_style),
            ]))
        })
        .collect();

    let list = List::new(items)
        .block(Block::default().borders(Borders::ALL).title(" Advanced Options (Tab to exit) "))
        .highlight_style(
            Style::default()
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD),
        )
        .highlight_symbol("> ");

    let mut state = ListState::default();
    state.select(Some(app.wizard.advanced_field_index));

    // Split for list and help
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Min(7),
            Constraint::Length(4),
        ])
        .split(area);

    f.render_stateful_widget(list, chunks[0], &mut state);

    // Help text
    let help = Paragraph::new(vec![
        Line::from(Span::styled("Space: toggle on/off", Style::default().fg(Color::DarkGray))),
        Line::from(Span::styled("Enter: edit value", Style::default().fg(Color::DarkGray))),
        Line::from(Span::styled("Tab: exit advanced mode", Style::default().fg(Color::DarkGray))),
    ])
    .block(Block::default().borders(Borders::TOP));
    f.render_widget(help, chunks[1]);
}
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors
- [ ] Wizard displays 7 steps in progress indicator

#### Manual Verification:
- [ ] Optimizations step appears after Binary selection
- [ ] Space key toggles optimization checkboxes
- [ ] Info panel updates to show enabled flags
- [ ] Can navigate back and forward through the step

---

## Phase 3: Event Handling Updates

### Overview
Update event handling to support presets, advanced mode, and field editing.

### Changes Required:

#### 1. Add Preset and Advanced Mode Handlers
**File**: `rust/crates/tronko-bench/src/tui/event.rs`

Find the wizard key handling section and add comprehensive handling for the Optimizations step:

```rust
// Handle Optimization step - depends on mode (preset selection vs advanced)
if app.wizard.step == WizardStep::SelectOptimizations {
    if app.wizard.editing_field {
        // Field editing mode
        match key.code {
            KeyCode::Enter => {
                // Apply the edited value
                apply_edited_value(app);
                app.wizard.editing_field = false;
                app.wizard.edit_buffer.clear();
            }
            KeyCode::Esc => {
                // Cancel editing
                app.wizard.editing_field = false;
                app.wizard.edit_buffer.clear();
            }
            KeyCode::Backspace => {
                app.wizard.edit_buffer.pop();
            }
            KeyCode::Char(c) => {
                // Only allow valid characters for numeric input
                if c.is_ascii_digit() || c == '.' {
                    app.wizard.edit_buffer.push(c);
                }
            }
            _ => {}
        }
    } else if app.wizard.advanced_mode {
        // Advanced mode navigation
        match key.code {
            KeyCode::Tab => {
                // Exit advanced mode
                app.wizard.advanced_mode = false;
                app.wizard.optimization_config.preset = OptimizationPreset::Custom;
            }
            KeyCode::Down | KeyCode::Char('j') => {
                if app.wizard.advanced_field_index < 4 {
                    app.wizard.advanced_field_index += 1;
                }
            }
            KeyCode::Up | KeyCode::Char('k') => {
                app.wizard.advanced_field_index = app.wizard.advanced_field_index.saturating_sub(1);
            }
            KeyCode::Char(' ') => {
                // Toggle boolean fields (indices 0 and 3)
                match app.wizard.advanced_field_index {
                    0 => {
                        app.wizard.optimization_config.early_termination =
                            !app.wizard.optimization_config.early_termination;
                        app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                    }
                    3 => {
                        app.wizard.optimization_config.enable_pruning =
                            !app.wizard.optimization_config.enable_pruning;
                        app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                    }
                    _ => {}
                }
            }
            KeyCode::Enter => {
                // Start editing numeric fields (indices 1, 2, 4)
                match app.wizard.advanced_field_index {
                    1 => {
                        if app.wizard.optimization_config.early_termination {
                            app.wizard.editing_field = true;
                            app.wizard.edit_buffer = format!("{:.2}", app.wizard.optimization_config.strike_box);
                        }
                    }
                    2 => {
                        if app.wizard.optimization_config.early_termination {
                            app.wizard.editing_field = true;
                            app.wizard.edit_buffer = format!("{}", app.wizard.optimization_config.max_strikes);
                        }
                    }
                    4 => {
                        if app.wizard.optimization_config.enable_pruning {
                            app.wizard.editing_field = true;
                            app.wizard.edit_buffer = format!("{:.2}", app.wizard.optimization_config.pruning_factor);
                        }
                    }
                    _ => {}
                }
            }
            KeyCode::Esc => {
                // Exit advanced mode on Esc
                app.wizard.advanced_mode = false;
            }
            _ => {}
        }
    } else {
        // Preset selection mode
        match key.code {
            KeyCode::Down | KeyCode::Char('j') => {
                if app.wizard.selected_preset_index < 3 {
                    app.wizard.selected_preset_index += 1;
                }
            }
            KeyCode::Up | KeyCode::Char('k') => {
                app.wizard.selected_preset_index = app.wizard.selected_preset_index.saturating_sub(1);
            }
            KeyCode::Enter | KeyCode::Char(' ') => {
                // Apply selected preset
                let preset = OptimizationPreset::all()[app.wizard.selected_preset_index];
                if preset == OptimizationPreset::Custom {
                    // Enter advanced mode for custom
                    app.wizard.advanced_mode = true;
                    app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                } else {
                    // Apply preset config
                    app.wizard.optimization_config = preset.to_config();
                }
            }
            KeyCode::Char('a') => {
                // Enter advanced mode
                app.wizard.advanced_mode = true;
                app.wizard.optimization_config.preset = OptimizationPreset::Custom;
            }
            KeyCode::Right | KeyCode::Tab => {
                // Advance to Confirm step
                app.wizard.step = WizardStep::Confirm;
                app.wizard.list_index = 0;
            }
            KeyCode::Left | KeyCode::Backspace => {
                // Go back to Binary step
                app.wizard.step = WizardStep::SelectBinary;
            }
            KeyCode::Esc => {
                // Cancel wizard
                app.mode = AppMode::Normal;
                app.wizard.reset();
            }
            _ => {}
        }
    }
    return; // Handled optimization step
}

/// Apply the value from edit_buffer to the appropriate field
fn apply_edited_value(app: &mut App) {
    match app.wizard.advanced_field_index {
        1 => {
            // Strike Box (f64)
            if let Ok(val) = app.wizard.edit_buffer.parse::<f64>() {
                if val > 0.0 && val <= 10.0 {
                    app.wizard.optimization_config.strike_box = val;
                    app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                }
            }
        }
        2 => {
            // Max Strikes (u32)
            if let Ok(val) = app.wizard.edit_buffer.parse::<u32>() {
                if val > 0 && val <= 100 {
                    app.wizard.optimization_config.max_strikes = val;
                    app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                }
            }
        }
        4 => {
            // Pruning Factor (f64)
            if let Ok(val) = app.wizard.edit_buffer.parse::<f64>() {
                if val > 0.0 && val <= 10.0 {
                    app.wizard.optimization_config.pruning_factor = val;
                    app.wizard.optimization_config.preset = OptimizationPreset::Custom;
                }
            }
        }
        _ => {}
    }
}
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors

#### Manual Verification:
- [ ] Preset selection: Up/Down navigates, Enter/Space applies preset
- [ ] Pressing 'a' enters advanced mode
- [ ] In advanced mode: Space toggles booleans, Enter edits numeric fields
- [ ] Tab exits advanced mode
- [ ] Field editing: type numbers, Enter confirms, Esc cancels
- [ ] Right/Tab advances to Confirm, Left goes back to Binary

---

## Phase 4: Command Builder Update

### Overview
Update the benchmark execution to include optimization flags in the tronko-assign command.

### Changes Required:

#### 1. Update build_command Function
**File**: `rust/crates/tronko-bench/src/benchmark/execution.rs`

Update `build_command()` (around line 78-122) to add optimization flags:

```rust
fn build_command(config: &BenchmarkConfig, output_path: &Path) -> Result<Command> {
    let mut cmd = Command::new(&config.params.tronko_assign_path);

    // Reference database
    cmd.arg("-f").arg(&config.database.reference_path);

    // FASTA index
    cmd.arg("-a").arg(&config.database.fasta_path);

    // Query files based on mode
    match (&config.queries.files, config.queries.read_mode) {
        (QueryFiles::Paired { forward, reverse }, ReadMode::Paired) => {
            cmd.arg("-p")
                .arg("-1").arg(forward)
                .arg("-2").arg(reverse);
        }
        (QueryFiles::Single { path }, ReadMode::UnpairedForward) => {
            cmd.arg("-s").arg("-g").arg(path);
        }
        (QueryFiles::Single { path }, ReadMode::UnpairedReverse) => {
            cmd.arg("-s").arg("-g").arg(path);
        }
        _ => {
            return Err(anyhow::anyhow!("Mismatched query files and read mode"));
        }
    }

    // Output file
    cmd.arg("-o").arg(output_path);

    // Additional flags
    cmd.arg("-r"); // Read taxonomic assignments

    if config.params.use_needleman_wunsch {
        cmd.arg("-w");
    }

    // Threads
    cmd.arg("-n").arg(config.params.cores.to_string());

    // Batch size
    cmd.arg("-b").arg(config.params.batch_size.to_string());

    // NEW: Add optimization flags
    for arg in config.params.optimizations.to_args() {
        cmd.arg(arg);
    }

    Ok(cmd)
}
```

#### 2. Update create_config_from_wizard
**File**: `rust/crates/tronko-bench/src/benchmark/execution.rs`

Update `create_config_from_wizard()` (around line 278-284) to include optimization config:

```rust
params: super::ExecutionParams {
    tronko_assign_path: binary.path.clone(),
    optimize_memory: binary.optimize_memory,
    cores: 1,
    use_needleman_wunsch: false,
    batch_size: 50000,
    optimizations: app.wizard.optimization_config.clone(),
},
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors

#### Manual Verification:
- [ ] Run benchmark with early termination enabled, verify `--early-termination` in logs
- [ ] Run benchmark with pruning enabled, verify `--enable-pruning` in logs
- [ ] Run baseline benchmark, verify no optimization flags in command

---

## Phase 5: Display Updates

### Overview
Update the confirmation step and comparison view to show optimization settings.

### Changes Required:

#### 1. Update Confirmation Display
**File**: `rust/crates/tronko-bench/src/tui/ui/wizard.rs`

Update `render_confirmation()` (around line 352-425) to show optimizations:

Add after the Binary section (around line 412):
```rust
// Optimizations
let opt_config = &app.wizard.optimization_config;
lines.push(Line::from(""));
lines.push(Line::from(vec![
    Span::styled("  Optimizations:", Style::default().fg(Color::Cyan)),
]));

if opt_config.any_enabled() {
    if opt_config.early_termination {
        lines.push(Line::from(vec![
            Span::raw("    "),
            Span::styled("✓", Style::default().fg(Color::Green)),
            Span::raw(" Early Termination"),
        ]));
    }
    if opt_config.enable_pruning {
        lines.push(Line::from(vec![
            Span::raw("    "),
            Span::styled("✓", Style::default().fg(Color::Green)),
            Span::raw(" Subtree Pruning"),
        ]));
    }
} else {
    lines.push(Line::from(vec![
        Span::raw("    "),
        Span::styled("None (baseline mode)", Style::default().fg(Color::DarkGray)),
    ]));
}

// Show command preview
lines.push(Line::from(""));
lines.push(Line::from(vec![
    Span::styled("Command preview:", Style::default().add_modifier(Modifier::BOLD)),
]));

// Build command args for preview
let mut cmd_parts = vec!["tronko-assign".to_string()];
// ... add key args ...
for arg in opt_config.to_args() {
    cmd_parts.push(arg);
}
lines.push(Line::from(vec![
    Span::styled(
        format!("  {}", cmd_parts.join(" ")),
        Style::default().fg(Color::DarkGray)
    ),
]));
```

#### 2. Update Benchmark Detail Display
**File**: `rust/crates/tronko-bench/src/tui/ui/detail.rs`

Add optimization info to the benchmark detail panel. Find where benchmark details are rendered and add:

```rust
// Optimizations section
lines.push(Line::from(""));
lines.push(Line::from(vec![
    Span::styled("Optimizations:", Style::default().fg(Color::Cyan)),
]));

if config.params.optimizations.any_enabled() {
    for arg in config.params.optimizations.to_args() {
        lines.push(Line::from(vec![
            Span::raw("  "),
            Span::styled(arg, Style::default().fg(Color::Green)),
        ]));
    }
} else {
    lines.push(Line::from(vec![
        Span::raw("  "),
        Span::styled("None (baseline)", Style::default().fg(Color::DarkGray)),
    ]));
}
```

#### 3. Update Comparison View
**File**: `rust/crates/tronko-bench/src/tui/ui/comparison.rs`

Add a section to compare optimization settings between benchmarks:

```rust
// After performance comparison section, add:

// Optimization Comparison
lines.push(Line::from(""));
lines.push(Line::from(vec![
    Span::styled("OPTIMIZATION SETTINGS", Style::default().add_modifier(Modifier::BOLD)),
]));
lines.push(Line::from("─".repeat(60)));

let opt_a = &config_a.params.optimizations;
let opt_b = &config_b.params.optimizations;

// Early Termination row
let et_a = if opt_a.early_termination { "ON" } else { "OFF" };
let et_b = if opt_b.early_termination { "ON" } else { "OFF" };
let et_diff = opt_a.early_termination != opt_b.early_termination;
lines.push(Line::from(vec![
    Span::styled(format!("{:<20}", "Early Termination"), Style::default()),
    Span::styled(format!("{:<15}", et_a), if opt_a.early_termination { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }),
    Span::styled(format!("{:<15}", et_b), if opt_b.early_termination { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }),
    if et_diff { Span::styled("*DIFF*", Style::default().fg(Color::Yellow)) } else { Span::raw("") },
]));

// Pruning row
let pr_a = if opt_a.enable_pruning { "ON" } else { "OFF" };
let pr_b = if opt_b.enable_pruning { "ON" } else { "OFF" };
let pr_diff = opt_a.enable_pruning != opt_b.enable_pruning;
lines.push(Line::from(vec![
    Span::styled(format!("{:<20}", "Subtree Pruning"), Style::default()),
    Span::styled(format!("{:<15}", pr_a), if opt_a.enable_pruning { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }),
    Span::styled(format!("{:<15}", pr_b), if opt_b.enable_pruning { Style::default().fg(Color::Green) } else { Style::default().fg(Color::DarkGray) }),
    if pr_diff { Span::styled("*DIFF*", Style::default().fg(Color::Yellow)) } else { Span::raw("") },
]));
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors

#### Manual Verification:
- [ ] Confirmation step shows optimization settings
- [ ] Benchmark detail panel shows optimization flags used
- [ ] Comparison view highlights optimization differences

---

## Phase 6: Help & Documentation Updates

### Overview
Update help text and status hints to document new optimization features.

### Changes Required:

#### 1. Update Help Overlay
**File**: `rust/crates/tronko-bench/src/tui/ui/help.rs`

Add new section for optimization-related keybindings:

```rust
// In the help text rendering:
lines.push(Line::from(""));
lines.push(Line::from(vec![
    Span::styled("OPTIMIZATION STEP", Style::default().add_modifier(Modifier::BOLD)),
]));
lines.push(Line::from("─".repeat(40)));
lines.push(Line::from("  Space        Toggle optimization on/off"));
lines.push(Line::from("  Enter        Continue to confirmation"));
lines.push(Line::from("  ←/Backspace  Go back to binary selection"));
```

#### 2. Update Status Bar Hints
**File**: `rust/crates/tronko-bench/src/tui/ui/status.rs`

Add context-specific hints for the optimization step:

```rust
// When in wizard mode, check for optimization step:
AppMode::BenchmarkWizard if app.wizard.step == WizardStep::SelectOptimizations => {
    "Space:toggle ↑↓:navigate Enter:continue ←:back Esc:cancel"
}
```

### Success Criteria:

#### Automated Verification:
- [ ] `cargo build -p tronko-bench` compiles without errors

#### Manual Verification:
- [ ] Help overlay shows optimization keybindings
- [ ] Status bar shows appropriate hints during optimization step

---

## Testing Strategy

### Unit Tests:
- Test `OptimizationPreset::to_config()` returns correct config for each preset
- Test `OptimizationConfig::to_args()` with various combinations
- Test `OptimizationConfig::short_label()` output for each preset and custom configs
- Test `OptimizationConfig::any_enabled()` logic
- Test `OptimizationConfig::has_custom_params()` logic
- Test backward compatibility: deserialize old benchmark JSON without optimization field

### Integration Tests:
```bash
# Create benchmark with Baseline preset
# Verify command doesn't include optimization flags

# Create benchmark with Balanced preset
# Verify --early-termination in command

# Create benchmark with Aggressive preset
# Verify both flags and non-default params in command

# Create benchmark with Custom config
# Verify specific flags in command

# Compare baseline vs aggressive
# Verify comparison view shows differences
```

### Manual Testing Checklist:

#### Preset Selection
- [ ] Start wizard, navigate to Optimizations step
- [ ] Verify Baseline is selected by default
- [ ] Navigate to Balanced preset, verify right panel shows settings
- [ ] Select Balanced preset, verify (●) indicator and command flags update
- [ ] Navigate to Aggressive, select it, verify all optimizations shown
- [ ] Navigate to Custom, verify it enters advanced mode

#### Advanced Mode
- [ ] Press 'a' from preset selection to enter advanced mode
- [ ] Verify Tab exits advanced mode
- [ ] Navigate to Early Termination, toggle with Space
- [ ] Navigate to Strike Box, press Enter to edit
- [ ] Type "0.75", press Enter to confirm
- [ ] Verify Strike Box shows 0.75
- [ ] Edit Pruning Factor, press Esc to cancel, verify value unchanged
- [ ] Verify preset shows "Custom" after manual changes

#### Navigation & Persistence
- [ ] Go back to Binary step, then forward - verify settings preserved
- [ ] Complete wizard with Balanced preset, verify name shows "+bal"
- [ ] Complete wizard with custom config, verify name shows "+ET+PR" etc.

#### Display & Comparison
- [ ] View benchmark details, verify preset name and optimization flags shown
- [ ] Create second benchmark with different preset
- [ ] Compare both benchmarks, verify optimization diff highlighted
- [ ] Restart app, verify old benchmarks without optimizations load correctly

---

## UI Mockups

### Optimization Step - Preset Selection (Default View)

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 6 of 7: Configure Optimizations                                   │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  DB ✓ → Project ✓ → Marker ✓ → Mode ✓ → Binary ✓ → [Opts] → Confirm     │ │
│  │                                                                         │ │
│  │  ┌─ Select Optimization Preset ──────┐ ┌─ Settings ─────────────────┐   │ │
│  │  │                                   │ │                            │   │ │
│  │  │ > (●) Baseline                    │ │ Preset: Baseline           │   │ │
│  │  │       No optimizations - use for  │ │                            │   │ │
│  │  │       baseline comparison         │ │ Early Termination: OFF     │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Balanced                    │ │ Subtree Pruning: OFF       │   │ │
│  │  │       Conservative optimizations, │ │                            │   │ │
│  │  │       maintains accuracy          │ │ ────────────────────────   │   │ │
│  │  │                                   │ │ Command flags:             │   │ │
│  │  │   ( ) Aggressive                  │ │   (none)                   │   │ │
│  │  │       Maximum speed, may have     │ │                            │   │ │
│  │  │       minor accuracy differences  │ │                            │   │ │
│  │  │                                   │ │ Press 'a' for advanced     │   │ │
│  │  │   ( ) Custom                      │ │ options                    │   │ │
│  │  │       Configure individual        │ │                            │   │ │
│  │  │       optimization settings       │ │                            │   │ │
│  │  │                                   │ │                            │   │ │
│  │  └───────────────────────────────────┘ └────────────────────────────┘   │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓:navigate  Enter/Space:select  a:advanced  →/Tab:continue  ←:back    │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Optimization Step - Balanced Preset Selected

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 6 of 7: Configure Optimizations                                   │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  DB ✓ → Project ✓ → Marker ✓ → Mode ✓ → Binary ✓ → [Opts] → Confirm     │ │
│  │                                                                         │ │
│  │  ┌─ Select Optimization Preset ──────┐ ┌─ Settings ─────────────────┐   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Baseline                    │ │ Preset: Balanced           │   │ │
│  │  │       No optimizations - use for  │ │                            │   │ │
│  │  │       baseline comparison         │ │ Early Termination: ON      │   │ │
│  │  │                                   │ │   Strike Box: 1.0          │   │ │
│  │  │ > (●) Balanced                    │ │   Max Strikes: 6           │   │ │
│  │  │       Conservative optimizations, │ │                            │   │ │
│  │  │       maintains accuracy          │ │ Subtree Pruning: OFF       │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Aggressive                  │ │ ────────────────────────   │   │ │
│  │  │       Maximum speed, may have     │ │ Command flags:             │   │ │
│  │  │       minor accuracy differences  │ │   --early-termination      │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Custom                      │ │ Press 'a' for advanced     │   │ │
│  │  │       Configure individual        │ │ options                    │   │ │
│  │  │       optimization settings       │ │                            │   │ │
│  │  │                                   │ │                            │   │ │
│  │  └───────────────────────────────────┘ └────────────────────────────┘   │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓:navigate  Enter/Space:select  a:advanced  →/Tab:continue  ←:back    │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Optimization Step - Advanced Mode

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 6 of 7: Configure Optimizations                                   │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  DB ✓ → Project ✓ → Marker ✓ → Mode ✓ → Binary ✓ → [Opts] → Confirm     │ │
│  │                                                                         │ │
│  │  ┌─ Presets (Tab to exit advanced) ──┐ ┌─ Advanced Options ─────────┐   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Baseline                    │ │ > Early Termination  ON    │   │ │
│  │  │   ( ) Balanced                    │ │     Strike Box       0.50  │   │ │
│  │  │   ( ) Aggressive                  │ │     Max Strikes      4     │   │ │
│  │  │   (○) Custom                      │ │   Subtree Pruning    ON    │   │ │
│  │  │                                   │ │     Pruning Factor   1.50  │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │                                   │ │ ────────────────────────   │   │ │
│  │  │                                   │ │ Space: toggle on/off       │   │ │
│  │  │                                   │ │ Enter: edit value          │   │ │
│  │  │                                   │ │ Tab: exit advanced mode    │   │ │
│  │  │                                   │ │                            │   │ │
│  │  └───────────────────────────────────┘ └────────────────────────────┘   │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  ↑↓:navigate  Space:toggle  Enter:edit  Tab:exit advanced  Esc:cancel   │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Optimization Step - Editing a Value

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 6 of 7: Configure Optimizations                                   │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  DB ✓ → Project ✓ → Marker ✓ → Mode ✓ → Binary ✓ → [Opts] → Confirm     │ │
│  │                                                                         │ │
│  │  ┌─ Presets (Tab to exit advanced) ──┐ ┌─ Advanced Options ─────────┐   │ │
│  │  │                                   │ │                            │   │ │
│  │  │   ( ) Baseline                    │ │   Early Termination  ON    │   │ │
│  │  │   ( ) Balanced                    │ │ >   Strike Box       0.75▌ │   │ │
│  │  │   ( ) Aggressive                  │ │     Max Strikes      4     │   │ │
│  │  │   (○) Custom                      │ │   Subtree Pruning    ON    │   │ │
│  │  │                                   │ │     Pruning Factor   1.50  │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │                                   │ │ ────────────────────────   │   │ │
│  │  │                                   │ │ Enter: confirm value       │   │ │
│  │  │                                   │ │ Esc: cancel edit           │   │ │
│  │  │                                   │ │                            │   │ │
│  │  │                                   │ │                            │   │ │
│  │  └───────────────────────────────────┘ └────────────────────────────┘   │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Type value, Enter:confirm, Esc:cancel                                  │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Updated Confirmation Step

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                                                                              │
│  ┌─ New Benchmark ─────────────────────────────────────────────────────────┐ │
│  │                                                                         │ │
│  │  Step 7 of 7: Confirm & Run                                             │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  DB ✓ → Project ✓ → Marker ✓ → Mode ✓ → Binary ✓ → Opts ✓ → [Confirm]   │ │
│  │                                                                         │ │
│  │  Review your benchmark configuration:                                   │ │
│  │                                                                         │ │
│  │  Database:     16S_Bacteria/reference_tree.trkb                         │ │
│  │  Project:      16S                                                      │ │
│  │  Marker:       16S                                                      │ │
│  │  Mode:         Paired-end                                               │ │
│  │  Binary:       tronko-assign (standard)                                 │ │
│  │                                                                         │ │
│  │  Optimizations:                                                         │ │
│  │    ✓ Early Termination                                                  │ │
│  │    ✓ Subtree Pruning                                                    │ │
│  │                                                                         │ │
│  │  Command:                                                               │ │
│  │    tronko-assign -f ... --early-termination --enable-pruning            │ │
│  │                                                                         │ │
│  │  ─────────────────────────────────────────────────────────────────────  │ │
│  │  Press Enter to start the benchmark, or Esc to cancel.                  │ │
│  └─────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### Benchmark List with Optimization Labels

```
┌─ Benchmarks ─────────────────────────────────────────────────────────────────┐
│                                                                              │
│ > [ ] 16S paired (binary+agg)          38.1s     ← Aggressive preset         │
│       42.3 MB   2 min ago         [OK]                                       │
│   [*] 16S paired (binary+bal)          41.2s     ← Balanced preset           │
│       45.1 MB   5 min ago         [OK]                                       │
│   [*] 16S paired (binary)              89.7s     ← Baseline (no suffix)      │
│       128.4 MB  1 hour ago        [OK]                                       │
│   [ ] 16S paired (binary+ET+PR)        36.5s     ← Custom config             │
│       40.1 MB   2 hours ago       [OK]                                       │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

**Label Key:**
- `(no suffix)` - Baseline preset (no optimizations)
- `+bal` - Balanced preset
- `+agg` - Aggressive preset
- `+ET`, `+PR`, `+ET+PR` - Custom with specific optimizations enabled

---

## Migration Notes

- **Backward Compatibility**: Existing benchmark JSON files will deserialize correctly due to `#[serde(default)]` on the new `optimizations` field
- **No Database Changes**: This affects only tronko-bench, not tronko-assign or tronko-build
- **Graceful Degradation**: If running against an older tronko-assign without optimization flags, the flags will be ignored (tronko-assign will error if flags are unknown - user should update)

---

## References

- Tier 1 Optimizations Plan: `thoughts/shared/plans/2026-01-01-tier1-optimizations-toggleable.md`
- Benchmark Workflow Plan: `thoughts/plans/2026-01-01-tronko-bench-unified-benchmarking-workflow.md`
- Current config: `rust/crates/tronko-bench/src/benchmark/config.rs`
- Current wizard: `rust/crates/tronko-bench/src/tui/ui/wizard.rs`
- Current execution: `rust/crates/tronko-bench/src/benchmark/execution.rs`
