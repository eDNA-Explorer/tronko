# Tronko Documentation Plan

This document outlines our comprehensive approach to documenting the Tronko codebase, with a particular focus on information flow between components, assignment generation processes, scoring mechanisms, and read filtering/discarding.

## 1. Core Components to Document

### 1.1. tronko-build
**Purpose:** Documents how reference databases are constructed for taxonomic assignment.

**Key files to analyze:**
- `tronko-build/tronko-build.c` - Main program logic
- `tronko-build/readfasta.c` - FASTA sequence parsing
- `tronko-build/readreference.c` - Reference tree handling
- `tronko-build/likelihood.c` - Likelihood calculations

**Documentation focus:**
- Database building workflow
- Tree partitioning algorithm
- Fractional likelihood calculations and storage
- Format of the output reference database

### 1.2. tronko-assign
**Purpose:** Documents how query sequences are assigned taxonomic classifications.

**Key files to analyze:**
- `tronko-assign/tronko-assign.c` - Main program logic
- `tronko-assign/assignment.c` - Core assignment logic
- `tronko-assign/alignment.c` and `alignment_scoring.c` - Alignment scoring
- `tronko-assign/placement.c` - Sequence placement on the tree

**Documentation focus:**
- Read processing workflow (both single and paired-end)
- BWA alignment to leaf nodes
- WFA/Needleman-Wunsch alignment details
- Scoring mechanisms (including when reads are discarded)
- LCA calculation methodology

### 1.3. External Libraries
**Purpose:** Document how Tronko integrates with external libraries.

**Key components:**
- BWA (Burrows-Wheeler Aligner)
- WFA2 (Wavefront Alignment Algorithm)
- RAxML/FastTree (for phylogenetic tree building)

## 2. Documentation Approach

### 2.1. Component Diagrams
- Create UML-style diagrams showing relationships between major modules
- Map data flow from input files through processing steps to outputs
- Highlight decision points for read filtering/discarding

### 2.2. Process Documentation
- Detailed flowcharts for database building process
- Step-by-step pipeline diagrams for taxonomic assignment
- Pseudocode representation of critical algorithms

### 2.3. API/Function Documentation
- Document key function signatures, inputs, and outputs
- Identify and explain important data structures
- Map memory allocation and deallocation patterns

### 2.4. Scoring and Filtering Deep Dive
- Detailed explanation of scoring formulas and thresholds
- Documentation of cases where reads are discarded
- Analysis of the LCA cutoff parameter and its effects

## 3. Implementation Plan

### 3.1. Code Analysis Phase
1. Systematic source code review focusing on main components
2. Trace execution flow through key functions
3. Identify critical data structures and their transformations
4. Document scoring calculations and cutoff thresholds

### 3.2. Documentation Development
1. Create component-level markdown files with core descriptions
2. Develop diagrams using Mermaid or PlantUML notation
3. Document core algorithms with annotated pseudo-code
4. Create cross-referenced documentation linking related components

### 3.3. Validation Phase
1. Verify documentation against the code implementation
2. Test documentation by following assignment processes end-to-end
3. Ensure all filtering/discarding conditions are documented
4. Cross-check the scoring mechanisms against publications

## 4. Documentation Structure

```
docs/
├── overview/
│   ├── architecture.md
│   └── workflow.md
├── tronko-build/
│   ├── database_format.md
│   ├── tree_partitioning.md
│   └── likelihood_calculation.md
├── tronko-assign/
│   ├── assignment_algorithm.md
│   ├── alignment_methods.md
│   ├── scoring_system.md
│   └── filtering_criteria.md
└── examples/
    ├── single_tree_workflow.md
    └── multiple_tree_workflow.md
```

## 5. Rationale

This documentation approach focuses specifically on transparency in the algorithmic processes that drive Tronko's taxonomic assignment. By documenting the information flow and decision points, we will provide researchers with:

1. **Algorithmic transparency** - Understanding exactly how assignments are made
2. **Scoring clarity** - Insight into how reads are scored and why some may be discarded
3. **Configuration guidance** - Better understanding of parameter effects (e.g., LCA cutoff)
4. **Implementation details** - For those who need to modify or extend the code

The emphasis on information flow and scoring addresses the specific need to understand how taxonomic assignments are generated, which is crucial for researchers evaluating results or comparing against other methods.