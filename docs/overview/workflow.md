# Tronko Workflow

This document describes the end-to-end workflow of Tronko, from database building to taxonomic assignment, with a focus on data flow and processing steps.

## Overall Workflow

```mermaid
flowchart TD
    %% Database Building Phase
    subgraph "Phase 1: Database Building"
        MSA[Multiple Sequence Alignment] --> TB[tronko-build]
        TAX[Taxonomy File] --> TB
        TREE[Phylogenetic Tree] --> TB
        TB --> DB[Reference Database]
    end
    
    %% Assignment Phase
    subgraph "Phase 2: Read Assignment"
        READS[Query Reads] --> TA[tronko-assign]
        DB --> TA
        REFFASTA[Reference Sequences FASTA] --> TA
        TA --> RESULTS[Taxonomic Assignments]
    end
```

## tronko-build Workflow

```mermaid
flowchart TD
    %% Input Processing
    MSA[Multiple Sequence Alignment] --> RF[Read FASTA Sequences]
    TREE[Phylogenetic Tree] --> RR[Parse Newick Tree]
    TAX[Taxonomy File] --> PT[Process Taxonomic Information]
    
    %% Core Processing
    RF --> |For each MSA| TP[Tree Partitioning]
    RR --> TP
    
    %% Decision Point
    TP --> DEC{Partition Tree?}
    DEC -->|Yes| SPLIT[Split Tree]
    SPLIT --> REBUILD[Rebuild Sub-trees]
    REBUILD --> CALC
    DEC -->|No| CALC[Calculate Node Likelihoods]
    
    %% Final Processing
    PT --> TAX_MAP[Map Taxonomy to Tree]
    CALC --> LIKE[Store Likelihood Values]
    TAX_MAP --> OUT[Write Reference Database]
    LIKE --> OUT
```

## tronko-assign Workflow

```mermaid
flowchart TD
    %% Input Processing
    READS[Query Reads] --> PROC[Process Reads]
    PROC -->|Single-end| SE[Single Read Processing]
    PROC -->|Paired-end| PE[Paired Read Processing]
    
    %% BWA Alignment
    SE --> BWA_S[BWA Alignment]
    PE --> BWA_P[BWA Alignment]
    
    %% Alignment and Assignment
    BWA_S --> |Choose alignment method| ALG_S{Alignment Method}
    BWA_P --> |Choose alignment method| ALG_P{Alignment Method}
    
    ALG_S -->|WFA| WFA_S[Wavefront Alignment]
    ALG_S -->|NW| NW_S[Needleman-Wunsch]
    ALG_P -->|WFA| WFA_P[Wavefront Alignment]
    ALG_P -->|NW| NW_P[Needleman-Wunsch]
    
    WFA_S --> SCORE_S[Score Alignment]
    NW_S --> SCORE_S
    WFA_P --> SCORE_P[Score Alignment]
    NW_P --> SCORE_P
    
    %% Placement and Scoring
    SCORE_S --> PLACE_S[Place on Phylogeny]
    SCORE_P --> PLACE_P[Place on Phylogeny]
    
    PLACE_S --> LCA_S[Calculate LCA]
    PLACE_P --> LCA_P[Calculate LCA]
    
    %% Filtering
    LCA_S --> FILTER_S{Score > Cutoff?}
    LCA_P --> FILTER_P{Score > Cutoff?}
    
    FILTER_S -->|Yes| ASSIGN_S[Assign Taxonomy]
    FILTER_S -->|No| DISCARD_S[Discard Read]
    FILTER_P -->|Yes| ASSIGN_P[Assign Taxonomy]
    FILTER_P -->|No| DISCARD_P[Discard Read]
    
    %% Output
    ASSIGN_S --> OUT[Output Results]
    ASSIGN_P --> OUT
```

## Information Flow Details

### 1. Database Building

1. **Input Files Processing**:
   - Multiple sequence alignment (MSA) in FASTA format is parsed
   - Phylogenetic tree in Newick format is processed
   - Taxonomy information is mapped to tree nodes

2. **Tree Partitioning** (optional):
   - Sum-of-pairs score or minimum leaf node count used to determine partitioning
   - Trees can be split into multiple sub-trees for more accurate alignment

3. **Likelihood Calculation**:
   - Fractional likelihoods calculated for all tree nodes
   - This is a key differentiator from traditional LCA methods that only use leaf nodes

4. **Output Generation**:
   - Reference database file containing tree structure, node likelihoods, and taxonomic information
   - Format optimized for efficient loading by tronko-assign

### 2. Read Assignment

1. **Read Processing**:
   - Handles both single-end and paired-end reads
   - Optional reverse complementing based on read orientation

2. **Initial Alignment**:
   - BWA used for fast mapping to leaf nodes in the reference tree
   - Provides initial positioning for more detailed alignment

3. **Detailed Alignment**:
   - Two options available:
     - Wavefront Alignment Algorithm (WFA2) - Faster, default option
     - Needleman-Wunsch - Traditional alignment algorithm

4. **Scoring and Placement**:
   - Alignment scores calculated based on matches/mismatches
   - Reads placed on the phylogenetic tree based on alignment scores

5. **LCA Calculation**:
   - Lowest Common Ancestor determined based on score cutoff
   - Critical filtering point: reads below threshold are discarded

6. **Output Generation**:
   - Tab-delimited results file with taxonomic assignments and scores
   - Includes mismatch information and tree/node placement details

## Read Filtering Conditions

Reads may be filtered out (discarded) at several points in the workflow:

1. **During BWA Alignment**:
   - Reads with no hits to reference sequences are discarded
   - Threshold determines minimum mapping quality

2. **During Score Evaluation**:
   - LCA cutoff parameter (-c) determines the score threshold
   - Reads with scores below the threshold are assigned to higher taxonomic levels or discarded

3. **During Taxonomic Assignment**:
   - Reads that cannot be confidently assigned to any taxonomic level are discarded

## Key Parameters Affecting Workflow

- **LCA Cutoff** (-c): Controls the stringency of taxonomic assignment
- **Alignment Method** (-w): Selects between WFA (default) or Needleman-Wunsch
- **Thread Count** (-C): Number of parallel processing threads
- **Score Constant** (-u): Affects likelihood calculations (default: 0.01)