# Tronko Documentation

This directory contains comprehensive documentation for the Tronko codebase, a rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets.

## Documentation Structure

```
docs/
├── overview/
│   ├── architecture.md - System architecture and component relationships
│   └── workflow.md - End-to-end workflow and information flow
├── tronko-build/
│   ├── database_format.md - Reference database structure and format
│   ├── likelihood_calculation.md - How fractional likelihoods are calculated
│   └── tree_partitioning.md - Tree partitioning algorithms and methods
├── tronko-assign/
│   ├── alignment_methods.md - Sequence alignment algorithms (WFA2 and Needleman-Wunsch)
│   ├── assignment_algorithm.md - Core taxonomic assignment process
│   ├── filtering_criteria.md - How reads are filtered and when they are discarded
│   └── scoring_system.md - Scoring mechanisms for alignment and taxonomic placement
├── examples/
│   ├── single_tree_workflow.md - Complete workflow using a single reference tree
│   └── multiple_tree_workflow.md - Complete workflow using multiple reference trees
└── README.md - This file
```

## Key Topics

### Information Flow

The documentation focuses on how information flows through the system:

1. From input files (MSA, tree, taxonomy) to the reference database
2. From the reference database to taxonomic assignments
3. Through various filtering stages where reads may be discarded

### Assignment Generation

The assignment process is documented in detail:

1. Initial BWA alignment to leaf nodes
2. Detailed alignment using WFA2 or Needleman-Wunsch
3. Scoring and tree placement
4. LCA determination based on score cutoffs
5. Taxonomic assignment output generation

### Scoring Mechanisms

Scoring is critical to taxonomic assignment:

1. Alignment scores based on matches, mismatches, and gaps
2. Tree placement scores incorporating fractional likelihoods
3. LCA cutoff thresholds for determining taxonomic resolution

### Read Filtering

The documentation explains when and why reads are discarded:

1. During BWA alignment (no hits or low mapping quality)
2. During alignment scoring (poor alignment quality)
3. During tree placement (ambiguous placement)
4. During LCA determination (below score threshold)

## How to Use This Documentation

- **New Users**: Start with `overview/workflow.md` for a high-level understanding
- **Users Running Analyses**: Follow example workflows in the `examples/` directory
- **Parameter Optimization**: See `tronko-assign/filtering_criteria.md` for tuning advice
- **Developers**: Refer to architecture and algorithm documentation for implementation details

## Additional Resources

- **Main README**: See the root `/README.md` file for installation and basic usage
- **Example Datasets**: Find example data in `example_datasets/` directory
- **Publication**: Pipes L, and Nielsen R (2022) A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. bioRXiv. https://www.biorxiv.org/content/10.1101/2022.12.06.519402v1