# Tronko Reference Database Format

This document describes the structure and format of the reference database created by `tronko-build` and used by `tronko-assign` for taxonomic classification.

## Overview

The reference database (`reference_tree.txt`) is a text file containing all necessary information for taxonomic assignment, including:

1. Tree structure information
2. Fractional likelihood values for each node
3. Taxonomic information for leaf nodes
4. MSA positions and other metadata

## File Structure

The file is organized into sections with a specific format:

```
<Number of trees>
<For each tree>
    <Number of nodes in tree>
    <For each node>
        <Node data>
        <Likelihood data>
        <Taxonomic data (if leaf node)>
    <Tree structure data>
```

## Node Entry Format

Each node in the tree is represented with the following fields:

```
<Node ID> <Parent ID> <Number of Children> [Child IDs...] <Leaf Flag> <Taxonomic Path (if leaf)>
```

Where:
- **Node ID**: Unique identifier for the node
- **Parent ID**: ID of the parent node (-1 for root)
- **Number of Children**: Count of child nodes
- **Child IDs**: List of child node IDs
- **Leaf Flag**: 1 if leaf node, 0 otherwise
- **Taxonomic Path**: Semicolon-delimited taxonomy (for leaf nodes only)

## Likelihood Data Format

Each node contains likelihood data that enables taxonomic assignment:

```
<Fractional likelihood value> <Cumulative likelihood> <Position data>
```

Where:
- **Fractional likelihood**: Probability value for this node
- **Cumulative likelihood**: Accumulated probability for this subtree
- **Position data**: Information about MSA positions relevant to this node

## File Generation Process

1. `tronko-build` parses the input phylogenetic tree and constructs an internal tree representation
2. Node relationships (parent-child) are established
3. Likelihood values are calculated for all nodes based on the MSA
4. Taxonomic information is attached to leaf nodes
5. The complete tree is serialized to the output file

## Example (Simplified)

```
1                       # Number of trees
1000                    # Number of nodes in tree 1
0 -1 2 1 2 0            # Root node (ID 0) with 2 children (1,2)
1 0 2 3 4 0             # Internal node (ID 1) with 2 children (3,4)
2 0 2 5 6 0             # Internal node (ID 2) with 2 children (5,6)
3 1 0 1 Eukaryota;Chordata;Aves;Passeriformes;Paridae;Parus;Parus major   # Leaf node with taxonomy
...
```

## Accessing Database Content

To read and parse the reference database:

1. Read the number of trees
2. For each tree, read the number of nodes
3. For each node, parse its data and store in appropriate data structures
4. Construct the tree hierarchy based on parent-child relationships
5. Associate likelihood values with each node

## Database Size Considerations

The database file size is determined by:
- Number of reference sequences
- Tree topology (balanced vs. unbalanced)
- Number of partitions (if using multiple trees)

A typical reference database for a medium-sized dataset (1,000-10,000 sequences) can range from a few MB to hundreds of MB.

## Compression

The reference database can be gzipped to reduce storage requirements. `tronko-assign` can read gzipped databases directly using the zlib library.

## Custom Database Generation

For specialized applications, researchers can customize the database by:
1. Selecting specific reference sequences
2. Using custom taxonomies
3. Adjusting likelihood calculations through partitioning parameters