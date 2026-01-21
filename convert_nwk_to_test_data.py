#!/usr/bin/env python3
"""
Convert Newick tree file to C++ test data format.

Usage:
    python convert_nwk_to_test_data.py <input.nwk> <output.txt>

The output format is:
    Line 1: number of vertices
    Line 2: number of edges
    Next n lines: vertex names (one per line)
    Next m lines: edges in format "u v weight"
"""

import sys
import os
from ete3 import Tree
from python_preprocessing.python_patristic_optimized import compute_patristic_distances_optimized


def convert_newick_to_test_data(newick_file: str, output_file: str):
    """
    Convert Newick tree to C++ test data format.
    
    Args:
        newick_file: Path to input Newick file
        output_file: Path to output test data file
    """
    print(f"Loading tree from: {newick_file}")
    
    # Try different format combinations
    tree = None
    for fmt in [0, 1, 2, 3, 5, 8, 9]:
        for quoted in [False, True]:
            try:
                tree = Tree(newick_file, format=fmt, quoted_node_names=quoted)
                print(f"  ✓ Loaded with format={fmt}, quoted={quoted}")
                break
            except:
                continue
        if tree is not None:
            break
    
    if tree is None:
        raise ValueError(f"Could not load tree from {newick_file}")
    
    # Get all leaf names
    leaves = [leaf.name for leaf in tree.get_leaves()]
    n = len(leaves)
    print(f"  Tree has {n} leaves")
    
    # Compute pairwise patristic distances
    print(f"  Computing {n}x{n} distance matrix...")
    edges_list = compute_patristic_distances_optimized(newick_file, leaves, verbose=False)
    
    print(f"  Computed {len(edges_list)} edges")
    
    # Write to output file
    print(f"Writing to: {output_file}")
    with open(output_file, 'w') as f:
        # Write number of vertices
        f.write(f"{n}\n")
        
        # Write number of edges
        f.write(f"{len(edges_list)}\n")
        
        # Write vertex names (one per line)
        for leaf in leaves:
            f.write(f"{leaf}\n")
        
        # Write edges (format: "u v weight")
        for leaf1, leaf2, dist in edges_list:
            f.write(f"{leaf1} {leaf2} {dist}\n")
    
    print(f"  ✓ Successfully wrote {n} vertices and {len(edges_list)} edges")
    print(f"  ✓ Output file: {output_file}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_nwk_to_test_data.py <input.nwk> <output.txt>")
        print("\nExample:")
        print("  python convert_nwk_to_test_data.py ../primates_genus.nwk data/primates_genus_test.txt")
        sys.exit(1)
    
    newick_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Handle relative paths
    if not os.path.exists(newick_file):
        # Try parent directory
        parent_path = os.path.join('..', newick_file)
        if os.path.exists(parent_path):
            newick_file = parent_path
        else:
            root_path = os.path.join('../..', newick_file)
            if os.path.exists(root_path):
                newick_file = root_path
    
    if not os.path.exists(newick_file):
        print(f"Error: File not found: {newick_file}")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        convert_newick_to_test_data(newick_file, output_file)
        print("\n✅ Conversion complete!")
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()





