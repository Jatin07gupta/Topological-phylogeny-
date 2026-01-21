# MLVMST C++ Input Format Explanation

## Overview

The MLVMST C++ executable (`main_mlvmst.cpp`) takes a **text file** as input that contains:
1. A list of vertices (species/taxa names)
2. A list of weighted edges (pairwise distances between vertices)

This input file represents a **complete weighted graph** where:
- **Vertices** = Species/taxa (e.g., "Homo_sapiens", "Pan_troglodytes")
- **Edges** = Pairwise distances between species (e.g., patristic distances from a phylogenetic tree)
- **Weights** = Distance values (floating-point numbers)

---

## Input File Format

The input file is a **plain text file** with the following structure:

```
<number_of_vertices>
<number_of_edges>
<vertex_name_1>
<vertex_name_2>
...
<vertex_name_n>
<vertex_u_1> <vertex_v_1> <weight_1>
<vertex_u_2> <vertex_v_2> <weight_2>
...
<vertex_u_m> <vertex_v_m> <weight_m>
```

### **Format Details:**

1. **Line 1**: Integer - Number of vertices (n)
2. **Line 2**: Integer - Number of edges (m)
3. **Lines 3 to n+2**: Vertex names (one per line, n lines total)
4. **Lines n+3 onwards**: Edges (format: `u v weight`, one per line, m lines total)

### **Edge Format:**
- Each edge is on a single line
- Format: `<vertex_u> <vertex_v> <weight>`
- Vertices are separated by spaces
- Weight is a floating-point number
- Edges are **undirected** (u-v is the same as v-u)

---

## Example Input File

Here's a simple example with 4 vertices and 6 edges:

```
4
6
Homo_sapiens
Pan_troglodytes
Gorilla_gorilla
Pongo_pygmaeus
Homo_sapiens Pan_troglodytes 0.0064
Homo_sapiens Gorilla_gorilla 0.0086
Homo_sapiens Pongo_pygmaeus 0.0152
Pan_troglodytes Gorilla_gorilla 0.0086
Pan_troglodytes Pongo_pygmaeus 0.0152
Gorilla_gorilla Pongo_pygmaeus 0.0152
```

**Explanation:**
- 4 vertices: Homo_sapiens, Pan_troglodytes, Gorilla_gorilla, Pongo_pygmaeus
- 6 edges: All pairwise combinations (4 choose 2 = 6)
- Weights: Patristic distances between species

---

## How the Input is Processed

### **Step 1: File Loading** (`load_test_data()`)

The `load_test_data()` function in `src/data_loader.cpp` reads the file:

```cpp
bool load_test_data(const std::string& filename, VertexList& vertices, EdgeList& edges) {
    // 1. Open file
    std::ifstream file(filename);
    
    // 2. Read number of vertices
    int n_vertices;
    file >> n_vertices;
    
    // 3. Read number of edges
    int n_edges;
    file >> n_edges;
    
    // 4. Read vertex names (one per line)
    for (int i = 0; i < n_vertices; ++i) {
        std::string vertex;
        std::getline(file, vertex);
        vertices.push_back(vertex);
    }
    
    // 5. Read edges (format: "u v weight")
    for (int i = 0; i < n_edges; ++i) {
        std::string u, v;
        double weight;
        file >> u >> v >> weight;
        edges.push_back({u, v, weight});
    }
    
    return true;
}
```

### **Step 2: Data Structures Created**

After loading, the following data structures are created:

**1. `VertexList`** (std::vector<std::string>):
```cpp
["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla", "Pongo_pygmaeus"]
```

**2. `EdgeList`** (std::vector<std::tuple<std::string, std::string, double>>):
```cpp
[
    ("Homo_sapiens", "Pan_troglodytes", 0.0064),
    ("Homo_sapiens", "Gorilla_gorilla", 0.0086),
    ("Homo_sapiens", "Pongo_pygmaeus", 0.0152),
    ("Pan_troglodytes", "Gorilla_gorilla", 0.0086),
    ("Pan_troglodytes", "Pongo_pygmaeus", 0.0152),
    ("Gorilla_gorilla", "Pongo_pygmaeus", 0.0152)
]
```

---

## Complete Pipeline Flow

Here's how the input flows through the MLVMST pipeline:

```
Input File (text)
    ↓
load_test_data()
    ↓
VertexList + EdgeList
    ↓
┌─────────────────────────────────────┐
│ Algorithm 2: construct_FC_and_GU()  │
│ Input: vertices, edges              │
│ Output: F_C, G_U_edges              │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│ Algorithm 3: build_mlvmst()         │
│ Input: vertices, edges, F_C, G_U    │
│ Output: M* edges, order, leaf_count │
└─────────────────────────────────────┘
    ↓
MLVMSTResult (M* spanning tree)
```

---

## What MLVMST Does with the Input

### **1. Algorithm 2 (F_C and G_U Construction)**

**Input:**
- `vertices`: List of all vertex names
- `edges`: List of all weighted edges

**Process:**
- Sorts edges by weight
- Runs Kruskal's algorithm to find MST components
- Builds Common Laminar Family (F_C)
- Builds MST-Union Graph (G_U)

**Output:**
- `F_C`: List of vertex sets (laminar family)
- `G_U_edges`: Set of edges in MST-union graph

### **2. Algorithm 3 (MLVMST Construction)**

**Input:**
- `vertices`: List of all vertex names
- `edges`: List of all weighted edges (same as Algorithm 2)
- `F_C`: Common Laminar Family from Algorithm 2
- `G_U_edges`: MST-Union Graph edges from Algorithm 2
- `ordering`: "increasing" or "decreasing" (default: "increasing")

**Process:**
1. **Compute delta_max**: For each vertex, compute how many F_C sets cover its neighbors
2. **Build total order**: Sort vertices by delta_max (increasing or decreasing)
3. **Build VMST**: Run Kruskal's algorithm with tie-breaking by vertex order
4. **Count leaves**: Count vertices with degree 1 in the resulting tree

**Output:**
- `M_star_edges`: Edges of the minimum leaves spanning tree
- `order`: Vertex ordering used
- `delta_max`: Delta_max scores for each vertex
- `leaf_count`: Number of leaves in M*

---

## Real-World Example

For the **Primates Species** dataset (449 nodes):

**Input File Structure:**
```
449
100576
Perodicticus_potto
Arctocebus_calabarensis
Arctocebus_aureus
...
[449 vertex names total]
...
Perodicticus_potto Arctocebus_calabarensis 22.09774
Perodicticus_potto Arctocebus_aureus 22.09774
...
[100,576 edge entries total]
```

**What Happens:**
1. File is loaded → 449 vertices, 100,576 edges
2. Algorithm 2 runs → F_C (860 components), G_U (11,595 edges)
3. Algorithm 3 runs → M* (448 edges), 173 leaves
4. Result: Minimum leaves spanning tree with 173 leaves

---

## Key Points

### **1. Complete Graph**
- The input represents a **complete weighted graph**
- All pairwise distances are provided (n choose 2 = n(n-1)/2 edges)
- This is typically generated from a phylogenetic tree (patristic distances)

### **2. Edge Weights**
- Weights are **floating-point numbers**
- Represent distances (e.g., patristic distances, evolutionary distances)
- Higher precision is preserved (`.17g` format in comparison scripts)

### **3. Vertex Names**
- Vertex names are **strings**
- Can contain underscores, hyphens, etc.
- Must be unique (no duplicates)

### **4. Undirected Edges**
- Edges are **undirected** (u-v = v-u)
- The algorithm handles canonicalization internally
- Only one direction needs to be provided in the file

---

## Command-Line Usage

```bash
# Basic usage (default: increasing ordering)
./build/bin/mlvmst_main data/primates_species_test.txt

# Explicit ordering
./build/bin/mlvmst_main data/primates_species_test.txt increasing
./build/bin/mlvmst_main data/primates_species_test.txt decreasing
```

**Arguments:**
1. **Required**: Input file path (text file with vertices and edges)
2. **Optional**: Ordering ("increasing" or "decreasing", default: "increasing")

---

## How Input Files are Generated

In practice, input files are generated from phylogenetic trees:

1. **Load Newick tree** (`.nwk` file)
2. **Compute patristic distances** (pairwise distances between all leaves)
3. **Generate input file** with:
   - Number of leaves
   - All pairwise edges with distances

**Example Python code** (from comparison scripts):
```python
# Load tree
tree = Tree("primates_species.nwk")

# Get leaf names
leaves = [leaf.name for leaf in tree.get_leaves()]

# Compute patristic distances
edges = compute_patristic_distances_optimized("primates_species.nwk", leaves)

# Write to file
with open("data/primates_species_test.txt", 'w') as f:
    f.write(f"{len(leaves)}\n")
    f.write(f"{len(edges)}\n")
    for leaf in leaves:
        f.write(f"{leaf}\n")
    for u, v, dist in edges:
        f.write(f"{u} {v} {dist:.17g}\n")
```

---

## Summary

**Input to MLVMST C++:**
- **Format**: Plain text file
- **Content**: 
  - Number of vertices
  - Number of edges
  - Vertex names (one per line)
  - Weighted edges (format: `u v weight`, one per line)
- **Represents**: Complete weighted graph (all pairwise distances)
- **Typical Source**: Patristic distances from phylogenetic tree

**What MLVMST Does:**
1. Loads vertices and edges from file
2. Runs Algorithm 2 to compute F_C and G_U
3. Runs Algorithm 3 to compute minimum leaves spanning tree (M*)
4. Outputs M* edges, vertex order, and leaf count

**Output:**
- M* edges: The minimum leaves spanning tree
- Order: Vertex ordering used for construction
- Leaf count: Number of leaves in M*

---

**See Also:**
- `src/data_loader.cpp` - File loading implementation
- `src/main_mlvmst.cpp` - Main executable
- `compare_mlvmst_python_cpp.py` - Example of input file generation

