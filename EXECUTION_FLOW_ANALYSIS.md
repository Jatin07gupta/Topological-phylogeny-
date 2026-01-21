# Complete Execution Flow Analysis: `main_clnj_boruvka.cpp`

## 📋 Table of Contents
1. [Input Format](#input-format)
2. [Execution Flow Overview](#execution-flow-overview)
3. [Step-by-Step Execution](#step-by-step-execution)
4. [Function Call Tree](#function-call-tree)
5. [Function Source Files](#function-source-files)
6. [Data Flow](#data-flow)

---

## 📥 Input Format

### Command Line Arguments
```bash
./clnj_boruvka_main <test_data_file> [ordering] [num_threads]
```

**Arguments:**
- **`test_data_file`** (required): Path to input data file
- **`ordering`** (optional, default: "increasing"): "increasing" or "decreasing"
- **`num_threads`** (optional, default: 4): Number of threads for parallel execution (0 = use OpenMP default)

### Input File Format (.txt file)

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

**Example:**
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

**Format Details:**
- **Line 1**: Integer - Number of vertices (n)
- **Line 2**: Integer - Number of edges (m)
- **Lines 3 to n+2**: Vertex names (one per line, n lines total)
- **Lines n+3 onwards**: Edges (format: `u v weight`, one per line, m lines total)
  - Each edge: `<vertex_u> <vertex_v> <weight>` (space-separated)
  - Edges are **undirected** (u-v = v-u)
  - Weight is a floating-point number

---

## 🔄 Execution Flow Overview

```
[Command Line Input]
        ↓
[Parse Arguments]
        ↓
[Load Test Data] → data_loader.cpp
        ↓
[Step 1: Algorithm 2] → algorithm2.cpp
  → construct_FC_and_GU()
        ↓
[Step 2: Delta_max & Order] → mlvmst.cpp
  → build_adjacency()
  → compute_delta_max()
  → build_total_order()
        ↓
[Step 3: MLVMST Borůvka Parallel] → mlvmst_boruvka.cpp
  → build_vmst_from_order_hybrid_parallel()
  → count_leaves()
        ↓
[Step 4: Prepare CLNJ Input]
  → build_distance_matrix() [LOCAL FUNCTION]
  → build_mlvmst_adjacency() [LOCAL FUNCTION]
        ↓
[Step 5: CLNJ] → clnj_fresh.cpp
  → CLNJ_clean()
    → neighbor_joining_local()
      → compute_distance_with_hidden()
        ↓
[Print Results]
  → print_clnj_result() [LOCAL FUNCTION]
```

---

## 📝 Step-by-Step Execution

### **Phase 0: Initialization** (lines 172-230)

#### **Function: `main(int argc, char* argv[])`**
- **Location**: `main_clnj_boruvka.cpp` (lines 172-344)
- **Purpose**: Entry point of the program

**Execution Steps:**
1. **Parse Command Line Arguments** (lines 173-203)
   - Extract filename from `argv[1]`
   - Extract ordering from `argv[2]` (default: "increasing")
   - Extract num_threads from `argv[3]` (default: 4)
   - Validate arguments

2. **OpenMP Setup** (lines 205-213)
   - If OpenMP available: Set thread count
   - If OpenMP not available: Sequential execution

3. **Load Input Data** (lines 215-230)
   - **Function Called**: `load_test_data(filename, vertices, edges)`
   - **Source**: `data_loader.cpp` (line 9)
   - **Input**: Filename (string)
   - **Output**: `VertexList vertices`, `EdgeList edges`
   - **What it does**:
     - Opens the input file
     - Reads number of vertices and edges
     - Reads vertex names into `vertices`
     - Reads edges (u, v, weight) into `edges`

---

### **Phase 1: Algorithm 2 - Construct F_C and G_U** (lines 232-244)

#### **Function: `construct_FC_and_GU(vertices, edges)`**
- **Location**: Called from `main_clnj_boruvka.cpp` (line 238)
- **Source**: `algorithm2.cpp` (declared in `algorithm2.h`)
- **Purpose**: Builds the Common Laminar Family (F_C) and MST-Union Graph (G_U)

**Input:**
- `VertexList vertices`: List of vertex names
- `EdgeList edges`: List of weighted edges (u, v, weight)

**Output:**
- `FC_GU_Result alg2_result`:
  - `F_C`: Common laminar family (list of vertex sets)
  - `G_U_edges`: MST-union graph edges (set of vertex pairs)

**What it does:**
1. Sorts edges by weight
2. Runs Kruskal's algorithm with DSU (Disjoint Set Union)
3. At each weight tier, adds components to F_C
4. Adds edges between different components to G_U

**Timing**: Measured in milliseconds

---

### **Phase 2: Compute Delta_max and Vertex Order** (lines 246-266)

#### **Step 2.1: Build G_U Adjacency List**
- **Function**: `build_adjacency(alg2_result.G_U_edges)`
- **Location**: Called from line 254
- **Source**: `mlvmst.cpp` (declared in `mlvmst.h` line 33)
- **Purpose**: Converts edge set to adjacency list representation
- **Input**: `std::set<std::pair<Vertex, Vertex>>` (edge set)
- **Output**: `std::unordered_map<Vertex, VertexSet>` (adjacency list)

#### **Step 2.2: Compute Delta_max**
- **Function**: `compute_delta_max(vertices, alg2_result.F_C, G_U_adj)`
- **Location**: Called from line 257
- **Source**: `mlvmst.cpp` (declared in `mlvmst.h` line 54)
- **Purpose**: Computes delta_max score for each vertex
- **Algorithm**:
  - For each vertex i:
    - Get neighbors N_i from G_U
    - For each component C in F_C (ordered by decreasing size):
      - If C contains neighbors of i but not i itself:
        - Increment delta_max(i)
        - Remove those neighbors from N_i
- **Input**:
  - `VertexList vertices`
  - `std::vector<VertexSet> F_C`
  - `std::unordered_map<Vertex, VertexSet> G_U_adj`
- **Output**: `std::map<Vertex, int> delta_max`

#### **Step 2.3: Build Vertex Order**
- **Function**: `build_total_order(vertices, delta_max, ordering)`
- **Location**: Called from line 260
- **Source**: `mlvmst.cpp` (declared in `mlvmst.h` line 72)
- **Purpose**: Sorts vertices by delta_max scores
- **Input**:
  - `VertexList vertices`
  - `std::map<Vertex, int> delta_max`
  - `std::string ordering` ("increasing" or "decreasing")
- **Output**: `VertexList vertex_order` (sorted vertices)

**Timing**: Measured in milliseconds

---

### **Phase 3: Build MLVMST using Borůvka Parallel** (lines 268-285)

#### **Step 3.1: Build VMST**
- **Function**: `build_vmst_from_order_hybrid_parallel(vertices, edges, vertex_order, num_threads)`
- **Location**: Called from line 274
- **Source**: `mlvmst_boruvka.cpp` (declared in `mlvmst_boruvka.h` line 58)
- **Purpose**: Constructs the Minimum Leaves VMST using parallel Borůvka algorithm
- **Algorithm**:
  - Borůvka-style phased approach with OpenMP parallelization
  - Each phase: Find best edge per component, merge components
  - Maintains same result as Kruskal's algorithm but faster
- **Input**:
  - `VertexList vertices`
  - `EdgeList edges`
  - `VertexList vertex_order`
  - `int num_threads`
- **Output**: `std::set<std::pair<Vertex, Vertex>> mlvmst_edges`

#### **Step 3.2: Count Leaves**
- **Function**: `count_leaves(mlvmst_edges, vertices)`
- **Location**: Called from line 281
- **Source**: `mlvmst.cpp` (declared in `mlvmst.h` line 103)
- **Purpose**: Counts vertices with degree 1 (leaves) in the MLVMST
- **Input**: 
  - `std::set<std::pair<Vertex, Vertex>> mlvmst_edges`
  - `VertexList vertices`
- **Output**: `int leaf_count`

**Timing**: Measured in milliseconds

---

### **Phase 4: Prepare CLNJ Input** (lines 287-311)

#### **Step 4.1: Build Distance Matrix**
- **Function**: `build_distance_matrix(vertices, edges)`
- **Location**: Defined in `main_clnj_boruvka.cpp` (lines 39-70)
- **Purpose**: Creates a distance matrix from edge list
- **Input**: `VertexList vertices`, `EdgeList edges`
- **Output**: `std::vector<std::vector<double>> distance_matrix_vec`
- **Algorithm**:
  - Creates m×m matrix
  - Fills matrix from edges (symmetric: D[i][j] = D[j][i])

#### **Step 4.2: Build MLVMST Adjacency Matrix**
- **Function**: `build_mlvmst_adjacency(vertices, mlvmst_edges)`
- **Location**: Defined in `main_clnj_boruvka.cpp` (lines 75-105)
- **Purpose**: Creates an adjacency matrix from MLVMST edges
- **Input**: `VertexList vertices`, `std::set<std::pair<Vertex, Vertex>> mlvmst_edges`
- **Output**: `std::vector<std::vector<int>> mst_adjacency_vec`
- **Algorithm**:
  - Creates m×m matrix (0/1 values)
  - Sets A[i][j] = 1 if edge exists between i and j

#### **Step 4.3: Convert to Eigen Matrices**
- **Location**: Lines 297-304
- **Purpose**: Convert std::vector matrices to Eigen::MatrixXd/MatrixXi
- **Creates**:
  - `Eigen::MatrixXd distance_matrix(m, m)`
  - `Eigen::MatrixXi mst_adjacency(m, m)`

**Timing**: Measured in milliseconds

---

### **Phase 5: Run CLNJ** (lines 313-326)

#### **Function: `CLNJ_clean(distance_matrix, mst_adjacency, threshold, true)`**
- **Location**: Called from line 322
- **Source**: `clnj_fresh.cpp` (declared in `clnj.h` line 114)
- **Purpose**: Reconstructs latent phylogenetic tree using Complete Latent Neighbor Joining
- **Input**:
  - `Eigen::MatrixXd distance_matrix`: m×m distance matrix
  - `Eigen::MatrixXi mst_adjacency`: m×m MLVMST adjacency matrix
  - `double threshold`: Threshold for hidden node creation (= -log(0.9))
  - `bool verbose`: Print progress messages (true)

**Algorithm (High Level):**
1. Initialize adjacency from MLVMST
2. Find internal nodes (degree > 1)
3. Sort internal nodes by degree (descending)
4. For each internal node:
   - Get its neighborhood S_v
   - Run local Neighbor Joining on S_v
   - Remove old edges within S_v
   - Add new edges from NJ result
   - Create hidden nodes if branch lengths exceed threshold

**Internal Functions Called:**
- `neighbor_joining_local()` - Local NJ algorithm
  - **Source**: `clnj_fresh.cpp` (line 79)
  - **Calls**: `compute_distance_with_hidden()` 
    - **Source**: `clnj_fresh.cpp` (line 23)

**Output**: `CLNJResult clnj_result`
- `adjacency`: Final tree adjacency
- `edge_weights`: Edge weights
- `hidden_info`: Hidden node information

**Timing**: Measured in milliseconds

---

### **Phase 6: Print Results** (lines 328-341)

#### **Function: `print_clnj_result(clnj_result, m)`**
- **Location**: Defined in `main_clnj_boruvka.cpp` (lines 107-170)
- **Purpose**: Pretty-prints CLNJ results
- **Output**:
  - Summary statistics (observed nodes, hidden nodes, total edges)
  - Adjacency list (parseable format)
  - Hidden node information

#### **Timing Summary**
- Prints timing for each phase
- Calculates total execution time

---

## 🌳 Function Call Tree

```
main()
├── load_test_data()                    [data_loader.cpp]
│
├── construct_FC_and_GU()               [algorithm2.cpp]
│   ├── (uses DisjointSetUnion class)
│   └── (uses ComponentTracker class)
│
├── build_adjacency()                   [mlvmst.cpp]
│
├── compute_delta_max()                 [mlvmst.cpp]
│
├── build_total_order()                 [mlvmst.cpp]
│
├── build_vmst_from_order_hybrid_parallel() [mlvmst_boruvka.cpp]
│   └── (uses OpenMP for parallelization)
│
├── count_leaves()                      [mlvmst.cpp]
│
├── build_distance_matrix()             [LOCAL - main_clnj_boruvka.cpp]
│
├── build_mlvmst_adjacency()            [LOCAL - main_clnj_boruvka.cpp]
│
├── CLNJ_clean()                        [clnj_fresh.cpp]
│   ├── neighbor_joining_local()        [clnj_fresh.cpp]
│   │   └── compute_distance_with_hidden() [clnj_fresh.cpp]
│   └── (processes internal nodes iteratively)
│
└── print_clnj_result()                 [LOCAL - main_clnj_boruvka.cpp]
```

---

## 📚 Function Source Files

### **Functions Defined in `main_clnj_boruvka.cpp` (Local Functions)**

| Function | Lines | Purpose |
|----------|-------|---------|
| `build_distance_matrix()` | 39-70 | Builds distance matrix from edge list |
| `build_mlvmst_adjacency()` | 75-105 | Builds adjacency matrix from MLVMST edges |
| `print_clnj_result()` | 107-170 | Prints CLNJ results in parseable format |
| `main()` | 172-344 | Main entry point |

### **Functions from `data_loader.cpp`**

| Function | Header File | Purpose |
|----------|-------------|---------|
| `load_test_data()` | `data_loader.h` | Loads vertices and edges from text file |

### **Functions from `algorithm2.cpp`**

| Function | Header File | Purpose |
|----------|-------------|---------|
| `construct_FC_and_GU()` | `algorithm2.h` | Algorithm 2: Constructs F_C and G_U |

### **Functions from `mlvmst.cpp`**

| Function | Header File | Purpose |
|----------|-------------|---------|
| `build_adjacency()` | `mlvmst.h` | Converts edge set to adjacency list |
| `compute_delta_max()` | `mlvmst.h` | Computes delta_max scores |
| `build_total_order()` | `mlvmst.h` | Builds vertex ordering from delta_max |
| `count_leaves()` | `mlvmst.h` | Counts leaves in tree |

### **Functions from `mlvmst_boruvka.cpp`**

| Function | Header File | Purpose |
|----------|-------------|---------|
| `build_vmst_from_order_hybrid_parallel()` | `mlvmst_boruvka.h` | Builds VMST using parallel Borůvka |

### **Functions from `clnj_fresh.cpp`**

| Function | Header File | Purpose |
|----------|-------------|---------|
| `CLNJ_clean()` | `clnj.h` | Main CLNJ algorithm |
| `neighbor_joining_local()` | `clnj.h` | Local NJ algorithm |
| `compute_distance_with_hidden()` | `clnj.h` | Computes distance with hidden nodes |

### **Standard Library Functions Used**

- `std::chrono::high_resolution_clock::now()` - Timing
- `std::chrono::duration_cast<std::chrono::milliseconds>()` - Time conversion
- `std::atoi()` - String to integer conversion
- `std::getline()` - File reading
- `std::cout`, `std::cerr` - Output streams

### **OpenMP Functions (if available)**

- `omp_get_max_threads()` - Get maximum threads
- OpenMP parallelization in `build_vmst_from_order_hybrid_parallel()`

---

## 🔀 Data Flow

### **Input Data Flow**
```
Input File (.txt)
    ↓
load_test_data()
    ↓
VertexList vertices  [std::vector<std::string>]
EdgeList edges       [std::vector<std::tuple<Vertex, Vertex, double>>]
```

### **Algorithm 2 Data Flow**
```
vertices + edges
    ↓
construct_FC_and_GU()
    ↓
FC_GU_Result:
  - F_C: std::vector<VertexSet>
  - G_U_edges: std::set<std::pair<Vertex, Vertex>>
```

### **Delta_max & Order Data Flow**
```
F_C + G_U_edges
    ↓
build_adjacency() → G_U_adj [std::unordered_map<Vertex, VertexSet>]
    ↓
compute_delta_max() → delta_max [std::map<Vertex, int>]
    ↓
build_total_order() → vertex_order [VertexList]
```

### **MLVMST Data Flow**
```
vertices + edges + vertex_order
    ↓
build_vmst_from_order_hybrid_parallel()
    ↓
mlvmst_edges [std::set<std::pair<Vertex, Vertex>>]
    ↓
count_leaves() → leaf_count [int]
```

### **CLNJ Preparation Data Flow**
```
vertices + edges + mlvmst_edges
    ↓
build_distance_matrix() → distance_matrix_vec [std::vector<std::vector<double>>]
build_mlvmst_adjacency() → mst_adjacency_vec [std::vector<std::vector<int>>]
    ↓
Convert to Eigen:
    - distance_matrix [Eigen::MatrixXd]
    - mst_adjacency [Eigen::MatrixXi]
```

### **CLNJ Data Flow**
```
distance_matrix + mst_adjacency + threshold
    ↓
CLNJ_clean()
    ↓ (iterates over internal nodes)
neighbor_joining_local()
    ↓
CLNJResult:
  - adjacency [std::unordered_map<int, std::unordered_set<int>>]
  - edge_weights [std::map<std::pair<int, int>, double>]
  - hidden_info [std::unordered_map<int, HiddenNodeInfo>]
    ↓
print_clnj_result()
```

---

## ✅ Summary

**Complete Pipeline:**
1. **Input**: Text file with vertices and edges
2. **Algorithm 2**: Constructs F_C (laminar family) and G_U (MST-union graph)
3. **Delta_max**: Computes vertex scores and ordering
4. **MLVMST**: Builds minimum leaves spanning tree using parallel Borůvka
5. **CLNJ**: Reconstructs latent tree with hidden nodes
6. **Output**: Final phylogenetic tree structure

**Files Executed:**
- ✅ `main_clnj_boruvka.cpp` - Main orchestration
- ✅ `data_loader.cpp` - Input loading
- ✅ `algorithm2.cpp` - Algorithm 2 implementation
- ✅ `mlvmst.cpp` - MLVMST utilities
- ✅ `mlvmst_boruvka.cpp` - Parallel Borůvka MLVMST
- ✅ `clnj_fresh.cpp` - CLNJ implementation

**All code executes sequentially** through the pipeline, with parallelization only in the Borůvka phase (Step 3).

---

**Last Updated**: Based on code analysis of `main_clnj_boruvka.cpp` and related files
