# Topological Phylogeny - CLNJ Tree Reconstruction

A high-performance C++ implementation of phylogenetic tree reconstruction using **Algorithm 2**, **MLVMST (Borůvka)**, and **CLNJ (Complete Latent Neighbor Joining)**.

## Overview

This project implements a complete pipeline for reconstructing phylogenetic trees from distance matrices:

1. **Algorithm 2**: Constructs F_C (forest of cliques) and G_U (union graph) from input edges
2. **MLVMST**: Builds Minimum Leaves Vertex MST using parallel Borůvka algorithm
3. **CLNJ**: Complete Latent Neighbor Joining for tree reconstruction

## Project Structure

```
cpp_cursor_conversion/
├── src/                              # C++ source files
│   ├── algorithm2.cpp                # Algorithm 2 implementation
│   ├── clnj_fresh.cpp                # CLNJ implementation
│   ├── data_loader.cpp               # Input file parser
│   ├── mlvmst_boruvka.cpp            # Parallel Borůvka MLVMST
│   ├── mlvmst.cpp                    # Basic MLVMST (Kruskal)
│   └── main_clnj_boruvka.cpp         # Main executable
├── include/                          # Header files
│   ├── algorithm2.h
│   ├── clnj.h
│   ├── common_types.h
│   ├── data_loader.h
│   ├── mlvmst_boruvka.h
│   └── mlvmst.h
├── data/                             # Sample test datasets
│   ├── primates_genus_test.txt
│   ├── primates_species_test.txt
│   ├── carnivora_species_test.txt
│   └── small_test_*.txt
├── python_preprocessing/             # Python utilities
│   └── python_patristic_optimized.py
├── convert_nwk_to_test_data.py       # Convert .nwk to input format
├── CMakeLists.txt                    # Build configuration
└── README.md
```

## Dependencies

### Required
- **C++17 compiler** (g++ 7.0+ or clang++ 5.0+)
- **CMake** 3.10+
- **Eigen3** (matrix operations)
- **OpenMP** (parallel processing)

### Installation (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install build-essential cmake libeigen3-dev libomp-dev
```

### Python (Optional - for data conversion)
```bash
pip install ete3 numpy
```

## Build Instructions

```bash
# Clone the repository
git clone https://github.com/Jatin07gupta/Topological-phylogeny-.git
cd Topological-phylogeny-

# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make clnj_boruvka_main

# The executable will be in build/bin/
```

## Usage

### Basic Usage
```bash
./build/bin/clnj_boruvka_main <input_file> [ordering] [num_threads]
```

### Arguments
| Argument | Description | Options | Default |
|----------|-------------|---------|---------|
| `input_file` | Path to test data file | `.txt` file | Required |
| `ordering` | Vertex ordering method | `increasing`, `decreasing` | `increasing` |
| `num_threads` | Number of parallel threads | `1`, `2`, `4`, `8`, etc. | `4` |

### Examples
```bash
# Run on primates dataset with 4 threads
./build/bin/clnj_boruvka_main data/primates_genus_test.txt increasing 4

# Run on carnivora dataset with 8 threads
./build/bin/clnj_boruvka_main data/carnivora_species_test.txt increasing 8
```

## Input File Format

The input file is a plain text file with the following structure:

```
<number_of_vertices>
<number_of_edges>
<vertex_name_1>
<vertex_name_2>
...
<vertex_name_n>
<vertex_u> <vertex_v> <weight>
<vertex_u> <vertex_v> <weight>
...
```

### Example
```
5
10
Species_A
Species_B
Species_C
Species_D
Species_E
Species_A Species_B 0.5
Species_A Species_C 0.8
...
```

## Converting Newick Files

To convert a Newick tree file (`.nwk`) to the required input format:

```bash
# Activate Python environment (if needed)
source env/bin/activate

# Convert .nwk to .txt
python convert_nwk_to_test_data.py input_tree.nwk data/output_test.txt
```

## Performance

| Dataset | Species | Edges | Time (4 threads) |
|---------|---------|-------|------------------|
| Primates Genus | 70 | 2,415 | ~2s |
| Primates Species | 449 | 100,576 | ~30s |
| Carnivora | 277 | 38,226 | ~15s |

## Pipeline Flow

```
Input (.txt)
    │
    ▼
┌─────────────────┐
│  Data Loader    │  ← Reads vertices and edges
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Algorithm 2    │  ← Constructs F_C and G_U
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│    MLVMST       │  ← Builds MST using parallel Borůvka
│   (Borůvka)     │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│     CLNJ        │  ← Reconstructs phylogenetic tree
└────────┬────────┘
         │
         ▼
Output (Newick tree)
```

## Sample Datasets

| File | Species | Description |
|------|---------|-------------|
| `primates_genus_test.txt` | 70 | Primate genera |
| `primates_species_test.txt` | 449 | Primate species |
| `carnivora_species_test.txt` | 277 | Carnivore species |
| `small_test_*.txt` | 5-10 | Quick test cases |

## Documentation

- `MLVMST_INPUT_EXPLANATION.md` - Detailed input format explanation
- `EXECUTION_FLOW_ANALYSIS.md` - Complete pipeline flow analysis

## License

This project is for academic/research purposes.

## Author

Jatin Gupta
