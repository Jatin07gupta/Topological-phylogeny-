/**
 * Main CLNJ executable using Borůvka Parallel for MLVMST construction
 * 
 * This version uses the faster Borůvka parallel implementation instead of
 * the Kruskal-based build_mlvmst() function.
 * 
 * Pipeline:
 * 1. Algorithm 2: Construct F_C and G_U
 * 2. Compute delta_max and vertex order
 * 3. Build MLVMST using Borůvka Parallel
 * 4. Run CLNJ on the MLVMST
 */

#include "clnj.h"
#include "mlvmst.h"
#include "mlvmst_boruvka.h"
#include "algorithm2.h"
#include "data_loader.h"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include <cmath>

// OpenMP support (optional - will compile without if not available)
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace phylo;

/**
 * Build distance matrix from edge list.
 * Assumes the edge list contains all pairwise distances.
 */
std::vector<std::vector<double>> build_distance_matrix(
    const VertexList& vertices,
    const EdgeList& edges
) {
    int m = static_cast<int>(vertices.size());
    std::vector<std::vector<double>> distance_matrix(m, std::vector<double>(m, 0.0));
    
    // Create vertex to index mapping
    std::unordered_map<Vertex, int> vertex_to_index;
    for (int i = 0; i < m; ++i) {
        vertex_to_index[vertices[i]] = i;
    }
    
    // Fill distance matrix from edges
    for (const auto& edge : edges) {
        const Vertex& u = std::get<0>(edge);
        const Vertex& v = std::get<1>(edge);
        double weight = std::get<2>(edge);
        
        auto it_u = vertex_to_index.find(u);
        auto it_v = vertex_to_index.find(v);
        
        if (it_u != vertex_to_index.end() && it_v != vertex_to_index.end()) {
            int i = it_u->second;
            int j = it_v->second;
            distance_matrix[i][j] = weight;
            distance_matrix[j][i] = weight;  // Symmetric
        }
    }
    
    return distance_matrix;
}

/**
 * Build MLVMST adjacency matrix from MLVMST edges.
 */
std::vector<std::vector<int>> build_mlvmst_adjacency(
    const VertexList& vertices,
    const std::set<std::pair<Vertex, Vertex>>& mlvmst_edges
) {
    int m = static_cast<int>(vertices.size());
    std::vector<std::vector<int>> adjacency(m, std::vector<int>(m, 0));
    
    // Create vertex to index mapping
    std::unordered_map<Vertex, int> vertex_to_index;
    for (int i = 0; i < m; ++i) {
        vertex_to_index[vertices[i]] = i;
    }
    
    // Fill adjacency matrix from MLVMST edges
    for (const auto& edge : mlvmst_edges) {
        const Vertex& u = edge.first;
        const Vertex& v = edge.second;
        
        auto it_u = vertex_to_index.find(u);
        auto it_v = vertex_to_index.find(v);
        
        if (it_u != vertex_to_index.end() && it_v != vertex_to_index.end()) {
            int i = it_u->second;
            int j = it_v->second;
            adjacency[i][j] = 1;
            adjacency[j][i] = 1;  // Symmetric
        }
    }
    
    return adjacency;
}

void print_clnj_result(const CLNJResult& result, int m) {
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "CLNJ RESULTS\n";
    std::cout << std::string(80, '=') << "\n";
    
    int total_nodes = static_cast<int>(result.adjacency.size());
    int num_hidden = static_cast<int>(result.hidden_info.size());
    
    // Count total edges
    int total_edges = 0;
    for (const auto& pair : result.adjacency) {
        total_edges += static_cast<int>(pair.second.size());
    }
    total_edges /= 2;  // Each edge counted twice
    
    std::cout << "\nSummary:\n";
    std::cout << "  Observed nodes: " << m << "\n";
    std::cout << "  Hidden nodes: " << num_hidden << "\n";
    std::cout << "  Total nodes: " << total_nodes << "\n";
    std::cout << "  Total edges: " << total_edges << "\n";
    std::cout << "  Expected edges (for tree): " << (total_nodes - 1) << "\n";
    std::cout << "  Valid tree: " << (total_edges == total_nodes - 1 ? "Yes" : "No") << "\n";
    
    // Output adjacency for verification (parseable format)
    std::cout << "\nAdjacency (for verification):\n";
    
    // Sort nodes for consistent output
    std::vector<int> sorted_nodes;
    for (const auto& pair : result.adjacency) {
        sorted_nodes.push_back(pair.first);
    }
    std::sort(sorted_nodes.begin(), sorted_nodes.end());
    
    for (int node : sorted_nodes) {
        auto adj_it = result.adjacency.find(node);
        if (adj_it == result.adjacency.end()) continue;
        
        std::vector<int> neighbors(adj_it->second.begin(), adj_it->second.end());
        std::sort(neighbors.begin(), neighbors.end());
        
        std::cout << "NODE:" << node << ":";
        for (size_t i = 0; i < neighbors.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << neighbors[i];
        }
        std::cout << "\n";
    }
    
    // Output hidden nodes info
    if (!result.hidden_info.empty()) {
        std::cout << "\nHidden nodes (for verification):\n";
        std::vector<int> sorted_hidden;
        for (const auto& pair : result.hidden_info) {
            sorted_hidden.push_back(pair.first);
        }
        std::sort(sorted_hidden.begin(), sorted_hidden.end());
        
        for (int hidden_id : sorted_hidden) {
            const HiddenNodeInfo& info = result.hidden_info.at(hidden_id);
            std::cout << "HIDDEN:" << hidden_id << ":anchor=" << info.anchor
                      << ":dist=" << std::fixed << std::setprecision(10) << info.dist_to_anchor << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <test_data_file> [ordering] [num_threads]\n";
        std::cerr << "\nExample:\n";
        std::cerr << "  " << argv[0] << " ../data/carnivora_test.txt increasing 4\n";
        std::cerr << "  " << argv[0] << " ../data/primates_test.txt decreasing 0\n";
        std::cerr << "\nArguments:\n";
        std::cerr << "  ordering: 'increasing' or 'decreasing' (default: 'increasing')\n";
        std::cerr << "  num_threads: Number of threads for Borůvka parallel (0 = use OpenMP default, default: 4)\n";
        return 1;
    }
    
    std::string filename = argv[1];
    
    // Get ordering from command line argument (default: "increasing")
    std::string ordering = "increasing";
    if (argc >= 3) {
        ordering = argv[2];
        if (ordering != "increasing" && ordering != "decreasing") {
            std::cerr << "Warning: ordering must be 'increasing' or 'decreasing', using 'increasing'\n";
            ordering = "increasing";
        }
    }
    
    // Get number of threads from command line (default: 4)
    int num_threads = 4;
    if (argc >= 4) {
        num_threads = std::atoi(argv[3]);
        if (num_threads < 0) {
            num_threads = 0;
        }
    }
    
    #ifdef _OPENMP
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    std::cout << "Using " << num_threads << " threads for Borůvka parallel\n";
    #else
    std::cout << "OpenMP not available - using sequential Borůvka\n";
    num_threads = 0;
    #endif
    
    std::cout << "Loading test data from: " << filename << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    VertexList vertices;
    EdgeList edges;
    
    if (!load_test_data(filename, vertices, edges)) {
        std::cerr << "Failed to load test data\n";
        return 1;
    }
    
    int m = static_cast<int>(vertices.size());
    std::cout << "\nLoaded:\n";
    std::cout << "  Vertices: " << m << "\n";
    std::cout << "  Edges: " << edges.size() << "\n";
    std::cout << "  Ordering: " << ordering << "\n";
    
    // Step 1: Run Algorithm 2 to get F_C and G_U
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "STEP 1: Algorithm 2 - Construct F_C and G_U\n";
    std::cout << std::string(80, '=') << "\n";
    
    auto start_alg2 = std::chrono::high_resolution_clock::now();
    FC_GU_Result alg2_result = construct_FC_and_GU(vertices, edges);
    auto end_alg2 = std::chrono::high_resolution_clock::now();
    auto alg2_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_alg2 - start_alg2);
    
    std::cout << "  F_C components: " << alg2_result.F_C.size() << "\n";
    std::cout << "  G_U edges: " << alg2_result.G_U_edges.size() << "\n";
    std::cout << "  Time: " << alg2_time.count() << " ms\n";
    
    // Step 2: Compute delta_max and vertex order
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "STEP 2: Compute Delta_max and Vertex Order\n";
    std::cout << std::string(80, '=') << "\n";
    
    auto start_delta = std::chrono::high_resolution_clock::now();
    
    // Build G_U adjacency list
    std::unordered_map<Vertex, VertexSet> G_U_adj = build_adjacency(alg2_result.G_U_edges);
    
    // Compute delta_max
    std::map<Vertex, int> delta_max = compute_delta_max(vertices, alg2_result.F_C, G_U_adj);
    
    // Build vertex order from delta_max
    VertexList vertex_order = build_total_order(vertices, delta_max, ordering);
    
    auto end_delta = std::chrono::high_resolution_clock::now();
    auto delta_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_delta - start_delta);
    
    std::cout << "  Vertex order size: " << vertex_order.size() << "\n";
    std::cout << "  Time: " << delta_time.count() << " ms\n";
    
    // Step 3: Build MLVMST using Borůvka Parallel
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "STEP 3: Build MLVMST using Borůvka Parallel\n";
    std::cout << std::string(80, '=') << "\n";
    
    auto start_mlvmst = std::chrono::high_resolution_clock::now();
    std::set<std::pair<Vertex, Vertex>> mlvmst_edges = build_vmst_from_order_hybrid_parallel(
        vertices, edges, vertex_order, num_threads
    );
    auto end_mlvmst = std::chrono::high_resolution_clock::now();
    auto mlvmst_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_mlvmst - start_mlvmst);
    
    // Count leaves
    int leaf_count = count_leaves(mlvmst_edges, vertices);
    
    std::cout << "  MLVMST edges: " << mlvmst_edges.size() << "\n";
    std::cout << "  Leaf count: " << leaf_count << "\n";
    std::cout << "  Time: " << mlvmst_time.count() << " ms\n";
    
    // Step 4: Build distance matrix and MLVMST adjacency
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "STEP 4: Prepare CLNJ Input\n";
    std::cout << std::string(80, '=') << "\n";
    
    auto start_prep = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> distance_matrix_vec = build_distance_matrix(vertices, edges);
    std::vector<std::vector<int>> mst_adjacency_vec = build_mlvmst_adjacency(vertices, mlvmst_edges);
    
    // Convert to Eigen matrices
    Eigen::MatrixXd distance_matrix(m, m);
    Eigen::MatrixXi mst_adjacency(m, m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            distance_matrix(i, j) = distance_matrix_vec[i][j];
            mst_adjacency(i, j) = mst_adjacency_vec[i][j];
        }
    }
    
    auto end_prep = std::chrono::high_resolution_clock::now();
    auto prep_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_prep - start_prep);
    
    std::cout << "  Distance matrix: " << m << "x" << m << "\n";
    std::cout << "  MLVMST adjacency: " << m << "x" << m << "\n";
    std::cout << "  Time: " << prep_time.count() << " ms\n";
    
    // Step 5: Run CLNJ
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "STEP 5: CLNJ - Reconstruct Latent Tree\n";
    std::cout << std::string(80, '=') << "\n";
    
    double threshold = -std::log(0.9);
    std::cout << "  Threshold: " << threshold << " (=-log(0.9))" << std::endl;
    
    auto start_clnj = std::chrono::high_resolution_clock::now();
    CLNJResult clnj_result = CLNJ_clean(distance_matrix, mst_adjacency, threshold, true);
    auto end_clnj = std::chrono::high_resolution_clock::now();
    auto clnj_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_clnj - start_clnj);
    
    std::cout << "  CLNJ time: " << clnj_time.count() << " ms\n";
    
    // Print results
    print_clnj_result(clnj_result, m);
    
    // Summary
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_clnj - start_alg2);
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "TIMING SUMMARY\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << "  Algorithm 2: " << alg2_time.count() << " ms\n";
    std::cout << "  Delta_max & Order: " << delta_time.count() << " ms\n";
    std::cout << "  MLVMST (Borůvka Parallel): " << mlvmst_time.count() << " ms\n";
    std::cout << "  Preparation: " << prep_time.count() << " ms\n";
    std::cout << "  CLNJ: " << clnj_time.count() << " ms\n";
    std::cout << "  Total: " << total_time.count() << " ms\n";
    
    return 0;
}

