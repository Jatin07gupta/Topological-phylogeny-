#include "mlvmst.h"
#include "algorithm2.h"
#include <algorithm>
#include <unordered_map>

namespace phylo {

//=============================================================================
// Helper Functions
//=============================================================================

std::unordered_map<Vertex, VertexSet> build_adjacency(
    const std::set<std::pair<Vertex, Vertex>>& edges
) {
    std::unordered_map<Vertex, VertexSet> adjacency;
    // Pre-allocate: estimate 2*edges.size() total insertions
    adjacency.reserve(edges.size() * 2);
    
    for (const auto& edge : edges) {
        const Vertex& u = edge.first;
        const Vertex& v = edge.second;
        adjacency[u].insert(v);
        adjacency[v].insert(u);
    }
    
    return adjacency;
}

//=============================================================================
// Delta_max Computation
//=============================================================================

std::map<Vertex, int> compute_delta_max(
    const VertexList& vertices,
    const std::vector<VertexSet>& F_C,
    const std::unordered_map<Vertex, VertexSet>& G_U_adj
) {
    // Step 1: Order F_C by decreasing set size (F_C^≥)
    // OPTIMIZED: Use vector of indices/pointers to avoid copying large sets
    std::vector<size_t> F_C_indices;
    F_C_indices.reserve(F_C.size());
    for (size_t idx = 0; idx < F_C.size(); ++idx) {
        if (!F_C[idx].empty()) {
            F_C_indices.push_back(idx);
        }
    }
    
    std::sort(F_C_indices.begin(), F_C_indices.end(),
              [&F_C](size_t a, size_t b) {
                  return F_C[a].size() > F_C[b].size();  // Decreasing order
              });
    
    std::map<Vertex, int> delta_max;
    
    // Step 2: Compute delta_max for each vertex
    for (const auto& i : vertices) {
        // Get neighbors of i in G_U
        VertexSet N_i;
        auto it = G_U_adj.find(i);
        if (it != G_U_adj.end()) {
            N_i.reserve(it->second.size());  // Pre-allocate
            N_i = it->second;  // Copy neighbors
        }
        
        // Initialize delta_max(i) = 0
        delta_max[i] = 0;
        
        // Greedy covering: process sets in decreasing size order
        for (size_t idx : F_C_indices) {
            const VertexSet& C = F_C[idx];
            
            // Check if i is in C (quick check first)
            if (C.find(i) != C.end()) {
                continue;  // Skip if i is in C
            }
            
            // OPTIMIZED: Check intersection - iterate smaller set for efficiency
            bool has_intersection = false;
            if (C.size() < N_i.size()) {
                // C is smaller, iterate C
                for (const auto& v : C) {
                    if (N_i.find(v) != N_i.end()) {
                        has_intersection = true;
                        break;
                    }
                }
            } else {
                // N_i is smaller or equal, iterate N_i
                for (const auto& v : N_i) {
                    if (C.find(v) != C.end()) {
                        has_intersection = true;
                        break;
                    }
                }
            }
            
            if (has_intersection) {
                // This set covers some neighbors and doesn't contain i
                delta_max[i]++;
                
                // OPTIMIZED: Set difference - iterate smaller set
                if (C.size() < N_i.size()) {
                    // C is smaller, iterate C and remove from N_i
                    for (const auto& v : C) {
                        N_i.erase(v);
                    }
                } else {
                    // N_i is smaller, iterate N_i and check C
                    for (auto nit = N_i.begin(); nit != N_i.end(); ) {
                        if (C.find(*nit) != C.end()) {
                            nit = N_i.erase(nit);
                        } else {
                            ++nit;
                        }
                    }
                }
            }
            
            // Early termination: if all neighbors covered, stop
            if (N_i.empty()) {
                break;
            }
        }
    }
    
    return delta_max;
}

//=============================================================================
// Total Order Construction
//=============================================================================

VertexList build_total_order(
    const VertexList& vertices,
    const std::map<Vertex, int>& delta_max,
    const std::string& ordering
) {
    // OPTIMIZED: Pre-compute delta values to avoid repeated map lookups
    std::vector<std::pair<Vertex, int>> vertex_delta;
    vertex_delta.reserve(vertices.size());
    
    for (const auto& v : vertices) {
        auto it = delta_max.find(v);
        int delta = (it != delta_max.end()) ? it->second : 0;
        vertex_delta.push_back({v, delta});
    }
    
    if (ordering == "increasing") {
        // Paper's stated ordering: nondecreasing (smaller delta_max first)
        std::sort(vertex_delta.begin(), vertex_delta.end(),
                  [](const std::pair<Vertex, int>& a, const std::pair<Vertex, int>& b) {
                      if (a.second != b.second) {
                          return a.second < b.second;
                      }
                      return a.first < b.first;  // Tie-break: lexicographic
                  });
    } else {
        // Empirically better: nonincreasing (larger delta_max first) - REVERSED
        std::sort(vertex_delta.begin(), vertex_delta.end(),
                  [](const std::pair<Vertex, int>& a, const std::pair<Vertex, int>& b) {
                      if (a.second != b.second) {
                          return a.second > b.second;  // Larger first
                      }
                      return a.first < b.first;  // Tie-break: lexicographic
                  });
    }
    
    VertexList result;
    result.reserve(vertex_delta.size());
    for (const auto& p : vertex_delta) {
        result.push_back(p.first);
    }
    
    return result;
}

//=============================================================================
// VMST Construction (Algorithm 1)
//=============================================================================

std::set<std::pair<Vertex, Vertex>> build_vmst_from_order(
    const VertexList& vertices,
    const EdgeList& edges,
    const VertexList& vertex_order
) {
    if (vertices.empty()) {
        return {};
    }
    
    // OPTIMIZED: Create rank mapping with pre-allocation
    std::unordered_map<Vertex, int> vertex_rank;
    vertex_rank.reserve(vertex_order.size());
    for (size_t i = 0; i < vertex_order.size(); ++i) {
        vertex_rank[vertex_order[i]] = static_cast<int>(i);
    }
    
    // OPTIMIZED: Use unordered_map instead of map for O(1) average lookup
    struct PairHash {
        std::size_t operator()(const std::pair<Vertex, Vertex>& p) const {
            return std::hash<Vertex>{}(p.first) ^ (std::hash<Vertex>{}(p.second) << 1);
        }
    };
    
    std::unordered_map<std::pair<Vertex, Vertex>, double, PairHash> edge_weights;
    edge_weights.reserve(edges.size());  // Pre-allocate
    for (const auto& edge : edges) {
        const Vertex& u = std::get<0>(edge);
        const Vertex& v = std::get<1>(edge);
        double w = std::get<2>(edge);
        
        if (u == v) {  // Skip self-loops
            continue;
        }
        
        std::pair<Vertex, Vertex> canonical = (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
        
        auto it = edge_weights.find(canonical);
        if (it == edge_weights.end() || w < it->second) {
            edge_weights[canonical] = w;
        }
    }
    
    // Convert to edge list with sorting
    EdgeList E_sorted;
    E_sorted.reserve(edge_weights.size());
    
    for (const auto& pair : edge_weights) {
        E_sorted.push_back(std::make_tuple(pair.first.first, pair.first.second, pair.second));
    }
    
    // OPTIMIZED: Pre-compute ranks for all edges to avoid repeated lookups during sort
    struct EdgeWithRank {
        Edge edge;
        int min_rank;
        int max_rank;
        
        EdgeWithRank(const Edge& e, const std::unordered_map<Vertex, int>& rank_map) 
            : edge(e) {
            Vertex u = std::get<0>(e);
            Vertex v = std::get<1>(e);
            int rank_u = rank_map.count(u) ? rank_map.at(u) : 1000000;
            int rank_v = rank_map.count(v) ? rank_map.at(v) : 1000000;
            min_rank = std::min(rank_u, rank_v);
            max_rank = std::max(rank_u, rank_v);
        }
    };
    
    std::vector<EdgeWithRank> edges_with_rank;
    edges_with_rank.reserve(E_sorted.size());
    for (const auto& e : E_sorted) {
        edges_with_rank.emplace_back(e, vertex_rank);
    }
    
    // Sort once with pre-computed ranks
    std::sort(edges_with_rank.begin(), edges_with_rank.end(),
              [](const EdgeWithRank& a, const EdgeWithRank& b) {
                  double w_a = std::get<2>(a.edge);
                  double w_b = std::get<2>(b.edge);
                  
                  if (w_a != w_b) {
                      return w_a < w_b;
                  }
                  
                  if (a.min_rank != b.min_rank) {
                      return a.min_rank < b.min_rank;
                  }
                  return a.max_rank < b.max_rank;
              });
    
    // Update E_sorted with sorted edges
    E_sorted.clear();
    E_sorted.reserve(edges_with_rank.size());
    for (const auto& ewr : edges_with_rank) {
        E_sorted.push_back(ewr.edge);
    }
    
    // Kruskal's algorithm
    DisjointSetUnion dsu(vertices);
    std::set<std::pair<Vertex, Vertex>> vmst_edges;
    
    for (const auto& edge : E_sorted) {
        const Vertex& u = std::get<0>(edge);
        const Vertex& v = std::get<1>(edge);
        
        Vertex root_u = dsu.find(u);
        Vertex root_v = dsu.find(v);
        
        if (root_u != root_v) {
            // Add edge to VMST
            std::pair<Vertex, Vertex> canonical = (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
            vmst_edges.insert(canonical);
            dsu.union_sets(u, v);
            
            // Early termination: MST has exactly n-1 edges
            if (vmst_edges.size() == vertices.size() - 1) {
                break;
            }
        }
    }
    
    return vmst_edges;
}

//=============================================================================
// Leaf Counting
//=============================================================================

int count_leaves(
    const std::set<std::pair<Vertex, Vertex>>& edges,
    const VertexList& vertices
) {
    if (edges.empty()) {
        return static_cast<int>(vertices.size());  // All vertices are leaves (disconnected)
    }
    
    // Count degree of each vertex
    std::unordered_map<Vertex, int> degree;
    degree.reserve(vertices.size());  // Pre-allocate
    
    for (const auto& edge : edges) {
        degree[edge.first]++;
        degree[edge.second]++;
    }
    
    // A leaf has degree 1
    int leaf_count = 0;
    for (const auto& v : vertices) {
        if (degree[v] == 1) {
            leaf_count++;
        }
    }
    
    return leaf_count;
}

//=============================================================================
// Main Algorithm 3: MLVMST
//=============================================================================

MLVMSTResult build_mlvmst(
    const VertexList& vertices,
    const EdgeList& edges,
    const std::vector<VertexSet>& F_C,
    const std::set<std::pair<Vertex, Vertex>>& G_U_edges,
    const std::string& ordering
) {
    // Edge case: empty graph
    if (vertices.empty()) {
        MLVMSTResult result;
        result.order = {};
        result.delta_max = {};
        result.M_star_edges = {};
        result.leaf_count = 0;
        return result;
    }
    
    // Edge case: single vertex
    if (vertices.size() == 1) {
        std::map<Vertex, int> delta_max;
        delta_max[vertices[0]] = 0;
        MLVMSTResult result;
        result.order = vertices;
        result.delta_max = delta_max;
        result.M_star_edges = {};
        result.leaf_count = 1;
        return result;
    }
    
    // Step 1: Build adjacency list for G_U
    std::unordered_map<Vertex, VertexSet> G_U_adj = build_adjacency(G_U_edges);
    
    // Step 2: Compute delta_max for each vertex
    std::map<Vertex, int> delta_max = compute_delta_max(vertices, F_C, G_U_adj);
    
    // Step 3: Build total order <*
    VertexList order = build_total_order(vertices, delta_max, ordering);
    
    // Step 4: Run VMST with order <* to obtain M*
    std::set<std::pair<Vertex, Vertex>> M_star_edges = build_vmst_from_order(vertices, edges, order);
    
    // Step 5: Count leaves in M*
    int leaf_count = count_leaves(M_star_edges, vertices);
    
    MLVMSTResult result;
    result.order = order;
    result.delta_max = delta_max;
    result.M_star_edges = M_star_edges;
    result.leaf_count = leaf_count;
    return result;
}

//=============================================================================
// Comparison Operator
//=============================================================================

bool MLVMSTResult::operator==(const MLVMSTResult& other) const {
    if (order.size() != other.order.size() ||
        delta_max.size() != other.delta_max.size() ||
        M_star_edges.size() != other.M_star_edges.size() ||
        leaf_count != other.leaf_count) {
        return false;
    }
    
    // Compare order (should be identical)
    for (size_t i = 0; i < order.size(); ++i) {
        if (order[i] != other.order[i]) {
            return false;
        }
    }
    
    // Compare delta_max
    for (const auto& pair : delta_max) {
        auto it = other.delta_max.find(pair.first);
        if (it == other.delta_max.end() || it->second != pair.second) {
            return false;
        }
    }
    
    // Compare M_star_edges
    if (M_star_edges != other.M_star_edges) {
        return false;
    }
    
    return true;
}

} // namespace phylo

