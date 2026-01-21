#include "mlvmst_boruvka.h"
#include "algorithm2.h"
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <vector>
#include <mutex>

// OpenMP support (optional - will compile without if not available)
#ifdef _OPENMP
#include <omp.h>
#endif

namespace phylo {

//=============================================================================
// Helper Functions
//=============================================================================

/**
 * Canonicalize an undirected edge (min, max) for set operations.
 */
std::pair<Vertex, Vertex> canonicalize_edge(const Vertex& u, const Vertex& v) {
    return (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
}

/**
 * Compute VMST tie-breaking key for an edge.
 * 
 * Returns: (weight, min(rank[u], rank[v]), max(rank[u], rank[v]))
 */
struct EdgeSortKey {
    double weight;
    int min_rank;
    int max_rank;
    
    EdgeSortKey(double w, int min_r, int max_r) 
        : weight(w), min_rank(min_r), max_rank(max_r) {}
    
    bool operator<(const EdgeSortKey& other) const {
        if (weight != other.weight) {
            return weight < other.weight;
        }
        if (min_rank != other.min_rank) {
            return min_rank < other.min_rank;
        }
        return max_rank < other.max_rank;
    }
};

/**
 * Compute edge sort key for an edge given vertex ranks.
 */
EdgeSortKey compute_edge_key(
    const Edge& edge,
    const std::unordered_map<Vertex, int>& vertex_rank
) {
    Vertex u = std::get<0>(edge);
    Vertex v = std::get<1>(edge);
    double w = std::get<2>(edge);
    
    int rank_u = (vertex_rank.count(u) > 0) ? vertex_rank.at(u) : std::numeric_limits<int>::max();
    int rank_v = (vertex_rank.count(v) > 0) ? vertex_rank.at(v) : std::numeric_limits<int>::max();
    
    int min_rank = std::min(rank_u, rank_v);
    int max_rank = std::max(rank_u, rank_v);
    
    return EdgeSortKey(w, min_rank, max_rank);
}

//=============================================================================
// Hybrid Borůvka VMST Construction
//=============================================================================

std::set<std::pair<Vertex, Vertex>> build_vmst_from_order_hybrid(
    const VertexList& vertices,
    const EdgeList& edges,
    const VertexList& vertex_order
) {
    // Edge case: no vertices
    if (vertices.empty()) {
        return {};
    }
    
    // Create rank mapping from vertex order
    std::unordered_map<Vertex, int> vertex_rank;
    vertex_rank.reserve(vertex_order.size());
    for (size_t i = 0; i < vertex_order.size(); ++i) {
        vertex_rank[vertex_order[i]] = static_cast<int>(i);
    }
    
    // Preprocess edges: remove self-loops and keep minimum weight for duplicates
    struct PairHash {
        std::size_t operator()(const std::pair<Vertex, Vertex>& p) const {
            return std::hash<Vertex>{}(p.first) ^ (std::hash<Vertex>{}(p.second) << 1);
        }
    };
    
    std::unordered_map<std::pair<Vertex, Vertex>, double, PairHash> edge_weights;
    edge_weights.reserve(edges.size());
    
    for (const auto& edge : edges) {
        const Vertex& u = std::get<0>(edge);
        const Vertex& v = std::get<1>(edge);
        double w = std::get<2>(edge);
        
        if (u == v) {  // Skip self-loops
            continue;
        }
        
        std::pair<Vertex, Vertex> canonical = canonicalize_edge(u, v);
        
        auto it = edge_weights.find(canonical);
        if (it == edge_weights.end() || w < it->second) {
            edge_weights[canonical] = w;
        }
    }
    
    // Convert to clean edge list
    EdgeList E_clean;
    E_clean.reserve(edge_weights.size());
    for (const auto& pair : edge_weights) {
        E_clean.push_back(std::make_tuple(pair.first.first, pair.first.second, pair.second));
    }
    
    // Initialize DSU and MST edge set
    DisjointSetUnion dsu(vertices);
    std::set<std::pair<Vertex, Vertex>> vmst_edges;
    
    size_t num_components = vertices.size();
    size_t target_edges = vertices.size() - 1;
    
    // Main Borůvka-style phased loop
    while (num_components > 1 && vmst_edges.size() < target_edges) {
        // (a) Find best outgoing edge for each component (Borůvka-style)
        // PARALLEL-FRIENDLY: This loop could be parallelized per component
        // in a future optimized implementation, but currently sequential.
        std::unordered_map<Vertex, Edge> best_edge;
        
        for (const auto& edge : E_clean) {
            const Vertex& u = std::get<0>(edge);
            const Vertex& v = std::get<1>(edge);
            
            Vertex root_u = dsu.find(u);
            Vertex root_v = dsu.find(v);
            
            // Skip edges within the same component
            if (root_u == root_v) {
                continue;
            }
            
            EdgeSortKey key_e = compute_edge_key(edge, vertex_rank);
            
            // Update best edge for component root_u
            auto it_u = best_edge.find(root_u);
            if (it_u == best_edge.end()) {
                best_edge[root_u] = edge;
            } else {
                EdgeSortKey key_old = compute_edge_key(it_u->second, vertex_rank);
                if (key_e < key_old) {
                    best_edge[root_u] = edge;
                }
            }
            
            // Update best edge for component root_v
            auto it_v = best_edge.find(root_v);
            if (it_v == best_edge.end()) {
                best_edge[root_v] = edge;
            } else {
                EdgeSortKey key_old = compute_edge_key(it_v->second, vertex_rank);
                if (key_e < key_old) {
                    best_edge[root_v] = edge;
                }
            }
        }
        
        // If no best edges found, graph is disconnected
        if (best_edge.empty()) {
            break;
        }
        
        // (b) Collect candidate edges and deduplicate by canonical undirected pair
        std::unordered_map<std::pair<Vertex, Vertex>, Edge, PairHash> candidate_by_canon;
        candidate_by_canon.reserve(best_edge.size());
        
        for (const auto& pair : best_edge) {
            const Edge& e = pair.second;
            const Vertex& u = std::get<0>(e);
            const Vertex& v = std::get<1>(e);
            
            std::pair<Vertex, Vertex> canon = canonicalize_edge(u, v);
            
            auto it = candidate_by_canon.find(canon);
            if (it == candidate_by_canon.end()) {
                candidate_by_canon[canon] = e;
            } else {
                // Keep the edge with the smaller key (shouldn't happen, but safe)
                EdgeSortKey key_new = compute_edge_key(e, vertex_rank);
                EdgeSortKey key_old = compute_edge_key(it->second, vertex_rank);
                if (key_new < key_old) {
                    candidate_by_canon[canon] = e;
                }
            }
        }
        
        // (c) Sort candidate edges globally by VMST key (Kruskal-style)
        std::vector<std::pair<EdgeSortKey, Edge>> candidate_list;
        candidate_list.reserve(candidate_by_canon.size());
        
        for (const auto& pair : candidate_by_canon) {
            EdgeSortKey key = compute_edge_key(pair.second, vertex_rank);
            candidate_list.push_back({key, pair.second});
        }
        
        std::sort(candidate_list.begin(), candidate_list.end(),
                  [](const std::pair<EdgeSortKey, Edge>& a, 
                     const std::pair<EdgeSortKey, Edge>& b) {
                      return a.first < b.first;
                  });
        
        // (d) Process candidates Kruskal-style with DSU
        for (const auto& pair : candidate_list) {
            if (vmst_edges.size() >= target_edges) {
                break;
            }
            
            const Edge& e = pair.second;
            const Vertex& u = std::get<0>(e);
            const Vertex& v = std::get<1>(e);
            
            Vertex root_u = dsu.find(u);
            Vertex root_v = dsu.find(v);
            
            // Skip if already in same component (shouldn't happen after dedup, but safe)
            if (root_u == root_v) {
                continue;
            }
            
            // Add edge to MST and merge components
            std::pair<Vertex, Vertex> canon = canonicalize_edge(u, v);
            vmst_edges.insert(canon);
            dsu.union_sets(root_u, root_v);
            num_components--;  // Each union merges two components into one
        }
    }
    
    return vmst_edges;
}

//=============================================================================
// Parallel Hybrid Borůvka VMST Construction
//=============================================================================

std::set<std::pair<Vertex, Vertex>> build_vmst_from_order_hybrid_parallel(
    const VertexList& vertices,
    const EdgeList& edges,
    const VertexList& vertex_order,
    int num_threads
) {
    // Edge case: no vertices
    if (vertices.empty()) {
        return {};
    }
    
    // Set number of threads if specified
    #ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    #endif
    
    // Create rank mapping from vertex order
    std::unordered_map<Vertex, int> vertex_rank;
    vertex_rank.reserve(vertex_order.size());
    for (size_t i = 0; i < vertex_order.size(); ++i) {
        vertex_rank[vertex_order[i]] = static_cast<int>(i);
    }
    
    // Preprocess edges: remove self-loops and keep minimum weight for duplicates
    struct PairHash {
        std::size_t operator()(const std::pair<Vertex, Vertex>& p) const {
            return std::hash<Vertex>{}(p.first) ^ (std::hash<Vertex>{}(p.second) << 1);
        }
    };
    
    std::unordered_map<std::pair<Vertex, Vertex>, double, PairHash> edge_weights;
    edge_weights.reserve(edges.size());
    
    for (const auto& edge : edges) {
        const Vertex& u = std::get<0>(edge);
        const Vertex& v = std::get<1>(edge);
        double w = std::get<2>(edge);
        
        if (u == v) {  // Skip self-loops
            continue;
        }
        
        std::pair<Vertex, Vertex> canonical = canonicalize_edge(u, v);
        
        auto it = edge_weights.find(canonical);
        if (it == edge_weights.end() || w < it->second) {
            edge_weights[canonical] = w;
        }
    }
    
    // Convert to clean edge list
    EdgeList E_clean;
    E_clean.reserve(edge_weights.size());
    for (const auto& pair : edge_weights) {
        E_clean.push_back(std::make_tuple(pair.first.first, pair.first.second, pair.second));
    }
    
    // Initialize DSU and MST edge set
    DisjointSetUnion dsu(vertices);
    std::set<std::pair<Vertex, Vertex>> vmst_edges;
    
    size_t num_components = vertices.size();
    size_t target_edges = vertices.size() - 1;
    
    // Main Borůvka-style phased loop
    while (num_components > 1 && vmst_edges.size() < target_edges) {
        // (a) Find best outgoing edge for each component (Borůvka-style) - PARALLELIZED
        // STEP 1: Take snapshot of all current component roots (thread-safe snapshot)
        // This avoids needing locks during parallel processing
        std::unordered_map<Vertex, Vertex> vertex_to_root;
        vertex_to_root.reserve(vertices.size());
        
        for (const auto& v : vertices) {
            Vertex root = dsu.find(v);
            vertex_to_root[v] = root;
        }
        
        // If only one component, we're done
        std::set<Vertex> component_roots;
        for (const auto& pair : vertex_to_root) {
            component_roots.insert(pair.second);
        }
        
        if (component_roots.size() <= 1) {
            break;
        }
        
        // STEP 2: Parallel edge scanning with thread-local best edges
        // Each thread maintains its own local_best map, then we merge at the end
        std::vector<std::unordered_map<Vertex, Edge>> thread_best_edges;
        
        #ifdef _OPENMP
        // Determine number of threads
        int actual_threads = (num_threads > 0) ? num_threads : omp_get_max_threads();
        thread_best_edges.resize(actual_threads);
        
        // Parallel edge scanning loop
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            auto& local_best = thread_best_edges[thread_id];
            
            // Each thread processes a portion of edges
            #pragma omp for nowait
            for (size_t i = 0; i < E_clean.size(); ++i) {
                const Edge& edge = E_clean[i];
                const Vertex& u = std::get<0>(edge);
                const Vertex& v = std::get<1>(edge);
                
                // Use snapshot mapping (read-only, no locks needed)
                auto it_u = vertex_to_root.find(u);
                auto it_v = vertex_to_root.find(v);
                
                if (it_u == vertex_to_root.end() || it_v == vertex_to_root.end()) {
                    continue;
                }
                
                Vertex root_u = it_u->second;
                Vertex root_v = it_v->second;
                
                // Skip edges within the same component
                if (root_u == root_v) {
                    continue;
                }
                
                EdgeSortKey key_e = compute_edge_key(edge, vertex_rank);
                
                // Update best edge for component root_u
                auto it_best_u = local_best.find(root_u);
                if (it_best_u == local_best.end()) {
                    local_best[root_u] = edge;
                } else {
                    EdgeSortKey key_old = compute_edge_key(it_best_u->second, vertex_rank);
                    if (key_e < key_old) {
                        local_best[root_u] = edge;
                    }
                }
                
                // Update best edge for component root_v
                auto it_best_v = local_best.find(root_v);
                if (it_best_v == local_best.end()) {
                    local_best[root_v] = edge;
                } else {
                    EdgeSortKey key_old = compute_edge_key(it_best_v->second, vertex_rank);
                    if (key_e < key_old) {
                        local_best[root_v] = edge;
                    }
                }
            }
        }
        
        // STEP 3: Merge thread-local results into global best_edge (sequential)
        std::unordered_map<Vertex, Edge> best_edge;
        for (const auto& local_best : thread_best_edges) {
            for (const auto& pair : local_best) {
                const Vertex& comp_root = pair.first;
                const Edge& e = pair.second;
                
                auto it = best_edge.find(comp_root);
                if (it == best_edge.end()) {
                    best_edge[comp_root] = e;
                } else {
                    EdgeSortKey key_new = compute_edge_key(e, vertex_rank);
                    EdgeSortKey key_old = compute_edge_key(it->second, vertex_rank);
                    if (key_new < key_old) {
                        best_edge[comp_root] = e;
                    }
                }
            }
        }
        #else
        // Fallback to sequential if OpenMP not available
        std::unordered_map<Vertex, Edge> best_edge;
        
        for (const auto& edge : E_clean) {
            const Vertex& u = std::get<0>(edge);
            const Vertex& v = std::get<1>(edge);
            
            auto it_u = vertex_to_root.find(u);
            auto it_v = vertex_to_root.find(v);
            
            if (it_u == vertex_to_root.end() || it_v == vertex_to_root.end()) {
                continue;
            }
            
            Vertex root_u = it_u->second;
            Vertex root_v = it_v->second;
            
            if (root_u == root_v) {
                continue;
            }
            
            EdgeSortKey key_e = compute_edge_key(edge, vertex_rank);
            
            auto it_best_u = best_edge.find(root_u);
            if (it_best_u == best_edge.end()) {
                best_edge[root_u] = edge;
            } else {
                EdgeSortKey key_old = compute_edge_key(it_best_u->second, vertex_rank);
                if (key_e < key_old) {
                    best_edge[root_u] = edge;
                }
            }
            
            auto it_best_v = best_edge.find(root_v);
            if (it_best_v == best_edge.end()) {
                best_edge[root_v] = edge;
            } else {
                EdgeSortKey key_old = compute_edge_key(it_best_v->second, vertex_rank);
                if (key_e < key_old) {
                    best_edge[root_v] = edge;
                }
            }
        }
        #endif
        
        // If no best edges found, graph is disconnected
        if (best_edge.empty()) {
            break;
        }
        
        // (b) Collect candidate edges and deduplicate by canonical undirected pair (SEQUENTIAL)
        std::unordered_map<std::pair<Vertex, Vertex>, Edge, PairHash> candidate_by_canon;
        candidate_by_canon.reserve(best_edge.size());
        
        for (const auto& pair : best_edge) {
            const Edge& e = pair.second;
            const Vertex& u = std::get<0>(e);
            const Vertex& v = std::get<1>(e);
            
            std::pair<Vertex, Vertex> canon = canonicalize_edge(u, v);
            
            auto it = candidate_by_canon.find(canon);
            if (it == candidate_by_canon.end()) {
                candidate_by_canon[canon] = e;
            } else {
                // Keep the edge with the smaller key (shouldn't happen, but safe)
                EdgeSortKey key_new = compute_edge_key(e, vertex_rank);
                EdgeSortKey key_old = compute_edge_key(it->second, vertex_rank);
                if (key_new < key_old) {
                    candidate_by_canon[canon] = e;
                }
            }
        }
        
        // (c) Sort candidate edges globally by VMST key (Kruskal-style) (SEQUENTIAL)
        std::vector<std::pair<EdgeSortKey, Edge>> candidate_list;
        candidate_list.reserve(candidate_by_canon.size());
        
        for (const auto& pair : candidate_by_canon) {
            EdgeSortKey key = compute_edge_key(pair.second, vertex_rank);
            candidate_list.push_back({key, pair.second});
        }
        
        std::sort(candidate_list.begin(), candidate_list.end(),
                  [](const std::pair<EdgeSortKey, Edge>& a, 
                     const std::pair<EdgeSortKey, Edge>& b) {
                      return a.first < b.first;
                  });
        
        // (d) Process candidates Kruskal-style with DSU (SEQUENTIAL)
        for (const auto& pair : candidate_list) {
            if (vmst_edges.size() >= target_edges) {
                break;
            }
            
            const Edge& e = pair.second;
            const Vertex& u = std::get<0>(e);
            const Vertex& v = std::get<1>(e);
            
            Vertex root_u = dsu.find(u);
            Vertex root_v = dsu.find(v);
            
            // Skip if already in same component (shouldn't happen after dedup, but safe)
            if (root_u == root_v) {
                continue;
            }
            
            // Add edge to MST and merge components
            std::pair<Vertex, Vertex> canon = canonicalize_edge(u, v);
            vmst_edges.insert(canon);
            dsu.union_sets(root_u, root_v);
            num_components--;  // Each union merges two components into one
        }
    }
    
    return vmst_edges;
}

} // namespace phylo

