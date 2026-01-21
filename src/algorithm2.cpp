#include "algorithm2.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <functional>

namespace phylo {

// Hash function for VertexSet (for fast deduplication)
struct VertexSetHash {
    std::size_t operator()(const VertexSet& s) const {
        std::size_t hash = 0;
        for (const auto& v : s) {
            hash ^= std::hash<Vertex>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

// Hash function for edge pairs
struct PairHash {
    std::size_t operator()(const std::pair<Vertex, Vertex>& p) const {
        return std::hash<Vertex>{}(p.first) ^ (std::hash<Vertex>{}(p.second) << 1);
    }
};

// Static member definition
const VertexSet ComponentTracker::empty_set_;

//=============================================================================
// DisjointSetUnion Implementation
//=============================================================================

DisjointSetUnion::DisjointSetUnion(const VertexList& vertices) {
    parent_.reserve(vertices.size());
    rank_.reserve(vertices.size());
    for (const auto& v : vertices) {
        parent_[v] = v;
        rank_[v] = 0;
    }
}

Vertex DisjointSetUnion::find(const Vertex& v) {
    if (parent_[v] != v) {
        parent_[v] = find(parent_[v]);  // Path compression
    }
    return parent_[v];
}

Vertex DisjointSetUnion::union_sets(const Vertex& u, const Vertex& v) {
    Vertex root_u = find(u);
    Vertex root_v = find(v);
    
    if (root_u == root_v) {
        return root_u;  // Already in same component
    }
    
    // Union by rank: attach smaller tree under larger tree
    if (rank_[root_u] < rank_[root_v]) {
        parent_[root_u] = root_v;
        return root_v;
    } else if (rank_[root_u] > rank_[root_v]) {
        parent_[root_v] = root_u;
        return root_u;
    } else {
        parent_[root_v] = root_u;
        rank_[root_u] += 1;
        return root_u;
    }
}

//=============================================================================
// ComponentTracker Implementation
//=============================================================================

ComponentTracker::ComponentTracker(const VertexList& vertices) {
    component_sets_.reserve(vertices.size());
    for (const auto& v : vertices) {
        component_sets_[v] = {v};
    }
}

void ComponentTracker::union_sets(const Vertex& old_root_u, const Vertex& old_root_v, const Vertex& new_root) {
    // Use move semantics to avoid copying
    VertexSet vertices_u = std::move(component_sets_[old_root_u]);
    VertexSet vertices_v = std::move(component_sets_[old_root_v]);
    
    component_sets_.erase(old_root_u);
    component_sets_.erase(old_root_v);
    
    // Merge sets efficiently - insert smaller into larger
    if (vertices_u.size() < vertices_v.size()) {
        vertices_v.insert(vertices_u.begin(), vertices_u.end());
        component_sets_[new_root] = std::move(vertices_v);
    } else {
        vertices_u.insert(vertices_v.begin(), vertices_v.end());
        component_sets_[new_root] = std::move(vertices_u);
    }
}

const VertexSet& ComponentTracker::component_vertices(const Vertex& root) const {
    auto it = component_sets_.find(root);
    if (it != component_sets_.end()) {
        return it->second;
    }
    return empty_set_;
}

//=============================================================================
// Helper Functions
//=============================================================================

namespace {
    /**
     * Add all singleton components {v} to F_C.
     * 
     * According to the paper, F_C should contain all components that appear
     * in every MST. Before processing any edges (w < min weight), every vertex
     * is its own singleton component.
     */
    void add_all_singletons(const VertexList& vertex_list, std::vector<VertexSet>& F_C) {
        F_C.reserve(F_C.size() + vertex_list.size());
        for (const auto& v : vertex_list) {
            F_C.push_back({v});
        }
    }
    
    /**
     * Process all edges at weight tier w (tier boundary flush).
     * 
     * Steps:
     * 1. Union all buffered edges E_w in the DSU
     * 2. Add resulting component vertex sets to F_C for all touched vertices V_w
     * 3. Clear buffers E_w and V_w
     */
    void flush_tier(
        std::vector<std::pair<Vertex, Vertex>>& E_w,
        VertexSet& V_w,
        DisjointSetUnion& dsu,
        ComponentTracker& tracker,
        std::vector<VertexSet>& F_C
    ) {
        // Step 1: Union all buffered edges
        for (const auto& edge : E_w) {
            Vertex u = edge.first;
            Vertex v = edge.second;
            
            Vertex root_u = dsu.find(u);
            Vertex root_v = dsu.find(v);
            
            if (root_u != root_v) {
                // Perform union in DSU
                Vertex new_root = dsu.union_sets(u, v);
                // Update component tracking
                tracker.union_sets(root_u, root_v, new_root);
            }
        }
        
        // Step 2: Add component vertex sets to F_C for all touched vertices
        // OPTIMIZED: Track by root first to avoid processing same component multiple times
        std::unordered_set<Vertex> processed_roots;
        processed_roots.reserve(V_w.size());
        // OPTIMIZED: Use normalized vectors for deduplication (like Python's frozenset approach)
        using NormalizedComponent = std::vector<Vertex>;
        struct NormalizedHash {
            std::size_t operator()(const NormalizedComponent& vec) const {
                std::size_t hash = 0;
                for (const auto& v : vec) {
                    hash ^= std::hash<Vertex>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
                }
                return hash;
            }
        };
        
        std::unordered_set<NormalizedComponent, NormalizedHash> added_components;
        added_components.reserve(V_w.size());
        added_components.max_load_factor(0.75);
        
        for (const auto& v : V_w) {
            Vertex root = dsu.find(v);
            
            // Skip if we've already processed this root (multiple vertices map to same root)
            if (processed_roots.find(root) != processed_roots.end()) {
                continue;
            }
            processed_roots.insert(root);
            
            const VertexSet& component = tracker.component_vertices(root);
            
            if (component.empty()) continue;
            
            // Convert to normalized form for reliable comparison
            NormalizedComponent normalized(component.begin(), component.end());
            std::sort(normalized.begin(), normalized.end());
            
            // O(1) hash lookup - check if we've already added this exact component
            if (added_components.find(normalized) == added_components.end()) {
                F_C.push_back(component);  // Store the component
                added_components.insert(normalized);
            }
        }
        
        // Step 3: Clear buffers for next tier
        E_w.clear();
        V_w.clear();
    }
    
    /**
     * Helper to canonicalize edge (ensure u < v for consistent ordering)
     */
    std::pair<Vertex, Vertex> canonicalize_edge(const Vertex& u, const Vertex& v) {
        if (u < v) {
            return {u, v};
        }
        return {v, u};
    }
}

//=============================================================================
// Main Algorithm Implementation
//=============================================================================

FC_GU_Result construct_FC_and_GU(
    const VertexList& vertex_list,
    const EdgeList& edge_list
) {
    // Step 1: Sort edges by weight (Kruskal's algorithm requirement)
    // CRITICAL: Use stable_sort WITHOUT secondary key to match Python's sorted() behavior
    // Python's sorted() preserves original order for equal weights, so we must do the same
    EdgeList E_sorted;
    E_sorted.reserve(edge_list.size());
    E_sorted = edge_list;  // Move or copy efficiently
    std::stable_sort(E_sorted.begin(), E_sorted.end(), 
                     [](const Edge& a, const Edge& b) {
                         return std::get<2>(a) < std::get<2>(b);
                     });
    
    // Edge case: no edges at all
    if (E_sorted.empty()) {
        // F_C contains V_G and all singletons
        std::vector<VertexSet> F_C;
        VertexSet V_G(vertex_list.begin(), vertex_list.end());
        F_C.push_back(V_G);
        add_all_singletons(vertex_list, F_C);
        return FC_GU_Result{F_C, std::set<std::pair<Vertex, Vertex>>()};
    }
    
    // Step 2: Initialize data structures
    DisjointSetUnion dsu(vertex_list);
    ComponentTracker tracker(vertex_list);
    
    // F_C starts with the complete vertex set V_G
    std::vector<VertexSet> F_C;
    // More aggressive pre-allocation for large datasets
    size_t estimated_fc_size = std::min(edge_list.size() / 5, size_t(1000000)) + vertex_list.size();
    F_C.reserve(estimated_fc_size);
    VertexSet V_G(vertex_list.begin(), vertex_list.end());
    F_C.push_back(std::move(V_G));
    
    // FIXED: Add all singleton components FIRST (before processing any edges)
    // This represents the state before the first tier: all vertices are isolated
    add_all_singletons(vertex_list, F_C);
    
    // G_U will accumulate edges where endpoints are in different components
    // OPTIMIZED: Use unordered_set for O(1) insertion instead of O(log n) with set
    std::unordered_set<std::pair<Vertex, Vertex>, PairHash> G_U_edges;
    // More aggressive estimate for large datasets
    size_t estimated_gu_size = std::min(edge_list.size() / 5, size_t(5000000));  // Cap at 5M
    G_U_edges.reserve(estimated_gu_size);
    
    // Buffers for tier-based processing
    std::vector<std::pair<Vertex, Vertex>> E_w;  // Edges at current weight tier
    E_w.reserve(1000);  // Larger reserve for datasets with many edges per tier
    VertexSet V_w;  // Vertices touched at current tier
    V_w.reserve(1000);
    
    double w_previous = std::get<2>(E_sorted[0]);  // Weight of first edge
    
    // Step 3: Process edges in weight-sorted order (Kruskal sweep)
    for (const auto& edge : E_sorted) {
        Vertex u = std::get<0>(edge);
        Vertex v = std::get<1>(edge);
        double w = std::get<2>(edge);
        
        // Tier boundary: weight has increased, flush previous tier
        if (w > w_previous) {
            flush_tier(E_w, V_w, dsu, tracker, F_C);
            w_previous = w;
            // Reserve for next tier - larger reserve for big datasets
            E_w.reserve(1000);
            V_w.reserve(1000);
        }
        
        // Check if u and v are in different components
        Vertex root_u = dsu.find(u);
        Vertex root_v = dsu.find(v);
        
        if (root_u != root_v) {
            // Edge connects different components → add to G_U
            G_U_edges.insert(canonicalize_edge(u, v));
            
            // Buffer edge for tier flush
            E_w.push_back({u, v});
            V_w.insert(u);
            V_w.insert(v);
        }
        // else: edge within same component, ignore (cycle in MST)
    }
    
    // Step 4: Final flush for the last weight tier
    flush_tier(E_w, V_w, dsu, tracker, F_C);
    
    // Step 5: OPTIMIZED deduplicate F_C using hash set (O(n) instead of O(n²))
    // CRITICAL FIX: Use sorted vector as key (like Python's frozenset) for reliable comparison
    using NormalizedComponent = std::vector<Vertex>;  // Sorted vertex list for comparison
    
    struct NormalizedHash {
        std::size_t operator()(const NormalizedComponent& vec) const {
            std::size_t hash = 0;
            for (const auto& v : vec) {
                hash ^= std::hash<Vertex>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };
    
    std::vector<VertexSet> F_C_deduplicated;
    std::unordered_set<NormalizedComponent, NormalizedHash> seen;
    seen.reserve(F_C.size());
    seen.max_load_factor(0.75);
    F_C_deduplicated.reserve(F_C.size());
    
    for (const auto& component : F_C) {
        // Convert to normalized form (sorted vector) like Python's frozenset
        NormalizedComponent normalized(component.begin(), component.end());
        std::sort(normalized.begin(), normalized.end());
        
        if (seen.find(normalized) == seen.end()) {
            F_C_deduplicated.push_back(component);
            seen.insert(normalized);
        }
    }
    
    // Sort by size for consistent output (laminarity is preserved by construction)
    std::sort(F_C_deduplicated.begin(), F_C_deduplicated.end(),
              [](const VertexSet& a, const VertexSet& b) {
                  return a.size() < b.size();
              });
    
    // Convert G_U_edges from unordered_set to set for return type compatibility
    std::set<std::pair<Vertex, Vertex>> G_U_set(G_U_edges.begin(), G_U_edges.end());
    
    return FC_GU_Result{F_C_deduplicated, G_U_set};
}

//=============================================================================
// Validation Function
//=============================================================================

bool validate_laminar(const std::vector<VertexSet>& F_C) {
    for (size_t i = 0; i < F_C.size(); ++i) {
        for (size_t j = i + 1; j < F_C.size(); ++j) {
            const VertexSet& X = F_C[i];
            const VertexSet& Y = F_C[j];
            
            // Check intersection
            VertexSet intersection;
            for (const auto& x : X) {
                if (Y.find(x) != Y.end()) {
                    intersection.insert(x);
                }
            }
            
            // If they intersect, one must be a subset of the other
            if (!intersection.empty()) {
                bool X_subset_Y = true;
                for (const auto& x : X) {
                    if (Y.find(x) == Y.end()) {
                        X_subset_Y = false;
                        break;
                    }
                }
                
                bool Y_subset_X = true;
                for (const auto& y : Y) {
                    if (X.find(y) == X.end()) {
                        Y_subset_X = false;
                        break;
                    }
                }
                
                if (!X_subset_Y && !Y_subset_X) {
                    std::cerr << "NOT LAMINAR: Sets intersect but neither contains the other\n";
                    return false;
                }
            }
        }
    }
    
    return true;
}

//=============================================================================
// FC_GU_Result Comparison (for testing)
//=============================================================================

bool FC_GU_Result::operator==(const FC_GU_Result& other) const {
    // Compare F_C sizes
    if (F_C.size() != other.F_C.size()) {
        return false;
    }
    
    // Compare F_C by converting to sorted vectors of sorted vertex lists
    auto normalize_F_C = [](const std::vector<VertexSet>& fc) {
        std::vector<std::vector<Vertex>> normalized;
        for (const auto& set : fc) {
            std::vector<Vertex> vec(set.begin(), set.end());
            std::sort(vec.begin(), vec.end());
            normalized.push_back(vec);
        }
        std::sort(normalized.begin(), normalized.end());
        return normalized;
    };
    
    auto F_C_norm = normalize_F_C(F_C);
    auto other_F_C_norm = normalize_F_C(other.F_C);
    
    if (F_C_norm != other_F_C_norm) {
        return false;
    }
    
    // Compare G_U_edges
    if (G_U_edges != other.G_U_edges) {
        return false;
    }
    
    return true;
}

} // namespace phylo

