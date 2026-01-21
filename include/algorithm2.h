#pragma once

#include "common_types.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <tuple>

namespace phylo {

/**
 * Result from Algorithm 2: Construct F_C and G_U
 * 
 * F_C: Common laminar family (list of vertex sets)
 * G_U_edges: MST-union graph edges (set of vertex pairs)
 */
struct FC_GU_Result {
    std::vector<VertexSet> F_C;  // Common laminar family
    std::set<std::pair<Vertex, Vertex>> G_U_edges;  // MST-union graph edges (canonicalized)
    
    // Comparison operators for testing
    bool operator==(const FC_GU_Result& other) const;
};

/**
 * Disjoint Set Union (Union-Find) with path compression and union by rank.
 * 
 * Tracks component membership for vertices during Kruskal's algorithm sweep.
 * O(α(n)) amortized complexity per operation.
 */
class DisjointSetUnion {
public:
    DisjointSetUnion(const VertexList& vertices);
    
    /**
     * Find the root representative of vertex v's component.
     * Uses path compression for amortized O(α(n)) complexity.
     */
    Vertex find(const Vertex& v);
    
    /**
     * Union the components containing u and v.
     * Uses union by rank for balanced trees.
     * Returns the new root representative.
     */
    Vertex union_sets(const Vertex& u, const Vertex& v);

private:
    std::unordered_map<Vertex, Vertex> parent_;
    std::unordered_map<Vertex, int> rank_;
};

/**
 * Track which vertices belong to each component in the DSU.
 * 
 * Maintains vertex sets for each root representative to avoid
 * repeatedly traversing the DSU to find component members.
 */
class ComponentTracker {
public:
    ComponentTracker(const VertexList& vertices);
    
    /**
     * Update component tracking after DSU union operation.
     * Merge vertex sets from old roots into new root.
     */
    void union_sets(const Vertex& old_root_u, const Vertex& old_root_v, const Vertex& new_root);
    
    /**
     * Get all vertices in the component represented by root.
     */
    const VertexSet& component_vertices(const Vertex& root) const;

private:
    std::unordered_map<Vertex, VertexSet> component_sets_;
    static const VertexSet empty_set_;  // For returning empty set when not found
};

/**
 * Algorithm 2: Construct F_C (common laminar family) and G_U (MST-union graph).
 * 
 * FIXED VERSION: Properly includes all singleton components in F_C.
 * 
 * F_C Structure (laminar family):
 * 1. V_G - the complete vertex set (always included)
 * 2. {v} - singleton for every vertex v (before any edges processed)
 * 3. All nontrivial components formed at each weight tier boundary
 * 
 * The algorithm performs a Kruskal-by-tiers sweep:
 * - Edges are processed in sorted order by weight
 * - At each weight tier boundary, components are flushed to F_C
 * - An edge (u,v) is added to G_U if u,v are in different components
 * 
 * Complexity: O(E log E) for sorting + O(E α(V)) for DSU operations
 * 
 * Args:
 *   vertex_list: List of all vertices in graph G
 *   edge_list: List of edges as (u, v, weight) tuples
 * 
 * Returns:
 *   FC_GU_Result with:
 *     - F_C: Common laminar family (list of vertex sets)
 *     - G_U_edges: MST-union graph edges (set of vertex pairs)
 */
FC_GU_Result construct_FC_and_GU(
    const VertexList& vertex_list,
    const EdgeList& edge_list
);

/**
 * Validate that F_C forms a laminar family:
 * For any two sets X, Y in F_C, either:
 * - X ⊆ Y, or
 * - Y ⊆ X, or
 * - X ∩ Y = ∅
 * 
 * Returns true if F_C is laminar, false otherwise.
 */
bool validate_laminar(const std::vector<VertexSet>& F_C);

} // namespace phylo

