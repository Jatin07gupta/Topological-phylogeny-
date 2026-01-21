#ifndef MLVMST_H
#define MLVMST_H

#include "common_types.h"
#include "algorithm2.h"  // For DisjointSetUnion
#include <vector>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <algorithm>

namespace phylo {

    /**
     * Result structure for Algorithm 3 (MLVMST)
     */
    struct MLVMSTResult {
        std::vector<Vertex> order;              // Total order <* used for VMST construction
        std::map<Vertex, int> delta_max;        // Per-vertex delta_max scores
        std::set<std::pair<Vertex, Vertex>> M_star_edges;  // Edges of M* (minimum leaves VMST)
        int leaf_count;                         // Number of leaves in M*
        
        bool operator==(const MLVMSTResult& other) const;
    };

    /**
     * Build adjacency list from undirected edge set.
     * 
     * @param edges Set of undirected edges (canonicalized)
     * @return Map from vertex to set of neighbors
     */
    std::unordered_map<Vertex, VertexSet> build_adjacency(
        const std::set<std::pair<Vertex, Vertex>>& edges
    );

    /**
     * Compute delta_max(i) for each vertex i exactly as in Algorithm 3.
     * 
     * Algorithm (from paper Lines L999-L1009):
     * For each vertex i:
     *   1. N_i ← neighbors of i in G_U
     *   2. delta_max(i) ← 0
     *   3. For each C in F_C^≥ (ordered by decreasing size):
     *      - If (C ∩ N_i ≠ ∅) and (i ∉ C):
     *        * delta_max(i) ← delta_max(i) + 1
     *        * N_i ← N_i \ C
     * 
     * @param vertices List of all vertices
     * @param F_C Common laminar family (from Algorithm 2)
     * @param G_U_adj Adjacency list of G_U (MST-union graph)
     * @return Map from vertex to delta_max score
     */
    std::map<Vertex, int> compute_delta_max(
        const VertexList& vertices,
        const std::vector<VertexSet>& F_C,
        const std::unordered_map<Vertex, VertexSet>& G_U_adj
    );

    /**
     * Build total order <* based on delta_max scores.
     * 
     * Note: Empirically, REVERSED ordering (larger delta_max first) gives better results
     * than the paper's stated nondecreasing order. Default is "decreasing" (larger first).
     * 
     * @param vertices List of all vertices
     * @param delta_max delta_max scores for each vertex
     * @param ordering Ordering direction: "decreasing" (larger first, default, best) 
     *                 or "increasing" (smaller first, paper's stated ordering)
     * @return List of vertices in total order <*
     */
    VertexList build_total_order(
        const VertexList& vertices,
        const std::map<Vertex, int>& delta_max,
        const std::string& ordering = "decreasing"
    );

    /**
     * Run Algorithm 1 (VMST) with given vertex order.
     * 
     * This is Kruskal's algorithm with stable tie-breaking using vertex order.
     * When two edges have the same weight, prefer the edge with lexicographically
     * smaller vertex pair according to the given order.
     * 
     * @param vertices List of vertex identifiers
     * @param edges List of weighted edges (u, v, w)
     * @param vertex_order Total order <* over vertices
     * @return Set of edges in the VMST (canonicalized, undirected)
     */
    std::set<std::pair<Vertex, Vertex>> build_vmst_from_order(
        const VertexList& vertices,
        const EdgeList& edges,
        const VertexList& vertex_order
    );

    /**
     * Count the number of leaves (degree-1 vertices) in a tree.
     * 
     * @param edges Set of tree edges
     * @param vertices List of all vertices
     * @return Number of leaves in the tree
     */
    int count_leaves(
        const std::set<std::pair<Vertex, Vertex>>& edges,
        const VertexList& vertices
    );

    /**
     * Implements Algorithm 3 (MLVMST) - Construct a minimum leaves VMST.
     * 
     * This algorithm finds a VMST with the minimum number of leaves among all
     * possible VMSTs by computing a special vertex ordering based on delta_max scores.
     * 
     * Algorithm 3 (Lines L991-L1011):
     * 1. Order F_C by decreasing set size → F_C^≥
     * 2. Build adjacency of G_U
     * 3. For each vertex i, compute delta_max(i)
     * 4. Build total order <* by delta_max (default: nonincreasing)
     * 5. Run VMST (Algorithm 1) with this order to obtain M*
     * 
     * @param vertices List of vertex identifiers
     * @param edges List of weighted edges (u, v, w)
     * @param F_C Common laminar family from Algorithm 2
     * @param G_U_edges MST-union graph edges from Algorithm 2
     * @param ordering Ordering direction: "decreasing" (larger first, default, best) 
     *                 or "increasing" (smaller first, paper's stated ordering)
     * @return MLVMSTResult containing order, delta_max, M* edges, and leaf count
     */
    MLVMSTResult build_mlvmst(
        const VertexList& vertices,
        const EdgeList& edges,
        const std::vector<VertexSet>& F_C,
        const std::set<std::pair<Vertex, Vertex>>& G_U_edges,
        const std::string& ordering = "decreasing"
    );

} // namespace phylo

#endif // MLVMST_H

