#ifndef MLVMST_BORUVKA_H
#define MLVMST_BORUVKA_H

#include "common_types.h"
#include "algorithm2.h"  // For DisjointSetUnion
#include <vector>
#include <set>
#include <unordered_map>

namespace phylo {

    /**
     * Hybrid Borůvka–VMST constructor with the same semantics as build_vmst_from_order.
     * 
     * This function constructs the VMST using a Borůvka-style phased approach while
     * maintaining logical equivalence to the Kruskal-based build_vmst_from_order.
     * It produces exactly the same VMST edges for the same input.
     * 
     * The algorithm works in phases:
     * 1. For each component (DSU root), find its best outgoing edge according to
     *    the VMST key (weight, min_rank, max_rank).
     * 2. Collect and deduplicate candidate edges from all components.
     * 3. Sort candidates globally by the VMST key.
     * 4. Process candidates Kruskal-style with DSU to merge components.
     * 5. Repeat until one component remains or MST is complete.
     * 
     * The per-component edge scanning step (step 1) is structured to be
     * parallel-friendly, though this implementation is sequential.
     * 
     * @param vertices List of vertex identifiers
     * @param edges List of weighted edges (u, v, w)
     * @param vertex_order Total order <* over vertices
     * @return Set of edges in the VMST (canonicalized, undirected)
     * 
     * Note:
     * This function is logically equivalent to build_vmst_from_order and
     * should produce identical results for the same input.
     */
    std::set<std::pair<Vertex, Vertex>> build_vmst_from_order_hybrid(
        const VertexList& vertices,
        const EdgeList& edges,
        const VertexList& vertex_order
    );

    /**
     * Parallelized hybrid Borůvka–VMST constructor with the same semantics as build_vmst_from_order.
     * 
     * This is the parallel version that uses OpenMP to parallelize the edge scanning phase.
     * The parallel part: finding best edge per component (edge scanning loop)
     * The sequential parts: DSU operations, sorting, candidate processing
     * 
     * @param vertices List of vertex identifiers
     * @param edges List of weighted edges (u, v, w)
     * @param vertex_order Total order <* over vertices
     * @param num_threads Number of threads to use (0 = use OpenMP default)
     * @return Set of edges in the VMST (canonicalized, undirected)
     */
    std::set<std::pair<Vertex, Vertex>> build_vmst_from_order_hybrid_parallel(
        const VertexList& vertices,
        const EdgeList& edges,
        const VertexList& vertex_order,
        int num_threads = 0
    );

} // namespace phylo

#endif // MLVMST_BORUVKA_H

