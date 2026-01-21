#ifndef CLNJ_H
#define CLNJ_H

#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <Eigen/Dense>

namespace phylo {

/**
 * Structure to track hidden node information.
 * Maps hidden_node_id -> {anchor: observed_node, dist_to_anchor: distance}
 */
struct HiddenNodeInfo {
    int anchor;           // Observed node that serves as anchor
    double dist_to_anchor; // Distance from hidden node to anchor
    
    HiddenNodeInfo() : anchor(-1), dist_to_anchor(0.0) {}
    HiddenNodeInfo(int a, double d) : anchor(a), dist_to_anchor(d) {}
};

/**
 * Structure for an edge in the local NJ result.
 */
struct LocalEdge {
    int node1;
    int node2;
    double distance;
    
    LocalEdge() : node1(-1), node2(-1), distance(0.0) {}
    LocalEdge(int n1, int n2, double d) : node1(n1), node2(n2), distance(d) {}
};

/**
 * Structure for a newly created hidden node.
 */
struct NewHiddenNode {
    int id;
    int anchor;
    double dist_to_anchor;
    
    NewHiddenNode() : id(-1), anchor(-1), dist_to_anchor(0.0) {}
    NewHiddenNode(int i, int a, double d) : id(i), anchor(a), dist_to_anchor(d) {}
};

/**
 * Result structure for CLNJ algorithm.
 */
struct CLNJResult {
    std::unordered_map<int, std::unordered_set<int>> adjacency;  // node_id -> set of neighbors
    std::map<std::pair<int, int>, double> edge_weights;          // Edge weights: (min, max) -> distance
    std::unordered_map<int, HiddenNodeInfo> hidden_info;         // hidden_id -> HiddenNodeInfo
    
    CLNJResult() = default;
};

/**
 * Compute distance between any two nodes (observed or hidden).
 * 
 * @param i Node index (observed if i < m, hidden otherwise)
 * @param j Node index (observed if j < m, hidden otherwise)
 * @param distance_matrix Original m×m observed distance matrix
 * @param hidden_info Maps hidden_node_id -> HiddenNodeInfo
 * @param m Number of observed nodes
 * @return Distance between i and j
 */
double compute_distance_with_hidden(
    int i,
    int j,
    const Eigen::MatrixXd& distance_matrix,
    const std::unordered_map<int, HiddenNodeInfo>& hidden_info,
    int m
);

/**
 * Local Neighbor Joining algorithm.
 * 
 * This performs NJ on a local set of nodes, potentially creating hidden nodes.
 * 
 * @param nodes List of node indices (can be observed or hidden)
 * @param distance_matrix Original m×m observed distance matrix
 * @param hidden_info Global hidden info (will be updated locally during NJ)
 * @param m Number of observed nodes
 * @param threshold Threshold for creating hidden nodes (default: -log(0.9))
 * @return Pair of (local_edges, new_hidden_nodes)
 */
std::pair<std::vector<LocalEdge>, std::vector<NewHiddenNode>> neighbor_joining_local(
    const std::vector<int>& nodes,
    const Eigen::MatrixXd& distance_matrix,
    const std::unordered_map<int, HiddenNodeInfo>& hidden_info,
    int m,
    int& next_hidden_id,  // Pass by reference to update globally
    double threshold = -std::log(0.9)
);

/**
 * Clean CLNJ implementation using adjacency list approach.
 * 
 * This is the main CLNJ algorithm that processes internal nodes in the MST,
 * runs local NJ on their neighborhoods, and updates the adjacency structure.
 * 
 * @param distance_matrix m×m observed distance matrix
 * @param mst_adjacency m×m MLVMST adjacency matrix (0/1 values)
 * @param verbose Print progress messages
 * @return CLNJResult containing final adjacency and hidden node info
 */
CLNJResult CLNJ_clean(
    const Eigen::MatrixXd& distance_matrix,
    const Eigen::MatrixXi& mst_adjacency,
    double threshold = -std::log(0.9),
    bool verbose = true
);

} // namespace phylo

#endif // CLNJ_H


