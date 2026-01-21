/**
 * CLNJ Fresh Implementation - Matching Python Exactly
 * Based on clnj_clean_implementation_sorted.py
 * Updated to use Eigen and include edge_weight tracking
 * 
 * To enable debug logging, compile with: -DCLNJ_DEBUG
 */

#include "clnj.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <iomanip>

namespace phylo {

/**
 * Compute distance between any two nodes (observed or hidden)
 * Matches Python: compute_distance_with_hidden()
 */
double compute_distance_with_hidden(
    int i,
    int j,
    const Eigen::MatrixXd& distance_matrix,
    const std::unordered_map<int, HiddenNodeInfo>& hidden_info,
    int m
) {
    // Both observed
    if (i < m && j < m) {
        return distance_matrix(i, j);
    }
    
    // i is hidden, j is observed
    if (i >= m && j < m) {
        auto it = hidden_info.find(i);
        if (it == hidden_info.end()) {
            std::cerr << "ERROR: Hidden node " << i << " not found in hidden_info" << std::endl;
            return 0.0;  // Return 0.0 instead of infinity to avoid breaking Q-matrix minimum finding
        }
        int anchor_i = it->second.anchor;
        double dist_i_anchor = it->second.dist_to_anchor;
        return distance_matrix(anchor_i, j) - dist_i_anchor;
    }
    
    // i is observed, j is hidden
    if (i < m && j >= m) {
        auto it = hidden_info.find(j);
        if (it == hidden_info.end()) {
            std::cerr << "ERROR: Hidden node " << j << " not found in hidden_info" << std::endl;
            return 0.0;  // Return 0.0 instead of infinity to avoid breaking Q-matrix minimum finding
        }
        int anchor_j = it->second.anchor;
        double dist_j_anchor = it->second.dist_to_anchor;
        return distance_matrix(i, anchor_j) - dist_j_anchor;
    }
    
    // Both hidden
    auto it_i = hidden_info.find(i);
    auto it_j = hidden_info.find(j);
    if (it_i == hidden_info.end() || it_j == hidden_info.end()) {
        std::cerr << "ERROR: Hidden node not found" << std::endl;
        return 0.0;  // Return 0.0 instead of infinity to avoid breaking Q-matrix minimum finding
    }
    
    int anchor_i = it_i->second.anchor;
    double dist_i_anchor = it_i->second.dist_to_anchor;
    int anchor_j = it_j->second.anchor;
    double dist_j_anchor = it_j->second.dist_to_anchor;
    
    return distance_matrix(anchor_i, anchor_j) - dist_i_anchor - dist_j_anchor;
}

/**
 * Local Neighbor Joining
 * Matches Python: neighbor_joining_local() exactly
 */
std::pair<std::vector<LocalEdge>, std::vector<NewHiddenNode>> neighbor_joining_local(
    const std::vector<int>& nodes,
    const Eigen::MatrixXd& distance_matrix,
    const std::unordered_map<int, HiddenNodeInfo>& hidden_info,
    int m,
    int& next_hidden_id,  // Pass by reference - will be updated
    double threshold,
    int clnj_iteration  // DEBUG: Track which CLNJ iteration this belongs to
) {
    int n_local = static_cast<int>(nodes.size());
    
    std::vector<LocalEdge> local_edges;
    std::vector<NewHiddenNode> new_hidden;
    
    if (n_local < 2) {
        return {local_edges, new_hidden};
    }
    
    if (n_local == 2) {
        double d = compute_distance_with_hidden(nodes[0], nodes[1], distance_matrix, hidden_info, m);
        local_edges.push_back(LocalEdge(nodes[0], nodes[1], d));
        return {local_edges, new_hidden};
    }
    
    // Build local distance matrix
    Eigen::MatrixXd D(n_local, n_local);
    for (int i = 0; i < n_local; i++) {
        for (int j = 0; j < n_local; j++) {
            D(i, j) = compute_distance_with_hidden(nodes[i], nodes[j], distance_matrix, hidden_info, m);
        }
    }
    
    // Track active nodes (indices in local matrix)
    std::vector<int> active;
    for (int i = 0; i < n_local; i++) {
        active.push_back(i);
    }
    
    // Map local index to global node ID
    std::unordered_map<int, int> node_map;
    for (int i = 0; i < n_local; i++) {
        node_map[i] = nodes[i];
    }
    
    // Create local copy of hidden_info that we'll update as we create new hidden nodes
    std::unordered_map<int, HiddenNodeInfo> local_hidden_info = hidden_info;
    
    // DEBUG: Track NJ iterations
    int nj_iteration = 0;
    
    // NJ iterations until 2 nodes remain
    while (static_cast<int>(active.size()) > 2) {
        nj_iteration++;
        int n_active = static_cast<int>(active.size());
        
        // Compute Q-matrix: Q(i,j) = (n-2)*d(i,j) - sum_k(d(i,k)) - sum_k(d(j,k))
        Eigen::MatrixXd Q = Eigen::MatrixXd::Constant(n_active, n_active, std::numeric_limits<double>::infinity());
        std::vector<double> row_sums(n_active, 0.0);
        
        // Compute row sums
        for (size_t idx_i = 0; idx_i < active.size(); idx_i++) {
            int i = active[idx_i];
            for (size_t idx_j = 0; idx_j < active.size(); idx_j++) {
                int j = active[idx_j];
                if (i != j) {
                    row_sums[idx_i] += D(i, j);
                }
            }
        }
        
        // Compute Q values
        for (size_t idx_i = 0; idx_i < active.size(); idx_i++) {
            int i = active[idx_i];
            for (size_t idx_j = 0; idx_j < active.size(); idx_j++) {
                int j = active[idx_j];
                if (i != j) {
                    Q(idx_i, idx_j) = (n_active - 2) * D(i, j) - row_sums[idx_i] - row_sums[idx_j];
                }
            }
        }
        
        // Find minimum Q(i,j) - matching Python's np.argmin behavior EXACTLY
        // np.argmin flattens in row-major (C-order) and returns FIRST occurrence of minimum
        // Step 1: Find the minimum value
        double min_q = std::numeric_limits<double>::infinity();
        for (int idx_i = 0; idx_i < n_active; idx_i++) {
            for (int idx_j = 0; idx_j < n_active; idx_j++) {
                if (Q(idx_i, idx_j) < min_q) {
                    min_q = Q(idx_i, idx_j);
                }
            }
        }
        
        // Step 2: Find FIRST occurrence of min_q in row-major order (matching np.argmin)
        // This is critical for tie-breaking: np.argmin returns first occurrence
        // CRITICAL: Use tighter epsilon for floating-point comparison
        // Python uses 64-bit floats (numpy.float64), so we need to match that precision
        // Using 1e-12 instead of 1e-10 for better precision matching with Python
        const double epsilon = 1e-12;
        int min_flat_idx = -1;
        for (int idx_i = 0; idx_i < n_active; idx_i++) {
            for (int idx_j = 0; idx_j < n_active; idx_j++) {
                // Use epsilon for floating-point comparison to handle precision issues
                if (std::abs(Q(idx_i, idx_j) - min_q) < epsilon) {
                    min_flat_idx = idx_i * n_active + idx_j;
                    goto found_first_min;  // Exit nested loops at first occurrence
                }
            }
        }
        found_first_min:
        
        // Convert flat index to row/column indices (matching Python line 108-109)
        int min_i_local = min_flat_idx / n_active;
        int min_j_local = min_flat_idx % n_active;
        
        // Ensure i < j (matching Python line 112-113)
        if (min_i_local > min_j_local) {
            std::swap(min_i_local, min_j_local);
        }
        
        int min_i = active[min_i_local];
        int min_j = active[min_j_local];
        
        // Compute branch lengths using NJ formula
        double d_ij = D(min_i, min_j);
        double sum_i = 0.0, sum_j = 0.0;
        
        // DEBUG: Check for inf/NaN values before summing
        bool has_inf_sum_i = false, has_inf_sum_j = false;
        for (int k : active) {
            double d_ik = D(min_i, k);
            double d_jk = D(min_j, k);
            if (std::isinf(d_ik) || std::isnan(d_ik)) {
                has_inf_sum_i = true;
                std::cerr << "ERROR: D(" << min_i << "," << k << ") = " << d_ik << " (inf/nan)" << std::endl;
            }
            if (std::isinf(d_jk) || std::isnan(d_jk)) {
                has_inf_sum_j = true;
                std::cerr << "ERROR: D(" << min_j << "," << k << ") = " << d_jk << " (inf/nan)" << std::endl;
            }
            sum_i += d_ik;
            sum_j += d_jk;
        }
        
        if (has_inf_sum_i || has_inf_sum_j) {
            std::cerr << "ERROR: sum_i or sum_j contains inf/nan. D matrix size: " << D.rows() << "x" << D.cols() << std::endl;
            std::cerr << "  active list: ";
            for (int k : active) std::cerr << k << " ";
            std::cerr << std::endl;
            std::cerr << "  min_i=" << min_i << ", min_j=" << min_j << std::endl;
        }
        
        double delta_i, delta_j;
        if (n_active > 2) {
            delta_i = 0.5 * d_ij + (sum_i - sum_j) / (2.0 * (n_active - 2));
        } else {
            delta_i = 0.5 * d_ij;
        }
        delta_j = d_ij - delta_i;
        
        // Ensure non-negative
        delta_i = std::max(0.0, delta_i);
        delta_j = std::max(0.0, delta_j);
        
        // DEBUG: Log threshold comparison
        // Note: clnj_iteration will be set by CLNJ_clean when calling neighbor_joining_local
        // DIVERGENCE ANALYSIS: Commented out - issue fixed with stable sort in Python
        std::cout << std::fixed << std::setprecision(17);
        std::cout << "[DEBUG] clnj_iter=" << clnj_iteration
                  << " nj_iter=" << nj_iteration
                  << " n_active=" << n_active
                  << " min_i=" << min_i << " min_j=" << min_j
                  << " min_q=" << min_q
                  << " d_ij=" << d_ij
                  << " sum_i=" << sum_i << " sum_j=" << sum_j
                  << " delta_i=" << delta_i 
                  << " delta_j=" << delta_j
                  << " threshold=" << threshold
                  << " delta_i_gt=" << (delta_i > threshold)
                  << " delta_j_gt=" << (delta_j > threshold)
                  << " create_hidden=" << (delta_i > threshold && delta_j > threshold) << std::endl;
        
        // Check if we should create hidden node
        if (delta_i > threshold && delta_j > threshold) {
            // Create hidden parent node
            int u_id = next_hidden_id;
            // DEBUG: Log hidden node creation with ID
            std::cout << "[HIDDEN_NODE_CREATED] clnj_iter=" << clnj_iteration
                      << " nj_iter=" << nj_iteration
                      << " u_id=" << u_id
                      << " next_hidden_id_before=" << next_hidden_id
                      << " min_i=" << node_map[min_i]
                      << " min_j=" << node_map[min_j] << std::endl;
            next_hidden_id += 1;
            
            // Add edges
            local_edges.push_back(LocalEdge(node_map[min_i], u_id, delta_i));
            local_edges.push_back(LocalEdge(node_map[min_j], u_id, delta_j));
            
            // Track hidden node - anchor must be an observed node (< m)
            // If min_i is observed, use it; otherwise use its anchor
            int anchor_node;
            double dist_to_anchor;
            
            if (node_map[min_i] < m) {
                anchor_node = node_map[min_i];
                dist_to_anchor = delta_i;
            } else if (node_map[min_j] < m) {
                anchor_node = node_map[min_j];
                dist_to_anchor = delta_j;
            } else {
                // Both are hidden, use min_i's anchor
                anchor_node = local_hidden_info[node_map[min_i]].anchor;
                dist_to_anchor = local_hidden_info[node_map[min_i]].dist_to_anchor + delta_i;
            }
            
            NewHiddenNode hidden_node(u_id, anchor_node, dist_to_anchor);
            new_hidden.push_back(hidden_node);
            
            // Add to local hidden info so we can use it for future nodes in this local NJ
            local_hidden_info[u_id] = HiddenNodeInfo(anchor_node, dist_to_anchor);
            
            // Add new node to mapping - use next available local index
            int u_local = static_cast<int>(D.rows());
            node_map[u_local] = u_id;
            
            // Expand distance matrix
            Eigen::MatrixXd D_new(D.rows() + 1, D.cols() + 1);
            D_new.topLeftCorner(D.rows(), D.cols()) = D;
            
            // CRITICAL FIX: Explicitly initialize diagonal element for new node to 0
            // Eigen should do this automatically, but we ensure it's set correctly
            D_new(u_local, u_local) = 0.0;
            
            // Compute distances from u to other active nodes
            for (int k : active) {
                if (k != min_i && k != min_j) {
                    double d_uk = 0.5 * (D(min_i, k) + D(min_j, k) - d_ij);
                    d_uk = std::max(0.0, d_uk);
                    D_new(u_local, k) = d_uk;
                    D_new(k, u_local) = d_uk;
                }
            }
            
            D = D_new;
            
            // Remove min_i and min_j, add u
            // Matching Python's list.remove() - removes first occurrence
            auto it_i = std::find(active.begin(), active.end(), min_i);
            if (it_i != active.end()) {
                active.erase(it_i);
            }
            auto it_j = std::find(active.begin(), active.end(), min_j);
            if (it_j != active.end()) {
                active.erase(it_j);
            }
            active.push_back(u_local);
        } else {
            // Don't create hidden node - connect directly
            local_edges.push_back(LocalEdge(node_map[min_i], node_map[min_j], d_ij));
            
            // Merge min_j into min_i
            for (int k : active) {
                if (k != min_i && k != min_j) {
                    double d_uk = 0.5 * (D(min_i, k) + D(min_j, k) - d_ij);
                    d_uk = std::max(0.0, d_uk);
                    D(min_i, k) = d_uk;
                    D(k, min_i) = d_uk;
                }
            }
            
            // Remove min_j from active (matching Python's list.remove())
            auto it_j = std::find(active.begin(), active.end(), min_j);
            if (it_j != active.end()) {
                active.erase(it_j);
            }
        }
    }
    
    // Final edge between last 2 nodes
    if (active.size() == 2) {
        int i = active[0];
        int j = active[1];
        double d = D(i, j);
        local_edges.push_back(LocalEdge(node_map[i], node_map[j], d));
    }
    
    return {local_edges, new_hidden};
}

/**
 * CLNJ Clean Implementation
 * Matches Python: CLNJ_clean() exactly
 * Includes edge_weight tracking
 */
CLNJResult CLNJ_clean(
    const Eigen::MatrixXd& distance_matrix,
    const Eigen::MatrixXi& mst_adjacency,
    double threshold,
    bool verbose
) {
    int m = static_cast<int>(distance_matrix.rows());
    
    if (verbose) {
        std::cout << "CLNJ Clean (SORTED): Processing " << m << " observed nodes" << std::endl;
    }
    
    // Initialize global adjacency list
    std::unordered_map<int, std::unordered_set<int>> adjacency;
    std::map<std::pair<int, int>, double> edge_weights;
    
    for (int i = 0; i < m; i++) {
        adjacency[i] = std::unordered_set<int>();
    }
    
    // Populate from MLVMST
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < m; j++) {
            if (mst_adjacency(i, j) > 0) {
                adjacency[i].insert(j);
                adjacency[j].insert(i);
                edge_weights[{std::min(i, j), std::max(i, j)}] = distance_matrix(i, j);
            }
        }
    }
    
    // Track hidden nodes
    std::unordered_map<int, HiddenNodeInfo> hidden_info;
    int next_hidden_id = m;
    
    // Find internal nodes in ORIGINAL MLVMST
    std::vector<std::pair<int, int>> internal_nodes;  // (degree, node)
    for (int i = 0; i < m; i++) {
        int degree = 0;
        for (int j = 0; j < m; j++) {
            if (mst_adjacency(i, j) > 0) {
                degree++;
            }
        }
        if (degree > 1) {
            internal_nodes.push_back({degree, i});
        }
    }
    
    // Sort by degree descending (matching Python's np.argsort(-degree[internal_nodes]))
    // np.argsort is stable and preserves order of equal elements
    // When degrees are equal, preserve original order (by node index)
    std::stable_sort(internal_nodes.begin(), internal_nodes.end(),
                     [](const auto& a, const auto& b) {
                         if (a.first != b.first) {
                             return a.first > b.first;  // Descending by degree
                         }
                         return a.second < b.second;  // Secondary: ascending by node index (preserves original order)
                     });
    
    if (verbose) {
        std::cout << "Found " << internal_nodes.size() << " internal nodes to process" << std::endl;
    }
    
    // Process each internal node
    int clnj_iteration = 0;
    for (size_t iter_num = 0; iter_num < internal_nodes.size(); iter_num++) {
        clnj_iteration++;
        int center = internal_nodes[iter_num].second;
        
        // Get CURRENT neighborhood from adjacency
        std::unordered_set<int> neighbors = adjacency[center];
        std::unordered_set<int> S_v = {center};
        S_v.insert(neighbors.begin(), neighbors.end());
        
        // MODIFIED: Use sorted order instead of hash-based order for deterministic matching
        // Matching Python: S_v_list = sorted(list(S_v))
        // CRITICAL: Python's sorted() is STABLE (Timsort), so we must use std::stable_sort
        std::vector<int> S_v_list(S_v.begin(), S_v.end());
        std::stable_sort(S_v_list.begin(), S_v_list.end());
        
        if (verbose) {
            std::cout << "  Iteration " << (iter_num + 1) << "/" << internal_nodes.size()
                      << ": Processing node " << center
                      << " (current degree " << neighbors.size()
                      << ", neighborhood size " << S_v.size() << ")" << std::endl;
        }
        
        // Calculate next_hidden_id (matching Python: next_hidden_id = m + len(hidden_info))
        next_hidden_id = m + static_cast<int>(hidden_info.size());
        
        // DEBUG: Log next_hidden_id initialization
        std::cout << "[NEXT_HIDDEN_ID] clnj_iter=" << clnj_iteration
                  << " center=" << center
                  << " hidden_info_size=" << hidden_info.size()
                  << " next_hidden_id=" << next_hidden_id << std::endl;
        
        // Run local NJ
        auto [local_edges, new_hidden] = neighbor_joining_local(
            S_v_list, distance_matrix, hidden_info, m, next_hidden_id, threshold, clnj_iteration
        );
        
        if (verbose && !new_hidden.empty()) {
            std::cout << "    Created " << new_hidden.size() << " new hidden node(s)" << std::endl;
        }
        
        // Update hidden_info
        for (const auto& h : new_hidden) {
            hidden_info[h.id] = HiddenNodeInfo(h.anchor, h.dist_to_anchor);
            adjacency[h.id] = std::unordered_set<int>();
        }
        
        // SPLICE: Remove old edges within S_v, add new edges
        // Step 1: Remove all edges within S_v
        // Matching Python: nested loop handles both directions naturally
        for (int node1 : S_v) {
            for (int node2 : S_v) {
                if (node1 != node2 && adjacency.count(node1) && adjacency[node1].count(node2)) {
                    adjacency[node1].erase(node2);
                    auto edge_key = std::make_pair(std::min(node1, node2), std::max(node1, node2));
                    edge_weights.erase(edge_key);
                }
            }
        }
        
        // Step 2: Add all edges from local NJ result
        for (const auto& edge : local_edges) {
            int node1 = edge.node1;
            int node2 = edge.node2;
            double dist = edge.distance;
            
            // Ensure both nodes are in adjacency
            if (adjacency.find(node1) == adjacency.end()) {
                adjacency[node1] = std::unordered_set<int>();
            }
            if (adjacency.find(node2) == adjacency.end()) {
                adjacency[node2] = std::unordered_set<int>();
            }
            
            adjacency[node1].insert(node2);
            adjacency[node2].insert(node1);
            auto edge_key = std::make_pair(std::min(node1, node2), std::max(node1, node2));
            edge_weights[edge_key] = dist;
        }
    }
    
    if (verbose) {
        int total_nodes = static_cast<int>(adjacency.size());
        int total_edges = static_cast<int>(edge_weights.size());
        int num_hidden = static_cast<int>(hidden_info.size());
        std::cout << "CLNJ complete: " << total_nodes << " total nodes ("
                  << m << " observed + " << num_hidden << " hidden)" << std::endl;
        std::cout << "  Total edges: " << total_edges << std::endl;
    }
    
    // Build result
    CLNJResult result;
    result.adjacency = adjacency;
    result.edge_weights = edge_weights;
    result.hidden_info = hidden_info;
    
    return result;
}

} // namespace phylo
