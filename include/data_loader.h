#pragma once

#include "common_types.h"
#include <string>

namespace phylo {

/**
 * Load test data from file.
 * 
 * File format:
 * Line 1: number of vertices (n)
 * Line 2: number of edges (m)
 * Lines 3 to n+2: vertex names (one per line)
 * Lines n+3 onwards: edges (format: "u v weight", one per line)
 * 
 * Returns true on success, false on error.
 */
bool load_test_data(const std::string& filename, VertexList& vertices, EdgeList& edges);

} // namespace phylo

