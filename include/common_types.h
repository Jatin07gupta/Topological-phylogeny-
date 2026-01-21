#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <set>

namespace phylo {

// Type aliases for clarity
using Vertex = std::string;
using Edge = std::tuple<Vertex, Vertex, double>;  // (u, v, weight)
using VertexSet = std::unordered_set<Vertex>;
using VertexList = std::vector<Vertex>;
using EdgeList = std::vector<Edge>;

} // namespace phylo

