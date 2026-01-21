#include "algorithm2.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>

namespace phylo {

bool load_test_data(const std::string& filename, VertexList& vertices, EdgeList& edges) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file: " << filename << std::endl;
        return false;
    }
    
    int n_vertices, n_edges;
    
    // Read number of vertices
    file >> n_vertices;
    if (n_vertices <= 0) {
        std::cerr << "ERROR: Invalid number of vertices: " << n_vertices << std::endl;
        return false;
    }
    
    // Read number of edges
    file >> n_edges;
    if (n_edges < 0) {
        std::cerr << "ERROR: Invalid number of edges: " << n_edges << std::endl;
        return false;
    }
    
    // Skip to next line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    // Read vertex names
    vertices.clear();
    vertices.reserve(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        std::string vertex;
        std::getline(file, vertex);
        if (vertex.empty()) {
            std::cerr << "ERROR: Empty vertex name at line " << (i + 3) << std::endl;
            return false;
        }
        vertices.push_back(vertex);
    }
    
    // Read edges
    edges.clear();
    edges.reserve(n_edges);
    
    // Optimize for large files - read in chunks
    std::string line;
    line.reserve(256);  // Pre-allocate line buffer
    
    for (int i = 0; i < n_edges; ++i) {
        if (!std::getline(file, line)) {
            break;  // EOF reached
        }
        
        if (line.empty()) {
            continue;
        }
        
        // Fast parsing - find spaces
        size_t first_space = line.find(' ');
        if (first_space == std::string::npos) continue;
        
        size_t second_space = line.find(' ', first_space + 1);
        if (second_space == std::string::npos) continue;
        
        std::string u = line.substr(0, first_space);
        std::string v = line.substr(first_space + 1, second_space - first_space - 1);
        double weight = std::stod(line.substr(second_space + 1));
        
        edges.push_back({std::move(u), std::move(v), weight});
    }
    
    file.close();
    
    std::cout << "Loaded: " << n_vertices << " vertices, " << edges.size() << " edges" << std::endl;
    return true;
}

} // namespace phylo

