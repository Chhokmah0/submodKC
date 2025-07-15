#pragma once

#include <vector>
#include <algorithm>
#include <fstream>

using Index = int;
using Graph = std::vector<std::vector<Index>>;

template <typename T>
struct LOCGraph {
    int n;
    int m;
    std::vector<std::vector<T>> gain;
};

struct INFGraph {
    int n,m;
    std::vector<std::vector<double>> prob;
};

struct INFGraphOneProb {
    int n,m;
    std::vector<double> prob;
    std::vector<std::vector<Index>> adj;
};

struct COVGraph {
    int n,m;
    std::vector<double> weights;
    std::vector<std::vector<Index>> cover; 
};

namespace utils {

inline Graph read_graph(const std::string& filename) {
    std::ifstream file(filename);
    int n, m;
    std::string _;
    file >> _ >> _ >> n >> m;
    Graph graph(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        file >> u >> v;
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    return graph;
}

inline Graph read_graph_edges(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<int, int>> edges;
    int u, v;
    int n = 0, m = 0;
    while (file >> u >> v) {
        edges.emplace_back(u, v);
        n = std::max({n, u + 1, v + 1});
        m++;
    }
    Graph graph(n);
    for (auto [u, v] : edges) {
        graph[u].push_back(v);
        graph[v].push_back(u);
    }
    return graph;
}

}  // namespace utils
