#pragma once

#include <algorithm>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <numbers>

#include "structures/submod_function.hpp"
#include "utils/graph.hpp"

namespace utils {

inline std::vector<Weight> random_uniform_weight(int n, int seed) {
    std::vector<Weight> weights(n);
    std::mt19937 rnd(seed);
    std::uniform_real_distribution<Weight> dis(0.4, 1.6);
    for (Weight& weight : weights) {
        weight = dis(rnd);
    }
    return weights;
}

inline std::vector<Weight> random_normal_weight(int n, int seed) {
    std::vector<Weight> weights(n);
    std::mt19937 rnd(seed);
    std::normal_distribution<Weight> normal_dis(1, 0.2);
    for (Weight& weight : weights) {
        weight = std::clamp(normal_dis(rnd), 0.1f, 1.9f);
    }
    return weights;
}

inline std::vector<Weight> dom_weight(int n, int seed) {
    std::vector<Weight> weights(n);
    std::mt19937 rnd(seed);
    std::uniform_real_distribution<Weight> dis(0.1, 1.0);
    for (Weight& weight : weights) {
        weight = dis(rnd);
    }
    return weights;
}

inline std::vector<Weight> unit_weight(int n) {
    return std::vector<Weight>(n, 1);
}

inline std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(str);
    std::string token;

    while (std::getline(stream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

inline std::vector<std::string> read_lines(const std::string& file_path) {
    std::ifstream file(file_path);
    if (file.fail()) {
        throw std::runtime_error("Fail when open file " + file_path);
    }

    std::string s;
    std::vector<std::string> lines;
    while (!file.eof()) {
        std::getline(file, s);
        if (s.empty()) {
            continue;
        }
        // trim
        s.erase(s.find_last_not_of(" \n\r\t") + 1);
        lines.emplace_back(std::move(s));
        s.clear();
    }
    return lines;
}

inline auto DOM_read_from(const std::string& file_path) {
    DomFunctionWithCache dom_func{read_graph_edges(file_path)};
    return dom_func;
}

inline auto INF_read_from(const std::string& file_path) {
    auto lines = read_lines(file_path);
    auto n = split(lines[0], ',').size();
    auto m = lines.size();
    std::vector<std::vector<double>> prob(n, std::vector<double>(m));
    for (int j = 0; j < m; j++) {
        auto tokens = split(lines[j], ',');
        for (int i = 0; i < n; i++) {
            prob[i][j] = std::stod(tokens[i]);
        }
    }
    INFGraph graph{.n = static_cast<int>(n), .m = static_cast<int>(m), .prob = std::move(prob)};
    return INFFunctionWithCache(std::move(graph));
}

inline auto LOC_read_from(const std::string& file_path) {
    auto lines = read_lines(file_path);
    // auto n = lines.size();
    // auto m = split(lines[0], ',').size();
    auto n = split(lines[0], ',').size();
    auto m = lines.size();
    std::vector<std::vector<double>> gain(n, std::vector<double>(m));
    for (int j = 0; j < m; j++) {
        auto tokens = split(lines[j], ',');
        for (int i = 0; i < n; i++) {
            gain[i][j] = std::stod(tokens[i]);
        }
    }
    LOCGraph graph{.n = static_cast<int>(n), .m = static_cast<int>(m), .gain = std::move(gain)};
    return LOCFunctionWithCache(std::move(graph));
}

inline auto COV_read_from(const std::string& file_path) {
    auto lines = read_lines(file_path);
    auto n = split(lines[1], ',').size();
    auto m = split(lines[0], ',').size();
    std::vector<double> weights(m);
    auto tokens = split(lines[0], ',');
    for (int i = 0; i < m; i++) {
        weights[i] = std::stod(tokens[i]);
    }
    std::vector<std::vector<int>> cover(n);
    for (int j = 0; j < m; j++) {
        auto tokens = split(lines[j + 1], ',');
        for (int i = 0; i < n; i++) {
            if (std::stod(tokens[i]) == 1) {
                cover[i].push_back(j);
            }
        }
    }
    COVGraph graph{.n = static_cast<int>(n), .m = static_cast<int>(m), .weights = std::move(weights), .cover = std::move(cover)};
    return COVFunctionWithCache(std::move(graph));
}

inline auto COV_read_from_Facebook_like(const std::string& file_path) {
    auto lines = read_lines(file_path);

    std::vector<std::tuple<Index, Index, int>> temp;

    // parse
    Index l = 0, r = 0;
    for (auto& line_str : lines) {
        std::istringstream line(line_str);
        Index u, v;
        int w;
        line >> u >> v >> w;
        l = std::max(l, u);
        r = std::max(r, v);
        if (u == 0 || v == 0) {
            throw "0-index in file " + file_path;
        }
        u--;
        v--;
        temp.push_back({u, v, w});
    }

    // construct graph
    COVGraph g;
    g.n = l;
    g.m = r;
    g.cover.resize(g.n);
    g.weights.resize(g.m);
    for (auto [u, v, w] : temp) {
        g.cover[u].emplace_back(v);
        g.weights[v] += w;
    }
    return COVFunctionWithCache(std::move(g));
}

constexpr float EARTH_RADIUS = 6371.0f;

inline float to_radians(float degree) { return degree * std::numbers::pi / 180.0f; }

inline float haversine_distance(float lat1, float lon1, float lat2, float lon2) {
    lat1 = to_radians(lat1);
    lon1 = to_radians(lon1);
    lat2 = to_radians(lat2);
    lon2 = to_radians(lon2);

    float dLat = lat2 - lat1;
    float dLon = lon2 - lon1;

    float a = std::sin(dLat / 2) * std::sin(dLat / 2) +
              std::cos(lat1) * std::cos(lat2) * std::sin(dLon / 2) * std::sin(dLon / 2);

    float c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));

    float distance = EARTH_RADIUS * c;

    return distance;
}

auto LOC_read_from_MTA_Swubway_Stations(const std::string& file_path) {
    std::vector<std::string> lines = read_lines(file_path);
    std::vector<std::pair<float, float>> stations;
    float minx = std::numeric_limits<float>::max();
    float maxx = std::numeric_limits<float>::lowest();
    float miny = std::numeric_limits<float>::max();
    float maxy = std::numeric_limits<float>::lowest();

    for (int i = 1; i < lines.size(); i++) {
        float x, y;
        if (sscanf(split(lines[i], ',').back().c_str(), "POINT (%f %f)", &x, &y) != 2) {
            throw std::runtime_error("Error when read line in file " + file_path);
        }
        stations.emplace_back(x, y);
        minx = std::min(minx, x);
        maxx = std::max(maxx, x);
        miny = std::min(miny, y);
        maxy = std::max(maxy, y);
    }

    LOCGraph<double> g;
    g.n = 100;
    g.m = lines.size() - 1;

    g.gain.resize(g.n);
    float dx = (maxx - minx) / 10;
    float dy = (maxx - minx) / 10;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            float x = minx + dx * i;
            float y = miny + dy * j;

            Index id = i * 10 + j;
            g.gain[id].reserve(g.m);
            for (const auto& [x_, y_] : stations) {
                float distance = haversine_distance(x, y, x_, y_);
                if (distance == 0) {
                    g.gain[id].emplace_back(10);
                    continue;
                }
                g.gain[id].emplace_back(1 / haversine_distance(x, y, x_, y_));
            }
        }
    }

    return LOCFunctionWithCache(std::move(g));
}

auto INF_read_from_Movie_100k(const std::string& file_path) {
    std::vector<std::string> lines = read_lines(file_path);

    std::vector<std::tuple<Index, Index, int>> temp;
    Index l = 0;
    Index r = 0;
    for (const std::string& s : lines) {
        std::istringstream line(s);
        Index u, v;
        // u: movie
        // v: people
        int rating;
        line >> v >> u >> rating;

        l = std::max(l, u);
        r = std::max(r, v);
        temp.emplace_back(u - 1, v - 1, rating);
    }
    std::vector<int> ratings(l);
    std::vector<int> cnt(l);
    for (const auto& [u, _, p] : temp) {
        ratings[u] += p;
        cnt[u] += 1;
    }

    INFGraphOneProb g;
    g.n = l;
    g.m = r;
    g.prob.resize(g.n);
    for (int i = 0; i < g.n; i++) {
        g.prob[i] = ((double)ratings[i] / cnt[i]) / 5;
    }

    g.adj.resize(g.n);
    for (const auto& [u, v, p] : temp) {
        g.adj[u].emplace_back(v);
    }

    return INFOneProbFunction(std::move(g));
}

}  // namespace utils
