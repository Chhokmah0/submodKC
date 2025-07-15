#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <execution>
#include <limits>
#include <numeric>

#include "utils/graph.hpp"

template <typename Val>
struct SimpleCache {
    Val cur_val;
};

template <typename Cache, typename Val>
concept CacheConcept = requires(Cache cache, Val val) {
    { cache.cur_val } -> std::convertible_to<Val>;
};

using Index = int;
using Weight = float;

template <typename Val, typename Cache = SimpleCache<Val>>
    requires CacheConcept<Cache, Val>
struct SubmodFunction {
    using EvalVal = Val;
    using EvalCache = Cache;

    Index max_index;
    EvalVal empty_set_val;
    EvalVal max_val_ub;

    SubmodFunction(Index max_index, EvalVal empty_set_val, EvalVal max_val_ub = std::numeric_limits<EvalVal>::max())
        : max_index(max_index), empty_set_val(empty_set_val), max_val_ub(max_val_ub) {}

    virtual EvalVal operator()(const std::vector<Index>& set) const = 0;

    [[nodiscard]] virtual EvalCache init_cache() const { return EvalCache{.cur_val = empty_set_val}; }

    virtual void init_cache(EvalCache& cache) const { cache.cur_val = empty_set_val; }

    virtual void update_cache(EvalCache& cache, Index i, EvalVal add_val) const { cache.cur_val += add_val; }

    [[nodiscard]] virtual EvalVal eval_margin(const EvalCache& cache, const std::vector<Index>& set, Index i) const = 0;
};

struct DOMFunction : SubmodFunction<int> {
    using Val = SubmodFunction<int>::EvalVal;
    using Cache = SubmodFunction<int>::EvalCache;

    Graph graph;

    DOMFunction(Graph&& _graph) : SubmodFunction(_graph.size(), 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    Val operator()(const std::vector<Index>& set) const override {
        std::vector<bool> is_covered(max_index, false);
        int cnt = 0;
        for (Index i : set) {
            if (!is_covered[i]) {
                cnt++;
                is_covered[i] = true;
            }
            for (Index j : graph[i]) {
                if (!is_covered[j]) {
                    cnt++;
                    is_covered[j] = true;
                }
            }
        }
        return cnt;
    }

    [[nodiscard]] Val eval_margin(const Cache& cache, const std::vector<Index>& set, Index i) const override {
        if (set.empty()) {
            return graph[i].size() + 1;
        }
        std::vector<bool> is_covered(max_index, false);
        int cnt = 0;
        for (Index j : set) {
            if (!is_covered[j]) {
                cnt++;
                is_covered[j] = true;
            }
            for (Index k : graph[j]) {
                if (!is_covered[k]) {
                    cnt++;
                    is_covered[k] = true;
                }
            }
        }
        if (!is_covered[i]) {
            cnt++;
            is_covered[i] = true;
        }
        for (Index j : graph[i]) {
            if (!is_covered[j]) {
                cnt++;
                is_covered[j] = true;
            }
        }
        return cnt - cache.cur_val;
    }
};

struct LOCFunction : SubmodFunction<double> {
    using Val = SubmodFunction<double>::EvalVal;
    using Cache = SubmodFunction<double>::EvalCache;

    LOCGraph<double> graph;

    LOCFunction(LOCGraph<double>&& _graph)
        : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    Val operator()(const std::vector<Index>& set) const override {
        std::vector<double> max_q(graph.m, 0);
        for (Index i : set) {
            for (size_t j = 0; j < graph.m; j++) {
                max_q[j] = std::max(max_q[j], graph.gain[i][j]);
            }
        }
        return std::reduce(max_q.cbegin(), max_q.cend());
    }

    [[nodiscard]] Val eval_margin(const Cache& cache, const std::vector<Index>& set, Index i) const override {
        if (set.empty()) {
            return std::reduce(graph.gain[i].cbegin(), graph.gain[i].cend());
        }
        std::vector<double> max_q = graph.gain[i];
        for (Index j : set) {
            for (size_t k = 0; k < graph.m; k++) {
                max_q[k] = std::max(max_q[k], graph.gain[j][k]);
            }
        }
        return std::reduce(max_q.cbegin(), max_q.cend()) - cache.cur_val;
    }
};

struct INFFunction : SubmodFunction<double> {
    using Val = SubmodFunction<double>::EvalVal;
    using Cache = SubmodFunction<double>::EvalCache;

    INFGraph graph;

    INFFunction(INFGraph&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    Val operator()(const std::vector<Index>& set) const override {
        std::vector<double> fail_p(graph.m, 1);
        for (Index i : set) {
            for (size_t j = 0; j < graph.m; j++) {
                fail_p[j] *= 1 - graph.prob[i][j];
            }
        }
        for (size_t j = 0; j < graph.m; j++) {
            fail_p[j] = 1 - fail_p[j];
        }
        return std::reduce(fail_p.cbegin(), fail_p.cend());
    }

    [[nodiscard]] Val eval_margin(const Cache& cache, const std::vector<Index>& set, Index i) const override {
        if (set.empty()) {
            return std::reduce(graph.prob[i].cbegin(), graph.prob[i].cend());
        }
        std::vector<double> fail_p(graph.m, 1);
        for (Index j : set) {
            for (size_t k = 0; k < graph.m; k++) {
                fail_p[k] *= 1 - graph.prob[j][k];
            }
        }
        for (size_t k = 0; k < graph.m; k++) {
            fail_p[k] *= 1 - graph.prob[i][k];
            fail_p[k] = 1 - fail_p[k];
        }
        return std::reduce(fail_p.cbegin(), fail_p.cend()) - cache.cur_val;
    }
};

struct INFOneProbFunction : SubmodFunction<double> {
    using Val = SubmodFunction<double>::EvalVal;
    using Cache = SubmodFunction<double>::EvalCache;

    INFGraphOneProb graph;

    INFOneProbFunction(INFGraphOneProb&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    Val operator()(const std::vector<Index>& set) const override {
        std::vector<double> fail_p(graph.m, 1);
        for (Index i : set) {
            for(auto j : graph.adj[i]) {
                fail_p[j] *= 1 - graph.prob[i];
            }
        }
        for (size_t j = 0; j < graph.m; j++) {
            fail_p[j] = 1 - fail_p[j];
        }
        return std::reduce(fail_p.cbegin(), fail_p.cend());
    }

    [[nodiscard]] Val eval_margin(const Cache& cache, const std::vector<Index>& set, Index i) const override {
        std::vector<double> fail_p(graph.m, 1);
        for (Index j : set) {
            for (auto k : graph.adj[j]) {
                fail_p[k] *= 1 - graph.prob[j];
            }
        }
        for (auto k : graph.adj[i]) {
            fail_p[k] *= 1 - graph.prob[i];
        }
        for (size_t k = 0; k < graph.m; k++) {
            fail_p[k] = 1 - fail_p[k];
        }
        return std::reduce(fail_p.cbegin(), fail_p.cend()) - cache.cur_val;
    }
};

struct COVFunction : SubmodFunction<double> {
    using Val = SubmodFunction<double>::EvalVal;
    using Cache = SubmodFunction<double>::EvalCache;

    COVGraph graph;

    COVFunction(COVGraph&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    double operator()(const std::vector<Index>& set) const override {
        std::vector<char> covered(graph.m, false);
        double sum = 0;
        for (Index i : set) {
            for (Index j : graph.cover[i]) {
                if (!covered[j]) {
                    sum += graph.weights[j];
                    covered[j] = true;
                }
            }
        }
        return sum;
    }

    [[nodiscard]] double eval_margin(const Cache& cache, const std::vector<Index>& set, Index i) const override {
        std::vector<char> covered(graph.m, false);
        double sum = 0;
        for (Index j : set) {
            for (Index k : graph.cover[j]) {
                if (!covered[k]) {
                    sum += graph.weights[k];
                    covered[k] = true;
                }
            }
        }
        for (Index j : graph.cover[i]) {
            if (!covered[j]) {
                sum += graph.weights[j];
                covered[j] = true;
            }
        }
        return sum - cache.cur_val;
    }
};

// -----------------------------------------------------------

struct DomCache {
    int cur_val;
    std::vector<char> is_covered;
};

struct DomFunctionWithCache : SubmodFunction<int, DomCache> {
    Graph graph;

    DomFunctionWithCache(Graph&& _graph) : SubmodFunction(_graph.size(), 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    int operator()(const std::vector<Index>& set) const override {
        std::vector<bool> is_covered(max_index, false);
        int cnt = 0;
        for (Index i : set) {
            if (!is_covered[i]) {
                cnt++;
                is_covered[i] = true;
            }
            for (Index j : graph[i]) {
                if (!is_covered[j]) {
                    cnt++;
                    is_covered[j] = true;
                }
            }
        }
        return cnt;
    }

    [[nodiscard]] DomCache init_cache() const override {
        return DomCache{.cur_val = 0, .is_covered = std::vector<char>(max_index, false)};
    }

    void init_cache(DomCache& cache) const override {
        cache.cur_val = 0;
        cache.is_covered = std::vector<char>(max_index, false);
    }

    void update_cache(DomCache& cache, Index i, int add_val) const override {
        cache.cur_val += add_val;
        cache.is_covered[i] = true;
        for (Index j : graph[i]) {
            cache.is_covered[j] = true;
        }
    }

    [[nodiscard]] int eval_margin(const DomCache& cache, const std::vector<Index>&, Index i) const override {
        int cnt = 0;
        if (!cache.is_covered[i]) {
            cnt++;
        }
        for (Index j : graph[i]) {
            if (!cache.is_covered[j]) {
                cnt++;
            }
        }
        return cnt;
    }
};

struct INFCache {
    double cur_val;
    std::vector<double> fail_p;
};

struct INFFunctionWithCache : SubmodFunction<double, INFCache> {
    INFGraph graph;

    INFFunctionWithCache(INFGraph&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    double operator()(const std::vector<Index>& set) const override {
        std::vector<double> fail_p(graph.m, 1);
        for (Index i : set) {
            for (size_t j = 0; j < graph.m; j++) {
                fail_p[j] *= 1 - graph.prob[i][j];
            }
        }
        for (size_t j = 0; j < graph.m; j++) {
            fail_p[j] = 1 - fail_p[j];
        }
        return std::reduce(fail_p.cbegin(), fail_p.cend());
    }
    [[nodiscard]] INFCache init_cache() const override {
        return INFCache{.cur_val = 0, .fail_p = std::vector<double>(graph.m, 1)};
    }
    void init_cache(INFCache& cache) const override {
        cache.cur_val = 0;
        cache.fail_p = std::vector<double>(graph.m, 1);
    }
    void update_cache(INFCache& cache, Index i, double add_val) const override {
        cache.cur_val += add_val;
        for (size_t j = 0; j < graph.m; j++) {
            cache.fail_p[j] *= 1 - graph.prob[i][j];
        }
    }
    [[nodiscard]] double eval_margin(const INFCache& cache, const std::vector<Index>& set, Index i) const override {
        if (set.empty()) {
            return std::reduce(graph.prob[i].cbegin(), graph.prob[i].cend());
        }
        return std::transform_reduce(graph.prob[i].cbegin(), graph.prob[i].cend(), cache.fail_p.cbegin(), 0.0);
    }
};

struct LOCCache {
    double cur_val;
    std::vector<double> max_gain;
};

struct LOCFunctionWithCache : SubmodFunction<double, LOCCache> {
    LOCGraph<double> graph;

    LOCFunctionWithCache(LOCGraph<double>&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    double operator()(const std::vector<Index>& set) const override {
        std::vector<double> max_gain(graph.m, 0);
        for (Index i : set) {
            for (size_t j = 0; j < graph.m; j++) {
                max_gain[j] = std::max(max_gain[j], graph.gain[i][j]);
            }
        }
        return std::reduce(max_gain.cbegin(), max_gain.cend());
    }
    [[nodiscard]] LOCCache init_cache() const override {
        return LOCCache{.cur_val = 0, .max_gain = std::vector<double>(graph.m, 0)};
    }
    void init_cache(LOCCache& cache) const override {
        cache.cur_val = 0;
        cache.max_gain = std::vector<double>(graph.m, 0);
    }
    void update_cache(LOCCache& cache, Index i, double add_val) const override {
        cache.cur_val += add_val;
        for (size_t j = 0; j < graph.m; j++) {
            cache.max_gain[j] = std::max(cache.max_gain[j], graph.gain[i][j]);
        }
    }
    [[nodiscard]] double eval_margin(const LOCCache& cache, const std::vector<Index>& set, Index i) const override {
        if (set.empty()) {
            return std::reduce(graph.gain[i].cbegin(), graph.gain[i].cend());
        }
        double ans = 0;
        for (size_t j = 0; j < graph.m; j++) {
            ans += std::max(graph.gain[i][j], cache.max_gain[j]);
        }
        return ans - cache.cur_val;
    }
};

struct COVCache {
    double cur_val;
    std::vector<char> covered;
};

struct COVFunctionWithCache : SubmodFunction<double, COVCache> {
    COVGraph graph;

    COVFunctionWithCache(COVGraph&& _graph) : SubmodFunction(_graph.n, 0), graph(std::move(_graph)) {
        // Only right when submod_func is monotonic
        std::vector<Index> indices(max_index);
        std::iota(indices.begin(), indices.end(), 0);
        max_val_ub = (*this)(indices);
    }

    double operator()(const std::vector<Index>& set) const override {
        std::vector<char> covered(graph.m, false);
        double sum = 0;
        for (Index i : set) {
            for (Index j : graph.cover[i]) {
                if (!covered[j]) {
                    sum += graph.weights[j];
                    covered[j] = true;
                }
            }
        }
        return sum;
    }

    [[nodiscard]] COVCache init_cache() const override {
        return COVCache{.cur_val = 0, .covered = std::vector<char>(graph.m, false)};
    }

    void init_cache(COVCache& cache) const override {
        cache.cur_val = 0;
        cache.covered = std::vector<char>(graph.m, false);
    }

    void update_cache(COVCache& cache, Index i, double add_val) const override {
        cache.cur_val += add_val;
        for (Index j : graph.cover[i]) {
            cache.covered[j] = true;
        }
    }

    [[nodiscard]] double eval_margin(const COVCache& cache, const std::vector<Index>& set, Index i) const override {
        double cnt = 0;
        for (Index j : graph.cover[i]) {
            if (!cache.covered[j]) {
                cnt += graph.weights[j];
            }
        }
        return cnt;
    }
};
