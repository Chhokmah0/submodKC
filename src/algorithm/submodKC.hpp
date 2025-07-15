#pragma once
#include <algorithm>
#include <atomic>
#include <cassert>
#include <compare>
#include <limits>
#include <nlohmann/json.hpp>
#include <numeric>
#include <queue>
#include <ranges>
#include <vector>

#include "structures/submod_function.hpp"
#include "utils/args.hpp"
#include "utils/timer.hpp"

template <typename SubmodFunc>
class SubmodKC {
   public:
    using Val = typename SubmodFunc::EvalVal;
    using Cache = typename SubmodFunc::EvalCache;

   private:
    SubmodFunc submod_func;
    Weight max_weight;
    std::vector<Weight> weights;

    size_t search_nodes = 0;
    // timer
    utils::Timer<> running_timer;
    bool early_stop = false;
    // time limit
    std::chrono::nanoseconds max_time = std::chrono::nanoseconds(0);

    double root_upper_bound = std::numeric_limits<double>::min();

   public:
    Val best_val;
    std::vector<Index> best_set;

    template <typename T>
    SubmodKC(SubmodFunc _submod_func, Weight _max_weight, std::vector<Weight> _weights, T _max_time = 0)
        : submod_func(std::move(_submod_func)),
          best_set(),
          best_val(submod_func.empty_set_val),
          max_weight(_max_weight),
          weights(std::move(_weights)),
          max_time(std::chrono::duration_cast<std::chrono::nanoseconds>(_max_time)) {
        assert(weights.size() == submod_func.max_index);
    }
    SubmodKC(SubmodFunc _submod_func, Weight _max_weight, std::vector<Weight> _weights)
        : SubmodKC(std::move(_submod_func), _max_weight, std::move(_weights), std::chrono::nanoseconds(0)) {}

    bool check() {
        bool ok = true;
        Weight weight = 0;
        for (auto i : best_set) {
            weight += weights[i];
        }
        if (weight > max_weight) {
            std::cout << "best_set weight is not correct\n";
            ok = false;
        }
        auto val = submod_func(best_set);
        if (std::abs(val - best_val) > 1e-10) {
            std::cout << "best_val is not correct\n";
            std::cout << "val: " << val << '\n';
            std::cout << "best_val: " << best_val << '\n';
            ok = false;
        }
        return ok;
    }

   private:
    struct Candidate {
        Val margin;
        Weight weight;
        Index id;
        bool acc;  // lazy update

        std::partial_ordering operator<=>(const Candidate& other) const {
            return margin * other.weight <=> other.margin * weight;
        }
    };

    static bool stable_cand_order(const Candidate& a, const Candidate& b) {
        if (a.margin * b.weight == b.margin * a.weight) {
            return a.acc < b.acc;
        }
        return a.margin * b.weight < b.margin * a.weight;
    }

    struct Candidates {
        std::vector<Candidate> data;
        std::vector<Val> suffix_margin_sum;
        std::vector<Weight> suffix_weight_sum;

        [[nodiscard]] bool empty() const { return data.empty(); }
        void clear() {
            data.clear();
            suffix_margin_sum.clear();
            suffix_weight_sum.clear();
        }
        [[nodiscard]] Weight tot_weight() const {
            return suffix_weight_sum.front() - suffix_weight_sum.back() + data.front().weight;
        }
        [[nodiscard]] Val tot_margin() const {
            return suffix_margin_sum.front() - suffix_margin_sum.back() + data.front().margin;
        }

        void maintain_suffix_sum() {
            suffix_margin_sum.resize(data.size());
            suffix_weight_sum.resize(data.size());
            if (data.empty()) {
                return;
            }
            suffix_margin_sum.back() = 0;
            suffix_weight_sum.back() = 0;
            for (size_t i = data.size() - 1; i > 0; --i) {
                suffix_margin_sum[i - 1] = suffix_margin_sum[i] + data[i].margin;
                suffix_weight_sum[i - 1] = suffix_weight_sum[i] + data[i].weight;
            }
        }
        void sort() {
            std::ranges::sort(data, stable_cand_order);
            maintain_suffix_sum();
        }
        void clear_acc() {
            for (auto& cand : data) {
                cand.acc = false;
            }
        }
        [[nodiscard]] double fractional_knapsack(Weight max_weight) const {
            if (data.empty() || max_weight == 0) {
                return 0;
            }
            max_weight += suffix_weight_sum.back();
            if (tot_weight() <= max_weight) {
                return tot_margin();
            }
            auto p =
                std::ranges::upper_bound(suffix_weight_sum, max_weight, std::greater<>()) - suffix_weight_sum.begin();
            return suffix_margin_sum[p] - suffix_margin_sum.back() +
                   (max_weight - suffix_weight_sum[p]) * data[p].margin / data[p].weight;
        }
        Candidate pop_back() {
            auto cand = std::move(data.back());
            data.pop_back();
            suffix_margin_sum.pop_back();
            suffix_weight_sum.pop_back();
            return cand;
        }
        Candidate erase_addable(Weight res_weight) {
            for (int i = (int)data.size() - 1; i >= 0; --i) {
                if (data[i].weight <= res_weight) {
                    auto cand = std::move(data[i]);
                    data.erase(data.begin() + i);
                    maintain_suffix_sum();
                    return cand;
                }
            }
            return Candidate{.margin = 0, .weight = 0, .id = -1, .acc = false};  // Return an invalid candidate
        }

        Weight reduction_weight(Val res_val) const {
            if (res_val <= 0) {
                return 0;
            }
            res_val += suffix_weight_sum.back();
            auto p = std::ranges::upper_bound(suffix_margin_sum, res_val, std::greater<>()) - suffix_margin_sum.begin();
            if (p == suffix_margin_sum.size()) {
                return 0;
            }
            return suffix_weight_sum[p] - suffix_weight_sum.back() +
                   (res_val - suffix_margin_sum[p]) * data[p].weight / data[p].margin;
        }
    };

    struct Node {
        std::vector<Index> set;
        Candidates candidates;
        Cache cache;
        Weight cur_weight;
    };

   public:
    Val solve() {
        std::cout << "submod_func max_index: " << submod_func.max_index << '\n';
        running_timer.start();
        Node root{
            .set = {},
            .candidates = Candidates{},
            .cache = submod_func.init_cache(),
            .cur_weight = 0,
        };
        root.candidates.data.reserve(submod_func.max_index);
        for (Index i = 0; i < submod_func.max_index; ++i) {
            root.candidates.data.push_back(Candidate{.margin = submod_func.eval_margin(root.cache, root.set, i),
                                                     .weight = weights[i],
                                                     .id = i,
                                                     .acc = true});
        }
        root.candidates.sort();
        // root_upper_bound = upper_bound_with_arg(root);
        // return 0;
        switch (args::algorithm) {
            case args::Algorithm::BASIC:
                basic_search(root);
                break;
            case args::Algorithm::DUAL:
                dual_search(root);
                break;
            case args::Algorithm::BFS:
                raw_bfs_dom_search(std::move(root));
                break;
            case args::Algorithm::GREEDY:
                seq_greedy(root);
                break;
            default:
                throw std::runtime_error("Invalid algorithm");
        }
        std::cout << "searched nodes: " << search_nodes << '\n';
        std::cout << "running time: " << (double)running_timer.stop().count() / 1e9 << "s\n";
        return best_val;
    }

    // recursive search
   private:
    void basic_search(Node& node) {
        search_nodes++;
        maintain_best_solution(node);
        if (node.candidates.data.empty()) {
            return;
        }
        if (node.cache.cur_val + node.candidates.tot_margin() <= best_val) {
            return;
        }
        if (upper_bound_fractional_knapsack(node) <= best_val) {
            return;
        }
        if (node.cache.cur_val >= 51.87) {
            int stop = 0;
        }
        reduction(node);
        if (node.candidates.data.empty()) {
            return;
        }
        if (node.cache.cur_val + node.candidates.tot_margin() <= best_val) {
            return;
        }
        if (upper_bound_fractional_knapsack(node) <= best_val) {
            return;
        }
        double ub = upper_bound_with_arg(node);
        if (root_upper_bound == std::numeric_limits<double>::min()) {
            root_upper_bound = ub;
        }
        if (ub <= best_val) {
            return;
        }

        while (!node.candidates.empty()) {
            if (best_val >= submod_func.max_val_ub) {
                return;
            }
            if (max_time > std::chrono::nanoseconds(0) && running_timer.cur_time() > max_time) {
                early_stop = true;
                return;
            }
            if (node.cache.cur_val + node.candidates.tot_margin() <= best_val) {
                return;
            }
            if (upper_bound_fractional_knapsack(node) <= best_val) {
                return;
            }
            if (upper_bound_with_arg(node) <= best_val) {
                return;
            }
            auto cand = node.candidates.pop_back();
            if (node.cur_weight + cand.weight > max_weight) {
                continue;
            }

            Node new_node = node;
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            choose(new_node, cand);
            update_margin_with_arg(new_node);
            basic_search(new_node);

            // reduction(node);
        }
    }

    void dual_search(Node& node) {
        search_nodes++;
        maintain_best_solution(node);
        if (node.candidates.data.empty()) {
            return;
        }
        if (node.cache.cur_val + node.candidates.tot_margin() <= best_val) {
            return;
        }
        if (upper_bound_fractional_knapsack(node) <= best_val) {
            return;
        }
        reduction(node);
        if (node.candidates.data.empty()) {
            return;
        }
        if (node.cache.cur_val + node.candidates.tot_margin() <= best_val) {
            return;
        }
        if (upper_bound_fractional_knapsack(node) <= best_val) {
            return;
        }

        bool cut = false;
        Weight res_weight = max_weight - node.cur_weight;

        std::vector<Node> children;
        std::vector<double> ub_rs;

        while (!node.candidates.empty()) {
            if (best_val >= submod_func.max_val_ub) {
                return;
            }
            if (max_time > std::chrono::nanoseconds(0) && running_timer.cur_time() > max_time) {
                early_stop = true;
                return;
            }
            // Here, we cannot immediately prune through ub_fk because we need to calculate ub_rs.
            if (!cut && upper_bound_fractional_knapsack(node) <= best_val) {
                cut = true;
            }
            auto cand = node.candidates.erase_addable(max_weight - node.cur_weight);
            if (cand.id == -1) {
                break;
            }
            if (!cut) {
                children.emplace_back(node);
                ub_rs.emplace_back(std::numeric_limits<double>::max());
            }

            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            choose(node, cand);
            // update_margin_threshold(node);
            if (args::set_bound_update_type == args::SetBoundUpdateType::GLOBAL) {
                update_margin_with_arg(node, res_weight);
            } else {
                update_margin_with_arg(node);
            }
            // update_margin(node);
            maintain_best_solution(node);

            for (int i = 0; i < children.size(); ++i) {
                ub_rs[i] =
                    std::min(ub_rs[i], upper_bound_fractional_knapsack(node, max_weight - children[i].cur_weight));
            }
            if (ub_rs[0] <= best_val) {
                break;
            }
        }
        if (root_upper_bound == std::numeric_limits<double>::min()) {
            root_upper_bound = upper_bound_fractional_knapsack(node);
        }
        if (ub_rs[0] <= best_val) {
            return;
        }

        for (int i = 0; i < children.size(); ++i) {
            int c = (int)children.size() - 1 - i;
            if (ub_rs[c] <= best_val) {
                continue;
            }
            dual_search(children[c]);
        }
    }

    // best first search
    void raw_bfs_dom_search(Node node) {
        struct BFSNode {
            double ub;
            Node node;

            bool operator<(const BFSNode& other) const { return ub < other.ub; }
        };
        auto cmp_id = [](const Candidate& a, const Candidate& b) { return a.id > b.id; };
        std::priority_queue<BFSNode> pq;
        double g_ub = upper_bound_dom(node);
        root_upper_bound = g_ub;
        pq.push({.ub = g_ub, .node = std::move(node)});
        // constexpr size_t max_bfs_node = (size_t)32 * 1024 * 1024;
        while (!pq.empty() && pq.top().ub > best_val) {
            // std::cout << "pq size: " << pq.size() << '\n';
            // std::cout << "upper bound: " << pq.top().ub << '\n';
            if (best_val >= submod_func.max_val_ub) {
                return;
            }
            if (max_time > std::chrono::nanoseconds(0) && running_timer.cur_time() > max_time) {
                early_stop = true;
                return;
            }
            search_nodes++;
            auto [ub, node] = std::move(pq.top());
            pq.pop();
            if (std::abs(ub - node.cache.cur_val) <= 1e-6) {
                return;
            }
            g_ub = std::min(g_ub, ub);
            std::ranges::sort(node.candidates.data, cmp_id);
            while (!node.candidates.empty()) {
                if (max_time > std::chrono::nanoseconds(0) && running_timer.cur_time() > max_time) {
                    early_stop = true;
                    return;
                }
                auto cand = node.candidates.pop_back();
                if (node.cur_weight + cand.weight > max_weight) {
                    continue;
                }
                auto new_node = node;
                choose(new_node, cand);
                maintain_best_solution(new_node);
                std::erase_if(new_node.candidates.data, [&](const Candidate& cand) {
                    // over weight reduction
                    if (cand.weight > max_weight - new_node.cur_weight) {
                        return true;
                    }
                    // non-positive margin reduction
                    if (cand.margin <= 0) {
                        return true;
                    }
                    return false;
                });
                update_margin(new_node);
                
                double ub = upper_bound_dom(new_node);
                if (best_val >= g_ub || best_val >= submod_func.max_val_ub) {
                    return;
                }
                if (ub > best_val) {
                    pq.push({.ub = ub, .node = std::move(new_node)});
                }
            }
        }
    }

    // lazy update
   private:
    void update_margin(Node& node) {
        for (auto& cand : node.candidates.data) {
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
        }
        node.candidates.sort();
    }

    void update_margin_threshold(Node& node, const Candidate& threshold) {
        bool updated = false;
        for (auto& cand : node.candidates.data | std::views::reverse) {
            if (stable_cand_order(cand, threshold)) {
                break;
            }
            if (!cand.acc) {
                updated = true;
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
        }
        if (updated) {
            node.candidates.sort();
        }
    }
    void update_margin_threshold(Node& node) {
        Candidate threshold{
            .margin = best_val - node.cache.cur_val, .weight = max_weight - node.cur_weight, .id = -1, .acc = true};
        update_margin_threshold(node, threshold);
    }
    void update_margin_threshold(Node& node, const Weight weight_threshold) {
        Candidate threshold{.margin = best_val - node.cache.cur_val, .weight = weight_threshold, .id = -1, .acc = true};
        update_margin_threshold(node, threshold);
    }

    void update_margin_max(Node& node, Candidate threshold) {
        bool updated = false;
        for (auto& cand : node.candidates.data | std::views::reverse) {
            if (stable_cand_order(cand, threshold)) {
                break;
            }
            if (!cand.acc) {
                updated = true;
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            threshold = std::max(threshold, cand);
        }
        if (updated) {
            node.candidates.sort();
        }
    }
    void update_margin_max(Node& node) {
        Candidate threshold{
            .margin = best_val - node.cache.cur_val, .weight = max_weight - node.cur_weight, .id = -1, .acc = true};
        update_margin_max(node, threshold);
    }
    void update_margin_max(Node& node, const Weight weight_threshold) {
        Candidate threshold{.margin = best_val - node.cache.cur_val, .weight = weight_threshold, .id = -1, .acc = true};
        update_margin_max(node, threshold);
    }

    void update_margin_with_arg(Node& node, const Weight weight_threshold) {
        switch (args::update_type) {
            case args::UpdateType::ALL:
                update_margin(node);
                break;
            case args::UpdateType::THRESHOLD:
                update_margin_threshold(node, weight_threshold);
                break;
            case args::UpdateType::MAX:
                update_margin_max(node, weight_threshold);
                break;
            default:
                throw std::runtime_error("Invalid update type");
        }
    }
    void update_margin_with_arg(Node& node) {
        Weight weight_threshold = max_weight - node.cur_weight;
        update_margin_with_arg(node, weight_threshold);
    }

    // reduction
   private:
    void reduction(Node& node) { reduction(node, max_weight - node.cur_weight); }

    void reduction(Node& node, Weight weight_threshold) {
        if (upper_bound_fractional_knapsack(node, weight_threshold) <= best_val) {
            node.candidates.clear();
            return;
        }
        int pre_size = (int)node.candidates.data.size();

        auto check_reduction = [&](const Candidate& cand) {
            // over weight reduction
            if (cand.weight > weight_threshold) {
                return true;
            }
            // non-positive margin reduction
            if (cand.margin <= 0) {
                return true;
            }
            // upper bound reduction (Candinate Reduction)
            double val = node.cache.cur_val + cand.margin;
            return val + node.candidates.fractional_knapsack(weight_threshold - cand.weight) <= best_val;
        };
        std::vector<char> is_reduced(node.candidates.data.size(), false);
        for (int i = 0; i < node.candidates.data.size(); i++) {
            is_reduced[i] = check_reduction(node.candidates.data[i]);
        }

        int p = 0;
        for (int i = 0; i < node.candidates.data.size(); i++) {
            if (!is_reduced[i]) {
                node.candidates.data[p] = std::move(node.candidates.data[i]);
                p++;
            }
        }
        node.candidates.data.resize(p);

        if (node.candidates.data.size() < pre_size) {
            node.candidates.maintain_suffix_sum();
        }
    }

    // maintain node
   private:
    void choose(Node& node, const Candidate& cand) {
        assert(cand.acc);
        node.set.push_back(cand.id);
        submod_func.update_cache(node.cache, cand.id, cand.margin);
        node.cur_weight += cand.weight;
        node.candidates.clear_acc();
    }

    // upper bounds
   private:
    double upper_bound_fractional_knapsack(const Node& node) {
        return node.cache.cur_val + node.candidates.fractional_knapsack(max_weight - node.cur_weight);
    }

    double upper_bound_fractional_knapsack(const Node& node, const Weight weight_threshold) {
        return node.cache.cur_val + node.candidates.fractional_knapsack(weight_threshold);
    }

    double upper_bound_knapsack(const Node& node, const double epsilon = 1) {
        const size_t n = node.candidates.data.size();
        const Weight res_weight = max_weight - node.cur_weight;

        Val max_margin = 0;
        for (const auto& cand : node.candidates.data) {
            max_margin = std::max(max_margin, cand.margin);
        }

        if (max_margin <= 0) {
            return node.cache.cur_val;
        }

        const double K = epsilon * max_margin / n;
        constexpr Weight MAX_WEIGHT = std::numeric_limits<Weight>::max() / 3;

        int max_val = 0;
        double sum_weight = 0;
        // Bound of the knapsack
        for (const auto& cand : node.candidates.data | std::views::reverse) {
            max_val += ceil(cand.margin / K);
            sum_weight += cand.weight;
            if (sum_weight > res_weight) {
                break;
            }
        }

        std::vector<Weight> dp(max_val + 1, MAX_WEIGHT);

        dp[0] = 0;

        for (const auto& cand : node.candidates.data | std::views::reverse) {
            const auto value = ceil(cand.margin / K);
            for (int val = dp.size() - 1; val >= value; val--) {
                dp[val] = std::min(dp[val], dp[val - value] + cand.weight);
            }
        }

        int dp_val = max_val;
        for (int val = 0; val < dp.size(); val++) {
            if (dp[val] != MAX_WEIGHT && node.cur_weight + dp[val] <= max_weight) {
                dp_val = val;
            }
        }
        return node.cache.cur_val + dp_val * K;
    }

    double upper_bound_refined_subset(Node node) {
        Weight res_weight = max_weight - node.cur_weight;
        double ub_rs = upper_bound_fractional_knapsack(node);
        while (!node.candidates.empty()) {
            // auto cand = node.candidates.pop_back();
            // if (node.cur_weight + cand.weight > res_weight) {
            //     break;
            // }
            auto cand = node.candidates.erase_addable(max_weight - node.cur_weight);
            if (cand.id == -1) {
                break;
            }
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            choose(node, cand);
            maintain_best_solution(node);
            if (args::set_bound_update_type == args::SetBoundUpdateType::GLOBAL) {
                update_margin_with_arg(node, res_weight);
            } else {
                update_margin_with_arg(node);
            }
            ub_rs = std::min(ub_rs, upper_bound_fractional_knapsack(node, res_weight));
            if (ub_rs <= best_val) {
                break;
            }
        }
        return ub_rs;
    }

    double upper_bound_dom(Node node) {
        Weight res_weight = max_weight - node.cur_weight;
        double ub_dom = upper_bound_fractional_knapsack(node);
        double beta = 1;
        Val pre_val = node.cache.cur_val;

        while (!node.candidates.empty()) {
            double s_ub = upper_bound_fractional_knapsack(node, res_weight) - node.cache.cur_val;
            // auto cand = node.candidates.pop_back();
            // if (node.cur_weight + cand.weight > res_weight) {
            //     break;
            // }
            auto cand = node.candidates.erase_addable(max_weight - node.cur_weight);
            if (cand.id == -1) {
                break;
            }
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            if (s_ub != 0) {
                double beta_i = (1 - cand.margin / s_ub);
                beta *= beta_i;
            } else {
                return node.cache.cur_val;
            }
            choose(node, cand);
            maintain_best_solution(node);
            if (args::set_bound_update_type == args::SetBoundUpdateType::GLOBAL) {
                update_margin_with_arg(node, res_weight);
            } else {
                update_margin_with_arg(node);
            }
        }
        if (beta != 1) {
            ub_dom = std::min(ub_dom, pre_val + (node.cache.cur_val - pre_val) / (1 - beta));
        }
        return ub_dom;
    }

    double upper_bound_with_arg(const Node& node) {
        switch (args::upper_bound_type) {
            case args::UpperBoundType::KNAPSACK:
                return upper_bound_knapsack(node);
            case args::UpperBoundType::FRAC_KNAPSACK:
                return upper_bound_fractional_knapsack(node);
            case args::UpperBoundType::REFINE_SUBSET:
                return upper_bound_refined_subset(node);
            case args::UpperBoundType::DOM:
                return upper_bound_dom(node);
            default:
                throw std::runtime_error("Invalid upper bound type");
        }
    }

    void seq_greedy(Node& node) {
        while (!node.candidates.data.empty()) {
            auto cand = node.candidates.pop_back();
            if (node.cur_weight + cand.weight > max_weight) {
                continue;
            }
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
            choose(node, cand);
            maintain_best_solution(node);
            update_margin(node);
            reduction(node);
        }
    }

    // maintain solution
   private:
    void maintain_best_solution(const Node& node) {
        if (node.cache.cur_val <= best_val) {
            return;
        }
        if (node.cache.cur_val > best_val) {
            std::cout << "update time: " << (double)running_timer.cur_time().count() / 1e9 << "s\n";
            std::cout << "update: " << node.cache.cur_val << '\n';
            best_val = node.cache.cur_val;
            best_set = node.set;
        }
    }

   public:
    void output_json(const std::string& file_path) const {
        nlohmann::json json;
        std::string algorithm_name;
        switch (args::algorithm) {
            case args::Algorithm::BASIC:
                algorithm_name = "basic";
                json["algorithm"] = "basic";
                break;
            case args::Algorithm::DUAL:
                algorithm_name = "dual";
                json["algorithm"] = "dual";
                break;
            case args::Algorithm::BFS:
                algorithm_name = "bfs";
                json["algorithm"] = "bfs";
                break;
            case args::Algorithm::GREEDY:
                algorithm_name = "greedy";
                json["algorithm"] = "greedy";
                break;
            default:
                throw std::runtime_error("Invalid algorithm");
        }
        switch (args::upper_bound_type) {
            case args::UpperBoundType::KNAPSACK:
                algorithm_name += "-k";
                json["upper_bound"] = "knapsack";
                break;
            case args::UpperBoundType::FRAC_KNAPSACK:
                algorithm_name += "-fk";
                json["upper_bound"] = "frac_knapsack";
                break;
            case args::UpperBoundType::REFINE_SUBSET:
                algorithm_name += "-rs";
                json["upper_bound"] = "refine_subset";
                break;
            case args::UpperBoundType::DOM:
                algorithm_name += "-dom";
                json["upper_bound"] = "dom";
                break;
            default:
                throw std::runtime_error("Invalid upper bound type");
        }
        if (args::algorithm == args::Algorithm::GREEDY) {
            json["algorithm_name"] = "greedy";
        } else if (args::algorithm == args::Algorithm::BFS) {
            json["algorithm_name"] = "BFSTC";
        } else {
            json["algorithm_name"] = algorithm_name;
        }
        switch (args::submod_func_type) {
            case args::SubmodFuncType::DOM:
                json["submod_func"] = "dom";
                break;
            case args::SubmodFuncType::INF:
                json["submod_func"] = "inf";
                break;
            case args::SubmodFuncType::LOC:
                json["submod_func"] = "loc";
                break;
            case args::SubmodFuncType::COV:
                json["submod_func"] = "cov";
                break;
            default:
                throw std::runtime_error("Invalid submod function type");
        }
        switch (args::weight_type) {
            case args::WeightType::UNIT:
                json["weight_type"] = "unit";
                break;
            case args::WeightType::UNIFORM:
                json["weight_type"] = "uniform";
                break;
            case args::WeightType::NORMAL:
                json["weight_type"] = "normal";
                break;
            default:
                throw std::runtime_error("Invalid weight type");
        }
        switch (args::update_type) {
            case args::UpdateType::ALL:
                json["update_type"] = "all";
                break;
            case args::UpdateType::THRESHOLD:
                json["update_type"] = "threshold";
                break;
            case args::UpdateType::MAX:
                json["update_type"] = "max";
                break;
            default:
                throw std::runtime_error("Invalid update type");
        }
        switch (args::set_bound_update_type) {
            case args::SetBoundUpdateType::GLOBAL:
                json["set_bound_update_type"] = "global";
                break;
            case args::SetBoundUpdateType::LOCAL:
                json["set_bound_update_type"] = "local";
                break;
            default:
                throw std::runtime_error("Invalid set bound update type");
        }
        json["input_file"] = args::input_file;

        json["root_upper_bound"] = root_upper_bound;
        json["search_nodes"] = search_nodes;
        json["max_weight"] = max_weight;
        json["running_time"] = (double)running_timer.elapsed().count() / 1e9;
        json["early_stop"] = early_stop;
        json["weights"] = weights;

        json["best_val"] = best_val;
        json["best_set"] = best_set;

        // write to file
        std::ofstream file(file_path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + file_path);
        }
        file << json.dump(4);
        file.close();
        std::cout << "Output to " << file_path << '\n';
    }
};
