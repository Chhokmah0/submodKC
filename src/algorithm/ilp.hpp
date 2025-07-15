#pragma once

#include <gurobi_c++.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <memory>
#include <nlohmann/json.hpp>

#include "gurobi_c.h"
#include "structures/submod_function.hpp"
#include "utils/args.hpp"
#include "utils/timer.hpp"

// using SubmodFunc = DOMFunction;
template <typename SubmodFunc>
class CGModel {
    using Val = SubmodFunc::EvalVal;
    using Cache = SubmodFunc::EvalCache;

    struct MyCallback : public GRBCallback {
        CGModel<SubmodFunc>& model;

        MyCallback(CGModel<SubmodFunc>& _cg_model) : model(_cg_model) {}

        void callback() override {
            if (where == GRB_CB_MIPSOL) {
                double z_val = getSolution(model.z);
                // assert(z_val == obj_val);

                auto x_val = std::unique_ptr<double[]>(getSolution(model.x.get(), model.submod_func.max_index));

                Node node{
                    .set = {},
                    .candidates = Candidates{},
                    .cache = model.submod_func.init_cache(),
                    .cur_weight = 0,
                };
                for (int i = 0; i < model.submod_func.max_index; ++i) {
                    if (x_val[i] > 0.5) {
                        auto margin = model.submod_func.eval_margin(node.cache, node.set, i);
                        node.set.push_back(i);
                        model.submod_func.update_cache(node.cache, i, margin);
                        node.cur_weight += model.weights[i];
                    }
                }
                if (node.cache.cur_val > model.best_val) {
                    model.best_val = node.cache.cur_val;
                    model.best_set = node.set;
                }
                if (model.best_val >= z_val) {
                    return;
                }

                GRBLinExpr submod_expr;
                submod_expr = node.cache.cur_val;
                for (int i = 0; i < model.submod_func.max_index; ++i) {
                    if (x_val[i] < 0.5) {
                        submod_expr += model.x[i] * model.submod_func.eval_margin(node.cache, node.set, i);
                    }
                }

                addLazy(model.z <= submod_expr);
            }
        }
    };

   public:
    SubmodFunc submod_func;
    std::vector<Weight> weights;
    Weight max_weight;

    GRBModel model;
    GRBVar z;
    std::unique_ptr<GRBVar[]> x;

    MyCallback callback;

    // 统计数据
    utils::Timer<> running_timer;
    double best_val;
    std::vector<Index> best_set;

    double root_upper_bound;

   public:
    CGModel(GRBEnv& env, SubmodFunc _submod_func, std::vector<Weight> _weights, Weight _max_weight,
            double time_limit_second = 0)
        : submod_func(std::move(_submod_func)),
          weights(std::move(_weights)),
          max_weight(_max_weight),
          model(env),
          z(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS, "z")),
          x(model.addVars(submod_func.max_index, GRB_BINARY)),
          callback(*this),
          best_val(submod_func.empty_set_val) {
        model.set(GRB_IntParam_LazyConstraints, 1);
        model.set(GRB_IntParam_Threads, 1);
        if (time_limit_second > 0) {
            model.set(GRB_DoubleParam_TimeLimit, time_limit_second);
        }
        model.set(GRB_DoubleParam_MemLimit, 32);
        model.setCallback(&callback);
        GRBLinExpr obj_expr = z;
        model.setObjective(obj_expr, GRB_MAXIMIZE);

        GRBLinExpr weight_expr;
        for (int i = 0; i < submod_func.max_index; ++i) {
            weight_expr += weights[i] * x[i];
        }
        model.addConstr(weight_expr <= max_weight, "weight");

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

        GRBLinExpr submod_expr;
        submod_expr = root.cache.cur_val;
        for (int i = 0; i < root.candidates.data.size(); ++i) {
            int id = root.candidates.data[i].id;
            submod_expr += x[id] * root.candidates.data[i].margin;
        }
        model.addConstr(z <= submod_expr);

        while (!root.candidates.empty()) {
            auto cand = root.candidates.pop_back();
            if (root.cur_weight + cand.weight > max_weight) {
                break;
            }
            choose(root, cand);
            update_margin(root);
            submod_expr = root.cache.cur_val;
            for (int i = 0; i < root.candidates.data.size(); ++i) {
                int id = root.candidates.data[i].id;
                submod_expr += x[id] * root.candidates.data[i].margin;
            }
            model.addConstr(z <= submod_expr);
        }
        model.update();
    }

    void solve() {
        auto relaxed_model = model.relax();
        relaxed_model.set(GRB_IntParam_LazyConstraints, 0);
        relaxed_model.optimize();
        root_upper_bound = relaxed_model.get(GRB_DoubleAttr_ObjVal);

        running_timer.start();
        model.optimize();
        running_timer.stop();
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
    };

    struct Node {
        std::vector<Index> set;
        Candidates candidates;
        Cache cache;
        Weight cur_weight;
    };

    void update_margin(Node& node) {
        for (auto& cand : node.candidates.data) {
            if (!cand.acc) {
                cand.margin = submod_func.eval_margin(node.cache, node.set, cand.id);
                cand.acc = true;
            }
        }
        node.candidates.sort();
    }

    void choose(Node& node, const Candidate& cand) {
        assert(cand.acc);
        node.set.push_back(cand.id);
        submod_func.update_cache(node.cache, cand.id, cand.margin);
        node.cur_weight += cand.weight;
        node.candidates.clear_acc();
    }

   public:
    void output_json(const std::string& file_path) const {
        nlohmann::json json;
        std::string algorithm_name;
        json["algorithm"] = "ilp";
        json["upper_bound"] = "ilp";
        json["algorithm_name"] = "ilp";
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
        json["input_file"] = args::input_file;

        json["root_upper_bound"] = root_upper_bound;
        json["max_weight"] = max_weight;
        json["running_time"] = (double)running_timer.elapsed().count() / 1e9;
        json["early_stop"] = model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT || model.get(GRB_IntAttr_Status) == GRB_MEM_LIMIT || 
                               model.get(GRB_IntAttr_Status) == GRB_INTERRUPTED;
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