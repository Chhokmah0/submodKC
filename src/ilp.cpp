#include <argparse/argparse.hpp>
#include <cstdint>
#include <iostream>

#include "algorithm/ilp.hpp"
#include "gurobi_c++.h"
#include "structures/submod_function.hpp"
#include "utils/args.hpp"
#include "utils/data.hpp"

template <typename T>
void run(T submod_func) {
    std::vector<Weight> weights;
    switch (args::weight_type) {
        case args::WeightType::UNIT:
            weights = utils::unit_weight(submod_func.max_index);
            break;
        case args::WeightType::UNIFORM:
            weights = utils::random_uniform_weight(submod_func.max_index, args::seed);
            break;
        case args::WeightType::NORMAL:
            weights = utils::random_normal_weight(submod_func.max_index, args::seed);
            break;
        default:
            throw std::runtime_error("Invalid weight type");
    }
    GRBEnv env;
    CGModel<decltype(submod_func)> cg_model(env, std::move(submod_func), std::move(weights), args::weight,
                                            args::time_limit_second);
    cg_model.solve();
    std::cout << "best_val: " << cg_model.best_val << '\n';
    std::cout << "best_set: ";
    for (auto i : cg_model.best_set) {
        std::cout << i << ' ';
    }
    std::cout << '\n';

    // check
    Weight weight = 0;
    for(int i : cg_model.best_set) {
        weight += cg_model.weights[i];
    }
    if (weight > args::weight) {
        std::cout << "best_set weight is not correct\n";
    }
    auto val = cg_model.submod_func(cg_model.best_set);
    if (std::abs(val - cg_model.best_val) > 1e-10) {
        std::cout << "best_set value is not correct\n";
        std::cout << "best_set value: " << val << '\n';
        std::cout << "model value: " << cg_model.model.get(GRB_DoubleAttr_ObjVal) << '\n';
    }
    if (args::output_file.empty()) {
        std::cout << "output file is empty\n";
        return;
    }
    cg_model.output_json(args::output_file);
}

int main(int argc, char** argv) {
    try {
        args::ilp_parse_args(argc, argv);
        switch (args::submod_func_type) {
            case args::SubmodFuncType::DOM:
                run(utils::DOM_read_from(args::input_file));
                break;
            case args::SubmodFuncType::INF:
                run(utils::INF_read_from(args::input_file));
                break;
            case args::SubmodFuncType::LOC:
                run(utils::LOC_read_from(args::input_file));
                break;
            case args::SubmodFuncType::COV:
                run(utils::COV_read_from(args::input_file));
                break;
            default:
                throw std::runtime_error("Invalid submod function type");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred.\n";
        return 1;
    }
    return 0;
}
