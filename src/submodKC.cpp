#include <argparse/argparse.hpp>
#include <cstdint>
#include <iostream>

#include "algorithm/submodKC.hpp"
#include "structures/submod_function.hpp"
#include "utils/args.hpp"
#include "utils/data.hpp"
#include "test/test_BFSTC.hpp"

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
    SubmodKC<decltype(submod_func)> solver(std::move(submod_func), args::weight, std::move(weights),
                                              std::chrono::nanoseconds((int64_t)(args::time_limit_second * 1e9)));
    auto val = solver.solve();
    std::cout << "best_val: " << val << '\n';
    std::cout << "best_set: ";
    for (auto i : solver.best_set) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    // debug
    solver.check();
    if (args::output_file.empty()) {
        std::cout << "output file is empty\n";
        return;
    }
    solver.output_json(args::output_file);
}

int main(int argc, char** argv) {
    try {
        args::seq_submod_parse_args(argc, argv);
        if (args::algorithm == args::Algorithm::BFS) {
            // close lazy update
            args::update_type = args::UpdateType::ALL;
        }
        // test_BFSTC_COV();
        // return 0;
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
