#include <argparse/argparse.hpp>
#include <cstdint>
#include <iostream>

#include "algorithm/submodKC.hpp"
#include "structures/submod_function.hpp"
#include "utils/args.hpp"
#include "utils/data.hpp"

void test_BFSTC_COV() {
    args::update_type = args::UpdateType::ALL;
    args::submod_func_type = args::SubmodFuncType::COV;
    args::weight_type = args::WeightType::UNIFORM;
    args::weight = 1.25;
    auto cov = utils::COV_read_from_Facebook_like("data/OF_two-mode_weightedmsg.txt");
    auto weights = utils::dom_weight(cov.max_index, 0);
    SubmodKC<decltype(cov)> solver(std::move(cov), args::weight, weights,
                                      std::chrono::nanoseconds((int64_t)(args::time_limit_second * 1e9)));
    auto val = solver.solve();
    std::cout << "best_val: " << val << '\n';
    std::cout << "best_set: ";
    for (auto i : solver.best_set) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    solver.output_json(args::output_file);
}

void test_BFSTC_LOC() {
    args::update_type = args::UpdateType::ALL;
    args::submod_func_type = args::SubmodFuncType::LOC;
    args::weight_type = args::WeightType::UNIFORM;
    args::weight = 1.5;
    auto loc = utils::LOC_read_from_MTA_Swubway_Stations("data/MTA_Subway_Stations_20241122.csv");
    auto weights = utils::dom_weight(loc.max_index, 0);
    SubmodKC<decltype(loc)> solver(std::move(loc), args::weight, weights,
                                      std::chrono::nanoseconds((int64_t)(args::time_limit_second * 1e9)));
    auto val = solver.solve();
    std::cout << "best_val: " << val << '\n';
    std::cout << "best_set: ";
    for (auto i : solver.best_set) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    solver.output_json(args::output_file);
}

void test_BFSTC_INF() {
    args::update_type = args::UpdateType::ALL;
    args::submod_func_type = args::SubmodFuncType::LOC;
    args::weight_type = args::WeightType::UNIFORM;
    args::weight = 1;
    auto inf = utils::INF_read_from_Movie_100k("data/ml-100k/u.data");
    auto weights = utils::dom_weight(inf.max_index, 0);
    SubmodKC<decltype(inf)> solver(std::move(inf), args::weight, weights,
                                      std::chrono::nanoseconds((int64_t)(args::time_limit_second * 1e9)));
    auto val = solver.solve();
    std::cout << "best_val: " << val << '\n';
    std::cout << "best_set: ";
    for (auto i : solver.best_set) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    solver.output_json(args::output_file);
}
