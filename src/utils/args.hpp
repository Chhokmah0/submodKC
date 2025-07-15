#pragma once
#include <argparse/argparse.hpp>

namespace args {

inline std::string input_file;
inline std::string output_file;
inline double weight = 5;

enum class WeightType : std::uint8_t { UNIT, UNIFORM, NORMAL };
inline WeightType weight_type = WeightType::UNIT;

enum class Algorithm : std::uint8_t {
    BASIC,
    DUAL,
    BFS,
    GREEDY,
};
inline Algorithm algorithm = Algorithm::DUAL;

enum class SubmodFuncType : std::uint8_t {
    DOM,
    INF,
    LOC,
    COV,
};
inline SubmodFuncType submod_func_type = SubmodFuncType::DOM;

inline double time_limit_second = 0;
inline int seed = 0;

enum class UpperBoundType : std::uint8_t {
    KNAPSACK,
    FRAC_KNAPSACK,
    REFINE_SUBSET,
    DOM,
};
inline UpperBoundType upper_bound_type = UpperBoundType::KNAPSACK;

enum class UpdateType : std::uint8_t {
    ALL,
    THRESHOLD,
    MAX,
};
inline UpdateType update_type = UpdateType::THRESHOLD;

enum class SetBoundUpdateType : std::uint8_t {
    GLOBAL,
    LOCAL,
};
inline SetBoundUpdateType set_bound_update_type = SetBoundUpdateType::GLOBAL;

void seq_submod_parse_args(int argc, char** argv) {
    argparse::ArgumentParser program("par_submod");

    program.add_argument("-i", "--input").help("input file").required().action([](const std::string& value) {
        input_file = value;
    });
    program.add_argument("-o", "--output").help("output file").action([](const std::string& value) {
        output_file = value;
    });
    program.add_argument("-w", "--max-weight").help("max weight").store_into(weight);
    program.add_argument("-g", "--gen-weight").help("gen weight type").action([](const std::string& value) {
        if (value == "unit") {
            weight_type = WeightType::UNIT;
        } else if (value == "uniform") {
            weight_type = WeightType::UNIFORM;
        } else if (value == "normal") {
            weight_type = WeightType::NORMAL;
        } else {
            throw std::runtime_error("Invalid weight type");
        }
    });
    program.add_argument("-a", "--algorithm")
        .help("algorithm")
        .choices("basic-k", "basic-fk", "basic-rs", "basic-dom", "dual-rs", "BFSTC", "greedy")
        .action([](const std::string& value) {
            if (value == "greedy") {
                algorithm = Algorithm::GREEDY;
                return;
            }
            if (value == "BFSTC") {
                algorithm = Algorithm::BFS;
                upper_bound_type = UpperBoundType::DOM;
                return;
            }
            // split at the first '-'
            std::string algorithm_name = value.substr(0, value.find('-'));
            if (algorithm_name == "basic") {
                algorithm = Algorithm::BASIC;
            } else if (algorithm_name == "dual") {
                algorithm = Algorithm::DUAL;
            }else if (algorithm_name == "greedy") {
                algorithm = Algorithm::GREEDY;
            } else {
                throw std::runtime_error("Invalid algorithm");
            }
            std::string algorithm_type = value.substr(value.find('-') + 1);
            if (algorithm_type == "k") {
                upper_bound_type = UpperBoundType::KNAPSACK;
            } else if (algorithm_type == "fk") {
                upper_bound_type = UpperBoundType::FRAC_KNAPSACK;
            } else if (algorithm_type == "rs") {
                upper_bound_type = UpperBoundType::REFINE_SUBSET;
            } else if (algorithm_type == "dom") {
                upper_bound_type = UpperBoundType::DOM;
            } else {
                throw std::runtime_error("Invalid upper bound type");
            }
        });

    program.add_argument("-f", "--submod_func")
        .help("submod function type")
        .required()
        .action([](const std::string& value) {
            if (value == "dom") {
                submod_func_type = SubmodFuncType::DOM;
            } else if (value == "inf") {
                submod_func_type = SubmodFuncType::INF;
            } else if (value == "loc") {
                submod_func_type = SubmodFuncType::LOC;
            } else if (value == "cov") {
                submod_func_type = SubmodFuncType::COV;
            } else {
                throw std::runtime_error("Invalid submod function type");
            }
        });
    program.add_argument("-u", "--update").help("update type").action([](const std::string& value) {
        if (value == "all") {
            update_type = UpdateType::ALL;
        } else if (value == "threshold") {
            update_type = UpdateType::THRESHOLD;
        } else if (value == "max") {
            update_type = UpdateType::MAX;
        } else {
            throw std::runtime_error("Invalid update type");
        }
    });
    program.add_argument("-s", "--seed").help("random seed").action([](const std::string& value) {
        seed = std::stoi(value);
    });
    program.add_argument("-t", "--time").help("time limit in seconds").store_into(time_limit_second);
    program.add_argument("-su", "--set-bound-update")
        .help("set bound update type")
        .action([](const std::string& value) {
            if (value == "global") {
                set_bound_update_type = SetBoundUpdateType::GLOBAL;
            } else if (value == "local") {
                set_bound_update_type = SetBoundUpdateType::LOCAL;
            } else {
                throw std::runtime_error("Invalid set bound update type");
            }
        });

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << '\n';
        std::cout << program;
        exit(0);
    }
}

void ilp_parse_args(int argc, char** argv) {
    argparse::ArgumentParser program("ilp");

    program.add_argument("-i", "--input").help("input file").required().action([](const std::string& value) {
        input_file = value;
    });
    program.add_argument("-o", "--output").help("output file").action([](const std::string& value) {
        output_file = value;
    });
    program.add_argument("-w", "--max-weight").help("max weight").store_into(weight);
    program.add_argument("-g", "--gen-weight").help("gen weight type").action([](const std::string& value) {
        if (value == "unit") {
            weight_type = WeightType::UNIT;
        } else if (value == "uniform") {
            weight_type = WeightType::UNIFORM;
        } else if (value == "normal") {
            weight_type = WeightType::NORMAL;
        } else {
            throw std::runtime_error("Invalid weight type");
        }
    });
    program.add_argument("-s", "--seed").help("random seed").action([](const std::string& value) {
        seed = std::stoi(value);
    });

    program.add_argument("-f", "--submod_func")
        .help("submod function type")
        .required()
        .action([](const std::string& value) {
            if (value == "dom") {
                submod_func_type = SubmodFuncType::DOM;
            } else if (value == "inf") {
                submod_func_type = SubmodFuncType::INF;
            } else if (value == "loc") {
                submod_func_type = SubmodFuncType::LOC;
            } else if (value == "cov") {
                submod_func_type = SubmodFuncType::COV;
            } else {
                throw std::runtime_error("Invalid submod function type");
            }
        });
    program.add_argument("-t", "--time").help("time limit in seconds").store_into(time_limit_second);

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << '\n';
        std::cout << program;
        exit(0);
    }
}
}  // namespace args
