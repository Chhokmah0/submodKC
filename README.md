# submodKC
Exact algorithm for the submodular knapsack problem.

## Prerequist 
1. Install Xmake
Xmake is a cross-platform build utility based on the Lua scripting language.

With cURL
```bash
curl -fsSL https://xmake.io/shget.text | bash
```

With Wget
```bash
curl -fsSL https://xmake.io/shget.text | bash
```

With PowerShell
```sh
irm https://xmake.io/psget.text | iex
```
2. Gurobi Version 12.01.

Gurobi is a commercial mixed integer linear solver, see https://www.gurobi.com/.
It is known that a free license can be applied for academic purpose.



## Build
The first build will require downloading the `argparse` and `nlohmann_json` libraries.
```bash
xmake f -m release
xmake b submodKC
```

You can change the compiler used using the following command line. 

```bash
xmake f --toolchain=clang
```

If you want the ILP, please set the Gurobi path in `xmake.lua` and build it with the following command:

```bash
xmake b ilp
```

## Command-line arguments
You can run the program with `xmake r submodKC`. Or find the binary in path like `./build/linux/x86_64/release`.

- `-i`, `--input`: The input file.
- `-o`, `--ouput`: The output `.json` file.
- `-a`, `--algorithm`: Algorithm. \[`basic-k`, `basic-fk`, `basic-rs`, `basic-dom`, `dual-rs`, `BFSTC`, `greedy`\]
- `-f`, `--submod_func`: Submod function type. \[`dom`,`inf`,`loc`,`cov`\].
- `-w`, `--max-weight`: $W$.
- `-g`, `--gen-weight`: Weight distribution. \[`normal`, `uniform`, `unit`\]
- `-s`, `--seed`: The seed for weight generation.
- `-t`, `--time`: Time limit in seconds.

# Citation
If you use our codes, please cite our paper
"Efficient Branch-and-Bound for Submodular Function Maximization under Knapsack Constraint."
