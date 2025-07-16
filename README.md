# submodKC
This is the code for paper "Efficient Branch-and-Bound for Submodular Function Maximization under Knapsack Constraint"

The work will be reported at the 28th European Conference on Artificial Intelligence (ECAI 2025). 

## Prerequist 
1. Install Xmake
Xmake is a cross-platform build utility based on the Lua scripting language.

With cURL
```bash
curl -fsSL https://xmake.io/shget.text | bash
```

With Wget
```bash
wget https://xmake.io/shget.text -O - | bash
```

With PowerShell
```sh
Invoke-Expression (Invoke-Webrequest 'https://xmake.io/psget.text' -UseBasicParsing).Content
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
```
@misc 
{hao2025efficient,
      title={Efficient Branch-and-Bound for Submodular Function Maximization under Knapsack Constraint}, 
      author={Yimin Hao and Yi Zhou and Chao Xu and Zhang-Hua Fu},
      year={2025},
      eprint={2507.11107},
      archivePrefix={arXiv},
      primaryClass={cs.DS},
      url={https://arxiv.org/abs/2507.11107}, 
}
```
