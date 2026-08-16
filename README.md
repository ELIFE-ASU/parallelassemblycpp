# AssemblyCpp v5

https://arxiv.org/abs/2410.09100#

This repository contains the C++ implementation of v5 of the assembly algorithm, in a "script" style version.

School of Chemistry, The University of Glasgow, University Avenue, Glasgow G12 8QQ, United Kingdom

Authors Ian Seet, Leroy Cronin


## Build with Conda

Install [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or another
Conda distribution, then run the following commands from the repository root:

```bash
conda env create --file environment.yml
conda activate assemblycpp-v5
mkdir -p build
c++ v5/main.cpp -std=c++23 -O3 \
  -I"${CONDA_PREFIX}/include" -o build/AssemblyCpp
```

The environment installs a C++ compiler, Boost (including the Boost Graph Library),
and Python for the test runner. Activating it places the Conda-provided compiler
on `PATH`. Recreate the environment after changing `environment.yml` with:

```bash
conda env remove --name assemblycpp-v5
conda env create --file environment.yml
```

To confirm that the executable was built successfully:

```bash
./build/AssemblyCpp --help
```

## Tests

The test runner can compile the executable and run a quick regression sample in
one command. Each molecule runs in an isolated temporary directory, so test
artifacts do not modify `unitTests/`:

```bash
python unitTests/unitTester.py --build --limit 20 --timeout 120
```

Run the complete 944-case battery with:

```bash
python unitTests/unitTester.py --build
```

Independent cases can run concurrently. Choose a worker count suitable for the
available CPU and memory:

```bash
python unitTests/unitTester.py --build --jobs 4
```

Use `python unitTests/unitTester.py --help` to see all options. The original
three positional arguments—executable, molecule list, and expected-results
file—remain supported.

GitHub Actions runs the quick regression sample for every push and pull request.
The complete battery can also be started manually from the Actions tab.
