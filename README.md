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
c++ v5/main.cpp -std=c++23 -O3 \
  -I"${CONDA_PREFIX}/include" -o AssemblyCpp
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
./AssemblyCpp --help
```

### Run the tests

From the repository root, with the Conda environment active:

```bash
cd unitTests
python unitTester.py ../AssemblyCpp batteryTest2 batteryTest2Base
```
