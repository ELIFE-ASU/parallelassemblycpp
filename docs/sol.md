# ASU Sol HPC setup

[Project overview](../README.md) · [Command line](cli.md) · [Build and tests](development.md)

Run source-tree examples from the repository root unless stated otherwise.

From the repository root on Sol, submit the environment installation job:

```bash
sbatch slurm/install-sol.sbatch
```

The job loads Sol's `mamba/latest` module and creates
`$HOME/.conda/envs/parallelassemblycpp` from `environment.yml`. Resubmitting updates
that environment to satisfy the file. It then builds the `release` preset in
`build/sol-release` and checks `ParallelAssemblyCpp --help`. The environment
includes the compiler, CMake, Ninja, Open MPI, Python, Matplotlib, and Ruff; the
release executable uses serial search. See the [parallel guide](parallel.md)
for parallel build presets.

The script requests one node, four CPUs, 16 GB RAM, and two hours in
`lightwork` with the `public` QoS, following ASU's guidance for
[environment creation and compilation](https://docs.rc.asu.edu/partitions-and-qos/#lightwork).
Output and errors go to `slurm-parallelassemblycpp-install-<job-id>.out` in the job's
working directory. An account can be selected with
`sbatch --account=<your-account> slurm/install-sol.sbatch`; use `myaccounts`
on Sol to list available accounts.

To choose another persistent environment location, pass an absolute prefix:

```bash
sbatch slurm/install-sol.sbatch /data/your_group/envs/parallelassemblycpp
```

An optional second argument supplies the absolute repository path when
submitting from another directory. Arguments are used because the script's
`--export=NONE` does not inherit custom variables from the submitting shell.
For example:

```bash
sbatch /path/to/parallelassemblycpp/slurm/install-sol.sbatch \
  /data/your_group/envs/parallelassemblycpp /path/to/parallelassemblycpp
```

After the setup job succeeds, activate the environment in subsequent jobs
using [ASU's supported activation syntax](https://docs.rc.asu.edu/mamba/):

```bash
module load mamba/latest
source activate "$HOME/.conda/envs/parallelassemblycpp"
./build/sol-release/ParallelAssemblyCpp unitTests/alanine.mol
```

Use your chosen prefix in `source activate` if you changed the default, and
run the executable from the repository root inside a compute allocation.
Finish jobs using the environment before resubmitting the installer, since
an update can replace their dependencies.

To run the benchmark suites and OpenMP plus MPI/OpenMP hybrid scaling sweeps on
a whole 128-core Sol node, submit `sbatch slurm/benchmark-sol.sbatch`. Results
include a combined `summary.csv`, raw JSON, and plots in
`build/sol-benchmarks/<job-id>/`. See the
[Sol benchmark job instructions](../benchmarks/README.md#asu-sol-batch-job) for
resource overrides, thread counts, and output details.
