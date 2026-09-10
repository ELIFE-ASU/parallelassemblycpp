#!/bin/bash
#SBATCH --job-name=parallelassemblycpp-install
#SBATCH --partition=lightwork
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --export=NONE

# Submit from the repository root on ASU Sol:
#   sbatch slurm/install-sol.sbatch [ENV_PREFIX [REPO_DIR]]
# Optional arguments must be absolute paths. The default environment is
# ~/.conda/envs/parallelassemblycpp; REPO_DIR defaults to SLURM_SUBMIT_DIR because
# Slurm executes a spooled copy of this script, not the file in the checkout.
# Site guidance: https://docs.rc.asu.edu/mamba/
# https://docs.rc.asu.edu/partitions-and-qos/#lightwork

# Module and Conda activation scripts may reference unset shell variables.
# Enable nounset only after activation has completed.
set -eo pipefail

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submit this script with sbatch on Sol; installation needs a job allocation." >&2
    exit 1
fi
if (( $# > 2 )); then
    echo "Usage: sbatch slurm/install-sol.sbatch [ENV_PREFIX [REPO_DIR]]" >&2
    exit 1
fi

env_prefix="${1:-${HOME}/.conda/envs/parallelassemblycpp}"
repo_dir="${2:-${SLURM_SUBMIT_DIR:?Missing submission directory}}"
if [[ "$env_prefix" != /* || "$repo_dir" != /* ]]; then
    echo "ENV_PREFIX and REPO_DIR must be absolute paths." >&2
    exit 1
fi
if [[ ! -f "$repo_dir/environment.yml" || ! -f "$repo_dir/CMakePresets.json" ]]; then
    echo "Repository files not found in $repo_dir; submit from the repository root." >&2
    exit 1
fi
cd "$repo_dir"

module purge
module load mamba/latest

# Keep Mamba 1.x and 2.x noninteractive without changing the user's config.
export CONDA_ALWAYS_YES=true
export MAMBA_ALWAYS_YES=true
export CONDA_CHANNEL_PRIORITY=strict
export MAMBA_CHANNEL_PRIORITY=strict

if [[ -f "$env_prefix/conda-meta/history" ]]; then
    echo "Updating environment: $env_prefix"
    mamba env update --prefix "$env_prefix" --file "$repo_dir/environment.yml"
elif [[ -e "$env_prefix" ]]; then
    echo "Refusing to overwrite a path that is not a Conda environment: $env_prefix" >&2
    exit 1
else
    mkdir -p "$(dirname "$env_prefix")"
    echo "Creating environment: $env_prefix"
    mamba env create --prefix "$env_prefix" --file "$repo_dir/environment.yml"
fi

# ASU supports source activate for both named environments and absolute prefixes.
# shellcheck disable=SC1091
source activate "$env_prefix"
set -u

# Use a separate build tree so existing local compiler settings do not leak in.
# Refresh the CMake cache in case a repeated installation updates the compiler.
build_dir="$repo_dir/build/sol-release"
cmake --preset release --fresh -B "$build_dir"
cmake --build "$build_dir" --parallel "${SLURM_CPUS_PER_TASK:-1}"
"$build_dir/ParallelAssemblyCpp" --help

printf '\nEnvironment ready. In future Sol jobs, run:\n'
printf 'module load mamba/latest\nsource activate %q\n' "$env_prefix"
printf 'Executable: %s/ParallelAssemblyCpp\n' "$build_dir"
