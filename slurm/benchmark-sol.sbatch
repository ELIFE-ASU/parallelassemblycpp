#!/bin/bash
#SBATCH --job-name=assemblycpp-benchmarks
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --export=NONE

# Create the environment with install-sol.sbatch first. Submit from the repo:
#   sbatch slurm/benchmark-sol.sbatch
# Reserve a standard Sol CPU node: 128 physical cores and all node memory.
# Hardware: https://docs.rc.asu.edu/supercomputer-hardware/
# Sol CPU jobs: https://docs.rc.asu.edu/partitions-and-qos/
# Keep nounset off while loading site modules and Conda activation scripts.
set -eo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch [SLURM OPTIONS] slurm/benchmark-sol.sbatch [OPTIONS]
  --env-prefix PATH   Environment (default: ~/.conda/envs/assemblycpp-v5)
  --repo-dir PATH     Repository (default: SLURM_SUBMIT_DIR)
  --output-dir PATH   New result directory (default: build/sol-benchmarks/JOB_ID)
  --suites LIST       Comma-separated suites (default: full,profile,scaling)
  --threads LIST     OpenMP thread counts, e.g. 2,4,8 (default: powers of two
                     up to allocated CPUs, also including the allocation size)
  --runs N           Measured rounds per suite/thread count (default: 6)
  --warmup N         Warm-up rounds (default: 1)
  --timeout SECONDS  Per-calculation timeout (default: 600)
  --help             Show this help
All paths must be absolute. Use a new output directory for every job.
The thread sweep measures Paclitaxel against paired serial runs.
Default allocation: one exclusive Sol node, 128 physical cores, all node memory.
EOF
}

env_prefix="${HOME}/.conda/envs/assemblycpp-v5"
repo_dir="${SLURM_SUBMIT_DIR:-$PWD}"
output_dir=""
suites="full,profile,scaling"
thread_list=""
runs=6
warmup=1
timeout=600
while (( $# )); do
    if [[ "$1" == --help ]]; then
        usage
        exit 0
    fi
    if (( $# < 2 )); then
        echo "Missing value for $1" >&2
        exit 1
    fi
    case "$1" in
        --env-prefix) env_prefix="$2" ;;
        --repo-dir) repo_dir="$2" ;;
        --output-dir) output_dir="$2" ;;
        --suites) suites="$2" ;;
        --threads) thread_list="$2" ;;
        --runs) runs="$2" ;;
        --warmup) warmup="$2" ;;
        --timeout) timeout="$2" ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
    shift 2
done

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "Submit this script with sbatch on Sol; benchmarks need a compute allocation." >&2
    exit 1
fi
cpus="${SLURM_CPUS_PER_TASK:?Missing allocated CPU count}"
if [[ ! "$cpus" =~ ^[1-9][0-9]*$ ]] || (( cpus < 2 )); then
    echo "Request at least two CPUs with sbatch --cpus-per-task=N." >&2
    exit 1
fi
if [[ "${SLURM_JOB_NUM_NODES:-1}" != 1 || "${SLURM_NTASKS:-1}" != 1 ]]; then
    echo "This job measures single-node OpenMP scaling; request one node and one task." >&2
    exit 1
fi
output_dir="${output_dir:-$repo_dir/build/sol-benchmarks/$SLURM_JOB_ID}"
for path in "$env_prefix" "$repo_dir" "$output_dir"; do
    if [[ "$path" != /* ]]; then
        echo "Environment, repository, and output paths must be absolute: $path" >&2
        exit 1
    fi
done
if [[ ! -f "$repo_dir/benchmarks/paclitaxel_scaling.py" ||
      ! -f "$repo_dir/benchmarks/summarize_sol.py" ||
      ! -f "$repo_dir/CMakePresets.json" ]]; then
    echo "Repository files not found in $repo_dir; submit from the root or use --repo-dir." >&2
    exit 1
fi
if [[ ! -f "$env_prefix/conda-meta/history" ]]; then
    echo "Environment not found: $env_prefix. Run slurm/install-sol.sbatch first." >&2
    exit 1
fi
if [[ ! "$suites" =~ ^(quick|full|profile|scaling)(,(quick|full|profile|scaling))*$ ]]; then
    echo "--suites must contain comma-separated quick, full, profile, or scaling." >&2
    exit 1
fi
IFS=, read -r -a suite_names <<< "$suites"
declare -A selected_suites=()
for suite in "${suite_names[@]}"; do
    if [[ -n "${selected_suites[$suite]:-}" ]]; then
        echo "Duplicate suite: $suite" >&2
        exit 1
    fi
    selected_suites[$suite]=1
done
if [[ -n "$thread_list" ]]; then
    if [[ ! "$thread_list" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]]; then
        echo "--threads must be a comma-separated list such as 2,4,8." >&2
        exit 1
    fi
    IFS=, read -r -a thread_counts <<< "$thread_list"
else
    thread_counts=()
    for (( count=2; count<=cpus; count*=2 )); do
        thread_counts+=("$count")
    done
    if (( thread_counts[${#thread_counts[@]}-1] != cpus )); then
        thread_counts+=("$cpus")
    fi
fi

cd "$repo_dir"
module purge
module load mamba/latest
# shellcheck disable=SC1091
source activate "$env_prefix"
set -u
# srun must export the activated environment despite the job's --export=NONE.
export SLURM_EXPORT_ENV=ALL
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

# Check dependencies, arguments, and actual step affinity before compiling or
# measuring. The existing sweep validates thread counts and repetition limits.
python -c 'import matplotlib; import matplotlib.backends.backend_agg'
build_dir="$output_dir/build"
scaling_args=(
    --build-dir "$build_dir" --output-dir "$output_dir/parallel"
    --threads "${thread_counts[@]}" --runs "$runs" --warmup "$warmup"
    --timeout "$timeout"
)
srun --ntasks=1 --cpus-per-task="$cpus" --cpu-bind=cores \
    python benchmarks/paclitaxel_scaling.py "${scaling_args[@]}" --dry-run

mkdir -p "$(dirname "$output_dir")"
if ! mkdir "$output_dir"; then
    echo "Choose a new --output-dir; existing results will not be overwritten." >&2
    exit 1
fi
mkdir "$output_dir/suites"
printf 'running\n' > "$output_dir/status.txt"
finish() {
    local result=$?
    trap - EXIT
    # The sweep saves each count's JSON immediately, including before a later
    # count fails. Refresh the CSV on either success or failure.
    if ! python benchmarks/summarize_sol.py "$output_dir"; then
        if (( result == 0 )); then
            result=1
        fi
    fi
    if (( result == 0 )); then
        printf 'completed\n' > "$output_dir/status.txt"
    else
        printf 'failed (exit %s); summary.csv contains completed reports only\n' \
            "$result" > "$output_dir/status.txt"
    fi
    exit "$result"
}
trap finish EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

{
    printf 'Repository: %s\nEnvironment: %s\nResults: %s\n' \
        "$repo_dir" "$env_prefix" "$output_dir"
    printf 'Suites: %s\nThreads: %s\nRuns: %s\nWarmup: %s\nTimeout: %s\n' \
        "$suites" "${thread_counts[*]}" "$runs" "$warmup" "$timeout"
    date --utc --iso-8601=seconds
    hostname
    git rev-parse HEAD || true
    git status --short || true
    (set +u; module list)
    python --version
    cmake --version
    lscpu
    scontrol show job "$SLURM_JOB_ID"
} > "$output_dir/metadata.txt" 2>&1

cmake --preset performance --fresh -B "$build_dir" \
    -DASSEMBLYCPP_BUILD_OPENMP=ON 2>&1 | tee "$output_dir/build.log"
cmake --build "$build_dir" --parallel "$cpus" \
    --target AssemblyCpp AssemblyCppOMP 2>&1 | tee -a "$output_dir/build.log"

# Each suite uses one pinned CPU. The sweep below runs inside the full allocated
# CPU set and pins its serial reference consistently across all thread counts.
for suite in "${suite_names[@]}"; do
    srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
        python benchmarks/benchmark.py --executable "$build_dir/AssemblyCpp" \
        --candidate-parallel off --suite "$suite" \
        --runs "$runs" --warmup "$warmup" --timeout "$timeout" \
        --json-output "$output_dir/suites/$suite.json" \
        2>&1 | tee "$output_dir/suites/$suite.txt"
    python benchmarks/summarize_sol.py "$output_dir"
done

srun --ntasks=1 --cpus-per-task="$cpus" --cpu-bind=cores \
    python benchmarks/paclitaxel_scaling.py "${scaling_args[@]}" \
    2>&1 | tee "$output_dir/parallel.txt"

printf '\nResults: %s\nQuick analysis: %s/summary.csv\n' "$output_dir" "$output_dir"
