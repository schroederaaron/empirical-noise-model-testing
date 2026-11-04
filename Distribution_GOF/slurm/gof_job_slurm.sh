#!/bin/bash
# gof_job_slurm.sh -- run ONE Distribution_GOF R script in the arch-gfortran container.
# Submitted by gof_yeast48_submit.sh; can also be submitted by hand from the
# Tensor-Omics root:
#
#   sbatch --job-name=gof_x --cpus-per-task=32 --mem=64G \
#          --output=analysis/Distribution_GOF/results/logs/%x_%j.out \
#          --error=analysis/Distribution_GOF/results/logs/%x_%j.err \
#          gof_job_slurm.sh scripts/02_run_gof.R --DATASET=... [--NAME=value ...]
#
# The container runs as YOUR uid (outputs stay yours), sees /media/BioNAS2 at the
# same path as the host, and uses the Tensor-Omics root as working dir, so
# external/docker_r_libs is the package library. N_CORES follows the allocation.
#SBATCH --ntasks=1
set -euo pipefail

GOF_DIR="$PWD/analysis/Distribution_GOF"
[[ -f "$GOF_DIR/config_gof.R" ]] || { echo "run from the Tensor-Omics root ($GOF_DIR not found)" >&2; exit 64; }
[[ $# -ge 1 ]] || { echo "usage: gof_job_slurm.sh scripts/<name>.R|tests [--NAME=value ...]" >&2; exit 64; }
target="$1"; shift

cpus="${SLURM_CPUS_PER_TASK:-1}"
mem_arg=(); [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && mem_arg=(--memory "${SLURM_MEM_PER_NODE}m")

if [[ "$target" == "tests" ]]; then
  cmd=(Rscript -e 'GOF_ROOT <- Sys.getenv("GOF_ROOT"); source(file.path(GOF_ROOT, "R", "setup.R")); load_or_install("testthat"); testthat::test_dir(file.path(GOF_ROOT, "tests", "testthat"), stop_on_failure = TRUE)')
else
  [[ "$target" == scripts/*.R && -f "$GOF_DIR/$target" ]] || { echo "not a Distribution_GOF script: $target" >&2; exit 64; }
  cmd=(Rscript "$GOF_DIR/$target" "$@" "--N_CORES=$cpus")   # last --N_CORES wins (config_gof.R)
fi

name="gof_${SLURM_JOB_ID:-local}_$$"
echo "gof_job: $(date -Is) job=${SLURM_JOB_ID:-local} node=$(hostname) cpus=$cpus target=$target args: $*" >&2
cleanup() { docker kill "$name" >/dev/null 2>&1 || true; }
trap cleanup TERM INT

docker run --rm --init --name "$name" \
  --user "$(id -u):$(id -g)" --cpus "$cpus" "${mem_arg[@]}" \
  -e HOME=/tmp -e OMP_NUM_THREADS=1 -e GOF_ROOT="$GOF_DIR" \
  -e GOF_GIT_HEAD="${GOF_GIT_HEAD:-analysis/ copy of empirical-noise-model-testing (no git checkout)}" \
  -v /media/BioNAS2:/media/BioNAS2 -w "$PWD" \
  arch-gfortran "${cmd[@]}" &
pid=$!
set +e; wait "$pid"; rc=$?; set -e
echo "gof_job: $(date -Is) exit=$rc" >&2
exit "$rc"
