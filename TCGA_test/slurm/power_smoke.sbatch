#!/bin/bash
# power_smoke.sbatch -- SMOKE TEST: P+ (1000 genes) + Part R on two cohorts, 2 rounds,
# then power_report.R on the smoke output. Writes ONLY to power_raw_smoke / power_out_smoke.
#
# Submit from the directory that holds logs/ (the job itself cd-s into $WORK):   mkdir -p logs && sbatch power_smoke.sbatch
#
# TUNE (guesses, not measured): --cpus-per-task, --mem, --time. After the run:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,AllocCPUS,State
#SBATCH --job-name=power_smoke
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=main        # the TOX library is built with -march=native on the main node; it
                                # dies with "Illegal instruction" on the second (non-AVX) node
##SBATCH --account=<account>

set -euo pipefail

ROOT=/media/BioNAS2/TCGA_TOX_TEST                                   # mounted at the same path (default outputs live here)
WORK=$ROOT/Tensor-Omics                                              # working dir: holds analysis/, rcpp/, power_cohort/, disp_params.rds
IMAGE=/media/BioNAS2/apptainer/images/arch-gfortran.sif
S=analysis/TCGA_test/scripts
CPUS=${SLURM_CPUS_PER_TASK:-8}

# --- fail fast on missing inputs (paths relative to $WORK) -------------------------------------
[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
for f in disp_params.rds power_cohort/kirc_healthy.rds power_cohort/kirc_stageiv.rds \
         $S/power_test.R $S/power_report.R; do
  [[ -f "$WORK/$f" ]] || { echo "power_smoke: missing $WORK/$f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"   # optional provenance

echo "power_smoke: $(date -Is) host=$(hostname) cpus=$CPUS mem=${SLURM_MEM_PER_NODE:-?}MB commit=${GIT_COMMIT:-unknown}" >&2

# Runs in Apptainer (group image, read-only), as you, inside this Slurm job: the job's
# cores and memory apply to it and scancel stops it. --cleanenv passes ONLY the --env
# variables below; --home gives R an empty writable home (your ~/.Rprofile and ~/R stay
# out); -B /media/BioNAS2 makes BioNAS2 visible at the same path, so every relative path
# resolves from $WORK exactly as before.
# OMP_NUM_THREADS=1: avoids threads x mclapply workers oversubscribing the allocation.
cd "$WORK"
mkdir -p .container-home

set +e
apptainer exec --cleanenv --home "$WORK/.container-home" -B /media/BioNAS2 \
  --env OMP_NUM_THREADS=1 \
  --env POWER_PARTS="P+,R" \
  --env POWER_ROUNDS=2 \
  --env POWER_GENES=1000 \
  --env POWER_CORES="$CPUS" \
  --env POWER_DISP_PARAMS=disp_params.rds \
  --env POWER_REAL_RDS=power_cohort/kirc_healthy.rds,power_cohort/kirc_stageiv.rds \
  --env POWER_RAW_DIR="$ROOT/power_raw_smoke" \
  --env POWER_OUT_DIR="$ROOT/power_out_smoke" \
  ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
  "$IMAGE" bash -c "
    Rscript -e 'str(readRDS(\"disp_params.rds\"), max.level = 1)' &&
    Rscript $S/power_test.R &&
    Rscript $S/power_report.R"
rc=$?
set -e
echo "power_smoke: finished rc=$rc at $(date -Is)" >&2
exit "$rc"
