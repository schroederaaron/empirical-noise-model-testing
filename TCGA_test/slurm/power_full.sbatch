#!/bin/bash
# power_full.sbatch -- FULL RUN: P+ and Part R (all six TCGA base cohorts), 20 rounds,
# then power_report.R. Writes to the DEFAULT locations ($ROOT/power_raw, $ROOT/power_out).
#
# Submit from the directory that holds logs/ (the job itself cd-s into $WORK), AFTER the smoke test looked sane:
#   mkdir -p logs && sbatch power_full.sbatch
#
# What overwrites what:
#   * power_test.R clears the raw directory of every part it runs (here: P+ and R) and nothing else.
#   * power_report.R rewrites $ROOT/power_out from ALL raw parts present in $ROOT/power_raw.
#   Copy $ROOT/power_out (and power_raw/P+) elsewhere first if you want to keep the placeholder-parameter run.
#
# TUNE (guesses, not measured): --mem and --time. Use the smoke test's sacct output to size them.
#SBATCH --job-name=power_full
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=main        # the TOX library is built with -march=native on the main node; it
                                # dies with "Illegal instruction" on the second (non-AVX) node
##SBATCH --account=<account>

set -euo pipefail

ROOT=/media/BioNAS2/TCGA_TOX_TEST
WORK=$ROOT/Tensor-Omics
IMAGE=/media/BioNAS2/apptainer/images/arch-gfortran.sif
S=analysis/TCGA_test/scripts
CPUS=${SLURM_CPUS_PER_TASK:-32}
PARTS="P+,R"                     # add more (e.g. "R,P,P+,Q,C,S,N,V,D") to regenerate other parts too

COHORTS="coad_healthy coad_stageiv luad_healthy luad_stageiv kirc_healthy kirc_stageiv"
REAL_RDS=""
for c in $COHORTS; do
  [[ -f "$WORK/power_cohort/$c.rds" ]] || { echo "power_full: missing $WORK/power_cohort/$c.rds" >&2; exit 2; }
  REAL_RDS+="${REAL_RDS:+,}power_cohort/$c.rds"
done
[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
for f in disp_params.rds $S/power_test.R $S/power_report.R; do
  [[ -f "$WORK/$f" ]] || { echo "power_full: missing $WORK/$f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"

echo "power_full: $(date -Is) host=$(hostname) cpus=$CPUS mem=${SLURM_MEM_PER_NODE:-?}MB parts=$PARTS commit=${GIT_COMMIT:-unknown}" >&2

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
  --env POWER_PARTS="$PARTS" \
  --env POWER_ROUNDS=20 \
  --env POWER_CORES="$CPUS" \
  --env POWER_DISP_PARAMS=disp_params.rds \
  --env POWER_REAL_RDS="$REAL_RDS" \
  ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
  "$IMAGE" bash -c "
    Rscript -e 'str(readRDS(\"disp_params.rds\"), max.level = 1)' &&
    Rscript $S/power_test.R &&
    Rscript $S/power_report.R"
rc=$?
set -e
echo "power_full: finished rc=$rc at $(date -Is)" >&2
exit "$rc"
