#!/bin/bash
# knn_config.sbatch -- neighbourhood configuration grid (knn_config_test.R): the four
# neighbourhoods of noise_model_variance (mean, band, adaptive, 2-D ellipse) x no
# scaling / cross-fit JS, at several sizes, on Part RS (six TCGA base cohorts) and
# Part P+ (simulated), n = 3, 5, 10, 20 rounds; then power_report.R and the kNN summary.
# Writes ONLY to $ROOT/knn_raw / $ROOT/knn_out (the normal power_out is not touched).
#
# Submit from the directory that holds logs/:   mkdir -p logs && sbatch knn_config.sbatch
# Subsets: sbatch --export=ALL,KNN_FAMILIES="mean,ellipse",KNN_PARTS=RS knn_config.sbatch
#
# TUNE (guesses): --time / --mem. About 50 arms per dataset; RS 720 + P+ 480 datasets.
#SBATCH --job-name=knn_config
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=main        # the TOX library is built with -march=native on the main node
##SBATCH --account=<account>

set -euo pipefail

ROOT=/media/BioNAS2/TCGA_TOX_TEST
WORK=$ROOT/Tensor-Omics
IMAGE=/media/BioNAS2/apptainer/images/arch-gfortran.sif
S=analysis/TCGA_test/scripts
CPUS=${SLURM_CPUS_PER_TASK:-32}
COHORTS=${RS_COHORTS:-"coad_healthy coad_stageiv luad_healthy luad_stageiv kirc_healthy kirc_stageiv"}

[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
REAL_RDS=""
for c in $COHORTS; do
  [[ -f "$WORK/power_cohort/$c.rds" ]] || { echo "knn_config: missing $WORK/power_cohort/$c.rds" >&2; exit 2; }
  REAL_RDS+="${REAL_RDS:+,}power_cohort/$c.rds"
done
for f in disp_params.rds $S/knn_config_test.R $S/power_test.R $S/power_report.R; do
  [[ -f "$WORK/$f" ]] || { echo "knn_config: missing $WORK/$f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"

echo "knn_config: $(date -Is) host=$(hostname) cpus=$CPUS parts=${KNN_PARTS:-RS,P+} families=${KNN_FAMILIES:-all} commit=${GIT_COMMIT:-unknown}"

# Apptainer as in power_full.sbatch: --cleanenv passes only the --env variables.
cd "$WORK"
mkdir -p .container-home

set +e
apptainer exec --cleanenv --home "$WORK/.container-home" -B /media/BioNAS2 \
  --env OMP_NUM_THREADS=1 \
  --env POWER_CORES="$CPUS" \
  --env POWER_DISP_PARAMS=disp_params.rds \
  --env POWER_REAL_RDS="$REAL_RDS" \
  --env KNN_PARTS="${KNN_PARTS:-RS,P+}" \
  --env KNN_ROUNDS="${KNN_ROUNDS:-20}" \
  --env KNN_N="${KNN_N:-3,5,10}" \
  --env KNN_FAMILIES="${KNN_FAMILIES:-mean,band,adaptive,ellipse}" \
  --env KNN_SCALINGS="${KNN_SCALINGS:-none,xfit_js}" \
  --env KNN_RAW_DIR="$ROOT/knn_raw" \
  --env KNN_OUT_DIR="$ROOT/knn_out" \
  ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
  "$IMAGE" Rscript $S/knn_config_test.R
rc=$?
set -e
echo "knn_config: finished rc=$rc at $(date -Is)"
exit "$rc"
