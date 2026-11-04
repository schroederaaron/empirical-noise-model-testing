#!/bin/bash
# yeast_knn.sbatch -- neighbourhood configuration grid (TCGA_test/scripts/knn_config_test.R)
# on the yeast 48x48 data: Part RS only (2n samples per round from the WT and the Snf2
# replicates), n = 3, 5, 10, 20 rounds, the four neighbourhoods x none / cross-fit JS,
# then power_report.R and the kNN summary. Simulated data (P+) is not data-specific
# and stays in TCGA_test/slurm/knn_config.sbatch.
# Writes ONLY to $ROOT/yeast_knn_raw and $ROOT/yeast_knn_out.
#
# Submit from the directory that holds logs/:   mkdir -p logs && sbatch yeast_knn.sbatch
# Options (sbatch --export=ALL,...): YEAST_DATASET, YEAST_GTF, REBUILD=1 (as yeast_power.sbatch),
#   KNN_FAMILIES (default mean,band,adaptive,ellipse), KNN_N (3,5,10), KNN_ROUNDS (20)
#
# TUNE (guesses): --time / --mem.
#SBATCH --job-name=yeast_knn
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
Y=analysis/Yeast_test/scripts
DATA=$ROOT/data/yeast48
DATASET=${YEAST_DATASET:-$DATA/yeast48_clean.rds}
COH=$DATA/cohorts_$(basename "$DATASET" .rds)
CPUS=${SLURM_CPUS_PER_TASK:-32}

[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
for f in "$DATASET" "$WORK/$S/knn_config_test.R" "$WORK/$S/power_test.R" "$WORK/$S/power_report.R" \
         "$WORK/$Y/make_yeast_cohorts.R"; do
  [[ -f "$f" ]] || { echo "yeast_knn: missing $f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"
echo "yeast_knn: $(date -Is) host=$(hostname) cpus=$CPUS dataset=$DATASET families=${KNN_FAMILIES:-all} commit=${GIT_COMMIT:-unknown}"

cd "$WORK"
mkdir -p .container-home
run() { apptainer exec --cleanenv --home "$WORK/.container-home" -B /media/BioNAS2 \
          --env OMP_NUM_THREADS=1 "$@"; }

if [[ "${REBUILD:-0}" == 1 || ! -f "$COH/yeast48_wt.rds" || ! -f "$COH/yeast48_snf2.rds" ]]; then
  run "$IMAGE" Rscript $Y/make_yeast_cohorts.R "--DATASET=$DATASET" "--OUT_DIR=$COH" ${YEAST_GTF:+"--GTF=$YEAST_GTF"}
fi

set +e
run --env POWER_CORES="$CPUS" \
    --env POWER_REAL_RDS="$COH/yeast48_wt.rds,$COH/yeast48_snf2.rds" \
    --env KNN_PARTS=RS \
    --env KNN_ROUNDS="${KNN_ROUNDS:-20}" \
    --env KNN_N="${KNN_N:-3,5,10}" \
    --env KNN_FAMILIES="${KNN_FAMILIES:-mean,band,adaptive,ellipse}" \
    --env KNN_SCALINGS="${KNN_SCALINGS:-none,xfit_js}" \
    --env KNN_RAW_DIR="$ROOT/yeast_knn_raw" \
    --env KNN_OUT_DIR="$ROOT/yeast_knn_out" \
    ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
    "$IMAGE" Rscript $S/knn_config_test.R
rc=$?
set -e
echo "yeast_knn: finished rc=$rc at $(date -Is)"
exit "$rc"
