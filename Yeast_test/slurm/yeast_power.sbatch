#!/bin/bash
# yeast_power.sbatch -- power benchmark on the yeast 48x48 data (Schurch et al. 2016):
# Parts RS (2n samples per round, n = 3, 5, 7, 10, 50 rounds) and R (whole-condition
# halves, 20 rounds) on the WT and the Snf2 replicates, then power_report.R.
# The engine is shared with TCGA (analysis/TCGA_test/scripts/power_test.R); the yeast
# cohorts come from Yeast_test/scripts/make_yeast_cohorts.R, built here once from the
# Distribution_GOF dataset when they do not exist yet.
# Writes ONLY to $ROOT/yeast_power_raw and $ROOT/yeast_power_out.
#
# Submit from the directory that holds logs/:   mkdir -p logs && sbatch yeast_power.sbatch
# Options (sbatch --export=ALL,...):
#   YEAST_DATASET   default $DATA/yeast48_clean.rds (bad replicates removed); yeast48_all.rds = all 96
#   YEAST_GTF       yeast GTF for gene lengths (otherwise TPM = CPM, lengths = 1)
#   PARTS           default "RS,R"
#   POWER_ROUNDS_RS default 50
#   REBUILD=1       rebuild the cohorts even if they exist (e.g. after changing the dataset or GTF)
#
# TUNE (guesses): --time / --mem.
#SBATCH --job-name=yeast_power
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
PARTS=${PARTS:-RS,R}

[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
for f in "$DATASET" "$WORK/$S/power_test.R" "$WORK/$S/power_report.R" "$WORK/$Y/make_yeast_cohorts.R"; do
  [[ -f "$f" ]] || { echo "yeast_power: missing $f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"
echo "yeast_power: $(date -Is) host=$(hostname) cpus=$CPUS parts=$PARTS dataset=$DATASET commit=${GIT_COMMIT:-unknown}"

cd "$WORK"
mkdir -p .container-home
run() { apptainer exec --cleanenv --home "$WORK/.container-home" -B /media/BioNAS2 \
          --env OMP_NUM_THREADS=1 "$@"; }

if [[ "${REBUILD:-0}" == 1 || ! -f "$COH/yeast48_wt.rds" || ! -f "$COH/yeast48_snf2.rds" ]]; then
  run "$IMAGE" Rscript $Y/make_yeast_cohorts.R "--DATASET=$DATASET" "--OUT_DIR=$COH" ${YEAST_GTF:+"--GTF=$YEAST_GTF"}
fi

set +e
run --env POWER_PARTS="$PARTS" \
    --env POWER_ROUNDS_RS="${POWER_ROUNDS_RS:-50}" \
    --env POWER_CORES="$CPUS" \
    --env POWER_REAL_RDS="$COH/yeast48_wt.rds,$COH/yeast48_snf2.rds" \
    --env POWER_RAW_DIR="$ROOT/yeast_power_raw" \
    --env POWER_OUT_DIR="$ROOT/yeast_power_out" \
    ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
    "$IMAGE" bash -c "Rscript $S/power_test.R && Rscript $S/power_report.R"
rc=$?
set -e
echo "yeast_power: finished rc=$rc at $(date -Is)"
exit "$rc"
