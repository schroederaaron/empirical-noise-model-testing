#!/bin/bash
# power_rs.sbatch -- Part RS (note D WP3): real TCGA noise at n = 3 and n = 5 per group,
# null arm (pi1 = 0) and power arm (pi1 = 0.1), then power_report.R.
# Copied from power_full.sbatch; PARTS = RS, no dispersion file (P+ only).
# Two log files: tables and titles in logs/%x_%j.out (stdout), messages in .err (stderr).
#
# Submit from the directory that holds logs/ (the job itself cd-s into $WORK):   mkdir -p logs && sbatch power_rs.sbatch
#
# SMOKE TEST first (note D section 4 step 2): 2 rounds, two cohorts (one healthy, one
# Stage IV), separate raw/out folders, runtime and storage extrapolated to the full
# 2400 datasets on 32 cores (printed at the end of power_test.R):
#   sbatch --export=ALL,POWER_ROUNDS_RS=2,RS_COHORTS="kirc_healthy kirc_stageiv",RS_TAG=_rs_smoke power_rs.sbatch
# FULL RUN only after the smoke run looked sane and its estimates were checked:
#   sbatch power_rs.sbatch
#
# What overwrites what:
#   * power_test.R clears the raw directory of every part it runs (here: RS) and nothing else.
#   * power_report.R rewrites $ROOT/power_out$RS_TAG from ALL raw parts present in
#     $ROOT/power_raw$RS_TAG. With the default (empty) tag that includes Part R, which the
#     sample-size curve needs. Copy $ROOT/power_out elsewhere first to keep the current report.
#
# TUNE (guesses, not measured): --mem and --time. Size them from the smoke run's
# extrapolation and sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,AllocCPUS,State
#SBATCH --job-name=power_rs
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
ROUNDS=${POWER_ROUNDS_RS:-100}
TAG=${RS_TAG:-}                  # "" = default power_raw / power_out
COHORTS=${RS_COHORTS:-"coad_healthy coad_stageiv luad_healthy luad_stageiv kirc_healthy kirc_stageiv"}

REAL_RDS=""
for c in $COHORTS; do
  [[ -f "$WORK/power_cohort/$c.rds" ]] || { echo "power_rs: missing $WORK/power_cohort/$c.rds" >&2; exit 2; }
  REAL_RDS+="${REAL_RDS:+,}power_cohort/$c.rds"
done
[[ -f "$IMAGE" ]] || { echo "missing image $IMAGE" >&2; exit 2; }
for f in $S/power_test.R $S/power_report.R; do
  [[ -f "$WORK/$f" ]] || { echo "power_rs: missing $WORK/$f" >&2; exit 2; }
done
GIT_COMMIT="$(git -C "$WORK/analysis" rev-parse HEAD 2>/dev/null || true)"

echo "power_rs: $(date -Is) host=$(hostname) cpus=$CPUS mem=${SLURM_MEM_PER_NODE:-?}MB rounds=$ROUNDS cohorts=[$COHORTS] tag=${TAG:-none} commit=${GIT_COMMIT:-unknown}"

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
  --env POWER_PARTS=RS \
  --env POWER_ROUNDS_RS="$ROUNDS" \
  --env POWER_CORES="$CPUS" \
  --env POWER_EXTRAP_JOBS=2400 \
  --env POWER_EXTRAP_CORES=32 \
  --env POWER_REAL_RDS="$REAL_RDS" \
  --env POWER_RAW_DIR="$ROOT/power_raw$TAG" \
  --env POWER_OUT_DIR="$ROOT/power_out$TAG" \
  ${GIT_COMMIT:+--env POWER_GIT_COMMIT="$GIT_COMMIT"} \
  "$IMAGE" bash -c "
    Rscript $S/power_test.R &&
    Rscript $S/power_report.R"
rc=$?
set -e
echo "power_rs: finished rc=$rc at $(date -Is)"
exit "$rc"
