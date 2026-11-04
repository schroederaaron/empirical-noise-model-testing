#!/bin/bash
# gof_yeast48_submit.sh -- submit the omnibus GOF analysis (issue #188) of the yeast
# 48x48 data as a chain of SLURM jobs. Run it BY HAND from the Tensor-Omics root.
#
#   ./gof_yeast48_submit.sh [--dry-run] <stage> [--after=<jobid>]
#
# stages (each stage's jobs start only if the previous job exited 0: afterok,
# and a failed dependency cancels the rest: --kill-on-invalid-dep=yes):
#   setup  one-time package install (gof_setup_slurm.sh; needs a user-owned
#          external/docker_r_libs and cmake in the image)
#          -> unit tests -> rebuild both datasets in the container (R 4.6.1)
#   dev    G1-G4 gates (00_validate_families.R, exits 1 on a G1-G3 failure)
#          -> dev_a: whole pipeline on 500 genes, B = 100, then 03_report.R
#          -> dev_b: all genes, no bootstrap / held-out (fit failures, fit times, G5 reference)
#   full   G5 gate (01_calibration_sim.R, reference = dev_b)
#          -> full run on yeast48_clean and on yeast48_all (B = 500), each then 03_report.R
#   all    setup -> dev -> full in one chain
#
# --after=<jobid> makes the first job of the stage wait for that job (afterok).
# "full" needs dev_b's output, so submit it after dev completed or with --after=<dev_b id>.
# --dry-run prints the sbatch commands without submitting (job ids shown as <name>).
#
# Resources below are PLACEHOLDERS until dev_a/dev_b timings exist (boot_timing.csv,
# fit_diag_gall_<model>.csv). LocalQ has no time limit, so no --time is requested.
set -euo pipefail

CPUS_BIG=32;  MEM_BIG=64G       # dev_a, dev_b, G5, full runs
CPUS_VAL=16;  MEM_VAL=32G       # 00_validate_families.R
CPUS_SMALL=4; MEM_SMALL=16G     # tests, dataset build, reports

ROOT="$PWD"
GOF="$ROOT/analysis/Distribution_GOF"
RES="$GOF/results/yeast48"
LOGS="$GOF/results/logs"
DATA=/media/BioNAS2/TCGA_TOX_TEST/data/yeast48
CLEAN="$DATA/yeast48_clean.rds"
ALL96="$DATA/yeast48_all.rds"
TAG=round-error_off-tmm              # gof_mode_tag() under the default ROUNDING / OFFSET_MODE

DRY=0; STAGE=""; AFTER=""
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=1 ;;
    --after=*) AFTER="${a#--after=}" ;;
    setup|dev|full|all) STAGE="$a" ;;
    *) echo "unknown argument '$a'" >&2; exit 64 ;;
  esac
done
[[ -n "$STAGE" ]] || { sed -n '2,24p' "$0"; exit 64; }
[[ -f "$GOF/config_gof.R" && -f "$ROOT/gof_job_slurm.sh" ]] || { echo "run from the Tensor-Omics root" >&2; exit 64; }
[[ $DRY -eq 1 ]] || mkdir -p "$LOGS"

# submit <name> <dependency-jobid|""> <cpus> <mem> <script> [args...]  -> sets JID
submit() {
  local name="$1" dep="$2" cpus="$3" mem="$4"; shift 4
  local opts=(--parsable --job-name "gof_$name" --cpus-per-task "$cpus" --mem "$mem"
              --output "$LOGS/%x_%j.out" --error "$LOGS/%x_%j.err")
  [[ -n "$dep" ]] && opts+=(--dependency "afterok:$dep" --kill-on-invalid-dep=yes)
  if [[ $DRY -eq 1 ]]; then
    echo "sbatch ${opts[*]} $*"; JID="<$name>"
  else
    JID=$(sbatch "${opts[@]}" "$@")
    echo "submitted gof_$name: job $JID${dep:+ (after $dep)}"
  fi
}
job() { local n="$1" d="$2" c="$3" m="$4"; shift 4; submit "$n" "$d" "$c" "$m" "$ROOT/gof_job_slurm.sh" "$@"; }

stage_setup() {
  submit setup "$AFTER" 8 32G "$ROOT/gof_setup_slurm.sh"; local s=$JID
  job tests "$s" $CPUS_SMALL $MEM_SMALL tests; local t=$JID
  local common=(--SOURCE=yeast48 "--COUNT_DIR=$DATA/counts_raw"
                "--EXCLUDE_LIST=$DATA/profDGE48/Bad_replicate_identification/exclude.lst")
  job data_clean "$t" $CPUS_SMALL $MEM_SMALL scripts/make_dataset.R "${common[@]}" --DROP_BAD=TRUE  --LABEL=yeast48_clean "--OUT=$CLEAN"
  local d1=$JID
  job data_all "$d1" $CPUS_SMALL $MEM_SMALL scripts/make_dataset.R "${common[@]}" --DROP_BAD=FALSE --LABEL=yeast48_all "--OUT=$ALL96"
  AFTER=$JID
}

stage_dev() {
  job val_g1g4 "$AFTER" $CPUS_VAL $MEM_VAL scripts/00_validate_families.R "--DATASET=$CLEAN" "--OUT_ROOT=$RES/validation"
  local v=$JID
  job dev_a "$v" $CPUS_BIG $MEM_BIG scripts/02_run_gof.R "--DATASET=$CLEAN" "--OUT_ROOT=$RES/dev_a" \
      --MAX_GENES_OBS=500 --N_GENES_BOOT=500 --B=100
  job dev_a_report "$JID" $CPUS_SMALL $MEM_SMALL scripts/03_report.R "--RUN_DIR=$RES/dev_a/yeast48_clean_$TAG"
  job dev_b "$v" $CPUS_BIG $MEM_BIG scripts/02_run_gof.R "--DATASET=$CLEAN" "--OUT_ROOT=$RES/dev_b" \
      --MAX_GENES_OBS=0 --RUN_BOOTSTRAP=FALSE --RUN_HELDOUT=FALSE
  AFTER=$JID
}

stage_full() {
  if [[ -z "$AFTER" && $DRY -eq 0 && ! -f "$RES/dev_b/yeast48_clean_$TAG/gboot_data.rds" ]]; then
    echo "no dev_b output in $RES/dev_b/yeast48_clean_$TAG -- run 'dev' first or pass --after=<dev_b id>" >&2; exit 64
  fi
  job g5 "$AFTER" $CPUS_BIG $MEM_BIG scripts/01_calibration_sim.R "--G5_REF_DIR=$RES/dev_b/yeast48_clean_$TAG" "--OUT_ROOT=$RES/g5"
  local g=$JID
  job full_clean "$g" $CPUS_BIG $MEM_BIG scripts/02_run_gof.R "--DATASET=$CLEAN" "--OUT_ROOT=$RES/full" --B=500
  job rep_clean "$JID" $CPUS_SMALL $MEM_SMALL scripts/03_report.R "--RUN_DIR=$RES/full/yeast48_clean_$TAG"
  job full_all96 "$g" $CPUS_BIG $MEM_BIG scripts/02_run_gof.R "--DATASET=$ALL96" "--OUT_ROOT=$RES/full" --B=500
  job rep_all96 "$JID" $CPUS_SMALL $MEM_SMALL scripts/03_report.R "--RUN_DIR=$RES/full/yeast48_all_$TAG"
}

case "$STAGE" in
  setup) stage_setup ;;
  dev)   stage_dev ;;
  full)  stage_full ;;
  all)   stage_setup; stage_dev; stage_full ;;
esac
