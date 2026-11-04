#!/bin/bash
# gof_setup_slurm.sh -- ONE-TIME package install for Distribution_GOF (issue #188).
#
# Runs as YOUR uid, like every GOF job. Prerequisites (done once, by hand):
#   * external/docker_r_libs is owned by you (not root), and has no stale 00LOCK-* dirs;
#   * the arch-gfortran image contains cmake: glmmTMB -> lme4 -> nloptr needs it at
#     BUILD time (see the Distribution_GOF README, "Packages").
#
#   sbatch gof_setup_slurm.sh        (from the Tensor-Omics root)
#SBATCH --job-name=gof_setup
#SBATCH --output=analysis/Distribution_GOF/results/logs/%x_%j.out
#SBATCH --error=analysis/Distribution_GOF/results/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
set -euo pipefail

LIB="$PWD/external/docker_r_libs"
[[ -w "$LIB" ]] || { echo "gof_setup: $LIB is not writable by $(id -un) -- change its owner first" >&2; exit 64; }
if compgen -G "$LIB/00LOCK-*" >/dev/null; then
  echo "gof_setup: stale lock dirs in $LIB (00LOCK-*) block reinstalls -- remove them first" >&2; exit 64
fi

docker run --rm \
  --user "$(id -u):$(id -g)" -e HOME=/tmp -e MAKEFLAGS=-j8 \
  -e GOF_ROOT="$PWD/analysis/Distribution_GOF" \
  -v /media/BioNAS2:/media/BioNAS2 -w "$PWD" \
  arch-gfortran bash -c '
set -e
command -v cmake >/dev/null || { echo "gof_setup: cmake not in the image (needed to build nloptr)" >&2; exit 64; }
Rscript -e "
  GOF_ROOT <- Sys.getenv(\"GOF_ROOT\")
  source(file.path(GOF_ROOT, \"R\", \"setup.R\"))      # installs glmmTMB, edgeR, data.table, statmod, ...
  for (p in c(\"testthat\", \"goftest\", \"ggplot2\", \"tximport\")) load_or_install(p)
  for (p in c(\"glmmTMB\", \"TMB\", \"edgeR\", \"tximport\", \"testthat\", \"goftest\", \"ggplot2\"))
    cat(sprintf(\"%-9s %s\n\", p, as.character(packageVersion(p))))
"'
