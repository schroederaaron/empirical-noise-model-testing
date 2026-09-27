#!/bin/bash -p
# gof_run.sh -- run ONE Distribution_GOF R script (or the unit tests) in the container.
#
#   gof_run.sh scripts/<name>.R [--NAME=value ...] [--dry-run]
#   gof_run.sh tests
#
# Safety properties (why this file lives in ~/.claude/bin, outside TCGA_TOX_TEST):
#   * the container gets exactly ONE host mount: ALLOWED_ROOT (from gof_env.conf);
#   * it runs as your uid, never as root; no --privileged, no extra mounts;
#   * only scripts under $GOF_REPO/Distribution_GOF/scripts/*.R can be run;
#   * the environment is rebuilt from scratch (env -i), PATH is fixed, and bash runs
#     with -p, so BASH_ENV / exported functions from the caller are ignored;
#   * inside SLURM, cores/memory of the container follow the allocation, and
#     SIGTERM (scancel, time limit) kills the container.
set -euo pipefail

if [[ -z "${GOF_CLEAN:-}" ]]; then
  pass=(GOF_CLEAN=1 HOME="${HOME:-$(/usr/bin/getent passwd "$(/usr/bin/id -u)" | /usr/bin/cut -d: -f6)}" USER="${USER:-}" LANG=C.UTF-8)
  for v in SLURM_JOB_ID SLURM_CPUS_PER_TASK SLURM_MEM_PER_NODE SLURM_JOB_NAME; do
    [[ -n "${!v:-}" ]] && pass+=("$v=${!v}")
  done
  exec /usr/bin/env -i "${pass[@]}" /bin/bash -p "$(/usr/bin/readlink -f "$0")" "$@"
fi
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

die() { echo "gof_run: $*" >&2; exit 64; }
BIN="$(dirname "$(readlink -f "$0")")"
CONF="$BIN/gof_env.conf"
[[ -f "$CONF" ]] || die "missing $CONF"
# shellcheck source=/dev/null
source "$CONF"

inside() { local rp; rp="$(realpath -m -- "$1")"; [[ "$rp" == "$ROOT_RP" || "$rp" == "$ROOT_RP"/* ]]; }
ROOT_RP="$(realpath -m -- "$ALLOWED_ROOT")"
for p in "$GOF_REPO" "$WORK_DIR" "$LOG_DIR"; do inside "$p" || die "$p is outside $ALLOWED_ROOT (fix gof_env.conf)"; done
[[ "$RUNTIME_BIN" == /* && -x "$RUNTIME_BIN" ]] || die "RUNTIME_BIN must be an absolute path to an executable"
GOF_DIR="$(realpath -m -- "$GOF_REPO/Distribution_GOF")"

[[ $# -ge 1 ]] || die "usage: gof_run.sh scripts/<name>.R [--NAME=value ...] | tests"
target="$1"; shift
for a in "$@"; do
  [[ "$a" == "--dry-run" || "$a" =~ ^--[A-Za-z0-9_]+=[^$'\n']*$ ]] || die "bad argument '$a' (allowed: --NAME=value, --dry-run)"
done

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  cpus="${SLURM_CPUS_PER_TASK:-1}"
  mem_arg=(); [[ -n "${SLURM_MEM_PER_NODE:-}" ]] && mem_arg=(--memory "${SLURM_MEM_PER_NODE}m")
else
  cpus="$LOCAL_CPUS"; mem_arg=()
fi
[[ "$cpus" =~ ^[0-9]+$ ]] || die "bad cpu count '$cpus'"

if [[ "$target" == "tests" ]]; then
  cmd=(Rscript -e 'GOF_ROOT <- Sys.getenv("GOF_ROOT"); source(file.path(GOF_ROOT, "R", "setup.R")); load_or_install("testthat"); testthat::test_dir(file.path(GOF_ROOT, "tests", "testthat"), stop_on_failure = TRUE)')
else
  script="$(realpath -m -- "$GOF_DIR/$target")"
  [[ "$script" == "$GOF_DIR/scripts/"*.R && -f "$script" ]] || die "'$target' is not an existing Distribution_GOF/scripts/*.R file"
  cmd=(Rscript "$script" "$@" "--N_CORES=$cpus")      # last --N_CORES wins (config_gof.R)
fi

name="gof_${SLURM_JOB_ID:-local}_$$"
mkdir -p "$LOG_DIR"
echo "gof_run: $(date -Is) host=$(hostname) runtime=$RUNTIME image=$IMAGE cpus=$cpus target=$target" >&2

case "$RUNTIME" in
  docker)
    cleanup() { "$RUNTIME_BIN" kill "$name" >/dev/null 2>&1 || true; }
    trap cleanup TERM INT
    "$RUNTIME_BIN" run --rm --init --name "$name" \
      --user "$(id -u):$(id -g)" --cpus "$cpus" "${mem_arg[@]}" \
      -e HOME=/tmp -e OMP_NUM_THREADS=1 -e GOF_ROOT="$GOF_DIR" \
      -v "$ROOT_RP:$ROOT_RP" -w "$WORK_DIR" \
      "$IMAGE" "${cmd[@]}" &
    ;;
  apptainer)
    APPTAINERENV_HOME=/tmp APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_GOF_ROOT="$GOF_DIR" \
    "$RUNTIME_BIN" exec --containall --cleanenv --bind "$ROOT_RP" --pwd "$WORK_DIR" \
      "$IMAGE" "${cmd[@]}" &
    ;;
  *) die "RUNTIME must be docker or apptainer" ;;
esac
pid=$!
set +e; wait "$pid"; rc=$?; set -e
exit "$rc"
