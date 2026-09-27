#!/bin/bash -p
# gof_submit.sh -- submit ONE gof_run.sh call as a SLURM job, with validated resources.
#
#   gof_submit.sh <job-name> <dependency> <cpus> <mem-GB> <time> scripts/<name>.R [--NAME=value ...]
#     dependency: none | afterok:<jobid>[:<jobid>...]
#     time:       [D-]HH:MM:SS
# Prints the job id (sbatch --parsable). Jobs whose dependency fails are cancelled
# (--kill-on-invalid-dep=yes) instead of pending forever. The job gets this script's
# own scrubbed environment (env -i below) plus GOF_BIN, never the caller's.
set -euo pipefail
if [[ -z "${GOF_CLEAN:-}" ]]; then
  exec /usr/bin/env -i GOF_CLEAN=1 HOME="${HOME:-$(/usr/bin/getent passwd "$(/usr/bin/id -u)" | /usr/bin/cut -d: -f6)}" USER="${USER:-}" LANG=C.UTF-8 \
    /bin/bash -p "$(/usr/bin/readlink -f "$0")" "$@"
fi
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

die() { echo "gof_submit: $*" >&2; exit 64; }
BIN="$(dirname "$(readlink -f "$0")")"
# shellcheck source=/dev/null
source "$BIN/gof_env.conf"
[[ "$SBATCH_BIN" == /* && -x "$SBATCH_BIN" ]] || die "SBATCH_BIN must be an absolute path to sbatch"

[[ $# -ge 6 ]] || die "usage: gof_submit.sh <job-name> <dependency> <cpus> <mem-GB> <time> scripts/<name>.R [args]"
jname="$1"; dep="$2"; cpus="$3"; mem="$4"; tlim="$5"; shift 5
[[ "$jname" =~ ^[A-Za-z0-9_.-]{1,40}$ ]]                  || die "bad job name '$jname'"
[[ "$dep" == none || "$dep" =~ ^afterok(:[0-9]+)+$ ]]      || die "bad dependency '$dep'"
[[ "$cpus" =~ ^[0-9]+$ ]] && (( cpus >= 1 && cpus <= MAX_CPUS ))     || die "cpus must be 1..$MAX_CPUS"
[[ "$mem" =~ ^[0-9]+$ ]]  && (( mem >= 1 && mem <= MAX_MEM_GB ))     || die "mem-GB must be 1..$MAX_MEM_GB"
[[ "$tlim" =~ ^([0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$ ]]  || die "bad time '$tlim'"
[[ "$1" == scripts/*.R ]]                                   || die "target must be scripts/<name>.R"
for a in "${@:2}"; do
  [[ "$a" == "--dry-run" || "$a" =~ ^--[A-Za-z0-9_]+=[^$'\n']*$ ]] || die "bad argument '$a'"
done

mkdir -p "$LOG_DIR"
opts=(--parsable --job-name "gof_$jname" --cpus-per-task "$cpus" --mem "${mem}G" --time "$tlim"
      --output "$LOG_DIR/%x_%j.out" --error "$LOG_DIR/%x_%j.err" --export "ALL,GOF_BIN=$BIN")
[[ -n "$SLURM_PARTITION" ]] && opts+=(--partition "$SLURM_PARTITION")
[[ -n "$SLURM_ACCOUNT" ]]   && opts+=(--account "$SLURM_ACCOUNT")
[[ "$dep" != none ]]        && opts+=(--dependency "$dep" --kill-on-invalid-dep=yes)
"$SBATCH_BIN" "${opts[@]}" "$BIN/gof_job.sbatch" "$@"
