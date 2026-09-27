#!/bin/bash -p
# gof_job.sbatch -- fixed SLURM job body. All #SBATCH options come from gof_submit.sh
# (command-line options override directives; there are none here on purpose).
# GOF_BIN is set by gof_submit.sh from its own location (--export ALL,GOF_BIN=... of its scrubbed env).
set -euo pipefail
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
[[ -n "${GOF_BIN:-}" && -x "$GOF_BIN/gof_run.sh" ]] || { echo "gof_job: GOF_BIN not set" >&2; exit 64; }
echo "gof_job: $(date -Is) job=${SLURM_JOB_ID:-?} node=$(hostname) args: $*" >&2
exec "$GOF_BIN/gof_run.sh" "$@"
