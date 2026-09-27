# Task: omnibus GOF test (issue #188) on the yeast 48×48 data — autonomous, end to end

**Version 2, 25.09.2026. Replaces `TASK_yeast48_slurm.md`.**

Read `CLAUDE.md` first, then this file completely. Then execute the phases in order,
without asking, until the task is done or a **STOP** condition occurs.

This file was checked against `schroederaaron/empirical-noise-model-testing@a35b966`
(25.09.2026). Where it names a config key, output file or script argument, that
name exists in that commit. If `git pull` brings a newer commit, re-check every name
you use against the code, and STOP if one no longer exists.

---

## A. Rules (apply to every phase)

1. **Git: pull only.** Never `git commit`, `git push`, `gh`, never create branches.
   Changes stay in the working tree. In Phase 8 you write a patch.
2. **Writes only below `/media/BioNAS2/TCGA_TOX_TEST/`.** The sandbox and the guard
   hook enforce this. If something is blocked, do not look for a way around it; STOP.
3. **R and SLURM only through the two wrappers**, always called with the absolute path,
   alone on the command line (no `&&`, `|`, `;`, redirections, `$(...)`, `VAR=...`
   prefixes):
   - `/home/schroder/.claude/bin/gof_run.sh <target> [--NAME=value ...]`
     Runs `Distribution_GOF/<target>` in the R container, synchronously, on this host.
     `<target>` is `scripts/<name>.R` or `tests`.
   - `/home/schroder/.claude/bin/gof_submit.sh <job-name> <dependency> <cpus> <mem-GB> <time> scripts/<name>.R [--NAME=value ...]`
     Submits the same thing as a SLURM job and prints the job id. `<dependency>` is
     `none` or `afterok:<id>[:<id>...]`. `<time>` is `[D-]HH:MM:SS`.
   Never call `docker`, `apptainer`, `sbatch`, `srun`, `salloc`, `scancel` or `sudo`
   directly (they are blocked). You may use `squeue`, `sacct`, `sinfo` and
   `scontrol show` to read.
4. **Every number you report must come from a file or command output you produced in
   this task.** If you could not check something, write "not checked".
5. **STOP** means: write the reason, with the numbers and file paths that show it, into
   the status file (rule 6). End with a short message to Aaron. Do not try a
   workaround, do not change a tolerance, gate, decision value or sample set to get
   past a STOP.
6. **Status file:** `Distribution_GOF/results/yeast48/STATUS.md`. Create it in Phase 0.
   After every phase, and whenever you submit jobs, append: time, phase, what was done,
   job ids, output directories, the next step. **On a resumed session, read it first and
   continue from the recorded next step.**
7. **Waiting for jobs:** poll with `sleep 900` followed by
   `squeue -u schroder -h -o "%i %j %T %M %l %R"` and, once a job has left the queue,
   `sacct -j <id> -X -n -o JobID,JobName%30,State,ExitCode,Elapsed`. If you have been
   polling for 12 hours in one session, update STATUS.md ("waiting for jobs …") and end
   the session; Aaron resumes it later.
8. **A job that ends in any state other than `COMPLETED` with exit code `0:0` is a STOP.**
   Quote the last 40 lines of its `.err` and `.out` log (in
   `Distribution_GOF/results/logs/`, named `gof_<job-name>_<jobid>.{out,err}`).
9. Decisions D1–D7 of the brief keep their `config_gof.R` defaults unless this file
   sets them explicitly. It sets only what is listed below.

### Fixed settings for all runs in this task

| key | value | why |
|---|---|---|
| `ROUNDING` | `error` (default) | the yeast counts are integers; a non-integer is a STOP (Phase 1) |
| `OFFSET_MODE` | `tmm` (default) | htseq counts, no transcript lengths |
| `OFFSETS_BOOT` | `reestimate` (default) | |
| `EXCLUDE_SAMPLES` | not used | the two datasets are separate files with separate labels (see Phase 1). The code does not change the label when samples are excluded, so both would otherwise write into the same output directory. |
| `BASE_SEED` | `188` (default) | |

Paths used below (define them in STATUS.md once):

```
ROOT   = /media/BioNAS2/TCGA_TOX_TEST
REPO   = the empirical-noise-model-testing clone you were started in
GOF    = $REPO/Distribution_GOF
RES    = $GOF/results/yeast48
DATA   = $ROOT/data/yeast48
CLEAN  = $DATA/yeast48_clean.rds
ALL96  = $DATA/yeast48_all.rds
TAG    = round-error_off-tmm        # from gof_mode_tag(); part of every output dir name
```

An `02_run_gof.R` run writes to `<OUT_ROOT>/<label>_<TAG>/`. Every profile below
therefore gets its own `OUT_ROOT`.

---

## Phase 0 — Preflight

1. `git pull` in `$REPO`. Record the HEAD commit in STATUS.md. Check that these exist:
   `$GOF/scripts/{00_validate_families,01_calibration_sim,02_run_gof,03_report,make_dataset}.R`,
   `$GOF/R/io.R`, `$GOF/config_gof.R`. Missing → **STOP**.
2. **Sandbox probe.** Run `touch "$HOME/.gof_sandbox_probe"`. It must **fail** (read-only
   file system / permission denied). If it succeeds, delete the file with
   `rm "$HOME/.gof_sandbox_probe"` and **STOP**: the sandbox is off, and this task may
   only run autonomously with the sandbox on.
3. **Container.** `/home/schroder/.claude/bin/gof_run.sh scripts/02_run_gof.R --dry-run`
   must print the resolved config and exit 0. This also loads every R package
   (`R/setup.R`). If a package has to be built and fails (the README names `nloptr`,
   which needs CMake at build time), **STOP** and quote the error. Aaron installs it.
4. **Unit tests.** `/home/schroder/.claude/bin/gof_run.sh tests` must pass. Failure →
   **STOP**.
5. **Cluster.** Record in STATUS.md:
   - `sinfo -o "%P %a %l %c %m %D"`;
   - `scontrol show partition` (default partition, `MaxTime`);
   - the largest CPU count and memory of a node in the default partition.
   
   Define `MAXTIME` = the default partition's `MaxTime`, `NODE_CPUS`, and `NODE_MEM_GB`.
   `UNLIMITED` counts as 7 days for this task. You may not request more than
   `NODE_CPUS` cores or `NODE_MEM_GB − 8` GB.
6. **Compute-node test.**
   `/home/schroder/.claude/bin/gof_submit.sh preflight none 2 4 00:20:00 scripts/02_run_gof.R --dry-run`.
   Wait for it (rule 7). It must end `COMPLETED 0:0`, and its `.out` log must contain
   the resolved config. This checks that compute nodes have the container runtime, the
   image and the `/media/BioNAS2` mount. Otherwise **STOP**.

---

## Phase 1 — Data

### 1.1 Get it

1. `git clone --depth 1 https://github.com/bartongroup/profDGE48.git $DATA/profDGE48`.
   Record the commit (expected at the time of writing:
   `5ff79149b1eb06629f6dea059c586ebd34086106`; a different one is not a STOP, but
   record it).
2. `sha256sum` of `Preprocessed_data/WT_countdata.tar.gz` and
   `Preprocessed_data/Snf2_countdata.tar.gz`, recorded in STATUS.md.
3. Extract both into `$DATA/counts_raw/`.

### 1.2 Check it (**STOP** on any deviation, with numbers)

Expected (checked by the brief's author on 25.09.2026):

- exactly 96 files `*.gbgout`, 48 starting with `WT_`, 48 with `Snf2_`, named like
  `WT_rep21_MID62_allLanes_tophat2.0.5.bam.gbgout`;
- per file: two tab-separated columns (feature, count). After removing the 5 htseq-count
  summary rows (`no_feature`, `ambiguous`, `too_low_aQual`, `not_aligned`,
  `alignment_not_unique`), 7,126 features remain;
- the feature list is identical, in the same order, in all 96 files;
- all counts are non-negative integers; no feature name contains `ERCC`;
- `$DATA/profDGE48/Bad_replicate_identification/exclude.lst` has 10 lines; each matches
  exactly one count file after appending `.gbgout`. They are 6 WT (reps 21, 22, 25, 28,
  34, 36) and 4 Snf2 (reps 06, 13, 25, 35).

### 1.3 Code (the only code changes in this task)

1. `R/io.R`: add `"htseq_counts"` to `DATASET_SOURCES`. `validate_dataset()` rejects any
   other source string, and none of the existing ones is correct for these data.
2. `R/io.R`: add
   `build_from_htseq_dir(count_dir, exclude_file, label, drop_bad)`. It:
   - reads all `*.gbgout` files and drops the 5 summary rows;
   - checks that the feature lists are identical;
   - builds `counts` (7,126 × n, numeric matrix, colnames `WT_rep21`-style,
     i.e. `<condition>_rep<NN>`);
   - builds `samples` with `sample_id`, `condition` (factor, levels `WT`, `Snf2`),
     `replicate` (integer), `mid`, `bad_replicate` (logical, from `exclude_file`), and
     `lib_size` (column sums);
   - drops the bad replicates if `drop_bad`;
   - returns `new_dataset(counts, samples, design = ~ condition, lengths = NULL,
     source = "htseq_counts", label = label)`.
3. `scripts/make_dataset.R`: add a `yeast48` branch to the `switch(src, ...)`, reading
   `--COUNT_DIR`, `--EXCLUDE_LIST`, `--DROP_BAD` (TRUE/FALSE) and `--LABEL` through the
   existing `arg()` helper. Update the usage comment and the `stop()` message listing the
   valid sources.
4. `tests/testthat/test-htseq.R`: a test that writes 3 small fixture `.gbgout` files and
   an exclude list to `tempdir()`, and checks the dimensions, the dropped summary rows,
   `bad_replicate`, and that mismatching feature lists are an error.
5. Re-run `/home/schroder/.claude/bin/gof_run.sh tests`. Failure → **STOP**.

### 1.4 Build

```
/home/schroder/.claude/bin/gof_run.sh scripts/make_dataset.R --SOURCE=yeast48 --COUNT_DIR=$DATA/counts_raw --EXCLUDE_LIST=$DATA/profDGE48/Bad_replicate_identification/exclude.lst --DROP_BAD=TRUE --LABEL=yeast48_clean --OUT=$CLEAN
/home/schroder/.claude/bin/gof_run.sh scripts/make_dataset.R --SOURCE=yeast48 --COUNT_DIR=$DATA/counts_raw --EXCLUDE_LIST=$DATA/profDGE48/Bad_replicate_identification/exclude.lst --DROP_BAD=FALSE --LABEL=yeast48_all --OUT=$ALL96
```

Write the paths out literally; rule 3 forbids variables on the command line.

**STOP** unless the printed output shows:
- `yeast48_clean`: 7126 genes × 86 samples (42 WT + 44 Snf2);
- `yeast48_all`: 7126 × 96;
- design `~condition`;
- integer audit: 0 non-integer entries.

### 1.5 Document

Write `$DATA/README.md`. It contains:
- source URL and commit;
- checksums;
- provenance per Gierliński et al. 2015 (*Bioinformatics* 31:3625, §2.1): TopHat2
  v2.0.5, htseq-count v0.5.3p9, Ensembl v68;
- the exclusion list, and the fact that it is the authors' own replicate QC, fixed before
  any GOF result;
- the two dataset files with their dimensions.

---

## Phase 2 — Validation gates G1–G4 (SLURM)

Submit:

```
/home/schroder/.claude/bin/gof_submit.sh val_g1g4 none <cpus> <mem> <time> scripts/00_validate_families.R --DATASET=<CLEAN> --OUT_ROOT=<RES>/validation
```

- **Resources:** `cpus = min(16, NODE_CPUS)`, `mem = min(32, NODE_MEM_GB − 8)`,
  `time = min(1-00:00:00, MAXTIME)`.
- **Why `--DATASET`:** G4 then uses the clean dataset's n (86).
- **Outputs:** `00_validate_families.R` exits 1 if G1, G2 or G3 fails, so rule 8 covers
  the gate. Its outputs go to `<RES>/validation/validation/`, because the script appends
  `validation`.
- **After it completes:** copy the `verdict.csv` rows and the G4 bias table into
  STATUS.md verbatim.

---

## Phase 3 — Development runs (SLURM, both submitted together, both `afterok` on Phase 2)

Only after Phase 2 completed. Use `none` as the dependency if you submit after it
finished.

**Resources for both:** `cpus = min(32, NODE_CPUS)`, `mem = min(64, NODE_MEM_GB − 8)`,
`time = min(2-00:00:00, MAXTIME)`.

**dev-a — whole pipeline, small:**

```
/home/schroder/.claude/bin/gof_submit.sh dev_a none <cpus> <mem> <time> scripts/02_run_gof.R --DATASET=<CLEAN> --OUT_ROOT=<RES>/dev_a --MAX_GENES_OBS=500 --N_GENES_BOOT=500 --B=100
/home/schroder/.claude/bin/gof_submit.sh dev_a_report afterok:<dev_a id> 2 8 02:00:00 scripts/03_report.R --RUN_DIR=<RES>/dev_a/yeast48_clean_round-error_off-tmm
```

With `N_GENES_BOOT` ≥ the pool, `select_boot_genes()` returns all genes where every
model fitted, so G_BOOT is at most 500 genes.

**dev-b — all genes, no bootstrap, no held-out score.** This gives the fit failures and
fit times across the whole expression range, and it is the G5 reference:

```
/home/schroder/.claude/bin/gof_submit.sh dev_b none <cpus> <mem> <time> scripts/02_run_gof.R --DATASET=<CLEAN> --OUT_ROOT=<RES>/dev_b --MAX_GENES_OBS=0 --RUN_BOOTSTRAP=FALSE --RUN_HELDOUT=FALSE
```

`N_GENES_BOOT` stays at its default (2000). Step 5 of `02_run_gof.R` still runs and
writes `gboot_data.rds` and `fits_gall_<model>.rds`, which `01_calibration_sim.R` reads.

**Record in STATUS.md** (copied from the files, not paraphrased):

| run | what to record |
|---|---|
| both | `sacct` `Elapsed` |
| both | `preprocess.csv`: genes before/after `filterByExpr` |
| both | `analysis_genes.txt` line count |
| both | `fit_failures.csv` (all rows) |
| both | `g6_pit_asserts.csv`: whether every row passed |
| dev-a | `boot_timing.csv` (all rows) |
| dev-a | `summary_table.md` (as is) |
| dev-b | median of `time_s` per model from `fit_diag_gall_<m>.csv` |

**STOP** if:
- the analysis set of dev-b is smaller than 50% of the filtered genes. Report
  `fit_failures.csv`: which model fails where. Deciding what to do about a failing
  model is Aaron's call.
- `dev_a_report` produced no `summary_table.md`.

---

## Phase 4 — Projection of the full runs (no job; arithmetic only, shown in STATUS.md)

Write every input number and every step of the arithmetic into STATUS.md. These are
**rough linear extrapolations**; label them as such.

1. **Full `02_run_gof.R` run** (`N_GENES_BOOT = 2000`, `B = 500`, K = 8 held-out folds):

   `T_full ≈ Elapsed(dev_b) + Elapsed(dev_a) × (2000 / n_GBOOT_dev_a) × (500 / 100)`

   - `n_GBOOT_dev_a` is the line count of `<RES>/dev_a/yeast48_clean_round-error_off-tmm/gboot_genes.txt`.
   - The second term over-counts the parts of dev-a that don't scale with B. That is
     acceptable because it errs on the long side.
2. **G5** (`01_calibration_sim.R`, defaults M = 40, G = 300, B = 99; truths NB and PLN,
   each fitted as the true model and as Poisson). From dev-a's `boot_timing.csv`, per
   model m, take `c_m = median_s_per_replicate / n_genes`, i.e. seconds per gene per
   replicate, serial.

   `T_G5 ≈ 40 × 300 × 100 × (c_nb + c_pln + 2 × c_poisson) / cpus`

   The factor 100 is 99 replicates plus the observed fit.
3. Multiply both by **1.5** (margin). Request `time = ceil(1.5 × T)` for each job.
4. **STOP** if either `1.5 × T` exceeds `MAXTIME`. Neither script can resume a
   half-finished bootstrap, so the job would be killed and the work lost. Report both
   projections. Changing `N_GENES_BOOT`, `B` or G5's M/G/B is decision D4/G5 design,
   i.e. Aaron's call.

---

## Phase 5 — Gate G5 and full runs (one SLURM chain, submitted together)

**Resources:** `cpus = min(32, NODE_CPUS)`, `mem = min(64, NODE_MEM_GB − 8)`, `time`
from Phase 4.

```
G5          : gof_submit.sh g5 none <cpus> <mem> <time_G5> scripts/01_calibration_sim.R --G5_REF_DIR=<RES>/dev_b/yeast48_clean_round-error_off-tmm --OUT_ROOT=<RES>/g5
full_clean  : gof_submit.sh full_clean afterok:<G5> <cpus> <mem> <time_full> scripts/02_run_gof.R --DATASET=<CLEAN> --OUT_ROOT=<RES>/full --B=500
rep_clean   : gof_submit.sh rep_clean afterok:<full_clean> 2 8 02:00:00 scripts/03_report.R --RUN_DIR=<RES>/full/yeast48_clean_round-error_off-tmm
full_all96  : gof_submit.sh full_all96 afterok:<G5> <cpus> <mem> <time_full> scripts/02_run_gof.R --DATASET=<ALL96> --OUT_ROOT=<RES>/full --B=500
rep_all96   : gof_submit.sh rep_all96 afterok:<full_all96> 2 8 02:00:00 scripts/03_report.R --RUN_DIR=<RES>/full/yeast48_all_round-error_off-tmm
```

(Each line is a full `/home/schroder/.claude/bin/gof_submit.sh …` call, with the literal
job ids and paths.)

- **How the gate works:** `01_calibration_sim.R` exits 1 when the full-size G5 gate
  fails. Dependent jobs are then cancelled automatically (`--kill-on-invalid-dep=yes` in
  the wrapper).
- **G5 is only a gate at full size.** With M = 40, G = 300 and B = 99 unchanged, the
  output goes to `<RES>/g5/G5_calibration/`. A `G5_calibration_SMOKE` directory must not
  appear.
- **Which dataset G5 uses:** the reference directory is dev-b on the **clean** data. G5
  simulates from the NB and PLN fits there.

Record all job ids in STATUS.md. Then wait (rule 7).

**After G5:** copy `G5_summary.csv` into STATUS.md. `State` ≠ `COMPLETED 0:0` →
**STOP** (rule 8). If the gate failed, say so explicitly, and report which rows have
`pass = FALSE`.

---

## Phase 6 — Salmon feasibility check (no download; do it while Phase 5 jobs run)

The Salmon arm is not run in this task. Collect the facts Aaron needs to decide.

1. Is `salmon` available on the host? (`command -v salmon`; inside the container it is
   not part of the image.)
2. ENA file report, via WebFetch or `curl` (`www.ebi.ac.uk` is allowed):
   `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB5348&result=read_run&fields=run_accession,sample_alias,read_count,base_count,fastq_bytes,fastq_ftp&format=tsv`
   Save it as `$DATA/ena_filereport_PRJEB5348.tsv`. Report the number of runs, the total
   read count, and the total `fastq_bytes` (sum over all files, in GB).
3. Do **not** download FASTQ files. Put the findings into STATUS.md under "Salmon arm —
   facts for a decision".

---

## Phase 7 — Results (after `rep_clean` and `rep_all96` completed)

Write `$RES/RUN_REPORT.md`. It contains **only** numbers and tables copied from output
files, each with its source path:

1. Data: the two datasets' dimensions and audit (Phase 1).
2. Gates: G1–G4 verdict, G4 bias table, G5 summary.
3. For `yeast48_clean` and `yeast48_all`, each:
   - `summary_table.md` as is;
   - the per-model row counts of `fit_failures.csv`;
   - `boot_strata_summary.csv` restricted to the mean-expression strata;
   - `heldout_summary.csv`;
   - the list of plot files in `plots/`.
4. The known limits, copied verbatim from `TABLE_CAPTION` in `R/report.R`.

**No interpretation.** Do not write which distribution "fits", is "best" or is
"rejected". Aaron and Asis interpret the results.

---

## Phase 8 — Hand-over

1. `git status --porcelain` in `$REPO`.
2. `git diff --no-color > Distribution_GOF/results/yeast48/yeast48_code.patch` for the
   tracked files changed in Phase 1.3.
3. List the new untracked files (at least `tests/testthat/test-htseq.R`) with their full
   paths in STATUS.md. Do not stage anything.
4. Final STATUS.md entry: "DONE", with pointers to `RUN_REPORT.md`, the patch and the
   new files.
