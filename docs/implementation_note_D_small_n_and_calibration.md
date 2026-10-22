# Implementation note D — real-data calibration and power at n = 3 and n = 5, wider band grid, readable output, parameter guard

**For:** Claude Code. **Builds on:** notes A and B (variance-matched neighbourhoods; already implemented as `engine = "variance"` with `variance_mode = "band" | "shrunk"`), and **note C (console output): implement note C in full as work package 1 of this note.**
**Repos:** `empirical-noise-model-testing` (`main`, commit 36b00df when this was written) and `asishallab-group/Tensor-Omics` (branch `145-empirical-noise-model`).
**Status:** written without executing any code. Everything marked **VERIFY** must be checked in the code, not assumed.

---

## 1. Why this note exists (findings so far, in short)

* **Part R (full-cohort halves, n = 13-41 per group, 6 TCGA base cohorts, 20 rounds):** the variance-matched arms are a large improvement over `TOX-log|TPM-centred`. With the wide band (`c = 2`):
  AUPRC excess +0.052 to +0.109 (6/6 cohorts), 82-94% of the gap to limma closed, TPR at nominal 0.05 from 0.14-0.53 to 0.36-0.72, observed FDR at nominal lower than TOX-log in 6/6 cohorts (0.007-0.049 vs 0.054-0.071).
  At equal FDR (TPR at achieved FDR) A c2 beats edgeR in 6/6 and DESeq2 in 3/3 Stage IV cohorts, and is 0.013-0.030 below limma.
  Model B equals model A there (own-variance weight w about 0.88-0.94, d0 3.4-5.7).
* **Simulation (P+, complete null pi1 = 0):** calibration is worse than TOX-log at small n. Share of null genes with raw p <= 0.05, A (c = 1.5): 0.081-0.096 at n = 3, 0.065-0.072 at n = 5, 0.057-0.061 at n = 10 (TOX-log 0.051-0.063; limma 0.050).
  Probability of at least one BH call under the null is about 1.0 for A even at sigma_d = 0 (TOX-log: 0.00 at n = 3, 5). A wider band improves calibration monotonically while ranking barely changes (c = 1 -> 1.5 -> 2).
  At n = 3, sigma_d = 1.05 (nb): AUPRC excess TOX-log 0.320, A 0.352, B 0.362, limma 0.378; FDR 0.222 / 0.419 / 0.430.
* **Part V (5% sparse hypervariable genes):** the variance arms are worse than TOX-log (lower AUPRC at het = 0.5, higher FDR).
* **Never measured:** calibration and power on REAL data at n = 3 and n = 5, which is the most common use case and the one other tools handle with shrinkage.

**Questions this note must answer (pre-specified, section 6):**
Q1 Is variance-matched TOX calibrated on real data at n = 3 and n = 5? Q2 How much power does it have there compared with TOX-log, limma, edgeR, DESeq2?
Q3 Which band constant `c` (grid 1-4) is needed? Q4 Does shrinkage (B) help at n = 3, 5 beyond A? Q5 (engineering) readable logs, and no more stale-parameter accidents.

---

## 2. Work packages

### WP1 — Console output (note C, complete)
Implement `implementation_note_C_console_output.md` exactly (rules R1-R7, `print_tab()`, `print_pivot()`, single Slurm log). Everything added below (Part RS tables) must use these helpers.

### WP2 — Dispersion-parameter guard (prevents the stale-file accident)
Background: the P+ run used `disp_params.rds` fitted BEFORE the depth fix (pooled `c = 22.9`, per-cohort 14.3 / 35.6 / 21.2); the corrected file has `c = 19.09` (11.96 / 22.87 / 21.98). Nothing stopped the stale file from being used.
1. `estimate_dispersion_trend.R`: save `mu_depth` (the depth at which `mu` is defined; currently the literal `4e7`) in the RDS. Put the constant in ONE place sourced by both the estimator and `power_test.R` (e.g. `common/power_constants.R`, `POWER_LIB_SIZE <- 4e7`); `PCFG$lib_size` must read it.
2. `power_test.R`: when `POWER_DISP_PARAMS` is set, require `abs(mu_depth / PCFG$lib_size - 1) < 1e-9`; otherwise `stop()` with "re-run estimate_dispersion_trend.R". A file without `mu_depth` counts as stale.
3. Log header and manifest: record the parameter file path, its `date`, `phi_inf`, `c`, `sigma_hat`, `mu_depth`.
4. Test: fixture without `mu_depth` -> error; fixture with a different depth -> error; matching fixture -> passes.
5. After this change the user re-runs the estimator once (about a minute) and uses that file. Do not hard-code any parameter values.

### WP3 — New part `RS` (real data, small n)
Purpose: calibration AND power on real TCGA noise at n = 3 and n = 5 per group.

* **Registration:** add `"RS"` to `ALL_PARTS` (after `"R"`); config `PCFG$RS <- list(n_reps = c(3L, 5L), pi1 = c(0, 0.10), n_rounds = 100L, ...)` with env override `POWER_ROUNDS_RS`. Same base cohorts as Part R (the six RDS files in `POWER_REAL_RDS`; same loader, same duplicate-patient handling).
* **Rounds per part:** `build_jobs()` currently uses the global `PCFG$n_rounds` inside `add()`; make rounds a per-part setting (default `PCFG$n_rounds`). **VERIFY** that seeds remain unique: `j$seed <- PCFG$seed + 7919 * j$job + 1000003 * match(part, ALL_PARTS)` is computed from `job`, and Part R overwrites `job` per cohort. For RS use one running job counter across (cohort, n, pi1, round) so that no two jobs share a seed. pi1 = 0 and pi1 = 0.1 jobs must have different seeds.
* **Sampler:** extend `simulate_real(cohort, pi1, n = NULL)`: `n <- n %||% (S %/% 2L)`; draw `2n` distinct samples at random (the existing `sample.int(S, 2L * n)`), first `n` to group A, next `n` to group B (disjoint by construction). If `2n > S` skip the cell and record it in the manifest. Everything else (composition balancing by thinning, realised-effect bookkeeping, `n_per_group`, `n_dropped`) stays.
  **VERIFY:** `pi_hat` (composition used for balancing) is now estimated from only 2n samples; check that the realised composition shift `delta_realised` stays near 0 at n = 3 and report its median per cell.
* **pi1 = 0 path:** **VERIFY** that `draw_effects(G, 0)`, `balance_effects(<empty>, pi_hat, thinning = TRUE)` and `effect_diag(<empty>, <empty>)` work with zero DE genes; if not, special-case: no thinning, `is_de = FALSE`, `true_lfc = 0`, `delta = 0`. Metrics that need DE genes (TPR, AUPRC, |LFC|80) are NA there; calibration metrics must be computed (WP4).
* **Gene filter:** the shared `common_filter()` per round (as in Part R); store `n_kept`.
* **Methods in RS (restricted):** reference arms `limma`, `edgeR`, `DESeq2`; TOX comparator `TOX-log` with input `tpm_centred` (label `TOX-log|TPM-centred` as in Part R); variance arms (WP5). Do NOT run the 23-arm TOX input grid that Part R/N run (the `part %in% c("N", "R")` branch near line 688): give RS its own short arm list.
* **Pairing:** all methods in a round must see the identical dataset (existing design). Store a hash of the count matrix in the manifest (`data_hash`) and test that it is identical across methods of a job.
* **Manifest fields for RS:** `part`, `cohort`, `S`, `n_per_group`, `pi1`, `round`, `seed`, `patient_check`, `n_dup_removed`, `n_kept`, `delta_realised`, `data_hash`.
* **Storage and runtime estimate:** jobs = 6 cohorts x 2 n x 2 pi1 x 100 = 2400; with about 12-14 methods and roughly 12,000 genes per job that is on the order of 3-4 x 10^8 raw gene rows (rough guess, **VERIFY with the smoke run**).
  The smoke run must print (a) seconds per job per method, (b) extrapolated total runtime on 32 cores, (c) extrapolated raw size. If raw size exceeds about 5 GB (proposal), add a lean mode for pi1 = 0 jobs: store only `job, method, gene, p` (padj is recomputable) and keep the full format for pi1 = 0.1. Ask the user before enabling lean mode.

### WP4 — Metrics (in `common/power_metrics.R` / `power_report.R`; no change to existing definitions)
Computed from the per-gene raw output, so no re-run is needed for existing parts.
* `any_call_05 = 1 - zero_disc_05` as an explicit column (under pi1 = 0 this is the probability of at least one false BH call, which equals the FDR under the complete null).
* Null-calibration set, per round over null genes (all genes when pi1 = 0): `p_le01_null`, `p_le05_null` (exists), `p_le10_null`, `ks_null` (Kolmogorov-Smirnov distance of the null p-values to U(0,1)). Aggregate mean +- SE over rounds; NA p-values excluded and counted.
* Paired differences (by round, same dataset): extend `power_paired_vs_tox.csv` to RS; new `power_paired_B_vs_A.csv` with the equal-c pairs (B c1 vs A c1.5, B c2 vs A c2, B c3 vs A c3, B c4 vs A c4) for `auprc_excess_pe`, `tpr_nom_05`, `fdr_nom_05`, `p_le05_null`, `any_call_05`, with SE and t.
* **Scan-radius diagnostic** for the variance arms: the scanned slots are the `n_scanned` genes nearest in mean to the target, so `scan_radius_log2 = |log2(m_j + 1) - log2(m_0 + 1)|` of the `n_scanned`-th nearest gene (**VERIFY** that `choose_index` really visits candidates in order of increasing mean distance). Report its median per expression decile and the share of genes whose radius exceeds 2 log2 units. Also report `frac_short_pool`, `frac_skipped`, `d0`, `w`, `tau_var`.
* Unit tests with synthetic p-values: uniform p -> `p_le05_null` about 0.05 and small `ks_null`; all p = 1 -> `p_le05_null = 0`; all p = 0 -> 1; NA handling.

### WP5 — Wider band grid
Add arms (same pattern as `vmatch_arms`): `TOX-log-vband-c3`, `TOX-log-vband-c4`, `TOX-log-vshrunk-c3`, `TOX-log-vshrunk-c4` (existing: A c1, 1.5, 2; B c0.5, 1, 2). Parts: P+, R, RS (V optional). `PCFG$vmatch$parts` must include `"RS"`.
Reason: in P+ and Part R calibration improves monotonically with c while ranking barely moves, and the optimum was not reached at c = 2. **c is chosen on calibration only (section 6).**
For RS the variance arms are: A c1.5, c2, c3, c4 and B c1, c2, c3, c4 (8 arms) plus the 4 comparators.

### WP6 — Reporting and plots (RS plus the sample-size curve)
* Pivot tables (note C `print_pivot`) per part/n: rows = cohort x n x pi1, columns = methods, one metric per table: `p_le05_null`, `any_call_05`, `ks_null` (pi1 = 0) and `fdr_nom_05`, `tpr_nom_05`, `tpr_ach05`, `auprc_excess_pe` (pi1 = 0.1).
* Plots in the existing palette: (a) null calibration (`p_le05_null` with SE bars and the 0.05 line) by method and n, facet by cohort; (b) **sample-size curve:** for each cohort, metric vs `n_per_group` (log x) combining RS (3, 5) and R (13-41), lines per method, for AUPRC excess, TPR at nominal, FDR at nominal; (c) paired differences with SE (A and B vs TOX-log; B vs A).
* Documentation: add Part RS, the criteria below and the results to `docs/power_benchmark_design.md`; add the variance engine to `docs/noise_model_current_state.md` if missing.

### WP7 — Optional experiments (implement the switches, run only when the user says so)
* **E1 zero-free universe (hypothesis: exact zeros at high expression contaminate pools):** env `POWER_R_ZERO_FREE=1` adds a cohort-level filter "no zero count in any sample" for ALL methods in Parts R and RS. Compare the TOX-minus-limma gap with and without it, per expression decile.
* **E2 cross-fitted selection (only if RS fails the calibration gate at n = 3 or 5):** select case-pool neighbours by the control-group variance and control-pool neighbours by the case-group variance, so the selection variable is independent of the residuals that form the pool. This needs two matching vectors instead of one and rests on an equal-variance assumption in SELECTION only. Do not implement before the user decides.
* **E3** `zero_diagnostic.R` (already written) is run by the user on the six cohort RDS files; no code needed.

---

## 3. Experimental design for RS (statistical requirements)

* Each round: fresh disjoint subsample, fresh DE set and effects (new random draws). 100 rounds per cell: planning estimate from the existing runs, per-round SD of AUPRC excess about 0.056-0.062 (P+ at n = 3, 5 and Part R), so the SE of a method mean is about 0.006 at 100 rounds (unpaired; paired differences are tighter). The SD for real small-n subsamples may differ: report the observed SD.
* The cohort is fixed, so SEs reflect subsampling variability only. **Report every result per base cohort; never pool SEs across cohorts as if independent; state the number of base cohorts (6) next to any cross-cohort statement.**
* The smallest cohort has S = 26 (LUAD Stage IV): with n = 5 a round uses 10 of 26 samples, so rounds overlap; mention this in the report (the SE is conditional on that cohort).
* pi1 = 0 rounds (used to choose c) and pi1 = 0.1 rounds (used to judge power) are different datasets with different seeds; keep it that way.

---

## 4. Runbook (done by the user; Claude Code prepares the files)

1. **P+ re-run with the corrected parameters** (independent of the code work above): estimator output `disp_params_4e7.rds` already exists (pooled `c = 19.089`, `phi_inf = 0.1244`, `sigma_hat = 1.052`). After WP2 regenerate it once with `mu_depth`, then `PARTS="P+"` with `-e POWER_DISP_PARAMS=<new file>`. Copy the old `power_out` and `power_raw/P+` to a dated folder first (the report rewrites `power_out`; each run clears the raw folder of the parts it runs).
2. **RS smoke test:** `POWER_PARTS=RS POWER_ROUNDS_RS=2`, two cohorts (one healthy, one Stage IV), separate `POWER_RAW_DIR`/`POWER_OUT_DIR`. Claude Code creates `power_rs.sbatch` by copying `power_full.sbatch` with `PARTS="RS"`, `POWER_ROUNDS_RS=100` and WITHOUT the `--error` line (note C R7).
3. **RS full run** only after the smoke test passes and the user has seen the runtime and storage estimates.
4. Run the report; read the pivot tables and plots from WP6.

---

## 5. Acceptance tests

1. Sampler: groups disjoint, sizes n, reproducible by seed, all samples from the cohort, skip when `2n > S`; seeds unique across all RS jobs; pi1 = 0 and 0.1 jobs differ.
2. pi1 = 0 path runs end to end on one cohort; DE-dependent metrics are NA; calibration metrics are finite.
3. `data_hash` identical across methods within a job.
4. WP4 metric tests pass; paired tables contain RS rows; `power_paired_B_vs_A.csv` exists.
5. WP2 guard tests pass; the log header shows parameter file, date, `c`, `mu_depth`.
6. WP1 acceptance of note C (titles with every table, no table chunk without keys, log at most about 1,500 lines).
7. Smoke run: 2 rounds finish, prints runtime and storage extrapolation, report runs.

---

## 6. Evaluation protocol and decision criteria (PROPOSAL: the project owner confirms or edits before any run)

1. **Calibration first (RS, pi1 = 0, per cohort, n = 3 and n = 5).** Gate: `p_le05_null` of the arm is at most `TOX-log|TPM-centred` + 0.01 in at least 5 of 6 cohorts. Report `any_call_05`, `ks_null`, `p_le01_null`, `p_le10_null` as supporting (no gate: `any_call_05` is a harsh metric under the complete null).
2. **Choose c on calibration only:** for each of A and B and each n in {3, 5}, the smallest `c` in the grid that passes step 1. If none passes, report "fails calibration at n" and stop for that arm.
3. **Then power (RS, pi1 = 0.1), using the chosen c:** an arm is an improvement over `TOX-log|TPM-centred` if AUPRC excess rises by at least 0.02 in at least 5 of 6 cohorts, no cohort loses more than 0.01, and observed FDR at nominal 0.05 is at most TOX-log's + 0.02. Report TPR at nominal and at achieved FDR, FDR, |LFC|80, and the gap to limma, edgeR and DESeq2.
4. **Does shrinkage help? (Q4)** B (its chosen c) beats A (its chosen c) if AUPRC excess is higher by at least 0.02 in at least 5 of 6 cohorts and `p_le05_null` is not higher by more than 0.01. Define "helps a lot" this way BEFORE looking at the results.
5. **Report all arms** that were run, including failures. Documentation text about shrinkage ("introduces a mild assumption, helps at low n") is decided only after step 4.

---

## 7. Risks, things to verify, decisions embedded

* **Embedded defaults (change if wrong):** n in {3, 5}; pi1 in {0, 0.1}; 100 rounds; c grids A {1, 1.5, 2, 3, 4} and B {0.5, 1, 2, 3, 4}; RS arm list above; thresholds in section 6; lean-mode threshold 5 GB; part name `RS`.
* **Expectations from the toys/simulation (not targets):** gains mainly at larger n; at n = 3 neither A nor B is expected to match TOX-log's calibration (in P+ at n = 3 shrinkage adds only about +0.01 AUPRC over A); calibration may fail at n = 3 and pass at n = 5. Whatever RS shows is the answer; do not tune to the expectation.
* **VERIFY list:** pi1 = 0 handling in `simulate_real`; unique seeds; `choose_index` scan order; filter/arm `parts` logic for RS; the report's cell key columns with `cohort` + `n_per_group` for RS (no collisions with R rows); `fdr_nom_05` definition when there are zero DE genes.
* **Out of scope:** changes to the Fortran models, the multi-axis module, per-side variance matching, E2 before a decision, any claim about n = 3 before RS has been run.
