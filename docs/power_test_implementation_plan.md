# `power_test.R` — implementation plan (v3)

**Repo:** `schroederaaron/empirical-noise-model-testing`, `main` @ `dd39f00`
**Fortran API checked:** `asishallab-group/Tensor-Omics`, branch `145-empirical-noise-model` @ `1136156`
**Script:** `TCGA_test/scripts/power_test.R` · **Results reviewed:** `Simulated_data/results/power_out/`
**Status:** proposal, nothing modified in either repo. v3 records all decisions (§0), including the Part R design (option A) confirmed after v2.

**Limits of this plan.** R is not installed in the environment this plan was written in, so no R code has been executed. Numbers quoted come from the committed CSVs, or from two labelled independent re-computations (§3 I3 and §4.9), both in Python. Everything in R must be tested on your side; §7 says how.

---

## 0. Decisions recorded

| # | Topic | Decision | Consequence |
|---|---|---|---|
| 1 | Dispersion (I2) | Trend depends on mean expression; test several magnitudes of gene-wise heterogeneity around it | §3 I2 |
| 2 | Effect-size truth (I3) | Balanced draw (option b) | §3 I3 |
| 3 | Contamination diagnostic (I5) | No "pool width" output exists; use the Fortran neighbourhood size | **Partly** — size is available, width is not; see §5 |
| 4 | `tpr_ach05` (I6) | Keep as is | I6 closed; caveat documented only |
| 5 | `bimodal` (I8) | Separate stress arm, out of the headline grid | Main evidence comes from real cohorts (Part R) |
| 6 | Real cohorts (I10) | TCGA-COAD, TCGA-LUAD, TCGA-KIRC; healthy and Stage IV; counts and TPM exist, loaders exist. **Option A:** each of the six sets is a base cohort; the **full** cohort is split into random halves, a known log2FC is injected, methods are evaluated; **20 repeats** | §3 I10 |
| 7 | AUPRC naming (A1) | Rename immediately, no aliases | `ap_*` → `auprc_*` |
| 8 | Architecture | One script analyses and saves CSV; a second script computes metrics and plots from those CSVs | §1 |

---

## 1. New architecture (decision 8)

Today `power_test.R` simulates, runs every method, computes metrics and plots in one pass. Only *aggregated* metrics are saved. Per-gene p-values are discarded, so any new metric (a PR curve, a different tie-break) forces a full re-run.

### 1.1 Scripts and data flow

| Component | Location | Role | Output |
|---|---|---|---|
| `power_test.R` (analysis) | `TCGA_test/scripts/` | Build truth / load real cohort, run every method, **save raw per-gene results. No metrics.** | `power_raw/` |
| `power_report.R` (interpretation) | `TCGA_test/scripts/` | Read `power_raw/`, compute **all** metrics, summaries, paired differences, oracle gap, **and all plots** | `power_out/*.csv`, `power_out/*.png` |
| `power_metrics.R` (library) | `common/` | Pure functions sourced by `power_report.R`; exists only so they can be unit-tested without the Fortran build | none |

The name `power_test.R` is kept for the analysis script so existing commands and docs stay valid. Adding a metric then means editing `power_metrics.R` / `power_report.R` and re-running only the report, which needs no Fortran and no reference-method runs.

### 1.2 Raw output (analysis script)

Per job (one simulated or thinned dataset), under `power_raw/<part>/`:

| File | Content |
|---|---|
| `results_<job>.csv.gz` | one row per (method, gene): `p`, `padj`, `lfc`, `tie`, `floor_th`, and for TOX arms `nbhd_own_case`, `nbhd_own_control` |
| `truth_<job>.csv.gz` | one row per gene, **method-independent** (stored once): `is_de`, `true_lfc`, `expr`, `kept` (post-filter) |
| `diag_<job>.csv.gz` | Part Q only: pool diagnostics (§5) on a gene subset |

Plus one manifest `power_jobs.csv`: `job_id, part, dist, n_rep, pi1, het, frac_up, sigma_d, cohort, round, seed, delta, delta_realised, n_de_total`, the effective-effect-size diagnostics from I3 (`eff_lfc_median_abs`, `frac_sign_flip`, `frac_below_floor`), and the **git commit hash and full `PCFG` dump** so every result is traceable to code and parameters.

- **Format.** `csv.gz` written with `data.table::fwrite` (already a dependency elsewhere in the repo). Raw files are **not committed** (add `power_raw/` to `.gitignore`); the metric CSVs and PNGs in `power_out/` stay committed as now.
- **Unavoidable re-run.** Because per-gene results were never saved, the committed `power_out` cannot feed the new report script. One full re-run of the new analysis script is required. It is also required anyway by I2, I3 and I4, so nothing extra is lost. After it, new metrics never need a re-run.

### 1.3 Report script behaviour

`power_report.R` reads `truth_*` + `results_*` (+ manifest), recomputes BH on each method's `p` where needed, and writes exactly the files the current script writes (`power_per_round`, `power_summary`, `power_paired_vs_ref`, `power_fdr_tpr_curve`, `power_stratified`, `power_oracle_gap`), plus the new AUPRC outputs (§4) and every plot. It can be run repeatedly, on a laptop, in minutes.

---

## 2. Issue register (updated)

| ID | Issue | Severity | Status after decisions |
|---|---|---|---|
| I1 | `fdr_nom_*` averaged over rounds with ≥1 call only | **High** | Fix in `power_metrics.R` |
| I2 | Constant dispersion → TOX's exchangeability holds by construction | **High** | Mean-dependent trend + heterogeneity magnitudes (decision 1) |
| I3 | Effective TPM-space effects ≠ documented range; ~23% sign flips | Medium | Balanced draw (decision 2) |
| I4 | Part C: Δ ≠ 0 at up-fraction 0.5, Δ range too narrow | Medium | Δ-targeted redesign |
| I5 | P2/P3 specified on a metric that cannot test them | Medium | Re-specify; diagnostics per §5 |
| I6 | `tpr_ach05` oracle-thresholded | Low–Medium | **Closed** — kept as is (decision 4); caveat documented |
| I7 | DESeq2 ≈2× nominal FDR on its own model, unexplained | Medium | Diagnostic run |
| I8 | `bimodal` anomalies | Medium | **Moved to stress arm S** (decision 5) |
| I9 | Documentation drift | Low | Fix |
| I10 | Real-data arms not run | **High** | **Promoted to primary evidence** (decisions 5, 6) |
| A1 | AUPRC | — | §4 |

---

## 3. Issue details

### I1 — FDR averaging bias (High)

**Evidence.** `metrics_one()` returns `fdr_nom_* = NA` when a method makes no calls, and `mean_se()` averages finite values only, giving E[V/R | R>0] instead of E[V/max(R,1)]. 29 of 200 (part, cell, method) rows at α = 0.05 are affected. Example: `TOX-log-boot`, `bimodal`, n=3 reports 0.570 (3 of 20 rounds have calls) vs **0.086** with zero-call rounds counted as 0.

**Fix** (in `power_metrics.R`): `fdr_nom = fp / max(1, tp + fp)`; keep `zero_disc_*`; add `fdr_cond_*` (the old conditional value) for anyone who wants it; print `frac_zero_disc` next to every FDR in the report and plot subtitles.

**Acceptance.** Unit test: a round with zero calls contributes 0. The 29 affected cells in the old results, recomputed from `power_per_round.csv`, match the zero-filled values (documents the size of the old error; no need to publish corrected old numbers).

### I2 — Heterogeneous, mean-dependent dispersion (High)

**Problem.** `simulate_power()` uses one scalar `disp = 0.2` for every gene. TOX pools residuals from mean-neighbouring genes and treats them as exchangeable; with identical dispersion and gene-independent noise that is exact by construction. Parity with limma on `nb`/`lnpois` does not test the assumption most likely to fail on real data (project evidence: edgeR is 2.2–2.8× anti-conservative on all twelve TCGA groups yet well calibrated on the synthetic arms).

**Model (decision 1).** Per-gene dispersion

`φ_g = φ_trend(μ_g) · exp(σ_d · z_g − σ_d² / 2)`, `z_g ~ N(0, 1)`,

- `φ_trend(μ)` decreases with mean expression, in the standard parametric form `φ(μ) = φ_∞ + c / μ`.
- The `−σ_d² / 2` term makes the multiplier **mean-one**, so raising σ_d increases *scatter* around the trend without also raising the average dispersion. Without it, σ_d would be confounded with overall noise level.
- σ_d = 0 gives the trend only (no gene-wise scatter); it is the control and is the closest analogue of today's run once `φ_trend` is flat.

**Anchoring the parameters in data (proposed default).** Estimate `φ_∞`, `c` and σ̂_d once from the TCGA healthy cohorts (COAD, LUAD, KIRC; edgeR trended vs tagwise dispersion). Caveat that must be recorded with the estimate: gene-wise dispersion estimates carry their own sampling noise, which inflates apparent scatter, so σ̂_d is an **upper bound** on the true heterogeneity. The magnitudes tested would be `σ_d ∈ {0, 0.5·σ̂_d, σ̂_d, 2·σ̂_d}`, bracketing the estimate instead of picking arbitrary values.

**Mechanics to verify before use.** `draw_counts()` takes `disp` as a scalar. For the `nb` and `lnpois` branches a per-gene vector expanded to the G×S layout should recycle column-major, so gene g uses its own value in every sample. This must be unit-tested (correct per-gene variance, checked empirically) before any run. `bimodal` has no dispersion parameter and is unaffected.

**Runs.** Both a **null arm (π₁ = 0, calibration)** and the power arm at every σ_d. Power without calibration under heterogeneity is not interpretable. Report each method's FDR at nominal α (zero-filled per I1), TPR at achieved FDR, and the paired difference against limma.

**Acceptance.** At σ_d = 0 with a flat trend the results reproduce the current run within Monte-Carlo error (regression check). Whether TOX degrades with σ_d is an open empirical question; the plan does not presuppose the answer.

### I3 — Balanced effect-size draw (decision 2)

**Evidence** (my Python re-implementation of the truth construction only, not the R code; G = 5000, π₁ = 0.10, 2000 replicates). Renormalising within the DE set forces `Σ π_a 2^(effective LFC) = Σ π_a` over that set, so the effective TPM-space LFC is `β + log2(cs)` with `log2(cs) < 0`:

| Quantity | Result |
|---|---|
| median `log2(cs)` | −0.89 (IQR −1.12 to −0.68) |
| median \|effective LFC\| | 1.34 (nominal geometric median 1.0) |
| DE genes with \|effective LFC\| < 0.25 | ≈10% (5th–95th pct 5–19%) |
| DE genes with sign flipped vs β | ≈23% (5th–95th pct 8–34%) |
| rounds with ≥10% of DE genes flipped | 93% |

**Implementation.**
1. Draw the *effective* signs and magnitudes first: sign symmetric, |LFC| log-uniform on `lfc_range`.
2. Enforce `Σ_{de} π_a · 2^e = Σ_{de} π_a` (nulls exactly null in TPM) by scaling the magnitudes of the side that would otherwise dominate, using the same `uniroot` pattern already in `simulate_real()` (Part R). The realised distribution then matches the target instead of being shifted.
3. The existing guard `stopifnot(max(abs(true_lfc[!is_de])) == 0)` stays.
4. **Known trade-off, must be reported, not hidden:** when a few very abundant genes dominate the DE set, one side may need substantial shrinking, so the realised |LFC| distribution can fall below the target range. Record `eff_lfc_median_abs`, `frac_sign_flip`, `frac_below_floor` per job (manifest) and report them beside every result.

**Acceptance.** Over ≥1000 simulated truths: sign split ≈ 50/50 in TPM space, `frac_sign_flip` ≈ 0 by construction, realised median |LFC| within a stated tolerance of the target, and the null guard holds in every round.

### I4 — Part C: Δ-targeted design

**Evidence.** Realised Δ (mean ± SD): 0.128 ± 0.095 / 0.181 ± 0.066 / 0.252 ± 0.041 / 0.288 ± 0.075 at up-fraction 0.5 / 0.7 / 0.9 / 1.0. It is never 0, and TOX-log on raw TPM is already at 2.4× nominal FDR at the first point, so the design doc's breakdown point (§4.2: smallest Δ with FDR > 2× nominal) cannot be located.

**Plan.**
1. Parameterise by **target Δ** ∈ {0, 0.25, 0.5, 1.0}: construct β, then solve for the scaling of one side's magnitudes that yields the target Δ (same root-finding pattern as I3). Δ = 0 is the reference, using the I3 balanced construction.
2. Add n_rep ∈ {3, 10} beside 5, and `nb` beside `lnpois` (currently `lnpois`, n=5 only).
3. Breakdown point defined explicitly: smallest Δ where zero-filled FDR at nominal 0.05 exceeds 0.10, with a bootstrap CI over rounds.

### I5 — Re-specify P2/P3 (see §5 for the diagnostics)

**Evidence.** Both predictions were stated on TPR at achieved FDR, which (a) rises with π₁ even without contamination, because FDR = π₀·FPR / (π₀·FPR + π₁·TPR), and (b) is invariant to a monotone inflation of p-values, which is how contamination acts. Observed: oracle − production ≈ 0 in achieved-FDR TPR (|Δ| ≤ 0.007, all |t| < 1.9) but +0.005 → +0.092 nominal TPR (SE ≤ 0.004) and +0.004 null FPR at het = 0.5, π₁ = 0.30.

**Plan.** Rewrite `docs/power_benchmark_design.md` §3.2 and §7:
- **P2′:** paired TOX − limma contrast in TPR at achieved FDR versus π₁, at het = 0.5 and het = 1. (Observed −0.027 → +0.161 at het = 0.5, opposite to the original P2. The mechanism behind TOX's advantage is an **untested hypothesis**, not a finding.)
- **P3′:** contamination is a *calibration* effect. Test it on paired oracle − production differences in zero-filled `fdr_nom_05`, `fpr_null05`, `n_called_05`, `tpr_nom_05`, and on the pool diagnostics of §5, with achieved-FDR TPR as the expected-null control.

### I6 — Closed (decision 4)

`tpr_ach05` is kept as is. Documentation only: state in the design doc that it is truth-thresholded and takes the maximum over all steps with realised FDR ≤ 0.05 (slightly optimistic on a non-monotone curve, and not user-realisable); `tpr_nom_*` is the user-facing number.

### I7 — DESeq2 diagnostic

`nb`, α = 0.05: DESeq2 observed FDR 0.107 / 0.103 / 0.092 at n = 3 / 5 / 10 vs limma 0.039–0.044 and edgeR 0.041–0.050. `run_deseq2()` uses `independentFiltering = FALSE`, otherwise defaults. Cause not established. **Diagnostic only; the reference arm is not changed until diagnosed:**
1. π₁ = 0 null run on `nb` with the current wrapper (p-value uniformity).
2. Variants: `independentFiltering = TRUE`; `cooksCutoff = FALSE`; record the number of Cook's-flagged NA p-values per round and `n_tested`, to rule out a different gene set being scored.

### I8 — `bimodal` becomes stress arm S (decision 5)

- Remove `bimodal` from the Part P grid. Add **Part S (stress)**: `bimodal` only, same n_rep grid, reported separately and labelled as a stress case in every table and plot. The simulator is not close to real data, and TOX's resampling is known to suffer from it (your assessment; it is consistent with the observed TOX-log ≈ 0 at n=3).
- The unexplained limma behaviour (power falling with n: 0.453 / 0.328 / 0.295) is **not investigated now**. It is recorded as a known open observation in the design doc so nobody interprets Part S as a ranking.
- Headline power claims are drawn from Part R (real cohorts), not from Parts P/Q/C/S.

### I9 — Documentation fixes

`docs/power_benchmark_design.md` and the script header: script path (`Simulated_data/scripts/` → `TCGA_test/scripts/`); `claude/power_benchmark_design.md` → `docs/power_benchmark_design.md`; "not yet run, no result exists" → point to `power_out` and its parameters; remove the phantom `power_fidelity.csv`; document the two-script architecture (§1) and the decisions in §0.

### I10 — Real cohorts as primary evidence (decisions 5, 6)

Part R (binomial thinning on real data: real correlation, real dispersion, real outliers, exactly known truth) becomes the main power evidence.

**Data.** `TCGA-COAD`, `TCGA-LUAD`, `TCGA-KIRC`; `healthy` (matched normals) and `Stage IV`. Counts via `load_counts_matrix(pid, stage)` (genes × samples), TPM via `load_stage_data(pid, ..., raw, apply_mean = FALSE, normalize = FALSE)` (samples × genes). `load_real_cohort()` already orients TPM by dimnames.

**Design (option A, confirmed).** For each of the six base cohorts (3 cancers × {healthy, Stage IV}):
1. Take the **full** cohort. If the sample count S is odd, drop one sample at random so both halves have `n = floor(S/2)`. (The number dropped is recorded in the manifest.)
2. Randomly split the samples into two halves A and B.
3. Draw the DE set and effective log2FCs with the balanced construction already used in `simulate_real()` (Part R; the `uniroot` step equalises the expected relative abundance removed from each side, so Δ = 0 and true nulls stay null in TPM), and inject them by binomial thinning.
4. Run every method on the thinned counts / TPM and save raw per-gene results (§1).
5. Repeat **20 times**, each with a fresh random split *and* a fresh DE set / effect draw.

**Consequences of using the full cohort (stated so they are not discovered later).**
- **n per group is cohort-specific** (≈ S/2), no longer the fixed grid {3, 5, 10}. The subsampled small-n arms are therefore *removed from Part R* by this decision; the sample-size curve now comes only from the synthetic Part P and from the different cohort sizes. Say if you want a subsampled arm back.
- **What 20 repeats do and do not reduce.** They average out the variance from the random split and the random injection. They do **not** reduce cohort-level uncertainty: the cohort is fixed and the 20 splits reuse the same samples, so the SE across repeats is the SE *conditional on that cohort*. Results are reported **per base cohort**; six cohorts is the effective sample for any cross-cohort statement, and SEs are not pooled across cohorts as if independent.
- **Large n changes method behaviour.** At larger n most moderate effects become detectable, so TPR approaches its ceiling and the small-|LFC| bins carry the discriminating information (the `lfc_bins` and `lfc80` outputs matter more here than the aggregate TPR). The gene-blocked bootstrap arms enumerate exactly only up to n_rep = 5 and fall back to sampling (p floor ≈ 1e-4, `m*` ≈ 10–25) above it; DESeq2 runtime grows with n. The `floor_blocks05` / `m_star` columns are kept for exactly this reason.
- The oracle arm and pool diagnostics (§5) are Part Q only and are not run on real cohorts.

**Checks that must be added before trusting any Part R number.**
1. **Length recovery.** Part R recovers gene lengths as the per-gene median of count/TPM after removing each sample's constant. That is only valid if TPM was derived from *these* counts. After per-sample scaling, the per-gene spread of count/TPM across samples must be ≈ 0; report the fraction of genes above a tolerance and drop them. Otherwise the TPM built from thinned counts is on a different scale from the TPM the null is calibrated on.
2. **Exchangeability of the two halves.** A random split makes true nulls null only if samples are independent. If several samples come from the same patient, halves are not exchangeable and variances are understated. Add a check for one sample per patient; the sample-ID format must be confirmed on the compute host before the check is written (it is not visible from the repositories). Cohorts that fail are reported, not silently used.
3. **Realised vs intended fold change.** Thinning can only remove reads, so a requested fold change may not be realised for low-count genes. For every injected gene compare the realised log2 ratio of group means to the intended value and store both; the report flags genes where they disagree instead of assuming the design was achieved. `delta_realised` (median null log-ratio, target 0) is reported per round; rounds where balancing failed (`delta = NA`) are counted and listed.
4. **Gene universe.** Keep the shared `common_filter()` universe for the method comparison, and report TPR over all true-DE genes (`n_de_total`) so TOX's stricter filter is charged honestly (design doc §4.3).

**Manifest fields added for Part R:** `cohort` (e.g. `TCGA-COAD/healthy`), `S`, `n_per_group`, `n_dropped`, `split_seed`.

**P5** (does edgeR's advantage shrink or reverse on real thinned data, as its calibration did between `nb` and TCGA?) becomes directly testable here.

---

## 4. A1 — AUPRC

### 4.1 Current state

`avg_precision()` already computes a **step-function AUPRC** (average precision): `Σ_k (R_k − R_{k−1}) · P_k` over tie-block ends. It is written as `ap_p` (p-value only, ties stay tied) and `ap_pe` (p-value, then |log₂FC of group means| as tie-break).

Independently validated (Python port of `roc_steps()` + `avg_precision()`, not the R code): tie-free scores match scikit-learn's `average_precision_score` (0.36897 vs 0.36897); all-tied scores give exactly the prevalence (0.09905); random ranking averages 0.1008 vs prevalence 0.0991; a perfect ranking gives 1.0.

Missing: a name that says what it is, a prevalence baseline, a PR curve, plots, tests, and inclusion in the oracle-gap and spread reports.

### 4.2 Definition (fixed)

- **Estimator:** step function, no interpolation. Linear interpolation is not valid in PR space and can be optimistic with sparse points or large tie blocks; no trapezoid variant is added.
- **Ranking:** `-p` for `_p` (ties are one step; no threshold splits a tie) and `score_pe` for `_pe`. The gap `_pe − _p` measures resolution loss from the discrete p grid, consistent with AUC/pAUC.
- **Ties:** precision at a block end is what a threshold can deliver. It is slightly pessimistic vs random tie-breaking (test: 0.3418 vs 0.3443 for a 200-gene tied block holding 103 DE genes). No third variant.
- **Untested genes** rank last, tied. **Recall denominator:** tested DE genes `P` for `auprc_*`, plus `auprc_all_*` with `n_de_total` (filtered-out DE genes never retrieved), mirroring `tpr_all_*`.

### 4.3 Columns (renamed immediately; no `ap_*` aliases)

| Column | Definition |
|---|---|
| `auprc_p`, `auprc_pe` | step AUPRC (numerically identical to the old `ap_p`, `ap_pe`) |
| `auprc_all_p`, `auprc_all_pe` | recall denominator `n_de_total` |
| `prevalence` | `n_de / n_tested` (random-ranking baseline) |
| `auprc_excess_pe` | `auprc_pe − prevalence` |
| `auprc_lift_pe` | `auprc_pe / prevalence` |
| `prec_top_nde` | precision at k = `n_de` (R-precision) |

**Why the baseline.** Prevalence spans 0.01 / 0.05 / 0.10 / 0.30 in Part Q (measured 0.0100 / 0.0498 / 0.0996 / 0.2977); a random ranking scores ≈ prevalence, so raw AUPRC is not comparable across those cells. Within one cell every method shares one prevalence, so raw and excess give the same order; the transform matters for **cross-cell** figures. `prec_top100/250/500` are capped at ≈ 0.48 / 0.19 / 0.096 at π₁ = 0.01 (~48 DE genes), so `prec_top_nde` is added.

### 4.4 PR curve

`curve_pr_one()`: precision on a recall grid `seq(0, 1, 0.05)`, using the interpolated-precision envelope (maximum precision over steps with recall ≥ r), for **display only**. The scalar AUPRC uses the raw step estimator and is never computed from the envelope. Defined everywhere because untested genes form a final tied block that brings recall to 1 at precision = prevalence. Averaged across rounds by vertical averaging → `power_pr_curve.csv` (mean ± SE per recall value), keyed by `(part, dist, n_rep, pi1, het, frac_up, sigma_d, cohort, method, recall)`.

### 4.5 Plots (all produced by `power_report.R`)

1. `pr_curves.png`: facets by cell, fixed method palette, dashed line at prevalence (Parts P and R; Part S separately, labelled stress).
2. Part Q: AUPRC **excess** vs π₁ (log-x, by het, oracle arms included), added to `power_vs_pi1.png`.
3. AUPRC excess added next to AUC and pAUC in `auc_vs_pauc.png` so the P6 question (does the metric separate methods?) covers AUPRC.

### 4.6 Uncertainty, comparison, reporting

Mean ± SE across rounds on every value; method comparisons use the **paired per-round** difference. Add `auprc_pe` and `auprc_p` to the paired-vs-reference and oracle-gap tables, and an `auprc_range` column to the P6 spread table. Cells with π₁ = 0.01 (~48 DE genes) are noisy for any AUPRC: SE beside every value, no ranking without paired differences.

### 4.7 Acceptance

1. `auprc_p` / `auprc_pe` equal the old `ap_p` / `ap_pe` on the same inputs (regression test against a Python-generated fixture).
2. All-tied input returns the prevalence exactly; perfect ranking returns 1.
3. Every `power_pr_curve.csv` recall grid reaches 1.0 with precision ≥ prevalence.
4. A random-score method has `auprc_excess_pe ≈ 0` within SE in a smoke run.

---

## 5. Contamination diagnostics (decision 3, refined)

**What the Fortran returns.** Each pipeline call returns three integer vectors per gene: `neighborhood_size_own_case`, `neighborhood_size_own_control`, `neighborhood_size_case`. In `tox_noise_model_exact.F90` these are the **pool size in residuals** (`n_pool_case`, `n_pool_control_own`). They are available today for the exact and bootstrap engines, and `power_test.R` already reads two of them for the p-value floor. They are **not** a measure of pool *width* (spread).

**What that means for P3′.** Pool size answers whether contamination changes the adaptive growth (does the pool stop at a different size?). It cannot show a widened null on its own.

**Width diagnostic.** `common/tox_null_reimpl.R::tox_diagnose()` — an R port of the exact model with its own `validate_against_fortran()` gate — returns per gene `sd_pool_case/control`, `sd_null`, `n_resid_pool_*`, `mean_span_*`, pool kurtosis and the spread of per-neighbour SDs.
- It applies to the **exact engine, pooled null only** (not the blocked arms).
- **The validation gate must pass (max |Δp| ≈ 0) on the simulated matrices before any diagnostic is used.**
- Cost is an R loop per gene, so run it on a fixed random gene subset per Part Q job, saved as `diag_<job>.csv.gz`.

**Confound to control for.** The oracle removes DE genes from the candidate set, so its k nearest neighbours lie further away in mean expression (at π₁ = 0.30, about 1/(1 − 0.30) ≈ 1.4× wider span in gene-rank terms). Oracle − production differences could partly reflect neighbour distance, not contamination. `mean_span_*` from `tox_diagnose()` measures this directly; report it alongside `sd_pool_*` and include it in the interpretation.

---

## 6. Parts overview after the plan

| Part | Purpose | Data | Role |
|---|---|---|---|
| **R** | Power with known truth on real noise | COAD / LUAD / KIRC × {healthy, Stage IV}; full cohort, random halves, 20 repeats | **Primary evidence** |
| P | Mechanism, effect-size and n_rep curves | `nb`, `lnpois`, `tpois` (balanced draw, I3) | Supporting |
| P+ | Dispersion heterogeneity (null + power) | σ_d grid (I2) | Supporting; tests TOX's key assumption |
| Q | π₁ × het, oracle arm + diagnostics | `lnpois` | Mechanism |
| C | Composition stress | Δ-targeted (I4) | Limitation, with a number attached |
| S | Stress | `bimodal` | Stress only; not for ranking |

---

## 7. Sequencing and testing

| Phase | Work | Needs |
|---|---|---|
| **1 — refactor** | Extract `common/power_metrics.R`; split `power_test.R` into analysis + `power_report.R`; implement I1 (zero-filled FDR) and A1 (AUPRC, PR curve, plots) | Tests only, no full run |
| **2 — simulator** | I3 balanced draw; I2 mean-dependent heterogeneous dispersion (with the recycling unit test); Part S; manifest fields | Smoke: `POWER_GENES=1000 POWER_ROUNDS=2 POWER_CORES=4` |
| **3 — Part R** | Full-cohort half-split logic (drop-one if S odd); length-recovery, one-sample-per-patient and realised-vs-intended-LFC checks; cohort RDS builder for the six base cohorts; smoke run, then full run | Access to `/media/BioNAS2/...` data on the compute host |
| **4 — diagnostics** | I7 (DESeq2), §5 diagnostics, I4 Δ-targeted Part C | Compute host |
| **5 — full run and report** | One full analysis run → `power_report.R` → update design doc (I5, I9, I6 caveat) | 20 rounds, 32 cores |

**Testing without R here.** Unit tests in R (`testthat`) for `power_metrics.R`, to be run on your side: perfect ranking → AUPRC 1; all tied → exactly the prevalence; random ≈ prevalence; a small tie-free case against an O(n²) brute-force reference; zero-call round contributes FDR 0; per-gene dispersion recycling; balanced-draw guard. Fixtures for the numerical checks will be generated with Python ports (scikit-learn for AP), so expected values are independent of the R implementation. The Python validation verifies the algorithm; the R tests verify the R code.

---

## 8. Open points

None blocking. The Part R interpretation is settled (option A, §0 row 6). Two things that depend on the compute host, not on a decision: the per-cohort sample counts S, and the sample-ID format used by the one-sample-per-patient check (§3 I10, check 2). Both are read from the data during Phase 3.

---

## 9. Not in scope

- Changes to the Fortran / TOX null models.
- Investigating why limma's ranking worsens with n on `bimodal` (recorded as an open observation only).
- Replacing the in-script tie-aware metrics with iCOBRA.
- Any conclusion about TOX vs limma on real data before Part R has run.
