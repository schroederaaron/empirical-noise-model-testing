# Implementation brief — Omnibus goodness-of-fit test for RNA-seq count distributions (Tensor-Omics issue #188)

**Audience:** Claude Code. **Author of the analysis:** Aaron. **Status:** 23.09.2026.
This brief supersedes `claude_omnibus_gof_sketch.md` (16.09.2026).

---

## 0. Ground rules — read first, follow throughout

1. **Git: pull only.** Never `git commit`, never `git push`, never create or delete branches,
   and never open PRs. Aaron does all commits himself. Deliver changes as new files in the
   working tree **plus** one patch file:
   `git diff --no-color > gof_188.patch` (tracked changes), with new files listed at the top
   of the patch or delivered alongside it.
2. **No results are claimed that were not produced.** Every number in any README, doc, log
   or summary must come from a file this code wrote. If a run was not done, write
   "not run".
3. **Validation gates are hard stops.** If a gate in §9 fails, stop. Report what failed,
   with numbers. Do not loosen tolerances, drop genes, or switch methods to make a gate
   pass.
4. **Do not use textbook p-values** of CvM/AD/KS as results. They may be logged only as
   a diagnostic, clearly labelled as uncalibrated.
5. **Do not silently round counts.** Non-integer input is an error unless the config
   explicitly selects a rounding mode (§4.2).
6. **Do not use `glmmTMB`'s Dunn–Smyth residuals** (`residuals(type = "dunn-smyth")`) as
   the PIT source. They return `qnorm(u)` and replace infinite values by 0, which turns a
   PIT of exactly 0 or 1 into 0.5 (read in `glmmTMB/R/methods.R`, `pit_norm_resids`).
   Families like `compois` are not covered by them either.
7. **Match the repository's conventions** (§2): `load_or_install()`, `source("config.R")`,
   `mclapply` with `OMP_NUM_THREADS=1`, results under `results/` which is git-ignored.
8. **Record provenance** in every output directory. That means `sessionInfo()`, package
   versions (`glmmTMB`, `TMB`, `edgeR`, `tximport`), the full resolved config, RNG kind
   and base seed, and the git HEAD of the repo.

---

## 1. What was checked for this brief (23.09.2026)

| source | state |
|---|---|
| Issue #188 | text unchanged since opening (`updatedAt 2026-09-09T16:35:41Z`), 0 comments, no linked branch or PR. The issue text is the specification. |
| `schroederaaron/empirical-noise-model-testing@main` (HEAD `b7f4df4`, 23.09.2026) | no Salmon/`tximport` code anywhere. The closest existing script is `TCGA_test/scripts/count_distribution.R`: per-gene intercept-only Poisson vs NB (`MASS::glm.nb`), offset `log(colSums(counts))` (raw library size, not TMM), cap of 1,500 genes, LRT/AIC. The new work extends that question; it does not replace the script. |
| `asishallab-group/Tensor-Omics` | no branch for #188. The analysis is pure R and does not need the Fortran build. |
| `glmmTMB` master (`44370fc`, 21.09.2026, version 1.1.15) | exports `dgenpois`/`pgenpois`/`rgenpois` (GP) — but **no** d/p/r functions for `compois`. Dunn–Smyth residuals behave as in rule 6. |

**Checks run on the installed `glmmTMB` 1.1.8** (Ubuntu package, R 4.3.3, 1 CPU core, simulated NB data, n = 200 for the logLik checks). Each check compares the sum of an own-written log-pmf, evaluated at the `glmmTMB` estimates, with `logLik(fit)`:

| family | own pmf parameterisation | own log-pmf sum vs `logLik` |
|---|---|---|
| `nbinom2` | `dnbinom(y, mu = fitted, size = sigma(fit))` | −789.6286 vs −789.6286 |
| `nbinom2` + `ziformula = ~1` | `π = predict(type = "zprob")`, `μ = predict(type = "conditional")` | −789.6286 vs −789.6286 (π̂ = 3·10⁻¹⁰ on data without excess zeros) |
| `genpois` | Joe–Zhu: `α = 1 − 1/√φ`, `λ1 = μ(1−α)`, `λ2 = α`, with **`φ = predict(type = "disp")`** | −804.1559 vs −804.1559 (with `φ²` instead: −914.38, i.e. wrong) |
| `compois` | mean-parameterised, `ν = 1/sigma(fit)`, `λ` solved numerically so that `E[Y] = μ` | −799.7514 vs −799.7514 (μ ≈ 20–35) |

**Fit time per gene** (n = 96, 2 conditions, NB-simulated counts, same 1-core container):

| family | μ = 5 | μ = 50 | μ = 500 | μ = 5000 |
|---|---|---|---|---|
| poisson / nbinom2 / zinb / genpois / PLN-Laplace | 0.15–0.24 s | 0.15–0.20 s | 0.15–0.20 s | 0.16–0.20 s |
| **compois** | 0.26 s | 0.90 s | **5.6 s** | **50 s, `convergence ≠ 0`** |

This is one gene per cell, so it is indicative only. But it means **`glmmTMB::compois` in
1.1.8 is not usable for highly expressed genes.** §5.6 specifies an own CMP fitter. Re-time
with whatever `glmmTMB` version is installed on the analysis machine; the dev version may
behave differently. This was not tested.

---

## 2. Where the code goes

In `empirical-noise-model-testing`, alongside the existing analyses:

```
Distribution_GOF/
  README.md                     # how to run; decisions; NO results unless produced
  config_gof.R                  # every knob, env-var overridable, printed to output
  R/
    io.R                        # input contract, tximport import, TCGA loader, integer audit
    preprocess.R                # filter, offsets, rounding modes, gene subset
    families.R                  # six family specs sharing one interface (§5)
    cmp_core.R                  # CMP normalising constant, λ(μ,ν), pmf/cdf/rng, own ML fit
    pln_core.R                  # PLN adaptive Gauss–Hermite pmf/cdf, own ML fit
    fit_dataset.R               # per-gene loop, error capture, parallel
    pit.R                       # PIT bounds, randomisation, re-randomisation
    gof_stats.R                 # CvM, AD, KS (+ N-normalised)
    strata.R                    # stratification rules, applied identically to every dataset
    bootstrap.R                 # parametric bootstrap
    heldout.R                   # K-fold predictive log score
    plots.R
    report.R
  scripts/
    00_validate_families.R      # gates G1–G4
    01_calibration_sim.R        # gate G5: does the whole test work when the truth is known?
    02_run_gof.R                # main analysis for one dataset
    03_report.R                 # table + plots from saved outputs (no refitting)
  tests/testthat/               # unit tests for gof_stats, pit, families, cmp_core, pln_core
  results/.gitkeep
```

Add to `.gitignore` (in the patch): `Distribution_GOF/results/*` and
`!Distribution_GOF/results/.gitkeep`, mirroring the existing entries.

Package loading: copy the `load_or_install()` / `.missing_pkg()` block from
`TCGA_test/scripts/count_distribution.R` verbatim, with the same
`LIB_DIR = external/docker_r_libs`. Packages needed: `glmmTMB`, `edgeR`, `tximport`,
`data.table`, `ggplot2`, `goftest` (tests only), `parallel`, `testthat`.

---

## 3. Input contract

`02_run_gof.R` consumes one **dataset object** (`.rds`):

```r
list(
  counts   = <numeric matrix, genes x samples, rownames = gene ids, colnames = sample ids>,
  samples  = <data.frame, one row per column of counts; includes the design columns>,
  design   = <formula over columns of `samples`, e.g. ~ condition>,
  lengths  = <optional numeric matrix like counts (tximport $length); NULL otherwise>,
  source   = <"salmon_tximport" | "featureCounts" | "tcga_counts_rds" | "simulated">,
  label    = <short id used in output paths>
)
```

Builders in `R/io.R`:

- `build_from_tximport(quant_files, tx2gene, samples, design, label)`. Calls
  `tximport(..., type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")` and keeps
  `$counts` and `$length`. **Before writing the import code, read the installed tximport
  vignette.** Copy the edgeR offset recipe from it rather than from memory.
- `build_from_featurecounts(file, samples, design, label)`.
- `build_from_tcga(project_id, stage, label)`. Reuses `load_counts_matrix()` from
  `count_distribution.R` (samples × genes `expression_vectors` → genes × samples). Design
  `~ 1` (single cohort), as in `count_distribution.R`. The README must state that
  unmodelled between-patient heterogeneity counts as lack of fit for every model.

**Integer audit** (always run, written to `input_audit.csv`):
- fraction of non-integer entries (`abs(y − round(y)) > 1e-8`);
- the distribution of the fractional part;
- both of the above per mean-expression decile;
- min/max count and library sizes.

Do **not** assume the TCGA `_counts.rds` files are integers. Report what the audit finds.

Quantifying FASTQs (Salmon / STAR + featureCounts) is **out of scope** for this code. The
prerequisites are listed in §11.

---

## 4. Preprocessing (identical for all six models)

### 4.1 Gene filter
`edgeR::filterByExpr(DGEList(counts), group = <design factor>)`. This is the downstream
criterion the issue asks for, not a stricter one. It is applied **once, to the observed
data**. The gene set is then fixed for all models, all bootstrap replicates and the
held-out score. Record the number of genes before and after filtering.

### 4.2 Fractional counts — `ROUNDING` (no silent default)

| value | behaviour |
|---|---|
| `"error"` (default) | stop if the integer audit finds non-integer entries |
| `"stochastic"` | `y* = floor(y) + Bernoulli(y − floor(y))`, applied **once**, before fitting, with its own recorded seed. Unbiased in the mean; adds ≤ 0.25 variance per entry. |
| `"floor"` | `y* = floor(y)`. Sensitivity arm only; biased downward. |

Background: R's `pnbinom`/`ppois` floor a non-integer `q`, so "not rounding" is implicit
flooring (verified: `pnbinom(3.4, …) == pnbinom(3, …)`). Which mode is primary is Aaron's
and Asis's decision (§12, D1). Write the mode into every output file name
(e.g. `…_round-stochastic`).

### 4.3 Offsets — `OFFSET_MODE`

| value | offset `o_gi` |
|---|---|
| `"tmm"` | `log(lib.size_i × norm.factors_i)` from `edgeR::normLibSizes` on the filtered counts |
| `"tmm_length"` | the tximport-vignette edgeR recipe (gene × sample offset from `$length`); only allowed if `lengths` is non-NULL |

The existing `count_distribution.R` uses `log(colSums(counts))`. That is not TMM and is
not reused.

### 4.4 Gene subsets
- `G_ALL`: every filtered gene. Used for the observed statistics and plots.
- `G_BOOT`: a fixed random subset, stratified by mean-expression decile, of size
  `N_GENES_BOOT` (default 2000; recorded seed). The bootstrap test and the held-out score
  run on `G_BOOT`. **`T_obs` for the bootstrap p-value is computed on `G_BOOT`**, never on
  `G_ALL`. Report the `G_ALL` statistics descriptively.

---

## 5. Families — one interface, six implementations

Every family is a list with this interface. It is used identically for fitting, the PIT,
the log score and simulation:

```r
fam <- list(
  name     = "nb",
  n_extra  = 1L,                         # gene-level parameters beyond the regression coefs
  fit      = function(y, X, off) -> list(beta, extra = list(...), mu = <per obs>,
                                          loglik, conv_ok = <lgl>, conv_msg, engine),
  logpmf   = function(y, par) -> numeric, # par = list(mu = <per obs>, extra...)
  cdf      = function(q, par) -> numeric, # F(q); must return 0 for q < 0, exact at integers
  rgen     = function(par) -> integer vector (one draw per obs)
)
```

`X` is the model matrix from `design` (built once). `off` is the offset row of the gene.
All `glmmTMB` calls use `control = glmmTMBControl(parallel = 1)` (parallelism comes from
`mclapply`).

### 5.1 Poisson
`glmmTMB(y ~ 0 + X + offset(off), family = poisson())`, or `glm.fit` (same MLE; check
agreement in G3). `cdf = ppois`, `rgen = rpois`.

### 5.2 NB
`family = nbinom2()`. `size = sigma(fit)`, `V = μ + μ²/size`. The issue's `φ` is
`1/size`; the output reports **both** under explicit names (`nb_size`, `nb_phi_issue`).
`cdf = pnbinom(q, mu, size)`, `rgen = rnbinom`.

### 5.3 ZINB
`family = nbinom2(), ziformula = ~1`. π is constant per gene and does **not** scale with
library size. Document this as a modelling choice.
`cdf(q) = π + (1−π)·pnbinom(q, …)` for `q ≥ 0`, 0 for `q < 0`.
`rgen`: `ifelse(runif < π, 0, rnbinom(…))`. π̂ at the boundary (≈ 0) is a legitimate
outcome, not a failure.

### 5.4 Generalised Poisson
`family = genpois()`. Own pmf with the mapping verified in §1:
`φ = predict(fit, type = "disp")`, `α = 1 − 1/√φ`, `λ1 = μ(1−α)`, `λ2 = α`,
`log p(y) = log λ1 + (y−1) log(λ1 + yλ2) − (λ1 + yλ2) − lgamma(y+1)`.
`cdf` = cumulative sum of the pmf over `0..q`.
**Underdispersion (φ < 1 ⇒ λ2 < 0):** terms with `λ1 + yλ2 ≤ 0` are zero, so the pmf need
not sum to 1. Compute the total mass `S_g` per gene and observation and store it. If
`S < 1 − 1e-8`, flag the gene. Do **not** renormalise silently; the flag count goes into
the report.
`rgen`: inversion from the cdf (or `glmmTMB::rgenpois` if the installed version exports
it; G2 must pass either way).

### 5.5 Poisson-lognormal (own fitter, `R/pln_core.R`)
Model: `Y | Z ~ Pois(exp(o + xβ + σZ))`, `Z ~ N(0,1)`. The regression coefficients are then
on the median scale; `E[Y] = exp(o + xβ + σ²/2)`. The offset still scales the rate
multiplicatively.
- pmf: **adaptive** Gauss–Hermite. Centre and scale the nodes on the mode of the integrand
  in `z` for each `(y, η, σ)`; default 30 nodes. Plain GH is inaccurate when σ is small and
  y is large, because the integrand is very sharply peaked.
- Fit: `nlminb` on `(β, log σ)`, starting from the Poisson fit and σ = 0.3.
- cdf: `Σ_{k} w_k · ppois(q, exp(η + σ z_k))`, using the same adapted nodes (see the
  G1 tolerance). `rgen`: `rpois(exp(η + σ·rnorm))`.
- Cross-check (reported, not a gate): σ̂ from `glmmTMB(y ~ … + (1|obs), poisson())`
  (Laplace approximation) on 200 genes.

### 5.6 Conway–Maxwell–Poisson (own fitter, `R/cmp_core.R`)
Mean-parameterised, as in `glmmTMB::compois` (Huang 2017). **Not** the
`COMPoissonReg`/Sellers–Shmueli λ-parameterisation: there an offset enters `log λ`, and
`E[Y] ≠ λ` unless ν = 1.
- `p(y) ∝ λ^y / (y!)^ν`, with `λ(μ, ν)` solved so that `E[Y] = μ`. Newton iterations on
  `log λ`, started from the asymptotic approximation `log λ ≈ ν·log(μ + (ν−1)/(2ν))`, with
  a bracketing fallback (`uniroot`).
- Normalising sums: log-sum-exp over a **window** `[max(0, m − K·s), m + K·s]` around the
  mode `m`, where `s ≈ sqrt(μ/ν)`. Default K = 40. Assert that the boundary terms are
  < 1e-15 of the maximum and widen the window otherwise. Cost is then O(s), not O(μ).
- Fit: `nlminb` on `(β, log ν)`, starting from the NB fit (`ν₀` chosen so the variance
  matches).
- Gate G3: on genes with `max(μ) ≤ 50`, logLik and estimates must agree with
  `glmmTMB::compois` (§9). `glmmTMB::compois` is used only for this cross-check.

---

## 6. Fitting a dataset (`R/fit_dataset.R`)

```r
fit_dataset(Y, X, O, fam, genes, n_cores) -> list(
  par   = <per gene: beta, extra, mu vector>,
  diag  = data.table(gene, loglik, conv_ok, conv_msg, engine, gp_mass_min, time_s),
  ok    = <logical per gene>
)
```

- `mclapply(genes, …, mc.cores = n_cores, mc.preschedule = TRUE)`. Every fit is wrapped in
  `tryCatch`; a failure yields `conv_ok = FALSE` with the message.
- **A fit counts as OK** when the optimiser converged (`convergence == 0`), `loglik` is
  finite, and all parameters are finite. `pdHess = FALSE` at a boundary (NB size → ∞,
  ZINB π → 0) is recorded, not treated as a failure.
- **Analysis gene set:** genes where **all six** models are OK on the observed data. Also
  report per model how many genes failed and why, and where the failures sit (their mean
  decile). As a sensitivity analysis, recompute the per-model observed statistics on each
  model's own OK set.
- In bootstrap replicates, genes whose refit fails are dropped from that replicate's
  statistic. Record the drop fraction per replicate; the report flags it if > 1%.

---

## 7. PIT and statistics

### 7.1 PIT (`R/pit.R`)
For every observation: `a = F(y−1)` and `b = F(y)`, where `F(−1) = 0`. Store `a` and `b`
(not only `u`). `u = runif(a, b)` with a recorded seed.
`rerandomise(a, b, R, seed)` redraws `u` R times (default R = 20) and reports the spread of
each statistic. The headline uses draw 1, stated as such.
Sanity asserts: `0 ≤ a ≤ b ≤ 1`, and `b − a` equal to `exp(logpmf(y))` to 1e-10 relative.

### 7.2 Statistics (`R/gof_stats.R`)
These formulas already agree with `goftest::cvm.test`, `goftest::ad.test` and
`stats::ks.test` to printed precision:

```r
cvm_u <- function(u) { u <- sort(u); N <- length(u)
  1/(12*N) + sum((u - (2*seq_len(N) - 1)/(2*N))^2) }
ad_u  <- function(u, eps = 1e-12) { u <- pmin(pmax(sort(u), eps), 1 - eps)
  N <- length(u); i <- seq_len(N)
  -N - mean((2*i - 1) * (log(u) + log1p(-rev(u)))) }
ks_u  <- function(u) { u <- sort(u); N <- length(u)
  max(seq_len(N)/N - u, u - (seq_len(N) - 1)/N) }
```

Report `W²`, `A²` and `D`, plus `W²/N` and `A²/N`. Under a fixed departure, `W²` and `A²`
grow roughly linearly with N, so the N-normalised forms are what can be compared across
datasets. Clamping in AD is needed because numerically `u` can equal 0 or 1 in the tails.
Record how many values were clamped.

### 7.3 Bootstrap reference
Per model: `excess = T_obs / median(T_b)` and `z = (T_obs − mean T_b) / sd(T_b)`.
These are **descriptive** cross-model comparisons. They correct for each model's own
in-sample estimation effect, which differs between models (Poisson has 0 extra
parameters per gene, ZINB has 2). They are not a formal model-selection criterion, and the
README must say so.

---

## 8. Parametric bootstrap, held-out score, strata, plots

### 8.1 Bootstrap (`R/bootstrap.R`)
For a model `m` fitted on `G_BOOT`:

```
for b in 1..B (parallel over b; genes serial inside each worker, or the reverse — pick one and document):
   Y_b   <- rgen(par_obs) for every gene in G_BOOT            # integer data
   O_b   <- OFFSETS_BOOT == "fixed" ? O : recompute offsets on Y_b (same OFFSET_MODE)
   fit_b <- fit_dataset(Y_b, X, O_b, fam, G_BOOT)             # identical procedure
   u_b   <- PIT(Y_b, fit_b)                                   # identical PIT code
   store: T_b = (W², A², D), per-stratum T_b (§8.3), 50-bin PIT histogram counts,
          per-sample mean of qnorm(clamp(u_b)), drop fraction
p_GOF = (1 + #{T_b ≥ T_obs}) / (B + 1)    # separately for W², A², D
```

- RNG: `RNGkind("L'Ecuyer-CMRG")`. Derive the stream of each `(model, b)` deterministically
  from `BASE_SEED`, so a single replicate can be rerun on its own.
- No re-filtering of simulated genes (the gene set stays fixed).
- `B_DEV = 100`, `B_FINAL = 500` (1000 if time allows), as in the issue.
- **Known limit (write it into the README and the table caption):** genes are simulated
  independently given their fitted parameters. Gene–gene correlation and shared
  sample-level effects are not in the null. If they are present in the real data,
  `p_GOF` is anti-conservative. At RNA-seq N, a floor `p_GOF = 1/(B+1)` for all models is
  the expected outcome, not a finding.

### 8.2 Held-out predictive score (`R/heldout.R`)
- K-fold over samples (default K = 8), stratified by the design factor so that every level
  stays in training. For `n < 16`, use leave-one-sample-out.
- Score: `mean log p(y_gi | θ̂_{−fold(i)})` per observation, per gene and model. The
  held-out sample's offset comes from the full-data TMM; note this small leakage in the
  README.
- Report per model: total and per-observation score, and the paired difference to NB per
  gene (mean, median, fraction of genes > 0). An interval from resampling genes is labelled
  approximate, because genes are not independent.
- By-product: the out-of-sample PITs (`a`, `b`), from which W², A² and D are computed as a
  second calibration view.
- Runs on `G_BOOT`.

### 8.3 Strata (`R/strata.R`)
Each rule is a function `(Y, O) -> factor per gene`. It is applied to the observed data
**and re-applied to every bootstrap dataset**. Stratifying by empirical dispersion or zero
fraction is selection on the same data that produce the PIT; re-applying the rule in each
bootstrap dataset is what calibrates that.

| rule | definition |
|---|---|
| `mean` | tertiles (and deciles for plots) of `log2` mean CPM using the offsets |
| `disp` | tertiles of edgeR tagwise dispersion (`estimateDisp` on the same DGEList) |
| `zero` | fraction of zeros: `0`, `(0, 0.1]`, `> 0.1` (the report shows the actual bin counts) |
| `condition` | observations split by design level (per observation, not per gene) |

For each stratum, compute `T_obs,s`, `p_GOF,s` against `T_b,s`, and the excess.

Per-sample diagnostic: the mean of `qnorm(clamp(u))` per sample, drawn against its
bootstrap envelope (2.5–97.5%). A sample outside the envelope across models points to a
sample-level effect.

### 8.4 Plots (`R/plots.R`, called only from `03_report.R`)
For each model:
- **PIT histogram** (50 bins), with a bootstrap 2.5–97.5% band per bin from the stored
  replicate histograms (calibrated, unlike a binomial band);
- **uniform QQ plot**, thinned to 10⁴ quantiles;
- **ECDF vs uniform CDF**, plus the difference plot `F̂(u) − u` with a bootstrap envelope;
- a mean-stratum panel for each.

Headless fonts: reuse `ensure_fonts()` from `null_calibration.R`.

### 8.5 Summary table (`03_report.R`)
One row per model:

`W², A², D` (on `G_BOOT`) · `W²/N` · `A²/N` · `p_GOF` (for each of W², A², D) ·
`excess(W²)` · held-out score/obs · Δ vs NB · `rank_heldout` · `rank_excess_W2` ·
failed fits · GP mass flags.

**Two rank columns, not one composite.** The issue asks for a single "Rank". Choosing
which criterion it uses is decision D5 (§12); the code shows both.

---

## 9. Validation gates (hard stops)

| gate | script | check | pass criterion |
|---|---|---|---|
| **G1** pmf/cdf | `00_validate_families.R` | on a parameter grid (μ ∈ {0.5, 5, 50, 500, 5·10³, 5·10⁴}, three dispersion levels per family): `Σ pmf` over the support; `cdf(q) − cdf(q−1) = pmf(q)` | `|Σ pmf − 1| < 1e-8` (except GP underdispersed: report the mass instead); increments agree to 1e-10 relative |
| **G2** rgen | same | 10⁵ draws per grid point, compared with the own cdf | max |ECDF − cdf| < 3·sqrt(1/10⁵) ≈ 0.0095; mean and variance within 4 standard errors |
| **G3** logLik consistency | same | own `logpmf` summed at the fitted parameters vs the engine's `logLik`; for CMP and PLN, own fitter vs `glmmTMB` on 200 simulated genes with μ ≤ 50 | glmmTMB families: |Δ| < 1e-6. CMP: |Δ logLik| < 1e-4 and relative Δ in β, ν < 1e-3. PLN: **reported only** (Laplace differs from quadrature) |
| **G4** parameter recovery | same | simulate 200 genes per family at n = the target dataset's n, then fit | median relative bias of the dispersion parameter reported; no gate threshold, but printed prominently |
| **G5** test calibration | `01_calibration_sim.R` | truth = NB, then truth = PLN; `M` datasets (default 40) with parameters drawn from the real-data NB/PLN fits; G = 300, n = target n, B = 99; run the full bootstrap for the **true** model and for Poisson | true-model `p_GOF` not concentrated near 0 (report its ECDF and the fraction < 0.05 with a binomial CI); Poisson rejected in ≥ 95% of datasets |
| **G6** PIT sanity | inside `02_run_gof.R` | the asserts of §7.1 | all pass |
| **G7** determinism | test | same seed ⇒ identical `T_obs` and `T_b[1:5]` | exact equality |

Unit tests (`testthat`) cover `gof_stats` against `goftest`/`ks.test`, the PIT bounds for
y = 0, `cmp_core` λ-solving (`E[Y] = μ` to 1e-10 relative), and the PLN quadrature against
`integrate()` on a grid.

The earlier pilot (sketch of 16.09, simulated NB, 300 genes) showed why G5 matters. With
**3 replicates per group, the true NB was rejected by the textbook CvM p-value
(p = 0.0011)**, with PIT sd 0.301 vs 0.289 for a uniform, i.e. excess mass at 0 and 1. That
is one simulation seed; G5 is its proper, replicated form.

---

## 10. Milestones and compute

| # | deliverable | done when |
|---|---|---|
| M0 | skeleton, config, `.gitignore`, provenance writer | `02_run_gof.R --dry-run` prints the resolved config |
| M1 | `families.R`, `cmp_core.R`, `pln_core.R` + unit tests | G1–G4 pass, report written |
| M2 | `pit.R`, `gof_stats.R`, `fit_dataset.R` | G6, G7 and unit tests pass; observed run on 500 genes of one dataset |
| M3 | `bootstrap.R` with `B_DEV = 100` on 500 genes | outputs complete; timing per model per replicate logged |
| M4 | `01_calibration_sim.R` | G5 report written |
| M5 | `heldout.R`, `strata.R` | outputs on the dev subset |
| M6 | `plots.R`, `report.R` | table + all plots produced from saved files, with no refitting |
| M7 | README + patch | Aaron reviews; **no commit** |

The full runs (`B_FINAL`, all datasets) are started by Aaron, not by Claude Code.

**Budget — measure first.** Use the per-fit times from M3 on the analysis machine. For
reference only, from the container timings in §1 (1 core, 0.16 s per fit, 5 non-CMP
models): `G_BOOT = 2000`, `B = 500` gives ≈ 2000 × 500 × 5 × 0.16 s ≈ 222 core-hours,
i.e. ≈ 7 h on 32 cores, **excluding CMP and the held-out score**. The held-out score adds
roughly K = 8 fits per gene and model, but no B-fold multiplier.

---

## 11. Data prerequisites (outside this code)

1. **Primary (proposed): yeast 48 × 48 WT vs Δsnf2, PRJEB5348**, quantified twice from the
   same FASTQs:
   - Salmon (+ `tx2gene`) → `build_from_tximport`;
   - STAR or HISAT2 + featureCounts → `build_from_featurecounts`.

   These two datasets differ only in quantification, which is the integer control the issue
   asks for. With 48 replicates per condition, the small-n estimation effect of §9 is
   minimal. Schurch et al. (*RNA* 2016) used 42 of 48 replicates per condition as
   high-quality, so a replicate-QC decision is needed (D6). The same design has been used
   for a distribution comparison before (Gierliński et al., *Bioinformatics* 2015), which
   gives a reference point for the result.
2. **Secondary:** TCGA cohorts via `build_from_tcga`. They are already on BioNAS, so they
   are the fastest dev data, but they are heterogeneous (§3).
3. Technical-replicate sets (MAQC/Quartet) answer a different question (technical noise
   only). If used, the dataset type goes into the table.

---

## 12. Decisions for Aaron / Asis (the code exposes each one as a config value; none is defaulted silently)

| id | decision | config | proposed |
|---|---|---|---|
| D1 | fractional counts | `ROUNDING` | `"stochastic"` primary, `"floor"` sensitivity. Needs Asis's OK, since the issue says not to round without justification |
| D2 | offset | `OFFSET_MODE` | `"tmm_length"` for tximport input, `"tmm"` for featureCounts/TCGA |
| D3 | offsets in the bootstrap | `OFFSETS_BOOT` | `"reestimate"` (the same procedure on simulated data) |
| D4 | bootstrap gene subset | `N_GENES_BOOT` | 2000, stratified by mean; revisit after M3 timing |
| D5 | the single "Rank" column | `RANK_BY` | held-out score (fair across parameter counts); excess-W² shown alongside |
| D6 | yeast replicate QC | `EXCLUDE_SAMPLES` | decide from QC before the run, **not** from GOF results |
| D7 | `glmmTMB` version | — | whatever is installed; record it. CMP runs on the own fitter regardless |

---

## 13. Out of scope

- No changes to the Tensor-Omics Fortran, C ABI, Rcpp or R wrappers.
- No shrinkage of dispersions. Per-gene MLE is used, as the issue specifies per-gene
  fitting. Shrinkage would be a different estimand; mention it in the README as a
  possible follow-up, not implemented.
- No conclusions about RNA-seq noise in any file until the full runs exist.
