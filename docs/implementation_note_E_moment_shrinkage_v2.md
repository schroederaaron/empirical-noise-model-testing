# Implementation note E (v2): distribution-free (moment-based) variance shrinkage, James-Stein type

**For:** Claude Code. **Status:** revised DRAFT for discussion (owner and Asis). **Start with this note's phase 1 (James-Stein shrinkage, level 2).**
**Scope of the first round (owner decision 27/10, evening): the cross-fitted moment-estimation shrinkage ONLY. The group-variance-ratio correction (section 2.5) and the level-2-plus-band arm are PARKED as ideas (phase 2). The shrinkage target is a switch, global median or LOESS trend, and both are run (section 2.3).**
**Builds on:** notes A (variance band), B (shrunk band, `variance_mode = 2`) and D (Part RS). **Repos:** Tensor-Omics (branch `145-empirical-noise-model`, module `noise_model_variance`) and `empirical-noise-model-testing`.
**Written without executing any Fortran or R.** Python toy simulations were run (section 7); everything marked **VERIFY** must be checked before relying on it.

---

## 0. What changed relative to v1 (review of 27/10/2026)

| # | v1 said | Finding | v2 |
|---|---|---|---|
| 1 | `V_g = n/(n-1)^3 * sum (w - wbar)^2` is an unbiased estimate of Var(s2) | **False at small n.** For normal noise E[V]/true = (n-2)/(n-1): 0.50 at n = 3, 0.75 at n = 5, 0.89 at n = 10; 0.43 (skewed) and 0.34 (heavy tail, t5) at n = 3. **At n = 3 it equals `s2^2/4` exactly** (verified to 1e-11): for 3 points `sum dev^4 = (sum dev^2)^2 / 2`, so a per-group data has no information on the fourth moment. | Replaced by the **cross-group contrast** estimator (section 2.2): unbiased for normal, skewed and heavy-tailed noise at n = 3, 5, 10 (ratio 0.997-1.006). |
| 2 | Level 2 scales each neighbour's residuals by its own pooled shrunk variance | **Not calibrated.** The residuals are internally studentised (part of their own denominator), so the null tail is too thin and the observed (externally studentised) t is compared against it. Toy: FDR 0.33 at n = 3 (nominal 0.05), null `p <= 0.001` rate 0.025 (25x nominal). The ideal lambda does not repair it. | **Cross-fitted scaling** (section 3): the case pool is scaled by the CONTROL-side shrunk variance and the control pool by the CASE-side one. Toy: FDR 0.06 / 0.04 / 0.04 at n = 3 / 5 / 10, tail rate 0.0016 / 0.0013 / 0.0014, ranking equal to a limma-like moderated t. |
| 3 | "Not needed from n around 7" | Not supported. With realistic heterogeneity (log-variance scatter about 1.05) the toy still gains at n = 10 (FDR 0.10 -> 0.04, TPR 0.60 -> 0.76), and the variance arms gained +0.05 to +0.11 AUPRC in Part R (n = 13-41). | **No expected threshold is assumed**; it comes from the tests (section 6). |
| 4 | Level 1 (selection only) and level 2 in the first implementation | Owner decision 27/10: start with the James-Stein shrinkage. | **Phase 1 = level 2 (cross-fitted moment shrinkage only); group-ratio correction and band arm parked.** Level 1 and per-bin lambda are phase 2 (section 8). |

**Verification log (what was and was not checked):**
code facts (module `noise_model_variance`, `compute_logvar_helper` returning pooled `s2`, `VARIANCE_MODE_BAND = 1` and `VARIANCE_MODE_SHRUNK = 2` so mode 3 is free, `loess_alloc(x, y, span, degree, fitted_values, mode, n_iters, ierr)` takes any `y`, `compute_pvalue` in `noise_model_exact`, `VAR_FLOOR = 1e-8`, RS registered in `power_test.R`): **confirmed**.
NOT verified (no R, no search access): the RS numbers quoted in section 1; the exact formulas in the cited papers and in `corpcor::var.shrink`; `CLAUDE.md` (not present in either repository as cloned; the documentation instruction refers to a local convention). Citations in this note are from memory and may be wrong: **check every one.**
Additional toy checks run for this version (Python, normal noise unless stated): unequal group variances (`toy_unequal_group_variance.py`, `toy_unequal_group_variance2.py`), global vs LOESS shrinkage target (`toy_const_target.py`, `toy_const_target2.py`). Not run: anything in Fortran or R.

---

## 1. Why this note exists

**Findings (Part RS, run of 27/10/2026; six TCGA base cohorts, 100 rounds per cell; as reported by the owner, not re-checked):**

* Null arm (pi1 = 0), share of null p <= 0.05: `TOX-log|TPM-centred` 0.054-0.067 at n = 3 and 0.047-0.063 at n = 5, i.e. close to nominal. limma 0.040-0.050.
* Power arm (pi1 = 0.1), observed FDR at nominal 0.05: `TOX-log|TPM-centred` 0.22-0.53 at n = 3 and 0.14-0.31 at n = 5. Variance band c3: 0.16-0.37 and 0.07-0.13. DESeq2: 0.20-0.43 at n = 3. limma: 0.005-0.071.
* The shrunk-selection variant (note B) does not change this at n = 3. A wider band helped more than shrinkage did.
* At n = 3 the variance-band arms rank genes nearly as well as limma, edgeR and DESeq2. The problem is the tail calibration of the p-values, not the ranking.

**Why a calibrated null can still give uncontrolled FDR.** `p <= 0.05` is a marginal statement, averaged over genes. BH works in the tail. The toy shows exactly the observed pattern: current TOX has null `p <= 0.05` at 0.05 but `p <= 0.001` at 0.0025-0.0070 (2.5x-7x nominal) and FDR 0.33-0.67 at n = 3. Noisy genes sit in the extreme tail of a pool built mostly from quieter genes.

**Diagnosis (hypothesis, supported by the toys, to be tested on real data):**

1. **The mean-matched model never uses the gene's own variance.** The noise scale comes entirely from the local trend: full pooling, the extreme of shrinkage (lambda = 1). A gene whose true noise is above its neighbours' is tested against too narrow a null.
2. **Note B shrinks only the selection, never the test.** limma, edgeR and DESeq2 shrink the variance (or dispersion) that enters the test statistic. That is what makes them work at small n.
3. **Note B is not distribution-free** (`d0` from trigamma moment matching assumes log-chi-square scatter). The estimator below uses only moments.

**Owner constraints (20-27/10/2026):**

* No distributional assumption for the data, the variances or the null. No externally pre-computed weights or parametric prior: every weight is estimated from the data inside the model.
* The shrinkage target is data-driven, never an assumed form: the median of the gene-wise variances (as in the original James-Stein / Opgen-Rhein and Strimmer shrinkage) or the LOESS trend of the variances on expression. Both are run in the first round (section 2.3).
* Shrinkage is a user switch (on/off). From which n it is no longer needed is decided by the tests.
* Goal: competitive with edgeR, limma and DESeq2 under most conditions, not to beat them. If the shrinkage does not help, it was a cheap test.

---

## 2. The estimator: moment-based shrinkage of gene-wise variances

James-Stein-type shrinkage with a data-driven intensity, as in the "shrinkage t" of Opgen-Rhein and Strimmer (2007), with two target options (global median as in the paper, or our LOESS trend) and with a **corrected, distribution-free estimate of the sampling variance** (section 2.2).

**References (from memory; VERIFY before citing):** James W, Stein C (1961), Proc. 4th Berkeley Symp. (the original shrinkage idea). Opgen-Rhein R, Strimmer K (2007), *Statistical Applications in Genetics and Molecular Biology* 6:9 (shrinkage t). Schaefer J, Strimmer K (2005), *SAGMB* 4:32 (analytic shrinkage intensity). Reference implementation for a cross-check: `corpcor::var.shrink` (**VERIFY** its formula, target and lambda clipping). Related work for the owner's reading: Jain N et al. (2003), local pooled error (LPE) test, *Bioinformatics* 19:1945; Baldi P, Long AD (2001), Cyber-T, *Bioinformatics* 17:509; Cui X et al. (2005), James-Stein variance shrinkage with permutation, *Biostatistics* 6:59.

### 2.1 Per-gene quantities

On the scale the model tests on: `y = log2(x + 1)` by default (`norm_method = 1`); the linear branch uses `y = x`. For gene g and group k in {c (case), t (control)} with `n_k` replicates:

```
s2_k,g  = sum_i (y_kgi - ybar_kg)^2 / (n_k - 1)                 unbiased group variance (df n_k - 1)
d       = n_c + n_t - 2
s2_g    = ((n_c - 1) s2_c,g + (n_t - 1) s2_t,g) / d              pooled variance (exactly `s2` of compute_logvar_helper)
delta_g = s2_c,g - s2_t,g                                         the cross-group contrast
```

### 2.2 Sampling variance of the variance estimates (distribution-free)

Assumptions: the two groups are independent; the true variance is equal in both groups (the same assumption the pooled `s2_g` already makes); fourth moments are finite. **No distribution is assumed.**
Then `E[s2_c] = E[s2_t]`, hence `E[delta_g^2] = Var(s2_c) + Var(s2_t)`, and `delta_g^2` is an unbiased estimator of that sum for every gene. Splitting it by group size:

```
Vhat_c,g    = ((n_t - 1) / d) * delta_g^2           estimate of Var(s2_c,g)
Vhat_t,g    = ((n_c - 1) / d) * delta_g^2           estimate of Var(s2_t,g)
Vhat_pool,g = ((n_c - 1)(n_t - 1) / d^2) * delta_g^2   estimate of Var(s2_g)          balanced design: delta^2 / 4
```

* **Balanced design (n_c = n_t): exact.** Unbalanced: the split between the groups assumes `Var(s2_k)` proportional to `1/(n_k - 1)`, which is exact for normal data and approximate otherwise; report the unbalanced case separately (**VERIFY** by simulation).
* Each `Vhat_g` is a single noisy contrast; only its SUM over genes is used (section 2.4).
* A gene with a true variance difference between the groups (e.g. a partly responding DE gene) inflates its own `delta_g^2`; the effect on the global sums is small when such genes are rare, and the direction is more shrinkage.
* **Do not use** the v1 formula `n/(n-1)^3 sum (w - wbar)^2` (biased low, exactly `s2^2/4` at n = 3).

### 2.3 Shrinkage target: global median (the paper's choice) or LOESS trend

```
shrink_target = "median":  T = median over tested genes of s2_g                      (one constant; robust)
shrink_target = "loess":   T(m) = LOESS of s2_g on m_g,   m_g = log2(mean of all samples of g + 1)   (same trend axis as notes A/B)
```

* Both are fit on the **variance scale**, not on `log s2` (no Jensen correction, no distribution needed). **Never use the global MEAN of `s2` as the target** (dominated by hypervariable genes; toy FDR 0.085 / 0.088 / 0.066 at n = 3 / 5 / 10 vs 0.047 / 0.049 / 0.035 for the median).
* `"loess"`: default non-robust (degree 1, span 0.3); `loess_robust` is the existing switch. Toy: robust and non-robust give the same results within noise (robust slightly more conservative in the tail). **VERIFY** that `tox_loess::loess_alloc` can fit `s2` with the same settings used for `l` (the signature takes any `y`).
* The same `T` is the target for the pooled and both group-wise shrinkages (`E[s2_c] = E[s2_t] = E[s2_g]`). Floor `T >= VAR_FLOOR`.
* **Evidence from the toy** (`toy_const_target.py`, `toy_const_target2.py`; normal noise, variance falls about 28x from low to high expression, sigma_v = 1.05, 10 % DE, 5-6 datasets per cell; differences below about 0.03 are noise). Median target vs LOESS target, n = 3 / 5 / 10:
  AUPRC excess 0.620 / 0.718 / 0.800 vs 0.631 / 0.728 / 0.808; FDR 0.047 / 0.049 / 0.035 vs 0.057 / 0.066 / 0.039; TPR 0.392 / 0.579 / 0.770 vs 0.394 / 0.600 / 0.789; lambda 0.37 / 0.23 / 0.09 vs 0.45 / 0.31 / 0.13.
  By expression tertile (low / mid / high) at n = 5, TPR 0.226 / 0.649 / 0.846 (median) vs 0.275 / 0.674 / 0.841 (LOESS); FDR in the low tertile 0.038 vs 0.128 (n = 3: 0.157 vs 0.251). So no loss at high expression, slightly less power but better calibration at low expression.
  Why the aggregate parity: neighbours are mean-matched, so the target and its neighbours are shifted by the same constant; what matters is the variance of a gene relative to its neighbours. With a constant target the variance trend counts as scatter, so lambda is smaller (less shrinkage).
* **Caveat:** real data may have steeper or irregular variance trends and zeros at high expression. Report TPR and FDR **by expression decile** for both targets, especially the top deciles where Part R showed the largest TOX deficit.

### 2.4 Shrinkage intensities and shrunk variances (global, three scalars)

(In the formulas below `T(m_g)` is the constant `T` when `shrink_target = "median"`.)

```
lambda_pool = clip( sum_g Vhat_pool,g / sum_g (s2_g   - T(m_g))^2 , 0, 1 )
lambda_c    = clip( sum_g Vhat_c,g    / sum_g (s2_c,g - T(m_g))^2 , 0, 1 )
lambda_t    = clip( sum_g Vhat_t,g    / sum_g (s2_t,g - T(m_g))^2 , 0, 1 )

s2s_pool,g  = lambda_pool * T(m_g) + (1 - lambda_pool) * s2_g
s2s_c,g     = lambda_c    * T(m_g) + (1 - lambda_c)    * s2_c,g
s2s_t,g     = lambda_t    * T(m_g) + (1 - lambda_t)    * s2_t,g
```

* **Meaning:** numerator = how noisy the gene-wise variance estimates are (sampling error); denominator = how far they scatter around the trend in total. All scatter is noise -> lambda -> 1 (every gene pulled onto the trend, today's TOX). Genes truly differ -> lambda small.
* One lambda per variance type (global), as in Opgen-Rhein and Strimmer; each `Vhat_g` is very noisy but lambda pools thousands of genes. `T` is estimated from thousands of genes, so its own variance is negligible next to `Vhat_g` (state this in the methods text).
* Toy: with normal noise the estimate reproduces the empirical-Bayes intensity of a limma-like moderated t (0.39 vs 0.40 at n = 3, 0.28 vs 0.25 at n = 5, 0.14 vs 0.13 at n = 10, log-variance scatter 1.05). With skewed or heavy-tailed noise it shrinks MORE (0.51 and 0.70 at n = 3) because the sampling variance is genuinely larger: the moment estimator adapts to the distribution, a normal-theory prior does not.
* Report `lambda_pool`, `lambda_c`, `lambda_t`, the unclipped ratios, the median `Vhat_pool`, the share of genes with `s2s_pool` within 10 % of `T`, and `n_floored`. Per-decile lambdas as a report-only diagnostic.
* Guard: `s2s >= VAR_FLOOR` before dividing; report the number of genes floored.

### 2.5 PARKED (phase 2, NOT in the first implementation): trend-level group-variance ratio

*Status 27/10: owner decision, keep as an idea, leave out for now. Everything below is a specification for later, with the toy evidence that motivated it.*

The cross-fitted scaling puts the case pool on the control scale and vice versa. If the groups differ systematically in noise (variance ratio `rho^2`), the null variance becomes `(rho^2 + rho^-2)/n` instead of `2/n` (conservative, power loss; toy, section 7). A **trend-level** correction needs no distribution and no per-gene ratio (a per-gene ratio is hopelessly noisy at n = 3):

```
Tc(m), Tt(m) = LOESS of s2_c and of s2_t on m (variance scale, same settings as T)
r(m)         = clip( Tc(m) / Tt(m), 0.05, 20 )                          trend-level case/control variance ratio
delta'_g     = s2_c,g - r(m_g) * s2_t,g                                  (E[delta'] = 0 when the ratio is right)
Vhat_c = delta'^2 / 2 ;  Vhat_t = delta'^2 / (2 r^2) ;  Vhat_pool = (delta'^2 / 2) * (1 + 1/r^2) / 4          (balanced design)
targets      T_c = Tc(m), T_t = Tt(m), T_pool = (Tc + Tt)/2 ; lambdas and s2s_* as in 2.4 with their own targets
pool scaling case residuals / sqrt(s2s_t,j * r(m_j)) ;  control residuals / sqrt(s2s_c,j / r(m_j))               (both on unit scale)
observed     t_g = delta_beta / sqrt(s2s_pool,g)                          s2s_pool estimates (sigma_c^2 + sigma_t^2)/2
```

* With `r = 1` everywhere this reduces to 2.2-2.4 (parity test: forcing `r = 1` gives the uncorrected p-values). Balanced designs only when enabled (**VERIFY** the unbalanced generalisation before enabling it).
* **Diagnostic (when implemented; reported with and without `group_ratio`):** the trend `r(m)` per contrast, as `rho_hat(m) = sqrt(r(m))`; the owner can measure it on real contrasts with `group_variance_ratio_diagnostic.R`.
* Assumption that remains: the ratio is a smooth function of the mean (trend level). Gene-specific ratios are not corrected.

---

## 3. Level 2: studentised exact null with cross-fitted scaling (PHASE 1)

The shrunk variance enters the **test**, as in the reference tools, while the null keeps its empirical shape. Neighbours: the plain mean-matched neighbourhood (adaptive kNN of the exact model); a variance band is not needed (whether it still helps is a later arm).

```
observed   t_g    = (ybar_c,g - ybar_t,g) / sqrt(s2s_pool,g)               scale factor sqrt(1/n_c + 1/n_t) dropped on both sides
residuals  z_c,j,i = eps_c,j,i / sqrt(s2s_t,j)   for every neighbour j in the CASE pool      (opposite-group variance)
           z_t,j,i = eps_t,j,i / sqrt(s2s_c,j)   for every neighbour j in the CONTROL pool   (opposite-group variance)
null pair  v      = z_c / sqrt(n_c) - z_t / sqrt(n_t)                       (eps with the Bessel factor, as now)
p_g       = (1 + #{ |v| >= |t_g| }) / (N + 1)                              exact pairwise count, unchanged `compute_pvalue`
```

* **Why the opposite group:** under H0 with equal variance the scale applied to a gene's case residuals is computed from the control samples only, hence independent of those residuals. The v1 construction (scaling by the pooled variance, which contains the residual being scaled) makes the null internally studentised: its tails are too thin (toy: tail rate 25x nominal). Cross-fitting removes this without any distributional assumption.
* **Implementation is small:** per-gene, per-side scale factors `1/sqrt(s2s_opposite,j)` applied when residuals are packed into the pools (the pools already carry the per-side `1/sqrt(n)`); pass `|beta_g| / sqrt(s2s_pool,g)` as the observed statistic; `compute_pvalue` is reused. The adaptive growth rule (`tau`) then acts on scaled residuals: **VERIFY** that it still behaves sensibly (the scaled pool is roughly unit-scale).
* **Switch:** `shrink = FALSE` (default until the tests decide) gives exactly the current exact model. `shrink = TRUE` gives this model. Benchmark-only negative control: `shrink_scaling = "same_group"` (the v1 construction) to confirm on real data what the toy predicts; not a user option.
* **Limits (built-in sanity checks):**
  * `lambda = 1` and a flat trend: `s2s` is the same constant for every gene, the scaling cancels, p-values equal the current exact model **exactly**.
  * `lambda = 1` and a non-flat trend: close to the current model; neighbours differ only by `T(m_j)/T(m_g)`.
  * `lambda = 0`: per-gene studentisation, unstable at n = 3 (why lambda is estimated).
  * So lambda interpolates between full pooling (today's TOX) and gene-wise studentisation, and the data choose where.
* **Assumptions, stated for the manual:** independent samples; equal variance in the two groups under H0 (already assumed by the pooled `s2`); finite fourth moments; one global lambda per variance type. **No distribution for the data, the variances or the null.**

---

## 4. Fortran / R changes (1-D)

* `noise_model_variance` (shares `compute_logvar_helper`, the LOESS and the gather): add `compute_moment_shrinkage_helper` (inputs: replicate matrices, `m`, target type and trend settings; the LOESS is called only for `shrink_target = "loess"`, the median needs a sort (**VERIFY** whether the median helper used for `delta_hat` in `noise_model_exact` can be reused); outputs: `s2_c`, `s2_t`, `s2`, `T`, `Vhat_*`, the three lambdas, `s2s_c`, `s2s_t`, `s2s_pool`). Pure, no allocation, as the other helpers.
* Level 2: a `shrink` argument in the variance pipeline (mode 0 = no band, the plain mean-matched gather; **VERIFY** whether reusing `noise_model_exact::gather_residuals_helper` is cleaner than a band with `tau_var` huge), the per-side residual scaling in the pool, and the scaled observed statistic. Outputs per gene: `s2s_pool`; scalar outputs: the three lambdas, `n_floored`.
* C wrapper, Rcpp and R wrapper: `tox_compute_noise_pvalues_pipeline_variance(..., variance_mode = c("band", "shrunk", "none"), shrink = FALSE, shrink_target = c("median", "loess"), loess_robust = FALSE)`; `norm_method = 1` stays the default.
* Documentation per the project convention (current_state, eval doc (new section), labbook); methods manuscript Section 11.2.5 ("in development") is rewritten once the tests decide.

---

## 5. Tests

### Unit tests (Fortran, `test/mod_test_noise_model_variance.f90`)

1. **`Vhat` (contrast) is unbiased:** mean of `Vhat_pool` over many simulated genes equals the empirical variance of `s2` within MC error, for normal, skewed (log-gamma) and heavy-tailed (t5) noise, at n = 3, 5, 10 (balanced); report the unbalanced bias. Include a documentation test that the v1 formula equals `s2^2/4` at n = 3.
2. **lambda:** homogeneous true variances (all genes on the trend) -> lambda_pool close to 1; strongly heterogeneous at n = 10 -> clearly below 1.
3. **Cross-check against `corpcor::var.shrink`** with the trend replaced by a constant (median) target and the v1 `V`, on a saved fixture (if available; tolerance 1e-8). **VERIFY** the target conventions first. (This checks the shrinkage algebra, not the corrected `V`.)
4. **Parity:** Level 2 with a constant target and lambda forced to 1: p-values identical to `noise_model_exact` (bit-exact or 1e-15).
5. **Scale invariance (linear branch):** multiplying all data of one gene's neighbourhood by a constant does not change that gene's p-value.
6. **Calibration regression test (fixed seeds, synthetic normal data, n = 3, heterogeneous variances):** null `p <= 0.001` at most 2x nominal for the cross-fitted scaling, and clearly above that for `shrink_scaling = "same_group"` (guards against reintroducing the internal studentisation).
7. Floor and validation paths (constant genes, `n_k = 2`, `n_k < 2` error, unbalanced designs).

### Benchmark (testing repo, `power_test.R`)

New arms (labels proposed):

```
TOX-log-stud-med        engine "variance", variance_mode = "none", shrink = TRUE, shrink_target = "median"    (level 2, global median target)   FIRST ROUND
TOX-log-stud-loess      same, shrink_target = "loess" (non-robust)                                             (level 2, LOESS-trend target)     FIRST ROUND
TOX-log-stud-nocf       same as -med, shrink_scaling = "same_group"                                            (negative control, v1)            FIRST ROUND
TOX-log-stud-rob        same as -loess, loess_robust = TRUE                                                    optional, only if -loess is promising

PARKED (phase 2; not implemented in the first round):
TOX-log-stud-gr         group_ratio = TRUE (section 2.5)
TOX-log-stud-vband-c3   variance_mode = "band", c = 3, shrink = TRUE    (level 2 + band; owner idea, keep in mind)
```

* Parts: **RS** (primary), R, P+ (incl. pi1 = 0), V (sparse hypervariable genes: does shrinking towards the trend hide true outliers?).
* **Extend RS to n = 7 and n = 10.** All six base cohorts have S >= 26 (`power_rs_checks.csv`), so both fit; at n = 10 a round uses up to 77 % of the smallest cohort (LUAD Stage IV, S = 26) and rounds overlap strongly. This is what decides "from which n is shrinkage not needed": compare `shrink = TRUE` vs `FALSE` per n.
* **P+ extension (equal-variance assumption):** case SD = rho x control SD, rho in {1, 1.25, 1.5, 2}, `TOX-log-stud-med`, `TOX-log-stud-loess` vs `TOX-log|TPM-centred`, limma, edgeR, DESeq2, pi1 in {0, 0.1}, n in {3, 5, 10}.
* **Side arm for the p-value-floor / neighbourhood hypothesis:** `TOX-log|TPM-centred` with a larger neighbourhood (k_max 100 and 200) and a stricter tau, at n = 3 and 5. Report FPR at 0.001 and 0.01 on the null arm.
* Report additionally: TPR and FDR by expression decile for `-med` and `-loess`; the three lambdas per cohort and n; FPR at 0.001 / 0.01 / 0.05 on the null arm; observed FDR at nominal 0.05 and 0.1; TPR at nominal and achieved FDR; AUPRC excess; gap to limma, edgeR, DESeq2.

---

## 6. Decision criteria (PROPOSAL: owner confirms before the run)

Reuse note D section 6 (calibration gate first, then power), with these additions:

1. **FDR gate:** observed FDR at nominal 0.05 at most 0.075 in at least 5 of 6 cohorts, per n.
2. **Tail gate (new):** on the null arm `p <= 0.001` at most 0.002 and `p <= 0.01` at most 0.015 in at least 5 of 6 cohorts.
3. **Shrinkage helps at n** if, with `shrink = TRUE` vs `FALSE` (same arm otherwise): the FDR gate passes with shrinkage and not without, or AUPRC excess rises by at least 0.02 in at least 5 of 6 cohorts with no loss in calibration.
4. **"Not needed from n on":** the smallest n in {3, 5, 7, 10} from which `shrink = FALSE` passes the calibration, tail and FDR gates and `shrink = TRUE` gives no AUPRC gain of 0.02 or more. No value is assumed in advance.
5. **Default of the switch:** decided after 1-4. If shrinkage helps at small n and costs nothing at large n, `shrink = TRUE` becomes the default.
6. Report every arm, including failures.

---

## 7. Expected behaviour from the Python toys (sanity checks, NOT targets)

`toy_note_E_level2_test.py` (single expression bin) and `toy_note_E2_trend_loess_cross_fit.py` (mean-variance trend, LOESS on the variance scale, 100-gene mean-bin pools). Normal / skewed / heavy-tailed noise, 10 % DE, 3000 genes, 5 datasets per cell: differences below about 0.03 are noise. Log-variance scatter sigma_v = 1.05. Cells: FDR / TPR at BH 0.05 | null `p <= 0.001` (nominal 0.001).

| noise | n | TOX exact (current) | Level 2 as written (v1) | **Level 2 revised (non-robust)** | limma-like moderated t |
|---|---|---|---|---|---|
| normal | 3 | 0.116 / 0.251 \| 0.0031 | 0.334 / 0.664 \| 0.0250 | **0.060 / 0.355 \| 0.0016** | 0.027 / 0.357 \| 0.0007 |
| normal | 5 | 0.098 / 0.404 \| 0.0027 | 0.194 / 0.741 \| 0.0096 | **0.041 / 0.597 \| 0.0013** | 0.041 / 0.604 \| 0.0010 |
| normal | 10 | 0.101 / 0.599 \| 0.0032 | 0.103 / 0.813 \| 0.0039 | **0.044 / 0.762 \| 0.0014** | 0.037 / 0.761 \| 0.0011 |
| skewed | 5 | 0.087 / 0.361 \| 0.0030 | 0.173 / 0.727 \| 0.0083 | **0.023 / 0.534 \| 0.0006** | 0.037 / 0.585 \| 0.0010 |
| heavy (t5) | 5 | 0.053 / 0.300 \| 0.0012 | 0.103 / 0.717 \| 0.0040 | **0.008 / 0.415 \| 0.0001** | 0.030 / 0.617 \| 0.0007 |

* With no heterogeneity (sigma_v = 0) the revised model equals the current model (lambda about 0.94-1.0): FDR 0.054 vs 0.047 at n = 3, AUPRC 0.603 vs 0.605.
* AUPRC excess at sigma_v = 1.05 (normal): n = 3 0.551 -> 0.627 (limma-like 0.635); n = 5 0.658 -> 0.737 (0.734); n = 10 0.748 -> 0.806 (0.803).
* **Heavy-tailed noise: the revised model is calibrated but conservative and loses power relative to the limma-like reference** (TPR 0.26 vs 0.41 at n = 3, 0.42 vs 0.62 at n = 5): the empirical null correctly carries the heavy tails, a normal-theory reference does not. Do not read this as a defect before comparing on real data.
* n = 3 stays hard for every method; none of these numbers say that n = 3 will be easy.
* Limits of the toys: normal/skewed/heavy noise at the observation level with equal variance in both groups, no zero inflation, no sample structure; pools are 100-gene mean bins, not per-gene adaptive kNN; the DE model is simple. TCGA with a known effect added is a stress test, not an ideal study (cohort heterogeneity).

---

## 8. Phase 2 (after phase 1 is evaluated) and risks

**Phase 2 (not in the first implementation):**
* The parked ideas: trend-level group-variance ratio correction (section 2.5, arm `-gr`) and level 2 plus the variance band (arm `-vband-c3`); robust LOESS target (`-rob`) unless the first round shows a reason to run it earlier.
* Level 1: selection-only successor of note B, `u_g = log(s2s_pool,g / T(m_g))`, `tau_var = c * 1.4826 * MAD_g(u_g)` (empirical spread, no trigamma), new `variance_mode = 3` ("mshrunk"); expected effect modest.
* Lambda per expression bin (deciles of m) if the per-decile diagnostics show mean-dependent heterogeneity.
* Multi-dimensional model: level 2 carries over per axis (`lambda_a`, `T_a`, `s2s` per axis, cross-fitted pools, product null). The MD model currently uses the mean-matched neighbourhood only, so level 2 is its natural first variance-aware step.
* Unbalanced-design handling beyond the proportional split; per-group n different by gene (missing values).

**Risks and open points:**
* **Equal variance across groups (for the same gene, under H0) is assumed** in three places: pooled `s2`, `Vhat`, cross-fit scaling. Toy test (`toy_unequal_group_variance.py`, case SD = rho x control SD for every gene): the cross-fitted null then has variance `(rho^2 + rho^-2)/n` instead of `2/n` (2.1x too wide at rho = 2), so it becomes conservative and loses power (n = 10, rho = 2: TPR 0.316 vs 0.596 for a pooled-t reference, FDR 0.006); at n = 3, rho = 2 FDR is 0.229 (anti-conservative) but still better than the current model (0.339). A pooled t is robust to this in balanced designs. **Part R and RS split ONE cohort at random, so rho = 1 by construction: the benchmark cannot test this assumption.** Second toy (`toy_unequal_group_variance2.py`, rho in {1, 1.25, 1.5, 2}): at the mild ratio rho = 1.25 (variance x1.56) the uncorrected cross-fit stays calibrated (FDR 0.06 / 0.04 / 0.04 at n = 3 / 5 / 10) and loses little relative to a pooled-t reference; the loss grows from rho = 1.5. The trend-level correction of section 2.5 restores most of the power without a cost at rho = 1 (n = 10, rho = 2: TPR 0.587 vs 0.301 uncorrected vs 0.602 pooled t, FDR 0.067). A variance-band pool WITHOUT scaling (crude version, tertiles of detrended pooled variance) does not help with the group ratio and is anti-conservative from the plug-in effect (FDR 0.41 / 0.22 / 0.12 at rho = 1); the group ratio is a systematic scale difference between the two pools, which selecting neighbours by pooled variance does not change. Follow-ups: (a) P+ extension with case SD = rho x control SD, rho in {1, 1.25, 1.5, 2}; (b) measure the real ratio with `group_variance_ratio_diagnostic.R`; (c) if the real ratios are not close to 1, un-park the group-ratio correction (section 2.5).
* **Constant (median) target:** shrinks genes at high expression (true variance far below the global median) towards a value that is too large. The toy shows no loss, but real data may differ; decided by the by-decile report (section 5).
* **Hypervariable genes:** shrinking towards the trend pulls genuinely noisy genes in. Part V answers whether that costs FDR; a global lambda may be too strong or too weak for them (per-bin lambda as the fallback).
* **Heavy tails:** conservative empirical null, power cost relative to a parametric reference (section 7).
* **The FDR problem may lie elsewhere.** If level 2 does not fix FDR at n = 3 on real data, candidates are the tail shape of the null (pool size, discreteness of log2(x + 1) at low counts) and the low-count genes. The side arm in section 5 tests the pool-size part.
* **Remaining parametric pieces elsewhere:** note A's `tau_var = c * sqrt(psi1(d/2))` assumes log-chi-square sampling; decide whether mode A should switch to an empirical width as well.
* **Out of scope here:** MD implementation, bootstrap and blocked engines, per-bin lambda before the global version has been tested.

---

## 9. Reading guide for the three papers (for the owner)

* **James-Stein / Opgen-Rhein and Strimmer:** how the shrinkage intensity is derived analytically (a ratio of sampling variance to total scatter around the target) and what they use as the target and for clipping; which estimator of the sampling variance they use and under which conditions they call it unbiased (compare with the contrast estimator in section 2.2); how they obtain significance (as far as I recall they rank with the shrinkage t and calibrate with an empirical null, check).
* **Local pooled error (LPE) test:** how the local pooled variance is defined (window over intensity, degrees of freedom), how p-values are obtained at 2-3 replicates, and how it compares with TOX's mean-matched neighbourhood (the closest relative).
* **Cyber-T:** how the prior degrees of freedom and the window size are chosen and how the prior variance is combined with the gene's own variance (compare with lambda and with limma's `d0 / (d0 + d)`).
