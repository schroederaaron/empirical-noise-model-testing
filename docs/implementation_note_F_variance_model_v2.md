# Implementation note F: variance-aware noise model v2 (James-Stein shrinkage, cross-fitting, adaptive growth)

**For:** Claude Code. **Status:** DRAFT for implementation after owner review. Written without executing any Fortran or R; Python toy checks were run (section 8).
**Repos:** Tensor-Omics, branch `145-empirical-noise-model`, module `noise_model_variance` (Fortran, C/Rcpp/R wrappers); `empirical-noise-model-testing` (benchmark).
**Companion:** note G (tests on the data structure: cross-fitting on real data, variance-ratio diagnostic).
**Supersedes:** note B (`variance_mode = 2`, "shrunk", prior df `d0` from trigamma) and note E v2 sections 4-5 (implementation and arm list). Note E v2 remains the reference for the estimator's derivation and the toy evidence.

---

## 0. Decisions this note implements (owner and Asis, meeting late October 2026)

1. Evaluate cross-fitting on real data (Asis had a comparable idea).
2. Replace the current shrinkage (note B: parametric `d0` from trigamma moment matching) by **James-Stein (moment-based, distribution-free) shrinkage**. Shrinkage stays a switch that can be turned off.
3. Compare three variance-aware approaches: **(a) the LOESS variance band as it is now**, **(b) cross-fitting without shrinkage**, **(c) cross-fitting with James-Stein shrinkage**, across replicate numbers and datasets.
4. Add an **adaptive, variance-aware neighbourhood growth** (Asis's "growth radius" idea) as a fourth neighbourhood option.
5. Asis's caution: adapting ever closer to limma/edgeR/DESeq2 risks losing TOX's standalone character. What stays distinct is listed in section 9.

**Assumptions made in this note (owner please confirm or correct):**
* "Cross-fitting without shrinkage" = the cross-fitted scaling of note E v2 with lambda = 0 (raw group variances), NOT cross-fitted neighbour *selection*. **Confirmed by the owner.**
* "Band as it is now" = the current band arm with the best calibration in RS, `TOX-log-vband-c3` (raw detrended log-variance, `tau_var = c * sqrt(psi1(d/2))`).
* The adaptive growth specification in section 5 is a proposal; Asis's exact intention is not documented.

---

## 1. Design: three orthogonal switches

| Switch | Values | Meaning |
|---|---|---|
| `variance_mode` (neighbour selection) | `"mean"` (plain mean-matched kNN of the exact model), `"band"` (note A, unchanged), `"adaptive"` (new, section 5) | which genes enter the pools |
| `match_variable` (for band/adaptive) | `"raw"` (note A, unchanged), `"js"` (James-Stein-shrunk, replaces note B) | the variance variable used to accept neighbours |
| `scaling` (null construction) | `"none"` (current), `"xfit_raw"` (cross-fitted, lambda = 0), `"xfit_js"` (cross-fitted, James-Stein) | how pool residuals and the observed statistic are scaled |

Plus `shrink_target` (`"median"` default, `"loess"`), `lambda_scale` (`"relative"` default, `"absolute"`), `loess_robust` (existing). The parametric note-B code path (`d0`, trigamma, inverse trigamma) is **removed** (alpha stage, no aliases); its old results stay in earlier runs for comparison.

---

## 2. James-Stein (moment) shrinkage of gene-wise variances

Scale `y = log2(x + 1)` (`norm_method = 1`) or `y = x` (linear branch). Groups c (case) and t (control), `n_c`, `n_t` replicates, `d = n_c + n_t - 2`.

```
s2_c,g, s2_t,g        unbiased group variances;   s2_g = ((n_c-1) s2_c,g + (n_t-1) s2_t,g) / d   (pooled; = s2 of compute_logvar_helper)
delta_g = s2_c,g - s2_t,g                          cross-group contrast
Vhat_c,g    = ((n_t-1)/d) * delta_g^2              unbiased for Var(s2_c,g)      (balanced: exact; unbalanced: approximate, see E v2 2.2)
Vhat_t,g    = ((n_c-1)/d) * delta_g^2
Vhat_pool,g = ((n_c-1)(n_t-1)/d^2) * delta_g^2     unbiased for Var(s2_g)        (balanced: delta^2 / 4)
T_g         = median_g(s2_g)  ("median", default)  |  LOESS of s2_g on m_g = log2(mean_g + 1) ("loess"); variance scale; floor VAR_FLOOR
lambda      = clip( sum_g w_g * Vhat_g / sum_g w_g * (s2_g - T_g)^2 , 0, 1 )      w_g = 1/T_g^2 ("relative", default) or 1 ("absolute")
s2s_g       = lambda * T_g + (1 - lambda) * s2_g                                    (pooled; same form with lambda_c, lambda_t for the groups)
```

* **Why the contrast estimator:** the textbook per-group estimator `n/(n-1)^3 sum (w - wbar)^2` is biased low at small n and at n = 3 equals `s2^2/4` exactly (no fourth-moment information). The contrast `delta_g^2` is unbiased for normal, skewed and heavy-tailed noise (verified, ratio 0.997-1.006; note E v2 section 0). Assumptions: independent groups, equal variance in both groups for a gene under H0, finite fourth moments. **No distribution.**
* **Why relative weights:** on the variance scale the sums are dominated by a few very noisy genes. The paper's effective dimension (below) is only about 4-9 of 3000 genes in the toy with realistic heterogeneity. Relative weights (each gene divided by its trend value) give the same lambda and the same performance in the toy but a 2-4x larger effective dimension (section 8). For the median target both are identical.
* **Distribution-free justification (to state in the methods text):** with `T` fixed and `s2_g` unbiased, `E sum_g (lambda T_g + (1-lambda) s2_g - sigma2_g)^2` is minimised at `lambda* = sum_g Var(s2_g) / sum_g E[(s2_g - T_g)^2]` (two moments only); `lambda` above is its plug-in estimate (the Ledoit-Wolf / Schaefer-Strimmer construction). The weighted version minimises the relative squared error.
* **Relation to the classical James-Stein estimator** (Gao and Hitchcock 2010, eq. (1), and the background they review): the positive-part James-Stein estimator towards a target `nu` is `nu + (1 - (p_hat - 2) / ||X - nu||^2_Q)^+ (X - nu)`. With `X` = the vector of gene-wise variances, `nu` = the target and `Q = diag(Var(s2_g))`, its shrinkage weight `(p_hat - 2) / ||X - nu||^2_Q` becomes, for thousands of genes, the moment ratio above; clipping to [0, 1] is the positive part (Baranchik 1964, as cited there). Note: the dominance theorems behind James-Stein (Stein; Bock 1975) are proven for **normal means**; we do not claim dominance, only the moment-optimal lambda.
* **Effective dimension diagnostic** (Bock 1975, as used by Gao and Hitchcock 2010): `p_hat = tr(Q) / lambda_max(Q)`. For diagonal `Q`: `p_hat = sum_g w_g Vhat_g / max_g w_g Vhat_g`; report also the participation ratio `(sum w Vhat)^2 / sum (w Vhat)^2`. These say how many genes effectively determine lambda. Bock's condition `p_hat > 2` must hold; a small value flags that lambda rests on a few genes.
* Report per run: `lambda_pool`, `lambda_c`, `lambda_t` (also unclipped), `p_hat`, participation ratio, `n_floored`.

---

## 3. Cross-fitted scaling (`scaling = "xfit_js"` or `"xfit_raw"`)

```
observed   t_g     = (ybar_c,g - ybar_t,g) / sqrt(s2s_g)                     (pooled shrunk variance)
case pool  z_c,j,i = eps_c,j,i / sqrt(s2s_t,j)    for neighbour j            (CONTROL-side variance)
ctrl pool  z_t,j,i = eps_t,j,i / sqrt(s2s_c,j)                               (CASE-side variance)
null pair  v       = z_c / sqrt(n_c) - z_t / sqrt(n_t);   p_g = (1 + #{|v| >= |t_g|}) / (N + 1)   (compute_pvalue unchanged)
xfit_raw:  lambda = lambda_c = lambda_t = 0 (raw variances, floor VAR_FLOOR);   xfit_js: lambdas from section 2
```

* Cross-fitting makes the scale of a neighbour's case residuals independent of those residuals (computed from the control samples), so the null is not internally studentised. Same-group scaling (note E v1) gave 25x too many null p <= 0.001 at n = 3 in the toy.
* Benchmark-only negative control: `shrink_scaling = "same_group"`.
* **Parity:** `xfit_js` with lambda forced to 1 and a constant target gives the current exact model's p-values exactly.
* The adaptive `tau` growth of the exact model then acts on scaled residuals (**VERIFY** behaviour).
* Equal variance between the groups is assumed for the scaling; under unequal group variances the cross-fitted null is too wide (conservative, power loss; E v2 section 8). The trend-level group-ratio correction stays **parked**.

---

## 4. Band with a James-Stein matching variable (`match_variable = "js"`, replaces note B)

```
u_g     = log( s2s_g / T_g )                         (James-Stein-shrunk, detrended; median or LOESS target as in section 2)
tau_var = c_js * 1.4826 * MAD_g(u_g)                 (empirical width; no trigamma)
```
* Acceptance as in note A: `|u_j - u_g| <= tau_var`, two-sided, fixed `k`, scan cap and short-pool rules unchanged.
* `match_variable = "raw"` keeps note A exactly (including `tau_var = c * sqrt(psi1(d/2))`); that is the "band as it is now" arm.
* No toy evidence for this arm yet; `c_js` grid in section 6.

---

## 5. Adaptive variance-aware growth (`variance_mode = "adaptive"`), PROPOSAL

Asis's idea: grow the neighbourhood in blocks, accept only neighbours of comparable variance, stop early when the neighbourhood stops being comparable. Per target gene and per side (case pool around the case mean, control pool around the control mean, as now):

1. Visit candidates in mean-closeness order (existing `choose_index`; the target itself first, self-inclusion as now).
2. Take the next **block of b candidates** (default b = 5). Candidates whose mean distance `|log2(mean_j + 1) - log2(mean_g + 1)|` exceeds the **growth radius `r_max`** are not visited: radius stop.
3. Accept the candidates of the block that pass the variance band (`match_variable` raw or js, width as in section 4 / note A).
4. Once at least `k_min` genes are accepted: compute the pool scale statistic (the quantity the current `tau` rule uses, mean |residual| of the pool) with and without the newly accepted genes; if the **absolute** relative change `|S_new / S_old - 1|` exceeds `tau_stab`, **discard the block and stop** (stability stop). Two-sided on purpose: the pool may become neither noticeably noisier nor noticeably quieter (the current exact-model `tau` rule only stops on increases).
5. Also stop when: `k_max` genes are accepted; the acceptance rate over the last `w` blocks is below `a_min` (the walk has left the comparable region); the scan cap (`20 * k_max` slots) is reached; candidates are exhausted.
6. Short pools: as note A (keep if at least `BAND_MIN_GENES` genes, otherwise p = -1).
7. **Resolution:** a smaller pool raises the gene's p-value floor `1/((k * n)^2 + 1)` (n = 3: 2.8e-4 at 20 genes vs 4.4e-5 at 50). Report each gene's floor and pool size; consider an n-dependent `k_min` (e.g. 50 at n = 3) if early stops cluster at small pools.

Defaults to start from (owner/Asis confirm): `b = 5, k_min = 20, k_max = 50, r_max = 0.5, tau_stab = 0.10, w = 2, a_min = 0.2`.
Per-gene diagnostics: genes accepted, slots scanned, final radius, p-value floor, **stop reason** (1 k_max, 2 radius, 3 stability, 4 acceptance rate, 5 scan cap, 6 exhausted).
**Parity test:** `b = 1, tau_var = Inf, r_max = Inf, tau_stab = Inf, a_min = 0, k_min = k_max = k` gives exactly `variance_mode = "band"` with an infinite band, which gives the exact model with `k_start = k_step = k_max = k`.
**Open alternative** (if "growth radius" was meant in two dimensions): accept `j` when `((m_j - m_g)/r_max)^2 + ((u_j - u_g)/tau_var)^2 <= 1`. Implement the one-dimensional version first.

---

## 6. Fortran / R changes and tests

**Code (Tensor-Omics, `noise_model_variance`):**
* Add `compute_moment_shrinkage_helper` (section 2): pure, no allocation; the LOESS is called only for `shrink_target = "loess"`; the median needs a sort (**VERIFY** whether the median helper used for `delta_hat` in `noise_model_exact` can be reused).
* Cross-fitted scaling at pool packing (section 3); the observed statistic passed as `|beta_g| / sqrt(s2s_g)`.
* Adaptive gather routine (section 5), next to the band gather (shares `choose_index`, the acceptance test, the scan cap).
* Remove the note-B `d0` path; `match_variable = "js"` replaces it.
* Wrapper: `tox_compute_noise_pvalues_pipeline_variance(..., variance_mode = c("mean", "band", "adaptive"), match_variable = c("raw", "js"), scaling = c("none", "xfit_raw", "xfit_js"), shrink_target = c("median", "loess"), lambda_scale = c("relative", "absolute"), c_tau, c_js, k, b, k_min, r_max, tau_stab, a_min, loess_robust)`. Scalar outputs: the lambdas, `p_hat`, participation ratio, `n_floored`; per-gene outputs: `s2s_g`, adaptive diagnostics.
* Documentation per the project convention (current_state, eval doc, labbook); methods text once the tests decide.

**Unit tests (`test/mod_test_noise_model_variance.f90`):**
1. Regression: `variance_mode = "band"`, `match_variable = "raw"`, `scaling = "none"` reproduces the current band results bit-identically (fixture stored before editing).
2. Parity: `variance_mode = "mean"`, `scaling = "none"` equals the exact model; `xfit_js` with lambda = 1 and constant target equals the exact model; adaptive parity (section 5).
3. Contrast `Vhat` unbiased (normal, skewed, t5; n = 3, 5, 10); documentation test that the v1 formula equals `s2^2/4` at n = 3.
4. lambda close to 1 for homogeneous variances; clearly below 1 for heterogeneous ones at n = 10; relative = absolute for a constant target.
5. Calibration regression (fixed seeds, n = 3, heterogeneous variances): `xfit_js` null `p <= 0.001` at most 2x nominal; `same_group` clearly above.
6. `xfit_raw` floor path (a group variance of 0); adaptive stop reasons each triggered by a constructed case; validation errors.

---

## 7. Benchmark arms and parts (first round)

```
references            limma, edgeR, DESeq2, TOX-log|TPM-centred
TOX-log-vband-c3      variance_mode "band", match "raw", scaling "none"                 (a) band as it is now
TOX-log-vband-js      variance_mode "band", match "js", c_js in {0.5, 1}, scaling "none"     replaces vshrunk
TOX-log-xfit-raw      variance_mode "mean", scaling "xfit_raw"                            (b) cross-fit, no shrinkage
TOX-log-xfit-js       variance_mode "mean", scaling "xfit_js", target "median"            (c) cross-fit + James-Stein   PRIMARY
TOX-log-xfit-js-loess same, target "loess", lambda "relative"
TOX-log-xfit-js-nocf  same as xfit-js, shrink_scaling "same_group"                        negative control (RS n = 3, 5 and P+ only)
TOX-log-adapt         variance_mode "adaptive", match "raw", scaling "none", r_max in {0.25, 0.5, 1.0}
TOX-log-adapt-xfit-js variance_mode "adaptive", match "js", scaling "xfit_js"             optional, if time permits
```

* All TOX arms on the input of the comparator (`tpm_centred`).
* **Parts:** RS (n = 3, 5, 7, 10; pi1 in {0, 0.1}; 100 rounds; six TCGA base cohorts) is primary; R (full-cohort halves); P+ (sigma_d in {0, 0.53, 1.05}, n = 3, 5, 10, pi1 in {0, 0.1}; plus the rho in {1, 1.25, 1.5, 2} extension for the `xfit` arms); V (hypervariable genes).
* **Choose the free parameters (`c_js`, `r_max`) on the pi1 = 0 rounds only**, then judge power on pi1 = 0.1 (note D protocol).
* Report, per cohort and n: null FPR at 0.001 / 0.01 / 0.05, observed FDR at nominal 0.05 and 0.1, TPR at nominal and achieved FDR, AUPRC excess, |LFC|80, gap to limma/edgeR/DESeq2, all **by expression decile** too; the lambdas, `p_hat`, participation ratio; adaptive stop-reason shares and radii.
* Decision criteria: note E v2 section 6 (calibration and tail gates first, then FDR and power); proposals until the owner confirms.

---

## 8. Expected behaviour from the toys (sanity checks, NOT targets)

Python toys, normal noise, mean-variance trend, gene scatter sigma_v = 1.05, 10 % DE, 3000 genes, 4-5 datasets per cell (differences below about 0.03 are noise). FDR / TPR at BH 0.05 (AUPRC excess):

| n | TOX exact (current) | cross-fit, no shrinkage | cross-fit + James-Stein (median target) |
|---|---|---|---|
| 3 | 0.115 / 0.274 (0.551) | 0.010 / 0.044 (0.574) | 0.047 / 0.368 (0.628) |
| 5 | 0.091 / 0.400 (0.636) | 0.011 / 0.378 (0.684) | 0.044 / 0.535 (0.696) |
| 10 | 0.078 / 0.623 (0.768) | 0.028 / 0.739 (0.802) | 0.050 / 0.768 (0.802) |

* **Cross-fitting without shrinkage is valid but very conservative at small n** (raw per-group variances with 2-4 degrees of freedom make a heavy-tailed null). Shrinkage is what makes cross-fitting usable at n <= 5; at n = 10 the two are close.
* LOESS target with relative lambda: same lambda (0.41 / 0.30 / 0.14) and performance as absolute, effective dimension 34 / 36 / 57 instead of 9 / 14 / 26.
* No toy evidence for `vband-js` and `adapt`: their behaviour is open.
* Limits: normal noise, equal group variances, mean-bin pools instead of per-gene kNN, no zeros. Real data decide.

---

## 9. What stays distinct from limma, edgeR, DESeq2 (Asis's caution)

As far as I know (please check): (1) the null is the empirical distribution of neighbouring genes' residual differences with exact pairwise counting, not a t or NB reference; (2) the shrinkage is distribution-free (moments only) and can be turned off, whereas the reference tools use parametric priors (scaled inverse chi-square / log-normal dispersion prior); (3) cross-fitted scaling of an empirical null; (4) neighbourhoods selected by mean and (optionally) variance. The tests should show whether the combination is competitive; being distinct is not by itself an argument.

**References:** Gao J, Hitchcock DB (2010), James-Stein shrinkage to improve k-means cluster analysis, *Computational Statistics and Data Analysis* 54:2113-2127, doi:10.1016/j.csda.2010.03.018 (owner-provided; formulas checked against the PDF). Cited therein: James and Stein (1961), Baranchik (1964), Bock (1975). Opgen-Rhein and Strimmer (2007) and Schaefer and Strimmer (2005): from memory, **verify**.
