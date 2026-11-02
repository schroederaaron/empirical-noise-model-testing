# Implementation note E: distribution-free (moment-based) variance shrinkage towards the LOESS trend

**For:** Claude Code. **Status:** DRAFT for discussion (owner and Asis). Written without executing code; everything marked **VERIFY** must be checked before relying on it.
**Builds on:** notes A (variance band), B (shrunk band, implemented as `variance_mode = 2`) and D (Part RS). **Repos:** Tensor-Omics (branch `145-empirical-noise-model`, module `noise_model_variance`) and `empirical-noise-model-testing`.

---

## 1. Why this note exists

**Findings (Part RS, run of 27/10/2026; six TCGA base cohorts, 100 rounds per cell):**

* Null arm (pi1 = 0), share of null p <= 0.05: `TOX-log|TPM-centred` 0.054-0.067 at n = 3 and 0.047-0.063 at n = 5, i.e. close to nominal. limma 0.040-0.050.
* Power arm (pi1 = 0.1), observed FDR at nominal 0.05: `TOX-log|TPM-centred` 0.22-0.53 at n = 3 and 0.14-0.31 at n = 5. Variance band c3: 0.16-0.37 and 0.07-0.13. DESeq2: 0.20-0.43 at n = 3. limma: 0.005-0.071.
* The current shrinkage variant (note B) does not change this at n = 3. On the null arm, B at c = 1 is 0.064-0.090 and A at c = 1.5 is 0.064-0.085; FDR is 0.28-0.47 for both. A wider band helped more than shrinkage did.
* At n = 3 the variance-band arms rank genes nearly as well as limma, edgeR and DESeq2 (AUPRC excess in `sample_size_curve.png`). The problem is the tail calibration of the p-values, not the ranking.

**Diagnosis (hypothesis, to be tested):**

1. **The mean-matched model never uses the gene's own variance.** The null is the pooled neighbour residuals, so the noise scale comes entirely from the local trend: this is full pooling, the extreme of shrinkage. A gene whose true noise is above that of its neighbours is tested against too narrow a null. At large n the variance band fixes this by choosing neighbours with similar noise. At n = 3 the selection variable is too noisy (sampling SD of the log variance about 0.8), so the selection is close to random.
2. **Note B shrinks only the selection, never the test.** limma, edgeR and DESeq2 shrink the variance (or dispersion) that enters the test statistic itself. That is what makes them work at small n.
3. **Note B is not distribution-free.** Its prior degrees of freedom `d0` come from moment matching with the trigamma function, which assumes that log variances scatter like a log-chi-square (normal errors). Note A's band width `c * sqrt(psi1(d/2))` rests on the same assumption.

**Owner constraints (20-27/10/2026):**

* No distributional assumption: no parametric prior, no chi-square sampling model.
* Shrinkage towards the data-driven LOESS trend is wanted.
* Shrinkage is a user switch (on/off). The manual should state from which n it is no longer needed; that n comes from the tests, not from a guess (expected around 7).
* Goal: competitive with edgeR, limma and DESeq2 under most conditions, not to beat them. The 1-D model is the working tool; the multi-dimensional model is the niche.

---

## 2. The estimator: moment-based shrinkage of gene-wise variances

James-Stein-type shrinkage with a data-driven intensity, as in the "shrinkage t" of Opgen-Rhein and Strimmer (2007). It uses only moments and assumes no distribution for the data or the variances. The target is our LOESS trend instead of their median.

**References (VERIFY before citing):** Opgen-Rhein R, Strimmer K (2007), *Statistical Applications in Genetics and Molecular Biology* 6:9 (shrinkage t). Schaefer J, Strimmer K (2005), *SAGMB* 4:32 (unbiased estimator of the variance of a sample variance, used below). Reference implementation for a cross-check: `corpcor::var.shrink` (target = median; **VERIFY** its exact formula and lambda clipping).

### 2.1 Per-gene quantities

All on the scale the model tests on: `y = log2(x + 1)` by default (`norm_method = 1`); the linear branch uses `y = x`.

For gene g and group k in {case, control}, with `n_k` replicates:

```
w_kgi   = (y_kgi - ybar_kg)^2                         i = 1..n_k
wbar_kg = mean_i w_kgi
s2_kg   = n_k / (n_k - 1) * wbar_kg                   unbiased group variance
V_kg    = n_k / (n_k - 1)^3 * sum_i (w_kgi - wbar_kg)^2   unbiased estimate of Var(s2_kg)   (Schaefer-Strimmer; VERIFY)
```

Pooled over the two groups (independent):

```
d       = n_case + n_ctrl - 2
s2_g    = ((n_case - 1) * s2_case,g + (n_ctrl - 1) * s2_ctrl,g) / d
V_g     = ((n_case - 1)^2 * V_case,g + (n_ctrl - 1)^2 * V_ctrl,g) / d^2     estimate of Var(s2_g)
```

`s2_g` is exactly the pooled variance `compute_logvar_helper` already computes (`s2`); only `V_g` is new.

### 2.2 Target: the LOESS trend on the variance scale

```
T(m) = LOESS of s2_g on m_g,   m_g = log2(mean of all samples of g + 1)   (the same trend axis as notes A/B)
```

* Fit on the **variance scale** (`s2`), not on `log s2`. `exp(LOESS of log s2)` underestimates the mean variance (Jensen), and correcting that needs a distribution (note B's digamma term). On the variance scale no correction is needed.
* **Robust or not (open, to test):** the bisquare-robust LOESS downweights hypervariable genes, which makes the target "typical" rather than "mean" variance. Default proposal: non-robust (degree 1, span 0.3), with `loess_robust` as the existing switch. Report both in the first runs.
* Floor: `T(m) >= VAR_FLOOR` (existing constant).

### 2.3 Shrinkage intensity and shrunk variance

```
lambda_hat = clip( sum_g V_g / sum_g (s2_g - T(m_g))^2 , 0, 1 )
s2s_g      = lambda_hat * T(m_g) + (1 - lambda_hat) * s2_g
```

* **What lambda means:** the numerator is how noisy the gene-wise variance estimates are (sampling error); the denominator is how far they scatter around the trend in total. If all scatter is sampling noise, lambda -> 1 and every gene is pulled onto the trend. If genes truly differ in noise, lambda stays small.
* **One lambda** for all genes (global), as in Opgen-Rhein and Strimmer. Each V_g is very noisy at n = 3, but lambda pools thousands of genes. **Optional, phase 2:** lambda per expression bin (e.g. deciles of m) if the noise heterogeneity clearly depends on the mean; report the per-decile values as a diagnostic from the start.
* The formula treats T as fixed. It is estimated from thousands of genes, so its own variance is negligible next to `V_g` (note this in the methods text).
* Report `lambda_hat`, the median `V_g`, and the share of genes with `s2s_g` within 10 % of `T(m_g)`.

---

## 3. Where the shrunk variance is used: two levels

### Level 1: selection only (distribution-free successor of note B)

Replace note B's matching variable and band width:

```
u_g     = log( s2s_g / T(m_g) )
tau_var = c * 1.4826 * MAD_g(u_g)          empirical spread of u, no trigamma
```

* Everything else as in notes A/B: two-sided band `|u_j - u_g| <= tau_var`, fixed k, scan cap, minimum genes.
* Degenerate case: if lambda_hat = 1 then u = 0 for all genes; pass `tau_var = 1` so the band code path accepts every candidate (as note B does for `d0 = Inf`). The model then equals the plain fixed-k model.
* New mode: `variance_mode = 3` ("mshrunk"), so A (1) and B (2) stay available for comparison.
* Expected effect: modest. Like B it only changes which neighbours are chosen. It removes B's parametric `d0` and makes the comparison A vs B vs mshrunk clean.

### Level 2: studentised exact null (the switch the owner asked for)

The shrunk variance enters the **test**, as in the reference tools, while the null keeps its empirical shape. Neighbours: the plain mean-matched neighbourhood (adaptive kNN of the exact model). The band is not needed, because each gene's own scale now accounts for noise differences; whether a band on top still helps is a test arm, not a default.

```
observed   t_g     = (ybar_case,g - ybar_ctrl,g) / sqrt(s2s_g)                     (the scale factor sqrt(1/n+ + 1/n-) cancels, see below)
residuals  z_kji   = eps_kji / sqrt(s2s_j)          for every neighbour gene j in the pool of side k (eps with the Bessel factor, as now)
null pair  v       = z_case / sqrt(n_case) - z_ctrl / sqrt(n_ctrl)
p_g        = (1 + #{ |v| >= |t_g| }) / (N + 1)       exact pairwise count, unchanged code (`compute_pvalue`)
```

* **Implementation is small:** scale each neighbour's residuals by `1/sqrt(s2s_j)` when they are added to the pool (the pools already carry per-side `1/sqrt(n)`), and pass `|beta_g| / sqrt(s2s_g)` as the observed statistic. `compute_pvalue` is reused as is. The moderated t would divide observed and null alike by `sqrt(1/n+ + 1/n-)`; dropping that constant on both sides leaves the comparison `|v| >= |t|` unchanged.
* **Switch:** `shrink = FALSE` (default until the tests decide) gives exactly the current exact model. `shrink = TRUE` gives the studentised null with moment shrinkage.
* **Limits, as built-in sanity checks:**
  * `lambda_hat = 1` and a flat trend: `s2s` is the same constant for every gene, the scaling cancels, and p-values equal the current exact model **exactly**.
  * `lambda_hat = 1` and a non-flat trend: close to the current model; neighbours differ only by `T(m_j) / T(m_g)`.
  * `lambda_hat = 0`: each gene is scaled by its own sample variance: a per-gene t with an empirical null of standardised residuals. Unstable at n = 3, which is exactly why lambda is estimated.
  * So lambda interpolates between full pooling (today's TOX) and gene-wise studentisation, and the data choose where.
* **Distributional assumptions:** none for the data, the variances or the null. The shrinkage needs finite fourth moments of the residuals; the null shape stays empirical. This is the property to state in the manual.
* **Guard:** `s2s_g >= VAR_FLOOR` before dividing. Report the number of genes floored.

### Multi-dimensional model (phase 2, after 1-D is decided)

Level 2 carries over per axis: compute `lambda_a`, `T_a`, `s2s_{g,a}` per axis from that axis' replicates, studentise `beta_{g,a}` and the axis pools, and keep the product null. `Var_null(a)` then follows from the standardised pools as now. Note that the MD model currently uses the mean-matched neighbourhood only (the variance band is 1-D), so level 2 is the natural first variance-aware step for MD as well. Out of scope for the first implementation.

---

## 4. Fortran / R changes (1-D)

* `noise_model_variance` (shares `compute_logvar_helper`, the LOESS and the gather): add `compute_moment_shrinkage_helper` (inputs: replicate matrices, `m`, trend settings; outputs: `s2`, `V`, `T`, `lambda_hat`, `s2s`). Pure, no allocation, as the other helpers. **VERIFY** that `tox_loess::loess_alloc` can fit `s2` (variance scale) with the same settings used for `l`.
* Level 1: `variance_mode = 3` in `compute_variance_match`.
* Level 2: a `shrink` argument in the variance pipeline (mode 0 = no band, i.e. the plain mean-matched gather, which the module must then offer; **VERIFY** whether reusing `noise_model_exact::gather_residuals_helper` is cleaner than a band with `tau_var = huge`), plus the residual scaling in the pool and the scaled observed statistic. Outputs per gene: `s2s`; scalar outputs: `lambda_hat`, `n_floored`.
* C wrapper, Rcpp and R wrapper (`tox_compute_noise_pvalues_pipeline_variance(..., variance_mode = c("band", "shrunk", "mshrunk", "none"), shrink = FALSE)`); `norm_method = 1` stays the default.
* Documentation per CLAUDE.md: current_state, eval doc (new section), labbook. Methods manuscript Section 11.2.5 ("in development") is rewritten once the tests decide.

---

## 5. Tests

### Unit tests (Fortran, `test/mod_test_noise_model_variance.f90`)

1. `V_g` against a brute-force computation on a small fixture; and its unbiasedness by simulation (mean of V over many draws equals the empirical variance of s2, within MC error; any residual distribution, e.g. normal and a skewed one).
2. lambda_hat: homogeneous true variances (all genes on the trend) -> lambda close to 1; strongly heterogeneous true variances at n = 10 -> lambda clearly below 1.
3. Cross-check against `corpcor::var.shrink` with the trend replaced by a constant (median) target, on a saved fixture (if available in the R environment; tolerance 1e-8). **VERIFY** the target conventions first.
4. Level 2 with a constant target and lambda forced to 1: p-values identical to `noise_model_exact` (bit-exact or 1e-15).
5. Level 2 invariance: multiplying all data of one gene's neighbourhood by a constant does not change that gene's p-value (scale-free statistic). On the log scale, test the linear branch.
6. Level 1 with lambda = 1: identical to the plain fixed-k gather (as note B's degradation test).
7. Floor and validation paths (constant genes, `n_k = 2`).

### Benchmark (testing repo, `power_test.R`)

New arms (labels proposed):

```
TOX-log-mshrunk        engine "variance", variance_mode = "mshrunk", c in {1.5, 2, 3}     (level 1)
TOX-log-stud           engine "variance", variance_mode = "none", shrink = TRUE           (level 2, mean-matched)
TOX-log-stud-vband-c3  engine "variance", variance_mode = "band", c = 3, shrink = TRUE    (level 2 + band, optional)
```

* Parts: **RS** (primary), R, P+ (incl. pi1 = 0), V (sparse hypervariable genes: does shrinking towards the trend hide true outliers?).
* **Extend RS to n = 7 and n = 10.** All six base cohorts have S >= 26 samples (`power_rs_checks.csv`), so both fit everywhere; note that at n = 10 a round uses up to 77 % of the smallest cohort (LUAD Stage IV, S = 26), so rounds overlap strongly there. This is what decides "from which n is shrinkage not needed": compare `shrink = TRUE` vs `FALSE` per n.
* **Side arm for the owner's p-value-floor / neighbourhood hypothesis:** the existing `TOX-log|TPM-centred` with a larger neighbourhood (e.g. k_max 100 and 200) and a stricter tau, at n = 3 and 5. Report FPR at 0.001 and 0.01 on the null arm, since BH works in the tail.
* Report additionally: `lambda_hat` per cohort and n; FPR at 0.001 / 0.01 / 0.05 on the null arm; observed FDR at nominal 0.05 and 0.1; TPR at nominal and achieved FDR; AUPRC excess; gap to limma, edgeR, DESeq2.

---

## 6. Decision criteria (PROPOSAL: owner confirms before the run)

Reuse note D section 6 (calibration gate first, then power), with these additions:

1. **FDR gate (new, because RS showed calibrated null p-values but uncontrolled FDR):** observed FDR at nominal 0.05 at most 0.075 in at least 5 of 6 cohorts. Run per n.
2. **Shrinkage helps at n** if, with `shrink = TRUE` vs `FALSE` (same arm otherwise): the FDR gate passes with shrinkage and not without, or AUPRC excess rises by at least 0.02 in at least 5 of 6 cohorts with no loss in calibration.
3. **"Not needed from n on":** the smallest n in {3, 5, 7, 10} from which `shrink = FALSE` passes the calibration and FDR gates and `shrink = TRUE` gives no AUPRC gain of 0.02 or more. That n goes into the manual.
4. **Default of the switch:** decided after 1-3. If shrinkage helps at small n and costs nothing at large n, `shrink = TRUE` becomes the default.
5. Report every arm, including failures.

---

## 7. Risks and open points

* **Hypervariable genes:** shrinking towards the trend pulls genuinely noisy genes in. Part V answers whether that costs FDR there; global lambda may be too strong or too weak for them (per-bin lambda as the fallback).
* **The FDR problem may lie elsewhere.** If level 2 does not fix FDR at n = 3, candidates are the tail shape of the null (pool size, discreteness of log2(x + 1) at low counts) and the low-count genes (the MD example). The side arm in section 5 tests the pool-size part.
* **Equal variance across groups** in `s2_g` (pooled), as in limma; the null keeps per-side pools. Per-side shrinkage is a possible extension.
* **Remaining parametric pieces elsewhere:** note A's `tau_var = c * sqrt(psi1(d/2))` assumes log-chi-square sampling; level 1 replaces it by the empirical MAD of u. Decide whether mode A should switch too.
* **Out of scope here:** MD implementation (section 3, phase 2), bootstrap and blocked engines, per-bin lambda before the global version has been tested.
