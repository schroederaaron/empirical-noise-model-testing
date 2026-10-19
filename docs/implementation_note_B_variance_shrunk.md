# Implementation note B — variance band on a shrunk matching variable ("vshrunk")

**For:** Claude Code. **Prerequisite:** note A is implemented (Fortran `match_mode = 1`, `match_logvar`, `tau_var`, `n_scanned_*`, `tox_match_logvar()` for model "A",
arm plumbing, report columns). This note adds ONLY: the shrunk matching variable (R side), its tau rule, one arm family, tests, and checks. **No further Fortran changes.**
**Status:** written without executing any R. Reference prototype: `toy_variance_moderation.py` (`fit_eb`, `moderate`) and `toy_variance_neighbourhood.py` (row "band two-sided, shrunk (EB)").

---

## 1. Goal and hypothesis

Model A matches on the **raw** detrended log-variance. With few replicates that variable is noisy (sampling SD of log s2 is 0.80 at n=3, 0.53 at n=5, 0.34 at n=10), and selecting neighbours on the
sample variance of the same residuals that form the null makes the null mimic the target's own estimate, which cost calibration in the toy when genes do **not** differ in noise.
Model B smooths the matching variable toward the data-driven LOESS trend before the band is applied.

* The **null stays purely empirical** (unchanged pairwise residual differences, exact tail count). Only the neighbour-selection variable is smoothed.
* The trend is the data-driven LOESS fit of note A, not an assumed one.
* **Graceful degradation:** if genes do not differ in noise beyond sampling error, the estimated prior df d0 becomes infinite, the matching variable is identically 0, every candidate is
  eligible, and Model B reduces exactly to current TOX with fixed k.
* **Hypothesis to test:** B keeps most of A's ranking gain at heterogeneity and removes A's calibration cost when there is none.

Note (for interpretation, do not "fix"): the one parametric element is the **estimate of d0** (moment matching with the trigamma sampling-variance correction, i.e. a log-F type prior). It sets only the smoothing
weight for neighbour SELECTION, never the null.

---

## 2. Definitions (the formula that was actually prototyped)

Use the quantities of note A: `L`, `d = n_c + n_t - 2`, `s2_g`, `l_g = log(s2_g + var_floor)`, `m_g`, the robust LOESS settings (span 0.3, degree 1, symmetric).
`psi` = digamma, `psi1` = trigamma, `uniroot` for the inverse.

1. Bias-corrected log variance: `e_g = l_g - psi(d/2) + log(d/2)`.
2. Trend: `T(m) =` robust LOESS of `e` on `m`. Detrended scatter: `r_g = e_g - T(m_g)`.
3. Excess scatter above sampling noise: `ev = V - psi1(d/2)` with `V = (1.4826 * mad(r))^2` (**default, robust**, matches the choice in `estimate_dispersion_trend.R`) or `V = var(r)` (option `robust_var = FALSE`).
4. Prior df: if `ev <= 0` then `d0 = Inf`; else `d0 = 2 * y` where `trigamma(y) = ev` (root of `psi1(y) - ev` on `[1e-6, 1e6]`; `psi1` is decreasing).
5. Prior variance at each gene's mean: `s0sq(m) = exp( T(m) + psi(d0/2) - log(d0/2) )`; if `d0 = Inf` then `s0sq(m) = exp(T(m))`.
6. **Moderated (shrunk) variance, on the VARIANCE scale (limma-style):** `s2m_g = ( d0 * s0sq(m_g) + d * s2_g ) / ( d0 + d )`; if `d0 = Inf` then `s2m_g = s0sq(m_g)`.
7. **Matching variable:** `u_g = log( s2m_g / s0sq(m_g) )`  (which is `log((d0 + d*q_g)/(d0 + d))` with `q_g = s2_g / s0sq(m_g)`; identically 0 when `d0 = Inf`).
8. Weight on the gene's own variance (diagnostic, reported): `w = d / (d + d0)` (0 when `d0 = Inf`).
9. **Band (two-sided)** exactly as in note A: eligible iff `|u_j - u_g| <= tau_var`, fixed `k = 50`.
10. **tau rule:** `tau_var = c * w * sqrt(trigamma(d/2))`, primary `c = 1`; sensitivity `c = 0.5` and `c = 2`.
    If `d0 = Inf` (so `w = 0`, `u = 0`): pass `u = 0` for all genes and `tau_var = 1` so the mode-1 code path is still exercised and accepts everything.

**Correction to earlier chat wording:** the shrinkage is applied on the variance scale (step 6) and the log is taken afterwards. It is NOT a log-scale convex combination. The two differ for very quiet genes:
on the variance scale `u` is bounded below by `log(d0/(d0+d))`, i.e. very quiet genes are treated as equivalent (the prior dominates their noisy estimate), while the upper side behaves like Model A shifted by `log(w)`.
Implement step 6 exactly as written.

---

## 3. R changes (benchmark repo)

Extend `common/tox_variance_match.R`: implement `model = "B"` in `tox_match_logvar()` returning, in addition to note A's fields: `d0`, `w`, `ev`, `scatter`, `u` (the shrunk variable), `s2m`, `s0sq`.
`tox_resolve_tau("B", c, d, w)` as in step 10.

### 3.1 Unit tests (testthat)

* **d0 recovery:** simulate `G = 15000` genes from the limma model: `sigma2_g = s0^2 * d0_true / rchisq(1, d0_true)`, then `s2_g = sigma2_g * rchisq(1, d)/d`, for `d0_true in {3, 8, 30}` and `d in {4, 8, 18}`.
  The classical-variance option must recover `d0_true` within 10%; report the robust option's bias (it is expected to differ, record by how much).
* **No heterogeneity:** constant `sigma2` -> `d0` is `Inf` or so large that `w < 0.05` in at least 95% of repeated simulations; `u` about 0.
* **Cross-language check (if Python is available):** with a constant trend (set `span` large / use `m` constant), `d0`, `s0sq` and `s2m` equal the Python `fit_eb` and `moderate` outputs for a saved fixture matrix (classical variance mode) to 1e-6.
* `u` is monotone in `s2_g` at fixed `m`; bounded below by `log(d0/(d0+d))` for finite `d0`.
* `d0 = Inf` returns `u == 0`, `w == 0`, and `tox_resolve_tau` returns `tau_var = 1`.

---

## 4. Benchmark arms

```r
list(label = "TOX-log-vshrunk",      engine = "exact", norm = 1L, null = 0L, input = "tpm", match = "B", c_tau = 1.0, k = 50L)
list(label = "TOX-log-vshrunk-c0.5", ..., c_tau = 0.5)     # sensitivity
list(label = "TOX-log-vshrunk-c2",   ..., c_tau = 2.0)     # sensitivity
```

Same parts, same `input` rule as note A (P+, D, R; in R use the comparator's input). Manifest additions: `d0`, `w`, `ev`, `scatter_robust`. The report adds `d0` and `w` columns to the summary so they can be read next to the results.

---

## 5. Checks specific to B (in addition to note A section 7)

1. **Exact degradation test (the key property):** on P+ with `sigma_d = 0` and on the pi1 = 0 null, `TOX-log-vshrunk` must reproduce `TOX-log` (fixed k = 50) within Monte-Carlo error; confirm `w` is about 0 there. If it does not, stop and investigate before reading any power result.
2. **Look at the estimated d0 on the real cohorts** (the six Part R cohorts, healthy and Stage IV) and on the P+ grid. Expected from the earlier dispersion estimate (scatter about 1.035): d0 about 2.7, i.e. `w` about 0.60 / 0.75 / 0.87 at n = 3 / 5 / 10. Report the realised values; large departures
   (very different d0 per cohort, or d0 near 0) must be discussed before interpreting the power numbers.
3. **Robust vs classical variance in the d0 step:** run Part R once with `robust_var = FALSE` as a sensitivity arm label (`TOX-log-vshrunk-classicalvar`) to see whether hypervariable genes drive d0.
4. **LOESS span sensitivity:** repeat the calibration step for `span in {0.15, 0.3, 0.5}`; the results should not change materially.
5. **Compare A vs B directly** (paired, same rounds): FDR at nominal on the null halves and P+ `sigma_d = 0`, AUPRC excess on P+ and Part R, and the pool diagnostics (`frac_short_pool`, `n_scanned_median`).

Select `c_tau` on calibration only (note A, section 7, step 1).

---

## 6. Expected behaviour from the Python toy (sanity checks, NOT targets)

Row "band two-sided, shrunk (EB), tau = 0.35", single mean bin, normal noise, 10% DE (6 datasets per cell; differences below about 0.03 are noise):

* sigma_v = 0: **identical to current TOX** (AUPRC 0.403 / 0.586 / 0.733 at n = 3 / 5 / 10; FDR 0.021 / 0.042 / 0.029) because `d0 = Inf`.
* sigma_v = 1.05: AUPRC 0.368 / 0.572 / 0.740 (current TOX: 0.287 / 0.413 / 0.662; moderated t: 0.442 / 0.618 / 0.748); FDR 0.47 / 0.25 / 0.11 (still above nominal at small n); TPR 0.50 / 0.57 / 0.69.
* So B is expected to **fix the no-heterogeneity cost of A** and keep most of A's ranking gain, but **calibration at n <= 5 under heterogeneity remains too liberal**, and n = 3 stays weak.

---

## 7. Risks and out of scope

* The weight rests on the d0 estimate; if the log-variance scatter is far from log-F-like (heavy tails, hypervariable genes), `w` can be mis-set. That is why robust and classical variants are both run, and why `c_tau` is chosen on calibration.
* Shrinkage of the matching variable is global (one `d0`); no gene-specific prior.
* Out of scope: any shrinkage of the test statistic or of the null, per-side variance, cross-fitting of the selection variable, bootstrap/blocked engines, a studentised null.
