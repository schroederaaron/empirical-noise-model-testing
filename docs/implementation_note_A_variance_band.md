# Implementation note A — variance-band neighbourhood ("vband")

**For:** Claude Code. **Do this note first; note B depends on the infrastructure built here.**
**Repos:** `asishallab-group/Tensor-Omics`, branch `145-empirical-noise-model` (Fortran/Rcpp/R) and
`schroederaaron/empirical-noise-model-testing`, `main` (benchmark).
**Status of this note:** written without executing any Fortran or R. Line numbers below are approximate; search by name.
**Reference prototypes (Python, normal-noise toys, NOT the product):** `toy_variance_moderation.py` (`fit_eb`, `moderate`) and
`toy_variance_neighbourhood.py` (neighbour selection + exact tail count). Put them under `docs/prototypes/` of the benchmark repo.

---

## 1. Goal and hypothesis

TOX builds each gene's null from the residuals of its k nearest neighbours **in mean expression only**. If genes at the same
mean differ in noise level (gene-wise dispersion scatter), the pool is a scale mixture: quiet genes are tested against a null
that is too wide (lost power) and noisy genes against one that is too narrow (false positives). In the P+ benchmark with the
TCGA-estimated parameters, the gap to limma grows with the scatter sigma_d (AUPRC excess gap about -0.01 at sigma_d=0, -0.10 at
sigma_d=1.05, -0.2 at 2.1), and on real cohorts (Part R) TOX trails limma by -0.06 to -0.12 AUPRC excess.

**Model A** keeps the empirical, distribution-free null exactly as it is and only changes WHICH genes become neighbours:
a candidate gene is accepted only if its detrended log-variance lies within a two-sided band around the target's.
**Hypothesis to test:** this recovers most of the ranking gap at moderate/large n without changing the null's construction.

**Non-negotiable (from the project owner):** the null stays empirical and non-parametric: pairwise differences of neighbours'
raw residuals, exact tail count. No parametric null, no change to how residuals are scaled or compared.

---

## 2. Definitions (fixed)

Let `L = log2(max(x,0) + 1)` per sample, on the **same normalised input matrix TOX receives** (TPM, TPM-centred or TMM-CPM;
centring does not change within-group variance). Groups have n_c and n_t samples, `d = n_c + n_t - 2`.

* Pooled within-group variance of gene g: `s2_g = [ sum_case (L - mean_case)^2 + sum_ctrl (L - mean_ctrl)^2 ] / d`
* `l_g = log(s2_g + var_floor)`, `var_floor = 1e-8` (numerical guard only; log2-scale units; report how many genes it touches).
* Mean axis for the trend: `m_g = log2( mean over ALL samples of the LINEAR input + 1 )`.
* Trend: `S(m)` = robust LOESS of `l` on `m` (`stats::loess`, `family = "symmetric"`, `degree = 1`, `span = 0.3` as defaults; expose all three).
* **Matching variable (Model A): `u_g = l_g - S(m_g)`** (raw, unshrunk, detrended log-variance ratio to the trend).
* **Band (two-sided):** candidate j is eligible for target g iff `u_j` is finite and `|u_j - u_g| <= tau_var`.
* **Fixed k:** the pool is the first `k` eligible genes in order of mean-closeness (the existing `choose_index` left/right expansion). k = 50
  (the current `k_max`; in the P+ runs the median pool already sat at 50 genes for sigma_d <= 1.05).
* **tau_var rule:** `tau_var = c * sqrt(trigamma(d/2))` (the sampling SD of log s2 at df d: 0.80 at n=3, 0.53 at n=5, 0.34 at n=10 per group).
  Primary `c = 1.5`; sensitivity `c = 1` and `c = 2`. Rationale: a fixed absolute band is narrower than the noise at n=3 and generous at n=10.

LOESS is used to **detrend**, not to smooth the matching variable: a LOESS fit over mean expression gives every gene at the same mean the same
value and carries no gene-specific information. The gene-specific signal is the deviation from the trend (`u_g`).

---

## 3. Fortran / ABI changes (Tensor-Omics)

Files: `src/tox/tox_noise_model_exact.F90` (main work), `src/tox/tox_noise_model.F90` (bootstrap: accept + reject only),
`rcpp/tensoromics_functions.cpp` and `.R`. Follow the **`null_method` precedent**: both entry points keep an identical ABI; the module
that does not implement a mode accepts the argument and returns `ierr = ERR_INVALID_INPUT` for any non-default value.

### 3.1 New arguments (add to `compute_noise_pvalue_pipeline`, `..._helper`, the `bind(C)` entry `compute_noise_pvalues_pipeline_exact_c`, the Rcpp function `tox_compute_noise_pvalues_pipeline_exact_rcpp`, and the R wrapper `tox_compute_noise_pvalues_pipeline_exact`)

| Argument | Type | Meaning |
|---|---|---|
| `match_mode` | int | 0 = current behaviour (default); 1 = variance band |
| `tau_var` | real64 | half-width of the band (log-variance units); used only in mode 1 |
| `match_logvar(n_genes)` | real64 | matching variable `u_g`, indexed by the ORIGINAL gene index; used only in mode 1 |

New outputs (both modes; `-1` where not computed, always `-1` in mode 0): `n_scanned_own_case(n_genes)`, `n_scanned_own_control(n_genes)` =
number of sorted candidate slots examined while building each pool.

**Mode 1 semantics:** `k_max` is the fixed neighbour count `k`; `k_start`, `k_step` and `tau` are ignored (they must still pass the existing
positive-integer validators, so the R wrapper passes `k_start = k_step = k_max = k`). The R wrapper defaults must make mode 0 the default
(`match_mode = 0L, tau_var = 0, match_logvar = NULL`, with `NULL` replaced by `numeric(n_genes)` before the call).

### 3.2 New routine `gather_residuals_band_helper`

Same contract as `gather_residuals_helper` (pure, no allocation, writes into the caller's `pooled_residuals`, hard cap `max_pool_size`), plus:
`target_u`, `match_u(:)` (original-index vector), `k`, `tau_var`, and outputs `n_genes_pool`, `n_scanned`.

```
pos = find_closest_helper(target_mean, means_sorted, n_genes)        ! closest slot is the FIRST candidate
candidates are then visited in mean-closeness order via choose_index (left_cand = pos-1, right_cand = pos+1)
for each candidate slot idx:
    n_scanned += 1
    u_j = match_u(original_indices(idx))
    if finite(u_j) and abs(u_j - target_u) <= tau_var:
        append sorted_data%residuals_packed(1:n_residuals(idx), idx) to the pool   (respect max_pool_size)
        n_genes_pool += 1
    stop when n_genes_pool == k, or n_scanned >= scan_max, or both sides exhausted, or pool is full
```

* `scan_max = BAND_MAX_SCAN_FACTOR * k`, `BAND_MAX_SCAN_FACTOR = 20` (module constant; document it). This bounds how far in mean the
  neighbourhood may drift when the band is strict.
* **Short pools:** do NOT silently relax tau. If `n_genes_pool < k` keep the shorter pool; if `n_genes_pool < BAND_MIN_GENES`
  (module constant, start with 10) treat the gene as not computable (p = -1, same as the existing `n_pool < 10` gate). The existing residual
  gate `n_pool < 10` stays.
* **Self-inclusion:** the existing code includes the target gene in its own pool (it is the closest gene to its own mean). Keep that
  (distance 0 always passes the band). The Python toys excluded self; the effect is about 1/k. Note it in the doc comment.
* **Equal-variance assumption in SELECTION only:** the same `u_g` is used for the case pool and the control pool (pooled variance, df = n_c+n_t-2).
  The null itself still uses separate per-side pools and per-side 1/sqrt(n) scaling, so the "no equal-variance assumption" property of the
  null is untouched. (A per-side `u` was worse in the toy because df = n-1 is much noisier. Out of scope.)
* Non-finite `target_u` -> gene skipped (p = -1), counted.
* The own-scaling by `1/sqrt(n_replicates)` and `compute_pvalue` are **unchanged**.

### 3.3 Validation (`tox_errors` patterns)

`match_mode in {0,1}`; in mode 1: `tau_var > 0`, `k_max >= 1`, array length `n_genes`. Bootstrap module: accept the new arguments, return
`ERR_INVALID_INPUT` if `match_mode /= 0`.

### 3.4 Docs

FORD comments (`!>`, `!|`) for every new argument/routine; extend `docs`/README text that lists the pipeline arguments and the `null_method`
precedent; add a short section to the benchmark repo's `docs/noise_model_current_state.md` ("Variance-band neighbourhood").

---

## 4. Fortran tests (`test/mod_test_noise_model.f90`, follow existing assertions)

1. **Backward compatibility:** `match_mode = 0` p-values are bit-identical to the pre-change build on an existing fixture (store the values before editing).
2. **Parity:** `match_mode = 1` with `tau_var` huge (all genes eligible), `k = K` equals `match_mode = 0` with `k_start = k_step = k_max = K` and any tau
   (identical pool, identical p-values) as long as `max_pool_size` does not bind.
3. **Hand-built example** (about 12 genes, known variances): assert exactly which genes enter the pool for a given target and band, and the p-value
   against an independent tail count (brute-force O(n^2) in the test).
4. **Edge cases:** non-finite `match_logvar` entries are never neighbours; non-finite target skipped; `scan_max` stops the search; short pool kept; pool below
   `BAND_MIN_GENES` skipped; `n_scanned_*` correct; bootstrap module rejects `match_mode = 1`.
5. **Dimension/range validation** failures return the right `ierr`.
6. R-level test in `rcpp/test` (if that is where wrapper tests live): the wrapper defaults reproduce mode 0 exactly.

---

## 5. R helper (benchmark repo): `common/tox_variance_match.R`

```r
tox_match_logvar(L_case, L_ctrl, mean_lin, model = c("A", "B"), span = 0.3, degree = 1,
                 robust_loess = TRUE, var_floor = 1e-8, ...)   # model "B" is note B
# L_case, L_ctrl : SAMPLES x GENES matrices of log2(max(x,0)+1) of the normalised input
# mean_lin       : per-gene mean of the LINEAR normalised input over all samples
# returns list(u, tau_unit = sqrt(trigamma(d/2)), d, trend = S(m) per gene, n_floored, ...)
tox_resolve_tau(model, c, d, w = NULL)    # A: c * sqrt(trigamma(d/2))
```

For model "A": `u = l - S(m)` as in section 2. Return diagnostics (d, tau_unit, number of genes at `var_floor`, LOESS span/degree).

**Unit tests (testthat, `TCGA_test/tests/`):**
* `u` is mean-free: median(u) about 0 and no residual dependence on `m` (rank correlation with `m` near 0).
* With no heterogeneity (constant sigma^2 across genes) the spread of `u` matches the sampling SD `sqrt(trigamma(d/2))` within 10%.
* Invariance: multiplying all expression of a gene by a constant leaves `u` unchanged only on the log2(x+1) scale where expected; document the actual behaviour, do not assert something false.
* `tox_resolve_tau` returns 0.80, 0.53, 0.34 (times c) for n = 3, 5, 10 per group (d = 4, 8, 18).

---

## 6. Benchmark plumbing

* `TCGA_test/scripts/power_test.R`: add arms next to `PCFG$tox_arms` (use the existing `grid_arms`/parts mechanism around the `for (a in PCFG$tox_arms)` loops, lines about 620-650,
  so the new arms do NOT run in every part):
  ```r
  list(label = "TOX-log-vband",    engine = "exact", norm = 1L, null = 0L, input = "tpm", match = "A", c_tau = 1.5, k = 50L)
  list(label = "TOX-log-vband-c1", ..., c_tau = 1.0)      # sensitivity
  list(label = "TOX-log-vband-c2", ..., c_tau = 2.0)      # sensitivity
  ```
  Parts: **P+, D, R** (in Part R use the same `input` as the comparator `TOX-log|TPM-centred`). `run_tox_arm()` (about line 506): if `arm$match` is set, compute `u` with
  `tox_match_logvar()` from the arm's own normalised input, resolve `tau_var`, and call the exact pipeline with `match_mode = 1L`, `k_max = k`.
* Raw output (`power_raw/`): add per-gene columns `u`, `n_scanned_case`, `n_scanned_control` for these arms; add to the job manifest: `tau_var`, `c_tau`, `span`, `n_floored`.
* `power_report.R` / `common/power_metrics.R`: add the labels to the palette/linetypes; add summary columns `frac_short_pool` (pool < k genes), `n_scanned_median`,
  `frac_skipped`. **No change to existing metric definitions.**
* `TCGA_test/scripts/null_calibration.R` (TCGA null halves): inspect how arms are defined there and add the same arms with the same helper. **Read the script first; do not assume its structure.**

---

## 7. Evaluation protocol (shared with note B; **decision criteria are a proposal: the project owner confirms or edits them before any run**)

Order matters; tune nothing on the power results.

1. **Calibration first:** TCGA null halves, then P+ with pi1 = 0 (`p_le05_null`, any-false-call probability, same statistics the existing calibration outputs use).
   Choose `c` (and for B the analogous constant) **on calibration only**.
2. **Then power/ranking:** P+ (sigma_d grid including 0, n = 3, 5, 10), Part D (hypervariable genes: `fpr_hv05`), and **Part R** (primary evidence).
3. Report per metric: FDR at nominal 0.05 (zero-call rounds = 0), TPR at nominal and achieved FDR, AUPRC excess, |LFC|80, with paired differences against `TOX-log` and `limma`.
4. Always show the diagnostics: `frac_short_pool`, `n_scanned_median`, realised pool size.

**Proposed acceptance (edit as needed):** a variant is *unacceptable* if its null-half `p_le05_null` exceeds `TOX-log`'s by more than 0.01 (absolute) or its Part R FDR at nominal
exceeds `TOX-log`'s by more than 0.02. Among acceptable variants, it is an *improvement* if Part R AUPRC excess rises by at least 0.02 in at least 5 of 6 cohorts with no cohort losing more than 0.01.

---

## 8. Expected behaviour from the Python toys (sanity checks, NOT targets)

Single mean bin, normal noise, 10% DE, sigma_v = log-variance scatter (6 datasets per cell; differences below about 0.03 are noise). Two-sided raw band, tau = 0.7:

* sigma_v = 1.05, n = 10: AUPRC excess 0.747 vs 0.662 for current TOX (moderated t: 0.748); FDR 0.074 vs 0.078; TPR 0.66 vs 0.23.
* sigma_v = 1.05, n = 5: AUPRC 0.576 vs 0.413; FDR 0.22 vs 0.46 (still too liberal).
* sigma_v = 1.05, n = 3: AUPRC 0.340 vs 0.287; FDR about 0.5. **n = 3 stays hard.**
* **sigma_v = 0 (no heterogeneity): Model A costs calibration at small n** (n = 5: FDR 0.19 vs 0.04; n = 3: 0.55 vs 0.02; n = 10: 0.05 vs 0.03). Selecting neighbours on the sample variance of the
  same residuals that build the null makes the null mimic the target's own estimated variance (a plug-in effect). If the real-data null halves show this, it is a finding, not a bug.
* One-sided bands are much worse (FDR 0.38 at n = 5, 0.19 at n = 10): keep the band two-sided.

---

## 9. Risks, decisions embedded, out of scope

* **Defaults chosen by the author (change if wrong):** k = 50; `c` in {1, 1.5, 2}; LOESS `span = 0.3`, `degree = 1`, symmetric; `var_floor = 1e-8`; `BAND_MAX_SCAN_FACTOR = 20`; `BAND_MIN_GENES = 10`; self-inclusion kept.
* **Neighbourhood drift in mean** when the band is strict: monitor `n_scanned_median` and the realised mean span; extend `tox_diagnose()`-style checks if needed.
* **Hypervariable genes** (real data): a few genes with extreme `u` get few eligible neighbours; they should end up with short pools or be skipped. Check how many.
* **Out of scope:** bootstrap/blocked engines (stub only), the multi-axis module, per-side variance matching, cross-fitting the selection variance (select by the other group's variance), a studentised null,
  adaptive k. Mention cross-fitting as the first follow-up if calibration at small n is poor.
