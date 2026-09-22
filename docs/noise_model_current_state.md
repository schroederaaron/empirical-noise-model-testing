# Current State
Current State as of 10.09.2026

The development of the noise model is mainly complete. Two interchangeable models are implemented behind one ABI — **`exact`** (sqrt-scaling, deterministic) and **`bootstrap`** (mean-level null) — and the bootstrap module additionally offers two **null constructions**, selected by `null_method`: the original **pooled** null (`0`) and the **gene-blocked** null (`1`, exactly enumerated where it fits). All are currently under test.

## Sqrt Scaling approach (model: "exact")
This model works as described in issue 145. For each gene it builds a residual pool in case and an independent residual pool in control from the mean-expression neighborhood (mean-centered, Bessel-corrected). The null is the distribution of the pairwise absolute differences `|r_case - r_control|` over every case×control residual pair — a direct measure of how much case-vs-control distance noise alone can produce. The p-value is the add-one-corrected fraction of those pairwise differences that reach or exceed the observed `|mean_case - mean_control|`.

The observed statistic is a difference of *means*, but the null is built from differences of *individual* residuals. An individual residual has standard deviation `sigma`, whereas a mean of `n_rep` replicates has standard deviation `sigma / sqrt(n_rep)`. So before the pairwise differences are formed, every residual is divided by `sqrt(n_rep)` (independently in case and control). This matches the null's variance to that of the observed mean-difference (`sigma^2_case / n_case + sigma^2_control / n_control`) exactly.

The scaling fixes the width but not the shape: the null still consists of differences of *individual* residuals, so it keeps the residual (heavy) tails, whereas a difference of means has its tails thinned by averaging (CLT). The null therefore sits a little wider in the tails than the statistic it scores, which makes the model mildly conservative — p-values biased slightly upward. Removing that residual bias is what the mean-bootstrap model below is for.

The advantage of this approach is that the null is computed exactly and deterministically — no sampling — via a sorted control pool and binary-search tail counting, so it is fast and reproducible.

## Mean-bootstrap approach (model: "bootstrap")

This is the model that resolves the over-conservatism of the sqrt-scaling approach. It builds the residual pools exactly the same way (a case pool and an independent control pool, both gathered from the mean-expression neighborhood, mean-centered and Bessel-corrected), but it constructs the null at the **mean level** instead of the individual-residual level.

For each gene we draw `N_BOOTSTRAP_DRAWS` (currently 10000) bootstrap samples. Each single draw:

1. resamples `n_rep` residuals **with replacement** from the case pool and averages them → one bootstrapped case mean,
2. resamples `n_rep` residuals with replacement from the control pool and averages them → one bootstrapped control mean,
3. takes the absolute difference of the two → one **null distance**.

The `n_rep` used per side is the actual per-gene replicate count (`sorted_*%max_resid_per_gene`, i.e. `n_samples`) — the same number of values the *observed* means were averaged over. The 10000 null distances form the null distribution, and the p-value is the add-one-corrected tail fraction:

```
p = (#{null_distance >= |observed|} + 1) / (n_boot + 1)
```

where the observed statistic is `|mean_case - mean_control|` on the residual scale.

Why this fixes the over-conservatism: the null is now a difference of two means-of-`n_rep`, so it carries the correct sampling variance `sigma^2_case / n_case + sigma^2_control / n_control` — matching the observed mean-difference directly, without the `sqrt(n_rep)` widening that the exact model has to live with. Each side is resampled from its own pool at its own `n_rep`, so there is **no equal-variance assumption** between case and control, and because the residuals are mean-centered the null is centered at zero shift (a true H0). Two boundaries fall out exactly: `|observed| = 0` gives `p = 1`, and an observed larger than any null distance gives `p = 1/(n_boot + 1)`.

The cost is that the null is sampled rather than computed in closed form, so it is slower than the exact model. The RNG is seeded once by `init_random(42)` at the start of the pipeline, so results are reproducible.

## Null construction: pooled vs gene-blocked (`null_method`)

The bootstrap module offers **two ways to build the mean-level null**, selected by the ABI parameter **`null_method`** (an integer, default `0`). The exact module accepts the argument so both C entry points keep an identical ABI, but **rejects any non-zero value** with `ierr = ERR_INVALID_INPUT` rather than silently ignoring it.

**`null_method = 0` — POOLED (the original behaviour).** Each draw resamples `n_rep` residuals **iid from the whole neighbourhood pool** and averages them. Because a mean-expression neighbourhood is mean-homogeneous but *not* variance-homogeneous, one draw can mix a residual from a quiet gene with one from a noisy gene — a combination no real gene's mean ever produces. The pooled residuals are therefore a **scale mixture**: the total variance is right, but the shape is wrong — a narrow core with heavy tails.

**`null_method = 1` — GENE-BLOCKED.** Each draw first picks one neighbour **gene**, then resamples `n_rep` residuals **within that gene**, so every null mean carries a **single coherent noise level**. This changes the null's *shape*, not its width, and it is the construction that addresses the raw arm's anti-conservatism (see `docs/raw_normalisation_diagnosis.md`).

Because resampling `n` of a gene's `n` residuals with replacement has only `C(2n-1, n)` distinct outcomes, the entire blocked null can be **written down exactly** as `n_genes_pool * C(2n-1, n)` weighted values per side and scored with the same sorted-array + binary-search tail count the exact model already uses — **no RNG and no `1/(n_boot+1)` floor**. Above the enumeration cap (`BLOCKED_MAX_ENUM_VALUES`, default 8000 — covers roughly `n_rep <= 5` at `k_max = 50`) it falls back to *sampling* the same blocked null.

Resolution (the p-value floor) therefore depends on the arm:

| construction | floor |
|---|---|
| exact | `1 / (n_pool_case * n_pool_control + 1)` |
| bootstrap, pooled (or blocked above the cap) | `1 / (N_BOOTSTRAP_DRAWS + 1)` |
| bootstrap, blocked & enumerated | `1 / (W_case * W_control + 1)`, `W = (n_pool / n_rep) * n_rep**n_rep` |

**Caveat:** `trim_frac` is **ignored** under `null_method = 1` — trimming sorts the pool in place, which destroys the per-gene block layout the blocked null reads, so the alloc layer forces it to 0. Trim + blocked is not a valid combination.

**RNG:** `random_number` is now called once per **chunk** of draws (a single array fill sized by `RNG_BUFFER_MAX`, allocated once and reused per gene) instead of once per draw — typically one RNG call per gene instead of `N_BOOTSTRAP_DRAWS`. The pooled arm is bit-for-bit unchanged by this (verified: `max |p_old - p_new| = 0`); the enumerated blocked path uses no RNG at all.

## The Bessel correction

The residual pool is built from mean-removed residuals `r_i = x_i - xbar` (linear) or `r_i = log2(x_i + c) - ghat` (log). Mean-removed residuals systematically **understate** the true spread: their variance is `sigma^2 * (n-1)/n`, not `sigma^2`, because one degree of freedom was spent estimating the mean. At small replicate counts this bias is large — at `n = 3` the residual standard deviation is only `sqrt(2/3) ~= 0.82` of the true `sigma`, an ~18% understatement.

If we sampled the null from an under-dispersed pool, the null would be too narrow and the model would be **anti-conservative** — it would call too many genes significant (inflated false-positive rate). To prevent this, every residual is scaled by the Bessel factor `sqrt(n / (n-1))` when the pool is built (see `prepare_sorted_data_helper`). This restores the pool to the unbiased `sigma`, so the bootstrap null has the correct width and the observed mean-difference is scored against a properly-scaled distribution. The correction matters most exactly where it is largest — the low-replicate datasets — which is the regime we care about.

## Variance stratification (removed)

Variance stratification has been **removed** from the code. The `own` null is now built directly from the whole gathered kNN pool (the "fallback" path that was already the only path in practice — see below). This section records what it was, why it never helped, and — since the code is gone — **how it was implemented**, so it can be rebuilt if the raw-normalization work ever needs it.

**The idea.** Subdivide the mean-expression neighborhood into **variance strata**: instead of comparing a gene against a null built from all nearby-mean genes, compare it against a null built only from genes that also share its *variance* regime. This targets the mixture case — same-mean genes that nevertheless have different variances (e.g. a subtype mixture within one cancer stage) — where a single pooled null over- or under-states the noise for individual genes.

**Why it made no measurable difference for the log normalization.** In log space the mean-variance relationship is largely stabilized (the whole point of the log/voom-style transform), so genes at the same mean already have approximately equal variance. A mean-expression neighborhood is therefore already variance-homogeneous, and there is nothing left to stratify — our tests confirmed no difference with stratification on vs. off under log normalization.

**Why it might still be useful for the raw normalization.** In raw/linear space the mean-variance trend is strong (variance grows with the mean) and heavy tails are common, so even within a narrow mean-neighborhood there can be genuine variance heterogeneity. Stratifying by variance could plausibly improve calibration there — the one scenario that would justify reintroducing it.

**Why its acceptance conditions were, in practice, impossible to meet.** Stratification only "activated" (used more than one bin) if a candidate binning passed three acceptance criteria. Two of them — `median(c_g) == 1` and `Pr(c_g > 2) < 0.1`, where `c_g` is the number of bins a gene spans — were trivially satisfied by construction: residuals were binned by their **source-gene mean**, and since all of a gene's residuals share that one mean, every gene fell entirely into a single bin, so `c_g == 1` always. That left the third as the only binding criterion: **every occupied bin must hold at least `STRATA_MIN_RESIDUALS_PER_BIN` (50) residuals**. But the binning variable is the source-gene mean while the neighborhood is a *k-nearest-neighbors set in mean space* — by construction mean-homogeneous, so the binning variable is nearly constant across the pool. Any split into ≥2 quantile bins either produced degenerate/near-empty bins or left almost all residuals in one bin, and the bins could not each independently reach the 50-residual floor. At low replicate counts the pool is also just too small: ~`k_max` (≈15) genes × `n_rep` residuals is ~45 at `n_rep = 3` — below the 50-floor even for a *single* bin. So the finer steps never passed and it always fell back to the single-bin whole-pool floor. Stratification-by-mean cannot separate variance regimes because the neighborhood it operates on has already been made mean-homogeneous by the neighborhood selection itself.

### How it was implemented (for future resurrection)

The layer lived, per side (case and control), between `gather_residuals_helper` and the p-value step, in **both** `tox_noise_model.F90` (bootstrap) and `tox_noise_model_exact.F90` (exact). It required a parallel `gene_id_per_residual` array from `gather_residuals_helper` recording, for each pooled residual, the sorted-gene slot it came from (built via `add_gene_id_to_pool_helper`).

- **Parameters.** `STRATA_BIN_COUNT_SCHEDULE` — candidate quantile bin counts, finest (5%/bin) to coarsest, with a **1-bin hard floor** last that is accepted unconditionally; `STRATA_N_SCHEDULE_STEPS` (its length); `STRATA_MAX_C_G_PROB_THRESHOLD` (= 0.1, criterion 2); `STRATA_MIN_RESIDUALS_PER_BIN` (= 50, criterion 3). Disabling was done by setting the schedule to `[1]`.
- **Five helper subroutines:**
  - `assign_residual_bins_helper` — assign each residual to an equal-*frequency* quantile bin, with cutpoints from `calc_percentile_helper` (linear-interpolation percentiles of the binning variable, i.e. the source-gene mean). Took a caller-supplied ascending `perm` so the O(n log n) sort was done once per pool, not per schedule step.
  - `locate_bin_helper` — binary-search the bin a value falls into, given the `n_bins+1` edges (values outside are clamped to the end bins).
  - `check_stratification_accepted_helper` — evaluate the three acceptance criteria for a candidate bin count, using work arrays `tmp_gene_min_bin` / `tmp_gene_max_bin` / `tmp_gene_seen` (an all-`.false.`-on-entry "seen" array it restored on return), `tmp_touched_gene_slots`, `tmp_c_g`, `tmp_bin_counts`.
  - `stratify_residuals_helper` — try each schedule step finest-first, calling assign + check, and return the coarsest-accepted (or the 1-bin floor) as `chosen_n_bins`, `chosen_bin_index_per_residual`, `chosen_bin_edges`, plus a `criteria_met` flag.
  - `select_stratum_for_target_helper` — given the accepted binning, copy out just the residuals in the bin containing the target gene's own mean.
- **Per-gene flow.** Stratify the case pool → select the target's stratum; same for control; then the p-value ran on the selected strata (bootstrap draws in the baseline model; sqrt-scale + exact tail count in the exact model), gated on each stratum having ≥ 10 residuals.
- **ABI diagnostic.** Two outputs, `chosen_n_bins_own_case` / `chosen_n_bins_own_control`, carried the chosen bin count per side **sign-encoded**: magnitude = bin count, sign = criteria met (+) vs coarse fallback (−), `-1` = not computed. These were surfaced through the Rcpp entry (`chosen_n_bins_own_case/control` list elements) and consumed by the R comparison/calibration scripts. They were removed with the layer.

## Residual-pool tail trimming (raw normalization)

A lighter-weight successor to stratification for the raw-normalization variance problem. Instead of trying to *partition* the neighborhood by variance (which the section above shows cannot work, because the mean-neighborhood is already mean-homogeneous), we simply **trim the tails of the pooled residuals**: after the kNN pool is gathered, sort it and drop the lower and upper `trim_frac` (default 5%) by value, keeping the central `1 − 2·trim_frac`. The null is then built from the trimmed pool exactly as before.

**The theory.** In raw/linear space the mean-variance trend is strong and heavy tails are common, so even a mean-homogeneous neighborhood can contain a handful of extreme residuals that widen (or, via a lone huge value, distort) the empirical null. If the pool *has* such artificial outliers, trimming removes them and tightens the null to the bulk of the noise; if it *doesn't*, the residuals are already a homogeneous spread, so the trimmed values sit close to the rest and almost nothing is lost. Either way the trimmed null is a more faithful estimate of typical noise. Under **log** normalization the mean-variance relationship is already stabilized and the tails are light, so trimming buys nothing there — it is therefore gated to raw only (`norm_method == 0`; the pipeline passes `trim_frac = 0` under log).

**How it is implemented.** A single shared helper `trim_pool_tails_helper(pool, n_pool, trim_frac)` in **both** `tox_noise_model.F90` and `tox_noise_model_exact.F90` (kept in sync). It is called on each side's pool immediately after `gather_residuals_helper` — before the `< 10`-residual gate and, in the exact model, before the sqrt-scaling — so both models score the trimmed central residuals. It sorts the pool (indirect `sort_real`) and keeps the central `n_pool − 2k` residuals, where `k = floor(n_pool · trim_frac)`. It no-ops when `trim_frac ≤ 0`, when `k` rounds to 0 (pool smaller than `1/trim_frac`, so low-`n_rep` pools are untouched), or when trimming would empty the pool (`trim_frac ≥ 0.5`); the reported `neighborhood_size_*` is the post-trim count. Unlike stratification this needs **no** `gene_id_per_residual` plumbing and **no** ABI diagnostics — just one real parameter, `trim_frac`, threaded after `tau` through the C entry, the Rcpp dispatcher, and the R wrappers (default `0.0`). Both calibration scripts expose it as an A/B arm (`null_calibration.R` → `TOX-raw-trim05`; `calibration_test.R` → a `_trim` companion for each raw arm), so trimmed-vs-untrimmed raw calibration appears in one sweep.

## Multidimensional (multi-axis) noise model — `noise_model_md`

A third module, `src/tox/tox_noise_model_md.F90` (`module noise_model_md`), generalises the model from a single case/control contrast to `d` **independent axes** (tissues, stages, conditions). A gene is **one vector**; the axes are not tested individually inside this test. It is separate from the two scalar modules because the ABI carries an extra dimension on almost every array.

### Statistic — two norms, two jobs

Per gene `g` and axis `a = 1…d`, `beta_hat(a) = mean_case(g,a) − mean_ctrl(g,a)` (with two groups and no covariates this is the OLS estimate, i.e. per axis identical to today's `obs_own`). Two norms are computed, and both are calibrated because the weights are applied to observed and null alike:

```
D     = sqrt( sum_a beta_hat(a)^2 )                    effect size, log2FC units
D_std = sqrt( sum_a beta_hat(a)^2 / Var_null(a) )      test statistic, unitless
```

They answer different questions — how large the response is, versus how confident we are that there is one — exactly as `logFC` and `t` do in limma. `D` is reported next to each gene; **`D_std` drives the p-value and the ranking**. With independent axes the null covariance is diagonal, so `D_std` **is** the full Mahalanobis statistic: no covariance matrix, no Cholesky, no ridge.

`D` is upward biased as a magnitude (`E[D] ≈ σ·√d ≠ 0` even under H0), so `d_sq_null_mean` is accumulated alongside it and the bias-corrected magnitude is `D²_adj = max(0, D² − d_sq_null_mean)`. Rank on the p-value regardless.

### Null — a product of residual vectors, per-side `1/√n`

Per gene, per axis, the residual pools are gathered by the existing adaptive kNN in mean space — `2d` gathers per gene, `gather_residuals_helper` unchanged. One null draw takes **one residual per axis per side**:

```
for a = 1..d:
    beta_null(a) = eps_case(a)/sqrt(n_case(a)) − eps_ctrl(a)/sqrt(n_ctrl(a))
D_std_null      = sqrt( sum_a beta_null(a)^2 / Var_null(a) )
```

**The scaling is per side, not per axis.** Each side is divided by the root of *its own* replicate count before the difference is taken; with `n_case = [4,6,7,6,5]` against `n_ctrl = [3,5,5,9,3]` each of the ten numbers enters its own term, and a single per-axis factor would be right only when both sides have equal `n`. The scaling is needed at all because the observed statistic is a difference of means over `n` replicates while a raw residual difference is at single-observation scale — without it the null is inflated by roughly `√n` and the test is strongly conservative. In the implementation the pool is scaled **once per gene** rather than each draw, which is exactly what `noise_model_exact` does in 1-D.

`Var_null(a)` then follows in closed form from the pools, with no sampling:

```
Var_null(a) = Var(pool_case(a))/n_case(a) + Var(pool_ctrl(a))/n_ctrl(a)
```

which on the already-scaled pools is just the sum of the two pool variances. It is returned per gene per axis as `var_null`, so the weights are inspectable rather than implicit. A degenerate axis (zero pool variance on both sides) gets **weight 0** rather than an infinity: it drops out of `D_std` on the observed and null sides alike, and the remaining axes decide the p-value.

Because the axes share no samples, drawing each axis independently reproduces the true joint null — that is the whole content of the scope assumption below. The `∏_a S_a` outcome space is never materialised; it is the size of an implicit measure, and it only decides which of the three paths runs.

### Three paths to the p-value, reported per gene

`S_a = n_pool_case(a) · n_pool_ctrl(a)` is every pairing of one case-side residual with one control-side residual. `method_used` records which path each gene took.

| `method_used` | path | when |
|---|---|---|
| 0 | **enumerated** — walk the joint space, exact tail count, no floor | `∏_a S_a ≤ enum_max_product` |
| 1 | **systematic** — fixed coprime stride per pool (the default sampler) | otherwise, `sampling_mode = 0` |
| 2 | **RNG** — classical Monte Carlo | otherwise, `sampling_mode = 1` |
| 3 | reserved for the convolution null | never produced yet |

**Enumeration is in practice the `d = 1` regression path and little else.** `S_a` grows *quadratically* in the neighbourhood size, so larger `k` makes enumeration worse, not better: at `k = 15`, `n_rep = 3` the joint space is 2 025 at `d = 1` but already 4.1 × 10⁶ at `d = 2` and 8.3 × 10⁹ at `d = 3`. The branch is kept because that is where `d = 1` lives — it reproduces `noise_model_exact`'s `pvalues_own` bit-for-bit. Do not plan around it for multidimensional use; the threshold is evaluated on the **product**, not per axis, because each axis looks harmless in isolation.

**Systematic sampling** walks each pool with a fixed stride chosen coprime to that pool's size, so the walk visits every entry before repeating and the marginals stay exactly balanced. It is deterministic and reproducible, needs no RNG and no seed, parallelises without stream splitting, and has lower variance than random draws at the same `M`. The caveat for the method text: `(1+B)/(M+1)` is justified as a *valid* Monte-Carlo p-value under **random** draws; under systematic sampling it is instead an approximation of the exact enumerated tail, and the error becomes a coverage error rather than a sampling error. Both samplers are implemented so they can be compared on the mock null.

`p = (1 + #{D_std_null ≥ D_std_obs}) / (M + 1)` on the sampled paths, and `(count + 1)/(∏_a S_a + 1)` on the enumerated one.

### Resolution

`M` (`n_draws_max`) is set by BH over `G` genes, **not by `d`**: the estimate is of a one-dimensional tail probability, `SE(p̂) = √(p(1−p)/M)`, with no dimensional surcharge. It is a user parameter defaulting to **2 × 10⁴**. A gene at the floor `1/(M+1)` is rejected only once its rank reaches `p_floor·G/q`, so at `G = 20 000`, `q = 0.05` roughly 20 genes must be tied at the floor before any is rejected; `M = 4 × 10⁵` is where the floor stops binding altogether.

**Besag–Clifford sequential stopping** (`n_exceed_target = c > 0`) samples until `c` exceedances (`p = c/m`) or `M` draws; expected draws `≈ min(M, c/p)`, so a 40× increase in `M` costs about 1.5× in average draws. It is a *different* estimator with granularity `c/m` and its interaction with BH is unverified, so it ships **off by default** and the fixed-`M` arm stays the reference.

### Scope assumption — enforced in R

Every sample must belong to exactly one axis: no plant, subject, extraction batch or control group shared between axes. That is precisely the condition under which the Cartesian-product null is correct, and it is why there is no `null_coupling` argument. It is metadata, not something the numbers reveal, so the guard lives in the R wrapper: `tox_compute_noise_pvalues_pipeline_md()` compares the sample identifiers (rownames) of **all** `2d` replicate matrices — a sample that is a case on one axis and a control on another violates the assumption just as badly — and **errors** if any identifier appears more than once, if a matrix has duplicate rownames, or if a matrix carries no usable identifiers at all. `check_axis_samples = FALSE` waives it explicitly. The helper is `tox_md_assert_disjoint_samples()`.

### Entry points and arguments

```fortran
subroutine compute_beta_md(replicates_case_packed, n_rep_case_per_axis, &
                           replicates_control_packed, n_rep_control_per_axis, &
                           n_genes, n_axes, norm_method, beta_obs, ierr)

subroutine compute_noise_pvalue_pipeline_md( &
    means_case, replicates_case_packed, n_rep_case_per_axis, &
    means_control, replicates_control_packed, n_rep_control_per_axis, &
    beta_obs, compute_pvalue_own, beta_mode, beta_centre, &
    n_genes, n_axes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
    null_method, sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
    seed, max_pool_size, &
    pvalues_own, d_obs, d_std_obs, d_sq_null_mean, method_used, n_draws_used, &
    neighborhood_size_case, neighborhood_size_control, var_null, &
    beta_check_cor, beta_check_mad, delta_hat, &
    n_genes_with_pvalue, ierr)
```

`compute_beta_md` takes the **packed replicates**, not the means: on the log2 branch the statistic that is scale-coherent with the residual pools is the difference of **Fréchet** means, `mean_i log2(x_i + c)`, and by Jensen's inequality that cannot be recovered from the linear-space mean. On the linear branch the two coincide. The `means_*` arrays keep their real job — the coordinate the kNN neighbourhood is matched in — and are not involved in the statistic.

**Data layout.** Ragged replicate counts per axis and per side are required, so there is no rectangular `(n_rep, n_genes, d)` array and no padding (padding would put sentinel values into a matrix whose consumers have no concept of sentinels):

```
replicates_packed(sum_a n_rep_a · n_genes)   ! axis blocks, each (n_rep_a, n_genes)
n_rep_per_axis(n_axes)
means(n_genes, n_axes)
```

**Samples-fastest within each block**, matching `prepare_sorted_data_helper`, which reads whole contiguous columns. Each axis is handed to that routine as an ordinary 2-D matrix by **pointer bounds remapping** — no copy, and not sequence association (these are module procedures with explicit interfaces, where a rank mismatch is a hard error). R matrices are already column-major in this layout, so packing on the Rcpp side is a `memcpy`.

| argument | shape | meaning |
|---|---|---|
| `beta_obs` | `(n_axes, n_genes)` | gene-contiguous; filled when `beta_mode = 0`, returned **centred** when `beta_centre = 1` |
| `beta_mode` | scalar | `0` compute internally (default, scale-coherent); `1` caller-supplied |
| `beta_centre` | scalar | `0` none (default); `1` subtract the per-axis median `beta_hat` across genes |
| `null_method` | scalar | ABI parity only; **must be 0** (see below) |
| `sampling_mode` | scalar | `0` systematic stride (default); `1` RNG |
| `n_draws_max` | scalar | `M`, user-settable, default 2×10⁴; floor `1/(M+1)` |
| `n_exceed_target` | scalar | Besag–Clifford stopping; `0` disables (default) |
| `enum_max_product` | scalar, int64 | enumerate while `∏_a S_a ≤` this; `0` forces sampling |
| `d_obs`, `d_std_obs` | `(n_genes)` | the two norms |
| `d_sq_null_mean` | `(n_genes)` | for the bias correction |
| `method_used` | `(n_genes)` | `0` enumerated, `1` systematic, `2` RNG, `3` convolution, `-1` not tested |
| `n_draws_used` | `(n_genes)` | achieved resolution |
| `var_null` | `(n_axes, n_genes)` | per-axis null variance; the weights, exposed |
| `neighborhood_size_*` | `(n_axes, n_genes)` | per-axis pool sizes, so a thin axis is visible |
| `beta_check_cor`, `beta_check_mad` | `(n_axes)` | supplied vs recomputed `beta_hat`; exactly `1` / `0` in internal mode |
| `delta_hat` | `(n_axes)` | estimated composition shift, reported even when centring is off |

`norm_method`, `k_start`, `k_step`, `k_max`, `tau`, `trim_frac` and `max_pool_size` keep their current meaning, applied per axis. `trim_frac` is raw-normalization only, as in the scalar modules.

**`null_method` is accepted only as 0.** Under this module's construction a draw takes one residual per side, and each neighbour gene contributes exactly `n_rep` residuals to the pool, so "pick a gene uniformly, then a residual within it" is the *same distribution* as "pick a residual uniformly from the pool". A gene-blocked argument would be a no-op that merely looked meaningful, so it is rejected with `ERR_INVALID_INPUT` rather than silently ignored — the same stance `noise_model_exact` takes.

**Thin axes drop the whole gene.** If any axis's pool falls below 10 residuals, or its `beta_hat` is non-finite, the gene is dropped entirely (`pvalues_own = -1`, `method_used = -1`). Dropping only the offending axis would make `D` incomparable between genes and break the magnitude ranking. All axes are gathered and their sizes recorded before that decision, so the failing axis stays identifiable. In practice the gate is not a live concern: `k_start`/`k_step`/`k_max` count **genes**, and 15 genes × 3 replicates already gives 45 residuals.

**External `beta_hat` is surfaced, not refused.** With `beta_mode = 1` the internal estimate is still computed and reported against the supplied one: correlation ≈ 1 with a large MAD is a **units** problem (edgeR's `coefficients` are on the natural-log scale, a factor of `ln 2`, and shrunk by default); degraded correlation is shrinkage or non-linearity. Only limma-voom is principled here.

### ABI and R entry

C entry `compute_noise_pvalues_pipeline_md_c` — deliberately **not** ABI-compatible with the two scalar entry points. It takes two extra arguments, `n_packed_case` / `n_packed_control`, the declared lengths of the packed buffers, cross-checked in Fortran against the value derived from `n_rep_*_per_axis` (`ERR_DIM_MISMATCH` on a mismatch, rather than reading past the caller's buffer). `enum_max_product` crosses as a `long long`, since the joint space passes 2³¹ at `d = 2` with quite ordinary neighbourhoods.

Rcpp binding `tox_compute_noise_pvalues_pipeline_md_rcpp`; R wrapper `tox_compute_noise_pvalues_pipeline_md()`, which takes **lists of per-axis matrices** (`samples × genes`), derives the per-axis means when not supplied, runs the sample-identifier guard, and returns `pvalues_own`, `d_obs`, `d_std_obs`, `d_sq_null_mean`, `method_used`, `n_draws_used`, `beta_obs`, the three `(n_axes × n_genes)` matrices (`neighborhood_size_case`, `neighborhood_size_control`, `var_null`), `beta_check_cor`, `beta_check_mad`, `delta_hat`, `n_rep_*_per_axis`, `n_success` and `ierr`.

Caller-side and still to be done in R: direction `beta_hat/D`, per-axis contribution shares `beta_a²/D²`, and BH over **tested genes only**.

### Normalisation

To be settled **empirically**, counts vs TPM, on the same dataset with everything else held fixed — not argued from the literature. The counts arm normalises in R (TMM `calcNormFactors`, or median-of-ratios); the TPM arm uses `beta_centre` on the Fortran side. Naming, so the comparison is labelled correctly: what `beta_centre` does is **median centring of log-ratios**, not TMM — TMM's weights come from count magnitudes, which do not exist on TPM. Same estimator family, more robust variant on counts. **`beta_centre` must be off whenever the R side has already normalised**, or the data are centred twice.

Why correct at all: on TPM `log2FC = beta − Δ` with `Δ` identical for every gene, and the null — built from *within-group* residuals — cannot contain information about a *between-group* offset. Running uncorrected is not the assumption-free option; it is the assumption `Δ = 0`. Median centring assumes only that the median gene is unchanged, which is strictly weaker. `delta_hat` is reported either way.

### Unit tests

`test/mod_test_noise_model_md.f90`, suite `noise_model_md`, 10 tests: the bit-exact `d = 1` enumerated reduction to `noise_model_exact` (both normalization branches, and with `n_case = 4` against `n_ctrl = 5` so it also pins the **per-side** scaling); `compute_beta_md` at `d = 1`; both norms against a hand-computed fixture; the closed-form `Var_null` against the enumerated null's own second moment; the plan's ragged `n_case = [4,6,7,6,5]` / `n_ctrl = [3,5,5,9,3]` fixture with an exactly-known residual variance; `enum_max_product` routing and `method_used`; systematic vs RNG agreement; `beta_centre` median recovery; degenerate input (thin axis, zero-variance axis, `D = 0`, non-finite `beta_hat`, and the validation error codes); and an off-scale supplied `beta_hat`.
