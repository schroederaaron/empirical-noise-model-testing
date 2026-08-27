# Current State
Current State as of 11.08.2026

The development of the noise model is mainly complete. Two different models were implemented and are currently being tested.

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
