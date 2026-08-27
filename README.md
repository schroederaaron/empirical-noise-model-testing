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
