# Why the raw normalisation is anti-conservative — what the existing results already prove, and what is still missing

Status as of 09.09.2026. Sections 1–3 are re-analyses of results already in this
repository (`TCGA_test/results/null_calibration_9414.out`, 2,688 summary rows:
8 cancers × 5 cohorts × 9 kNN configs × 4 replicate arms × 100 splits). Section 4
is a controlled simulation run for this note. Section 5 is what still has to be
measured on TCGA. Nothing here is asserted beyond what the cited numbers support.

---

## 0. Short answer

The raw arm's anti-conservatism is **not** a neighbourhood-size problem and **not**
a null-width problem. Both of those are ruled out by the data we already have. The
evidence points instead at the **shape** of the raw null: in linear space a
mean-neighbourhood is mean-homogeneous but not variance-homogeneous, so the pooled
residuals are a *scale mixture* — correct total variance, but a **narrow core with
heavy tails**. A moderate observed distance lands too far out in that core
(anti-conservative at α = 0.05) while an extreme one is still covered
(conservative at α = 0.01). That is exactly the pattern in the data, and it is the
one pattern a width error cannot produce.

This is a well-supported hypothesis, not yet a confirmed one. Section 5 lists the
per-gene measurements that would confirm or refute it, and the scripts that make them.

---

## 1. Three facts from the existing run

Reproduce with `Rscript TCGA_test/scripts/reanalyse_null_calibration_sweep.R`.

### 1.1 `k` cannot be the cause

Paired within (cancer, stage, replicate arm), moving from the smallest to the
largest kNN config:

| method | k8_30 → k20_50 | median excess (inflation − 1) |
|---|---|---|
| TOX-raw | **−0.093** (median −0.086, 128/128 negative) | **+0.641** |
| TOX-log | −0.109 | +0.052 |

Larger `k` helps, consistently, and by about **one seventh** of what would need
explaining. `k_max` 30 → 50 moves raw inflation by 0.005. Retuning the kNN
parameters cannot close a gap of 0.64. **Hypothesis (3) is excluded as the primary
cause.** (It is not excluded as a marginal improvement.)

### 1.2 The null is not too narrow — it is the wrong shape

Median over all cells:

| method | FPR@.05 | inflation | FPR@.01 | inflation | median p |
|---|---|---|---|---|---|
| TOX-raw | 0.0821 | **1.64** | 0.0120 | **1.20** | 0.423 |
| TOX-log | 0.0526 | 1.05 | 0.0132 | 1.32 | 0.505 |
| limma | 0.0486 | 0.97 | 0.0090 | 0.91 | 0.493 |
| edgeR | 0.1020 | 2.04 | 0.0335 | 3.36 | 0.420 |

and restricted to `k20_50`, by replicates per group:

| n/group | raw infl@.05 | raw infl@.01 | log infl@.05 | log infl@.01 |
|---|---|---|---|---|
| ≤12 | 1.45 | **0.75** | 1.06 | 1.40 |
| 13–25 | 1.61 | 1.04 | 1.01 | 1.29 |
| 26–50 | 1.70 | 1.31 | 0.98 | 1.19 |
| 51–100 | 1.79 | 1.45 | 0.99 | 1.20 |
| >100 | **1.90** | 1.74 | 0.97 | 1.15 |

At small n the raw model is **anti-conservative at 5% and conservative at 1% at the
same time**. No error in the null's width can do that — a null that is too narrow
inflates both levels, and inflates the 1% level *more*. Compare edgeR, whose defect
*is* width-like: 2.04 at 5% and 3.36 at 1%, the ordering a scale error produces.

Fitting the null as `λ · (standardised t_ν)` against a Gaussian observed statistic,
identified from those two FPRs:

| n/group | λ (scale) | ν | implied excess kurtosis | quantile ratio |
|---|---|---|---|---|
| ≤12 | 0.899 | 6.7 | 2.21 | 1.489 |
| 26–50 | 0.863 | 8.7 | 1.27 | 1.441 |
| >100 | 0.837 | 9.7 | 1.04 | 1.426 |
| TOX-log (all n) | 0.99–1.01 | — | — | 1.27–1.28 |

Gaussian reference: λ = 1, ratio 1.314, excess kurtosis 0. So the raw null is
**~10–16 % too WIDE overall** (λ < 1) and clearly **leptokurtic**; TOX-log sits at
λ ≈ 1 with a ratio slightly *below* Gaussian. The defect is shape, not scale.

> This is a two-parameter descriptive summary of where a p-value distribution
> departs from uniform, not a fitted generative model. Its use is that scale and
> shape leave opposite fingerprints, so it identifies which one is broken.

### 1.3 The failure grows with replicates — and only for raw

Raw inflation rises 1.45 → 1.90 from n ≤ 12 to n > 100; `ad_A2` rises 386 → 725;
median p falls 0.443 → 0.388. TOX-log is flat to slightly improving (1.06 → 0.97).

A **fixed** scale error is n-invariant: observed and null both shrink as
`1/√n`, so the ratio does not move. Something that gets *worse* with n means the
two are converging to different shapes. That is what the model does by
construction: the observed statistic is a difference of **means** and Gaussianises
with n (CLT), while the null is built from **individual** residuals, whose shape is
fixed however large n gets. The `1/√n` scaling matches the width and leaves the
shape mismatch untouched — as `noise_model_current_state.md` already states. What
that note gets the wrong way round is the *sign*: it reasons only about the far
tail, where a heavy-tailed null is indeed conservative. At α = 0.05 we are in the
**core**, where the same leptokurtic null is too narrow.

---

## 2. The hypothesis

In raw (linear) space genes at the same mean still have very different sds. The
kNN neighbourhood conditions on the **mean** only, so the pooled residual pool
mixes genes of different scale. A scale mixture has the right total variance but a
peaked core and heavy tails. Under log the mean–variance relation is stabilised,
the pool is near single-scale, and the mismatch disappears.

This is the same structural fact recorded in `noise_model_current_state.md` under
"Variance stratification (removed)" — that a mean-neighbourhood is mean-homogeneous
by construction and therefore cannot be stratified by mean. The note treats that as
the reason stratification *failed*; it is also the reason the raw null is shaped
wrongly.

---

## 3. Verdict on the three candidate causes

| # | hypothesis | verdict | basis |
|---|---|---|---|
| 1 | τ too low → neighbourhoods too small → null too narrow | **Not supported as stated.** The null is ~10–16 % too **wide** (λ = 0.84–0.90), not too narrow. The *variance-structure* half of the hypothesis is right; the *width* half has the sign backwards. | §1.2 |
| 2 | observed distances too large relative to residuals | **Untested.** Nothing measured so far separates it. Part 3 of the new diagnostic script is built for exactly this. | — |
| 3 | mis-calibrated kNN parameters | **Excluded as primary cause.** Full k range moves inflation by 0.09 of a 0.64 excess. | §1.1 |

---

## 4. Controlled check

`simulate_mechanism()` (Part 5 of `raw_anticonservative_diagnosis.R`) builds an H0
cohort with right-skewed means and per-gene CV drawn with dispersion `cv_het` —
i.e. `cv_het` *is* "how much variance heterogeneity survives at a fixed mean".
Everything else held fixed, 1,500 genes, exact null:

| cv_het | n | norm | pool excess kurtosis | infl@.05 | infl@.01 |
|---|---|---|---|---|---|
| 0.00 | 10 | raw | 0.65 | **0.83** | 0.33 |
| 0.00 | 40 | raw | 0.94 | **0.63** | 0.67 |
| 0.35 | 40 | raw | 4.73 | 1.25 | 0.80 |
| 0.70 | 10 | raw | 9.64 | 1.39 | 1.13 |
| 0.70 | 40 | raw | **29.5** | **1.63** | 1.20 |
| 0.70 | 40 | log | 3.64 | 1.21 | 2.00 |

With heterogeneity switched **off**, the raw exact model is *conservative* — the
sqrt-scaling result the design note predicts. Switching it on, and nothing else,
takes pool kurtosis from 0.9 to 29 and inflation@.05 from 0.63 to 1.63 while
inflation@.01 stays at 1.2, and the log arm stays far closer to 1. That is the
real-data signature, produced from a single controlled cause.

A mean-**bootstrap** null on the same pools gives infl@.05 = 1.20 (vs 1.63 exact)
at cv_het = 0.7, n = 40, but infl@.01 = 2.2 (vs 1.2). It repairs the core and not
the tail — as expected, since resampling cannot add tail mass the pool does not
contain. It is a partial fix, not a fix.

**This is a simulation, not TCGA.** It shows the mechanism is *sufficient* to
produce the observed pattern. It does not show it is what is happening in TCGA.

---

## 5. What is missing, and what now measures it

Everything above is derived from aggregate p-value distributions. The mechanism is
a statement about **individual residual pools**, and the Fortran returns only
`pvalues_own` plus three neighbourhood *sizes* — the pool's width, shape and
composition, and each gene's own variance, are computed inside
`gather_residuals_helper` and discarded.

`common/tox_null_reimpl.R` is a faithful R port of the exact `own` path
(`prepare_sorted_data` → `gather_residuals` → trim → `1/√n` → exact pairwise
p-value), instrumented with those quantities. It needs no Fortran change and no
rebuild. It ships with `validate_against_fortran()`, which runs both
implementations on the same input; both new scripts run that gate first and label
their results provisional if it fails. Verified here against a naive
re-implementation: exact p-values match brute force to 1e-14, and pool membership
and stop reason match on 600/600 genes across τ-binding and non-binding regimes.

### New per-gene quantities

- `rho = sd_own / sd_pool` — does the pool describe *this* gene's noise? Its spread
  across genes is the variance heterogeneity a mean-neighbourhood failed to remove.
  Reported against the estimation-noise floor `1/√(2(n−1))`, so the excess is
  separable from sampling noise in `sd_own`.
- `kurt_pool` — the pool's excess kurtosis: the scale mixture, measured directly.
- `sd_log_nb_sd` — spread of the per-neighbour sds inside one pool.
- `z_pool` vs `z_own` — the observed distance scored against the pooled null, and
  against the gene's own replicates. **This is the test that separates hypothesis
  (2).** If `z_own` is calibrated and `z_pool` is not, the observed distances are
  fine and the neighbourhood is at fault. If both are miscalibrated, the observed
  statistic itself is the problem and no neighbourhood change will help.
- `stop_reason`, `n_genes_pool`, `mean_span`, `frac_below` — whether τ is even the
  binding constraint, and whether the neighbourhood is one-sided in mean.
- `p_floor = 1/(n_a·n_b + 1)` — the smallest p a gene can receive.

### Predictions, stated so they can fail

| | prediction | falsifies what if wrong |
|---|---|---|
| P1 | `kurt_pool` (raw) ≫ (log), rising with the per-pool sd spread | the mechanism |
| P2 | inflation@.05 rises with `kurt_pool` across gene quintiles | the mechanism (weakly — see caveat in the script: neighbouring genes share pools, so rows are far from independent, and within the raw arm kurtosis varies little) |
| P3 | `z_own` far better calibrated than `z_pool` under raw | if both fail → hypothesis (2), and this whole line is wrong |
| P4 | `frac_stop_tau` small at τ = 0.1; no τ in {0.02 … 1.0} fixes it | hypothesis (1) revived if some τ does fix it |
| P5 | **tail trimming makes α = 0.05 worse, not better** | the mechanism, and it settles the open ToDo directly |

P5 is worth stating plainly: the open ToDo *"check whether a [5,95] percentile cap
improves the raw model"* is predicted to come back **negative**. Trimming removes
the pool's tails, which narrows an already-too-narrow core. It treats the one part
of the distribution that is already conservative and worsens the part that is not.

### Scripts

| script | what it does |
|---|---|
| `common/tox_null_reimpl.R` | instrumented R port + `validate_against_fortran()` + `tox_neighbours()` |
| `TCGA_test/scripts/reanalyse_null_calibration_sweep.R` | reproduces §1 from committed output; seconds, no data needed |
| `TCGA_test/scripts/raw_anticonservative_diagnosis.R` | Parts 0–5: fidelity gate, pool composition, scale/shape, `z_own` vs `z_pool`, τ/k/trim response, controlled simulation |
| `TCGA_test/scripts/pvalue_tail_methods.R` | analytic vs counted tails: saddlepoint / GPD / pairs / bootstrap, and the pool-uncertainty limit (§7b) |
| `TCGA_test/scripts/null_construction_comparison.R` | five null constructions on the same pools: calibration, resolution, BH false positives (§7c) |
| `TCGA_test/scripts/null_calibration.R` | *(modified, +31/−10)* arms now carry a `model` field; adds `TOX-raw-boot`, `TOX-log-boot`, `TOX-raw-trim02`, `TOX-raw-trim05` |

The two new calibration arms are the cheapest decisive test available: the
bootstrap null is the *same pools* with the CLT shape. If `TOX-raw-boot` recovers
α = 0.05 and not α = 0.01, the diagnosis holds. If it fixes both, the diagnosis is
incomplete. If it fixes neither, it is wrong.

---

## 6. If the diagnosis holds

Ordered by how much they change:

1. **Keep log as the default for raw-scale data.** Already the evidence-based
   choice (12/12 calibrated); this note explains *why* rather than changing it.
2. **Give the null the observed statistic's shape** — the mean-bootstrap null,
   raw only. Partial by the simulation: fixes the core, not the tail.
3. **Remove the mixture instead of the tails.** Rescale each neighbour gene's
   residuals to unit sd before pooling, then re-scale the pool by the target
   gene's own sd. This pools *shape* across the neighbourhood while taking *scale*
   locally, which is what the neighbourhood can and cannot estimate respectively.
   Untested — it is a design proposal, and it reintroduces dependence on a
   per-gene sd estimated from few replicates, which at n = 3 is itself poor.
4. **Do not adopt trimming for calibration reasons** unless P5 fails.

---

## 7. Separate deliverable — BH survivors with full statistics

`TCGA_test/scripts/tox_significant_genes_report.R` reports every gene TOX-log still
calls after BH, with p-values from the **production Fortran pipeline** (via
`compute_noise_pvalues()`, so the calls are the project's own) and every added
statistic from the validated R port:

- **expression**: mean / median / sd / var / CV / quartiles / range in both groups,
  on the TPM scale and the log2 scale, plus log2FC, the observed statistic and
  direction, and the sample counts.
- **the null it was scored against**: `sd_pool`, the gene's own residual sd, `rho`,
  `kurt_pool`, `sd_null`, `z_pool`, `z_own` — both sides.
- **neighbourhood**: pool sizes in residuals and in genes, mean span, `frac_below`,
  stop reason; plus a companion long table with one row per (survivor, neighbour
  gene) carrying the neighbour's id, mean, mean-distance and residual sd.
- **inference**: raw p, BH q, rank, and the resolution floor.

### Read this before interpreting any q-value

The exact p-value is a pairwise tail count, `p = (count + 1)/(n_a·n_b + 1)`, so the
**smallest p a gene can receive is `p_floor = 1/(n_a·n_b + 1)`**, set by its own
pool sizes. A gene at that floor is **censored** — its true p may be far smaller —
and is flagged `at_floor`.

**The floor does not make BH impossible — it makes rejection quantised.** If `r`
genes are tied at the floor, BH rejects the whole block as soon as
`p_floor ≤ q·r/G`, i.e.

> **r ≥ p_floor · G / q**

Verified directly against `p.adjust(..., "BH")`: at n_rep = 3, k = 50 genes
(150 × 150 = 22,500 pairs, `p_floor` = 4.44e-5), G = 12,000, q = 0.05 the predicted
threshold is r = 11, and R rejects 0 genes at r = 10 and exactly 11 at r = 11.
So a *cohort* of ≥ 11 simultaneous floor genes fires; fewer fire nothing. What is
needed for a **single** gene to be rejectable on its own:

| n_rep | q = 0.05 | q = 0.01 |
|---|---|---|
| 3 | k ≥ **164** genes (492 residuals) | k ≥ 366 |
| 5 | k ≥ 98 | k ≥ 220 |
| 10 | k ≥ 49 | k ≥ 110 |

Under H0 the expected number of floor genes is `G · p_floor` ≈ 0.53 at k = 50,
n_rep = 3, so r ≥ 11 essentially never occurs by chance. `hits_FDR05` = 0.0000
across all 12 calibration groups is therefore the **expected and conservative**
consequence of the floor, not evidence of a bug — and not evidence of good
calibration either, which is why it should not be read as reassurance. In the real
analysis pool sizes are far larger, the floor drops, and BH fires. The script
prints the median floor, the BH threshold, a `floor_blocks_BH` flag and the number
of floor genes per cohort, so an empty survivor list can be told apart from a
genuinely null cohort. In a synthetic check with 20 planted DE genes, all 20 were
recovered and **all 20 sat at the floor** — detected, but unrankable among
themselves.

### The pair count overstates the real resolution by 40–160×

`p_floor` is a *hard* limit. There is a second, softer one that binds much earlier:
the 22,500 pairs are not 22,500 independent draws. They come from 150 + 150
residuals, and at n_rep = 3 each gene contributes only 2 degrees of freedom.
Resampling whole neighbourhoods 500 times and reading off the effective sample size
as `N_eff = p(1−p)/var(p̂)` at a true tail probability of 1e-2:

| pool heterogeneity (sd of log per-gene sd) | N_eff | % of 22,500 pairs | honest resolution 1/N_eff | overstated by |
|---|---|---|---|---|
| 0.00 (single scale) | 595 | 2.6 % | 1.7e-3 | 38× |
| 0.25 | 345 | 1.5 % | 2.9e-3 | 65× |
| 0.50 | 191 | 0.85 % | 5.2e-3 | 118× |
| 0.75 | 145 | 0.64 % | 6.9e-3 | 156× |

and the scaling, at a fixed true tail p = 1e-2 (het = 0.5):

| k | n_rep | pool | pairs | N_eff | N_eff / pairs |
|---|---|---|---|---|---|
| 50 | 3 | 150 | 22,500 | 216 | 1.0 % |
| 100 | 3 | 300 | 90,000 | 353 | 0.4 % |
| 30 | 5 | 150 | 22,500 | 681 | 3.0 % |
| 15 | 10 | 150 | 22,500 | 3,514 | 15.6 % |

Two consequences:

1. **Replicates buy resolution; neighbours barely do.** At a fixed pool of 150
   residuals, going from n_rep = 3 to 10 raises `N_eff` from 216 to 3,514 — 16×.
   Doubling k at n_rep = 3 raises it 1.6×. The floor improves as `k²` while the
   real precision improves roughly as `k`, so **growing k widens the gap between
   nominal and honest resolution rather than closing it.** Reaching the BH-relevant
   region (~4e-6) on precision, not just on the floor, is out of reach at n_rep = 3
   for any k a transcriptome can supply.
2. **The two problems share a cause.** The same neighbourhood heterogeneity that
   makes the raw null leptokurtic (§2) also cuts `N_eff` by 4× (595 → 145). Fixing
   the mixture would improve both the calibration and the resolution.

Practical reading: below roughly `1/N_eff`, a p-value is reporting neighbourhood
sampling noise, not evidence about the gene. At n_rep = 3, k = 50 that boundary is
around 5e-3, three orders of magnitude above the printed floor. This does not make
the p-values invalid — the model is still, on average, scoring the right quantity —
but it means the extreme tail is driven by which genes happened to draw a quiet
neighbourhood, which is exactly the region BH reads. `p_floor` is already reported
per gene; `N_eff` is not yet estimated per gene, and should be.

---

## 7b. A derived tail instead of a counted one

### How edgeR / limma / DESeq2 avoid a floor

All three share one structure:

1. **Borrow noise information across genes.** limma: empirical-Bayes shrinkage of
   gene-wise residual variances toward a prior fitted across all genes — with
   `trend = TRUE` / voom, an *expression-dependent* prior from a LOESS
   mean–variance trend. edgeR: NB dispersion shrunk toward a mean-dependent trend,
   then a quasi-likelihood dispersion shrunk the same way. DESeq2: dispersion
   shrunk toward a mean-dependent trend.
2. **Standardise the effect** by that borrowed noise estimate.
3. **Read the tail off a theoretical distribution** whose df is inflated by the
   borrowing — moderated *t* on `d_g + d_0` df (limma), `F(1, d_g + d_0)` (edgeR
   QL), Wald normal (DESeq2).

Step 3 is analytic, so `p = 1e-30` is arithmetic rather than resolution. Nothing is
counted, so nothing has a floor.

**TOX already does step 1.** The kNN mean-neighbourhood is the same device as
limma's expression-dependent prior and edgeR's mean-dependent trend. What TOX does
instead of steps 2–3 is count pairs.

### Keeping the tail analytic without assuming a distribution

A theoretical tail need not be an *assumed* one. The saddlepoint (Lugannani–Rice)
approximation reads the tail off the pool's **own empirical cumulant generating
function** — the reference distribution stays the data's, but it is evaluated
analytically. Written for the **mean-difference** statistic, it also supplies the
CLT shape the observed statistic has, i.e. the shape fix from §2. One change,
both problems. It is deterministic (no RNG), and costs one 1-D root find per gene.

Implemented as `tox_pvalue_saddlepoint()`; evaluated by
`TCGA_test/scripts/pvalue_tail_methods.R`.

**A1 — validated against a null with a closed-form answer** (large N(0,1) pools, so
`D ~ N(0, 2/n)`):

| true p | 1e-2 | 1e-4 | 1e-6 | 1e-8 | 1e-10 | 1e-12 |
|---|---|---|---|---|---|---|
| saddlepoint / truth | 1.15 | 1.30 | 1.41 | 1.41 | 1.21 | 0.81 |

Bounded relative error twelve orders down. (An Edgeworth or normal expansion drifts
without limit instead.)

**A2 — the four methods on a heterogeneous n_rep = 3 pool** (k = 50, 22,500 pairs,
floor 4.44e-5; reference = 2e7 bootstrap draws from the same pools):

| true p | exact pairs | boot 1e4 | saddlepoint | GPD(1e4) |
|---|---|---|---|---|
| 1e-3 | 3.9e-3 | 9.0e-4 | 2.2e-3 | 1.0e-3 |
| 1e-4 | 4.4e-4 | 3.0e-4 | 2.3e-4 | 1.1e-4 |
| 1e-5 | **1.3e-4 (floored)** | **2.0e-4 (floored)** | 2.3e-5 | 9.4e-6 |

Both counted methods bottom out; the saddlepoint is still tracking, within ~2×.

### B — what it does *not* fix, and this is the deciding result

The saddlepoint evaluates the tail of a **given pool** to arbitrary precision. It
cannot reduce uncertainty in the pool. Across 200 independently redrawn
neighbourhoods (n_rep = 3, k = 50, het = 0.5), at the same true threshold:

| true p | median p̂ | IQR | 5–95 % span | within 2× of truth |
|---|---|---|---|---|
| 1e-2 | 7.7e-3 | 3.8e-3 – 1.5e-2 | 33× | 48 % |
| 1e-3 | 1.5e-4 | 3.9e-5 – 8.2e-4 | 1,451× | 22 % |
| 1e-4 | 3.0e-6 | 9.7e-8 – 3.0e-5 | 377,867× | 9 % |

So at low replicate counts an analytic tail **replaces a visible, honest floor with
a number that looks precise and is not.** The floor at least marks where the
information stops.

### Recommendation

- **Adopt the analytic tail where the pool is well estimated** — the real
  case-vs-control analysis, tens to hundreds of samples. There it removes the floor,
  fixes the shape mismatch, and removes the bootstrap's RNG dependence at once.
- **At n_rep = 3, report a bound** (`p < 1/N_eff`) rather than a small number,
  whichever tail method is used. No tail method turns n = 3 into n = 30.
- **Report `N_eff` per gene** alongside the p-value, so the two limits — floor and
  pool uncertainty — are both visible.
- GPD peaks-over-threshold (Knijnenburg et al. 2009, *Bioinformatics*) is the
  simpler, better-precedented alternative and performed comparably here; it needs a
  threshold choice and was less stable across the two draws tested.
- The most limma-like option — shrink `sd_own` toward the neighbourhood sd (that
  *is* empirical Bayes) and use a moderated-*t* tail — is simplest and fastest, and
  would not be as wrong as it sounds: the observed statistic is a difference of
  means and Gaussianises, which is precisely why the individual-residual null
  mismatches it. Untested here; worth a direct comparison.

---

## 7c. Choosing the null construction for publication

### Two doubts about the bootstrap, checked against the code

**"Does the bootstrap null already contain the signal?"** — No. The gene means are
never added back. `prepare_sorted_data_helper` stores
`bessel * (x_ig − mu_g)` (raw) or `bessel * (log2(x_ig + c) − ghat_g)` (log), both
centred *within each group*, and
`compute_pvalue_bootstrap_mean_helper` resamples from `pool_case` / `pool_control`
only — the draw loop touches nothing but those two arrays. The null is centred at
zero by construction and the observed mean-difference is scored against it. Had the
means been added, the null would be centred on `obs` itself and every p-value would
sit near 0.5 — a one-line diagnostic if you ever want to confirm it.

There *is* a smaller real version of the concern: the target gene is its own nearest
neighbour, so its own residuals form ~1/k of its own pool. That is within-group
noise, not signal, but it does make the null slightly self-dependent — a
leave-one-out option would settle it cheaply.

**"Residuals cancel when averaged, so the null understates the noise."** — The
cancellation is required, not a defect. The observed statistic is *also* a mean and
cancels identically; a mean of `n` values has sd `σ/√n`, and the null must too.
Averaging `|residuals|` instead would give a null far too wide and a hopelessly
conservative test — this is the same point the `√n`-scaling section of
`noise_model_current_state.md` makes.

But the intuition is aimed at something real, just one step over: an iid draw of
`n_rep` residuals from the *pooled* neighbourhood can combine one residual from a
quiet gene with one from a noisy gene. **No real gene's mean is formed that way** —
its `n_rep` values all come from one noise level. The flaw is incoherence in the
pooling, not cancellation in the averaging.

### Blocking the draw, and enumerating it

Fix the incoherence by **blocking**: pick a neighbour gene, then resample within
that gene. Each null draw then carries a single coherent noise level, like the
observed statistic.

And once blocked, the null can be written down exactly. Resampling `n` of a gene's
`n` residuals with replacement has only `C(2n−1, n)` distinct outcomes — 10 at
n_rep = 3, 126 at 5, 2.8e6 at 10 — so the whole blocked null is
`k · C(2n−1, n)` weighted values per side, scored with the same sorted-pool +
binary-search tail count the exact model already implements. No draws, no RNG, no
reproducibility caveat. Verified: multiset enumeration reproduces full `n^n`
enumeration to 2.5e-16, and reproduces the blocked bootstrap as B → ∞
(0.1002 vs 0.1001, 0.01130 vs 0.01120, 9.25e-4 vs 9.18e-4 at B = 2e6).

### The resolution arithmetic (k = 30, n_rep = 3, G = 12,000, q = 0.05)

BH needs the best p-value to reach `q/G` = **4.2e-6**.

| construction | support | floor |
|---|---|---|
| exact pairs | 90 × 90 = 8,100 | 1.2e-4 |
| bootstrap B = 25,000 | 25,001 | 4.0e-5 |
| bootstrap matching 4.2e-6 | — | **needs B ≈ 240,000** |
| **blocked enumerated** | 810 × 810 = 656,100 | **1.5e-6** |

The combinatorial ceiling is not what binds a bootstrap: distinct multisets of 3
from 90 residuals number 125,580 per side, ~1.6e10 pairs of them, so 25,000 draws
is nowhere near exhausting the space — it is simply 25,000 draws, and `1/(B+1)` is
the floor regardless of how much room is left above it.

### Calibration under a complete null

Split-half H0, raw normalisation, k_start = 30 / k_max = 50, B = 5,000,
6,000 gene-observations per cell, ± 2 Monte-Carlo SE:

| n_rep | construction | infl@.05 | infl@.01 | median p |
|---|---|---|---|---|
| 3 | exact pairs | 1.15 ± 0.12 | 0.57 ± 0.19 | 0.561 |
| 3 | bootstrap iid | 1.26 ± 0.13 | **1.42 ± 0.31** | 0.605 |
| 3 | saddlepoint (iid, B→∞) | 1.15 ± 0.12 | **1.88 ± 0.35** | 0.631 |
| 3 | bootstrap blocked | 0.83 ± 0.10 | 0.10 ± 0.08 | 0.481 |
| 10 | exact pairs | 1.40 ± 0.13 | 0.95 ± 0.25 | 0.550 |
| 10 | bootstrap iid | 1.24 ± 0.12 | **1.97 ± 0.36** | 0.636 |
| 10 | saddlepoint (iid, B→∞) | 1.00 ± 0.11 | **2.50 ± 0.40** | 0.694 |
| 10 | **bootstrap blocked** | **1.16 ± 0.12** | **1.07 ± 0.27** | 0.545 |

(cv_het = 0.7 throughout.) Three things fall out:

1. **The iid bootstrap is anti-conservative in the tail** — infl@.01 of 1.4 at
   n_rep = 3 and 2.0 at n_rep = 10, worse than the exact model's 0.57 / 0.95. The
   tail is exactly where BH reads. That is a serious mark against switching to the
   current bootstrap as the published default.
2. **The saddlepoint confirms this is structural, not Monte-Carlo noise** — as the
   B → ∞ limit of the same construction it gives the best α = 0.05 calibration in
   the table (1.00) and the worst tail (2.50).
3. **Blocking fixes the tail.** At n_rep = 10 the blocked null is the only
   construction calibrated at *both* levels. At n_rep = 3 it is over-conservative
   (0.83 / 0.10) — the safe direction, at a cost in power.

Enumerated blocked vs the current exact model at n_rep = 3, 5,000 observations per
cell:

| cv_het | exact pairs | blocked enumerated |
|---|---|---|
| 0.0 | 0.90 / 0.54 | 0.62 / 0.08 |
| 0.4 | 0.88 / 0.64 | 0.70 / 0.26 |
| 0.7 | 1.16 / 0.72 | 0.91 / 0.28 |

### BH false positives under a complete null — the metric that decides it

G = 2,500, n_rep = 3, k_start = 30, 16 runs, no gene DE anywhere, so every
rejection is false. Target 0.

| cv_het | construction | mean BH hits/run | max | runs with ≥1 | median min p |
|---|---|---|---|---|---|
| 0.4 | exact pairs | **0.81** | 5 | 3/16 | 4.4e-5 |
| 0.4 | blocked enumerated | 0.06 | 1 | 1/16 | 4.1e-3 |
| 0.7 | exact pairs | **0.75** | 5 | 3/16 | 4.4e-5 |
| 0.7 | blocked enumerated | 0.06 | 1 | 1/16 | 3.8e-3 |

This **corrects §7's reading of the floor.** A floor is not automatically safe. Genes
*pile up* at it — the median minimum p-value for the exact model is exactly the
floor — and BH rejects a tied block of `r` genes as soon as `r ≥ p_floor·G/q`. At
G = 2,500 that is `r ≥ 3`, which happens. So the coarse null *manufactures* BH
rejections at n_rep = 3 rather than preventing them. `hits_FDR05 = 0` in the TCGA
calibration runs reflects the larger G and larger pools there (`r ≥ 11`), not a
general property. The blocked-enumerated null has no pile-up and produces 13×
fewer false discoveries.

### What LFCseq and NOISeqBIO do

**LFCseq** is close to TOX by design. For each gene it takes a neighbourhood of
genes with similar expression strength (default **50 genes**), and builds the null
by splitting the samples of *each condition* into two subsets, computing log fold
changes between those within-condition subsets for every neighbourhood gene, and
pooling them; the probability of non-differential expression is the fraction of that
pool exceeding the observed |LFC|. Two things to take from it:

- It is independent validation of the mean-neighbourhood device, at a comparable k.
- **Its null lives at the difference-of-group-means scale by construction.** There
  is no individual-residual null and therefore no `√n` correction and no shape
  mismatch — the defect diagnosed in §1.2/§1.3 simply does not arise. It reaches the
  same place the blocked/mean-level null does, by a different route.
- It is still an empirical fraction, so it has the same resolution floor.

**NOISeqBIO** (verified from the package documentation; the primary paper was not
reachable — PMC is captcha-gated and the Bioconductor vignette robots-disallowed, so
the θ definition and the exact probability formula below are *not* verified):

- The noise distribution is generated by **`r` permutations** (default 50) by
  resampling.
- M and D are **divided by `S + a0`**, where S is the statistic's standard error and
  `a0` is the `a0per` percentile (default 0.9) of S across all features — a
  shrinkage/fudge-factor exactly analogous to limma's `s0` and SAM's, and something
  TOX has no equivalent of. It stabilises genes whose own S is near zero.
- The noise distribution is **smoothed by kernel density estimation** (`adj`) before
  probabilities are read off it.
- The DE probability comes from comparing a mixture density `f` against the noise
  density `f0` — an Efron-style local-fdr construction rather than a p-value, so
  NOISeqBIO sidesteps BH entirely.

The transferable point: **NOISeqBIO does not count raw draws either.** It smooths
the null first. That is a third route past resolution granularity, alongside
enumeration (§7c) and the saddlepoint (§7b) — though KDE tails are the least
trustworthy of the three far out. Its `S + a0` shrinkage is worth considering
independently of the tail question: it is cheap, well-precedented, and targets the
same per-gene variance instability that `rho` measures in §5.

### Recommendation for the published tool

1. **Do not adopt the iid bootstrap as the default.** It is anti-conservative at
   α = 0.01 (1.4–2.0×) — the region BH reads — and B = 25,000 does not reach
   BH-relevant resolution anyway (4.0e-5 vs the 4.2e-6 needed).
2. **Adopt the gene-blocked mean-level null, enumerated.** It combines what each
   existing model gets right: the bootstrap's mean-level (CLT-shaped) null, the
   exact model's determinism and enumeration machinery, and gene coherence on top.
   At n_rep = 3, k = 30 it gives a floor of 1.5e-6 with zero draws, and 13× fewer
   BH false positives than the current exact model.
3. **Enumerate up to n_rep ≈ 10** (`C(2n−1, n)·k` = 2.8e6 values, ~22 MB per side at
   k = 30) and fall back to the sampled blocked null above that.
4. **Consider a `S + a0`-style shrinkage** (NOISeqBIO, limma, SAM all use one).
5. **Keep reporting the floor and an `N_eff` estimate per gene.** Enumeration removes
   an artificial floor that was binding *above* the information limit; it does not
   create information. At n_rep = 3 the honest resolution is still ~1e-2.

All of §7c is simulation under one cohort model. `null_construction_comparison.R`
reproduces every table; the same comparison should be run on TCGA split-halves
before any of it is treated as settled.

---

## 8. Running

```bash
cd TCGA_test/scripts
Rscript reanalyse_null_calibration_sweep.R            # seconds, no data needed
Rscript raw_anticonservative_diagnosis.R              # 3 cancers x 2 stages x 2 arms
Rscript tox_significant_genes_report.R                # all 8 cancers, TOX-log
Rscript null_calibration.R                            # now incl. boot + trim arms
```

`raw_anticonservative_diagnosis.R` is deliberately scoped to 3 cancers × 2 stages ×
10 splits: `null_calibration.R` already establishes *that* raw is anti-conservative
across all 8 cancers, so the job here is per-gene detail, not more cells.

**Caveat carried through all of the above:** every number in §1 and §4 comes either
from a re-analysis of an existing run or from simulation. Nothing in §5 has been
run against TCGA yet — the per-gene confirmation is the next step, not a result.
