# Multidimensional TOX noise model — implementation plan

Issue #185. Status 17.09.2026, rev. 4a. Written in English to match the other
planning notes and because it is intended for the issue thread.

This supersedes `md_noise_model_plan_independent_axes.md` in full. It is
rewritten on the design agreed with Asis on 16.09.2026: **residual vectors with
per-axis √n scaling**, rather than the bootstrap-of-means construction in earlier
revisions. `md_noise_model_implementation_sketch.md` remains the evidence base
for the statistical claims referenced here by section number.

---

## 1. Decisions taken

| decision | source |
|---|---|
| Null is built from **residual vectors**: one residual drawn per axis per side, each scaled by `1/√n` of its own side. Not a bootstrap of means | Asis, 16.09. |
| Axes are treated as **independent**; the Cartesian product of per-axis residuals is the null | Asis, 16.09. |
| A gene is **one vector**; axes are not analysed individually within this test | Aaron, 16.09. |
| Report **both norms**: unweighted as effect size, per-axis standardised as test statistic | agreed 17.09. |
| **Normalisation is tested empirically**, counts vs TPM, rather than argued from the literature | Aaron, 16.09. |
| Enumerate when the outcome space allows it, sample otherwise, **report which was used** | Aaron, 16.09. |
| Systematic ("pseudo-Monte-Carlo") sampling instead of RNG where possible | Aaron, 16.09. |
| `< 10 residuals` gate is **not** a live concern (see §5.4) — closed | Aaron, 17.09. |
| Multidimensional expression filter **deferred**, not implemented in this draft | Aaron, 17.09. |
| Trajectories, velocity, acceleration, integrals: **out of scope** | Aaron, 17.09. |

**Scope assumption, enforced in code:** every sample belongs to exactly one axis.
No plant, subject, extraction batch or control group is shared between two axes.
This is what makes the Cartesian-product null correct. §9 gives the check and the
guard.

---

## 2. Statistic

Per gene `g`, per axis `a = 1…d`:

```
beta_hat(a)  =  mean_case(g,a)  -  mean_ctrl(g,a)          (#185 §1)
```

With two groups and no covariates this is the OLS estimate of `β_a` in
`x_as = α_a + β_a z_s + ε_as`, i.e. per axis identical to today's `obs_own`. So
`d = 1` reduces exactly to the current pipeline.

Two norms, both computed, both calibrated (the weights are applied to observed and
null alike):

```
D        = sqrt( sum_a beta_hat(a)^2 )                      effect size, log2FC units
D_std    = sqrt( sum_a beta_hat(a)^2 / Var_null(a) )        test statistic, unitless
```

They answer different questions — how large the response is, versus how confident
we are that there is one — exactly as `logFC` and `t` do in limma. `D` is reported
next to each gene; `D_std` drives the p-value and the ranking. With independent
axes the null covariance is diagonal, so `D_std` **is** the full Mahalanobis
statistic; no covariance matrix, no Cholesky, no ridge.

`Var_null(a)` is available in closed form from the pools, with no sampling:

```
Var_null(a) = Var(pool_case(a)) / n_case(a)  +  Var(pool_ctrl(a)) / n_ctrl(a)
```

`D` is upward-biased as a magnitude (`E[D] > ‖β‖`, and `E[D] ≈ σ_β√d ≠ 0` under
H0). Accumulate `mean(D²_null)` and also report
`D²_adj = max(0, D² − mean(D²_null))`. Rank on the p-value regardless.

---

## 3. Null

### 3.1 Construction

Per gene, per axis, gather the residual pools with the existing adaptive kNN in
mean space — one pool per side per axis, `2d` gathers per gene, using
`gather_residuals_helper` unchanged.

One null draw:

```
for a = 1..d:
    eps_plus   <- one residual from pool_case(a)
    eps_minus  <- one residual from pool_ctrl(a)
    beta_null(a) = eps_plus / sqrt(n_case(a))  -  eps_minus / sqrt(n_ctrl(a))
D_null      = sqrt( sum_a beta_null(a)^2 )
D_std_null  = sqrt( sum_a beta_null(a)^2 / Var_null(a) )
```

**Note the scaling is per side, not per axis.** Each side is divided by the root
of *its own* replicate count before the difference is taken. That is what makes
ragged designs correct — with `n_case = [4,6,7,6,5]` and
`n_ctrl = [3,5,5,9,3]` each of the ten numbers enters its own term. A single
per-axis factor would be right only when both sides have equal `n`.

Why the scaling is needed at all: the observed statistic is a difference of means
over `n` replicates, while a raw residual difference is at single-observation
scale. Without the `1/√n` the null would be inflated by roughly `√n` and the test
strongly conservative. This is the same correction the 1-D implementation already
carries and the original `Empirical_Noise-Based_Significance_of_Distances.tex`
does not.

### 3.2 On the Cartesian product

Asis's observation that we do not know which sample in axis `a₁` corresponds to
which in axis `a₂` is exactly the condition under which the product null is
correct. Because the axes share no samples, drawing each axis independently
reproduces the true joint null.

It also answers the `200⁵` worry directly: **that number is never materialised.**
It is the size of an implicit measure, not a set to be built. And it does not
represent `200⁵` units of information — the information is bounded by the `200`
residuals per axis. The product size decides only whether we enumerate or sample
(§4); it is not itself a quantity to be capped.

### 3.3 p-value

```
p = (1 + #{D_std_null >= D_std_obs}) / (M + 1)
```

for the sampled paths, and the exact weighted tail count for the enumerated one.
The formula and its interpretation are unchanged from 1-D — the null values are
now distances of `d`-vectors instead of scalars, and nothing about the tail count
cares which. Aaron's expectation here is correct.

---

## 4. Three paths to the p-value

The outcome space per axis has

```
S_a  =  (k_case(a) · n_case(a))  ×  (k_ctrl(a) · n_ctrl(a))     atoms
```

— every pairing of one case-side residual with one control-side residual. The
joint space is `∏_a S_a`. Which path is used is decided per gene and **reported
per gene** in the output.

### 4.1 Enumeration — exact, rarely applicable

Walk the joint space, weighted tail count. Exact, no floor.

The arithmetic is discouraging and should be stated plainly rather than assumed
away. At the defaults (`k = 15` genes, `n_rep = 3`) each pool holds 45 residuals,
so `S_a = 45 × 45 = 2025`:

| d | joint space |
|---|---|
| 1 | 2 025 |
| 2 | 4.1 × 10⁶ |
| 3 | 8.3 × 10⁹ |
| 5 | 3.4 × 10¹⁶ |

**Correction to what the existing constants are.** There is no 10,000-*pair*
routing threshold; earlier revisions of this note said there was. Reading the
source:

- `N_BOOTSTRAP_DRAWS = 10000` is the Monte-Carlo **draw count** in the pooled path
  of `tox_noise_model.F90`, not a routing threshold.
- `BLOCKED_MAX_ENUM_VALUES = 8000` caps the *blocked* multiset enumeration
  (`C(2n−1,n)·k` per side).
- `tox_noise_model_exact.F90` has **no threshold at all** — its own header states
  it needs neither constant.

The reason is in `compute_pvalue_helper`: it sorts the control pool once, then for
each case residual runs two binary searches for `#{b < a+t}` and `#{b ≤ a−t}`. No
pair is ever materialised, so the cost is `O(k log k)` and exactness in 1-D is
effectively free at any pool size. There is no "gap" to widen.

That also sharpens why `d ≥ 2` breaks: it is not a cost limit anyone chose. The
binary-search trick does not generalise, so the product would have to be built,
and no threshold value reaches 10¹⁰.

Note also that **larger neighbourhoods make enumeration worse, not better** —
`S_a` grows quadratically in `k`:

| k (genes) | n_rep | pool | S_a | d = 2 | d = 3 |
|---|---|---|---|---|---|
| 15 | 3 | 45 | 2 025 | 4.1 × 10⁶ | 8.3 × 10⁹ |
| 20 | 3 | 60 | 3 600 | 1.3 × 10⁷ | 4.7 × 10¹⁰ |
| 20 | 5 | 100 | 10 000 | 1.0 × 10⁸ | 1.0 × 10¹² |
| 50 | 5 | 250 | 62 500 | 3.9 × 10⁹ | 2.4 × 10¹⁴ |

So with the `k = 20+` settings normally used, the enumerated branch is in practice
the `d = 1` regression path and nothing more. Keep it for that purpose; do not
plan around it for multidimensional use. Aaron's `3⁵` example corresponds to pools
of a single gene.

So: keep the branch, because it is where `d = 1` lives and it is the regression
test against the current module, but do not expect it to run in multidimensional
use. The threshold must be evaluated on the **product** `∏_a S_a`, not per axis —
each axis looks harmless in isolation.

### 4.2 Sampling — the working default

`M` draws, `p = (1 + B)/(M + 1)`, floor `1/(M+1)`.

**Systematic sampling instead of RNG** (Aaron's pseudo-Monte-Carlo). Rather than
drawing random indices, walk the product space with a fixed stride per axis,
chosen coprime to that axis's pool size so the marginals stay balanced. Benefits:
deterministic and reproducible, no RNG overhead, no seed dependence, parallelises
without stream-splitting, and lower variance than random draws at the same `M`.

One caveat for the method text: `(1+B)/(M+1)` is justified as a *valid* Monte-Carlo
p-value under random draws. Under systematic sampling it is instead an
**approximation of the exact enumerated tail** — usually better in practice, but
the justification changes, and the error becomes a coverage error rather than a
sampling error. Both variants should be implemented and compared on the mock null
before systematic becomes the default.

**Resolution.** `M` is set by BH over `G` genes, not by `d`. A gene at the floor is
rejected only if `r ≥ p_floor·G/q`, so at least that many genes must be tied at
the floor before any of them is rejected, with no ranking inside the block. The
dimension does not enter: we estimate a one-dimensional tail probability,
`SE(p̂) = √(p(1−p)/M)`, independent of `d`. This answers Asis's under-sampling
question — there is no dimensional surcharge on `M`.

`M` is therefore a **user parameter** (`n_draws_max`), not a compile-time
constant, defaulting to **20 000**:

| M | p_floor | genes tied at floor before BH rejects (G = 20 000, q = 0.05) |
|---|---|---|
| 10⁴ | 1e-4 | 40 |
| **2×10⁴ (default)** | 5e-5 | 20 |
| 4×10⁵ | 2.5e-6 | 1 |

`M = 4×10⁵` is the point at which the floor stops binding altogether (the
top-ranked gene becomes rejectable). As a flat cost that is 4.8×10¹⁰ draws at
`d = 3`; with sequential stopping it is close to free (§8).

### 4.3 Convolution — exact tail without enumeration, no floor

`D² = Σ_a β_a²` is a **sum of independent terms**, and sums convolve. Per axis:
square the `S_a` atoms, bin onto a shared grid, FFT-convolve the `d`
distributions, `O(d N log N)`.

This fits Asis's design better than it fitted any earlier revision, because the
per-axis atom set *is* the existing pairwise enumeration — 2 025 values at the
defaults, trivially small. What was infeasible was only the product.

Measured against Monte Carlo (`md_toy6_convolution.py`, `d = 3`, 2 500 atoms per
axis, heterogeneous scales, `M = 2×10⁷`):

| target p | MC | convolution | ratio |
|---|---|---|---|
| 1e-3 | 1.000e-03 | 9.991e-04 | 0.999 |
| 1e-4 | 1.000e-04 | 1.020e-04 | 1.020 |
| 1e-5 | 1.000e-05 | 1.028e-05 | 1.028 |
| 1e-6 | 1.000e-06 | 1.260e-06 | 1.260 |

Within 3 % through `p = 10⁻⁵`. At `10⁻⁶` the MC reference rests on ~20
exceedances (relative SE ~22 %), so there the *reference* is the unreliable one.
The error changes character: discretisation error instead of sampling error, and
discretisation error has no floor.

Not in the first draft. Staged as Phase 6, conditional on the floor actually
binding (validation item 13). Open: per-gene cost is unmeasured, and the grid
range must be set per gene from the atom supports or the tail of interest is
truncated.

---

## 5. Module

### 5.1 File and reuse

New file `src/tox/tox_noise_model_md.F90`, module `noise_model_md`, following the
existing helper/main/alloc layering and FORD conventions. Separate from the two
existing modules because the ABI carries an extra dimension; `d = 1` then gives a
hard regression test.

Reuse by `use`, do not copy: `sorted_data_t`, `prepare_sorted_data`,
`gather_residuals_helper`, `trim_pool_tails_helper`, `find_closest_helper`,
`choose_index`. Either make the needed ones public (one line, no behavioural
change) or lift them into a shared `tox_noise_pools` module — the latter removes
an existing verbatim duplication between `tox_noise_model.F90` and
`tox_noise_model_exact.F90` but touches validated code. Asis's call; Q3.

### 5.2 Data layout

Ragged replicate counts per axis and per side are required. No rectangular
`(n_rep, n_genes, d)` array and no padding — padding wastes memory but, worse,
puts sentinel values into a matrix whose consumers have no concept of sentinels
(`prepare_sorted_data_helper` would average over them silently).

Flat contiguous buffer plus an offset table, per side:

```
replicates_packed(sum_a n_rep_a * n_genes)    ! axis blocks, each (n_rep_a, n_genes)
n_rep_per_axis(n_axes)
axis_offset(n_axes + 1)                       ! derived prefix sums
means(n_genes, n_axes)
```

**Samples-fastest within each block**, matching
`prepare_sorted_data_helper(replicates(n_samples, n_genes))`, which reads whole
contiguous columns `replicates(:, orig_idx)`. Genes-fastest would stride a gene's
replicates by `n_genes`.

Each axis is handed to the existing routine as an ordinary 2-D matrix by pointer
bounds remapping, with no copy:

```fortran
real(real64), dimension(:), intent(in), target :: replicates_packed
real(real64), pointer, contiguous              :: axis_reps(:,:)

do a = 1, n_axes
    axis_reps(1:n_rep_per_axis(a), 1:n_genes) => &
        replicates_packed(axis_offset(a)+1 : axis_offset(a+1))
    call prepare_sorted_data_helper(means(:, a), axis_reps, &
                                    n_rep_per_axis(a), n_genes, &
                                    norm_method, sorted(a), ...)
end do
```

Do **not** pass the rank-1 section straight to the rank-2 dummy: sequence
association would allow it for an external procedure, but these are module
procedures with explicit interfaces and rank mismatch is then a hard error.

Memory at `d = 3`, `n_rep = 5`, `G = 20 000`: ≈ 4.8 MB packed residuals. Not a
constraint.

### 5.3 Entry points

```fortran
!> beta_hat per axis, per #185 §1.
subroutine compute_beta_md(means_case, means_control, n_genes, n_axes, &
                           norm_method, beta_obs, ierr)

!> The test. Public; accepts beta_obs from compute_beta_md or an external estimator.
subroutine compute_noise_pvalue_pipeline_md( &
    means_case, replicates_case_packed, n_rep_case_per_axis, &
    means_control, replicates_control_packed, n_rep_control_per_axis, &
    beta_obs, compute_pvalue_own, beta_mode, beta_centre, &
    n_genes, n_axes, norm_method, k_start, k_step, k_max, tau, trim_frac, &
    null_method, sampling_mode, n_draws_max, n_exceed_target, enum_max_product, &
    seed, max_pool_size, &
    pvalues_own, d_obs, d_std_obs, d_sq_null_mean, &
    method_used, n_draws_used, &
    neighborhood_size_case, neighborhood_size_control, var_null, &
    beta_check_cor, beta_check_mad, delta_hat, &
    n_genes_with_pvalue, ierr)
```

| argument | shape | meaning |
|---|---|---|
| `beta_obs` | `(n_axes, n_genes)` | gene-contiguous; ignored when `beta_mode = 0` |
| `beta_mode` | scalar | `0` compute internally (default, scale-coherent); `1` caller-supplied |
| `beta_centre` | scalar | `0` none; `1` subtract per-axis median `β̂_a` across genes (§6) |
| `sampling_mode` | scalar | `0` systematic stride; `1` RNG |
| `n_draws_max` | scalar | `M`, user-settable, **default 2×10⁴**; p-value floor is `1/(M+1)` |
| `n_exceed_target` | scalar | Besag–Clifford stopping; `0` disables (default) |
| `enum_max_product` | scalar | enumerate while `∏_a S_a ≤` this; else sample |
| `d_obs`, `d_std_obs` | `(n_genes)` | the two norms of §2 |
| `d_sq_null_mean` | `(n_genes)` | for the bias correction |
| `method_used` | `(n_genes)` | `0` enumerated, `1` systematic, `2` RNG, `3` convolution |
| `n_draws_used` | `(n_genes)` | achieved resolution |
| `var_null` | `(n_axes, n_genes)` | per-axis null variance; the weights, exposed |
| `neighborhood_size_*` | `(n_axes, n_genes)` | per-axis pool sizes |
| `beta_check_cor/mad` | `(n_axes)` | supplied vs recomputed `β̂`; `beta_mode = 1` only |
| `delta_hat` | `(n_axes)` | estimated composition shift, reported even when centring is off |

`norm_method`, `k_start`, `k_step`, `k_max`, `tau`, `trim_frac`, `null_method`,
`max_pool_size` keep their current meaning, applied per axis.

There is **no `null_coupling` argument** — the scope assumption makes it constant;
§9 gives the guard instead.

### 5.4 The `< 10 residuals` gate

Confirmed not a live concern: `k_start`/`k_step`/`k_max` count **genes**, not
residuals, and each gene contributes `n_replicates`. At 15 genes × 3 replicates
the pool holds 45, and the gate still clears after trimming. With no
stratification nothing further is removed. The branch stays for safety; the
question of what to do when it fires is closed.

### 5.5 Caller-side, in R

Direction `β̂/D`, per-axis contribution shares `β̂_a²/D²`, and BH — over **tested
genes only**, as previously fixed. See §7 on when direction may be reported.

---

## 6. Normalisation

Tested empirically rather than argued. Two arms on the same dataset, everything
else held fixed:

| arm | input | composition correction | where |
|---|---|---|---|
| counts | raw counts | TMM (`calcNormFactors`) or median-of-ratios (`estimateSizeFactorsForMatrix`) | R, existing library functions |
| TPM | TPM | `beta_centre` — per-axis median of `β̂_a` across genes | Fortran |

Naming, so the comparison is labelled correctly: what runs on the Fortran side is
**median centring of log-ratios**, not TMM. TMM's weights come from count
magnitudes, which do not exist on TPM. Same estimator family, more robust variant
on counts.

**`beta_centre` must be off whenever the R side has already normalised**, or the
data are centred twice.

Why correct at all: on TPM, `log2FC = β − Δ` with `Δ` identical for every gene,
and the null — built from *within-group* residuals — cannot contain information
about a *between-group* offset. Running uncorrected is not the assumption-free
option; it is the assumption `Δ = 0`. Median centring assumes only that the median
gene is unchanged, which is strictly weaker. Report `Δ̂_a/σ̂_a` per axis either way.

Related and worth measuring in the same run: **TPM removes the count magnitude,
which is what predicts technical precision**, and its length correction actively
works against the kNN matching (at equal TPM, a longer gene carries more counts
and is more precise — roughly a 4–5× CV spread across a 500 bp–10 kb range among
genes the matcher treats as equivalent). This should show up in **power**, not
calibration, which is why the existing calibration results do not speak to it.

---

## 7. Filtering — deferred, open

Not implemented in this draft. Recorded so it is not lost.

The `< 10` gate closes, but the underlying question returns in a different form in
`d` dimensions. A gene strongly expressed in leaf and silent in root passes a
global filter but carries a **degenerate root axis**: its kNN neighbours there are
also near-zero genes with tiny residuals, so that axis's null is extremely narrow,
while a few stray reads still give `β̂_root ≠ 0`. A narrow null against a
non-narrow observation means that axis dominates `D` and drives the significance —
for a gene that is not expressed in root at all.

Aaron's proposal, to be worked out: something like **TPM > 1 in ≥ 70 % of samples
in at least one axis**, on the grounds that being expressed in one tissue and not
another is itself the information we want to keep, so requiring expression in all
`d` axes would discard exactly the tissue-specific genes of interest.

Cheap and decisive preliminary: count how many genes pass the filter in **all**
axes versus in **at least one**. For four TCGA cohorts of related tissue the gap
should be small; for Arabidopsis leaf/stem/root probably large. That number
decides whether this is a footnote or a design question.

---

## 8. Cost

Per gene: `2d` pool gathers (existing code), then `M · d` residual pairs. At the
default `M = 2×10⁴`, `d = 3`, `G = 20 000`: ≈ 2.4×10⁹ draws, against 2×10⁹ for the
current 1-D bootstrap path. So **three dimensions at double the draw count costs
roughly what one dimension costs today** — Asis's construction draws two residuals
per axis rather than resampling `n_rep` of them, which pays for the extra axes.

Further reductions, in order of how much they are worth:

- **Systematic sampling** removes the RNG entirely (§4.2).
- **Besag–Clifford sequential stopping**, which is what makes a large `M`
  affordable. Sample until `c` exceedances (`p = c/m`) or `M` draws; expected
  draws `≈ min(M, c/p)`, so genes with large `p` stop almost immediately whatever
  `M` is, and only the tail pays the cap. Averaging over `p̂ ~ U(0,1)` with
  `c = 20`:

  | M | expected draws/gene |
  |---|---|
  | 10⁴ | ≈ 144 |
  | 2×10⁴ | ≈ 158 |
  | 4×10⁵ | ≈ 218 |

  A 40× increase in `M` costs about 1.5× in average draws. That is the
  combination to aim for — **large `M` plus sequential stopping** — rather than a
  modest flat `M`. Caveats stand: the sequential p-value is a different estimator
  with granularity `c/m`, and its interaction with BH is unverified. Ship it off
  by default (`n_exceed_target = 0`) and enable once the mock null agrees with the
  fixed-`M` arm. If it validates, the floor question largely dissolves and Phase 6
  becomes optional rather than the plan's answer to resolution.
- **Per-gene deterministic streams** — with systematic sampling this comes free:
  each gene's walk is a pure function of its index, so the gene loop parallelises
  under OpenMP with no stream splitting and no order dependence.

---

## 9. The scope assumption, and how it is kept honest

Two checks on the intended dataset, before relying on this plan:

1. **Read the sample sheet.** Does a plant, subject or extraction/library batch ID
   recur across axes? Is one control group shared by several axes (the TCGA
   Stage I–IV shape)? If yes to either, the product null is not correct as written.
2. **Measure the across-axis residual correlation**, as a backstop for structure
   the metadata does not record. Report it per run alongside the other diagnostics.

Reference points from §4.1 of the sketch: `ρ ≈ 0.3` gives roughly 1.6× inflation
at `α = 0.01`; `ρ ≈ 0.1` is worth documenting and ignoring. Blocking on the
replicate index alone does **not** fix it when the correlation is gene-specific —
only a joint neighbourhood does, which is out of scope here.

**Guard in code.** The R wrapper compares sample identifiers across the per-axis
matrices and **errors** if any identifier appears in more than one axis, naming
this assumption in the message. A few lines, and it converts a silent validity
failure into a loud one.

---

## 10. Validation

**Fortran unit tests** (`test/mod_test_noise_model_md.f90`):

1. `d = 1` reproduces `noise_model`'s `pvalues_own` exactly under the enumerated
   path (where both are deterministic, this is an exact-match test, not a
   tolerance test).
2. `compute_beta_md` at `d = 1` reproduces today's `obs_own` exactly, on both the
   linear and the log2 branch.
3. Norm correctness for both `D` and `D_std` against a hand-computed fixture.
4. `Var_null` closed form agrees with the sample variance of a large null draw.
5. Ragged `n_case = [4,6,7,6,5]`, `n_ctrl = [3,5,5,9,3]`: pools have the right
   sizes and each side is scaled by its own `√n`. Deliberately asymmetric, since a
   single per-axis factor would pass a symmetric test.
6. `enum_max_product` routes correctly: a tiny fixture enumerates, a realistic one
   samples, and `method_used` reports it.
7. Systematic and RNG sampling agree within Monte-Carlo error on the same fixture.
8. `beta_centre = 1` subtracts the per-axis median; with a synthetic constant
   offset added to every gene on one axis, `delta_hat` recovers it.
9. Degenerate input: zero-variance axis, `D = 0`, non-finite `β̂` — each hits the
   intended branch and error code.
10. `beta_mode = 1` with a deliberately off-scale `β̂` gives a large
    `beta_check_mad` while the run completes.

**Statistical validation** (testing repo):

11. **Mock-null calibration in `d` dimensions** — split one homogeneous condition
    into two fake groups across *all* `d` axes; every gene null by construction.
    Inflation at α = 0.05 and 0.01 for `d = 1, 2, 3, 4`. Primary acceptance gate.
    `d = 1` must match the existing TPM results. Natural dataset: healthy samples
    of TCGA-KIRC / LUSC / LUAD / STAD split in half per cohort.
12. **Counts vs TPM**, both arms of §6, same dataset, both calibration and power.
    This is the measurement that settles the normalisation question.
13. **Resolution**: `n_draws_used`, `method_used` and `1/(M+1)` per gene; fraction
    of BH rejections sitting at the floor, per run. Run at `M = 2×10⁴` and at
    `M = 4×10⁵` to see whether the floor binds in practice, and compare the
    sequential-stopping arm against the fixed-`M` arm for agreement. Decides
    whether Phase 6 is needed at all.
14. **Composition**: `Δ̂_a/σ̂_a` per axis, `‖Δ̂‖` against `E[D_null]`, rank
    correlation between centred and uncentred orderings, and the mean cosine of
    reported directions against `−Δ̂`. Calibrate the `beta_centre` default on the
    mock-null arm, where "most genes unchanged" holds by construction — not on
    tumour-vs-normal, where it does not.
15. **Power**, with matched resolution across arms. This matters: a per-axis TOX
    comparator would otherwise run through the exact path with no floor while the
    omnibus is MC-floored, biasing the comparison toward per-axis testing
    precisely in the tail where BH reads. Run both through the same machinery at
    the same `M`, or restrict to the region above both floors and say so.
    Comparators: per-axis TOX combined by **Šidák** and by **Simes**; limma's
    moderated F across the same `d` contrasts; `glmQLFTest` with a multi-column
    contrast; DESeq2 LRT against a reduced design.

---

## 11. Phases

| phase | content | blocked by |
|---|---|---|
| 0 | the two checks of §9 on the intended dataset | — |
| 1 | `tox_noise_model_md.F90`: `compute_beta_md`, product null with √n scaling, both norms, enumeration + systematic sampling, unit tests 1–10 | 0 |
| 2 | C ABI, Rcpp, R wrapper incl. the §9 guard; `d = 1` end-to-end reproduction | 1 |
| 3 | mock-null calibration (11), counts vs TPM (12), resolution (13), composition (14) | 2 |
| 4 | power benchmark (15); decide the `beta_centre` default; decide the filter (§7) | 3 |
| 5 | sequential stopping, OpenMP | 3 |
| 6 | convolution null (§4.3), **if** item 13 shows the floor binds | 3 |

Phases 1–2 are a few days. Phase 3 is where the real answer is.

---

## 12. Open questions

**Q1.** The multidimensional expression filter (§7). Deferred, but needed before
any biological interpretation.

**Q2.** `beta_centre` default on or off? Item 14 answers it, but a default must be
chosen before that runs.

**Q3.** Is the `tox_noise_pools` refactor acceptable in the same PR?

**Q4.** For external `β̂` (`beta_mode = 1`): warn, refuse, or surface
`beta_check_*` only? Concrete failure modes: edgeR's coefficients are on the
**natural log** scale and shrunk by default (`prior.count = 0.125`, with
`unshrunk.coefficients` separate). The two-number diagnostic distinguishes these —
correlation ≈ 1 with large MAD is a units problem, degraded correlation is
shrinkage. Only limma-voom is principled (β linear in the data, `E` matrix the
matching residual scale).

---

## 13. Limitations to carry into the write-up

1. **Exactness by enumeration is gone at `d ≥ 2`.** Both 1-D devices are
   intrinsically one-dimensional: the sorted-pool binary search exploits that on a
   line "far from a" is two intervals, i.e. contiguous index ranges, and no
   ordering of ℝᵈ makes a sphere's exterior contiguous; the multiset enumeration
   exploits a small known outcome space, and the joint space is a product.
   A k-d tree does not help — it indexes an existing point set, and these points
   are never materialised. The compact intuition: the sum of `d` dice is a single
   number, and naive enumeration is still `6^d`. Exactness needs an *enumerable
   null*, not a scalar statistic. Partly recoverable by convolution (§4.3).
   *Precedent:* limma hit this with `roast` (`(b+1)/(nrot+1)`, floor
   `1/(nrot+1)`) and answered it with `fry`, an analytic approximation not limited
   by rotation count.
2. **Composition shift on TPM** (§6). Two harms, behaving differently. False-positive
   inflation is governed by `Δ/σ`, not `d` — the relative inflation of `E[D²]` is
   `1 + (Δ/σ)²` with no `d` in it; what grows with `d` is only the concentration of
   the `χ²_d` null, so FPR@.01 roughly doubles from `d = 1` to `d = 12`. **Direction
   bias is the genuinely new `d`-dimensional harm**: because `Δ` is common to all
   genes, every null gene's direction is pulled toward `−Δ` (mean cosine 0.44 at
   `Δ/σ = 0.5`, essentially independent of `d`), and it does **not** shrink as α is
   tightened.
3. **Direction is descriptive, not inferential.** With equal per-axis noise the
   null of `D` is rotation-invariant, so power depends only on `‖β‖` and not on
   which tissues respond (0.643 / 0.642 / 0.643 at `d = 6` for effects on 1, 3, 6
   axes). #185 §12 reads as though direction is tested; it is not. Combined with
   (2), **I would gate any direction or coordinated-response output on item 14**:
   a family "responding coherently across tissues" is exactly the signature a
   composition shift manufactures, and it is the most persuasive-looking wrong
   figure this pipeline can produce.
4. **The omnibus does not replace per-tissue DE.** At fixed `‖β‖` it loses power as
   `d` grows (0.81 at `d = 2` → 0.50 at `d = 12`); it wins for coordinated
   responses (`d = 6`, `r = 3`: 0.642 vs 0.518) and loses for tissue-specific ones
   (`r = 1`: 0.643 vs 0.717). Both should be reported. `diffSplice` faces the same
   choice over a gene's per-exon coefficients and the edgeR docs report the same
   crossover — independent confirmation, and the reason Simes belongs in item 15.
5. **Not assumption-free.** Distribution-free, yes. Composition is an assumption
   in every method including edgeR and DESeq2, and `Δ` is not identifiable from
   the data alone. Worth reporting the robustness margin `D − c` per gene — the
   smallest `‖Δ‖` that would render the gene non-significant — which asserts no
   value of `Δ` at all and is free to compute.
6. **All simulation figures quoted here are a Gaussian surrogate** — no
   mean-variance trend, no zeros, no kNN matching, no TPM structure. They size
   properties of the statistic and predict nothing about real data.
7. **Nothing has been compiled.** Claims about existing code come from reading
   `tox_noise_model.F90`, `tox_noise_model_exact.F90`, `tensoromics_functions.cpp/.R`
   and the testing repo on the branches cited in the sketch.