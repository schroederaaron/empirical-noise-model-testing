# The gene-blocked mean null — what it is, and how it differs from the current null

Reference for `null_method` in `tox_noise_model.F90`. Everything below describes code
as implemented; the measured numbers are from the simulations named at each point.

---

## 0. The difference in one paragraph

All three null constructions in TOX answer the same question — *how large a
case-vs-control distance can this gene's neighbourhood produce from noise alone?* —
and all three build the answer from the same residual pool. They differ only in **how
a single null "distance" is assembled from those residuals**:

| construction | one null distance is… |
|---|---|
| `exact` (tox_noise_model_exact) | \|one case residual − one control residual\|, each pool pre-scaled by 1/√n |
| `null_method = 0` "pooled" | \|mean of n residuals drawn iid from the case pool − same from the control pool\| |
| `null_method = 1` "blocked" | \|mean of n residuals drawn from **one** case gene − same from **one** control gene\| |

The blocked version changes one thing: the `n` residuals that get averaged into a null
mean now all come from **the same neighbour gene**, instead of being picked
independently from the pooled residuals of all of them. Everything else — the
neighbourhood, the residuals, the observed statistic, the p-value formula — is
untouched.

That one change turns out to matter for the *shape* of the null (§4) and, because the
blocked null has only finitely many outcomes, it lets the whole null be **enumerated
exactly** instead of sampled (§5), which removes the Monte-Carlo p-value floor (§7).

---

## 1. What is identical in all three

Unchanged by this work, and worth stating so the comparison is clean.

**Residuals** (`prepare_sorted_data_helper`). For every gene, in each group
separately, replicates are centred and Bessel-corrected:

```
raw :  r_ig = sqrt(n/(n-1)) * ( x_ig - mu_g )
log :  r_ig = sqrt(n/(n-1)) * ( log2(x_ig + c) - ghat_g ),   ghat_g = mean_i log2(x_ig + c)
```

Centring is **within group**, so the residuals carry no case-vs-control signal — a
gene's differential expression does not enter its own or anyone else's null. Genes are
then sorted by mean expression.

**Neighbourhood** (`gather_residuals_helper`). Starting at the gene closest in mean,
neighbour genes are added outward, nearest first, until `k_start` genes are in hand,
then in rounds of `k_step` while the relative change in the pool's mean absolute
residual stays ≤ `tau`, capped at `k_max` genes / `max_pool_size` residuals. Each gene
contributes exactly `n_rep` residuals, so **the pool is a sequence of contiguous
per-gene blocks** — which is what makes blocking possible without any extra
bookkeeping.

**Observed statistic.** `obs_g = mean_case − mean_control` on the residual scale
(linear for raw, log2 Fréchet means for log). Never modified by any null construction.

**p-value form.** Add-one corrected tail fraction,
`p = (#{null ≥ |obs|} + 1) / (#null + 1)`, so `p = 1` when `obs = 0` and the smallest
attainable p is `1/(#null + 1)`.

**BH.** Applied downstream in R, over tested genes only. Unchanged.

---

## 2. The current pooled null (`null_method = 0`), precisely

Let the case neighbourhood contain `k` genes, so the flat pool is
`P = (r_11 … r_1n | r_21 … r_2n | … | r_k1 … r_kn)`, `k·n` values, and similarly for
control. For each of `n_boot` draws:

1. draw `n` indices uniformly from `1 … k·n`, **independently**, and average those
   residuals → one bootstrapped case mean;
2. same, independently, from the control pool → one bootstrapped control mean;
3. the absolute difference is one null distance.

Because the draws in step 1 are independent, a single null mean can be built from a
residual of a very quiet gene and a residual of a very noisy one. **No real gene's mean
is formed that way**: a real `mean_case` for gene `g` averages `n` values that all come
from gene `g`, i.e. all share one noise level.

---

## 3. The blocked null (`null_method = 1`), precisely

Same, except step 1 becomes two steps:

1a. draw one neighbour **gene** `j` uniformly from `1 … k`;
1b. draw `n` residuals with replacement **from gene `j` only**, and average them.

Independently for the control side. So each null mean carries a single, coherent noise
level, drawn from the neighbourhood's distribution of noise levels — structurally the
same object as the observed mean it will be compared against.

Two consequences, one statistical (§4) and one computational (§5).

---

## 3b. Why the blocked null needs no `1/√n` correction

The `exact` model needs one because its null unit is a **single residual**, which has
variance `s²`, while the statistic it scores is a **mean of `n` replicates**, which has
variance `s²/n`. The null is built one scale up from the thing it describes, so it is
divided by `√n` after the fact — a correction that fixes the width but not the shape
(the documented source of the raw-normalisation problem).

Both bootstrap nulls avoid that entirely, because **the averaging is performed
explicitly**: a draw already *is* a mean of `n` values, so it lands on the observed
statistic's scale by construction. Concretely, for one gene with Bessel-corrected
residuals `r_1 … r_n` (which satisfy `Σ r_i = 0` and `(1/n)Σ r_i² = s²`, the unbiased
variance), resampling `n` of them with replacement gives a mean `M` with

```
E[M]   = mean(r) = 0                      (the null is centred, exactly)
Var[M] = ((1/n) Σ r_i²) / n = s² / n      (the variance of the observed mean, exactly)
```

Checked on a real example — replicates (4.1, 9.7, 6.0), so `s² = 8.11`:

| quantity | value |
|---|---|
| `Var[M]` from the enumeration | 2.70333333 |
| `s²/n` | 2.70333333 |
| ratio | **1.000000000000** |

and across neighbourhoods at `n_rep = 3, 5, 10`, the realised null sd divided by its
target is 0.9997 / 1.0005 / 0.9997. No correction is applied anywhere, and none is
needed.

### The "wouldn't the mean always be zero?" trap

It would — **if the resampling were without replacement.** Drawing all `n` residuals
without replacement is just a permutation of the same `n` values, so the mean is the
mean of all of them, which is exactly 0 by construction. The null would collapse to a
point mass at zero and every p-value would be 1. Verified: for
`r = (−3.0619, 3.7967, −0.7348)`, all 6 orderings give a mean of exactly 0.

**With** replacement the multiset can repeat values — `(r_1, r_1, r_1)` is a legal draw
— so the same three residuals produce 10 distinct means spanning −3.0619 … +3.7967. The
number of values averaged is still `n`, matching the observed mean; what varies is how
many times each one is used.

### What it does cost at low replicate counts

Only one balanced multiset (each residual used exactly once) has mean 0, but it carries
the largest weight, so a noticeable lump of the null sits exactly at zero:

| n_rep | P(M = 0 exactly) | excess kurtosis of the null, single-scale neighbourhood | sd / target |
|---|---|---|---|
| 3 | 0.222 | 1.30 | 0.9997 |
| 5 | 0.038 | 0.79 | 1.0005 |
| 10 | ~0 | 0.11 | 0.9997 |

These neighbourhoods have **one** noise level, so the excess kurtosis is not scale
mixing — it is the coarseness of resampling `n` values from `n` values. It is a genuine
approximation error at `n_rep = 3` and it decays quickly with replicates. Note it moves
the *shape*, never the width.

Splitting the two effects on the same neighbourhood structure and seed, `n_rep = 3`,
`k = 30`:

| neighbourhood | excess kurtosis | sd / target |
|---|---|---|
| single scale — discreteness only | 1.30 | 0.9997 |
| sd of log σ = 0.7 — both effects | 10.06 | 0.9999 |

So most of the blocked null's heavy tail is the scale mixture it is *supposed* to
reproduce (§4), with roughly 1.3 of excess kurtosis contributed by bootstrap
coarseness. Both push the same way — toward a slightly conservative test at
`n_rep = 3`, which is what the calibration simulation showed (`infl@.05` 0.62–0.91,
`infl@.01` 0.08–0.28). (This quantity is highly variable between neighbourhoods: it
depends on which σ's a given neighbourhood happens to contain — an earlier draw of the
same setup gave 3.96 rather than 10.06.)

---

## 4. Why the structure matters: it changes the shape, not the width

The two nulls have **the same variance**. Marginally, a residual drawn from the pooled
pool has variance `E[σ²]` over the neighbourhood, so a mean of `n` of them has variance
`E[σ²]/n` — exactly what a blocked draw has on average. What differs is what happens
*within* a draw:

- **pooled**: the `n` values averaged have `n` *different* scales, so averaging mixes
  the scales away. The result is pushed toward a normal with variance `E[σ²]/n`.
- **blocked**: the scale is constant within a draw and varies between draws — a genuine
  scale mixture.

And the thing the null has to describe *is* a scale mixture. Gene `g`'s observed
statistic has sd `σ_g·√(2/n)` for its own `σ_g`; across the genes being tested, `σ_g`
varies, so the ensemble of observed statistics is a mixture over `σ`.

Measured (one neighbourhood, `k = 30`, `n_rep = 3`, per-gene sd spread `sd(log σ) = 0.7`,
observed sds ranging 0.217–3.314, i.e. 15×; 4×10⁶ draws each):

| distribution | sd | excess kurtosis | P(\|D\| > 4·sd) |
|---|---|---|---|
| **(A)** ensemble of observed statistics — what must be described | 1.3330 | **17.69** | 8.0e-3 |
| **(B)** pooled null (current) | 1.2407 | **0.60** | 1.9e-4 |
| **(C)** blocked null (new) | 1.2396 | **3.96** | 3.4e-3 |

The widths agree to under 1 %. The shapes do not: the pooled null is almost Gaussian
(excess kurtosis 0.6) while the thing it is scoring has excess kurtosis 17.7, and in
the far tail the pooled null is **42× too light** where the blocked null is 2.4× too
light.

Scoring (A) against each null:

| null used | FPR @ 0.05 | FPR @ 0.01 |
|---|---|---|
| pooled | 0.0587 | **0.0277** (2.8× target) |
| blocked | 0.0438 | 0.0154 (1.5× target) |
| target | 0.0500 | 0.0100 |

This is the mechanism behind the calibration result reported in
`raw_normalisation_diagnosis.md` §7c — the pooled bootstrap running `infl@.01` of
1.4–2.0 while blocking restores it. **The tail is the region BH reads**, which is why
this matters more than the α = 0.05 column.

**Blocking improves the shape; it does not perfect it** (3.96 vs 17.69). Two honest
reasons: the mixing distribution available to the null is only the `k = 30` empirical
`σ_j` of that neighbourhood, not the continuum; and each `σ_j` is itself estimated from
`n_rep` residuals, which at `n_rep = 3` is 2 degrees of freedom. The blocked null is a
better model of the right object, not a correct one.

---

## 5. Why the blocked null can be enumerated

Resampling `n` values with replacement from `n` values has `n^n` ordered outcomes, but
the *mean* only depends on **how many times each residual was drawn**, not on the order.
So the distinct outcomes are the multisets — `C(2n−1, n)` of them — each carrying a
multinomial weight `n! / ∏_j m_j!`:

| n_rep | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|
| `n^n` | 4 | 27 | 256 | 3,125 | 46,656 | 16.8 M | 10¹⁰ |
| `C(2n−1, n)` | 3 | **10** | 35 | 126 | 462 | 6,435 | 92,378 |

Worked, for one gene with residuals `r = (−1.2, 0.1, 1.1)` and `n_rep = 3`:

| multiset | draw counts | mean | weight |
|---|---|---|---|
| 111 | 3,0,0 | −1.2000 | 1 |
| 112 | 2,1,0 | −0.7667 | 3 |
| 113 | 2,0,1 | −0.4333 | 3 |
| 122 | 1,2,0 | −0.3333 | 3 |
| 123 | 1,1,1 | 0.0000 | 6 |
| 222 | 0,3,0 | 0.1000 | 1 |
| 133 | 1,0,2 | 0.3333 | 3 |
| 223 | 0,2,1 | 0.4333 | 3 |
| 233 | 0,1,2 | 0.7667 | 3 |
| 333 | 0,0,3 | 1.1000 | 1 |

Weights sum to 27 = `n^n`; the weighted mean is 0 (the residuals are centred) and the
weighted variance is 0.2956 = `σ̂²/n` exactly — i.e. the enumeration reproduces the
sampling variance of a mean of `n` draws, as it must.

Ten values instead of 27, and 462 instead of 46,656 at `n_rep = 6`. `C(2n−1, n)` is
computed as the running product `C(n+j, j)`, integer at every step, with an int64
overflow guard (`n_multisets_helper`); the tuples are walked as a non-decreasing
odometer and the weights built as a product of binomials so every intermediate stays an
exact integer (`enumerate_gene_means_helper`).

Doing this for every gene block in the pool gives the **complete** blocked null as
`k · C(2n−1, n)` weighted values per side (`build_blocked_means_helper`). Only complete
blocks are used — a partial trailing block can only come from the `max_pool_size` cap
and is dropped.

---

## 6. Scoring it

With case values `a_i` (weights `w^a_i`) and control values `b_j` (weights `w^b_j`):

```
p = ( Σ_ij w^a_i w^b_j · 1[ |a_i − b_j| ≥ |obs| ]  +  1 ) / ( W_a · W_b + 1 ),
        W = Σ w = n_genes_pool · n_rep^n_rep
```

Computed exactly as the individual-residual null already is — sort the control values
once, prefix-sum their weights, then per case value two binary searches give the weight
of control values inside `(a_i − t, a_i + t)`; the rest is tail weight
(`compute_pvalue_blocked_exact_helper`). Cost is one sort plus `N` binary searches over
`N = k·C(2n−1,n)` values, not the `N²` pairs it scores.

No RNG, no draw count, no seed. The result is the `n_boot → ∞` limit of the sampled
blocked bootstrap — confirmed against it at B = 2×10⁶ (0.1002 vs 0.1001, 0.01130 vs
0.01120, 9.25e-4 vs 9.18e-4), and the Fortran agrees with an independent R
implementation to 1.65e-15 relative over 400 genes.

---

## 7. What changes operationally

**The p-value floor.** Now differs per arm — compare `min_p` against the arm's own:

| construction | floor | at `k = 30`, `n_rep = 3` |
|---|---|---|
| exact | `1/(n_pool_case · n_pool_control + 1)` | 1.2e-4 |
| sampled (pooled, or blocked above the cap) | `1/(n_boot + 1)` | 1.0e-4 |
| blocked & enumerated | `1/(W_a · W_b + 1)` | **1.5e-6** |

BH at `q` needs the smallest p to reach `q/G` — 4.2e-6 at `G = 12,000`, `q = 0.05`. Only
the enumerated arm can clear that from a single gene at `n_rep = 3`. `W` is derivable in
R from the reported neighbourhood sizes and the replicate count, so no new output array
was added.

**Determinism.** The enumerated path uses no RNG at all, so it is reproducible
independently of gene order and could be parallelised over genes as-is (it is left
sequential so one code path serves both).

**RNG cost, where sampling is still used.** `random_number` is now called once per
*chunk* of draws rather than once per draw — one array fill covering `chunk · per_iter`
values, buffer allocated once and reused across genes. At typical replicate counts that
is one RNG call per gene instead of 10,000. Verified behaviour-preserving: the pooled
arm reproduces the pre-change module **bit for bit** (`max |Δp| = 0`), because filling
one large array consumes the same stream in the same order as many small fills.

**Runtime** (2,000 genes, `k_max = 50`, gfortran -O2, single core, vs this module's own
10,000-draw pooled bootstrap on the same data):

| n_rep | values/side | blocked | pooled | ratio |
|---|---|---|---|---|
| 3 | 500 | **0.41 s** | 1.45 s | **0.28×** |
| 4 | 1,750 | 1.71 s | 1.79 s | 0.95× |
| 5 | 6,300 | 7.34 s | 2.18 s | 3.4× |
| 6 | 23,100 | 30.6 s | 2.63 s | 12× |
| 7 | 85,800 | 131 s | 3.03 s | 43× |

At `n_rep = 3` — the case where the floor actually bites — the blocked null is **3.5×
faster than the bootstrap it replaces and removes the floor**. Cost grows like
`C(2n−1, n)`, ~4× per extra replicate, so `BLOCKED_MAX_ENUM_VALUES = 8000` enumerates
through `n_rep = 5` at `k_max = 50` and falls back to the **sampled** blocked null above
that (~15 % more expensive than pooled, and back to the `1/(n_boot+1)` floor). The cap
is a runtime knob, not a memory one — 8,000 values is 64 KB per side.

**`trim_frac` is ignored under `null_method = 1`.** Trimming sorts the pool in place,
which destroys the per-gene block layout the blocked null reads. The alloc layer forces
it to 0 there; trim + blocked is not a valid combination.

---

## 8. What this does *not* do

- **It does not add information.** The enumerated floor of 1.5e-6 is not 1.5e-6 worth of
  evidence. The effective sample size behind a `k·n_rep = 90`-residual pool at
  `n_rep = 3` is of order 10² (see `raw_normalisation_diagnosis.md` §7b), so a p-value
  far below ~1e-2 there reflects which neighbourhood the gene happened to draw. What
  enumeration removes is an *artificial* floor that was binding well above the
  information limit — worth doing, and not the same as resolving finer.
- **It does not fix the raw normalisation.** The anti-conservatism diagnosed in §1–2 of
  `raw_normalisation_diagnosis.md` is about the *exact* model's individual-residual
  null under raw normalisation. Blocking addresses the mean-level null's tail; the
  raw-vs-log question is separate and still points at log.
- **It is not yet confirmed on real data.** Every calibration claim about blocking so
  far is simulation (`null_construction_comparison.R`). The TCGA split-half run through
  `TOX-raw-blocked` / `TOX-log-blocked` is what tests it.

---

## 9. Summary table

| | exact | pooled (`null_method = 0`) | blocked (`null_method = 1`) |
|---|---|---|---|
| null unit | individual residual pair, ÷√n | mean of `n` iid pooled residuals | mean of `n` residuals from **one gene** |
| scale within a draw | one residual | mixed across genes | constant, from one gene |
| shape vs observed ensemble | too heavy (individual residuals) | too light (excess kurtosis 0.6 vs 17.7) | closer (3.96 vs 17.7) |
| measured FPR @ .01 | — | 0.0277 | 0.0154 |
| computation | exact pairwise count | `n_boot` draws | exact enumeration, else `n_boot` draws |
| RNG | none | yes | none (enumerated) / yes (fallback) |
| floor at `k=30`, `n_rep=3` | 1.2e-4 | 1.0e-4 | **1.5e-6** |
| runtime at `n_rep = 3` | fastest | 1.45 s | 0.41 s |
| `trim_frac` honoured | yes | yes | no |
