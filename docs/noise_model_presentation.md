# The TOX Noise Model — Presentation Reference

A speaker's reference for a **20–30 min talk**, pitched at a **bioinformatics undergrad**.
Not slides — a "don't forget anything" cheat sheet. Each section has the **idea**, the
**plain-language version**, and (where useful) the **exact detail** if someone asks.

> Golden thread of the talk: *"We built a data-driven way to decide whether a gene's
> expression difference is real or just noise — and we showed it's as well-behaved as the
> best standard tools, without their failure mode."*

---

## 0. One-sentence pitch (open + close with this)

> The noise model detects **outlier genes** — genes whose expression differs between two
> groups (e.g. cancer vs healthy) **more than measurement + biological noise alone would
> produce** — by comparing each gene against an **empirical noise distribution built from
> the data itself**, instead of assuming a fixed statistical model.

---

## 1. The problem & why it matters (≈2 min)

- **Goal:** given expression for many genes in two conditions (e.g. tumour vs normal),
  find the genes that are *genuinely* different — the biologically interesting "outliers".
- **The catch:** every measurement is noisy. Two random halves of the *same* group already
  differ gene-by-gene just by chance. So "A ≠ B" is not enough — we need "A differs from B
  **more than noise would**".
- **Why not just use a t-test / edgeR / DESeq2?** Those assume a *parametric* model of the
  noise (a specific distribution). When that assumption is slightly wrong — which it is on
  real human data — they call too many false positives (we show this later).
- **Our angle (Tensor-Omics philosophy):** stay **data-driven and geometric** — measure how
  far a gene's group-difference is from the spread of "noise differences" seen in comparable
  genes, with **no assumed distribution**.

**Analogy to say out loud:** *Is a 5-point exam gap between two classes meaningful? Don't
assume a formula for "normal" gaps — look at the gaps you actually see between random halves
of students with similar ability. If 5 points is bigger than almost all of those, it's real.*

---

## 2. Background concepts the audience needs (≈3–4 min)

Define these up front; the rest of the talk leans on them.

- **RNA-seq counts:** per gene, per sample, an integer count ∝ how much that gene was
  expressed × how deeply that sample was sequenced (**library size**).
- **Two kinds of noise:**
  - *technical / shot noise* — from finite sampling of reads (Poisson).
  - *biological noise* — real variation between patients (adds **overdispersion**).
- **Replicates (`n_rep`):** how many samples per group. More replicates = more certainty.
- **The null hypothesis H0:** "this gene is NOT differentially expressed." A good method,
  run on data where H0 is true everywhere, should almost never reject it.
- **p-value:** probability of seeing a difference this big *if H0 were true*. Under H0, a
  well-behaved test gives p-values **uniform on [0,1]** — every value equally likely.
- **FPR (false-positive rate) @ α:** fraction of genes with p < α when H0 is true. Should
  ≈ α (so ~5% below 0.05). **Too many = anti-conservative** (bad); too few = conservative.
- **FDR (false discovery rate) & BH:** after testing thousands of genes you correct for
  multiple testing (Benjamini–Hochberg). Key fact: **under a complete null, controlling FDR
  becomes controlling FWER — so ANY discovery at all is a false positive.**
- **Parametric vs empirical null:**
  - *parametric* = "assume the noise follows distribution X, compute the tail." (edgeR/DESeq2)
  - *empirical* = "estimate the noise distribution from the data, read the tail off that."
    (our noise model) — robust when no clean formula fits.

---

## 3. How the noise model works (≈6–8 min — the core)

Walk this in five beats.

### Beat 1 — The observed statistic
For a gene, take the **difference of group means** (case vs control), on the residual scale:
- raw normalization: `mean_case − mean_control` (linear).
- log normalization: difference of **log-space (Fréchet) means**, i.e. means of `log2(x+1)`.

That number is "how far apart the two groups are for this gene". Question: **is it big
compared to noise?**

### Beat 2 — Residuals = the shape of noise
For each gene, subtract its own mean → **residuals** (`x_i − mean`). Residuals are the
gene's sample-to-sample wiggle = its noise. We collect residuals to learn "what does noise
look like at this expression level?"

- **Bessel correction** (detail): mean-subtracted residuals *understate* the true spread
  (you spent one degree of freedom estimating the mean). At `n=3` the residual SD is only
  ~82% of the truth. We scale every residual by `√(n/(n−1))` to fix this. Matters most at
  **low replicate counts** — exactly our regime. *(Say: "small technical correction so we
  don't underestimate noise with few samples.")*

### Beat 3 — The neighbourhood (kNN in mean-expression space)
Noise depends on expression level (low-count genes are noisier relatively; high-count genes
have big absolute variance). So for a target gene we pool residuals from **genes with a
similar mean expression** — its *k nearest neighbours in mean space*.

- Controlled by `k_start` (initial # genes), `k_step`, `k_max` (cap), `tau` (stop when the
  pooled noise stops changing). Counted in **genes**, not residuals.
- We build **two independent pools**: one from the case side, one from the control side.

**Analogy:** *compare a gene only against other genes "in its weight class".*

### Beat 4 — The empirical null + the p-value
Now: how big a case-vs-control difference does *noise alone* produce? Two interchangeable
engines (same interface, pick one):

- **`exact` model** (used for the TCGA results): scale each pool by `1/√n_rep` (because the
  observed statistic is a difference of *means*, but the pool is *individual* residuals —
  means are `√n_rep` less variable), then count the fraction of all pairwise
  `|r_case − r_control|` that reach or exceed the observed difference. Deterministic, fast
  (sorted + binary search). *Mildly conservative by design* (keeps heavy individual-residual
  tails).
- **`bootstrap` model:** for 10,000 draws, resample `n_rep` residuals per side, average each,
  take `|Δmean|` → a null distance. The 10,000 distances are the null. Correct tail *shape*.

**p-value** = add-one-corrected tail fraction:
`p = (#{null ≥ observed} + 1) / (N + 1)`. Small p ⇒ the gene's difference is bigger than
noise ⇒ candidate outlier.

### Beat 5 — Two normalizations + an optional knob
- **raw** vs **log** are the two main modes (log stabilises the mean–variance relationship).
- **Tail trimming (raw only, optional):** drop the most extreme residuals (e.g. top/bottom
  5%) from the pool before scoring, so a few outliers don't inflate the noise estimate.
  (Off by default; a knob we're studying.)

**Reproducibility note:** the bootstrap RNG is seeded once (`init_random(42)`); the exact
model is fully deterministic.

**Slide-worthy diagram to draw:** genes sorted by mean → highlight target gene → bracket its
kNN neighbours → pool their residuals (case & control) → histogram = the empirical null →
mark the observed Δmean → shaded tail = p-value.

---

## 4. The competitors we compare against (≈2 min)

All three are the field-standard differential-expression (DE) tools, run at their **current
recommended settings** (important for credibility — reviewers can't say we used them wrong):

- **edgeR** — Negative-Binomial GLM, quasi-likelihood F-test (`glmQLFit`/`glmQLFTest`, robust).
- **DESeq2** — Negative-Binomial GLM, Wald test, median-of-ratios normalization.
- **limma-voom** — log-transform + **empirical per-gene variance** + moderated t/F
  (`voomLmFit` + robust `eBayes`).

**One structural distinction to plant early (it explains the results):**
- edgeR & DESeq2 assume a **parametric NB variance** `Var = μ + φ·μ²` (a *dispersion* φ).
- limma estimates variance **empirically per gene** (no rigid parametric tail).
- TOX estimates the whole noise distribution **empirically** (most nonparametric of all).

---

## 5. How we test it (≈3–4 min)

Three complementary tests. Emphasise the first.

### 5a. Null calibration on real TCGA data (the decisive one)
- Take **one homogeneous group** (e.g. all Stage I tumours of a cancer). Randomly split its
  samples into fake groups **A vs B**. Same population, random labels ⇒ **H0 is exactly true,
  no gene is really DE.**
- Run every method A-vs-B. A correct method must give: **p-values uniform**, **FPR ≈ α**,
  **~0 discoveries at FDR<0.05**.
- Repeat **100 times** (random splits) per cell and average — so it's not one lucky split.
- Swept across: **8 cancers**, **4 stages + healthy**, replicate counts **10/20/40 + full
  half**, and a **kNN-parameter grid** (neighbourhood sizes), raw & log.
- Key reading: it isolates **calibration** (false positives), *not* power.

### 5b. Simulation (calibration_test.R)
- Generate data with a *known* truth. **Part A** = mock null (no DE) → p-values must be flat.
  **Part B** = with DE genes → check observed vs nominal FDR and ranking (AUC).
- Stress arms: `nb`, `lnpois` (lognormal-Poisson), `tpois` (heavy-tailed), `bimodal`
  (genes switch on/off — violates every tool's assumptions on purpose).

### 5c. Count-distribution diagnostic (count_distribution.R)
- Sanity check on the *data*: is it Poisson or Negative Binomial? Per gene, fit Poisson vs NB
  (with a library-size offset), likelihood-ratio test → is it overdispersed?
- Purpose: justify *why* the NB tools are the right comparison, and set the stage.

**Metrics glossary for the audience:** `FPR@0.05` (raw tail, target 0.05); `hits_FDR05`
(fraction of genes called DE after BH — target **0**, because H0); `ks_D / ad_A2 / cvm_W2`
(distances of the p-value distribution to uniform — 0 = perfect); `median_p` (should be 0.5).

---

## 6. Results (≈5–6 min — the payoff)

### 6a. The data is Negative Binomial (context)
- Across TCGA cohorts, **~100% of genes reject Poisson**; typical dispersion **φ ≈ 0.25–0.42**,
  i.e. **BCV (biological CV) ≈ 0.5–0.65**. Right in the expected range for human patient
  cohorts. ⇒ Poisson is far too narrow; **NB is the correct data model** — which is exactly
  why edgeR/DESeq2 use it. Good sanity check.

### 6b. Null calibration — the headline
Representative per-method calibration (H0, so target FPR ≈ 0.05, hits_FDR05 ≈ 0):

| method   | FPR@0.05 | hits_FDR05 | verdict |
|----------|---------:|-----------:|---------|
| limma    | ~0.049   | 0.000      | **calibrated** |
| TOX-log  | ~0.053   | ~0.0006    | **calibrated** (≈ limma) |
| TOX-raw  | ~0.082   | 0.000      | anti-conservative (p not uniform) — *open* |
| edgeR    | ~0.102   | 0.0036     | **anti-conservative (~2×)** |
| DESeq2   | ~0.097   | 0.0071     | **anti-conservative (~2×)** |

**The three big findings:**
1. **A two-family split.** The NB-dispersion tools (**edgeR, DESeq2**) are ~2× anti-conservative
   and make real false discoveries. The empirical-variance methods (**limma, TOX-log**) are
   well-calibrated. Mechanism = whether variance is forced into a parametric form or estimated
   from data.
2. **edgeR/DESeq2 get *worse* with more replicates** — p-value histograms drift further from
   uniform and they reject H0 in more and more of the 100 runs as n grows. **limma and TOX
   stay flat** at 10/20/40/150 samples. Under a true null, FPR growing with n = the model is
   misspecified and larger n just gives more power to detect that it's wrong. *(This is a
   documented effect — Li et al., "Exaggerated false positives…", Genome Biology 2022.)*
   → **This is the key figure: calibration vs. sample size, per method.**
3. **Consequence for the field:** because the "reference" DE tools over-call on real data, the
   gene lists they produce contain false positives — so they are **not clean ground truth**,
   especially in large studies.

### 6c. The honest open item — TOX-raw
- **TOX-log** is essentially as calibrated as limma. **TOX-raw** is **anti-conservative**
  (FPR ≈ 0.08, p-values inflated/left-shifted — *not* uniform). `hits_FDR05 = 0` for raw only
  means the excess isn't extreme enough to cross BH at ~15k genes — it's **not** a calibration
  pass.
- **Cause — now diagnosed (this is a nice "we chased it down" story):** in raw space a
  mean-neighbourhood is mean-homogeneous but **not variance-homogeneous**, so the pooled
  residuals are a **scale mixture**: the total variance is *correct*, but the **shape** is wrong —
  a **narrow core with heavy tails**. A moderate observed distance lands too far out in that core
  (anti-conservative at α = 0.05) while an extreme one is still covered (conservative at α = 0.01).
  That two-sided pattern is exactly what the data shows, and it is the one pattern a *width* error
  cannot produce. Log space stabilises the mean–variance trend, so TOX-log is unaffected.
- **Two candidates were ruled OUT by the data** (say this — it's the strongest part):
  - **Neighbourhood size (`k`) is not the cause.** Larger `k` helps consistently (128/128 paired
    comparisons improve, so `k20_50` is the best config) but only closes about **one seventh** of
    the gap — the raw excess is +0.64 while `k_max` 30→50 moves it by 0.005.
  - **The null is not too narrow** — a width error can't be anti-conservative at 0.05 *and*
    conservative at 0.01 simultaneously.
- **The fix — the gene-blocked null:** draw all `n_rep` residuals from **one neighbour gene**
  instead of from the pooled mixture, so every null mean carries a single coherent noise level.
  That changes the null's **shape, not its width**. In simulation it moves FPR@0.01 from 0.0277
  (pooled) to 0.0154, and — because resampling `n` of a gene's `n` residuals has only `C(2n−1, n)`
  outcomes — it can be **enumerated exactly**: no RNG and no p-value floor.
- *Frame as: "known, diagnosed, fix implemented — real-data validation is the next run."*

### 6d. Bottom line
> **Proof of concept achieved.** TOX-log matches the best-behaved standard tool (limma) on
> real-data calibration and beats the widely-used NB tools (edgeR/DESeq2), which degrade with
> sample size. It does this **nonparametrically**, which is the whole point — and the
> foundation we can now extend.

*(Caveat to state: this is a **calibration** result — controlling false positives under H0 —
not yet a power/precision claim. Power is the complementary evaluation, next.)*

---

## 7. Open questions & next steps (≈1–2 min)

- **Validate the gene-blocked null on real data.** It is implemented and verified in simulation
  (and against an independent R implementation to 1.65e-15); the TCGA split-half run is what
  tests it. Watch `FPR_0.01` for `TOX-raw-boot` vs `TOX-raw-blocked` — simulation says pooled is
  anti-conservative there (1.4–2.0×) and blocking restores it.
- **kNN retuning is closed** — `k20_50` fixed as the best config; it is not the raw fix.
- **Power / FDR-with-signal** evaluation (Part B) — show TOX keeps sensitivity, so "calibrated"
  isn't just "toothless".
- **Planned noise-model extensions** — a well-calibrated null is the base they build on (name
  your specific extensions here).
- Confirm the **edgeR/DESeq2 degrade-with-n** curve across all 8 cancers (final headline figure).

---

## 8. Takeaway soundbites (memorise 2–3)

- "We ask *is this bigger than noise?* — and we learn what noise looks like **from the data**,
  not from an assumption."
- "On real human data, the standard NB tools invent false positives, and it gets **worse** the
  more samples you have. Our empirical null doesn't."
- "TOX-log is as calibrated as the best standard method (limma) and beats edgeR/DESeq2 — with
  no distributional assumptions."
- "Calibration first, power next — you can't trust discoveries from a miscalibrated null."

---

## 9. Glossary / quick numbers (keep on a backup slide)

- **Residual:** `x_i − mean` (per gene). The noise signal.
- **Empirical null:** noise distribution estimated from data (vs assumed).
- **kNN neighbourhood:** genes with similar mean expression, pooled to estimate local noise.
- **Bessel factor:** `√(n/(n−1))`, fixes residual under-dispersion (big at n=3: ~82%→100%).
- **Three null constructions:** **exact** = deterministic pairwise tail, ÷√n (mildly
  conservative); **pooled bootstrap** (`null_method=0`) = mean of `n` iid residuals drawn from the
  whole pool (10,000 draws); **gene-blocked** (`null_method=1`) = mean of `n` residuals from **one
  gene** — right shape, exactly enumerable (no RNG, no floor). Floors at k=30, n_rep=3:
  1.2e-4 / 1.0e-4 / **1.5e-6**. `trim_frac` is ignored under blocked.
- **Overdispersion / φ / BCV:** `Var = μ + φ·μ²`; `BCV = √φ`. TCGA: φ≈0.25–0.42, BCV≈0.5–0.65.
- **FPR@0.05:** target 0.05. **hits_FDR05:** target 0 (any = false discovery under H0).
- **KS / AD / CvM (`ks_D/ad_A2/cvm_W2`):** distance of p-values to uniform; 0 = perfect.
- **Reference tools:** edgeR (NB QL-F), DESeq2 (NB Wald), limma-voom (log + empirical variance).
- **Datasets:** TCGA — 8 cancers, 4 stages + matched healthy.

---

## 10. Likely questions (prep answers)

- **"Why not just use edgeR/DESeq2?"** → they over-call under H0 on real data and worsen with
  n (Li et al. 2022); we want a method that doesn't rely on the NB tail.
- **"Isn't 'well-calibrated' just being conservative?"** → No — p-values are *uniform* and
  median ≈ 0.5, not skewed high; and power is the next evaluation.
- **"Why is TOX-raw off?"** → the raw mean–variance trend makes the neighbourhood
  variance-heterogeneous, so the pooled null is a **scale mixture** — right width, wrong
  *shape* (narrow core, heavy tails). Not a `k` problem (ruled out: k closes ~1/7 of the gap)
  and not a width problem. Fixed by the **gene-blocked null**; log is unaffected.
- **"Is the half-split really H0?"** → yes: same population, random labels → no true DE by
  construction; the cleanest possible null, no effect-size threshold needed.
- **"Poisson or NB?"** → strongly NB (100% of genes reject Poisson; BCV≈0.5) — as expected
  for human cohorts.
- **"Exact vs bootstrap — which and why?"** → results use exact (deterministic, fast);
  bootstrap fixes the exact model's mild conservatism if needed.

---

## 11. Suggested slide plan & timing (~25 min, ~14 slides)

1. Title + one-line pitch — 0:30
2. The problem: real vs noise difference (exam analogy) — 2:00
3. Background: counts, noise types, H0, p-values, FPR/FDR, empirical vs parametric — 3:30
4. Core idea in one picture (observed Δ vs empirical null) — 1:30
5. Residuals + Bessel — 2:00
6. kNN neighbourhood (the "weight class" diagram) — 2:00
7. The null + p-value; exact vs bootstrap — 2:30
8. raw vs log (+ trimming knob) — 1:30
9. The reference tools & the parametric-vs-empirical distinction — 2:00
10. How we test: the H0 half-split (+ simulation, count-dist) — 3:00
11. Result: data is NB (φ/BCV) — 1:30
12. Result: calibration table + two-family split — 2:30
13. **Result: calibration-vs-n (the money figure)** — 2:30
14. TOX-raw open item + next steps + takeaways — 2:00

Trim slides 3/5/8 first if running long; slides 12–13 are the payoff — protect their time.

---

## 12. What NOT to overclaim (keep yourself honest)

- It's a **calibration** result, not yet a power/precision benchmark.
- **TOX-raw is not calibrated yet** — say so; it's the open, understood item.
- Don't say limma is *immune* to inflation — it's *more robust*; at extreme n it can inflate
  mildly too (just far less). Claim only your tested range (≤150 samples).
- The TCGA calibration numbers are from the **exact** model; the bootstrap model is an
  alternative engine.
- Pooled median FPRs understate the story — the **with-n trend** is the real evidence.

---

# APPENDIX A — Figure specifications (what to draw, where, how)

Detailed, buildable specs for every figure. Each gives **Purpose / Type / Layout / Axes /
Elements & positions / Annotations / Colour / Source**. Schematics (A1–A7) are hand-built
(Illustrator / draw.io / PowerPoint / TikZ); data figures (A8–A12) come from the scripts'
outputs or the CSVs.

## A0. One visual system first (apply to ALL figures)

Fix these once so the deck reads as a single system:

- **Fixed method colours** (colour-blind-safe, Okabe–Ito), used identically in *every* figure
  and in the text:
  - **TOX-log** = blue `#0072B2` (our hero — always most saturated)
  - **TOX-raw** = sky blue `#56B4E9`
  - **limma** = green `#009E73` (the calibrated reference)
  - **edgeR** = vermilion `#D55E00` (the "bad" one — warm = alarm)
  - **DESeq2** = orange `#E69F00` (also NB family → warm)
- **Semantic colours** reused across schematics: *case* = a warm grey-blue, *control* = a warm
  grey-orange; *noise / null* = light grey fill; *observed statistic* = black; *"good" region*
  = faint green band, *"bad"* = faint red band.
- **Type:** one sans-serif (e.g. Inter/Helvetica), ≥ 24 pt on slides; axis titles bold.
- **Conventions:** dashed line = a *target/reference* value; solid = data; shaded tail = a
  p-value / rejection region; always label axes with units and keep the same axis ranges when
  two panels are compared.
- **Theme:** white background, minimal gridlines (light grey), no chartjunk, no 3-D.

---

## A1. Pipeline / functionality overview (the "how it works in one picture" schematic)

- **Purpose:** the whole method end-to-end on one slide; you'll return to it as a roadmap.
- **Type:** left-to-right flow diagram, 5 stages in rounded boxes joined by arrows.
- **Layout (left → right), each box stacked with a tiny inset sketch above the label:**
  1. **Input** — box: "Expression matrix (genes × samples), two groups: case / control".
     Inset: a small heatmap grid, columns tinted case-blue vs control-orange.
  2. **Residuals** — box: "Per gene: subtract mean → residuals; ×√(n/(n−1)) (Bessel)".
     Inset: a dot cloud around a horizontal mean line, arrows showing deviations.
  3. **Neighbourhood** — box: "Pool residuals from genes with similar mean (kNN in mean space)".
     Inset: genes as dots on a horizontal "mean expression" axis, a bracket around the target's
     neighbours.
  4. **Empirical null** — box: "Case vs control noise-difference distribution (exact or
     bootstrap)". Inset: a grey histogram.
  5. **p-value** — box: "Observed Δmean vs null tail → p". Inset: same histogram with a black
     vline and a shaded right tail.
- **Elements & positions:** boxes equal size, evenly spaced; arrows left→right; under the whole
  row a thin caption strip: "no distributional assumption — the null is learned from the data."
- **Annotations:** above stage 2–3 a small brace labelled "learn what noise looks like *here*";
  above stage 5 "small p ⇒ bigger than noise ⇒ candidate outlier".
- **Colour:** boxes white with grey borders; insets use the semantic palette (case-blue,
  control-orange, null-grey, observed-black).
- **Source:** hand-drawn schematic.
- **Undergrad touch:** number the boxes 1–5; you literally walk the arrow.

## A2. Core concept — observed difference vs the empirical null (THE key conceptual figure)

- **Purpose:** the single idea the whole talk rests on. Spend time here.
- **Type:** one histogram + a vertical marker.
- **Layout:** single centred panel, wide (16:9-friendly), lots of whitespace.
- **Axes:** x = "case−control noise difference (|Δ| on residual scale)", starting at 0;
  y = "frequency (from noise)". Keep y unlabelled-numbers minimal (it's illustrative).
- **Elements & positions:**
  - A right-skewed grey histogram/density = **the empirical null** (noise-only differences),
    centred near 0, tapering right.
  - A **black vertical line** at the observed |Δmean|, placed **out in the right tail**
    (~90th–97th percentile) so the shaded area is visibly small.
  - **Shade the area to the right** of the black line in faint red = the p-value.
- **Annotations:** label the black line "observed difference (this gene)"; a callout on the red
  tail "p = fraction of noise ≥ observed"; a small note under the bulk "noise alone rarely gets
  this far". Optionally a second, greyed-out black line near the centre labelled "a typical
  (non-DE) gene" for contrast.
- **Colour:** null = grey fill; observed = black; tail = faint red (`#D55E00` at ~20% alpha).
- **Source:** schematic (or generate a real null with the bootstrap for one gene and overlay).
- **Build tip:** make TWO versions on click — first the histogram alone ("what noise looks
  like"), then add the observed line + tail ("now score a gene").

## A3. Residual construction (+ Bessel)

- **Purpose:** define "residual" concretely and motivate the Bessel factor.
- **Type:** small scatter/strip plot for one gene.
- **Layout:** single panel; samples on the x-axis, expression on y.
- **Axes:** x = "samples (1..n)"; y = "expression (raw or log2)".
- **Elements & positions:** n dots (say n=3–5) scattered vertically; a **horizontal solid line
  at the mean**; **vertical dashed segments** from each dot to the mean line = the residuals,
  labelled "residual = xᵢ − mean".
- **Annotations:** a boxed note: "with few samples, residuals *understate* the true spread
  (n=3 → ~82%). Fix: ×√(n/(n−1))." Show a faint wider band (the corrected spread) behind the
  raw residuals.
- **Colour:** dots case-blue; mean line black; correction band light green.
- **Source:** schematic with a toy 3–5 point example.

## A4. kNN neighbourhood in mean space (the "weight class" figure)

- **Purpose:** show that noise is estimated *locally* in expression.
- **Type:** 1-D dot strip + a pulled-out residual pool.
- **Layout:** top strip spanning the width = all genes; a callout box below.
- **Axes:** top strip x = "mean expression (sorted)"; no y (jitter dots for visibility).
- **Elements & positions:**
  - Many small grey dots along the strip = genes by mean.
  - **Highlight the target gene** (black, larger) somewhere mid-strip.
  - Draw a **bracket / shaded window** around its `k` nearest neighbours.
  - An arrow from the window down to a **pool box** showing a little histogram of the pooled
    residuals (built from those neighbours).
- **Annotations:** window labelled "k nearest genes in mean (k_start…k_max, stop at τ)";
  note "two independent pools: case & control"; caption "compare a gene only against its
  'weight class'."
- **Colour:** target black; window faint blue; pool histogram grey.
- **Source:** schematic.
- **Optional real inset:** a real mean-vs-variance scatter (from a cohort) with the window
  overlaid to show neighbours really do share a noise level.

## A5. Building the null — exact vs bootstrap (two engines)

- **Purpose:** show the two interchangeable ways to turn pools into a null; keep it light.
- **Type:** two side-by-side mini-schematics sharing one output histogram style.
- **Layout:** left half "exact", right half "bootstrap", a divider line; both feed a null
  histogram at the bottom.
- **Elements & positions:**
  - **Left (exact):** two small stacks of residuals (case, control) → "scale by 1/√n_rep" →
    "count all pairwise |r_case − r_control| ≥ observed". Tag: "deterministic, fast, mildly
    conservative".
  - **Right (bootstrap):** "resample n_rep per side → average → |Δ|", drawn as a loop arrow
    with "×10,000" → collect into the null. Tag: "correct tail shape, uses RNG (seed 42)".
- **Annotations:** a brace under both: "same interface → pick either; TCGA results use exact".
- **Colour:** case-blue / control-orange stacks; null-grey output.
- **Source:** schematic.
- **If short on time:** cut this to a single line on the previous slide.

## A6. raw vs log — mean–variance stabilisation

- **Purpose:** explain *why* log behaves and raw is trickier (sets up the TOX-raw result).
- **Type:** two scatter panels, shared layout.
- **Layout:** side by side, **identical axis style**, left "raw", right "log2".
- **Axes:** left x = "gene mean (raw)", y = "gene variance (raw)"; right x = "mean log2",
  y = "variance log2".
- **Elements & positions:** each panel = a cloud of genes (dots). Left: a steep upward
  **curved trend** (variance grows with mean, `Var = μ + φμ²`) — draw the fitted curve.
  Right: a **roughly flat** cloud (variance stabilised) — draw a near-horizontal trend.
- **Annotations:** left: "steep trend ⇒ a mean-neighbourhood mixes very different variances
  ⇒ null mis-scaled". Right: "flat ⇒ neighbours share variance ⇒ null well-scaled".
- **Colour:** dots grey; trend lines black; a faint red arrow on the left panel pointing at the
  steep part ("this is what hurts TOX-raw").
- **Source:** real data — per-gene mean vs variance from one cohort (raw counts/TPM and log2).
  Trivial to make in R (`plot(rowMeans, rowVars)`); a LOESS line for the trend.

## A7. Parametric NB vs empirical variance (why the methods differ)

- **Purpose:** the conceptual reason edgeR/DESeq2 fail where limma/TOX don't.
- **Type:** 3-column concept strip (one column per "family").
- **Layout:** three equal columns: **edgeR/DESeq2**, **limma**, **TOX**.
- **Elements & positions (top icon + one line each):**
  - edgeR/DESeq2: icon = a *fixed formula* `Var = μ + φμ²` with a rigid curve; line
    "assumes the tail's shape → wrong on real data → over-calls".
  - limma: icon = a per-gene empirical variance (a little box-plot per gene); line
    "estimates variance from the data → robust".
  - TOX: icon = the grey empirical-null histogram; line "estimates the *whole* noise
    distribution → most assumption-free".
- **Annotations:** bottom brace grouping limma+TOX as "empirical / robust", edgeR+DESeq2 as
  "parametric / fragile". Foreshadow: "remember this grouping — it predicts the results."
- **Colour:** use the method colours on each column header.
- **Source:** schematic.

## A8. The data is Negative Binomial (context result)

- **Purpose:** justify NB comparators and set up "Poisson too narrow".
- **Type:** mean–variance scatter with two reference curves (preferred) OR a dispersion bar.
- **Layout:** single panel.
- **Axes:** log–log; x = "gene mean count", y = "gene variance". Log–log makes the power laws
  straight lines (easier to read).
- **Elements & positions:**
  - Genes = grey dots.
  - **Poisson line** `Var = μ` (slope 1) — dashed grey.
  - **NB curve** `Var = μ + φμ²` with φ≈0.3 — solid black; it sits **above** the dots' bulk at
    high mean.
  - The dot cloud clearly tracks **above** the Poisson line.
- **Annotations:** "100% of genes reject Poisson (LRT, FDR<0.05)"; "BCV = √φ ≈ 0.5–0.65 —
  typical for human cohorts"; label the two lines.
- **Colour:** Poisson dashed grey, NB solid black, dots grey.
- **Source:** `count_distribution_poisson_vs_nb.csv` for the φ/BCV numbers; the scatter from a
  cohort's per-gene mean & variance (raw counts).
- **Alt (simpler):** a 1-row table/bar of `median_bcv` per cohort, all ≈ 0.5 — but the
  scatter teaches more.

## A9. What calibration looks like — p-value histograms (H0)

- **Purpose:** teach "uniform = good" and let the audience *see* mis-calibration.
- **Type:** small-multiples of p-value histograms, one panel per method.
- **Layout:** a single row of 5 panels (edgeR, DESeq2, limma, TOX-log, TOX-raw), shared x/y.
- **Axes:** x = "p-value" 0→1 (20 bins); y = "count". Same y-scale across panels.
- **Elements & positions:**
  - Each panel = the method's null p-value histogram at a fixed, large n.
  - A **dashed horizontal line at the uniform expectation** (flat reference) across every panel.
  - Order panels worst→best or group NB-family | empirical-family with a faint divider.
- **Annotations:** limma & TOX-log = "flat ✓"; edgeR/DESeq2 = "left spike ✗ (excess small p)";
  TOX-raw = "mild left shift — not flat". A single callout "flat under H0 = calibrated".
- **Colour:** bars in each method's colour; uniform line dashed grey.
- **Source:** `null_pvalue_histograms.png` (script already makes this) — for the talk, **subset
  to one large n** and one representative cohort so it's 5 clean panels, not a wall. Re-plot
  from `null_calibration_per_run.csv` if you want full control.

## A10. THE money figure — calibration vs sample size

- **Purpose:** the headline: NB tools **degrade with n**, limma & TOX **stay flat**.
- **Type:** line plot (one line per method).
- **Layout:** single large panel (this is the climax slide — give it the whole frame).
- **Axes:** x = "replicates per group (n_per_group: 10, 20, 40, full≈…)" — ordered, treat as
  ordinal ticks; y = **either** "FPR @ 0.05" (with a dashed line at 0.05 = target) **or**
  "fraction of 100 runs that stay calibrated" (0→1, dashed at 1.0). Pick ONE; FPR is more
  intuitive.
- **Elements & positions:**
  - 5 lines with points, method colours; markers at each n.
  - **edgeR & DESeq2 rise** away from the target as n increases (the story).
  - **limma & TOX-log stay flat** near the target; **TOX-raw flat but elevated** (~0.08).
  - Dashed horizontal reference (0.05 target, or 1.0).
  - Optional faint green band around the target (e.g. 0.03–0.07 = "well-calibrated zone").
- **Annotations:** arrow along the edgeR line "more samples → more false positives"; a bracket
  on the right grouping the two flat families; one-line title "Under a true null, a correct
  method stays at α for any n."
- **Colour:** the fixed palette; make TOX-log the boldest line.
- **Source:** aggregate `null_calibration_per_run.csv` → mean FPR (or calibrated-fraction) by
  method × n_per_group. If you use the calibrated-fraction version, define "calibrated" as e.g.
  `hits_FDR05 == 0` or `FPR ∈ [0.03,0.07]` per run and average the indicator over the 100 runs.
- **Undergrad touch:** say "flat is the whole point — nothing to detect, so nothing should
  change with n."

## A11. The honest open item — TOX-raw is not uniform yet

- **Purpose:** show you know exactly where TOX-raw stands (credibility).
- **Type:** two small panels: (left) TOX-raw p-value histogram; (right) its QQ vs uniform.
- **Layout:** side by side.
- **Axes:** left as A9 (p 0→1 vs count); right x = "expected uniform quantile", y = "observed
  p", with the **y=x diagonal** drawn.
- **Elements & positions:** left histogram shows a **left lean / excess near 0** (not flat);
  right QQ **bows above the diagonal at small expected p** (= too many small p). Optionally
  overlay TOX-log (flat / on-diagonal) in its colour as the contrast.
- **Annotations:** "raw: FPR≈0.08, p left-shifted — anti-conservative"; "hits_FDR05=0 only
  because the excess doesn't reach BH — NOT a pass"; "suspected: raw mean–variance trend →
  fix via kNN sweep".
- **Colour:** TOX-raw sky blue; TOX-log overlay blue; diagonal dashed grey.
- **Source:** `null_pvalue_histograms.png` / `null_pvalue_qqplots.png`, subset to TOX-raw
  (+ TOX-log). The scripts already emit both.

## A12. (Optional) Scorecard table — styled

- **Purpose:** a clean summary to land on.
- **Type:** a small formatted table (not a raw R dump).
- **Layout:** rows = methods (in palette colour swatches), columns = `FPR@0.05`,
  `hits_FDR05`, `stable with n?`, `assumption`.
- **Elements:** cell **background shading** green→red by how close to target (conditional
  formatting); a ✓/✗ column for "stays calibrated as n grows".
- **Annotations:** footnote "H0 = random half-split; target FPR 0.05, target hits 0".
- **Colour:** method swatches in the palette; heatmap cells green (good) → red (bad).
- **Source:** the per-method summary you already have.

---

### Figure priority (if you build only a few)
1. **A2** (core concept) — non-negotiable.
2. **A10** (calibration vs n) — the result that wins the talk.
3. **A1** (pipeline) — the roadmap.
4. **A9** (p-value histograms) — teaches "uniform = good".
5. **A7** (parametric vs empirical) — explains *why* A10 looks like it does.

Everything else is supporting. A6 and A8 are worth it if you have time, because together they
explain the TOX-raw open item mechanistically.
</content>
