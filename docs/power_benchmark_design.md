# Power benchmarking for TOX — design document

**Status:** implemented as `Simulated_data/scripts/power_test.R` (23.09.2026) — Parts P, Q,
C and R below map onto §7. Not yet run on the Tensor-Omics build; the script was
smoke-tested only against a mock of the Fortran wrappers (the R port), so no power result
exists yet. Every number quoted from existing results is cited to its source file.

**Correction (23.09.2026) to §3.2.** The original mechanism — true DE genes inflating their
neighbours' null pools — is wrong as stated for *homogeneous* effects: TOX's pools are
*within-group* residuals, and within-group centring removes a constant DE shift completely.
Contamination can only enter through DE that inflates **within-group variance**
(heterogeneous response). §3.2 and predictions P2/P3 are rewritten accordingly.

---

## 0. Summary

Calibration (Type I) is well covered: `Simulated_data/scripts/calibration_test.R` Part A
and `TCGA_test/scripts/null_calibration.R` both test it, the latter without any simulator
assumption. Power is *partially* covered — Part B already computes an FDR sweep, a TPR and
a global AUC — but as currently configured it cannot separate methods, cannot attribute a
power difference to a cause, and has no real-data arm at all.

The three things this document proposes:

1. **A metric layer** that separates *ranking quality* from *threshold behaviour*, and that
   is robust to TOX's discrete p-values (§2–§3).
2. **A truth layer** — how to define ground truth for a TPM-native method without
   accidentally destroying the null (§4). This is not a detail; it decides whether the
   benchmark measures anything.
3. **A data layer**, tiered: keep the current π-space simulator, add one signal-injection
   arm on real data, and add real-truth datasets that need no simulator at all (§5–§6).

---

## 1. What the existing Part B does and does not show

`part_b_fdr()` computes, per (distribution × n_rep × library size × round × method):
`observed_fdr`, `tpr`, `n_called`, `zero_disc` at nominal FDR ∈ {0.01, 0.05, 0.10, 0.20},
plus one global `auc`.

Three problems are visible directly in `Simulated_data/results/calibration_test_9179.out`:

**(a) The AUC does not discriminate.** On the `nb` arm, AUC is 0.7878–0.8063 for *every*
method including DESeq2, edgeR, limma and all twelve TOX arms; on `lnpois` it is
0.7928–0.8054. Across the three k-configs the AUC moves in the fourth decimal
(0.7928 / 0.7932 / 0.7932). A metric with that little spread is measuring the effect-size
distribution of the simulation, not the methods. Meanwhile the TPRs at nominal FDR 0.05 on
the same `nb` rows span 0.1199 (`TOX_bootstrap_n0_k20_50`) → 0.3283 (DESeq2) — a 2.7×
difference the AUC does not see.

**(b) TPR and observed FDR are reported, but not jointly as a curve.** The pairing matters:
on `lnpois` at nominal 0.05, `TOX_exact_n1_k20_50` reaches TPR 0.3619 at observed FDR
0.0742 while limma reaches TPR 0.3517 at observed FDR 0.0512. Those are not comparable
points — one is a method operating above its nominal level. The only fair comparison is
TPR read at a *common achieved* FDR.

**(c) Zero-discovery rounds are a resolution artefact, not conservatism, and are not
separated as such.** `TOX_bootstrap_n0_*` on `lnpois` at nominal 0.01 shows
`tpr = 0.0000`, `frac_zero_disc = 1.00` with `auc = 0.7929` — the ranking is intact and the
threshold produces nothing. §3 makes that ceiling an explicit, computable quantity instead
of an observed failure.

---

## 2. What "power" should mean here — the metric set

Split the metrics into three groups. Report all three; they answer different questions and
they can disagree.

### 2.1 Threshold behaviour (primary)

| metric | definition | why |
|---|---|---|
| **FDR–TPR curve** | for every threshold on the adjusted p-value, plot observed FDR (x) vs TPR (y); overlay the points at nominal α ∈ {0.01, 0.05, 0.10} | the standard DE-benchmark display; shows power *and* whether the nominal level is honoured, in one object |
| `TPR @ achieved FDR = 0.05` | TPR read off that curve at observed FDR 0.05 | method-fair power, independent of whether the method's own α is calibrated |
| `TPR @ nominal α`, `observed FDR @ nominal α` | as now | what a user actually gets |
| `n_called`, `frac_zero_disc` | as now | keep; instability ≠ conservatism |

`iCOBRA` (Soneson & Robinson, *Nat Methods* 2016; Bioconductor) computes exactly these:
`calculate_performance()` → `plot_fdrtprcurve()` draws the curve with the nominal-threshold
points marked, plus `plot_roc()`, `plot_fdrnbrcurve()`, `plot_tpr()`. Using it means the
plots are directly comparable to the published DE-benchmarking literature, and it removes a
whole class of hand-rolled-metric bugs. Recommended as the metric backend.

### 2.2 Ranking quality (threshold-free)

Full ROC-AUC is retained but demoted, for the reason in §1(a) and because with π₁ ≈ 10% DE
genes the ROC is dominated by the region above FPR = 0.5, which nobody operates in.
Add:

- **pAUC at FPR ≤ 0.05 and ≤ 0.10**, standardised (McClish) so it is on a 0.5–1 scale. This
  is the region BH actually reads.
- **PR-AUC / average precision.** Under class imbalance the precision–recall plot is the
  more informative of the two (Saito & Rehmsmeier, *PLoS ONE* 2015); with π₁ = 10% the
  imbalance is real.
- **Precision @ top-k** for k ∈ {100, 250, 500} — what a biologist reading the top of the
  list gets.
- **Concordance across rounds (CAT curve)** — fraction of the top-k shared between two
  independent simulated datasets. Measures stability of the ranking, which is a distinct
  failure mode from low power and is not captured by any of the above.

### 2.3 Stratified power (where the mechanism shows)

Aggregate TPR hides everything interesting about TOX, because TOX's null is built from a
mean-expression neighbourhood. Report TPR at fixed achieved FDR, broken down by:

- **|true log2FC| bin** (e.g. 0.25/0.5/1/2/4) → the power curve proper;
- **mean-expression decile** → tests whether the kNN neighbourhood costs power at the ends
  of the expression range, where the neighbourhood is necessarily one-sided in mean;
- **replicates per group** n ∈ {2, 3, 5, 10, 20} → the sample-size curve;
- **DE fraction π₁** ∈ {0.001, 0.01, 0.05, 0.10, 0.30} → see §3.2, this is TOX-specific and
  a mock null can never test it.

**Single-number summary — the detectable effect size.** For each (method, n, expression
decile), the smallest |log2FC| at which TPR ≥ 0.8 at nominal FDR 0.05. This is the
`erccdashboard` LODR idea (Munro et al., *Nat Commun* 2014) applied to the whole
transcriptome; it is the number a user of TOX would actually want, and it is directly
comparable across methods.

### 2.4 Direction and effect size

TOX's statistic is `|mean_case − mean_control|` — two-sided by construction, with no
inferential direction. So report separately:

- **sign error rate** among called genes (using the reported `lfc`), i.e. the Type-S error
  rate;
- **RMSE / MAD of the estimated log2FC against the truth**, over called genes and over all
  genes. Note this is a property of the estimator TOX reports alongside the p-value, not of
  the test.

### 2.5 Monte-Carlo error — mandatory, not optional

Every cell must carry an uncertainty. The design is naturally paired (all methods see the
same simulated dataset), so:

- report `mean ± SD/√R` over R rounds for every metric;
- for method-vs-method claims, use the **paired per-round difference** and its bootstrap CI,
  not the difference of means. Paired comparison removes the dataset-to-dataset variance,
  which is the dominant term.
- R = 5 (current `n_rounds`) is enough for a qualitative read but too few for a ranking
  claim. R ≥ 20 for anything that goes into the thesis; state R and the SE next to every
  number.

---

## 3. Two TOX-specific effects that a power test must isolate

### 3.1 The p-value floor puts a hard ceiling on TPR — compute it, don't discover it

For the exact variant the null is the enumeration of all `n_a · n_b` pairwise differences,
so the attainable p-values are `(j+1)/(N+1)`, `N = pool²`, `pool = k_max · n_rep`, and
`p_min = 1/(N+1)`.

Because *all* genes at the floor are tied, BH either rejects all of them or none: with `m`
genes at the floor, BH rejects them iff `m ≥ G · p_min / α`. So the number of discoveries
is quantised and there is a **minimum non-zero discovery set size**:

```
m*  =  ceil( G · p_min / α )        # nothing at all can be called below this
TPR_ceiling  =  m* / n_de           # the smallest non-zero TPR the method can report
```

Worked, at the two configurations in the repo:

| config | pool | p_min | G | α | m* |
|---|---|---|---|---|---|
| sim: `k_max = 50`, `n_rep = 3` | 150 | 4.4e-5 | 5 000 | 0.05 | 5 |
| TCGA (old): `K_MAX = 20`, `n_rep = 3` | 60 | 2.8e-4 | 12 000 | 0.05 | 67 |
| TCGA (current `k20_50`): `k_max = 50`, `n_rep = 3` | 150 | 4.4e-5 | 12 000 | 0.05 | 11 |

The old TCGA row says TOX at that configuration could not return a result between 1 and ~66
discoveries; `null_calibration.R` has since moved to `k20_50`, which brings that to ~11. The
pooled bootstrap engine has a fixed floor of `1/10001` regardless of `k` (m* = 24 at
G = 12 000), and the gene-blocked null, enumerated for `n_rep ≤ 5`, has
`p_min = 1/(W² + 1)` with `W = k_max · n_rep^n_rep` (`changes_summary.md` A2). Any power comparison run at that configuration is partly measuring arithmetic.
`resolution_report()` already prints this quantity; the power script must print it **next
to every TPR** and flag cells where `TPR_observed` is at or near the ceiling.

Levers, all already touched by existing work: raise `k_max` (pool grows linearly, `p_min`
quadratically), or replace the counted tail with an analytic one — the GPD / saddlepoint
options in `TCGA_test/scripts/pvalue_tail_methods.R`. A power benchmark is the right place
to quantify what that buys, because on the null it buys nothing visible.

### 3.2 Null-pool contamination by true DE genes — only through variance

TOX builds a gene's null from the **within-group** residuals of its mean-expression
neighbours: `x_ig − mean_g` in the case group and, separately, in the control group
(`prepare_sorted_data_helper`). A DE gene with a *constant* effect has its shift removed by
exactly that centring — its case residuals are ordinary noise at its case mean, and it is
placed in the case ordering by that same case mean. So with homogeneous effects there is
**no** contamination mechanism, and TOX's power should not depend on π₁.

Contamination is real only when DE also changes the **within-group variance** — the usual
case in tumours, where a subset of samples responds. Such a gene's case residuals are
inflated, they enter its neighbours' case pools, those nulls widen, and power for the
*neighbours* falls; the effect grows with π₁ and cannot appear at π₁ = 0.

Prescribed test (`power_test.R` Part Q), stated so it can fail:

- sweep π₁ ∈ {0.01, 0.05, 0.10, 0.30} × responder fraction `het` ∈ {1, 0.5};
- **prediction:** at `het = 1`, TOX's TPR at achieved FDR 0.05 is flat in π₁ within MC
  error; at `het = 0.5` it declines with π₁, relative to limma;
- **the isolating arm:** the production Fortran re-run with the pool restricted to true
  nulls (`TOX-log-oracle`, `TOX-log-blocked-oracle`). The Fortran pools every gene it is
  handed, so this is done through the input: one call with the null genes only (their
  p-values), then one call per DE gene with *null genes + that gene*, testing it alone — its
  pool is its null neighbours plus its own residuals, exactly as in production. The paired
  oracle-minus-production gap (`power_oracle_gap.csv`) is contamination and nothing else. If
  it is ~0 at `het = 0.5` too, this mechanism is not operating and the line is closed.

The oracle arm is a diagnostic that uses the truth — legitimate precisely because it is
never a candidate method.

### 3.3 Ties, when ranking

With `m` genes pinned at `p_min`, any AUC/pAUC computed on p-values alone treats them as one
undifferentiated block. `auc_from_stat()` uses `rank()`, which assigns mid-ranks — that is
the correct tie-handling for the Mann–Whitney AUC, so the number is not wrong; but it *is*
an upper bound on what a user gets from the ranked list, and it flattens exactly the top of
the list that pAUC and precision@k look at.

Fix: rank by the pair `(p, −|statistic|)` — p-value first, observed effect size as the
tie-break — and report AUC/pAUC under **both** rankings. The gap between them is a direct,
interpretable measurement of resolution loss.

---

## 4. Defining the truth for a TPM-native method

This is where a power benchmark for TOX is easiest to get silently wrong.

### 4.1 The composition shift

Already derived in this project: `log2FC_TPM = β − Δ`, with

```
Δ = log2( Σ_g π_a,g · 2^{β_g} )
```

the global composition shift (β = the true per-gene log2 ratio of absolute abundance,
π_a = relative abundances in the control group).

Consequences for benchmarking:

- **Any count-space simulator with Δ ≠ 0 makes every "null" gene non-null in TPM**, by
  exactly −Δ. Scoring TOX against the count-space truth then charges it for false positives
  that are, in TPM space, true. Scoring the count-based tools against the TPM truth does the
  mirror-image damage.
- The current `make_truth()` avoids this deliberately — it renormalises *within* the DE set
  (`cs <- sum(pi_a[de]) / sum(pi_a[de] * 2^beta[de])`) and asserts
  `max(abs(true_lfc[!is_de])) == 0`. That assertion is the guard rail. **Keep it.**

### 4.2 The rule

- **Main power comparison: run at Δ = 0.** Balanced up/down DE, renormalised so the nulls
  are exactly null on both scales. Then the two estimands coincide and the comparison is
  fair for everyone.
- **Composition stress test: a separate, explicitly labelled experiment.** Sweep the up
  fraction (0.5 → 1.0) so Δ grows, and report each method against *both* truths
  (β and β − Δ). Expected and worth stating: TOX and every TPM-native method track β − Δ;
  TMM/median-of-ratios-normalised count methods track β; neither is wrong, they estimate
  different things. The deliverable is the **breakdown point** — the Δ at which the FPR of
  each method against its *own* estimand exceeds 2×. This is the empirical counterpart of
  the identifiability result already derived (the technical noise floor is not identifiable
  from TPM alone), and it belongs in the thesis as a limitation with a number attached.

### 4.3 A shared gene universe is required

`docs/noise_model_results_summary.md` records that TOX's filter is 15–20% stricter than
edgeR/limma's (12 287–12 727 vs 15 309–16 276 genes). TPR denominators must not differ
between methods. Report both, and say which is which:

- **TPR over the common universe** (the `*-degfilt` shared gene set already implemented) —
  the method comparison;
- **TPR over all true-DE genes**, counting a filtered-out truth gene as a miss — the honest
  end-to-end number, since a gene a method never tests is a gene the user never sees.

---

## 5. How well can count/TPM data be simulated?

### 5.1 The honest answer

**Counts:** marginal properties are reproduced well. NB-based simulators fitted to real data
match the mean–variance relationship, the dispersion trend, library-size variation and the
zero fraction closely. **Joint properties are not.** Gene–gene correlation, batch and
unwanted-variation structure, and outlier patterns are where parametric simulators diverge
from real data. A published assessment of exactly the five relevant bulk simulators
(compcodeR, powsimR, SimSeq, seqgendiff, SPsimSeq; *Genes* 2022, 13:2362) evaluated Q–Q,
mean–variance, dispersion–BCV, **feature–feature correlation**, DEG counts and PCA
separability, and found no universal best — with **SimSeq and seqgendiff best preserving
feature–feature correlation**, and a caution that seqgendiff can produce artificially
enhanced class separability because signal is added by construction.

**This project already has its own evidence for the gap, and it is strong.** edgeR was
well-calibrated on the `nb` and `lnpois` synthetic arms and is 2.19–2.77× anti-conservative
on all twelve TCGA groups (`docs/noise_model_results_summary.md`). The simulation did not
reproduce whatever property of real data breaks edgeR. That is a direct, in-project
demonstration that **a synthetic-only power ranking cannot be the headline claim.**

**TPM:** no published simulator emits TPM natively — every one of them emits counts. TPM
realism is therefore inherited plus two extra assumptions:

1. the gene-length model (the current script draws `L ~ lognormal(log 2000, 0.8)`, which is
   a stand-in; using the real length distribution of the annotation is free and strictly
   better);
2. the composition constraint of §4.

There is a third thing no count simulator models at all: **quantification uncertainty**.
Real TPM comes out of an EM over multi-mapping reads with effective-length correction, and
that adds gene-specific error, especially for short genes and large paralogue families.
Only read-level simulators put the quantifier in the loop (`polyester`, Frazee et al. 2015,
Bioconductor; or BEERS2, *Brief Bioinform* 2024). Worth one arm eventually; not worth
blocking the power work on.

### 5.2 Tool table

| tool | source | mechanism | what it buys here | caveat |
|---|---|---|---|---|
| **current π-space simulator** | this repo | truth defined on π, renormalised within the DE set | the only design that guarantees **exactly-null nulls in TPM**; full control of every stress axis | marginals only; no gene–gene correlation |
| **seqgendiff** | CRAN; Gerard, *BMC Bioinformatics* 2020 | **binomial thinning** of real counts: `ỹ|y ~ Bin(y, 2^q)`, so `log2 E[Ỹ] = log2 E[Y] + q` | real correlation, real dispersion, real outliers, real batch structure **plus an exactly known log2FC** | thinning only removes reads (small depth loss, extra binomial noise); a requested fold change may not be realised for near-zero genes |
| **SPsimSeq** | Bioconductor; Assefa et al., *Bioinformatics* 2020 | log-linear density estimation of per-gene marginals from a real source dataset + **Gaussian copula** for between-gene correlation | an independent second simulator that is real-anchored and correlation-aware; handles multimodal genes | needs a real source dataset with the contrast of interest |
| **compcodeR** | Bioconductor; Soneson 2014 | NB with parameters from real data; outlier/dispersion-shift modes; built-in DE runners and an HTML comparison report | the *lingua franca* — results directly comparable to published benchmarks | parametric; no correlation structure |
| **SimSeq** | CRAN; Benidt & Nettleton, *Bioinformatics* 2015 | non-parametric resampling of real samples from a large source dataset with a real contrast | fully distribution-free, like TOX itself — good philosophical match | needs a large source dataset; truth is inherited, not designed |
| **powsimR** | GitHub (`bvieth/powsimR`); Vieth et al., *Bioinformatics* 2017 | NB parameters estimated from real data, purpose-built for power/sample-size | ready-made sample-size planning curves | GitHub-only, heavy dependency chain — an install risk, budget time for it |
| **polyester** | Bioconductor; Frazee et al. 2015 | simulates **reads** from a reference transcriptome | puts alignment + quantification in the loop; the only way to test TPM end-to-end | expensive; needs the full pipeline |

### 5.3 Recommended tiering

- **Tier 0 — keep.** The current π-space simulator, for mechanism and stress axes (§3.2,
  §4.2). Two cheap upgrades: use the real gene-length distribution instead of the lognormal
  stand-in; widen the true-LFC distribution (currently `|N(1.0, 0.6)|`) so the effect-size
  power curve has a low end to resolve — an effect-size range that never gets hard is why
  all AUCs sit at 0.79.
- **Tier 1 — add first, highest value per hour: seqgendiff thinning on real data.** Take a
  real homogeneous cohort (a TCGA stage — the same data the mock null already uses, or the
  yeast WT set), split it into two fake groups, thin one side by a known per-gene `q`. The
  result has real correlation structure *and* exact truth. It is the missing arm: the mock
  null with the signal switched on. Since the two groups come from one population, Δ is
  controllable and the nulls are genuinely null by construction, not by assumption.
- **Tier 2 — a second, independent simulator: SPsimSeq** (correlation-aware) and/or
  **compcodeR** (comparability). Two simulators disagreeing is information; one simulator
  agreeing with itself is not.
- **Tier 3 — real truth, no simulator at all.** §6. This carries the headline claim.
- **Tier 4 — optional, later: polyester/BEERS2**, if quantification uncertainty in TPM turns
  out to matter.

### 5.4 The TPM truth under thinning — closed form

Thinning gene *g* by `2^{q_g}` multiplies its expected count by `2^{q_g}`. TPM is
`π ∝ count / L` renormalised per sample, so a *sample-wide* constant thinning cancels
exactly and only the across-gene contrast in `q` survives:

```
log2FC_TPM,g  =  q_g − Δ ,        Δ = log2( Σ_h π_a,h · 2^{q_h} )
```

So: **choose `q` centred at the abundance-weighted mean** (`q'_g = q_g − Δ`) and the null
genes are exactly null in TPM, with the DE genes carrying exactly the intended
contrast — the same guarantee `make_truth()` gives, now on top of real data. To run the
composition stress test instead, deliberately leave Δ ≠ 0 and score against both truths per
§4.2. Either way the truth is known in closed form, not estimated.

---

## 6. Real-data power designs (no simulator)

These carry the headline claim, for the reason in §5.1. All four datasets are already in the
project's dataset list.

1. **Everaert et al. 2017 (GSE83402), MAQC A/B with whole-transcriptome PrimePCR RT-qPCR
   truth (~13 045 protein-coding genes).** The strongest external truth available: it is
   independent of every DE tool being compared. Define positives as `|log2FC_qPCR| ≥ 1`,
   negatives as `|log2FC_qPCR| ≤ 0.2`, and **exclude the band between** — then report the
   whole analysis again at (0.5, 0.1) to show threshold sensitivity. Yields a genuine ROC,
   pAUC, PR curve and FDR–TPR curve on real data. Caveat to state: A vs B is a huge
   contrast, so this measures the easy regime and will compress method differences.

2. **SEQC/MAQC-III (GSE47792) with ERCC spike-in ratio pools.** The ERCC design has four
   subpools of 23 transcripts at Mix1:Mix2 ratios of **4:1, 1:1, 1:1.5 and 1:2**, spanning a
   ~2²⁰ dynamic range, with the 1:1 pool as designed true negatives (Munro et al., *Nat
   Commun* 2014). This is truth by construction inside a real experiment, at *small* fold
   changes — precisely the regime the qPCR contrast misses — and it gives the LODR of §2.3
   directly. `erccdashboard` (Bioconductor) implements the ROC/AUC and LODR computation.
   Caveat: 92 spike-ins is a small truth set and they are not subject to biological
   variation, so treat this as a technical-performance floor, not a biological power
   estimate.

3. **Quartet (D5/D6/F7/M8) ratio-based reference.** Ratio truth across ~10k genes per pair
   at *small* inter-sample differences — the clinically realistic regime — with n = 3, plus
   built-in mixture truths (T1 = 3:1, T2 = 1:3) and ERCC spike-ins (Nat Commun 2024,
   s41467-024-50420-y). Two independent things to measure: DE performance against the ratio
   reference, and the **titration/monotonicity test** on the mixtures, which needs no truth
   at all — a gene DE between D5 and D6 must be intermediate in the mixtures, and violation
   rate is a real-data validity metric that no simulator can fake.

4. **Yeast 48×48 WT vs Δsnf2 (PRJEB5348) — the subsampling power curve.** Schurch et al.
   (*RNA* 2016) used exactly this design: define a gold standard from the full high-quality
   replicate set (they used 42/48 per condition), then subsample to n = 2…20 and measure
   recovery. Their headline — 8 of 11 tools recover only 20–40% of the gold-standard DE
   genes at n = 3, >85% for >4-fold genes, and >20 replicates needed for >85% at all fold
   changes — is the reference curve TOX should be plotted against. **Circularity caveat, to
   be stated, not hidden:** a gold standard defined by a tool favours that tool. Mitigate by
   defining it as the intersection of edgeR ∩ DESeq2 ∩ limma at full n with an effect-size
   cut, reporting results against each single-tool standard separately, and noting that
   this arm measures *agreement at low n with the same family's answer at high n*, not truth.

The same cohort-halving machinery already in `null_calibration.R` supplies the negative
control for all of these — π₁ = 0 on the same data, same pipeline.

---

## 7. Implementation plan

**Implemented:** `Simulated_data/scripts/power_test.R` — Part P (main grid), Part Q (π₁ ×
heterogeneity + oracle), Part C (composition stress with TPM / median-centred TPM / TMM-CPM
inputs), Part R (binomial thinning on a user-supplied real cohort). The real-truth datasets
of §6 are not in it. Differences from the plan below: the metrics are computed in-script
(tie-aware step ROC, so p-value ties are never split) rather than via iCOBRA; and the
thinning is done directly with `rbinom` with the larger-loss side scaled down so both groups
lose the same expected relative abundance, instead of `thin_diff()` + centring — the centring
formula in the sketch below does not account for the per-gene max subtraction `thin_diff()`
applies, so the sketch is superseded.

The original plan: a new script, `Simulated_data/scripts/power_test.R` (plus a real-data sibling under
`TCGA_test/scripts/`), reusing `run_all_methods()` and the method wrappers unchanged.

**Shared metric module** — one file, used by both, so the numbers mean the same thing
everywhere. Skeleton (untested; the intent is the definitions, not the code):

```r
# ---- ranking metrics -------------------------------------------------------
# A ranking score, higher = more significant.
#   key = "p"     : p-value only (what the current auc_from_stat does; mid-ranks the ties)
#   key = "p_eff" : p-value first, |effect size| as the tie-break inside a tied block
# The GAP between the two AUCs is the resolution loss from the discrete p-value grid.
rank_score <- function(res, key = c("p_eff", "p")) {
  key <- match.arg(key)
  if (key == "p") return(-res$stat)
  ord <- order(res$stat, -abs(res$lfc), na.last = TRUE)    # best first
  sc  <- numeric(length(ord)); sc[ord] <- rev(seq_along(ord))
  sc
}

# standardised partial AUC over FPR in [0, fmax]  (McClish); 0.5 = random, 1 = perfect
pauc_std <- function(score, is_de, fmax = 0.05) {
  o  <- order(score, decreasing = TRUE)
  y  <- is_de[o]
  tp <- cumsum(y) / sum(y); fp <- cumsum(!y) / sum(!y)
  k  <- fp <= fmax
  a  <- sum(diff(c(0, fp[k])) * tp[k])                    # area under ROC up to fmax
  0.5 * (1 + (a - fmax^2 / 2) / (fmax - fmax^2 / 2))      # McClish standardisation
}

pr_auc <- function(score, is_de) {                        # average precision
  o <- order(score, decreasing = TRUE); y <- is_de[o]
  prec <- cumsum(y) / seq_along(y)
  sum(prec[y]) / sum(y)
}

prec_at_k <- function(score, is_de, k) {
  o <- order(score, decreasing = TRUE)[seq_len(min(k, length(score)))]
  mean(is_de[o])
}

# ---- TPR at a COMMON ACHIEVED FDR, not at a nominal level ------------------
tpr_at_achieved_fdr <- function(padj, is_de, target = 0.05) {
  o <- order(padj, na.last = NA); y <- is_de[o]
  fdr <- cumsum(!y) / seq_along(y)
  k   <- which(fdr <= target)
  if (!length(k)) return(0)
  sum(y[seq_len(max(k))]) / sum(is_de)
}

# ---- the discreteness ceiling (§3.1) ---------------------------------------
tpr_ceiling <- function(G, n_de, pool, alpha = 0.05) {
  p_min <- 1 / (pool^2 + 1)
  m_min <- ceiling(G * p_min / alpha)     # smallest non-zero discovery set
  list(p_min = p_min, min_discoveries = m_min, tpr_floor = m_min / n_de)
}
```

Then hand the same objects to `iCOBRA` (`COBRAData(pval=, padj=, truth=)` →
`calculate_performance()` → `plot_fdrtprcurve()`) for the curves, so the primary figures are
the standard ones.

**Signal-injection harness (Tier 1)** — the shape of it:

Use `thin_diff()`, not `thin_2group()`: `thin_2group()` draws the signal itself from
`signal_fun` and picks the DE genes internally, whereas `thin_diff()` takes an explicit
`design_fixed` (N × P, no intercept) and `coef_fixed` (G × P) — which is what lets us impose
a *prespecified, abundance-centred* per-gene `q`.

```r
# one real homogeneous cohort -> two fake groups -> known signal in group B
library(seqgendiff)

grp <- rep(0:1, length.out = ncol(counts))            # fake group assignment
q   <- rep(0, G); q[de] <- sgn * abs(rnorm(n_de, lfc_mean, lfc_sd))

# relative abundances in the (unthinned) cohort, using the REAL annotation lengths L
rate <- counts / L
pi_a <- rowMeans(sweep(rate, 2, colSums(rate), "/"))

q <- q - log2(sum(pi_a * 2^q))        # centre => Delta = 0 => nulls EXACTLY null in TPM (§5.4)

th <- thin_diff(mat         = counts,
                design_fixed = matrix(grp, ncol = 1),
                coef_fixed   = matrix(q,   ncol = 1))

# th$mat = thinned counts; truth in TPM space is exactly q, by construction.
# Build TPM from th$mat with the SAME L. Verify the guard rail before scoring:
#   stopifnot(max(abs(q[-de])) < 1e-12)
```

Note the caveat from the method paper: thinning can only remove reads, so a requested
fold change may not be realised for genes with near-zero counts — check the realised
`log2` ratio of group means against `q` and drop or flag genes where the two disagree,
rather than assuming the design was achieved.

**Run grid** (the full cross is large — stage it):

| axis | values | note |
|---|---|---|
| data source | π-space sim, seqgendiff-on-TCGA, seqgendiff-on-yeast, SPsimSeq | stage 1 = first two |
| observation model (sim only) | nb, lnpois, tpois, bimodal | as now |
| n per group | 2, 3, 5, 10, 20 | the sample-size curve |
| π₁ | 0.001, 0.01, 0.05, 0.10, 0.30 | §3.2 — the new axis |
| up-fraction | 0.5 (main), 0.7, 0.9, 1.0 (stress) | §4.2 |
| true LFC | wider than now: mixture over 0.25–4 | §5.3 |
| TOX arms | exact/bootstrap × norm 0/1 × k grid, + oracle-pool arm | §3.2 |
| rounds | ≥ 20 for anything reported | §2.5 |

**Prespecified predictions** (in the style of `claude/raw_normalisation_diagnosis.md` — each
stated so it can fail):

| | prediction | falsifies what if wrong |
|---|---|---|
| P1 | TOX-log TPR at *achieved* FDR 0.05 is within ~10% of limma's on `nb`/`lnpois`; the visible gap at *nominal* levels is threshold placement, not ranking | that the current TPR gap is a calibration artefact |
| P2 | at `het = 1` TOX's TPR is flat in π₁; at `het = 0.5` it declines with π₁ relative to limma (§3.2, corrected) | that within-group centring neutralises constant effects / that variance-inflating DE contaminates pools |
| P3 | `TOX-log-oracle` − `TOX-log` ≈ 0 at `het = 1` and > 0 at `het = 0.5`, growing with π₁ | the contamination mechanism, decisively |
| P4 | TPR is markedly lower in the lowest and highest expression deciles, where the kNN neighbourhood is one-sided in mean | that the neighbourhood is expression-neutral |
| P5 | on real thinned data, edgeR's power advantage over TOX shrinks or reverses relative to the synthetic arms — mirroring its calibration reversal between `nb` and TCGA | the claim that synthetic benchmarks flatter parametric tools |
| P6 | pAUC@0.05 separates the methods by more than full AUC does (which spans < 0.03 today) | that the AUC's flatness is an imbalance artefact rather than genuine equivalence |

---

## 8. Suggested order of work

1. Metric module + `iCOBRA` wiring; re-score the **existing** Part B output. Cheap, no new
   runs, and it immediately tells you whether P1 and P6 hold.
2. Add the π₁ axis and the oracle-pool arm to the existing simulator (P2, P3). Highest
   information per unit of compute, and it is the gap a mock null structurally cannot cover.
3. seqgendiff arm on a TCGA stage already used for the mock null (P5). Same data, same
   pipeline, signal switched on.
4. Real-truth arms: Everaert/qPCR first (largest truth set), then SEQC ERCC for the small
   fold-change regime and the LODR.
5. Yeast subsampling curve, with the circularity caveat stated.
6. Only if 1–5 leave the TPM question open: composition stress test at Δ ≠ 0, and/or a
   read-level arm.

---

## 9. How to run `power_test.R` (Phase 6 of `run_order.md` — independent of Phases 3–5)

**`Simulated_data/scripts/power_test.R`** — pure simulation, no TCGA data. Needs the
    Tensor-Omics build carrying `null_method` (the `*-blocked` arms). Run from the
    Tensor-Omics root next to `calibration_test.R` and `config.R`; every TOX p-value,
    including the oracle arms, comes from the compiled Fortran.
    - Smoke test first: `POWER_GENES=1000 POWER_ROUNDS=2 POWER_CORES=4 Rscript power_test.R`.
    - Full: `POWER_PARTS=P,Q,C` (default), 20 rounds, 32 cores. The oracle arms make
      one Fortran call per DE gene (up to 1,500 at π₁ = 0.3), each a full prepare+gather
      pass; DESeq2 is likely still the slowest arm.
    - Part R (real-data thinning): `POWER_PARTS=R POWER_REAL_RDS=<cohort.rds>`, where the RDS
      is `list(counts = genes × samples, tpm = … or lengths = …)` for ONE homogeneous cohort
      (e.g. a TCGA stage: `load_counts_matrix(pid, stage)` + its raw TPM, same samples).
    → `power_out/power_{per_round,summary,paired_vs_ref,fdr_tpr_curve,stratified,fidelity}.csv` + plots.

The real-truth datasets of `power_benchmark_design.md` §6 (Everaert qPCR, SEQC ERCC,
Quartet, yeast subsampling) are still not scripted.

---

## Sources

- Soneson & Robinson, iCOBRA: open, reproducible, standardized and live method benchmarking, *Nat Methods* 2016 — https://www.nature.com/articles/nmeth.3805; `plot_fdrtprcurve` — https://rdrr.io/bioc/iCOBRA/man/plot_fdrtprcurve.html
- Gerard, Data-based RNA-seq simulations by binomial thinning, *BMC Bioinformatics* 2020 — https://link.springer.com/article/10.1186/s12859-020-3450-9; package — https://cran.r-project.org/package=seqgendiff, https://dcgerard.github.io/seqgendiff/
- Assefa et al., SPsimSeq, *Bioinformatics* 2020 — https://academic.oup.com/bioinformatics/article/36/10/3276/5739438
- compcodeR — https://bioconductor.org/packages/release/bioc/html/compcodeR.html
- Benidt & Nettleton, SimSeq, *Bioinformatics* 2015 — https://academic.oup.com/bioinformatics/article/31/13/2131/196386; https://cran.r-project.org/package=SimSeq
- Vieth et al., powsimR, *Bioinformatics* 2017 — https://academic.oup.com/bioinformatics/article/33/21/3486/3952669; https://bvieth.github.io/powsimR/
- Frazee et al., polyester — https://bioconductor.org/packages/polyester
- Assessment of synthetic RNA-seq data generators, *Genes* 2022, 13:2362 — https://www.mdpi.com/2073-4425/13/12/2362
- Munro et al., ERCC spike-in ratio mixtures / LODR, *Nat Commun* 2014 — https://www.nature.com/articles/ncomms6125; erccdashboard — https://www.bioconductor.org/packages/release/bioc/vignettes/erccdashboard/inst/doc/erccdashboard.html
- Quartet/MAQC multi-centre RNA-seq benchmark, *Nat Commun* 2024 — https://www.nature.com/articles/s41467-024-50420-y
- Schurch et al., How many biological replicates…, *RNA* 2016 — https://rnajournal.cshlp.org/content/early/2016/03/28/rna.053959.115
- Saito & Rehmsmeier, precision–recall vs ROC on imbalanced data, *PLoS ONE* 2015 — https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118432
- In-project results quoted: `Simulated_data/results/calibration_test_9179.out`, `docs/noise_model_results_summary.md`, `Simulated_data/scripts/calibration_test.R`, `TCGA_test/scripts/null_calibration.R`