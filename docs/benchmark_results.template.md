# Benchmark results: variance-aware noise model (TCGA + yeast)

*Curated summary of the power benchmark and the kNN configuration test, as of 05/11/2026. The
raw and per-round outputs are not in git; every table below is generated from the run's
summary CSVs by `common/benchmark_results_tables.R` (template: `docs/benchmark_results.template.md`),
so no number was copied by hand. Edit the template, not `benchmark_results.md`.*

## 1. Runs covered

| run | data | parts | rounds per cell | output folder (local copy) |
|---|---|---|---|---|
| TCGA power | 6 TCGA base cohorts (COAD, KIRC, LUAD × healthy, Stage IV) + simulation | RS, R, P+, U (+ N, P, Q, C, S, V, D) | RS 50, others 20 | `TCGA_test/results/power_out` (04/11/2026) |
| TCGA kNN grid | same cohorts; P+ simulated | RS, P+ | 20 | `TCGA_test/results/knn_out` (04/11/2026) |
| Yeast power | yeast 48×48 clean (`yeast48_clean.rds`): WT (42), Δsnf2 (44) | RS, R | RS 50, R 20 | `Yeast_test/results/yeast_power_out` (05/11/2026) |
| Yeast kNN grid | same two cohorts | RS | 20 | `Yeast_test/results/yeast_knn_out` (05/11/2026) |

- **Neighbourhood in the power runs:** the `xfit-*` arms use the plain TOX neighbourhood (adaptive
  k 20→50 genes, tau 0.1). Band/adaptive arms use `vmatch$k = 50`. Fixed larger k was tested
  only in the kNN grid (section 6).
- **Ellipse** (`variance_mode = 4`) is present only in the yeast kNN grid. The TCGA kNN run
  predates the rebuilt library and has no ellipse rows.
- **Yeast gene lengths:** none were available, so yeast "TPM" is CPM (all lengths 1).
- **Cohort counts:** yeast has 2 base cohorts, TCGA 6. Cross-cohort statements on yeast are
  "1 of 2" at best and are weaker evidence.

## 2. How to read the tables

- **null p≤0.001 / null p≤0.01:** share of genes below the threshold in the pi1 = 0 rounds (no
  injected effect). Nominal = 0.001 / 0.01. *Tail fold* = null p≤0.001 / 0.001.
- **FDR@0.05:** observed FDR among genes with BH-adjusted p ≤ 0.05 in the pi1 = 0.1 rounds
  (rounds without calls count as 0). *Worst cohort* = the maximum over the base cohorts.
- **TPR@nominal 0.05:** power at BH ≤ 0.05, i.e. what a user gets, but it rewards
  anti-conservative methods.
- **TPR@achieved 0.05:** power at the threshold where the *observed* FDR is 0.05. This compares
  rankings fairly across methods with different calibration.
- **AUPRC excess:** area under the precision-recall curve above the prevalence. Ranking only.
- **Real-data tables:** mean over base cohorts. Part R has no null rounds, so no null columns.
- **Arm names:** `TOX-log|TPM-centred` = plain exact model on centred log TPM (the comparator).
  `xfit-js` = cross-fitted, James-Stein-shrunk scaling. `xfit-raw` = cross-fitted, unshrunk.
  `xfit-js-nocf` = same-group (pooled) scaling, the negative control. `vband-*` = variance band.
  `adapt-*` = adaptive variance-aware growth. `|<input>` = input matrix of the Part R grid.

Plots to compare against are named under each table.

How the variance-aware model works (scaling, cross-fitting, shrinkage) is explained in
**Appendix A** at the end.

## 3. Conclusions

1. **The recommended configuration is calibrated on real data from n = 5, and at n = 3 its
   p ≤ 0.01 rate is about right.** `TOX-log-xfit-js` (mean neighbourhood, cross-fitted JS
   scaling, TPM-centred) passes the calibration, tail and FDR gates in 6/6 TCGA cohorts at
   n = 5, 7 and 10 (section 9). Its null p ≤ 0.01 rate is at or below limma/edgeR's at every n
   on TCGA, and within 0.0005 of them on yeast (tables 4.1, 4.2).
2. **n = 3 is still not solved.** The tail (p ≤ 0.001) is 3.5× nominal on TCGA and 5.7× on
   yeast with the power-run neighbourhood. FDR@0.05 is 0.167 on TCGA (worst cohort 0.248) and
   0.096 on yeast. A larger fixed neighbourhood (k = 200) roughly halves the tail excess on both
   datasets (section 6). On yeast it then matches edgeR/limma on tail and FDR. On TCGA the tail
   comes close to edgeR's (fold 1.56 vs 1.33), but FDR stays about twice edgeR's (0.129 vs 0.065).
3. **On TCGA the n = 3 problem is concentrated in the Stage IV cohorts** (xfit-js FDR 0.22–0.25,
   vs 0.09–0.11 in the healthy cohorts, table 4.4). Yeast, with clean replicates, does not show
   it to the same degree. Biological heterogeneity of the cohort is the likely driver. Not
   proven.
4. **Cross-fitting is required.** Every arm that scales with a variance containing the residual
   itself is anti-conservative at n = 3–5: `xfit-js-nocf`, `vband-js-*`, `adapt-xfit-js`
   (tables 4.3, 8.1).
5. **Shrinkage pays at small n, and it costs nothing at larger n.** xfit-js vs xfit-raw at n = 3:
   TPR@nominal 0.147 vs 0.053 (TCGA), 0.595 vs 0.308 (yeast). AUPRC excess gains ≥ 0.02 in 6/6
   TCGA cohorts. The formal "shrinkage helps" rule still says FALSE, because at n = 3 xfit-js
   fails the tail gate (section 9).
6. **The remaining gap to the reference tools is in the ranking.** At an achieved FDR of 0.05,
   xfit-js trails limma by about 5–8 TPR points on TCGA (n = 5–10) and edgeR/limma by 4–14 points
   on yeast (largest at n = 3). In Part R (large n) the gap is 2 points on both datasets.
7. **limma is the only reference tool calibrated on TCGA at n ≥ 5 and in Part R.** edgeR and
   DESeq2 are anti-conservative there (edgeR Part R FDR 0.206; RS n = 10 0.147). On yeast all three references are close to
   nominal at n ≥ 5.
8. **Neighbourhood:** tau has no effect at all (identical results with and without it on both
   datasets). Fixed k = 200 beats the adaptive 20→50 at n ≤ 5. At n = 10 the size hardly
   matters. Ellipse is clearly anti-conservative (yeast). The adaptive family loses heavily at
   large n (Part R TPR@achieved 0.424 vs 0.638 TCGA, 0.775 vs 0.930 yeast).
9. **Simulation agrees.** Under dispersion scatter (P+, sigma_d = 1.05) plain TOX-log loses most of its power.
   xfit-js restores it to limma's level at n = 10, and to about 40–60 % of limma's TPR@achieved at
   n = 3. With unequal group variances (Part U) it holds to rho = 1.5 at n ≥ 5
   (`nb`, n = 5: FDR 0.078, just above the 0.075 gate) and breaks at rho = 2, especially for `nb`
   at small n (section 8).
10. **Input normalisation:** among plain TOX-log inputs, TMM ranks best on TCGA and raw counts
   worst. On yeast, TMM and MoR are the least calibrated (FDR 0.08). The variance arm has not
   been run on the input grid yet (open item).

**Recommended configuration (proposal):** `variance_mode = mean`, `scaling = xfit_js`,
`shrink_target = median`, a fixed neighbourhood of about 200 genes (`k_start = k_max = 200`), no
tau, TPM input, log, centred. The fixed k = 200 is supported by the kNN grid only. The power
benchmark has not been run with it yet.

## 4. Part RS: real data at small n

### 4.1 TCGA (6 cohorts)
{{RS_T}}

*Plots: `power_out/rs_null_calibration.png`, `fdr_tpr_curves_vmatch_vs_ref_RS.png`, `sample_size_curve.png`.*

### 4.2 Yeast (2 cohorts)
{{RS_Y}}

*Plots: `yeast_power_out/rs_null_calibration.png`, `fdr_tpr_curves_vmatch_vs_ref_RS.png`.*

### 4.3 All variance arms at n = 3 and 5

TCGA:
{{RSVAR_T}}

Yeast:
{{RSVAR_Y}}

*Plots: `fdr_tpr_curves_vmatch_all_RS.png`, `variance_match_RS.png`.*

### 4.4 Per cohort, edgeR vs xfit-js (n = 3, 5)

TCGA:
{{RSCOH_T}}

Yeast:
{{RSCOH_Y}}

## 5. Part R: whole-cohort halves (n = 13–41 TCGA, 21–22 yeast)

### 5.1 Variance arms
TCGA:
{{R_T}}

Yeast:
{{R_Y}}

*Plots: `variance_match_R.png`, `fdr_tpr_curves_vmatch_vs_ref_R.png`.*

### 5.2 Input grid (plain TOX arms)
TCGA:
{{RGRID_T}}

Yeast:
{{RGRID_Y}}

*Plots: `grid_R_fdr_nom_05.png`, `grid_R_tpr_nom_05.png`, `grid_R_tpr_ach05.png`.*

### 5.3 Per cohort
TCGA:
{{RCOH_T}}

Yeast (WT is harder for TOX than Δsnf2; the references do not show this):
{{RCOH_Y}}

## 6. kNN configuration grid (Part RS; 20 rounds)

Configs: `k20-50` = adaptive 20→50 with tau 0.1 (the power-run default). `fixK` = fixed K genes.
`k50-200-notau` / `-t0.25` = adaptive 50→200 without tau / with tau 0.25. `band k200-c5`,
`adaptive k200-r0.5` = the best members of those families. *p floor* = median smallest
attainable p.

### 6.1 TCGA
{{KNN_T}}

### 6.2 Yeast
{{KNN_Y}}

### 6.3 Ellipse (yeast only; range over k ∈ {50, 100, 200} × r ∈ {0.5, 1})
{{ELL_Y}}

### 6.4 TCGA Part P+ (simulated): no dispersion scatter, and sigma_d = 1.05
sigma_d = 0:
{{KNNP0}}

sigma_d = 1.05:
{{KNNP1}}

*Plots: `knn_out/knn_tradeoff_{RS,Pplus}.png`, `knn_by_size_{RS,Pplus}.png`, `yeast_knn_out/knn_tradeoff_RS.png`.*

## 7. Reading of the kNN grid

- **Tau:** `k50-200-notau`, `k50-200-t0.25` and `fix200` give the same numbers in every row, on
  both datasets. Tau is not a useful parameter.
- **Size:** at n = 3 going from `k20-50` to `fix200` halves the tail fold (TCGA 3.38 → 1.56, yeast
  2.59 → 1.08) and lowers FDR (TCGA 0.175 → 0.129, yeast 0.100 → 0.076). `fix400` adds little. At
  n = 10 all sizes are within noise.
- **Why size matters at n = 3:** with k20-50 the smallest attainable p is about 4e-5. Many
  true positives tie at that floor (yeast power run: about 250 genes per round), which caps
  TPR@achieved. k = 200 lowers the floor to about 3e-6.

## 8. Simulation: Parts P+ and U (TCGA power run)

### 8.1 Part P+ (gene-wise dispersion scatter sigma_d; n = 3 and 10; sigma_d 0 and 1.05)
{{PPLUS}}

*Plots: `variance_match_Pplus.png`, `dispersion_heterogeneity.png`, `fdr_tpr_curves_vmatch_vs_ref_Pplus.png`.*

### 8.2 Part U (case-group dispersion × rho²; rho 1, 1.5, 2)
{{U}}

*Plot: `unequal_variance_rho.png`, `variance_match_U.png`.*

## 9. Decision protocol (note E v2 §6), `power_rs_gates.csv`

Gates per arm and n: calibration (null p≤0.05 ≤ comparator + 0.01), tail (null p≤0.001 ≤ 0.002
and p≤0.01 ≤ 0.015), FDR (FDR@0.05 ≤ 0.075). Each must hold in at least (cohorts − 1) cohorts:
TCGA 5 of 6, yeast 1 of 2.

| n | arm | TCGA cal / tail / FDR (of 6) | yeast cal / tail / FDR (of 2) |
|---|---|---|---|
| 3 | TOX-log\|TPM-centred | 6 / 0 / 0 | 2 / 0 / 0 |
| 3 | xfit-raw | 6 / 6 / 3 | 2 / 2 / 2 |
| 3 | xfit-js | 6 / 0 / 0 | 2 / 0 / 0 |
| 5 | xfit-raw | 6 / 6 / 6 | 2 / 1 / 2 |
| 5 | xfit-js | 6 / 6 / 6 | 2 / 1 / 1 |
| 7 | xfit-js | 6 / 6 / 6 | 2 / 1 / 2 |
| 10 | xfit-js | 6 / 6 / 6 | 2 / 1 / 2 |

- **Shrinkage rule:** at n = 3 the AUPRC excess gain (xfit-js − xfit-raw) is ≥ 0.02 in 6/6 TCGA
  cohorts (mean +0.037) and 2/2 yeast cohorts (+0.071). The rule still returns "helps = FALSE",
  because it also requires xfit-js to pass the tail gate.
- **"Not needed":** the rule says shrinkage is "not needed" from n = 5 on both datasets. This
  rule looks at AUPRC only. The TPR@nominal gain at n = 5 is still large (tables 4.1/4.2), which
  is why amending it to include TPR@nominal was proposed.

## 10. Open items and caveats

- **Power benchmark with fixed k = 200:** not run yet. Every power-run number above uses k 20→50.
- **Variance arm on the input grid** (`TOX-log-xfit-js|counts`, `|TMM`, `|MoR`, …): implemented
  05/11/2026, not run yet.
- **n = 3 on TCGA Stage IV:** FDR around 0.25 even with cross-fitting. Cause not established.
- **TCGA n = 3, xfit-js (resolved 09/11/2026, approximately):** null p ≤ 0.01 is nominal
  (0.0100), yet FDR@0.05 is 0.167. These two numbers do not contradict each other. At n = 3 on
  TCGA, BH at 0.05 makes about 240 calls out of 13 900 genes, so the BH threshold sits at
  p ≈ 0.0009, where the null rate is 3.5× nominal (p ≤ 0.001: 0.0035). The resulting estimate is
  12 500 null genes × about 0.003 ≈ 39 false calls out of 243, i.e. FDR ≈ 0.16. The same
  estimate gives about 0.10 for yeast (observed 0.096) and about 0.06 for edgeR on TCGA (0.054).
  FDR is set by the null tail *at the BH threshold*. On TCGA, low power puts that threshold deep
  in the tail. The cause of the tail excess itself is still open.
- **Yeast WT vs Δsnf2:** TOX is less calibrated on WT in both RS and R, while the references
  are not. A replicate-heterogeneity check on WT is open.
- **Unequal group variances (rho = 2):** cross-fitting assumes equal variances under H0 and
  breaks there. Leave-one-out pooled scaling was proposed as an alternative (not implemented).
- **Ellipse on TCGA** needs a rerun with the rebuilt library (`KNN_FAMILIES=ellipse`). Yeast
  already rules it out.
- **Base cohorts are fixed:** SEs reflect subsampling only. Yeast is one experiment with two
  conditions.

## Appendix A. How the variance-aware model works

This appendix describes the concepts only. The authoritative description of the current design
is `noise_model_current_state.md`; the rationale and history are in the eval doc (§21, §23).

### A.1 Starting point: the exact model

For each gene the test runs on log2(x + 1). Per group, the model takes the gene's centre (the
mean of the log values), and the observed statistic is the difference of the case and control
centres. With centring on, the median of that difference over all genes is subtracted first.

The null distribution is built from other genes:

1. **Neighbourhood:** the genes closest to the tested gene in mean expression (default:
   adaptive, 20 to 50 genes).
2. **Residual pools:** for every neighbour, each sample minus its group mean. These are corrected
   for the degree of freedom lost to the mean and divided by √n, so that they are on the scale of
   a difference of means. This gives one case pool and one control pool.
3. **Null and p-value:** every case residual is paired with every control residual. The null
   values are their differences. The p-value is the fraction of null values at least as large as
   the observed statistic, with an add-one correction. The smallest possible p-value is therefore
   1 / (case pool size × control pool size + 1). That is the *p floor* in the kNN tables.

**The limitation:** neighbours are alike in mean but not in noise. A gene that is noisier than
its neighbours gets a null that is too narrow, so it is called too often. A quieter gene gets a
null that is too wide, so it loses power. Plain TOX-log loses most of its power under dispersion
scatter for this reason (Part P+, table 8.1).

### A.2 Scaling: putting the gene and its neighbours in the same units

The variance-aware arms divide every neighbour's residuals by that neighbour's standard
deviation, and the tested gene's statistic by the tested gene's standard deviation. The pool then
holds residuals in standard-deviation units, and the statistic becomes t-like. The neighbourhood
still supplies the *shape* of the noise (tails, skew), and each gene's own variance supplies the
*width*.

This raises the question of which standard deviation to divide by.

### A.3 Cross-fitting: why each group is scaled by the other group's variance

**The problem with a gene's own variance.** At small n, a variance estimated from a group's
residuals depends strongly on those very residuals. Dividing a residual by a standard deviation
computed from the same residuals ("internal studentisation") puts a hard cap on the result. With
n samples, a residual can be at most (n − 1)/√n of its own group's standard deviation away from
the group mean, which is about 1.15 at n = 3. Large standardised residuals become impossible, so
the null has tails that are far too short. The tested statistic has no such cap, so ordinary
genes look extreme too often.

The negative-control arm `xfit-js-nocf` does exactly this, with the pooled variance of both
groups (which still contains the residual). It is several-fold anti-conservative at n = 3–5
(tables 4.3, 8.1). In the toy example it gave 25× too many null p ≤ 0.001.

**The cross-fit.** For each neighbour gene:

- its **case** residuals are divided by the standard deviation of its **control** group;
- its **control** residuals are divided by the standard deviation of its **case** group.

Under the null hypothesis both groups of a gene share one true variance. The other group's
variance is then an estimate of the same quantity that is **independent** of the residual being
scaled. There is no cap, and the standardised residuals have the heavy tails of a real
t-variable, matching the tested statistic.

**The tested gene's statistic** is divided by the standard deviation from its pooled (both
groups) variance. That is legitimate for the statistic, because the null pool comes from other
genes.

**How it is applied:**

- The scaling is applied once, to the residuals of every gene, before any neighbourhood is
  gathered.
- Neighbour selection, pairing and tail counting are unchanged. Only the numbers in the pools
  and the statistic differ.
- With `scaling = xfit_raw` the raw variances are used. With `xfit_js` the shrunk variances of
  A.4 are used.

**The assumption:** cross-fitting relies on equal case and control variances under the null.
When they differ, a group is scaled by the wrong variance. Part U shows this is tolerable up to
an SD ratio of about 1.5 and breaks at 2 (table 8.2).

### A.4 James–Stein shrinkage: stabilising the variances

**The problem.** At n = 3 a group variance rests on 2 degrees of freedom and is very noisy. Some
genes get a variance far too small by chance, and dividing by it inflates their statistic. Others
get one far too large and lose power. Cross-fitting with raw variances (`xfit-raw`) is therefore
calibrated but very conservative at small n (TPR@nominal 0.05 at n = 3: 0.053 on TCGA, 0.308 on yeast).

**The shrinkage.** Each variance is pulled part of the way toward a common target:

> shrunk variance = λ · T + (1 − λ) · raw variance

- **The target T** is the median of the pooled variances over all genes (`shrink_target =
  median`, recommended). Alternatively it is a LOESS trend of the variance against mean
  expression (`loess`). The same T serves the pooled, case and control variances.
- **The intensity λ** is one number per kind of variance (pooled, case, control), shared by all
  genes and estimated from the data. Conceptually, λ is the share of the observed spread of the
  variances around T that is pure estimation noise:

  > λ = (estimated sampling noise of the variances) / (observed spread of the variances around T)

  Both parts are summed over all genes, weighted by 1 / T² so that the comparison is relative,
  and λ is clipped to [0, 1]. The genuine gene-to-gene differences, i.e. the remainder, are kept.
- **Estimating the sampling noise without a distributional assumption.** Under the null, a
  gene's case and control variances estimate the same quantity. The squared difference between
  them therefore measures how noisy a variance estimate is: its expectation is the sum of the two
  estimates' sampling variances. Splitting it by degrees of freedom gives the noise of the case,
  control and pooled variances. No normal, negative-binomial or chi-square model is assumed; only
  moments are used.
- **λ falls with n on its own.** More replicates mean less sampling noise and a smaller λ, so the
  shrinkage fades without a switch. Measured on the benchmark data (xfit-js, power rounds):

{{LAMBDA}}

  The case and control λ are larger than the pooled one, because a single group's variance has
  only about half the degrees of freedom.

**How the shrunk variances are used** (`xfit_js`):

1. Case residuals are scaled by the shrunk **control** variance, and control residuals by the
   shrunk **case** variance (cross-fit, A.3).
2. The tested statistic is divided by the square root of the gene's shrunk **pooled** variance.
3. In the `vband-js` arms, the shrunk variance is also used to *select* neighbours of similar
   variance. Those arms do not cross-fit, and they are anti-conservative.

**What shrinkage does to the test.** With λ about 0.45 at n = 3 (table above), about half of each scale
comes from the stable global target, so the extreme small variances at 2 degrees of freedom
largely disappear. This is what turns `xfit-raw`'s
conservative test into one with roughly limma-like calibration and much higher power at n = 3–5
(tables 4.1, 4.2). Cross-fitting is not undone by the shrinkage: the shrunk control variance mixes
the control variance with global quantities (T, λ), so it stays essentially independent of the
case residual it scales.

