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
| n | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| 3 | limma | 0.0009 | 0.0091 | 0.034 | 0.048 | 0.105 | 0.178 | 0.388 |
| 3 | edgeR | 0.0013 | 0.0129 | 0.054 | 0.080 | 0.122 | 0.171 | 0.380 |
| 3 | DESeq2 | 0.0124 | 0.0330 | 0.294 | 0.437 | 0.290 | 0.152 | 0.370 |
| 3 | TOX-log\|TPM-centred | 0.0081 | 0.0213 | 0.372 | 0.542 | 0.115 | 0.014 | 0.256 |
| 3 | TOX-log-xfit-js | 0.0035 | 0.0100 | 0.167 | 0.248 | 0.147 | 0.091 | 0.350 |
| 3 | TOX-log-xfit-raw | 0.0013 | 0.0032 | 0.075 | 0.095 | 0.053 | 0.049 | 0.312 |
| 5 | limma | 0.0008 | 0.0088 | 0.036 | 0.049 | 0.263 | 0.321 | 0.485 |
| 5 | edgeR | 0.0017 | 0.0161 | 0.069 | 0.105 | 0.288 | 0.301 | 0.470 |
| 5 | DESeq2 | 0.0090 | 0.0286 | 0.215 | 0.329 | 0.386 | 0.292 | 0.468 |
| 5 | TOX-log\|TPM-centred | 0.0048 | 0.0164 | 0.197 | 0.295 | 0.141 | 0.058 | 0.348 |
| 5 | TOX-log-xfit-js | 0.0011 | 0.0059 | 0.044 | 0.064 | 0.203 | 0.246 | 0.460 |
| 5 | TOX-log-xfit-raw | 0.0005 | 0.0025 | 0.028 | 0.040 | 0.145 | 0.223 | 0.448 |
| 7 | limma | 0.0007 | 0.0077 | 0.031 | 0.053 | 0.356 | 0.411 | 0.545 |
| 7 | edgeR | 0.0024 | 0.0180 | 0.105 | 0.205 | 0.385 | 0.365 | 0.523 |
| 7 | DESeq2 | 0.0072 | 0.0263 | 0.198 | 0.322 | 0.453 | 0.366 | 0.525 |
| 7 | TOX-log\|TPM-centred | 0.0031 | 0.0131 | 0.136 | 0.172 | 0.172 | 0.114 | 0.413 |
| 7 | TOX-log-xfit-js | 0.0005 | 0.0041 | 0.024 | 0.038 | 0.245 | 0.346 | 0.524 |
| 7 | TOX-log-xfit-raw | 0.0003 | 0.0024 | 0.015 | 0.022 | 0.206 | 0.335 | 0.519 |
| 10 | limma | 0.0006 | 0.0076 | 0.036 | 0.051 | 0.444 | 0.489 | 0.600 |
| 10 | edgeR | 0.0041 | 0.0216 | 0.147 | 0.279 | 0.478 | 0.421 | 0.573 |
| 10 | DESeq2 | 0.0062 | 0.0246 | 0.183 | 0.296 | 0.517 | 0.442 | 0.581 |
| 10 | TOX-log\|TPM-centred | 0.0024 | 0.0121 | 0.107 | 0.158 | 0.223 | 0.185 | 0.478 |
| 10 | TOX-log-xfit-js | 0.0002 | 0.0035 | 0.013 | 0.018 | 0.314 | 0.439 | 0.584 |
| 10 | TOX-log-xfit-raw | 0.0002 | 0.0023 | 0.009 | 0.013 | 0.287 | 0.431 | 0.581 |

*Plots: `power_out/rs_null_calibration.png`, `fdr_tpr_curves_vmatch_vs_ref_RS.png`, `sample_size_curve.png`.*

### 4.2 Yeast (2 cohorts)
| n | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| 3 | limma | 0.0043 | 0.0145 | 0.065 | 0.070 | 0.616 | 0.680 | 0.743 |
| 3 | edgeR | 0.0042 | 0.0145 | 0.065 | 0.070 | 0.618 | 0.678 | 0.742 |
| 3 | DESeq2 | 0.0073 | 0.0176 | 0.103 | 0.109 | 0.658 | 0.667 | 0.737 |
| 3 | TOX-log\|TPM-centred | 0.0072 | 0.0184 | 0.118 | 0.119 | 0.591 | 0.484 | 0.704 |
| 3 | TOX-log-xfit-js | 0.0057 | 0.0143 | 0.096 | 0.099 | 0.595 | 0.541 | 0.710 |
| 3 | TOX-log-xfit-raw | 0.0019 | 0.0044 | 0.058 | 0.060 | 0.308 | 0.330 | 0.639 |
| 5 | limma | 0.0021 | 0.0127 | 0.048 | 0.069 | 0.743 | 0.804 | 0.804 |
| 5 | edgeR | 0.0022 | 0.0129 | 0.048 | 0.069 | 0.744 | 0.805 | 0.804 |
| 5 | DESeq2 | 0.0058 | 0.0185 | 0.079 | 0.106 | 0.781 | 0.801 | 0.802 |
| 5 | TOX-log\|TPM-centred | 0.0042 | 0.0145 | 0.088 | 0.121 | 0.666 | 0.656 | 0.771 |
| 5 | TOX-log-xfit-js | 0.0031 | 0.0126 | 0.061 | 0.082 | 0.701 | 0.729 | 0.783 |
| 5 | TOX-log-xfit-raw | 0.0014 | 0.0064 | 0.032 | 0.042 | 0.573 | 0.674 | 0.767 |
| 7 | limma | 0.0017 | 0.0110 | 0.051 | 0.068 | 0.795 | 0.838 | 0.823 |
| 7 | edgeR | 0.0018 | 0.0114 | 0.053 | 0.070 | 0.797 | 0.838 | 0.824 |
| 7 | DESeq2 | 0.0039 | 0.0156 | 0.082 | 0.093 | 0.820 | 0.836 | 0.822 |
| 7 | TOX-log\|TPM-centred | 0.0033 | 0.0131 | 0.076 | 0.080 | 0.707 | 0.718 | 0.800 |
| 7 | TOX-log-xfit-js | 0.0028 | 0.0119 | 0.061 | 0.066 | 0.751 | 0.777 | 0.808 |
| 7 | TOX-log-xfit-raw | 0.0016 | 0.0078 | 0.040 | 0.042 | 0.688 | 0.752 | 0.800 |
| 10 | limma | 0.0008 | 0.0108 | 0.045 | 0.049 | 0.855 | 0.884 | 0.843 |
| 10 | edgeR | 0.0009 | 0.0112 | 0.047 | 0.051 | 0.857 | 0.883 | 0.844 |
| 10 | DESeq2 | 0.0026 | 0.0148 | 0.071 | 0.077 | 0.870 | 0.883 | 0.843 |
| 10 | TOX-log\|TPM-centred | 0.0031 | 0.0123 | 0.067 | 0.080 | 0.775 | 0.795 | 0.825 |
| 10 | TOX-log-xfit-js | 0.0021 | 0.0111 | 0.061 | 0.071 | 0.818 | 0.844 | 0.831 |
| 10 | TOX-log-xfit-raw | 0.0014 | 0.0081 | 0.047 | 0.056 | 0.786 | 0.828 | 0.827 |

*Plots: `yeast_power_out/rs_null_calibration.png`, `fdr_tpr_curves_vmatch_vs_ref_RS.png`.*

### 4.3 All variance arms at n = 3 and 5

TCGA:
| n | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| 3 | TOX-log\|TPM-centred | 0.0081 | 0.0213 | 0.372 | 0.542 | 0.115 | 0.014 | 0.256 |
| 3 | TOX-log-xfit-js | 0.0035 | 0.0100 | 0.167 | 0.248 | 0.147 | 0.091 | 0.350 |
| 3 | TOX-log-xfit-raw | 0.0013 | 0.0032 | 0.075 | 0.095 | 0.053 | 0.049 | 0.312 |
| 3 | TOX-log-xfit-js-loess | 0.0051 | 0.0147 | 0.217 | 0.328 | 0.158 | 0.072 | 0.332 |
| 3 | TOX-log-xfit-js-nocf | 0.0217 | 0.0410 | 0.377 | 0.480 | 0.333 | 0.013 | 0.337 |
| 3 | TOX-log-vband-c3 | 0.0073 | 0.0194 | 0.265 | 0.375 | 0.186 | 0.050 | 0.317 |
| 3 | TOX-log-vband-js-c0.5 | 0.0282 | 0.0500 | 0.411 | 0.513 | 0.365 | 0.008 | 0.350 |
| 3 | TOX-log-vband-js-c1 | 0.0170 | 0.0355 | 0.326 | 0.416 | 0.315 | 0.040 | 0.352 |
| 3 | TOX-log-adapt-r0.25 | 0.0073 | 0.0194 | 0.262 | 0.373 | 0.184 | 0.063 | 0.316 |
| 3 | TOX-log-adapt-r0.5 | 0.0073 | 0.0195 | 0.263 | 0.373 | 0.185 | 0.062 | 0.316 |
| 3 | TOX-log-adapt-r1 | 0.0073 | 0.0195 | 0.263 | 0.374 | 0.185 | 0.061 | 0.316 |
| 3 | TOX-log-adapt-xfit-js | 0.0148 | 0.0332 | 0.291 | 0.376 | 0.256 | 0.055 | 0.311 |
| 5 | TOX-log\|TPM-centred | 0.0048 | 0.0164 | 0.197 | 0.295 | 0.141 | 0.058 | 0.348 |
| 5 | TOX-log-xfit-js | 0.0011 | 0.0059 | 0.044 | 0.064 | 0.203 | 0.246 | 0.460 |
| 5 | TOX-log-xfit-raw | 0.0005 | 0.0025 | 0.028 | 0.040 | 0.145 | 0.223 | 0.448 |
| 5 | TOX-log-xfit-js-loess | 0.0018 | 0.0090 | 0.067 | 0.080 | 0.216 | 0.222 | 0.449 |
| 5 | TOX-log-xfit-js-nocf | 0.0091 | 0.0237 | 0.200 | 0.243 | 0.368 | 0.139 | 0.455 |
| 5 | TOX-log-vband-c3 | 0.0032 | 0.0125 | 0.099 | 0.120 | 0.266 | 0.208 | 0.446 |
| 5 | TOX-log-vband-js-c0.5 | 0.0092 | 0.0239 | 0.202 | 0.263 | 0.376 | 0.162 | 0.461 |
| 5 | TOX-log-vband-js-c1 | 0.0056 | 0.0174 | 0.140 | 0.182 | 0.334 | 0.225 | 0.463 |
| 5 | TOX-log-adapt-r0.25 | 0.0028 | 0.0112 | 0.091 | 0.116 | 0.249 | 0.218 | 0.428 |
| 5 | TOX-log-adapt-r0.5 | 0.0028 | 0.0113 | 0.091 | 0.116 | 0.249 | 0.217 | 0.429 |
| 5 | TOX-log-adapt-r1 | 0.0028 | 0.0113 | 0.091 | 0.115 | 0.250 | 0.217 | 0.429 |
| 5 | TOX-log-adapt-xfit-js | 0.0040 | 0.0142 | 0.117 | 0.170 | 0.268 | 0.217 | 0.392 |

Yeast:
| n | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| 3 | TOX-log\|TPM-centred | 0.0072 | 0.0184 | 0.118 | 0.119 | 0.591 | 0.484 | 0.704 |
| 3 | TOX-log-xfit-js | 0.0057 | 0.0143 | 0.096 | 0.099 | 0.595 | 0.541 | 0.710 |
| 3 | TOX-log-xfit-raw | 0.0019 | 0.0044 | 0.058 | 0.060 | 0.308 | 0.330 | 0.639 |
| 3 | TOX-log-xfit-js-loess | 0.0064 | 0.0162 | 0.107 | 0.113 | 0.623 | 0.539 | 0.719 |
| 3 | TOX-log-xfit-js-nocf | 0.0203 | 0.0381 | 0.255 | 0.265 | 0.736 | 0.308 | 0.722 |
| 3 | TOX-log-vband-c3 | 0.0088 | 0.0202 | 0.137 | 0.139 | 0.635 | 0.462 | 0.716 |
| 3 | TOX-log-vband-js-c0.5 | 0.0297 | 0.0506 | 0.318 | 0.326 | 0.744 | 0.151 | 0.707 |
| 3 | TOX-log-vband-js-c1 | 0.0176 | 0.0343 | 0.226 | 0.232 | 0.719 | 0.342 | 0.721 |
| 3 | TOX-log-adapt-r0.25 | 0.0084 | 0.0196 | 0.132 | 0.132 | 0.620 | 0.507 | 0.696 |
| 3 | TOX-log-adapt-r0.5 | 0.0084 | 0.0197 | 0.132 | 0.132 | 0.628 | 0.511 | 0.704 |
| 3 | TOX-log-adapt-r1 | 0.0085 | 0.0198 | 0.132 | 0.132 | 0.632 | 0.513 | 0.708 |
| 3 | TOX-log-adapt-xfit-js | 0.0167 | 0.0350 | 0.230 | 0.236 | 0.572 | 0.304 | 0.571 |
| 5 | TOX-log\|TPM-centred | 0.0042 | 0.0145 | 0.088 | 0.121 | 0.666 | 0.656 | 0.771 |
| 5 | TOX-log-xfit-js | 0.0031 | 0.0126 | 0.061 | 0.082 | 0.701 | 0.729 | 0.783 |
| 5 | TOX-log-xfit-raw | 0.0014 | 0.0064 | 0.032 | 0.042 | 0.573 | 0.674 | 0.767 |
| 5 | TOX-log-xfit-js-loess | 0.0036 | 0.0141 | 0.067 | 0.094 | 0.724 | 0.742 | 0.790 |
| 5 | TOX-log-xfit-js-nocf | 0.0133 | 0.0293 | 0.149 | 0.179 | 0.798 | 0.626 | 0.790 |
| 5 | TOX-log-vband-c3 | 0.0051 | 0.0168 | 0.085 | 0.110 | 0.738 | 0.725 | 0.791 |
| 5 | TOX-log-vband-js-c0.5 | 0.0142 | 0.0313 | 0.158 | 0.181 | 0.791 | 0.631 | 0.782 |
| 5 | TOX-log-vband-js-c1 | 0.0083 | 0.0222 | 0.107 | 0.130 | 0.775 | 0.700 | 0.791 |
| 5 | TOX-log-adapt-r0.25 | 0.0045 | 0.0153 | 0.076 | 0.097 | 0.704 | 0.706 | 0.759 |
| 5 | TOX-log-adapt-r0.5 | 0.0046 | 0.0154 | 0.076 | 0.098 | 0.714 | 0.715 | 0.768 |
| 5 | TOX-log-adapt-r1 | 0.0046 | 0.0155 | 0.076 | 0.098 | 0.718 | 0.719 | 0.772 |
| 5 | TOX-log-adapt-xfit-js | 0.0088 | 0.0237 | 0.120 | 0.142 | 0.633 | 0.553 | 0.654 |

*Plots: `fdr_tpr_curves_vmatch_all_RS.png`, `variance_match_RS.png`.*

### 4.4 Per cohort, edgeR vs xfit-js (n = 3, 5)

TCGA:
| n | cohort | method | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | null p≤0.001 |
|---|---|---|---|---|---|---|---|
| 3 | TCGA-COAD/healthy | edgeR | 0.056 | 0.229 | 0.295 | 0.496 | 0.0011 |
| 3 | TCGA-COAD/healthy | TOX-log-xfit-js | 0.089 | 0.227 | 0.192 | 0.466 | 0.0025 |
| 3 | TCGA-COAD/Stage IV | edgeR | 0.080 | 0.042 | 0.048 | 0.279 | 0.0012 |
| 3 | TCGA-COAD/Stage IV | TOX-log-xfit-js | 0.217 | 0.094 | 0.004 | 0.264 | 0.0035 |
| 3 | TCGA-KIRC/healthy | edgeR | 0.042 | 0.184 | 0.281 | 0.468 | 0.0021 |
| 3 | TCGA-KIRC/healthy | TOX-log-xfit-js | 0.112 | 0.167 | 0.142 | 0.409 | 0.0051 |
| 3 | TCGA-KIRC/Stage IV | edgeR | 0.044 | 0.045 | 0.092 | 0.307 | 0.0014 |
| 3 | TCGA-KIRC/Stage IV | TOX-log-xfit-js | 0.248 | 0.095 | 0.018 | 0.265 | 0.0040 |
| 3 | TCGA-LUAD/healthy | edgeR | 0.042 | 0.225 | 0.273 | 0.486 | 0.0008 |
| 3 | TCGA-LUAD/healthy | TOX-log-xfit-js | 0.088 | 0.234 | 0.185 | 0.461 | 0.0024 |
| 3 | TCGA-LUAD/Stage IV | edgeR | 0.059 | 0.006 | 0.037 | 0.243 | 0.0012 |
| 3 | TCGA-LUAD/Stage IV | TOX-log-xfit-js | 0.246 | 0.067 | 0.002 | 0.233 | 0.0035 |
| 5 | TCGA-COAD/healthy | edgeR | 0.070 | 0.359 | 0.389 | 0.533 | 0.0007 |
| 5 | TCGA-COAD/healthy | TOX-log-xfit-js | 0.032 | 0.266 | 0.335 | 0.521 | 0.0008 |
| 5 | TCGA-COAD/Stage IV | edgeR | 0.090 | 0.214 | 0.195 | 0.393 | 0.0026 |
| 5 | TCGA-COAD/Stage IV | TOX-log-xfit-js | 0.050 | 0.150 | 0.167 | 0.395 | 0.0013 |
| 5 | TCGA-KIRC/healthy | edgeR | 0.037 | 0.375 | 0.429 | 0.563 | 0.0017 |
| 5 | TCGA-KIRC/healthy | TOX-log-xfit-js | 0.037 | 0.238 | 0.313 | 0.525 | 0.0012 |
| 5 | TCGA-KIRC/Stage IV | edgeR | 0.074 | 0.221 | 0.228 | 0.409 | 0.0015 |
| 5 | TCGA-KIRC/Stage IV | TOX-log-xfit-js | 0.064 | 0.129 | 0.153 | 0.386 | 0.0012 |
| 5 | TCGA-LUAD/healthy | edgeR | 0.038 | 0.388 | 0.427 | 0.578 | 0.0011 |
| 5 | TCGA-LUAD/healthy | TOX-log-xfit-js | 0.031 | 0.331 | 0.392 | 0.567 | 0.0010 |
| 5 | TCGA-LUAD/Stage IV | edgeR | 0.105 | 0.170 | 0.137 | 0.346 | 0.0025 |
| 5 | TCGA-LUAD/Stage IV | TOX-log-xfit-js | 0.049 | 0.102 | 0.115 | 0.364 | 0.0009 |

Yeast:
| n | cohort | method | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | null p≤0.001 |
|---|---|---|---|---|---|---|---|
| 3 | YEAST48/Snf2 | edgeR | 0.070 | 0.679 | 0.717 | 0.769 | 0.0018 |
| 3 | YEAST48/Snf2 | TOX-log-xfit-js | 0.094 | 0.648 | 0.589 | 0.738 | 0.0039 |
| 3 | YEAST48/WT | edgeR | 0.061 | 0.556 | 0.639 | 0.716 | 0.0065 |
| 3 | YEAST48/WT | TOX-log-xfit-js | 0.099 | 0.542 | 0.492 | 0.682 | 0.0075 |
| 5 | YEAST48/Snf2 | edgeR | 0.028 | 0.788 | 0.839 | 0.825 | 0.0008 |
| 5 | YEAST48/Snf2 | TOX-log-xfit-js | 0.041 | 0.740 | 0.792 | 0.808 | 0.0015 |
| 5 | YEAST48/WT | edgeR | 0.069 | 0.700 | 0.770 | 0.783 | 0.0035 |
| 5 | YEAST48/WT | TOX-log-xfit-js | 0.082 | 0.661 | 0.667 | 0.758 | 0.0047 |

## 5. Part R: whole-cohort halves (n = 13–41 TCGA, 21–22 yeast)

### 5.1 Variance arms
TCGA:
| method | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|
| limma | 0.036 | 0.055 | 0.633 | 0.662 | 0.712 |
| edgeR | 0.206 | 0.335 | 0.667 | 0.517 | 0.673 |
| DESeq2 | 0.163 | 0.309 | 0.674 | 0.608 | 0.693 |
| TOX-log\|TPM-centred | 0.060 | 0.071 | 0.395 | 0.399 | 0.632 |
| TOX-log-xfit-js | 0.010 | 0.031 | 0.522 | 0.638 | 0.707 |
| TOX-log-xfit-raw | 0.009 | 0.029 | 0.512 | 0.634 | 0.705 |
| TOX-log-xfit-js-loess | 0.013 | 0.034 | 0.531 | 0.635 | 0.707 |
| TOX-log-vband-c3 | 0.020 | 0.038 | 0.565 | 0.640 | 0.706 |
| TOX-log-vband-js-c0.5 | 0.026 | 0.052 | 0.579 | 0.642 | 0.706 |
| TOX-log-vband-js-c1 | 0.014 | 0.024 | 0.553 | 0.637 | 0.708 |
| TOX-log-adapt-r0.25 | 0.017 | 0.032 | 0.367 | 0.422 | 0.462 |
| TOX-log-adapt-r0.5 | 0.017 | 0.032 | 0.368 | 0.424 | 0.464 |
| TOX-log-adapt-r1 | 0.017 | 0.032 | 0.369 | 0.424 | 0.464 |
| TOX-log-adapt-xfit-js | 0.017 | 0.034 | 0.458 | 0.532 | 0.582 |

Yeast:
| method | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|
| limma | 0.032 | 0.034 | 0.932 | 0.949 | 0.877 |
| edgeR | 0.034 | 0.037 | 0.934 | 0.949 | 0.878 |
| DESeq2 | 0.044 | 0.050 | 0.937 | 0.949 | 0.878 |
| TOX-log\|TPM-centred | 0.062 | 0.093 | 0.885 | 0.899 | 0.869 |
| TOX-log-xfit-js | 0.048 | 0.068 | 0.914 | 0.930 | 0.872 |
| TOX-log-xfit-raw | 0.043 | 0.060 | 0.909 | 0.926 | 0.871 |
| TOX-log-xfit-js-loess | 0.050 | 0.071 | 0.918 | 0.933 | 0.873 |
| TOX-log-vband-c3 | 0.050 | 0.069 | 0.924 | 0.937 | 0.874 |
| TOX-log-vband-js-c0.5 | 0.062 | 0.087 | 0.922 | 0.931 | 0.869 |
| TOX-log-vband-js-c1 | 0.043 | 0.059 | 0.922 | 0.940 | 0.875 |
| TOX-log-adapt-r0.25 | 0.040 | 0.054 | 0.750 | 0.763 | 0.708 |
| TOX-log-adapt-r0.5 | 0.041 | 0.055 | 0.763 | 0.775 | 0.719 |
| TOX-log-adapt-r1 | 0.041 | 0.055 | 0.768 | 0.780 | 0.724 |
| TOX-log-adapt-xfit-js | 0.051 | 0.067 | 0.794 | 0.805 | 0.753 |

*Plots: `variance_match_R.png`, `fdr_tpr_curves_vmatch_vs_ref_R.png`.*

### 5.2 Input grid (plain TOX arms)
TCGA:
| method | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|
| limma | 0.036 | 0.055 | 0.633 | 0.662 | 0.712 |
| edgeR | 0.206 | 0.335 | 0.667 | 0.517 | 0.673 |
| DESeq2 | 0.163 | 0.309 | 0.674 | 0.608 | 0.693 |
| TOX-log\|TPM | 0.061 | 0.073 | 0.396 | 0.402 | 0.625 |
| TOX-log\|TPM-centred | 0.060 | 0.071 | 0.395 | 0.399 | 0.632 |
| TOX-log\|counts | 0.063 | 0.099 | 0.346 | 0.354 | 0.599 |
| TOX-log\|counts-centred | 0.057 | 0.088 | 0.343 | 0.359 | 0.619 |
| TOX-log\|TMM | 0.062 | 0.074 | 0.418 | 0.421 | 0.640 |
| TOX-log\|MoR | 0.073 | 0.109 | 0.377 | 0.367 | 0.621 |
| TOX-log-blocked\|TPM-centred | 0.061 | 0.070 | 0.411 | 0.403 | 0.634 |
| TOX-log-blocked\|TMM | 0.062 | 0.072 | 0.431 | 0.426 | 0.642 |
| TOX-raw\|TPM | 0.008 | 0.019 | 0.062 | 0.114 | 0.505 |
| TOX-raw\|TMM | 0.010 | 0.031 | 0.115 | 0.178 | 0.543 |

Yeast:
| method | FDR@0.05 | FDR@0.05 worst cohort | TPR@nominal 0.05 | TPR@achieved 0.05 | AUPRC excess |
|---|---|---|---|---|---|
| limma | 0.032 | 0.034 | 0.932 | 0.949 | 0.877 |
| edgeR | 0.034 | 0.037 | 0.934 | 0.949 | 0.878 |
| DESeq2 | 0.044 | 0.050 | 0.937 | 0.949 | 0.878 |
| TOX-log\|TPM | 0.068 | 0.101 | 0.884 | 0.894 | 0.864 |
| TOX-log\|TPM-centred | 0.062 | 0.093 | 0.885 | 0.899 | 0.869 |
| TOX-log\|counts | 0.069 | 0.128 | 0.787 | 0.838 | 0.824 |
| TOX-log\|counts-centred | 0.016 | 0.024 | 0.783 | 0.895 | 0.870 |
| TOX-log\|TMM | 0.079 | 0.116 | 0.904 | 0.903 | 0.869 |
| TOX-log\|MoR | 0.077 | 0.114 | 0.896 | 0.899 | 0.869 |
| TOX-log-blocked\|TPM-centred | 0.051 | 0.073 | 0.873 | 0.899 | 0.869 |
| TOX-log-blocked\|TMM | 0.064 | 0.090 | 0.889 | 0.902 | 0.870 |
| TOX-raw\|TPM | 0.076 | 0.108 | 0.891 | 0.888 | 0.864 |
| TOX-raw\|TMM | 0.075 | 0.107 | 0.901 | 0.896 | 0.869 |

*Plots: `grid_R_fdr_nom_05.png`, `grid_R_tpr_nom_05.png`, `grid_R_tpr_ach05.png`.*

### 5.3 Per cohort
TCGA:
| n | cohort | method | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess |
|---|---|---|---|---|---|---|
| 13 | TCGA-LUAD/Stage IV | limma | 0.041 | 0.366 | 0.389 | 0.527 |
| 13 | TCGA-LUAD/Stage IV | edgeR | 0.335 | 0.435 | 0.222 | 0.461 |
| 13 | TCGA-LUAD/Stage IV | TOX-log-xfit-js | 0.008 | 0.207 | 0.358 | 0.519 |
| 13 | TCGA-LUAD/Stage IV | TOX-log\|TPM-centred | 0.071 | 0.141 | 0.132 | 0.412 |
| 20 | TCGA-COAD/healthy | limma | 0.055 | 0.717 | 0.745 | 0.768 |
| 20 | TCGA-COAD/healthy | edgeR | 0.149 | 0.730 | 0.637 | 0.746 |
| 20 | TCGA-COAD/healthy | TOX-log-xfit-js | 0.010 | 0.601 | 0.718 | 0.764 |
| 20 | TCGA-COAD/healthy | TOX-log\|TPM-centred | 0.061 | 0.461 | 0.478 | 0.701 |
| 29 | TCGA-LUAD/healthy | limma | 0.032 | 0.748 | 0.782 | 0.788 |
| 29 | TCGA-LUAD/healthy | edgeR | 0.085 | 0.758 | 0.750 | 0.779 |
| 29 | TCGA-LUAD/healthy | TOX-log-xfit-js | 0.031 | 0.701 | 0.758 | 0.781 |
| 29 | TCGA-LUAD/healthy | TOX-log\|TPM-centred | 0.067 | 0.532 | 0.526 | 0.724 |
| 32 | TCGA-COAD/Stage IV | limma | 0.024 | 0.606 | 0.643 | 0.713 |
| 32 | TCGA-COAD/Stage IV | edgeR | 0.250 | 0.654 | 0.394 | 0.655 |
| 32 | TCGA-COAD/Stage IV | TOX-log-xfit-js | 0.002 | 0.502 | 0.633 | 0.711 |
| 32 | TCGA-COAD/Stage IV | TOX-log\|TPM-centred | 0.054 | 0.393 | 0.397 | 0.644 |
| 36 | TCGA-KIRC/healthy | limma | 0.029 | 0.735 | 0.761 | 0.770 |
| 36 | TCGA-KIRC/healthy | edgeR | 0.127 | 0.756 | 0.692 | 0.753 |
| 36 | TCGA-KIRC/healthy | TOX-log-xfit-js | 0.004 | 0.617 | 0.735 | 0.768 |
| 36 | TCGA-KIRC/healthy | TOX-log\|TPM-centred | 0.055 | 0.456 | 0.470 | 0.687 |
| 41 | TCGA-KIRC/Stage IV | limma | 0.035 | 0.626 | 0.652 | 0.705 |
| 41 | TCGA-KIRC/Stage IV | edgeR | 0.288 | 0.670 | 0.407 | 0.646 |
| 41 | TCGA-KIRC/Stage IV | TOX-log-xfit-js | 0.004 | 0.501 | 0.623 | 0.700 |
| 41 | TCGA-KIRC/Stage IV | TOX-log\|TPM-centred | 0.056 | 0.386 | 0.393 | 0.626 |

Yeast (WT is harder for TOX than Δsnf2; the references do not show this):
| n | cohort | method | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess |
|---|---|---|---|---|---|---|
| 21 | YEAST48/WT | limma | 0.034 | 0.913 | 0.931 | 0.868 |
| 21 | YEAST48/WT | edgeR | 0.037 | 0.914 | 0.930 | 0.868 |
| 21 | YEAST48/WT | DESeq2 | 0.050 | 0.919 | 0.930 | 0.868 |
| 21 | YEAST48/WT | TOX-log-xfit-js | 0.068 | 0.896 | 0.905 | 0.860 |
| 21 | YEAST48/WT | TOX-log\|TPM-centred | 0.093 | 0.854 | 0.855 | 0.854 |
| 21 | YEAST48/WT | TOX-log\|TMM | 0.116 | 0.877 | 0.863 | 0.855 |
| 22 | YEAST48/Snf2 | limma | 0.030 | 0.952 | 0.968 | 0.886 |
| 22 | YEAST48/Snf2 | edgeR | 0.031 | 0.954 | 0.968 | 0.887 |
| 22 | YEAST48/Snf2 | DESeq2 | 0.039 | 0.955 | 0.968 | 0.887 |
| 22 | YEAST48/Snf2 | TOX-log-xfit-js | 0.029 | 0.932 | 0.955 | 0.884 |
| 22 | YEAST48/Snf2 | TOX-log\|TPM-centred | 0.031 | 0.916 | 0.943 | 0.883 |
| 22 | YEAST48/Snf2 | TOX-log\|TMM | 0.041 | 0.930 | 0.944 | 0.884 |

## 6. kNN configuration grid (Part RS; 20 rounds)

Configs: `k20-50` = adaptive 20→50 with tau 0.1 (the power-run default). `fixK` = fixed K genes.
`k50-200-notau` / `-t0.25` = adaptive 50→200 without tau / with tau 0.25. `band k200-c5`,
`adaptive k200-r0.5` = the best members of those families. *p floor* = median smallest
attainable p.

### 6.1 TCGA
| n | family | scaling | config | null p≤0.001 | tail fold | null p≤0.01 | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | p floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | reference |  | limma | 0.0009 | 0.95 | 0.0084 | 0.029 | 0.105 | 0.177 | 0.389 | – |
| 3 | reference |  | edgeR | 0.0013 | 1.33 | 0.0122 | 0.065 | 0.122 | 0.172 | 0.381 | – |
| 3 | reference |  | DESeq2 | 0.0115 | 11.55 | 0.0308 | 0.295 | 0.293 | 0.154 | 0.371 | – |
| 3 | mean | none | k20-50 | 0.0079 | 7.86 | 0.0204 | 0.368 | 0.117 | 0.021 | 0.256 | 5.0e-05 |
| 3 | mean | none | fix200 | 0.0046 | 4.61 | 0.0171 | 0.352 | 0.056 | 0.017 | 0.250 | 2.8e-06 |
| 3 | mean | xfit_js | k20-50 | 0.0034 | 3.38 | 0.0093 | 0.175 | 0.149 | 0.088 | 0.346 | 4.5e-05 |
| 3 | mean | xfit_js | fix100 | 0.0020 | 1.98 | 0.0070 | 0.144 | 0.107 | 0.099 | 0.349 | 1.1e-05 |
| 3 | mean | xfit_js | fix200 | 0.0016 | 1.56 | 0.0066 | 0.129 | 0.090 | 0.104 | 0.351 | 2.8e-06 |
| 3 | mean | xfit_js | fix400 | 0.0013 | 1.31 | 0.0062 | 0.114 | 0.079 | 0.110 | 0.353 | 6.9e-07 |
| 3 | mean | xfit_js | k50-200-notau | 0.0016 | 1.56 | 0.0066 | 0.129 | 0.090 | 0.104 | 0.351 | 2.8e-06 |
| 3 | mean | xfit_js | k50-200-t0.25 | 0.0016 | 1.56 | 0.0066 | 0.129 | 0.090 | 0.104 | 0.351 | 2.8e-06 |
| 3 | band | xfit_js | k200-c5 | 0.0018 | 1.82 | 0.0072 | 0.128 | 0.109 | 0.105 | 0.353 | 2.8e-06 |
| 3 | adaptive | xfit_js | k200-r0.5 | 0.0042 | 4.18 | 0.0128 | 0.178 | 0.177 | 0.100 | 0.346 | 2.8e-06 |
| 3 | adaptive | none | k200-r0.5 | 0.0057 | 5.67 | 0.0170 | 0.236 | 0.174 | 0.088 | 0.324 | 2.8e-06 |
| 5 | reference |  | limma | 0.0006 | 0.56 | 0.0067 | 0.036 | 0.266 | 0.331 | 0.494 | – |
| 5 | reference |  | edgeR | 0.0012 | 1.22 | 0.0137 | 0.077 | 0.292 | 0.306 | 0.479 | – |
| 5 | reference |  | DESeq2 | 0.0076 | 7.57 | 0.0253 | 0.217 | 0.393 | 0.294 | 0.476 | – |
| 5 | mean | none | k20-50 | 0.0041 | 4.13 | 0.0146 | 0.219 | 0.146 | 0.066 | 0.350 | 1.6e-05 |
| 5 | mean | none | fix200 | 0.0027 | 2.72 | 0.0134 | 0.191 | 0.099 | 0.058 | 0.346 | 1.0e-06 |
| 5 | mean | xfit_js | k20-50 | 0.0008 | 0.82 | 0.0046 | 0.050 | 0.205 | 0.247 | 0.467 | 1.6e-05 |
| 5 | mean | xfit_js | fix100 | 0.0005 | 0.53 | 0.0039 | 0.034 | 0.171 | 0.256 | 0.468 | 4.0e-06 |
| 5 | mean | xfit_js | fix200 | 0.0004 | 0.37 | 0.0037 | 0.027 | 0.159 | 0.264 | 0.470 | 1.0e-06 |
| 5 | mean | xfit_js | fix400 | 0.0003 | 0.33 | 0.0035 | 0.023 | 0.153 | 0.269 | 0.471 | 2.5e-07 |
| 5 | mean | xfit_js | k50-200-notau | 0.0004 | 0.37 | 0.0037 | 0.027 | 0.159 | 0.264 | 0.470 | 1.0e-06 |
| 5 | mean | xfit_js | k50-200-t0.25 | 0.0004 | 0.37 | 0.0037 | 0.027 | 0.159 | 0.264 | 0.470 | 1.0e-06 |
| 5 | band | xfit_js | k200-c5 | 0.0005 | 0.52 | 0.0044 | 0.036 | 0.199 | 0.275 | 0.472 | 1.0e-06 |
| 5 | adaptive | xfit_js | k200-r0.5 | 0.0014 | 1.36 | 0.0072 | 0.068 | 0.245 | 0.258 | 0.448 | 1.0e-06 |
| 5 | adaptive | none | k200-r0.5 | 0.0017 | 1.69 | 0.0087 | 0.087 | 0.245 | 0.233 | 0.439 | 1.0e-06 |
| 10 | reference |  | limma | 0.0007 | 0.69 | 0.0079 | 0.027 | 0.430 | 0.481 | 0.594 | – |
| 10 | reference |  | edgeR | 0.0047 | 4.74 | 0.0231 | 0.129 | 0.462 | 0.410 | 0.567 | – |
| 10 | reference |  | DESeq2 | 0.0070 | 7.03 | 0.0258 | 0.158 | 0.502 | 0.438 | 0.576 | – |
| 10 | mean | none | k20-50 | 0.0026 | 2.62 | 0.0130 | 0.083 | 0.208 | 0.192 | 0.479 | 4.0e-06 |
| 10 | mean | none | fix200 | 0.0021 | 2.11 | 0.0123 | 0.067 | 0.181 | 0.185 | 0.475 | 2.5e-07 |
| 10 | mean | xfit_js | k20-50 | 0.0003 | 0.32 | 0.0038 | 0.012 | 0.300 | 0.437 | 0.581 | 4.0e-06 |
| 10 | mean | xfit_js | fix100 | 0.0002 | 0.23 | 0.0035 | 0.009 | 0.285 | 0.440 | 0.582 | 1.0e-06 |
| 10 | mean | xfit_js | fix200 | 0.0002 | 0.18 | 0.0034 | 0.007 | 0.282 | 0.442 | 0.582 | 2.5e-07 |
| 10 | mean | xfit_js | fix400 | 0.0001 | 0.15 | 0.0033 | 0.006 | 0.279 | 0.444 | 0.583 | 6.2e-08 |
| 10 | mean | xfit_js | k50-200-notau | 0.0002 | 0.18 | 0.0034 | 0.007 | 0.282 | 0.442 | 0.582 | 2.5e-07 |
| 10 | mean | xfit_js | k50-200-t0.25 | 0.0002 | 0.18 | 0.0034 | 0.007 | 0.282 | 0.442 | 0.582 | 2.5e-07 |
| 10 | band | xfit_js | k200-c5 | 0.0003 | 0.26 | 0.0039 | 0.013 | 0.328 | 0.444 | 0.581 | 2.5e-07 |
| 10 | adaptive | xfit_js | k200-r0.5 | 0.0005 | 0.54 | 0.0051 | 0.023 | 0.314 | 0.386 | 0.505 | 2.5e-07 |
| 10 | adaptive | none | k200-r0.5 | 0.0006 | 0.60 | 0.0056 | 0.022 | 0.323 | 0.390 | 0.507 | 2.5e-07 |

### 6.2 Yeast
| n | family | scaling | config | null p≤0.001 | tail fold | null p≤0.01 | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | p floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | reference |  | limma | 0.0010 | 1.03 | 0.0089 | 0.075 | 0.599 | 0.668 | 0.736 | – |
| 3 | reference |  | edgeR | 0.0009 | 0.91 | 0.0088 | 0.075 | 0.600 | 0.669 | 0.736 | – |
| 3 | reference |  | DESeq2 | 0.0030 | 3.00 | 0.0117 | 0.112 | 0.650 | 0.655 | 0.731 | – |
| 3 | mean | none | k20-50 | 0.0043 | 4.31 | 0.0145 | 0.126 | 0.574 | 0.499 | 0.698 | 4.4e-05 |
| 3 | mean | none | fix200 | 0.0022 | 2.20 | 0.0123 | 0.097 | 0.524 | 0.520 | 0.701 | 2.8e-06 |
| 3 | mean | xfit_js | k20-50 | 0.0026 | 2.59 | 0.0100 | 0.100 | 0.577 | 0.550 | 0.704 | 4.4e-05 |
| 3 | mean | xfit_js | fix100 | 0.0015 | 1.46 | 0.0082 | 0.083 | 0.543 | 0.561 | 0.708 | 1.1e-05 |
| 3 | mean | xfit_js | fix200 | 0.0011 | 1.08 | 0.0077 | 0.076 | 0.536 | 0.588 | 0.710 | 2.8e-06 |
| 3 | mean | xfit_js | fix400 | 0.0008 | 0.80 | 0.0076 | 0.074 | 0.535 | 0.590 | 0.711 | 6.9e-07 |
| 3 | mean | xfit_js | k50-200-notau | 0.0011 | 1.08 | 0.0077 | 0.076 | 0.536 | 0.588 | 0.710 | 2.8e-06 |
| 3 | mean | xfit_js | k50-200-t0.25 | 0.0011 | 1.08 | 0.0077 | 0.076 | 0.536 | 0.588 | 0.710 | 2.8e-06 |
| 3 | band | xfit_js | k200-c5 | 0.0015 | 1.45 | 0.0087 | 0.083 | 0.550 | 0.572 | 0.708 | 2.8e-06 |
| 3 | adaptive | xfit_js | k200-r0.5 | 0.0039 | 3.88 | 0.0136 | 0.127 | 0.586 | 0.536 | 0.683 | 2.8e-06 |
| 3 | adaptive | none | k200-r0.5 | 0.0035 | 3.46 | 0.0141 | 0.119 | 0.602 | 0.568 | 0.700 | 2.8e-06 |
| 5 | reference |  | limma | 0.0005 | 0.45 | 0.0067 | 0.037 | 0.742 | 0.802 | 0.803 | – |
| 5 | reference |  | edgeR | 0.0005 | 0.52 | 0.0068 | 0.037 | 0.744 | 0.802 | 0.804 | – |
| 5 | reference |  | DESeq2 | 0.0021 | 2.09 | 0.0111 | 0.072 | 0.777 | 0.802 | 0.802 | – |
| 5 | mean | none | k20-50 | 0.0025 | 2.49 | 0.0108 | 0.084 | 0.662 | 0.638 | 0.772 | 1.6e-05 |
| 5 | mean | none | fix200 | 0.0016 | 1.56 | 0.0097 | 0.075 | 0.635 | 0.660 | 0.770 | 1.0e-06 |
| 5 | mean | xfit_js | k20-50 | 0.0013 | 1.30 | 0.0073 | 0.057 | 0.691 | 0.728 | 0.782 | 1.6e-05 |
| 5 | mean | xfit_js | fix100 | 0.0010 | 0.96 | 0.0067 | 0.048 | 0.681 | 0.749 | 0.783 | 4.0e-06 |
| 5 | mean | xfit_js | fix200 | 0.0007 | 0.69 | 0.0063 | 0.045 | 0.679 | 0.753 | 0.784 | 1.0e-06 |
| 5 | mean | xfit_js | fix400 | 0.0006 | 0.61 | 0.0061 | 0.044 | 0.679 | 0.755 | 0.784 | 2.5e-07 |
| 5 | mean | xfit_js | k50-200-notau | 0.0007 | 0.69 | 0.0063 | 0.045 | 0.679 | 0.753 | 0.784 | 1.0e-06 |
| 5 | mean | xfit_js | k50-200-t0.25 | 0.0007 | 0.69 | 0.0063 | 0.045 | 0.679 | 0.753 | 0.784 | 1.0e-06 |
| 5 | band | xfit_js | k200-c5 | 0.0013 | 1.33 | 0.0080 | 0.054 | 0.692 | 0.746 | 0.781 | 1.0e-06 |
| 5 | adaptive | xfit_js | k200-r0.5 | 0.0031 | 3.10 | 0.0118 | 0.084 | 0.689 | 0.703 | 0.751 | 1.0e-06 |
| 5 | adaptive | none | k200-r0.5 | 0.0017 | 1.70 | 0.0100 | 0.067 | 0.707 | 0.744 | 0.768 | 1.0e-06 |
| 10 | reference |  | limma | 0.0016 | 1.63 | 0.0127 | 0.036 | 0.849 | 0.888 | 0.845 | – |
| 10 | reference |  | edgeR | 0.0017 | 1.70 | 0.0130 | 0.038 | 0.851 | 0.888 | 0.845 | – |
| 10 | reference |  | DESeq2 | 0.0037 | 3.75 | 0.0163 | 0.059 | 0.863 | 0.888 | 0.844 | – |
| 10 | mean | none | k20-50 | 0.0029 | 2.87 | 0.0136 | 0.061 | 0.771 | 0.804 | 0.824 | 4.0e-06 |
| 10 | mean | none | fix200 | 0.0021 | 2.06 | 0.0129 | 0.060 | 0.762 | 0.794 | 0.822 | 2.5e-07 |
| 10 | mean | xfit_js | k20-50 | 0.0021 | 2.15 | 0.0108 | 0.055 | 0.812 | 0.827 | 0.830 | 4.0e-06 |
| 10 | mean | xfit_js | fix100 | 0.0019 | 1.93 | 0.0105 | 0.051 | 0.811 | 0.845 | 0.830 | 1.0e-06 |
| 10 | mean | xfit_js | fix200 | 0.0019 | 1.94 | 0.0106 | 0.050 | 0.811 | 0.847 | 0.830 | 2.5e-07 |
| 10 | mean | xfit_js | fix400 | 0.0019 | 1.94 | 0.0106 | 0.049 | 0.811 | 0.849 | 0.831 | 6.2e-08 |
| 10 | mean | xfit_js | k50-200-notau | 0.0019 | 1.94 | 0.0106 | 0.050 | 0.811 | 0.847 | 0.830 | 2.5e-07 |
| 10 | mean | xfit_js | k50-200-t0.25 | 0.0019 | 1.94 | 0.0106 | 0.050 | 0.811 | 0.847 | 0.830 | 2.5e-07 |
| 10 | band | xfit_js | k200-c5 | 0.0024 | 2.41 | 0.0118 | 0.069 | 0.815 | 0.837 | 0.827 | 2.5e-07 |
| 10 | adaptive | xfit_js | k200-r0.5 | 0.0030 | 2.99 | 0.0125 | 0.077 | 0.757 | 0.767 | 0.760 | 2.5e-07 |
| 10 | adaptive | none | k200-r0.5 | 0.0025 | 2.55 | 0.0115 | 0.055 | 0.770 | 0.798 | 0.771 | 2.5e-07 |

### 6.3 Ellipse (yeast only; range over k ∈ {50, 100, 200} × r ∈ {0.5, 1})
| n_rep | scaling | null p≤0.01 | FDR@0.05 | TPR@achieved |
|---|---|---|---|---|
| 3 | none | 0.028–0.050 | 0.210–0.333 | 0.138–0.462 |
| 5 | none | 0.017–0.028 | 0.121–0.190 | 0.506–0.720 |
| 10 | none | 0.016–0.020 | 0.082–0.108 | 0.806–0.836 |
| 3 | xfit_js | 0.028–0.047 | 0.220–0.334 | 0.109–0.394 |
| 5 | xfit_js | 0.019–0.027 | 0.133–0.191 | 0.433–0.670 |
| 10 | xfit_js | 0.015–0.017 | 0.099–0.118 | 0.734–0.801 |

### 6.4 TCGA Part P+ (simulated): no dispersion scatter, and sigma_d = 1.05
sigma_d = 0:
| n | family | scaling | config | null p≤0.001 | tail fold | null p≤0.01 | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | p floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | reference |  | limma | 0.0008 | 0.78 | 0.0089 | 0.034 | 0.189 | 0.209 | 0.377 | – |
| 3 | reference |  | edgeR | 0.0015 | 1.55 | 0.0122 | 0.067 | 0.212 | 0.201 | 0.377 | – |
| 3 | reference |  | DESeq2 | 0.0033 | 3.30 | 0.0192 | 0.140 | 0.244 | 0.194 | 0.372 | – |
| 3 | mean | none | k20-50 | 0.0024 | 2.37 | 0.0126 | 0.103 | 0.199 | 0.128 | 0.359 | 4.4e-05 |
| 3 | mean | xfit_js | k20-50 | 0.0023 | 2.27 | 0.0121 | 0.100 | 0.201 | 0.133 | 0.360 | 4.4e-05 |
| 3 | mean | xfit_js | fix200 | 0.0014 | 1.40 | 0.0103 | 0.061 | 0.183 | 0.178 | 0.365 | 2.8e-06 |
| 3 | mean | xfit_js | fix400 | 0.0012 | 1.21 | 0.0097 | 0.054 | 0.178 | 0.179 | 0.366 | 6.9e-07 |
| 5 | reference |  | limma | 0.0008 | 0.85 | 0.0092 | 0.041 | 0.280 | 0.289 | 0.453 | – |
| 5 | reference |  | edgeR | 0.0015 | 1.46 | 0.0125 | 0.071 | 0.301 | 0.285 | 0.453 | – |
| 5 | reference |  | DESeq2 | 0.0033 | 3.30 | 0.0186 | 0.136 | 0.329 | 0.272 | 0.446 | – |
| 5 | mean | none | k20-50 | 0.0017 | 1.69 | 0.0112 | 0.079 | 0.274 | 0.251 | 0.439 | 1.6e-05 |
| 5 | mean | xfit_js | k20-50 | 0.0014 | 1.37 | 0.0105 | 0.073 | 0.276 | 0.256 | 0.441 | 1.6e-05 |
| 5 | mean | xfit_js | fix200 | 0.0010 | 1.04 | 0.0096 | 0.057 | 0.268 | 0.264 | 0.444 | 1.0e-06 |
| 5 | mean | xfit_js | fix400 | 0.0009 | 0.94 | 0.0094 | 0.051 | 0.266 | 0.266 | 0.444 | 2.5e-07 |
| 10 | reference |  | limma | 0.0010 | 0.99 | 0.0098 | 0.040 | 0.414 | 0.424 | 0.559 | – |
| 10 | reference |  | edgeR | 0.0016 | 1.56 | 0.0128 | 0.068 | 0.434 | 0.420 | 0.561 | – |
| 10 | reference |  | DESeq2 | 0.0028 | 2.81 | 0.0174 | 0.100 | 0.448 | 0.411 | 0.555 | – |
| 10 | mean | none | k20-50 | 0.0016 | 1.55 | 0.0112 | 0.052 | 0.402 | 0.400 | 0.550 | 4.0e-06 |
| 10 | mean | xfit_js | k20-50 | 0.0012 | 1.21 | 0.0101 | 0.043 | 0.405 | 0.414 | 0.553 | 4.0e-06 |
| 10 | mean | xfit_js | fix200 | 0.0010 | 0.98 | 0.0098 | 0.040 | 0.401 | 0.418 | 0.554 | 2.5e-07 |
| 10 | mean | xfit_js | fix400 | 0.0010 | 0.95 | 0.0096 | 0.038 | 0.401 | 0.418 | 0.554 | 6.2e-08 |

sigma_d = 1.05:
| n | family | scaling | config | null p≤0.001 | tail fold | null p≤0.01 | FDR@0.05 | TPR@nominal | TPR@achieved | AUPRC excess | p floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | reference |  | limma | 0.0008 | 0.80 | 0.0095 | 0.034 | 0.132 | 0.156 | 0.388 | – |
| 3 | reference |  | edgeR | 0.0011 | 1.06 | 0.0112 | 0.045 | 0.146 | 0.157 | 0.390 | – |
| 3 | reference |  | DESeq2 | 0.0094 | 9.43 | 0.0277 | 0.275 | 0.303 | 0.159 | 0.390 | – |
| 3 | mean | none | k20-50 | 0.0060 | 6.00 | 0.0194 | 0.263 | 0.162 | 0.009 | 0.300 | 4.6e-05 |
| 3 | mean | xfit_js | k20-50 | 0.0032 | 3.16 | 0.0119 | 0.131 | 0.198 | 0.095 | 0.375 | 4.4e-05 |
| 3 | mean | xfit_js | fix200 | 0.0012 | 1.22 | 0.0077 | 0.062 | 0.145 | 0.135 | 0.381 | 2.8e-06 |
| 3 | mean | xfit_js | fix400 | 0.0010 | 0.97 | 0.0073 | 0.049 | 0.132 | 0.141 | 0.383 | 6.9e-07 |
| 5 | reference |  | limma | 0.0010 | 1.02 | 0.0095 | 0.042 | 0.295 | 0.308 | 0.491 | – |
| 5 | reference |  | edgeR | 0.0013 | 1.31 | 0.0110 | 0.051 | 0.307 | 0.310 | 0.494 | – |
| 5 | reference |  | DESeq2 | 0.0063 | 6.25 | 0.0232 | 0.186 | 0.404 | 0.310 | 0.493 | – |
| 5 | mean | none | k20-50 | 0.0040 | 3.97 | 0.0161 | 0.166 | 0.206 | 0.088 | 0.395 | 1.6e-05 |
| 5 | mean | xfit_js | k20-50 | 0.0017 | 1.66 | 0.0087 | 0.059 | 0.291 | 0.284 | 0.485 | 1.6e-05 |
| 5 | mean | xfit_js | fix200 | 0.0008 | 0.83 | 0.0072 | 0.037 | 0.267 | 0.294 | 0.488 | 1.0e-06 |
| 5 | mean | xfit_js | fix400 | 0.0007 | 0.72 | 0.0069 | 0.031 | 0.264 | 0.301 | 0.489 | 2.5e-07 |
| 10 | reference |  | limma | 0.0009 | 0.94 | 0.0097 | 0.040 | 0.462 | 0.475 | 0.603 | – |
| 10 | reference |  | edgeR | 0.0013 | 1.34 | 0.0119 | 0.053 | 0.479 | 0.475 | 0.605 | – |
| 10 | reference |  | DESeq2 | 0.0037 | 3.67 | 0.0179 | 0.116 | 0.524 | 0.475 | 0.604 | – |
| 10 | mean | none | k20-50 | 0.0030 | 3.02 | 0.0145 | 0.103 | 0.299 | 0.236 | 0.512 | 4.0e-06 |
| 10 | mean | xfit_js | k20-50 | 0.0010 | 1.00 | 0.0085 | 0.040 | 0.444 | 0.461 | 0.600 | 4.0e-06 |
| 10 | mean | xfit_js | fix200 | 0.0007 | 0.73 | 0.0081 | 0.030 | 0.439 | 0.468 | 0.601 | 2.5e-07 |
| 10 | mean | xfit_js | fix400 | 0.0007 | 0.67 | 0.0079 | 0.028 | 0.439 | 0.469 | 0.601 | 6.2e-08 |

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
| dist | n_rep | sigma_d | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | TPR@achieved | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| lnpois | 3 | 0 | limma | 0.0008 | 0.0093 | 0.035 | 0.205 | 0.370 |
| lnpois | 3 | 0 | edgeR | 0.0019 | 0.0146 | 0.090 | 0.187 | 0.358 |
| lnpois | 3 | 0 | TOX-log | 0.0022 | 0.0130 | 0.111 | 0.134 | 0.352 |
| lnpois | 3 | 0 | TOX-log-xfit-js | 0.0020 | 0.0124 | 0.106 | 0.133 | 0.354 |
| lnpois | 3 | 0 | TOX-log-xfit-raw | 0.0018 | 0.0044 | 0.147 | 0.009 | 0.238 |
| lnpois | 3 | 0 | TOX-log-xfit-js-nocf | 0.0072 | 0.0249 | 0.255 | 0.000 | 0.355 |
| lnpois | 3 | 1.05 | limma | 0.0008 | 0.0091 | 0.031 | 0.151 | 0.384 |
| lnpois | 3 | 1.05 | edgeR | 0.0011 | 0.0113 | 0.041 | 0.153 | 0.381 |
| lnpois | 3 | 1.05 | TOX-log | 0.0052 | 0.0185 | 0.240 | 0.000 | 0.306 |
| lnpois | 3 | 1.05 | TOX-log-xfit-js | 0.0032 | 0.0116 | 0.132 | 0.091 | 0.372 |
| lnpois | 3 | 1.05 | TOX-log-xfit-raw | 0.0018 | 0.0046 | 0.113 | 0.040 | 0.317 |
| lnpois | 3 | 1.05 | TOX-log-xfit-js-nocf | 0.0223 | 0.0447 | 0.418 | 0.000 | 0.362 |
| lnpois | 10 | 0 | limma | 0.0008 | 0.0095 | 0.044 | 0.447 | 0.567 |
| lnpois | 10 | 0 | edgeR | 0.0020 | 0.0153 | 0.088 | 0.422 | 0.554 |
| lnpois | 10 | 0 | TOX-log | 0.0014 | 0.0103 | 0.056 | 0.414 | 0.557 |
| lnpois | 10 | 0 | TOX-log-xfit-js | 0.0013 | 0.0098 | 0.052 | 0.422 | 0.560 |
| lnpois | 10 | 0 | TOX-log-xfit-raw | 0.0006 | 0.0065 | 0.025 | 0.390 | 0.547 |
| lnpois | 10 | 0 | TOX-log-xfit-js-nocf | 0.0029 | 0.0152 | 0.105 | 0.426 | 0.561 |
| lnpois | 10 | 1.05 | limma | 0.0011 | 0.0095 | 0.039 | 0.502 | 0.620 |
| lnpois | 10 | 1.05 | edgeR | 0.0016 | 0.0132 | 0.061 | 0.492 | 0.615 |
| lnpois | 10 | 1.05 | TOX-log | 0.0027 | 0.0135 | 0.092 | 0.307 | 0.538 |
| lnpois | 10 | 1.05 | TOX-log-xfit-js | 0.0010 | 0.0082 | 0.036 | 0.489 | 0.615 |
| lnpois | 10 | 1.05 | TOX-log-xfit-raw | 0.0007 | 0.0061 | 0.027 | 0.486 | 0.610 |
| lnpois | 10 | 1.05 | TOX-log-xfit-js-nocf | 0.0045 | 0.0184 | 0.128 | 0.494 | 0.615 |
| nb | 3 | 0 | limma | 0.0008 | 0.0089 | 0.040 | 0.167 | 0.329 |
| nb | 3 | 0 | edgeR | 0.0010 | 0.0104 | 0.051 | 0.180 | 0.339 |
| nb | 3 | 0 | TOX-log | 0.0024 | 0.0124 | 0.114 | 0.104 | 0.314 |
| nb | 3 | 0 | TOX-log-xfit-js | 0.0024 | 0.0120 | 0.111 | 0.105 | 0.316 |
| nb | 3 | 0 | TOX-log-xfit-raw | 0.0017 | 0.0047 | 0.118 | 0.011 | 0.216 |
| nb | 3 | 0 | TOX-log-xfit-js-nocf | 0.0084 | 0.0257 | 0.297 | 0.000 | 0.318 |
| nb | 3 | 1.05 | limma | 0.0007 | 0.0097 | 0.052 | 0.150 | 0.378 |
| nb | 3 | 1.05 | edgeR | 0.0009 | 0.0111 | 0.066 | 0.149 | 0.384 |
| nb | 3 | 1.05 | TOX-log | 0.0065 | 0.0199 | 0.299 | 0.000 | 0.281 |
| nb | 3 | 1.05 | TOX-log-xfit-js | 0.0033 | 0.0116 | 0.146 | 0.061 | 0.359 |
| nb | 3 | 1.05 | TOX-log-xfit-raw | 0.0020 | 0.0052 | 0.140 | 0.030 | 0.308 |
| nb | 3 | 1.05 | TOX-log-xfit-js-nocf | 0.0259 | 0.0480 | 0.460 | 0.000 | 0.348 |
| nb | 10 | 0 | limma | 0.0008 | 0.0090 | 0.044 | 0.396 | 0.539 |
| nb | 10 | 0 | edgeR | 0.0013 | 0.0103 | 0.054 | 0.411 | 0.552 |
| nb | 10 | 0 | TOX-log | 0.0015 | 0.0104 | 0.056 | 0.370 | 0.531 |
| nb | 10 | 0 | TOX-log-xfit-js | 0.0013 | 0.0102 | 0.054 | 0.378 | 0.534 |
| nb | 10 | 0 | TOX-log-xfit-raw | 0.0006 | 0.0065 | 0.032 | 0.340 | 0.521 |
| nb | 10 | 0 | TOX-log-xfit-js-nocf | 0.0033 | 0.0159 | 0.110 | 0.375 | 0.535 |
| nb | 10 | 1.05 | limma | 0.0009 | 0.0095 | 0.041 | 0.444 | 0.578 |
| nb | 10 | 1.05 | edgeR | 0.0010 | 0.0100 | 0.047 | 0.458 | 0.589 |
| nb | 10 | 1.05 | TOX-log | 0.0032 | 0.0143 | 0.116 | 0.190 | 0.477 |
| nb | 10 | 1.05 | TOX-log-xfit-js | 0.0008 | 0.0081 | 0.034 | 0.441 | 0.576 |
| nb | 10 | 1.05 | TOX-log-xfit-raw | 0.0007 | 0.0063 | 0.024 | 0.438 | 0.570 |
| nb | 10 | 1.05 | TOX-log-xfit-js-nocf | 0.0047 | 0.0198 | 0.138 | 0.447 | 0.576 |

*Plots: `variance_match_Pplus.png`, `dispersion_heterogeneity.png`, `fdr_tpr_curves_vmatch_vs_ref_Pplus.png`.*

### 8.2 Part U (case-group dispersion × rho²; rho 1, 1.5, 2)
| dist | n_rep | rho | method | null p≤0.001 | null p≤0.01 | FDR@0.05 | TPR@achieved | AUPRC excess |
|---|---|---|---|---|---|---|---|---|
| lnpois | 3 | 1 | limma | 0.0008 | 0.0096 | 0.034 | 0.176 | 0.407 |
| lnpois | 3 | 1 | edgeR | 0.0010 | 0.0114 | 0.041 | 0.182 | 0.406 |
| lnpois | 3 | 1 | TOX-log | 0.0059 | 0.0189 | 0.231 | 0.023 | 0.336 |
| lnpois | 3 | 1 | TOX-log-xfit-js | 0.0032 | 0.0118 | 0.139 | 0.115 | 0.393 |
| lnpois | 3 | 1 | TOX-log-xfit-raw | 0.0016 | 0.0044 | 0.131 | 0.038 | 0.336 |
| lnpois | 3 | 1.5 | limma | 0.0009 | 0.0106 | 0.042 | 0.112 | 0.345 |
| lnpois | 3 | 1.5 | edgeR | 0.0010 | 0.0128 | 0.048 | 0.105 | 0.340 |
| lnpois | 3 | 1.5 | TOX-log | 0.0058 | 0.0196 | 0.275 | 0.008 | 0.277 |
| lnpois | 3 | 1.5 | TOX-log-xfit-js | 0.0029 | 0.0110 | 0.145 | 0.053 | 0.325 |
| lnpois | 3 | 1.5 | TOX-log-xfit-raw | 0.0019 | 0.0047 | 0.121 | 0.028 | 0.268 |
| lnpois | 3 | 2 | limma | 0.0007 | 0.0109 | 0.025 | 0.074 | 0.282 |
| lnpois | 3 | 2 | edgeR | 0.0011 | 0.0148 | 0.040 | 0.074 | 0.272 |
| lnpois | 3 | 2 | TOX-log | 0.0066 | 0.0219 | 0.308 | 0.000 | 0.222 |
| lnpois | 3 | 2 | TOX-log-xfit-js | 0.0045 | 0.0152 | 0.228 | 0.005 | 0.252 |
| lnpois | 3 | 2 | TOX-log-xfit-raw | 0.0020 | 0.0042 | 0.180 | 0.011 | 0.200 |
| lnpois | 5 | 1 | limma | 0.0010 | 0.0097 | 0.045 | 0.317 | 0.514 |
| lnpois | 5 | 1 | edgeR | 0.0012 | 0.0120 | 0.051 | 0.317 | 0.510 |
| lnpois | 5 | 1 | TOX-log | 0.0037 | 0.0151 | 0.141 | 0.147 | 0.429 |
| lnpois | 5 | 1 | TOX-log-xfit-js | 0.0017 | 0.0090 | 0.056 | 0.301 | 0.508 |
| lnpois | 5 | 1 | TOX-log-xfit-raw | 0.0008 | 0.0049 | 0.031 | 0.267 | 0.490 |
| lnpois | 5 | 1.5 | limma | 0.0010 | 0.0103 | 0.043 | 0.249 | 0.446 |
| lnpois | 5 | 1.5 | edgeR | 0.0014 | 0.0138 | 0.052 | 0.244 | 0.436 |
| lnpois | 5 | 1.5 | TOX-log | 0.0037 | 0.0158 | 0.165 | 0.082 | 0.363 |
| lnpois | 5 | 1.5 | TOX-log-xfit-js | 0.0011 | 0.0062 | 0.052 | 0.214 | 0.422 |
| lnpois | 5 | 1.5 | TOX-log-xfit-raw | 0.0009 | 0.0039 | 0.039 | 0.176 | 0.402 |
| lnpois | 5 | 2 | limma | 0.0011 | 0.0114 | 0.060 | 0.177 | 0.371 |
| lnpois | 5 | 2 | edgeR | 0.0016 | 0.0163 | 0.078 | 0.167 | 0.355 |
| lnpois | 5 | 2 | TOX-log | 0.0045 | 0.0183 | 0.209 | 0.048 | 0.296 |
| lnpois | 5 | 2 | TOX-log-xfit-js | 0.0023 | 0.0111 | 0.106 | 0.119 | 0.332 |
| lnpois | 5 | 2 | TOX-log-xfit-raw | 0.0009 | 0.0035 | 0.052 | 0.105 | 0.303 |
| lnpois | 10 | 1 | limma | 0.0010 | 0.0094 | 0.046 | 0.499 | 0.623 |
| lnpois | 10 | 1 | edgeR | 0.0015 | 0.0130 | 0.069 | 0.487 | 0.618 |
| lnpois | 10 | 1 | TOX-log | 0.0025 | 0.0137 | 0.090 | 0.293 | 0.550 |
| lnpois | 10 | 1 | TOX-log-xfit-js | 0.0009 | 0.0083 | 0.044 | 0.483 | 0.619 |
| lnpois | 10 | 1 | TOX-log-xfit-raw | 0.0007 | 0.0059 | 0.033 | 0.473 | 0.611 |
| lnpois | 10 | 1.5 | limma | 0.0011 | 0.0105 | 0.045 | 0.416 | 0.560 |
| lnpois | 10 | 1.5 | edgeR | 0.0020 | 0.0177 | 0.085 | 0.396 | 0.547 |
| lnpois | 10 | 1.5 | TOX-log | 0.0026 | 0.0153 | 0.111 | 0.222 | 0.470 |
| lnpois | 10 | 1.5 | TOX-log-xfit-js | 0.0004 | 0.0051 | 0.017 | 0.381 | 0.530 |
| lnpois | 10 | 1.5 | TOX-log-xfit-raw | 0.0008 | 0.0058 | 0.019 | 0.371 | 0.517 |
| lnpois | 10 | 2 | limma | 0.0011 | 0.0111 | 0.056 | 0.348 | 0.507 |
| lnpois | 10 | 2 | edgeR | 0.0032 | 0.0262 | 0.145 | 0.311 | 0.478 |
| lnpois | 10 | 2 | TOX-log | 0.0028 | 0.0171 | 0.122 | 0.189 | 0.402 |
| lnpois | 10 | 2 | TOX-log-xfit-js | 0.0008 | 0.0086 | 0.042 | 0.275 | 0.438 |
| lnpois | 10 | 2 | TOX-log-xfit-raw | 0.0006 | 0.0045 | 0.028 | 0.249 | 0.410 |
| nb | 3 | 1 | limma | 0.0009 | 0.0099 | 0.034 | 0.155 | 0.387 |
| nb | 3 | 1 | edgeR | 0.0012 | 0.0109 | 0.042 | 0.162 | 0.394 |
| nb | 3 | 1 | TOX-log | 0.0066 | 0.0197 | 0.304 | 0.000 | 0.286 |
| nb | 3 | 1 | TOX-log-xfit-js | 0.0033 | 0.0121 | 0.129 | 0.063 | 0.374 |
| nb | 3 | 1 | TOX-log-xfit-raw | 0.0021 | 0.0049 | 0.124 | 0.038 | 0.320 |
| nb | 3 | 1.5 | limma | 0.0009 | 0.0105 | 0.043 | 0.072 | 0.302 |
| nb | 3 | 1.5 | edgeR | 0.0014 | 0.0124 | 0.069 | 0.057 | 0.304 |
| nb | 3 | 1.5 | TOX-log | 0.0095 | 0.0254 | 0.435 | 0.000 | 0.195 |
| nb | 3 | 1.5 | TOX-log-xfit-js | 0.0041 | 0.0109 | 0.218 | 0.007 | 0.268 |
| nb | 3 | 1.5 | TOX-log-xfit-raw | 0.0020 | 0.0047 | 0.182 | 0.015 | 0.246 |
| nb | 3 | 2 | limma | 0.0013 | 0.0132 | 0.019 | 0.030 | 0.234 |
| nb | 3 | 2 | edgeR | 0.0022 | 0.0155 | 0.198 | 0.009 | 0.224 |
| nb | 3 | 2 | TOX-log | 0.0193 | 0.0381 | 0.643 | 0.000 | 0.124 |
| nb | 3 | 2 | TOX-log-xfit-js | 0.0119 | 0.0227 | 0.510 | 0.000 | 0.164 |
| nb | 3 | 2 | TOX-log-xfit-raw | 0.0029 | 0.0052 | 0.360 | 0.004 | 0.164 |
| nb | 5 | 1 | limma | 0.0010 | 0.0103 | 0.040 | 0.310 | 0.494 |
| nb | 5 | 1 | edgeR | 0.0012 | 0.0113 | 0.051 | 0.307 | 0.503 |
| nb | 5 | 1 | TOX-log | 0.0043 | 0.0169 | 0.192 | 0.026 | 0.382 |
| nb | 5 | 1 | TOX-log-xfit-js | 0.0016 | 0.0091 | 0.051 | 0.284 | 0.488 |
| nb | 5 | 1 | TOX-log-xfit-raw | 0.0010 | 0.0051 | 0.037 | 0.251 | 0.473 |
| nb | 5 | 1.5 | limma | 0.0011 | 0.0116 | 0.053 | 0.210 | 0.403 |
| nb | 5 | 1.5 | edgeR | 0.0015 | 0.0126 | 0.066 | 0.195 | 0.410 |
| nb | 5 | 1.5 | TOX-log | 0.0080 | 0.0238 | 0.350 | 0.000 | 0.258 |
| nb | 5 | 1.5 | TOX-log-xfit-js | 0.0016 | 0.0066 | 0.078 | 0.118 | 0.367 |
| nb | 5 | 1.5 | TOX-log-xfit-raw | 0.0010 | 0.0042 | 0.050 | 0.142 | 0.375 |
| nb | 5 | 2 | limma | 0.0021 | 0.0165 | 0.091 | 0.126 | 0.322 |
| nb | 5 | 2 | edgeR | 0.0025 | 0.0156 | 0.121 | 0.099 | 0.324 |
| nb | 5 | 2 | TOX-log | 0.0174 | 0.0390 | 0.572 | 0.000 | 0.166 |
| nb | 5 | 2 | TOX-log-xfit-js | 0.0072 | 0.0185 | 0.348 | 0.000 | 0.238 |
| nb | 5 | 2 | TOX-log-xfit-raw | 0.0016 | 0.0046 | 0.145 | 0.034 | 0.281 |
| nb | 10 | 1 | limma | 0.0010 | 0.0098 | 0.038 | 0.478 | 0.605 |
| nb | 10 | 1 | edgeR | 0.0012 | 0.0106 | 0.049 | 0.489 | 0.615 |
| nb | 10 | 1 | TOX-log | 0.0032 | 0.0153 | 0.105 | 0.222 | 0.505 |
| nb | 10 | 1 | TOX-log-xfit-js | 0.0009 | 0.0085 | 0.036 | 0.472 | 0.603 |
| nb | 10 | 1 | TOX-log-xfit-raw | 0.0006 | 0.0065 | 0.029 | 0.467 | 0.597 |
| nb | 10 | 1.5 | limma | 0.0014 | 0.0138 | 0.077 | 0.351 | 0.509 |
| nb | 10 | 1.5 | edgeR | 0.0017 | 0.0132 | 0.082 | 0.362 | 0.527 |
| nb | 10 | 1.5 | TOX-log | 0.0081 | 0.0269 | 0.281 | 0.026 | 0.365 |
| nb | 10 | 1.5 | TOX-log-xfit-js | 0.0005 | 0.0049 | 0.027 | 0.285 | 0.484 |
| nb | 10 | 1.5 | TOX-log-xfit-raw | 0.0007 | 0.0044 | 0.030 | 0.320 | 0.496 |
| nb | 10 | 2 | limma | 0.0043 | 0.0277 | 0.179 | 0.266 | 0.436 |
| nb | 10 | 2 | edgeR | 0.0030 | 0.0190 | 0.131 | 0.281 | 0.460 |
| nb | 10 | 2 | TOX-log | 0.0217 | 0.0514 | 0.535 | 0.000 | 0.243 |
| nb | 10 | 2 | TOX-log-xfit-js | 0.0050 | 0.0189 | 0.199 | 0.066 | 0.347 |
| nb | 10 | 2 | TOX-log-xfit-raw | 0.0010 | 0.0046 | 0.052 | 0.186 | 0.409 |

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
- **TCGA n = 3, xfit-js:** null p ≤ 0.01 is nominal (0.0100), but FDR@0.05 in the power rounds
  is 0.167. The pi1 = 0 and pi1 = 0.1 rounds disagree, and why is not explained yet.
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

| data | part | n | lambda pooled | lambda case | lambda control |
|---|---|---|---|---|---|
| TCGA | RS | 3 | 0.453 | 0.642 | 0.634 |
| TCGA | RS | 5 | 0.302 | 0.478 | 0.501 |
| TCGA | RS | 7 | 0.240 | 0.428 | 0.401 |
| TCGA | RS | 10 | 0.167 | 0.297 | 0.308 |
| TCGA | R | 13 | 0.152 | 0.258 | 0.276 |
| TCGA | R | 20 | 0.110 | 0.221 | 0.205 |
| TCGA | R | 29 | 0.044 | 0.084 | 0.085 |
| TCGA | R | 32 | 0.067 | 0.120 | 0.135 |
| TCGA | R | 36 | 0.032 | 0.068 | 0.058 |
| TCGA | R | 41 | 0.053 | 0.106 | 0.098 |
| yeast | RS | 3 | 0.457 | 0.638 | 0.609 |
| yeast | RS | 5 | 0.290 | 0.466 | 0.510 |
| yeast | RS | 7 | 0.213 | 0.384 | 0.405 |
| yeast | RS | 10 | 0.152 | 0.292 | 0.275 |
| yeast | R | 21 | 0.056 | 0.101 | 0.113 |
| yeast | R | 22 | 0.086 | 0.180 | 0.147 |

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

