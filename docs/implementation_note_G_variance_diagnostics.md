# Implementation note G: tests on the data structure (variance-ratio diagnostic and cross-fitting on real data)

**For:** Claude Code. **Status:** DRAFT. **Companion:** note F (model changes and arms). Written without executing R; Python toy checks were run (section 5).
**Repo:** `empirical-noise-model-testing` (`power_test.R` stores what is needed, `power_report.R` / `common/` compute the diagnostics; no Fortran change).

---

## 0. Purpose (meeting with Asis, late October 2026)

1. **Evaluate cross-fitting on real data:** the arms of note F on Parts RS and R. This note fixes how the results are read (section 4, A5).
2. **Test the hypothesis before implementing anything:** genes whose noise changes between the two conditions (relative to the expected noise at their expression) are unreliable and produce false positives. Asis's proposed quantity:

```
r_g = ( s2_A,g / T_A(m_A,g) ) / ( s2_B,g / T_B(m_B,g) )        should be about 1
```

   collect it empirically, call the upper 5 % unreliable, and check whether its rank correlates with the false-positive rank.
3. Owner's additions (kept as specified): Spearman correlation between the LOESS-normalised variance ratio and the TOX p-values ("NNH p-values" in the meeting notes; **ASSUMPTION:** the p-values of each TOX arm, confirm) among known true negatives, in a **signed** and an **absolute** version, plus the false-positive rate within bins of the ratio. A correlation would indicate variance-dependent miscalibration; no correlation would not rule it out.
4. **Added (section 5 explains why):** a third quantity, the **overall** relative variance `v_g = log2(s2_g / T(m_g))`. In the toy it, not the imbalance, carries the false-positive signal of the current model.

---

## 1. Definitions (computed per round, per input type)

On `y = log2(x + 1)` of the input the TOX arm uses (TPM, TPM-centred, TMM-CPM; centring does not change variances):

```
group k in {A, B}:  s2_k,g (unbiased),  m_k,g = log2(mean_k,g(x) + 1)
T_k(m)              LOESS of s2_k,g on m_k,g over all tested genes (variance scale, degree 1, span 0.3, non-robust; floor VAR_FLOOR)
L_g = log2 r_g      = log2( s2_A,g / T_A(m_A,g) ) - log2( s2_B,g / T_B(m_B,g) )     signed imbalance
|L_g|                                                                                  absolute imbalance
v_g                 = log2( s2_g / T(m_g) ),  s2_g pooled, T = LOESS of s2_g on m_g = log2(mean over all samples + 1)   overall relative noise
```

* Each group is normalised by its own trend **at its own group mean**: variance differences explained by a mean difference (DE genes, mean-variance trend) are removed; a systematic group difference (e.g. tumour noisier overall) is absorbed by `T_A / T_B`. `L` therefore measures **gene-specific** imbalance only. The systematic part is reported separately (A4).
* **Chance thresholds (normal data, equal variances):** `s2_A / s2_B` follows F(n-1, n-1); the central 95 % ends at 39x (n = 3, |L| > 5.29), 9.6x (n = 5, 3.26), 4.0x (n = 10, 2.01), 2.5x (n = 20, 1.34). At n = 3 an "upper 5 %" flag therefore consists mostly of chance imbalance; real imbalance must be very large to be detectable.

---

## 2. Data and what must be stored

* **Primary:** Part RS (n = 3, 5, 7, 10; 100 rounds; six TCGA base cohorts). pi1 = 0 rounds: every gene is a true negative. pi1 = 0.1 rounds: the non-DE genes are the known negatives.
* Part R (full-cohort halves, n = 13-41): the same analyses at large n.
* **Reference ("null of the diagnostic"):** P+ with pi1 = 0 (simulated, equal variances, no outliers) for the same arms. Every real-data correlation is read against this baseline, because some arms can correlate with these quantities by construction.
* P+ rho extension (note F section 7): known systematic group imbalance.
* **Real contrasts** (healthy vs Stage IV per cancer; no truth): only the distribution of `L` and the group-trend ratio (A4); `group_variance_ratio_diagnostic.R` already does the trend part.
* **Storage:** `power_test.R` writes, per job and per input type, a method-independent file `vdiag_<job>.csv.gz` with `gene, input, s2_A, s2_B, mean_A, mean_B` (on the y-scale of that input). The report computes `T_A`, `T_B`, `T`, `L`, `v` from it. Per-gene p-values already exist in the raw results. Random halves: A and B are the two halves, kept as labelled by the sampler.

---

## 3. Arms analysed

`TOX-log|TPM-centred`, `TOX-log-vband-c3`, `TOX-log-xfit-raw`, `TOX-log-xfit-js` (and the other note F arms when present), plus limma, edgeR and DESeq2 as references (their p-values against the same quantities).

---

## 4. Analyses

**A1 Rank correlation (owner's specification, extended).** Per round, among true negatives:
`rho_signed = Spearman(L_g, p_g)`, `rho_abs = Spearman(|L_g|, p_g)`, `rho_overall = Spearman(v_g, p_g)` (average ranks for tied p-values at the floor).
Aggregate over rounds: mean +- SE (rounds are the replication unit, not genes: with thousands of genes even |rho| = 0.02 is "significant"). Report per cohort x n x arm, next to the P+ baseline. Sign convention: negative means a larger value goes with smaller p, i.e. more false positives.
* The signed version is meaningful when imbalance has a direction (real contrasts; in random halves an outlier sample lands in A or B at random, so signed effects cancel and the absolute version is the relevant one).

**A2 False-positive rate by bins.** Deciles of `|L|`, of `L`, and of `v`: share of true negatives with p <= 0.001 / 0.01 / 0.05, and for pi1 = 0.1 the share of BH calls (0.05) that are false, per bin. Enrichment = FPR in the top 5 % vs the rest, with a bootstrap CI over rounds.

**A3 Flag utility (only if A1/A2 show an effect).** Two flags: empirical (top 5 % of `|L|` in that round) and **null-calibrated** (threshold = 95th percentile of `|L|` in the pi1 = 0 rounds of the same cohort and n). On pi1 = 0.1 rounds: FDR and TPR with flagged genes removed vs kept, and the share of true DE genes flagged (the cost).

**A4 Real contrasts (no truth).** Distribution of `L` and the systematic ratio `T_A(m)/T_B(m)` by expression decile; share of genes beyond the null-calibrated threshold.

**A5 Cross-fitting on real data.** Note F's arms on RS and R with the note D / E v2 criteria (calibration and tail gates on pi1 = 0, then FDR and power), per n and cohort and by expression decile; plus A1-A2 per arm, to see whether the cross-fit removes the variance dependence of the false positives (toy: yes, section 5).

**Outputs:** `variance_diagnostic_corr.csv` (A1), `variance_diagnostic_bins.csv` (A2), `variance_flag_utility.csv` (A3), `variance_ratio_real_contrasts.csv` (A4); plots: FPR vs decile of `|L|`, `L` and `v` per n (one panel per arm), rho summary per cohort and n.

**Tests:** synthetic check that equal variances give |rho| < 0.03 for all three quantities; an injected gene-specific imbalance gives a negative `rho_signed` for the current model; `vdiag` files contain one row per gene and input; the LOESS fits are finite.

---

## 5. What the toy predicts (sanity, NOT targets)

`toy_variance_ratio.py`, `toy_variance_ratio2.py` (pure null, normal noise, mean-variance trend, sigma_v = 1.05, 3000 genes, 4 datasets per cell):

| scenario | current TOX: rho_signed / rho_abs | FPR(p <= 0.01) top-5 % vs rest |
|---|---|---|
| equal variances | about 0 / about 0 (|rho| <= 0.02) | no enrichment |
| gene-specific imbalance (case SD x exp(N(0, 0.35))) | -0.06 to -0.09 / about 0 | no enrichment (direction cancels in \|L\|) |
| outlier sample in 5 % of genes (one group) | -0.02 to -0.06 / about 0 | none at n = 3, up to 1.8x at n = 10 |

**Overall relative variance `v` (equal variances):** current TOX rho_overall -0.36 / -0.41 / -0.45 at n = 3 / 5 / 10; FPR(p <= 0.01) in the noisiest 5 % 0.10-0.13 vs 0.006-0.009 in the rest (11-22x). Cross-fit + James-Stein: rho_overall -0.01 to -0.04, no enrichment. Cross-fit without shrinkage: reversed sign (conservative).

**Reading rules (pre-registered):**
* Weak imbalance correlations and a strong overall-variance correlation (the toy pattern) mean that the main false-positive driver is gene-wise heterogeneity, which cross-fitting addresses, and that an imbalance flag would add little.
* Strong real-data imbalance correlations, clearly above the P+ baseline, especially at n >= 10 where the ratio is informative, would support an imbalance flag; then A3 decides whether it helps.
* At n = 3 the ratio is mostly chance (39x threshold): an absent correlation there says little.

---

## 6. Decision rule for an imbalance flag (PROPOSAL)

Implement a flag (as an annotation, not a filter, in the first version) only if, at a given n: FPR enrichment in the null-calibrated top 5 % is at least 2x with a bootstrap CI excluding 1 in at least 4 of 6 cohorts, AND removing flagged genes lowers observed FDR without lowering TPR by more than 0.02 (A3). Otherwise document the diagnostic as a negative result.
