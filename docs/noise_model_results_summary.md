# Noise model evaluation
As of 23/08/2026

----
![p-value distributions](./calibration_out/null_pvalue_qqplots.png)

pvalue qq-plots of all tools under different distributions (n_rep = {3,5})

----
### Simulated data
#### Negative Binomial (`nb`)
This section represents the exact distibution assumed behind edgeR and DESeq2.
All models/parameters are well calibrated for the nb distribution and show no visible deviation from the expected diagonal line.

#### Lognormal multiplier (`lnpois`)
All models are well calibrated and show no visible deviation from the expected diagonal line

#### Bernoulli (`bimodal`)
DESeq2 shows highly conservative behaviour, identifying large amounts of negatives where positives should be. The assumed negative binomal is no longer truly given, causing the expected failure of the tool. edgeR results in a comparable behaviour.
Limma performs slightly conservative, with slightly too high p values, but performances with `n_rep` = 5 is better than with `n_rep` = 3.

##### Different TOX models
The performance of tensor omics varies between the chosen model, parameters and the number of replicates. 
The bootstrap model performs slightly different, with the raw normalisation exhibiting less deviation from the diagonal than the log normalisation, while still showing a visible S-Shape. The S-Shape becomes more clear with more replicates and with larger neighborhoods, indicating that the neighborhood growing causes some instability for this underlying distribution.

The `exact` model on the other hand shows a well calibrated null under the raw normalisation, and only slight deviations under the log normalisation. Neighborhood size does not seem to affect the result. THe deviation from the diagonal increases with `n_rep = 5`, and becomes comparable between the raw and log normalisation.

### Data driven results summary
#### Negative Binomial

```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
                  DESeq2      nb     3       0.0651       0.0155     4.2708   0.4847 0.0000       0.60         NaN         NaN        NaN         NaN
                   edgeR      nb     3       0.0507       0.0092     5.5683   0.4986 0.0003       0.00         NaN         NaN        NaN         NaN
                   limma      nb     3       0.0485       0.0082     6.0157   0.4941 0.0004       0.00         NaN         NaN        NaN         NaN
 TOX_bootstrap_n0_k20_50      nb     3       0.0437       0.0049     9.0716   0.4830 0.0011       0.00           1           0    46.1511       150.0
 TOX_bootstrap_n0_k30_70      nb     3       0.0418       0.0040    10.5440   0.4858 0.0010       0.00           1           0    49.2420       210.0
 TOX_bootstrap_n0_k40_70      nb     3       0.0416       0.0038    10.9914   0.4859 0.0010       0.00           1           0    49.2420       210.0
 TOX_bootstrap_n1_k20_50      nb     3       0.0568       0.0122     4.7686   0.4924 0.0002       0.00           1           0    77.2331       150.0
 TOX_bootstrap_n1_k30_70      nb     3       0.0538       0.0113     4.8736   0.4964 0.0002       0.00           1           0    69.8272       210.0
 TOX_bootstrap_n1_k40_70      nb     3       0.0532       0.0110     4.9629   0.4971 0.0002       0.00           1           0    69.8272       210.0
     TOX_exact_n0_k20_50      nb     3       0.0451       0.0061     7.6779   0.4854 0.0001       0.00           1           0   107.6956       150.0
     TOX_exact_n0_k30_70      nb     3       0.0431       0.0049     8.9023   0.4871 0.0007       0.00           1           0    88.3001       210.0
     TOX_exact_n0_k40_70      nb     3       0.0429       0.0047     9.2116   0.4874 0.0007       0.00           1           0    88.3001       210.0
     TOX_exact_n1_k20_50      nb     3       0.0566       0.0118     4.9854   0.4916 0.0001       0.00           1           0   126.6531       150.0
     TOX_exact_n1_k30_70      nb     3       0.0539       0.0111     4.9501   0.4952 0.0001       0.00           1           0   152.8891       210.0
     TOX_exact_n1_k40_70      nb     3       0.0532       0.0108     5.0365   0.4966 0.0001       0.00           1           0   144.1757       210.0
     
```

with:
`frac_p_lt_05` being the fraction of p-values below 0.05. Target is 0.05 -> uniform distribution
`frac_p_lt_01` being the fraction of p-values below 0.01. Target is 0.01 -> uniform distribution
`tail_ration` being `frac_p_lt_05 / frac_p_lt_01` to answer if most values below 0.05 are also below 0.01 -> spike near zero. Target is 5
`median_p` being the median p value
`min_p` being the minimum p value
`n_rejected` average number of "positive" outliers across 5 runs.

Result: Under a negative binomial edgeR has the best calibration, with limma being slightly conservative and DESeq2 being slightly anti-conservative both `frac_p` below/above expected value. Under a log normalisation, the bootstrap shows very well calibration, with `frac_p` and `tail_ratio` being both very close to the expected value. The raw normalisation shows an overall acceptable calibration, but with `frac_p_lt_01` being only at ~0.004, while it should be at 0.01. The model therefore performs more conservative

The exact model shows nearly identical results and similar normalisation-dependent behaviour. Results are comparable at `n_rep = 5`

#### Lognormal multiplier

```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
                  DESeq2  lnpois     3       0.0702       0.0184     3.8134   0.4713 0.0000       0.80         NaN         NaN        NaN         NaN
                   edgeR  lnpois     3       0.0573       0.0124     4.6559   0.4872 0.0002       0.00         NaN         NaN        NaN         NaN
                   limma  lnpois     3       0.0487       0.0095     5.2188   0.5061 0.0003       0.00         NaN         NaN        NaN         NaN
 TOX_bootstrap_n0_k20_50  lnpois     3       0.0482       0.0063     7.7619   0.4777 0.0008       0.00           1           0    46.3676       150.0
 TOX_bootstrap_n0_k30_70  lnpois     3       0.0463       0.0056     8.3527   0.4804 0.0011       0.00           1           0    43.8424       210.0
 TOX_bootstrap_n0_k40_70  lnpois     3       0.0461       0.0055     8.4210   0.4807 0.0011       0.00           1           0    43.8424       210.0
 TOX_bootstrap_n1_k20_50  lnpois     3       0.0520       0.0121     4.3262   0.5017 0.0002       0.00           1           0    85.6873       150.0
 TOX_bootstrap_n1_k30_70  lnpois     3       0.0504       0.0109     4.6253   0.5039 0.0002       0.00           1           0    82.4243       210.0
 TOX_bootstrap_n1_k40_70  lnpois     3       0.0501       0.0109     4.6176   0.5049 0.0001       0.00           1           0    88.2829       210.0
     TOX_exact_n0_k20_50  lnpois     3       0.0484       0.0067     7.3144   0.4788 0.0002       0.00           1           0   100.2209       150.0
     TOX_exact_n0_k30_70  lnpois     3       0.0474       0.0060     8.0316   0.4806 0.0007       0.00           1           0    72.8840       210.0
     TOX_exact_n0_k40_70  lnpois     3       0.0472       0.0059     8.1570   0.4811 0.0009       0.00           1           0    66.7121       210.0
     TOX_exact_n1_k20_50  lnpois     3       0.0524       0.0128     4.1071   0.5028 0.0001       0.00           1           0   125.0697       150.0
     TOX_exact_n1_k30_70  lnpois     3       0.0505       0.0115     4.4117   0.5055 0.0001       0.00           1           0   135.5096       210.0
     TOX_exact_n1_k40_70  lnpois     3       0.0501       0.0113     4.4431   0.5064 0.0001       0.00           1           0   146.1783       210.0
     
```

with:
`frac_p_lt_05` being the fraction of p-values below 0.05. Target is 0.05 -> uniform distribution
`frac_p_lt_01` being the fraction of p-values below 0.01. Target is 0.01 -> uniform distribution
`tail_ration` being `frac_p_lt_05 / frac_p_lt_01` to answer if most values below 0.05 are also below 0.01 -> spike near zero. Target is 5
`median_p` being the median p value
`min_p` being the minimum p value
`n_rejected` average number of "positive" outliers across 5 runs.

Overall all models are well calibrated, with the raw normalisation showing slight issues in the `frac_p_lt_01`.

#### Bernoulli

```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
                  DESeq2 bimodal     3       0.0017       0.0016     1.0650   0.9311 0.0000       6.80         NaN         NaN        NaN         NaN
                   edgeR bimodal     3       0.0020       0.0014     1.4157   0.7650 0.0000       6.60         NaN         NaN        NaN         NaN
                   limma bimodal     3       0.0384       0.0025    15.8246   0.3151 0.0000       6.40         NaN         NaN        NaN         NaN
 TOX_bootstrap_n0_k20_50 bimodal     3       0.0305       0.0085     3.6161   0.7298 0.0001       0.00           1           0   100.0000        90.0
 TOX_bootstrap_n0_k30_70 bimodal     3       0.0305       0.0066     4.7478   0.7453 0.0001       0.00           1           0    94.1414       156.0
 TOX_bootstrap_n0_k40_70 bimodal     3       0.0295       0.0051     5.9734   0.7576 0.0002       0.00           1           0    82.3044       210.0
 TOX_bootstrap_n1_k20_50 bimodal     3       0.0511       0.0235     2.1859   0.7429 0.0001       0.00           1           0   100.0000        90.0
 TOX_bootstrap_n1_k30_70 bimodal     3       0.0436       0.0226     1.9468   0.7572 0.0001       0.00           1           0   100.0000       150.0
 TOX_bootstrap_n1_k40_70 bimodal     3       0.0405       0.0213     1.9204   0.7685 0.0001       0.00           1           0   100.0000       210.0
     TOX_exact_n0_k20_50 bimodal     3       0.0126       0.0017     7.7622   0.5176 0.0001       0.00           1           0   127.0197        90.0
     TOX_exact_n0_k30_70 bimodal     3       0.0094       0.0017     5.7144   0.5273 0.0000       0.00           1           0   186.6530       156.0
     TOX_exact_n0_k40_70 bimodal     3       0.0079       0.0017     4.6895   0.5370 0.0000       0.80           1           0   193.3498       210.0
     TOX_exact_n1_k20_50 bimodal     3       0.0382       0.0130     2.9850   0.5256 0.0000       0.00           1           0   148.4605        90.0
     TOX_exact_n1_k30_70 bimodal     3       0.0373       0.0110     3.4536   0.5350 0.0000       3.60           1           0   210.0000       150.0
     TOX_exact_n1_k40_70 bimodal     3       0.0370       0.0087     4.4135   0.5460 0.0000       5.80           1           0   210.0000       210.0
     
```

DESeq2 shows a clear miscalibration with `frac_p_lt_05 = frac_p_lt_01 = 0.0017`. So almost every pvalue below 0.05 is also below 0.01, but the fraction for both is way too small. DESeq2 therefore results in extremly large p values with a clear spike near 1 (median_p = 0.9311). But it also rejects the null ("true" positive) for an average of 6.8 genes.

EdgeR shows a misscalibration comparable to DESeq2 with a non-uniform p value distribution and a spike in 1 direction (median_p = 0.765). It also identifies "true positives".

Limma shows a far better calibration with `frac_p_lt_05` being close to 0.04, so only slightly off-chart. `frac_p_lt_01` sits around 0.0025, a quarter of what it should be. the median p on the other hand is at 0.31, suggesting that a lot of p-values sit in lower ranges, while only few could be considered significant.

TOX performance largly depends on model and normalisation. 

----
Bootstrap:
```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
 TOX_bootstrap_n0_k20_50 bimodal     3       0.0305       0.0085     3.6161   0.7298 0.0001       0.00           1           0   100.0000        90.0
 TOX_bootstrap_n0_k30_70 bimodal     3       0.0305       0.0066     4.7478   0.7453 0.0001       0.00           1           0    94.1414       156.0
 TOX_bootstrap_n0_k40_70 bimodal     3       0.0295       0.0051     5.9734   0.7576 0.0002       0.00           1           0    82.3044       210.0
 TOX_bootstrap_n1_k20_50 bimodal     3       0.0511       0.0235     2.1859   0.7429 0.0001       0.00           1           0   100.0000        90.0
 TOX_bootstrap_n1_k30_70 bimodal     3       0.0436       0.0226     1.9468   0.7572 0.0001       0.00           1           0   100.0000       150.0
 TOX_bootstrap_n1_k40_70 bimodal     3       0.0405       0.0213     1.9204   0.7685 0.0001       0.00           1           0   100.0000       210.0
 
```

The raw normalisation shows `frac_p_lt_05` comparable to limma at around 0.03 and `frac_p_lt_01` at around 0.005 - 0.008 with values decreasing with growing neighborhoods. 
The log normalisation on the other hand shows `frac_p_lt_05` at around 0.04 - 0.05 and `frac_p_lt_01` at around 0.022 with values decreasing with growing neighborhoods. Interestingly, the log normalisation exhibits more pvalues below 0.01 than we would expect, while all other methods show a different behaviour. The fraction below 0.05 is around the exepcted value, indicating a slight peak towards the 0.01, visible in the tail ratio of around 2. Under both normalisations the median_p sits around 0.7, indicating that the majority of pvalues is larger than expected.

---
Exact:
```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
     TOX_exact_n0_k20_50 bimodal     3       0.0126       0.0017     7.7622   0.5176 0.0001       0.00           1           0   127.0197        90.0
     TOX_exact_n0_k30_70 bimodal     3       0.0094       0.0017     5.7144   0.5273 0.0000       0.00           1           0   186.6530       156.0
     TOX_exact_n0_k40_70 bimodal     3       0.0079       0.0017     4.6895   0.5370 0.0000       0.80           1           0   193.3498       210.0
     TOX_exact_n1_k20_50 bimodal     3       0.0382       0.0130     2.9850   0.5256 0.0000       0.00           1           0   148.4605        90.0
     TOX_exact_n1_k30_70 bimodal     3       0.0373       0.0110     3.4536   0.5350 0.0000       3.60           1           0   210.0000       150.0
     TOX_exact_n1_k40_70 bimodal     3       0.0370       0.0087     4.4135   0.5460 0.0000       5.80           1           0   210.0000       210.0
     
```
The raw normalisation shows values in the `frac_p_lt_05` of around 0.008 - 0.012 and a value of 0.0017 in `frac_p_lt_01`. In both cases, the values are far below the expected value. On the other hand, the median p sits around 0.5, indicating that the mis-calibration may only affect smaller pvalues, leading to a conservative behaviour. The largest knn configuration identified a "true positive" in 4/5 runs, indicating a sligthly data depended behaviour.

The log normalisation on the other hand seems to be slightly betetr calibrated with `frac_p_lt_05` being around 0.037 and `frac_p_lt_01` being around 0.0087 - 0.013. Median p-values are around 0.5, while growing neighborhoods tend to identify a growing number of "outliers"

---

The differences become more clear when looking at the n_rep = 5 results.
```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
                  DESeq2 bimodal     5       0.0001       0.0000     1.0000   0.8003 0.0713       0.20         NaN         NaN        NaN         NaN
                   edgeR bimodal     5       0.0005       0.0003     1.5000   0.7849 0.0157       0.00         NaN         NaN        NaN         NaN
                   limma bimodal     5       0.0210       0.0024     9.3909   0.3390 0.0017       0.00         NaN         NaN        NaN         NaN
 TOX_bootstrap_n0_k20_50 bimodal     5       0.0506       0.0112     4.5876   0.4913 0.0001       0.00           1           0    91.5459       210.0
 TOX_bootstrap_n0_k30_70 bimodal     5       0.0452       0.0099     4.5826   0.5009 0.0002       0.00           1           0    79.8287       350.0
 TOX_bootstrap_n0_k40_70 bimodal     5       0.0403       0.0094     4.3204   0.5073 0.0002       0.00           1           0    79.8287       350.0
 TOX_bootstrap_n1_k20_50 bimodal     5       0.0796       0.0117     6.9380   0.3775 0.0001       0.00           1           0   100.0000       195.0
 TOX_bootstrap_n1_k30_70 bimodal     5       0.0769       0.0110     7.0970   0.3819 0.0001       0.00           1           0    94.1414       350.0
 TOX_bootstrap_n1_k40_70 bimodal     5       0.0747       0.0106     7.1377   0.3873 0.0001       0.00           1           0    94.1414       350.0
     TOX_exact_n0_k20_50 bimodal     5       0.0172       0.0089     1.9571   0.3268 0.0002       0.00           1           0   100.8979       210.0
     TOX_exact_n0_k30_70 bimodal     5       0.0146       0.0083     1.8079   0.3311 0.0002       0.20           1           0   128.2933       350.0
     TOX_exact_n0_k40_70 bimodal     5       0.0139       0.0080     1.7569   0.3373 0.0002       0.20           1           0   128.2933       350.0
     TOX_exact_n1_k20_50 bimodal     5       0.0201       0.0105     1.9343   0.2794 0.0002       0.00           1           0   103.1104       195.0
     TOX_exact_n1_k30_70 bimodal     5       0.0154       0.0105     1.4826   0.2793 0.0001       0.20           1           0   130.2112       350.0
     TOX_exact_n1_k40_70 bimodal     5       0.0135       0.0105     1.2978   0.2823 0.0002       0.20           1           0   128.9252       350.0
```

While DESeq2 and EdgeR show very highly conservative behaviour, limma shows a slightly less conservative behvaiour, but still a mis-calibration

---
Bootstrap:
```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
 TOX_bootstrap_n0_k20_50 bimodal     5       0.0506       0.0112     4.5876   0.4913 0.0001       0.00           1           0    91.5459       210.0
 TOX_bootstrap_n0_k30_70 bimodal     5       0.0452       0.0099     4.5826   0.5009 0.0002       0.00           1           0    79.8287       350.0
 TOX_bootstrap_n0_k40_70 bimodal     5       0.0403       0.0094     4.3204   0.5073 0.0002       0.00           1           0    79.8287       350.0
 TOX_bootstrap_n1_k20_50 bimodal     5       0.0796       0.0117     6.9380   0.3775 0.0001       0.00           1           0   100.0000       195.0
 TOX_bootstrap_n1_k30_70 bimodal     5       0.0769       0.0110     7.0970   0.3819 0.0001       0.00           1           0    94.1414       350.0
 TOX_bootstrap_n1_k40_70 bimodal     5       0.0747       0.0106     7.1377   0.3873 0.0001       0.00           1           0    94.1414       350.0
 ```

 The raw normalisation shows a well-calibration in every part, with all values being close to the expected ones, while the log normalisation shows a clear anti-conservative behaviour with `frac_p_lt_05` being around 0.75 and ad a median p of ~0.38, indicating overall low p values. On the other hand, `frac_p_lt_01` stays at almost 0.01
 
---
Exact:
```
                  method    dist n_rep frac_p_lt_05 frac_p_lt_01 tail_ratio median_p  min_p n_rejected frac_tested frac_fallbk pool_infer nbhd_median
     TOX_exact_n0_k20_50 bimodal     5       0.0172       0.0089     1.9571   0.3268 0.0002       0.00           1           0   100.8979       210.0
     TOX_exact_n0_k30_70 bimodal     5       0.0146       0.0083     1.8079   0.3311 0.0002       0.20           1           0   128.2933       350.0
     TOX_exact_n0_k40_70 bimodal     5       0.0139       0.0080     1.7569   0.3373 0.0002       0.20           1           0   128.2933       350.0
     TOX_exact_n1_k20_50 bimodal     5       0.0201       0.0105     1.9343   0.2794 0.0002       0.00           1           0   103.1104       195.0
     TOX_exact_n1_k30_70 bimodal     5       0.0154       0.0105     1.4826   0.2793 0.0001       0.20           1           0   130.2112       350.0
     TOX_exact_n1_k40_70 bimodal     5       0.0135       0.0105     1.2978   0.2823 0.0002       0.20           1           0   128.9252       350.0
```

The exact models results indicate an anti-conservative behaviour with p values sitting around 0.3. Under both normalisations, the median p is around 0.3, with the raw normalisation at around 0.33 and log around 0.27. The `frac_p_lt_05` is around 0.015 - 0.02 in all cases. So in general, the fraction of pvalues below 0.05 is lower than expected, but the median p is also lower than expected, indicating that pvalues are smaller than they should be, while they do not pass the significance threshold in most cases.



# Mock-null uniformity diagnostics — TOX vs edgeR / limma / DESeq2
(Note: This section was generated and aggregated by Claude and verified by me)

TOX rows collapse the three neighbourhood configs (`k20_50`, `k30_70`, `k40_70`) into a min–max range, since the earlier k-sweep showed negligible sensitivity to k — the range width itself is evidence of that. `edgeR` / `limma` / `DESeq2` have one config each and stay as single rows. `ks_p`, `ad_p`, `cvm_p` are omitted: per the script's own docstring they saturate near the resolution floor at this gene count and carry no information beyond the statistics themselves.

## Column reference

Every column below tests the same question — **is this p-value vector Uniform(0,1)?** — under the mock null, where every gene is truly non-DE and every rejection is a false positive.

| column | what it is | target |
|---|---|---|
| `median_p` | median of the raw p-values | **0.5** |
| `ks_D` | Kolmogorov–Smirnov statistic: the largest gap between the empirical CDF and the uniform CDF, `sup_x \|F_n(x) − x\|` | **0**, smaller = better |
| `ad_A2` | Anderson–Darling statistic: like KS but weighted `1/[F(x)(1−F(x))]`, so far more sensitive to deviations near 0 and 1 | **0**, smaller = better |
| `cvm_W2` | Cramér–von Mises statistic: integrated squared deviation across the whole range, `∫[F_n(x)−x]² dx` | **0**, smaller = better |
| `wilx_p` | one-sample Wilcoxon signed-rank p-value testing only whether the median equals 0.5 — a pure **location** check, not a shape check | **large** (not rejected) |

Three statistics, three sensitivities: KS catches the single worst deviation anywhere; CvM catches deviations spread across the range; AD catches deviations concentrated at the extremes — which is where BH reads and therefore the one most relevant to false-discovery behaviour. `wilx_p` is the odd one out: it only checks centring, so a distribution that is badly shaped but happens to be centred near 0.5 would pass it while failing all three others.

---

## Table 1 — `nb` (negative binomial: edgeR / DESeq2's exact generative model)

| n_rep | group | median_p | ks_D | ad_A2 | cvm_W2 | wilx_p |
|---|---|---|---|---|---|---|
| 3 | DESeq2 | 0.480 | 0.0264 | 8.27 | 1.47 | 1.3e-04 |
| 3 | edgeR | 0.494 | 0.0131 | 1.26 | 0.206 | 2.7e-01 |
| 3 | limma | 0.491 | 0.0221 | 3.46 | 0.641 | 2.2e-02 |
| 3 | TOX_exact_n0 | 0.479–0.480 | 0.0244–0.0278 | 5.48–5.93 | 0.995–1.16 | 3.0e-03–9.4e-03 |
| 3 | **TOX_exact_n1** | 0.479–0.484 | 0.0178–0.0232 | 2.93–5.02 | 0.569–0.944 | 2.9e-03–2.4e-02 |
| 3 | TOX_bootstrap_n0 | 0.475–0.476 | 0.0263–0.0301 | 6.54–7.50 | 1.18–1.46 | 9.1e-04–4.6e-03 |
| 3 | TOX_bootstrap_n1 | 0.479–0.485 | 0.0171–0.0226 | 2.45–4.56 | 0.465–0.854 | 4.5e-03–4.1e-02 |
| 5 | DESeq2 | 0.496 | 0.0214 | 4.20 | 0.314 | 1.7e-01 |
| 5 | edgeR | 0.505 | 0.0138 | 0.72 | 0.128 | 4.5e-01 |
| 5 | limma | 0.516 | 0.0177 | 1.11 | 0.248 | 2.7e-01 |
| 5 | TOX_exact_n0 | 0.481–0.482 | 0.0198–0.0206 | 2.13–2.43 | 0.398–0.451 | 8.2e-02–1.1e-01 |
| 5 | **TOX_exact_n1** | 0.508 | 0.0123–0.0130 | 0.81–0.85 | 0.155–0.175 | 4.8e-01–5.6e-01 |
| 5 | TOX_bootstrap_n0 | 0.490–0.491 | 0.0109–0.0128 | 0.73–0.78 | 0.104–0.117 | 5.5e-01–6.0e-01 |
| 5 | TOX_bootstrap_n1 | 0.511–0.513 | 0.0150–0.0159 | 1.19–1.30 | 0.247–0.281 | 2.7e-01–3.2e-01 |

## Table 2 — `lnpois` (identical mean–variance to `nb`, lognormal instead of gamma shape)

| n_rep | group | median_p | ks_D | ad_A2 | cvm_W2 | wilx_p |
|---|---|---|---|---|---|---|
| 3 | DESeq2 | 0.478 | 0.0296 | 14.3 | 1.92 | 1.0e-05 |
| 3 | edgeR | 0.494 | 0.0148 | 2.47 | 0.314 | 8.6e-02 |
| 3 | limma | 0.514 | 0.0163 | 2.05 | 0.361 | 8.2e-02 |
| 3 | TOX_exact_n0 | 0.485–0.487 | 0.0164–0.0173 | 2.06–2.79 | 0.380–0.533 | 2.5e-02–6.7e-02 |
| 3 | **TOX_exact_n1** | 0.504–0.509 | 0.0131–0.0149 | 1.56–1.98 | 0.227–0.342 | 1.3e-01–3.4e-01 |
| 3 | TOX_bootstrap_n0 | 0.487 | 0.0153–0.0174 | 2.04–3.09 | 0.380–0.599 | 1.9e-02–6.9e-02 |
| 3 | TOX_bootstrap_n1 | 0.503–0.508 | 0.0124–0.0139 | 1.22–1.54 | 0.166–0.262 | 1.9e-01–5.0e-01 |
| 5 | DESeq2 | 0.470 | 0.0376 | 14.9 | 2.58 | 4.4e-07 |
| 5 | edgeR | 0.482 | 0.0227 | 4.55 | 0.846 | 4.1e-03 |
| 5 | limma | 0.505 | 0.0161 | 1.77 | 0.250 | 1.6e-01 |
| 5 | TOX_exact_n0 | 0.466–0.467 | 0.0373–0.0380 | 11.8–12.3 | 2.43–2.53 | 5.1e-06–8.0e-06 |
| 5 | **TOX_exact_n1** | 0.503 | 0.0159–0.0162 | 1.50–1.57 | 0.236–0.249 | 1.8e-01–2.1e-01 |
| 5 | TOX_bootstrap_n0 | 0.480–0.481 | 0.0225–0.0226 | 3.55–3.57 | 0.716–0.717 | 2.4e-02 |
| 5 | TOX_bootstrap_n1 | 0.501–0.502 | 0.0145–0.0159 | 1.21–1.31 | 0.183–0.200 | 2.7e-01–3.3e-01 |

## Table 3 — `bimodal` (structural zeros, negative skew, no method's home turf)

| n_rep | group | median_p | ks_D | ad_A2 | cvm_W2 | wilx_p |
|---|---|---|---|---|---|---|
| 3 | DESeq2 | 0.927 | 0.5790 | 4570 | 699 | 0.0e+00 |
| 3 | edgeR | 0.766 | 0.4810 | 2100 | 446 | 0.0e+00 |
| 3 | limma | 0.316 | 0.1960 | 253 | 49.4 | 4.1e-40 |
| 3 | **TOX_exact_n0** | 0.439–0.461 | **0.0839–0.0920** | 88.6–101 | 13.2–17.1 | 1.1e-08–2.3e-22 |
| 3 | TOX_exact_n1 | 0.440–0.462 | 0.154–0.160 | 120–153 | 21.2–28.2 | 9.1e-31–3.3e-48 |
| 3 | TOX_bootstrap_n0 | 0.672–0.711 | 0.175–0.212 | 193–300 | 37.8–57.1 | 0.0e+00 |
| 3 | TOX_bootstrap_n1 | 0.686–0.723 | 0.188–0.227 | 268–334 | 51.5–64.4 | 7.1e-04–1.2e-14 |
| 5 | DESeq2 | 0.804 | 0.5490 | 3370 | 590 | 0.0e+00 |
| 5 | edgeR | 0.785 | 0.5110 | 2540 | 515 | 0.0e+00 |
| 5 | limma | 0.340 | 0.1670 | 250 | 30.3 | 5.7e-06 |
| 5 | **TOX_exact_n0** | 0.342–0.361 | **0.163–0.169** | 162–163 | 28.4–29.9 | 1.7e-04–5.1e-01 |
| 5 | TOX_exact_n1 | 0.289–0.295 | 0.207–0.212 | 228 | 41.6–42.9 | 1.3e-03–9.0e-01 |
| 5 | TOX_bootstrap_n0 | 0.502–0.518 | 0.238–0.245 | 589–664 | 67.3–76.4 | 0.0e+00 |
| 5 | TOX_bootstrap_n1 | 0.394–0.408 | 0.296–0.307 | 915–984 | 97.1–101 | 0.0e+00 |

---

## Main findings

**TOX matches edgeR on `nb` — edgeR's own generative model — at n=5.**  `TOX_exact_n1`'s `ks_D` (0.012–0.013) and `wilx_p` (0.48–0.56) sit essentially on top of edgeR's (0.0138, 0.447). At n=3 TOX is somewhat worse (`ks_D` 0.018–0.023 vs edgeR's 0.013), so the parity is an n=5 result specifically.

**`norm_method` is the dominant factor on `lnpois`, not `k` or `variant`.** `TOX_exact_n0`'s `ks_D` roughly doubles from n=3 to n=5 (0.017 → 0.038), tracking DESeq2's own degradation over the same change (0.030 → 0.038). `TOX_exact_n1` stays flat (0.013–0.016 across both) and sits closest to limma throughout. Since `nb` and `lnpois` share the identical mean–variance relationship and differ only in mixing-distribution shape, this says TOX's `norm_method = 1` mode is not sensitive to that shape difference, while `norm_method = 0` is.

**On `bimodal`, every TOX arm beats DESeq2 and edgeR by roughly an order of magnitude, and `TOX_exact_n0` is the best arm in the table.** `ks_D` for `TOX_exact_n0` is 0.084–0.17 against DESeq2's 0.55 and edgeR's 0.48–0.51 across both replicate counts. `TOX_exact_n0` at n=5 is also the only group with `wilx_p` reaching non-rejection (up to 0.51 and 0.90 across its k-configs) — every other group in this table, TOX included, rejects the location test outright.

**`TOX_exact_n0` and `TOX_exact_n1` trade places depending on the arm.** `n1` is better calibrated on `nb`/`lnpois`; `n0` is markedly better on `bimodal` (`ks_D` 0.08–0.17 vs `n1`'s 0.15–0.21, at both replicate counts). No single `norm_method` dominates across all three distributions in this table.

**The bootstrap variant is consistently worse than exact wherever they can be compared.** On every distribution and replicate count, `TOX_bootstrap_*`'s `ks_D`, `ad_A2`, and `cvm_W2` are equal to or larger than the matching `TOX_exact_*` row — most sharply on `bimodal`, where bootstrap's `ad_A2` (193–984) runs 2–8× exact's (89–228) depending on cell.

**The three shape statistics agree with each other far more than any of them agrees with `wilx_p`, and this matters most on `bimodal`.** DESeq2 and edgeR both fail `ks_D`, `ad_A2`, and `cvm_W2` badly there, yet neither `wilx_p` nor `median_p` would reveal *why* on their own — `wilx_p = 0` just says "not centred at 0.5," not "spiked" versus "depleted." Reading the shape statistics alongside the earlier `tail_ratio` diagnostic (not reproduced here) is what distinguishes DESeq2/edgeR's spike-at-zero pattern from limma's tail-depletion pattern on the same data — a distinction none of these five columns makes on its own.

**Caveat on the ranges.** Because k has little effect, a narrow range (e.g. `TOX_exact_n1` `ad_A2` on `nb` n=5: 0.81–0.85) is safe to read as a single number. A wide range (e.g. `TOX_bootstrap_n0` `ad_A2` on `bimodal` n=3: 193–300) means the three k-configs disagree enough that the per-k breakdown is worth pulling before drawing a conclusion at that level of detail.

## Real Null Data

![Null p values histograms](./null_calibration/plots/null_pvalue_histograms.png)

![Null p values qq plots](./null_calibration/plots/null_pvalue_qqplots.png)

# Real-data mock null — TCGA cancer cohorts

Unlike everything in the prior tables, this null is not simulated. Each row splits the samples from **one clinical group** (a single cancer type at a single stage) into two random halves, runs the full method, and repeats 100 times. Since both halves come from the same population, every gene is truly null by construction — this is the strongest calibration test available, because it uses no simulator assumptions at all, only real expression heterogeneity.

## Column reference

| column | what it is | target |
|---|---|---|
| `n_samples` | total samples in the group, before the 100 random halvings | — |
| `n_genes` | genes surviving each method's own filter | — (see note below) |
| `FPR_0.05` | mean, over 100 runs, of the fraction of genes with raw p < 0.05 | **0.05** |
| `FPR_0.01` | same at the 1% level | **0.01** |
| `median_p` | mean, over 100 runs, of the median raw p-value | **0.5** |
| `ks_D` | Kolmogorov–Smirnov distance of the p-values to Uniform(0,1), averaged over runs | **0** |
| `ad_A2` | Anderson–Darling statistic (as before, weighted toward the tails) | **0** |
| `hits_FDR05` | mean genes surviving BH q < 0.05 **per run** | **0** |
| `inflation_0.05` | `FPR_0.05 / 0.05` | **≈1** |
| `verdict` | categorical label from the script | — |

**On `n_genes`:** TOX tests 12,287–12,727 genes per cancer type depending on cancer, versus 15,309–16,276 for edgeR/limma — TOX's filter is consistently ~15–20% stricter. This is a real difference in gene universe, not just labelling, and it matters for interpreting the comparison (see Findings).

**On `hits_FDR05`:** confirmed from source (`round(mean(hits_fdr05), 4)`, averaged directly over the `N_SPLITS` runs) — a genuine per-run mean, not a fraction of runs. E.g. edgeR's 0.0343 means roughly one run in thirty produces a single BH-significant gene purely from chance. `TOX-raw`'s value is exactly 0.0000 in **all twelve** groups (`hits_FDR05_sd` is non-zero in two of them — 0.0001, 0.0001 — meaning a handful of individual runs did produce one hit, but rarely enough that the mean rounds to zero).

**On `verdict`:** confirmed from `null_calibration.R`: a literal fixed cutoff, `inflation_0.05 > 1.5 → "ANTI-CONSERVATIVE"`, `inflation_0.05 < 0.67 → "conservative"`, otherwise `"calibrated"`. The `conservative` category exists but never appears in this data — no method landed below 0.67. `TOX-raw` is the only method that crosses the 1.5 boundary, at NSCLC stage II (1.49, just under).

---

## Table 1 — Kidney cancer

| method | stage | n_samples | FPR_0.05 | FPR_0.01 | median_p | ks_D | ad_A2 | hits_FDR05 | inflation | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| edgeR | I | 272 | 0.1308 | 0.0606 | 0.406 | 0.1074 | 913.1 | 0.0343 | 2.62 | ANTI-CONSERVATIVE |
| limma | I | 272 | 0.0441 | 0.0083 | 0.508 | 0.0656 | 247.5 | 0.0000 | 0.88 | calibrated |
| **TOX-log** | I | 272 | 0.0471 | 0.0096 | 0.484 | 0.1090 | 568.9 | 0.0001 | 0.94 | calibrated |
| TOX-raw | I | 272 | 0.0837 | 0.0133 | 0.442 | 0.0960 | 524.9 | 0.0000 | 1.67 | ANTI-CONSERVATIVE |
| edgeR | II | 60 | 0.1315 | 0.0517 | 0.386 | 0.1224 | 1050.4 | 0.0165 | 2.63 | ANTI-CONSERVATIVE |
| limma | II | 60 | 0.0401 | 0.0070 | 0.513 | 0.0798 | 356.3 | 0.0001 | 0.80 | calibrated |
| **TOX-log** | II | 60 | 0.0486 | 0.0115 | 0.472 | 0.1055 | 486.7 | 0.0010 | 0.97 | calibrated |
| TOX-raw | II | 60 | 0.0782 | 0.0096 | 0.424 | 0.1016 | 519.4 | 0.0000 | 1.56 | ANTI-CONSERVATIVE |
| edgeR | III | 123 | 0.1227 | 0.0493 | 0.400 | 0.1087 | 806.1 | 0.0182 | 2.45 | ANTI-CONSERVATIVE |
| limma | III | 123 | 0.0507 | 0.0100 | 0.490 | 0.0626 | 235.8 | 0.0000 | 1.01 | calibrated |
| **TOX-log** | III | 123 | 0.0552 | 0.0135 | 0.487 | 0.1060 | 528.8 | 0.0003 | 1.10 | calibrated |
| TOX-raw | III | 123 | 0.0846 | 0.0130 | 0.435 | 0.0916 | 507.8 | 0.0000 | 1.69 | ANTI-CONSERVATIVE |
| edgeR | IV | 83 | 0.1205 | 0.0471 | 0.403 | 0.1061 | 750.6 | 0.0153 | 2.41 | ANTI-CONSERVATIVE |
| limma | IV | 83 | 0.0513 | 0.0096 | 0.487 | 0.0560 | 190.5 | 0.0000 | 1.03 | calibrated |
| **TOX-log** | IV | 83 | 0.0549 | 0.0130 | 0.487 | 0.1079 | 534.9 | 0.0003 | 1.10 | calibrated |
| TOX-raw | IV | 83 | 0.0866 | 0.0119 | 0.411 | 0.1084 | 639.7 | 0.0000 | 1.73 | ANTI-CONSERVATIVE |

## Table 2 — Non-small-cell lung cancer (NSCLC)

| method | stage | n_samples | FPR_0.05 | FPR_0.01 | median_p | ks_D | ad_A2 | hits_FDR05 | inflation | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| edgeR | I | 295 | 0.1304 | 0.0560 | 0.398 | 0.1122 | 889.5 | 0.0250 | 2.61 | ANTI-CONSERVATIVE |
| limma | I | 295 | 0.0508 | 0.0105 | 0.495 | 0.0514 | 185.2 | 0.0001 | 1.02 | calibrated |
| **TOX-log** | I | 295 | 0.0528 | 0.0133 | 0.491 | 0.0630 | 184.2 | 0.0004 | 1.06 | calibrated |
| TOX-raw | I | 295 | 0.0814 | 0.0100 | 0.379 | 0.1243 | 608.9 | 0.0000 | 1.63 | ANTI-CONSERVATIVE |
| edgeR | II | 126 | 0.1216 | 0.0473 | 0.403 | 0.1059 | 741.2 | 0.0147 | 2.43 | ANTI-CONSERVATIVE |
| limma | II | 126 | 0.0524 | 0.0109 | 0.492 | 0.0506 | 161.8 | 0.0000 | 1.05 | calibrated |
| **TOX-log** | II | 126 | 0.0505 | 0.0128 | 0.505 | 0.0603 | 169.5 | 0.0004 | 1.01 | calibrated |
| TOX-raw | II | 126 | 0.0746 | 0.0081 | 0.393 | 0.1108 | 485.2 | 0.0000 | **1.49** | **calibrated** |
| edgeR | III | 85 | 0.1193 | 0.0464 | 0.405 | 0.1030 | 679.0 | 0.0134 | 2.39 | ANTI-CONSERVATIVE |
| limma | III | 85 | 0.0503 | 0.0102 | 0.497 | 0.0480 | 150.3 | 0.0000 | 1.01 | calibrated |
| **TOX-log** | III | 85 | 0.0559 | 0.0145 | 0.493 | 0.0703 | 247.3 | 0.0005 | 1.12 | calibrated |
| TOX-raw | III | 85 | 0.0804 | 0.0097 | 0.405 | 0.1011 | 454.6 | 0.0000 | 1.61 | ANTI-CONSERVATIVE |
| edgeR | IV | 26 | 0.1095 | 0.0350 | 0.406 | 0.1006 | 584.2 | 0.0004 | 2.19 | ANTI-CONSERVATIVE |
| limma | IV | 26 | 0.0492 | 0.0092 | 0.493 | 0.0471 | 128.9 | 0.0000 | 0.98 | calibrated |
| **TOX-log** | IV | 26 | 0.0611 | 0.0183 | 0.509 | 0.0511 | 116.6 | 0.0015 | 1.22 | calibrated |
| TOX-raw | IV | 26 | 0.0843 | 0.0074 | 0.421 | 0.0866 | 324.8 | 0.0000 | 1.69 | ANTI-CONSERVATIVE |

## Table 3 — Stomach cancer

| method | stage | n_samples | FPR_0.05 | FPR_0.01 | median_p | ks_D | ad_A2 | hits_FDR05 | inflation | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| edgeR | I | 59 | 0.1258 | 0.0516 | 0.399 | 0.1097 | 814.8 | 0.0180 | 2.52 | ANTI-CONSERVATIVE |
| limma | I | 59 | 0.0494 | 0.0095 | 0.493 | 0.0471 | 147.2 | 0.0000 | 0.99 | calibrated |
| **TOX-log** | I | 59 | 0.0518 | 0.0129 | 0.514 | 0.0766 | 239.3 | 0.0004 | 1.04 | calibrated |
| TOX-raw | I | 59 | 0.0806 | 0.0088 | 0.414 | 0.0944 | 429.9 | 0.0000 | 1.61 | ANTI-CONSERVATIVE |
| edgeR | II | 122 | 0.1234 | 0.0512 | 0.405 | 0.1050 | 767.8 | 0.0202 | 2.47 | ANTI-CONSERVATIVE |
| limma | II | 122 | 0.0526 | 0.0104 | 0.489 | 0.0515 | 177.4 | 0.0000 | 1.05 | calibrated |
| **TOX-log** | II | 122 | 0.0588 | 0.0161 | 0.505 | 0.0848 | 353.9 | 0.0008 | 1.18 | calibrated |
| TOX-raw | II | 122 | 0.0855 | 0.0119 | 0.405 | 0.1038 | 561.8 | 0.0000 | 1.71 | ANTI-CONSERVATIVE |
| edgeR | III | 168 | 0.1385 | 0.0624 | 0.390 | 0.1209 | 1011.1 | 0.0311 | 2.77 | ANTI-CONSERVATIVE |
| limma | III | 168 | 0.0499 | 0.0095 | 0.493 | 0.0530 | 163.8 | 0.0000 | 1.00 | calibrated |
| **TOX-log** | III | 168 | 0.0513 | 0.0125 | 0.507 | 0.0792 | 268.9 | 0.0002 | 1.03 | calibrated |
| TOX-raw | III | 168 | 0.0865 | 0.0117 | 0.395 | 0.1125 | 566.4 | 0.0000 | 1.73 | ANTI-CONSERVATIVE |
| edgeR | IV | 39 | 0.1109 | 0.0375 | 0.409 | 0.0986 | 617.4 | 0.0026 | 2.22 | ANTI-CONSERVATIVE |
| limma | IV | 39 | 0.0492 | 0.0091 | 0.492 | 0.0515 | 165.1 | 0.0000 | 0.98 | calibrated |
| **TOX-log** | IV | 39 | 0.0532 | 0.0136 | 0.516 | 0.0938 | 399.2 | 0.0004 | 1.06 | calibrated |
| TOX-raw | IV | 39 | 0.0762 | 0.0072 | 0.430 | 0.0858 | 411.7 | 0.0000 | 1.52 | ANTI-CONSERVATIVE |

---

## Headline summary

| method | verdict record (of 12) | inflation range | mean inflation |
|---|---|---|---|
| **edgeR** | 12/12 ANTI-CONSERVATIVE | 2.19 – 2.77 | **2.48** |
| **limma** | 12/12 calibrated | 0.80 – 1.05 | **0.98** |
| **TOX-log** | 12/12 calibrated | 0.94 – 1.22 | **1.07** |
| **TOX-raw** | 11/12 ANTI-CONSERVATIVE | 1.49 – 1.73 | **1.64** |

---

## Findings

**edgeR is badly miscalibrated on real data, across every single cancer type and stage — and this is the opposite of what the simulated benchmark showed.** edgeR was well-calibrated on `nb`, `lnpois`, and even mostly on `bimodal` in the synthetic tests. On TCGA it over-calls by a consistent 2.2–2.8× at nominal 5% in all 12 groups, with `hits_FDR05` up to 0.034 — meaningfully more BH-level false discoveries than any other method. This is not sampling noise: with `FPR_0.05_sd = 0.0421` on 100 runs, the standard error of the reported mean is `0.0421/√100 ≈ 0.0042`, so the kidney-stage-I deviation from nominal (0.1308 − 0.05 = 0.0808) is about **19 standard errors** from zero. Real cancer-cohort expression violates something edgeR's negative-binomial GLM relies on, in a way none of the synthetic arms captured — most plausibly correlation structure between genes (co-regulation, pathway effects) that a per-gene independent dispersion model doesn't see, though I don't have direct evidence of that specific mechanism from this table alone.

**limma is the best-calibrated method in every single group, on the continuous statistics as well as the pass/fail verdict.** `ks_D` for limma is smaller than `TOX-log`'s in **all twelve** comparisons — verified directly, not a near-tie. So while both are labelled "calibrated," limma's p-values are measurably closer to uniform throughout.

**`TOX-log` passes the practical calibration check in all 12 groups** — inflation 0.94–1.22, tightly banded, no verdict failures. Given the two synthetic-data findings so far — TOX matching edgeR on `nb`/`lnpois` and beating both edgeR and DESeq2 on `bimodal` — this is the third and most important confirmation, because it's the one dataset with no simulator assumptions at all.

**`TOX-raw` fails calibration almost everywhere on real data — this is the opposite of the `bimodal` synthetic result, where `norm_method = 0` was the *better*-calibrated mode.** 11 of 12 groups are ANTI-CONSERVATIVE at 1.49–1.73× inflation, with the sole exception (NSCLC stage II, 1.49) sitting just under the confirmed 1.5 verdict boundary. Real TPM is far more right-skewed than the synthetic `bimodal` arm's two-point mixture; the log transform likely does real work here that the synthetic multimodal test didn't need it to do. This is a case where the simulated and real-data results actively disagree, and the real data should be weighted more heavily for choosing a default `norm_method`.

**`TOX-raw`'s inflation is real at the raw-p level but essentially never reaches actual BH discoveries.** `hits_FDR05 = 0.0000` in all 12 groups despite FPR inflation up to 1.73×, while `edgeR` — with higher raw inflation (2.2–2.8×) — produces real BH hits in every single group. Consistent with the `ad_A2` values: `TOX-raw`'s tail statistic (325–640) is well below `edgeR`'s (584–1050) even where their raw `FPR_0.05` numbers are closer. `TOX-raw`'s miscalibration looks like a broad, shallow shift across the p-value distribution rather than a concentrated spike in the extreme tail — the same regime BH actually reads. Practically: `TOX-raw` would rarely hand a user an outright false discovery under BH in a homogeneous cohort, even though its raw p-values are measurably non-uniform. This doesn't make `TOX-raw` safe to use — the raw inflation is a genuine defect — but it bounds the practical damage differently than the raw number alone suggests.

**Kidney cancer is the hardest cancer type for `TOX-log`'s tail behaviour, and the effect is disproportionate to `TOX-log` specifically.** Mean `ad_A2` by cancer type:

| method | kidney | NSCLC | stomach | kidney / NSCLC ratio |
|---|---|---|---|---|
| TOX-log | 529.8 | 179.4 | 315.3 | **3.0×** |
| limma | 257.5 | 156.6 | 163.4 | 1.6× |
| edgeR | 880.0 | 723.5 | 802.8 | 1.2× |

All three methods find kidney harder than NSCLC, so this is partly a property of the kidney data itself (plausibly more within-stage biological heterogeneity, e.g. subtype mixture within a TCGA stage label). But the ratio is roughly double for `TOX-log` versus limma, and the aggregate `FPR_0.05` check doesn't show it — kidney's inflation values (0.94–1.10) look as good as anywhere else. This is a real gap between what the pass/fail verdict reports and what the tail-sensitive statistic reports, and it's a concrete argument for reporting `ad_A2` alongside `FPR_0.05` rather than the verdict alone.

**The gene-count gap has a confirmed mechanism, and it sharpens the concern rather than resolving it.** `edgeR`/`limma` filter with `filterByExpr()`, recomputed **fresh on each random half's counts** every split. TOX filters with `compute_gene_keep_mask()`, computed **once per cancer type from the full cohort**, before the stage/split loop even begins — which is why TOX's gene count is constant across all four stages of a cancer type while edgeR/limma's drifts slightly split-to-split. This has one good property: TOX's filter can't be responding to the particular random partition, so there's no split-specific cherry-picking. But it also means the two filters answer different questions — "expressible given this design" (edgeR/limma, standard and well-validated) vs. "meets TOX's threshold given the whole cohort" (TOX-specific, criteria not reviewed here). Whether the excluded ~15–20% skews toward the genes hardest to calibrate is still open, and I'd want `compute_gene_keep_mask()`'s actual criteria before treating TOX's real-data calibration advantage as settled rather than partly a product of an easier gene set.

**A related, unverified risk: silently dropped genes could be diluting `FPR_0.05`'s denominator.** In `run_tox_once()`, p-values carrying the not-computed sentinel are removed with no count kept: `p[p < 0 | p > 1] <- NA; p[!is.na(p)]`. If the `n_genes` reported in the summary table is taken **before** this drop rather than after, `FPR_0.05` for TOX would be computed against a numerator from fewer genes than its own denominator implies — which would bias `FPR_0.05` **downward** and make TOX look artificially better calibrated. Worth a direct check: `length(res$pvalues_own)` after the drop vs. `ncol(m_tpm)` for a single run, and confirming which count feeds `n_genes` in `make_row()`.


# Open ToDo's:
- Check wether a [5,95] percentile cap improves the performance of the raw model
- Implement multi-dimenionality handling
- Run on all 8 cancer cohorts
- Run on mammal dataset
- 