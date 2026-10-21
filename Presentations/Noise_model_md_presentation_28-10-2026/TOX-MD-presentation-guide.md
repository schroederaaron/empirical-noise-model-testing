# Presentation guide: Evaluating the multi-dimensional TOX noise model

Student project kick-off, scientific presentation (30–50 min). Source document for a LaTeX/Beamer first draft.

Presentation date: **Wednesday, 28 October 2026**.

Status of this document: written from the Tensor-Omics repository (branch `145-empirical-noise-model`, commit `613a36d`, 2026-10-19, not merged into `main`) and from the project discussion. Items marked **[VERIFY]** could not be checked when this was written. Those that concern repository contents can be checked by the drafting assistant, who has both repositories locally; the others need the presenter. Items marked **[TODO]** need an input from the presenter.

---

## 0. Instructions for the drafting assistant (Claude Code)

1. Build a Beamer deck (16:9) following Section 4, one frame per slide, in the order given. Put the speaker notes in `\note{...}`.
2. Keep slide text short: at most 6 lines per slide, no full sentences where a keyword works. Put the explanation into the speaker notes.
3. **Do not invent numbers, results, citations or dataset names.** Where a number is needed and not given here, insert `\todo{...}` (a visible red marker). Add references only after the presenter confirms them.
4. Do not call the method "distribution-free" or "assumption-free". Say "no parametric count distribution; assumes genes with similar group means have exchangeable within-group residuals".
5. Draw diagrams in TikZ. Create the plots listed in Section 5 as placeholder figures (clearly labelled) unless real data is provided.
6. Use a colour-blind-safe palette. Use the same colour for the same object throughout: case, control, null, observed.
7. Backup slides go after `\appendix` and are not counted in the time.
8. Build with `latexmk -pdf`. Keep figures in `figures/`.
9. Both repositories are available locally: `Tensor-Omics` (branch `145-empirical-noise-model`) and `empirical-noise-model-testing`. Before writing, check every **[VERIFY]** that concerns repository content against them (paths, function names, argument names, defaults) and correct this document where it is wrong.
10. The presenter will provide a link to Issue #185 (the multi-dimensional model). Read it and reconcile slide 11 with it. Report discrepancies to the presenter instead of silently choosing one version.
11. Read the testing repository for the TCGA scripts and current results, but do not put result numbers on slides without the presenter's confirmation.
12. Default parameter values are **subject to change** (the presenter will decide them later). Show them only as "current default, may change", or avoid them.

---

## 1. Requirements from the organizers

The organizers asked for 30–50 minutes covering the six points below. Audience: students with basic knowledge of statistics, data science and machine learning, but only rudimentary life-science knowledge. The biology must be introduced.

| # | Requirement | Slides |
|---|---|---|
| 1 | Scientific question | 3 |
| 2 | Material (biological material and data) | 4, 5, 15, 16 |
| 3 | Methods the students should apply | 6–14, 17 |
| 4 | Data availability (available vs expected) | 15, 16 |
| 5 | First steps (first two weeks) | 19 |
| 6 | Outlook (coming months) | 20 |

Context that changed during planning: the project is supervised for a **team of students**, so work can be split into parallel tracks. The proposal text has already been submitted. It lists TCGA as available and the self-assembled multi-tissue dataset as optional; in this talk the assembly is presented as **highly recommended**: not a strict requirement, but expected as part of the students' work, and there is enough time in the semester for it. The number of students is not yet known, so the split into tracks is decided on the go.

---

## 2. Message of the talk

One sentence: *We have a statistical test that estimates its own noise distribution from the data instead of assuming one; your project is to find out when it can be trusted and how it compares with the standard tools, now for genes measured in several tissues at once.*

Three take-homes for the students:
1. Why the noise model matters for every differential-expression result.
2. What is new in TOX and where its assumptions can fail.
3. What they will build: a calibration and power evaluation, in parallel tracks.

Tone: honest about open questions. Students will find the weak spots anyway, and naming them makes the project credible.

---

## 3. Timing

Target 40 minutes: about 36 minutes of content plus a buffer for questions.

| Block | Slides | Min |
|---|---|---|
| Opening | 1–3 | 2.5 |
| Biology and statistics primers | 4–6 | 7 |
| Existing tools and TOX (1-D) | 7–9 | 6.5 |
| Multi-dimensional model | 10–12 | 6.5 |
| Evaluation design | 13–14 | 3.5 |
| Data | 15–16 | 3 |
| Students' work | 17–20 | 6.5 |
| Close | 21 | 0.5 |

**To shorten to 30 min:** merge slides 4 and 5, move slide 9 to the backup, drop slide 14, merge 15 and 16.
**To extend to 50 min:** add backup slides B1–B3 to the main line and include a short live demo of the R call from Section 6.

---

## 4. Slide-by-slide specification

Each entry: purpose, on-slide content, speaker notes, visual, source.

### Slide 1: Title (0.5 min)
- Title: "Evaluating a multi-dimensional noise model for differential expression in multi-tissue RNA-seq data".
- Names, affiliation, date (28 October 2026). **[TODO: affiliation and team name]**

### Slide 2: Agenda (0.5 min)
- Six items matching the organizers' list in Section 1: question, material, methods, data, first steps, outlook.

### Slide 3: Scientific question (1.5 min)
- On slide: "Does the multi-dimensional TOX model control false positives, and how does its power compare with limma, edgeR and DESeq2?"
- Notes: two parts: calibration (are p-values honest?) and power (does it find true effects?). Result for the students: a reproducible evaluation showing where the model is and is not reliable.
- Source: proposal text.

### Slide 4: Biology primer, expression and RNA-seq (2.5 min)
- On slide: gene → RNA → protein; RNA-seq gives one read count per gene per sample; replicates; library size differs between samples.
- Visual: simple pipeline diagram: tissue sample → sequencing → count matrix (genes × samples).
- Notes: counts are not directly comparable between samples because of different sequencing depth, hence normalization (TMM and median-of-ratios for counts; TPM as a within-sample scale). Replicates are independent biological samples of the same condition and give an estimate of natural variability. Photosynthesis genes are high in leaf and near zero in root: expression is tissue specific.
- Keep it to what a data-science student needs; do not go into library preparation.

### Slide 5: Case vs control and tissues (1.5 min)
- On slide: contrast case vs control; effect = log2 fold change; one gene across tissues = a vector.
- Visual: a gene's log2FC in three tissues as a small bar plot, one tissue responding and two not.
- Notes: this is why a multi-dimensional test is needed: the question is whether a gene responds anywhere, not in each tissue separately.

### Slide 6: Statistics primer (3 min)
- On slide: null distribution, p-value, calibration, FDR, power.
- Visual: two small plots: (a) histogram of p-values under the null, flat when calibrated; (b) observed FDR vs nominal.
- Notes: under the null, p-values should be uniform. Calibration means the false-positive rate at a threshold matches that threshold. Benjamini–Hochberg controls the FDR, the expected share of false discoveries among the calls. Power is the share of true effects found (TPR). A method can look powerful because it is anti-conservative, so power is only meaningful together with calibration.

### Slide 7: Existing tools (2 min)
- On slide: table: edgeR, DESeq2: negative binomial model with estimated dispersion; limma: linear model on transformed values with empirical Bayes variance moderation.
- Notes: each one assumes a noise model. When the assumption fails, error control can suffer. Do not overstate: these tools are well validated; the question is the multi-tissue setting and non-parametric assumptions.
- Multi-tissue note: there is no single obvious reference method for a joint test across tissues. This is deliberately left open: finding and justifying a comparison is part of the students' task, to be researched and discussed in the meetings. Present candidate ideas, not recommendations:
  1. Omnibus tests in a tissue-by-condition design (moderated F-test in limma, quasi-likelihood F-test in edgeR, likelihood-ratio test in DESeq2). They ask "any effect in any tissue", whereas TOX uses the norm of the effect vector.
  2. Combining per-tissue p-values.
  3. The presenter's suggestions `limma::fry()` and `edgeR::diffSpliceDGE()` (note: `diffSpliceDGE` is an edgeR function; limma's counterpart is `diffSplice`). They answer related but different questions (gene-set testing; differential effect across features within a gene), so the students must work out what each one tests. `fry()` is a gene-set test with set-level p-values only (no per-gene effects); `diffSpliceDGE()` tests relative exon-versus-gene effects. Assessing whether these approaches and comparisons are feasible per gene across tissues, including their limitations, is a student task.
  - Tell the students to look for the limitations of each approach: there will be no 1:1 comparability.

### Slide 8: TOX idea, 1-D (2.5 min)
- On slide: estimate the noise empirically; match genes by mean expression; pool their residuals; read off the null.
- Visual: scatter plot of replicate values against gene mean, with a window around the target gene's mean highlighted.
- Notes: for a gene, find neighbor genes with similar mean expression, collect their within-group residuals (replicate minus group mean), and use these as the empirical noise distribution. No negative-binomial assumption. The assumption is that genes with similar means share comparable noise.
- Source: `misc/Tensor_Omics_Methods/Empirical_Noise-Based_Significance_of_Distances.tex`.

### Slide 9: 1-D mechanics (2 min) *(move to backup in the 30-min version)*
- On slide: five steps: group means → residual pool (adaptive kNN, counted in genes) → scale by 1/√n → null for the mean difference → p-value from the tail.
- Notes: the residuals are scaled by 1/√n of their own group so that a single residual lands on the scale of a mean difference. The neighborhood grows outward in sorted-mean order until k_start genes, then adaptively until the pool changes too much (tau) or k_max is reached. Residuals carry a Bessel correction.
- Source: the 1/√n scaling and the exact pairwise tail count are in `src/tox/tox_noise_model_exact.F90` (`tox_noise_model.F90` bootstraps the means instead); the shared helpers `gather_residuals_helper` and `prepare_sorted_data_helper` live in `src/tox/tox_noise_model.F90`.

### Slide 10: Multi-tissue problem (1.5 min)
- On slide: a gene is one vector over d tissues (axes); axes treated as independent.
- Visual: gene vectors in a 3-axis space; case and control points.
- Notes: independence means no sample is shared between axes (no plant, patient or batch used in two tissues). This is metadata, not something the numbers reveal, so the R wrapper checks that no sample identifier appears on more than one axis (`check_axis_samples`). That check cannot detect a shared individual under different identifiers. In the TCGA design the axes are cancer cohorts (e.g. TCGA-KIRC, TCGA-LUAD, TCGA-STAD) and different cohorts contain different patients, which fits the independence assumption.

### Slide 11: Multi-dimensional model (3 min)
- On slide: per axis a: β̂ₐ = mean_case − mean_ctrl; D = √Σβ̂ₐ² (effect size, log2FC units with the default `norm_method = 1`); D_std = √Σ β̂ₐ²/Var_null(a) (test statistic, unitless).
- Reconciled with Issue #185 (presenter, 20/10/2026): the issue is the original spec, the code is current. The test statistic is D_std and D is the effect size; β̂ is computed internally from the replicates (supplying β̂ from edgeR, limma or DESeq2 is planned, but not supported yet). The direction of β̂ is not tested: that needs an angle-based noise model, a separate construction, and is left out.
- Visual (TikZ): for each axis a separate neighbor pool for case and control; one residual drawn from each; differences per axis; combined into one score. See Section 5, figure F4.
- Notes:
  - Neighbors are found separately per axis and per group, in 1-D, on that group's own mean. It is **not** a joint nearest-neighbor search in d dimensions. This is deliberate: a gene high in leaf and near zero in root has neighbors at a similar level in each tissue separately, whereas a joint search would find genes that match poorly on at least one axis.
  - Null draw: one case residual and one control residual per axis, both scaled by 1/√n; per-axis difference; weighted sum of squares. The joint null is the product of the per-axis nulls.
  - With independent axes the null covariance is diagonal, so D_std is the full Mahalanobis statistic.
  - p = (number of null draws ≥ observed + 1)/(m + 1). Three computation paths: exact enumeration (mostly d = 1), systematic stride sampling (default), random sampling.
- Source: `src/tox/tox_noise_model_md.F90`, `compute_noise_pvalue_pipeline_md_helper`, `compute_pvalue_md_helper`.

### Slide 12: Assumptions, open questions and what to look for (2 min)
- On slide:
  1. Genes with similar group means have exchangeable within-group residuals.
  2. Axes share no samples.
  3. Open: noise scale shared across tissues is not captured.
  4. Open: sensitivity to k, tau, normalization.
- Notes for 3: neighbors are matched on the mean only. If a gene has an overdispersion that is shared across tissues, its per-axis effects scale together, whereas the null draws scale each axis independently. This could make the p-values anti-conservative, more so with more axes. This is a hypothesis to test, not a result. It is shown deliberately as an outlook for the students: their task is to find out where problems occur (error control, difficulties with parameters, data issues). Present it as "what could go wrong and what you should look for", not as a known flaw.
- Notes for 4: the students explore parameters on a pre-specified grid and report all settings, not only the best one.

### Slide 13: Evaluation design (2.5 min)
- On slide: random split of a large dataset into two groups; (a) no effect → calibration; (b) known effect injected into one group → power.
- Visual: sample matrix split into halves, with a "known effect" arrow on one half.
- Notes: real data has no known truth, so truth is created. A random split gives a mock null. Adding a known log2FC to selected genes in one half gives a ground truth. The same approach is already used with TCGA.
- Pitfalls worth saying out loud: (1) repeated splits of one finite dataset are not independent replicates, so variability across splits understates uncertainty; (2) samples within a tissue must be exchangeable, otherwise the split inherits the heterogeneity; (3) the injected effect should be realistic (see the existing scripts); (4) average observed FDR over all rounds, counting rounds without calls as zero false discoveries.

### Slide 14: Metrics and plots (1 min) *(drop in the 30-min version)*
- On slide: QQ plot of p-values, false-positive rate at 0.001/0.01/0.05, observed vs nominal FDR, TPR at observed FDR.
- Visual: the existing TPR-versus-achieved-FDR figure from the proposal, if the presenter wants to show the scalar results. The proposal figure itself was not available; the draft uses the current scalar-model plots from the testing repository instead (Part R: observed FDR at nominal 0.05 and TPR at achieved FDR 0.05, run of 19 Oct 2026) on an added slide 14b, plus backups B7/B8. All are scalar (1-D) results; the presenter confirms them before showing.

### Slide 15: Material and data (1.5 min)
- On slide: available: TCGA RNA-seq, healthy samples from several cohorts as axes (e.g. TCGA-KIRC, TCGA-LUAD, TCGA-STAD), existing scripts; highly recommended: a self-assembled multi-tissue dataset.
- Notes: TCGA is the guaranteed core and is available now. Assembling a new multi-tissue dataset is highly recommended: not a strict requirement, but expected as part of the students' work. In the TCGA design each cohort is one axis and only the healthy (normal-tissue) samples are used; one sentence is enough, no detail needed. A caveat worth a short note: normal samples in TCGA are often taken adjacent to tumors and may not be perfectly healthy, and cohorts differ in how many normals they have. In the existing null-calibration run the healthy groups held 72 (TCGA-KIRC), 59 (TCGA-LUAD) and 36 (TCGA-STAD) samples. How the normal files were built is not needed for the talk.
- Not recommended as the main data: the multi-tissue mammal dataset in `material/` (`kallisto_sex_data_<Tissue>.tsv`). It has few replicates (6 samples per tissue, 8 in testis, so few per group once split), so calibration by random splits is hard, and it contains real signal (sex differences), so a split is not a fully known null. It can still be used for exploration.

### Slide 16: Criteria for a new dataset (1.5 min)
- On slide: checklist: replicates per tissue enough to split (roughly at least 6; fewer leaves little to split); independent samples across tissues; homogeneous samples within a tissue (same condition, comparable study); count data available; if studies are combined, the same protocols (sequencing platform and quantification).
- Notes: independence has to be read from the methods sections of the publications: check whether tissues come from the same individuals. Different studies per tissue can be acceptable because each axis is only compared within itself, provided the protocols are the same (same sequencing platform and the same quantification). Studies using different platforms or quantifiers should not be mixed. Some plant atlases use standardized extraction and processing protocols for exactly this reason (not essential for the talk). Within-tissue homogeneity must still be checked. Do not name a specific dataset unless verified. **[TODO: optional candidate list after checking]**

### Slide 17: What the students do, in tracks (2 min)
- On slide: Track A: dataset survey and assembly; Track B: TOX on TCGA, calibration scripts; Track C: reference tools and metrics, including ideas for a multi-tissue comparison (e.g. omnibus tests, `fry()`, `diffSpliceDGE()`; to be researched and discussed in the meetings); shared: onboarding, final report.
- Notes: tracks are a suggestion. The number of students is not yet known, so the split is decided on the go. Choosing and justifying the comparison method is left to the students; limitations are part of the task. Students use R and a provided Claude Code license; generated code must be checked by the student. No Fortran knowledge is needed: TOX is called through its R interface.

### Slide 18: Repository tour (1.5 min)
- On slide: the paths from Section 6 (README, manuscript, R interface, tests, docs).
- Notes: show where to start reading and where to run something. No live coding.

### Slide 19: First two weeks (2 min)
- On slide: (1) set up R, build or install TOX, run the example call; (2) read: RNA-seq normalization and the three tools; (3) explore TCGA: library sizes, filtering, mean-variance plot; (4) Track A: one-page comparison of candidate datasets against the criteria; (5) reproduce one existing calibration plot of the scalar model.
- Notes: the aim of the first two weeks is that everyone can run the pipeline and agrees on the dataset decision date.

### Slide 20: Outlook (1 min)
- On slide: months 1–2: onboarding and data; then calibration; then power; then parameter exploration; final report. About 180 h per student as a university guideline (6 ECTS), not a hard target.
- Notes: do as much as possible; calibration comes first and power as time allows.

### Slide 21: Take-home and questions (0.5 min)
- The three take-homes from Section 2.

---

### Backup slides (after `\appendix`)
- **B1:** full null construction for d axes with a worked example for d = 2 (pool sizes, number of atoms, weights). Sanity check: the expected value of the weighted null statistic per gene is the number of axes with positive weight.
- **B2:** the three p-value paths and the resolution: p-value floor 1/(M+1) on the sampled paths, default M = 20,000 (the enumerated path also uses (count+1)/(m+1), so its floor is 1/(∏ S_a + 1), tiny but not zero); why M is set by the number of genes under BH and not by d.
- **B3:** parameters and their current defaults, marked as subject to change (Section 6).
- **B4:** 1-D scaling argument: why residuals are scaled by 1/√n, and the caveat that this fixes the variance but not the shape of the distribution of a mean.
- **B5:** the joint-neighborhood alternative and why it was not chosen (dimensionality, tissue-specific genes).
- **B6:** glossary (Section 8).

---

## 5. Figures to produce

| ID | Content | How |
|---|---|---|
| F1 | RNA-seq pipeline: sample → counts matrix | TikZ |
| F2 | One gene's log2FC across tissues | R/ggplot or TikZ, illustrative data, labelled as such |
| F3 | Mean–residual scatter with the neighborhood window | R/ggplot from a real dataset, or schematic |
| F4 | MD null: d axes, case and control pools, one draw per axis, weighted sum | TikZ |
| F5 | Random-split design with injected effect | TikZ |
| F6 | p-value histogram (flat vs anti-conservative), observed vs nominal FDR | R/ggplot with simulated illustrations, clearly labelled as illustrative |
| F7 | Parallel tracks timeline, first two weeks | TikZ |

Existing figures in the repo: `misc/figures/tox_methods.svg`, `tox_methods_1a.png`, `tox_methods_1b.png`. Checked: they show the general Tensor-Omics geometry (shifts, angles, fields), not the noise model, so they are not used. Existing result plots used instead: `TCGA_test/results/power_out/grid_R_*.png`, `fdr_tpr_curves.png`, `TCGA_test/results/null_calibration/plots/null_fpr_inflation.png`.

---

## 6. Where to find things in the repository

Repository: `https://github.com/asishallab-group/Tensor-Omics`. Noise-model work is on branch `145-empirical-noise-model` (commit `613a36d` at time of writing). The branch is **not merged into `main`**.

| What | Path |
|---|---|
| Project overview, folder structure, build and test instructions | `README.md`, `concise_project_info.md` |
| Methods manuscript | `misc/Tensor_Omics_Methods.pdf`, chapters in `misc/Tensor_Omics_Methods/*.tex` |
| Chapter on the empirical noise model | `misc/Tensor_Omics_Methods/Empirical_Noise-Based_Significance_of_Distances.tex` |
| Manual | web documentation in `doc/` (there is no `misc/tox_manual.pdf`) |
| 1-D noise model (Fortran) | `src/tox/tox_noise_model.F90` (bootstrap of means; `null_method` pooled or gene-blocked), `tox_noise_model_exact.F90` (exact, 1/√n scaling) |
| Multi-dimensional noise model (Fortran) | `src/tox/tox_noise_model_md.F90` |
| Variance-aware variant (work in progress, not for this talk) | `src/tox/tox_noise_model_variance.F90` |
| R interface | `rcpp/tensoromics_functions.R`: `tox_compute_noise_pvalues_pipeline`, `..._exact`, `..._variance`, `..._md` |
| Fortran unit tests | `test/mod_test_noise_model.f90`, `mod_test_noise_model_md.f90`, `mod_test_noise_model_variance.f90` |
| Build and test scripts | `build.sh`, `test_runner.sh`, `run_all_tests.sh`; compiler container: `misc/gfortran.docker` |
| Example data | `material/` (e.g. `kallisto_sex_data_<Tissue>.tsv`, `tpm.tsv`); the mammal data is not recommended for calibration, see slide 15 |
| TCGA test scripts, power-test script and current results | separate repository `empirical-noise-model-testing`, checked out as the submodule `experiments/Noise_Model_Test`: `TCGA_test/scripts/null_calibration.R`, `power_test.R`, `power_report.R`; `Simulated_data/scripts/calibration_test.R`, `calibration_test_md.R` (multi-dimensional; smoke run only, no results yet); results in `*/results/` |
| Issue describing the multi-dimensional model | Issue #185 in the Tensor-Omics repository (link provided by the presenter; read it and reconcile slide 11) |

Note: paths and function names were read directly from the repository. Line numbers are not given because they change.

### Facts verified from the code (safe to put on slides)

- The R function for the multi-dimensional model is `tox_compute_noise_pvalues_pipeline_md`. It takes **lists of per-axis replicate matrices** (samples × genes), one for case and one for control. Replicate counts may differ between axes; the gene columns must be identical and in the same order; row names must be sample IDs.
- Defaults (all subject to change): `norm_method = 1` (log2), `k_start = 20`, `k_step = 1`, `k_max = 50` (counted in genes, not residuals), `tau = 0.1`, `sampling_mode = 0` (systematic stride), `n_draws_max = 20000`, `enum_max_product = 1e5`, `n_exceed_target = 0` (off).
- `null_method` exists only for interface parity and must be 0 in the multi-dimensional function.
- Outputs per gene include: `pvalues_own`, `d_obs` (D), `d_std_obs` (D_std), `d_sq_null_mean`, `method_used`, `var_null` per axis, neighborhood sizes per axis. Genes that were not tested have `pvalues_own = -1`.
- A gene is dropped entirely if any axis has fewer than 10 residuals in its pool or a non-finite effect.
- `beta_centre`: 0 = none (default), 1 = subtract the per-axis median beta. It must be 0 if the data were already normalized on the R side, otherwise the data are centred twice.
- The manuscript chapter now describes the adaptive neighbourhood (`k_start`, `k_step`, `k_max`, `tau`, sorted means with binary search); the old k = 100 and k-d-tree text was corrected on 20/10/2026. On slides still say "k is a parameter" until the defaults are decided. On slides say "k is a parameter" and do not quote a default; the presenter will decide the defaults later.

### Minimal call to show on a slide or in a demo (illustrative; check argument names against the current file)

```r
res <- tox_compute_noise_pvalues_pipeline_md(
  case_replicates    = list(leaf = case_leaf, stem = case_stem, root = case_root),
  control_replicates = list(leaf = ctrl_leaf, stem = ctrl_stem, root = ctrl_root),
  max_pool_size      = 50000L
)
```

`max_pool_size` is deliberately left off the slides: per the presenter it is set high (about 50,000) and is never reached. The R function has no default for it, so it must be passed (checked against the local file).

---

## 7. Claims to avoid and results to verify

- Avoid: "distribution-free", "reliable regardless of the data's true distribution", "outperforms".
- Do not show any benchmark numbers against limma, edgeR or DESeq2 until they are taken from the current results in the local `empirical-noise-model-testing` repository and confirmed by the presenter.
- Show the scalar-model results as they are, including where TOX is worse than the established tools. Students will see this anyway.
- Do not present the LOESS/variance-aware variant as a result; it is under development.
- No literature citations unless verified by the presenter.

---

## 8. Likely questions and short answers

| Question | Answer |
|---|---|
| Why not simply use DESeq2? | It is a strong baseline and part of the comparison. The question is how a non-parametric empirical null behaves, especially for a joint test across tissues. |
| How are tissue-specific genes handled? | Neighbors are found per axis on that group's own mean, so a gene high in leaf and near zero in root gets suitable neighbors in each tissue. |
| What if a tissue has few replicates? | Each axis needs at least 2 replicates per side. Pools are built from many genes, so the residual pool is larger than the replicate number; the effect of small n on the null shape is part of the evaluation. |
| Why scale residuals by 1/√n? | The observed statistic is a difference of means; a single residual is at single-observation scale. Without scaling the null is too wide (conservative). |
| What does independence of axes mean in practice? | No sample, individual or batch appears in more than one axis. The code checks sample identifiers, not biology. |
| What is the difference between D and D_std? | D is an effect size in log2FC units; D_std divides each axis by its null variance and is the test statistic that drives the p-value. |
| Why a p-value floor? | The p-value is estimated from M draws, so its smallest value is 1/(M+1). With many genes and BH, a floor that is too high can prevent rejections. |
| Do the students need Fortran? | No. TOX is called from R. |
| Is ground truth available? | Not in real data. It is created by random splits (null) and injected effects (power). |
| What if we find that TOX is worse? | That is a valid result. The deliverable is an honest evaluation. |
| Which method is the reference for the multi-tissue test? | There is no obvious one-to-one reference. Finding candidates and understanding their limitations is part of the project; we will discuss it in the meetings. |

---

## 9. Decisions and remaining open points

Decided:
1. Dataset assembly is presented as **highly recommended**: not strict, but expected as part of the students' work.
2. The comparison method for the multi-tissue test is left open for the students, who research ideas (omnibus tests, `fry()`, `diffSpliceDGE()`, others) and discuss limitations with the presenter.
3. TCGA axes are cohorts (e.g. TCGA-KIRC, TCGA-LUAD, TCGA-STAD), using the healthy samples.
4. Team size and track split are decided on the go.
5. Different studies per tissue are acceptable if the protocols are the same (platform and quantifier); not essential for the talk.
6. The shared-noise-scale concern is shown on slide 12 as an outlook. The mammal dataset is not recommended for calibration.

7. Slide 11 follows the code, not Issue #185: D_std is the test statistic, D the effect size, β̂ computed internally; direction left out (angle-based, separate construction).
8. `norm_method = 1` (log2) is the default of the multi-dimensional R wrapper (changed 20/10/2026).
9. The feasibility of `fry()`, `diffSpliceDGE()` and other comparisons is a student task (slide 7).
10. How the TCGA normal files were built is not relevant for the talk.

Remaining open:
- Default parameter values other than `norm_method` (to be decided later; marked as subject to change).

---

## 10. Preparation plan (talk on Wednesday, 28 October 2026)

| Day | Task |
|---|---|
| Tue 20 – Wed 21 Oct | Give the drafting assistant this document, both local repositories and the Issue #185 link; generate the first draft; talk to the supervisor about the open points. |
| Thu 22 Oct | Check all **[VERIFY]** items against the repositories; produce figures F1–F7. |
| Fri 23 Oct | Revise content; run the minimal R call once on a small example; first timed run-through. |
| Sat 24 – Sun 25 Oct | Buffer. |
| Mon 26 Oct | Second timed run-through; trim against the table in Section 3. |
| Tue 27 Oct | Final pass: remove all `\todo` markers, check Section 7, prepare the backup slides. |
| Wed 28 Oct | Talk. |

---

## 11. Prompt for Claude Code

> Read `TOX-MD-presentation-guide.md`. Create a Beamer (16:9) presentation `talk.tex` following Section 4 slide by slide. Use Section 5 for the figures (TikZ for diagrams, placeholder figures labelled "illustrative" for plots), and Sections 0 and 7 for the rules on claims. Put speaker notes into `\note{}`. Put all unverified items into visible `\todo{}` markers. Add the backup slides after `\appendix`. Make sure the document compiles with `latexmk -pdf talk.tex` and report any warnings. Do not invent numbers, results or references.
