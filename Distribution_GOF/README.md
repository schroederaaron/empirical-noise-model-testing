# Distribution_GOF — omnibus goodness-of-fit for RNA-seq count distributions (issue #188)

Which count distribution describes RNA-seq replicate noise? Six per-gene models —
**Poisson, NB, ZINB, generalised Poisson, Poisson-lognormal, Conway–Maxwell–Poisson** —
are fitted to every filtered gene, and each is judged by the uniformity of its
randomised PIT (CvM `W²`, AD `A²`, KS `D`) against a **parametric bootstrap**, plus a
K-fold **held-out predictive log score**. Built to
`docs/claude_omnibus_gof_implementation.md` (23.09.2026).

**This README contains no results.** Every number about data lives in a file a
script wrote (see *Outputs*). Nothing here draws conclusions about RNA-seq noise.

## Running

Run from the directory that holds `external/docker_r_libs` (the Tensor-Omics root on
the cluster), in the `arch-gfortran` image. Every knob is in `config_gof.R`; override
with `--NAME=value` or the env var `GOF_NAME`.

```
Rscript Distribution_GOF/scripts/02_run_gof.R --dry-run                    # M0: print resolved config
Rscript Distribution_GOF/scripts/00_validate_families.R                    # G1-G4 (hard stop on G1-G3)
Rscript Distribution_GOF/scripts/02_run_gof.R --DATASET=ds.rds             # main analysis, one dataset
Rscript Distribution_GOF/scripts/01_calibration_sim.R --G5_REF_DIR=<02 output dir>   # G5
Rscript Distribution_GOF/scripts/03_report.R --RUN_DIR=<02 output dir>     # table + plots, no refitting
Rscript -e 'GOF_ROOT <- "Distribution_GOF"; source(file.path(GOF_ROOT, "R/setup.R")); testthat::test_dir(file.path(GOF_ROOT, "tests/testthat"))'
```

Datasets follow the contract of brief §3. For TCGA and simulated data use
`scripts/make_dataset.R` (it prints the integer audit, so the ROUNDING decision can be made
before the run):

```
Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=tcga --PROJECT=TCGA-KIRC --STAGE="Stage I"
Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=tcga --PROJECT=TCGA-KIRC --STAGE=healthy
Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=simulated --N_GENES=2000 --N_PER_GROUP=24
```

For Salmon/tximport or featureCounts input (which need a sample sheet) call
`build_from_tximport()` / `build_from_featurecounts()` in `R/io.R` and `saveRDS()` the result.

### Packages

`R/setup.R` self-installs into `external/docker_r_libs` (the repo's `load_or_install()`,
copied verbatim). Two build-time traps in the arch-gfortran image, both hit on
23.09.2026:

* `fs` (a testthat dependency) needs libuv headers — solved in `setup.R` by
  `USE_BUNDLED_LIBUV=1` (static libuv, no runtime dependency);
* `nloptr` (glmmTMB → lme4 → nloptr) needs **CMake at build time** to compile its
  bundled static nlopt. The image has none: install it once in the container that does
  the first install (`pacman -Sy --noconfirm cmake`). Afterwards it is not needed.

## Decisions (brief §12) — each a config value, printed and recorded in provenance

| id | config | default here | note |
|---|---|---|---|
| D1 | `ROUNDING` | `error` | `stochastic` / `floor` must be chosen explicitly; written into every output path |
| D2 | `OFFSET_MODE` | `tmm` | `tmm_length` for tximport input (see *Offsets*) |
| D3 | `OFFSETS_BOOT` | `reestimate` | |
| D4 | `N_GENES_BOOT` | 2000 | stratified by mean-expression decile |
| D5 | `RANK_BY` | `heldout` | the table always shows both ranks |
| D6 | `EXCLUDE_SAMPLES` | none | decide from QC, not from GOF results |
| D7 | glmmTMB version | whatever is installed | recorded in `provenance.txt` |

## Method notes and limits (read before interpreting anything)

* **Bootstrap null ignores gene–gene correlation and shared sample-level effects.**
  Genes are simulated independently given their fitted parameters. If such structure
  is present, `p_GOF` is anti-conservative. At RNA-seq N, `p_GOF = 1/(B+1)` for every
  model is the expected outcome, **not a finding**.
* **`excess = T_obs / median(T_b)` and `z` are descriptive**, cross-model comparisons
  that correct for each model's own in-sample estimation effect (Poisson has 0 extra
  parameters per gene, ZINB 2). They are **not a formal model-selection criterion**.
* **TCGA cohorts** are fitted as a single group (`~ 1`): unmodelled between-patient
  heterogeneity counts as lack of fit for **every** model.
* **Held-out score:** the held-out sample's offset comes from the full-data offsets — a
  small leakage. The gene-resampling interval for the difference to NB is labelled
  approximate: genes are not independent.
* **ZINB:** π is constant per gene and does not scale with library size (a modelling choice).
* **No dispersion shrinkage:** per-gene MLE, as the issue specifies. Shrinkage would be a
  different estimand — a possible follow-up, not implemented.
* **Bootstrap parallelism:** over replicates; genes run serially inside a worker. Each
  replicate's streams are `seed_for(BASE_SEED, "boot", model, b)` and
  `seed_for(BASE_SEED, "bootpit", model, b)` under L'Ecuyer-CMRG, so one replicate can be
  rerun alone (`run_one_replicate()`).
* **Observed bootstrap reference:** `T_obs` is computed on `G_BOOT` from a refit with
  offsets estimated on the `G_BOOT` counts — exactly what every replicate does — not from
  the `G_ALL` fit (whose offsets use all filtered genes).

## Deviations from the brief (and why)

1. **`tmm_length` offsets use `edgeR::DGEListFromTximport()`.** The brief says to copy
   the edgeR recipe from the installed tximport vignette. The installed vignette
   (tximport 1.40.0) no longer carries the manual `normMat`/`scaleOffset` recipe; it
   delegates to `DGEListFromTximport(txi)` (edgeR ≥ 4.10.0). That function stores
   `offset.prior = log(length) − rowMeans(log(length))`; after `normLibSizes()`,
   `getOffset()` returns `offset.prior + log(lib.size × norm.factors)`. This is what is used.
2. **PLN cdf is not "Σ wₖ ppois on the same adapted nodes".** Nodes adapted to the pmf
   peak span only ±~8τ in z; when τ is small (large y, moderate σ) the cdf integrand
   φ(z)·ppois(q,·) has most of its mass outside that span (y = 1000, σ = 0.5: nodes
   within |z| < 0.61 while F(999) ≈ 0.5 lives on the whole line). Instead
   `F(q) = Φ(z*) − ∫_{−∞}^{z*} φ(1−S) + ∫_{z*}^{∞} φS` with composite Gauss–Legendre on
   graded panels, and exact pmf summation in the deep lower tail (F < 1e-25). Verified
   against `integrate()`, the summed pmf and a fine trapezoid rule.
3. **`PLN_NODES` default is 60, not 30.** With 30 adaptive Gauss–Hermite nodes the pmf
   has a 2.5e-10 relative error at σ = 1, μ = 0.5 (against a trapezoid reference), which
   fails G1's 1e-10 increment check; 60 passes the whole G1 grid. This raises precision;
   no tolerance was changed.
4. **Own PLN and CMP fitters use analytic gradients.** With finite differences `nlminb`
   stopped with "singular convergence" on most genes (convergence ≠ 0 ⇒ not OK). The
   scores come from the same quadrature/window sums (PLN: `y − E[λ|y]`, `E[z(y−λ)|y]`;
   CMP: `(y−μ)μ/V` and `(y−μ)Cov(Y, log Y!)/V − log y! + E[log Y!]`), checked against
   finite differences.
5. **G6 / G1 increment tolerance has an absolute floor.** "`b − a` equal to
   `exp(logpmf(y))` to 1e-10 relative" cannot hold in double precision once `b ≈ 1`: `b − a`
   is then a difference of two numbers near 1 with rounding ~ε·b, for any family (it
   already fails for R's own `pnbinom`/`dnbinom`). Implemented as
   `|(b−a) − p| ≤ 1e-10·p + 64·ε·b`; the number of observations where the floor binds is
   reported (`g6_pit_asserts.csv`, `G1_pmf_cdf.csv`). The 1e-10 relative criterion is
   unchanged wherever it is representable.
6. **Interface additions:** `mu_fn` (means for held-out samples) and `moments` (true
   mean/variance for G2) on every family.
7. **`R/setup.R`** (not in the brief's file list) holds the single copy of the
   package bootstrap every script needs; **`scripts/make_dataset.R`** builds dataset files.
8. **G3, CMP cross-check: genes where the glmmTMB reference fails are excluded**
   (decision 24.09.2026). `glmmTMB::compois` (1.1.14) failed to converge or errored on some
   simulated genes where the own fitter converged near the truth. Such genes — reference
   `fit$convergence != 0` or an error — are excluded from the gate and listed in
   `results/validation/G3_excluded_reference_failures.csv`, counted in `verdict.csv`. A
   failure of the OWN fitter is never excluded, and a glmmTMB fit that reports convergence
   is always compared.

## Outputs (`results/<label>_round-<mode>_off-<mode>/`)

`input_audit.csv`, `preprocess.csv`, `fit_diag_gall_<m>.csv`, `params_gall_<m>.csv`,
`fit_failures.csv`, `analysis_genes.txt`, `g6_pit_asserts.csv`, `stats_gall_obs.csv`
(analysis set + each model's own OK set), `rerandomisation_gall.csv`,
`gboot_genes.txt`, `stats_gboot_obs.csv`, `strata_gboot_obs.csv`, `boot_<m>.rds`,
`boot_summary.csv`, `boot_strata_summary.csv`, `boot_timing.csv`, `heldout_<m>.rds`,
`heldout_summary.csv`, `heldout_pit_stats.csv`, `provenance.txt`, `sessionInfo.txt`,
`config_resolved.rds`; after `03_report.R`: `summary_table.{csv,md}`, `plots/`.
Validation: `results/validation/` (G1–G4), `results/G5_calibration[_SMOKE]/`.

## Validation status

See `results/validation/verdict.csv` and `results/G5_calibration*/G5_summary.csv` for
whatever has been run. Status at delivery is reported in the delivery note, not here.
