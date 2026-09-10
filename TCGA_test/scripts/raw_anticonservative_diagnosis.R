#!/usr/bin/env Rscript
# raw_anticonservative_diagnosis.R
# -----------------------------------------------------------------------------
# WHY the raw normalisation is anti-conservative -- a decomposition, not a sweep.
#
# THE OBSERVATION TO EXPLAIN (null_calibration_9414.out, 8 cancers x 5 cohorts
# x 9 kNN configs x 4 replicate arms, 100 splits each):
#
#   TOX-raw   FPR@.05 / .05 = 1.64 (median)   FPR@.01 / .01 = 1.20   median_p = 0.42
#   TOX-log   FPR@.05 / .05 = 1.05            FPR@.01 / .01 = 1.32   median_p = 0.51
#
# and, restricted to k20_50, as a function of replicates per group:
#
#   n<=12  13-25  26-50  51-100  >100        n<=12  13-25  26-50  51-100  >100
#   1.45   1.61   1.70   1.79    1.90        1.06   1.01   0.98   0.99    0.97
#   0.75   1.04   1.31   1.45    1.74        1.40   1.29   1.19   1.20    1.15
#   ^ TOX-raw infl@.05 / infl@.01            ^ TOX-log infl@.05 / infl@.01
#
# Three facts constrain the explanation before any new code runs:
#
#   (i)  k barely matters. Across all cells, TOX-raw's inflation moves from 1.69
#        (k8_30) to 1.61 (k20_50) -- 0.08 out of a 0.64 excess -- and k_max 30 ->
#        50 moves it by 0.005. A mis-set kNN parameter cannot be the cause.
#   (ii) At small n the raw model is anti-conservative at alpha = 0.05 (1.45) and
#        CONSERVATIVE at alpha = 0.01 (0.75) AT THE SAME TIME. No error in the
#        null's WIDTH can do that: a null that is too narrow inflates both levels,
#        and inflates the 1% level MORE. Fitting D_null = lambda * (standardised
#        t_nu) against a Gaussian observed gives lambda = 0.84-0.90 (null ~10-15%
#        too WIDE) with excess kurtosis 1.0-2.2. The defect is the null's SHAPE.
#   (iii) The inflation GROWS with n for raw (1.45 -> 1.90) and is flat for log.
#        A fixed scale error is n-invariant, because observed and null both shrink
#        as 1/sqrt(n). Something that gets worse with n is a shape mismatch: the
#        observed statistic is a difference of MEANS and Gaussianises with n (CLT),
#        while the null keeps the shape of INDIVIDUAL residuals however large n gets.
#
# THE HYPOTHESIS THIS SCRIPT TESTS
#   In raw (linear) space, genes at the SAME mean still have very different sds.
#   A mean-neighbourhood is mean-homogeneous by construction but NOT variance-
#   homogeneous, so the pooled residual pool is a SCALE MIXTURE. A scale mixture
#   has the correct total variance but a NARROW CORE and HEAVY TAILS. Scored
#   against it, a moderate observed distance lands too far out (anti-conservative
#   at alpha = 0.05) while an extreme one is still covered (conservative at
#   alpha = 0.01). Under log the mean-variance trend is stabilised, the pool is
#   near-single-scale, and the mismatch disappears.
#
#   Controlled check (simulate_mechanism(), Part 5): with per-gene CV
#   heterogeneity switched OFF the raw exact model is CONSERVATIVE (infl@.05 =
#   0.83 / 0.63 at n = 10 / 40, matching the sqrt-scaling theory in
#   noise_model_current_state.md). Switching heterogeneity on (sd of log CV = 0.7)
#   takes pool excess kurtosis from 0.9 to 29 and infl@.05 from 0.63 to 1.63 at
#   n = 40, while infl@.01 stays at 1.2 -- the real-data pattern, reproduced from
#   a single controlled cause.
#
# WHAT THIS SCRIPT PRODUCES (per cancer x stage x replicate arm x normalisation)
#   Part 1  pool composition   sd_pool vs the gene's own sd (rho), the spread of
#                              per-neighbour sds inside one pool, pool excess
#                              kurtosis -- the mixture, measured directly
#   Part 2  scale vs shape     the (lambda, nu) decomposition above, computed from
#                              this run's own p-values, per cell
#   Part 3  is it the data or the pool?  calibration of z_own (observed scored
#                              against the gene's OWN sd) vs z_pool (against the
#                              pooled null). This is what separates the user's
#                              hypothesis (2) from (1)/(3): if z_own is calibrated
#                              and z_pool is not, the observed distances are fine
#                              and the neighbourhood is at fault
#   Part 4  parameter response tau sweep + kNN sweep + trim sweep, with the share
#                              of genes whose pool actually STOPPED on tau. If tau
#                              almost never binds, "tau too low" is excluded
#                              directly rather than by inference
#   Part 5  simulate_mechanism() the controlled cv-heterogeneity experiment
#
# PREDICTIONS, stated before the run so the result can falsify them:
#   P1  kurt_pool (raw) >> kurt_pool (log), and rises with the per-pool spread of
#       neighbour sds.
#   P2  Splitting genes by kurt_pool, inflation@.05 rises monotonically with it;
#       the low-kurtosis genes are calibrated or conservative.
#   P3  z_own is far better calibrated than z_pool under raw. If BOTH are
#       miscalibrated, the cause is in the data/observed statistic (hypothesis 2),
#       not the neighbourhood, and this whole hypothesis is wrong.
#   P4  frac_stop_tau is small at tau = 0.1. Lowering tau shrinks pools and does
#       NOT fix the inflation; raising it does not either.
#   P5  Tail trimming makes alpha = 0.05 WORSE, not better. It removes the pool's
#       tails, which narrows the core further -- it treats the symptom that is
#       already conservative and worsens the one that is not. (This is the open
#       ToDo "check whether a [5,95] percentile cap improves the raw model": the
#       mechanism above says it will not, and this is the test.)
#
# RUN:  Rscript raw_anticonservative_diagnosis.R          (from TCGA_test/scripts)
# -----------------------------------------------------------------------------

suppressMessages({library(dplyr); library(data.table)})
source("outlier_significance_analysis.R")   # loaders, preprocess_replicates, wrappers
source("config.R")
# tox_null_reimpl.R lives in common/; the repo's scripts are run from more than one
# working directory, so look for it rather than assuming one.
.source_reimpl <- function() {
  cands <- c("tox_null_reimpl.R", "common/tox_null_reimpl.R", "../../common/tox_null_reimpl.R",
             "../common/tox_null_reimpl.R", "../../../common/tox_null_reimpl.R")
  for (p in cands) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Could not locate tox_null_reimpl.R (looked in: ", paste(cands, collapse = ", "), ")")
}
.source_reimpl()
options(width = 200)

OUT_DIR <- file.path(dirname(TOX_TEST_DIR), "raw_diagnosis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==================== CONFIG ====================
# Deliberately SMALL: this is a mechanism experiment, not a calibration sweep.
# null_calibration.R already establishes THAT raw is anti-conservative over all 8
# cancers; the job here is to find out WHY, which needs per-gene detail on a few
# representative cells, not more cells.
DIAG_CANCERS <- CANCER_TYPES[c("kidney cancer", "non small cell lung cancer",
                               "stomach cancer")]
DIAG_STAGES  <- c("Stage I", "Stage III")
N_SPLITS     <- 10L        # per cell; per-gene tables make each split informative
N_PER_GROUP  <- c(10L, 40L)  # the replicate arms where the n-trend is visible
K_CFG        <- list(k_start = 20L, k_step = 1L, k_max = 50L)
TAU_GRID     <- c(0.02, 0.05, 0.10, 0.25, 1.00)   # 1.00 == effectively never stops
K_GRID_DIAG  <- list(list(k_start = 8L,  k_step = 1L, k_max = 30L, name = "k8_30"),
                     list(k_start = 20L, k_step = 1L, k_max = 50L, name = "k20_50"),
                     list(k_start = 40L, k_step = 1L, k_max = 120L, name = "k40_120"))
TRIM_GRID    <- c(0.00, 0.02, 0.05)
MAX_POOL     <- 70000L
N_CORES      <- max(1L, min(16L, parallel::detectCores() - 1L))
set.seed(42)

# ==================== helpers ====================

#' Scale/shape decomposition of a p-value vector.
#'
#' Models the null distance as lambda * S, S standardised symmetric with the tail
#' weight of a t_nu, and the observed as N(0,1). Then
#'   FPR(a) = P(|Z| > lambda * s_a(nu)),  s_a(nu) = (1-a) quantile of standardised |t_nu|
#' so the RATIO of the two implied normal quantiles identifies nu (shape) and
#' their level identifies lambda (scale). lambda > 1 = null too NARROW;
#' nu small / excess kurtosis > 0 = null leptokurtic (narrow core, heavy tail).
#' Gaussian reference: ratio 1.314, excess kurtosis 0, lambda 1.
#'
#' This is a DESCRIPTIVE two-number summary of where a p-value distribution
#' departs from uniform, not a fitted generative model. Its value is that scale
#' and shape have opposite fingerprints, so it says which one is broken.
fit_scale_shape <- function(p) {
  f05 <- mean(p < 0.05); f01 <- mean(p < 0.01)
  if (!is.finite(f05) || !is.finite(f01) || f05 <= 0 || f01 <= 0 || f05 >= 1 || f01 >= 1)
    return(c(lambda = NA, nu = NA, excess_kurt = NA, tail_ratio_z = NA))
  z1 <- qnorm(1 - f05 / 2); z2 <- qnorm(1 - f01 / 2); r <- z2 / z1
  sq <- function(nu, a) if (nu > 200) qnorm(1 - a / 2) else qt(1 - a / 2, nu) * sqrt((nu - 2) / nu)
  g  <- function(nu) sq(nu, 0.01) / sq(nu, 0.05) - r
  if (g(2.5) * g(400) > 0)   # outside the t family: lighter than normal, or extreme
    return(c(lambda = z1 / qnorm(0.975), nu = NA, excess_kurt = NA, tail_ratio_z = r))
  nu <- uniroot(g, c(2.5, 400))$root
  c(lambda = z1 / sq(nu, 0.05), nu = nu,
    excess_kurt = if (nu > 4) 6 / (nu - 4) else Inf, tail_ratio_z = r)
}

#' Calibration of a two-sided z-statistic against the standard normal.
#' Used to score z_own and z_pool on the same footing (Part 3).
z_calib <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 50) return(c(fpr05 = NA, fpr01 = NA, med_p = NA))
  p <- 2 * pnorm(-abs(z))
  c(fpr05 = mean(p < 0.05), fpr01 = mean(p < 0.01), med_p = median(p))
}

#' One A/B split of a samples x genes matrix at a target per-group size.
pick_split <- function(n, size) {
  if (size < 3L || n < 2L * size) return(NULL)
  ix <- sample.int(n, 2L * size)
  list(a = ix[seq_len(size)], b = ix[size + seq_len(size)])
}

#' Load one cohort's raw TPM matrix (samples x genes), TOX's own gene filter applied.
load_cohort <- function(pid, stage, kept_ids) {
  d <- tryCatch(load_stage_data(pid, stage, "cancer", use_constant_healthy = FALSE,
                                norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !is.matrix(d$expression_vectors)) return(NULL)
  m <- d$expression_vectors
  if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
  m[, intersect(kept_ids, colnames(m)), drop = FALSE]
}

# ==================== Part 0: fidelity gate ====================
# Every number below comes from the R port, so the port must be shown to agree
# with the Fortran FIRST. A failure here voids the whole run -- do not "fix" it
# by loosening the tolerance; the two implementations compute the same integer
# pair count, so they should agree exactly.
cat("\n=== Part 0: R port vs Fortran ===\n")
gate_ok <- TRUE
for (nm in c(0L, 1L)) {
  pid <- unname(DIAG_CANCERS[1]); km <- compute_gene_keep_mask(pid)
  m <- load_cohort(pid, DIAG_STAGES[1], names(km)[km])
  if (is.null(m)) { cat("  cohort unavailable -- SKIPPING GATE (results are unvalidated)\n"); gate_ok <- FALSE; break }
  sp <- pick_split(nrow(m), 20L)
  v <- validate_against_fortran(m[sp$a, , drop = FALSE], m[sp$b, , drop = FALSE], nm,
                                K_CFG$k_start, K_CFG$k_step, K_CFG$k_max, 0.1, 0.0, MAX_POOL)
  cat(sprintf("  norm=%d  n=%d genes  max|dp|=%.3g  max|d nbhd|=%d  ok=%s\n",
              nm, v$n_compared, v$max_abs_diff, v$max_abs_diff_nbhd, v$ok))
  if (!isTRUE(v$ok)) gate_ok <- FALSE
}
if (!gate_ok) warning("R port did NOT reproduce the Fortran exactly -- every result below is provisional.")

# ==================== Parts 1-3: per-gene diagnostics ====================

diag_one <- function(m, norm, n_per, kcfg = K_CFG, tau = 0.1, trim = 0.0) {
  sp <- pick_split(nrow(m), n_per); if (is.null(sp)) return(NULL)
  tox_diagnose(m[sp$a, , drop = FALSE], m[sp$b, , drop = FALSE], norm,
               kcfg$k_start, kcfg$k_step, kcfg$k_max, tau, trim, MAX_POOL, verbose = FALSE)
}

cat("\n=== Parts 1-3: per-gene pool diagnostics ===\n")
per_gene <- list(); cells <- list()
for (cname in names(DIAG_CANCERS)) {
  pid <- unname(DIAG_CANCERS[cname])
  km <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL); if (is.null(km)) next
  kept <- names(km)[km]
  for (stage in DIAG_STAGES) {
    m <- load_cohort(pid, stage, kept); if (is.null(m)) next
    for (n_per in N_PER_GROUP) {
      if (nrow(m) < 2L * n_per) next
      for (norm in c(0L, 1L)) {
        d <- do.call(rbind, parallel::mclapply(seq_len(N_SPLITS), function(s) {
          x <- diag_one(m, norm, n_per); if (is.null(x)) return(NULL); x$split <- s; x
        }, mc.cores = N_CORES))
        if (is.null(d) || !nrow(d)) next
        d$cancer <- cname; d$stage <- stage; d$n_per_group <- n_per
        d$norm <- ifelse(norm == 0L, "raw", "log")
        per_gene[[length(per_gene) + 1L]] <- d
        cat(sprintf("  %-26s %-10s n=%-3d %-4s  genes=%d\n", cname, stage, n_per,
                    ifelse(norm == 0L, "raw", "log"), nrow(d)))
      }
    }
  }
}
if (!length(per_gene)) stop("No cohorts produced diagnostics -- check BASE_DATA_DIR.")
PG <- as.data.frame(rbindlist(per_gene, fill = TRUE))
fwrite(PG, file.path(OUT_DIR, "per_gene_diagnostics.csv.gz"))

# ---- Part 1: what the pools are made of -------------------------------------
# rho = sd_own / sd_pool. Its spread across genes is variance heterogeneity
# WITHIN a mean-neighbourhood -- but part of that spread is just estimation noise
# in sd_own at n replicates, whose floor is ~1/sqrt(2(n-1)) on the log scale.
# `sd_log_rho_excess` subtracts that floor, so what remains is heterogeneity the
# neighbourhood genuinely failed to remove.
cat("\n--- Part 1: pool composition (median over genes; mean over splits) ---\n")
p1 <- PG %>% group_by(cancer, stage, n_per_group, norm) %>%
  summarise(genes = dplyr::n(),
            rho_med          = round(median(rho_case, na.rm = TRUE), 3),
            sd_log_rho       = round(sd(log(pmax(rho_case, 1e-9)), na.rm = TRUE), 3),
            sd_log_rho_floor = round(1 / sqrt(2 * (first(n_case) - 1)), 3),
            kurt_pool        = round(median(kurt_pool_case, na.rm = TRUE), 2),
            nb_sd_spread     = round(median(sd_log_nb_sd_case, na.rm = TRUE), 3),
            n_genes_pool     = round(median(n_genes_pool_case), 1),
            frac_below       = round(median(frac_below_case), 3),
            frac_stop_tau    = round(mean(stop_case == "tau"), 3),
            .groups = "drop") %>%
  mutate(sd_log_rho_excess = round(pmax(0, sd_log_rho - sd_log_rho_floor), 3)) %>%
  as.data.frame()
print(p1, row.names = FALSE)
write.csv(p1, file.path(OUT_DIR, "part1_pool_composition.csv"), row.names = FALSE)

# P2: does inflation track pool kurtosis GENE BY GENE? This is the strongest
# single test, because it is a within-cell contrast -- the genes being compared
# share cohort, split, replicate count and every parameter, and differ only in
# how heterogeneous their own neighbourhood happens to be.
# CAVEAT on the row count: neighbouring genes share most of their pool, so rows
# are far from independent and the effective n per quintile is much smaller than
# `genes` suggests. Read a MONOTONE trend across quintiles, not a single contrast,
# and do not attach a p-value to it. A flat profile does not by itself refute the
# mechanism either: within the raw arm every pool is leptokurtic, so the spread of
# kurtosis BETWEEN raw genes is small next to the raw-vs-log gap. The raw-vs-log
# contrast in Part 1 is the primary evidence; this is corroboration.
cat("\n--- Part 1b: FPR by pool-kurtosis quintile (within cell; P2) ---\n")
p1b <- PG %>% filter(is.finite(kurt_pool_case)) %>%
  group_by(cancer, stage, n_per_group, norm) %>%
  mutate(q = ntile(kurt_pool_case, 5)) %>% ungroup() %>%
  group_by(norm, n_per_group, q) %>%
  summarise(kurt_med = round(median(kurt_pool_case), 2),
            rho_spread = round(sd(log(pmax(rho_case, 1e-9))), 3),
            infl_05 = round(mean(p < 0.05) / 0.05, 3),
            infl_01 = round(mean(p < 0.01) / 0.01, 3),
            median_p = round(median(p), 3), genes = dplyr::n(), .groups = "drop") %>%
  as.data.frame()
print(p1b, row.names = FALSE)
write.csv(p1b, file.path(OUT_DIR, "part1b_fpr_by_kurtosis.csv"), row.names = FALSE)

# ---- Part 2: scale vs shape, per cell ---------------------------------------
cat("\n--- Part 2: scale/shape decomposition (lambda>1 = null too NARROW) ---\n")
p2 <- PG %>% group_by(cancer, stage, n_per_group, norm) %>%
  summarise(FPR_05 = round(mean(p < 0.05), 4), FPR_01 = round(mean(p < 0.01), 4),
            median_p = round(median(p), 3),
            lambda = round(fit_scale_shape(p)["lambda"], 3),
            nu_eff = round(fit_scale_shape(p)["nu"], 1),
            excess_kurt_implied = round(fit_scale_shape(p)["excess_kurt"], 2),
            tail_ratio_z = round(fit_scale_shape(p)["tail_ratio_z"], 3),
            .groups = "drop") %>% as.data.frame()
print(p2, row.names = FALSE)
write.csv(p2, file.path(OUT_DIR, "part2_scale_shape.csv"), row.names = FALSE)

# ---- Part 3: the pool, or the data? (P3) ------------------------------------
# z_own scores each observed distance against the SAME gene's own replicate sd;
# z_pool scores it against the neighbourhood null the model actually uses. Both
# are two-sided normal scores, so they are directly comparable.
#   z_own calibrated, z_pool not  -> the neighbourhood is at fault  (hyp. 1 / 3)
#   both miscalibrated            -> the observed distances are genuinely large
#                                    relative to within-group noise (hyp. 2)
# CAVEAT: at small n, sd_own is itself noisy and z_own is t-like, not normal, so
# z_own will look mildly anti-conservative even when correct. Read the raw-vs-log
# CONTRAST in z_own, not its absolute level -- the bias applies equally to both.
cat("\n--- Part 3: observed vs own sd (z_own) and vs pooled null (z_pool); P3 ---\n")
p3 <- PG %>% group_by(cancer, stage, n_per_group, norm) %>%
  summarise(own_FPR05  = round(z_calib(z_own)["fpr05"], 4),
            own_FPR01  = round(z_calib(z_own)["fpr01"], 4),
            pool_FPR05 = round(z_calib(z_pool)["fpr05"], 4),
            pool_FPR01 = round(z_calib(z_pool)["fpr01"], 4),
            exact_FPR05 = round(mean(p < 0.05), 4),
            sd_null_over_own = round(median(sd_null / sd_own_null, na.rm = TRUE), 3),
            .groups = "drop") %>% as.data.frame()
print(p3, row.names = FALSE)
write.csv(p3, file.path(OUT_DIR, "part3_own_vs_pool.csv"), row.names = FALSE)

# ==================== Part 4: parameter response ====================
# tau, k and trim, on the SAME splits, so differences are attributable to the
# parameter and not to the partition. `frac_stop_tau` says how often the tau rule
# is the thing that ends pool growth: if it is near zero at tau = 0.1, then tau is
# not binding and cannot be the cause, whatever the calibration does.
cat("\n=== Part 4: parameter response (same splits throughout) ===\n")
sweep_rows <- list()
for (cname in names(DIAG_CANCERS)[1]) {          # one cancer is enough for a response curve
  pid <- unname(DIAG_CANCERS[cname])
  km <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL); if (is.null(km)) next
  for (stage in DIAG_STAGES[1]) {
    m <- load_cohort(pid, stage, names(km)[km]); if (is.null(m)) next
    for (n_per in N_PER_GROUP) {
      if (nrow(m) < 2L * n_per) next
      grid <- c(
        lapply(TAU_GRID,   function(t) list(kind = "tau",  tau = t,   kcfg = K_CFG, trim = 0.0, lab = sprintf("tau=%.2f", t))),
        lapply(K_GRID_DIAG,function(k) list(kind = "k",    tau = 0.1, kcfg = k,     trim = 0.0, lab = k$name)),
        lapply(TRIM_GRID,  function(f) list(kind = "trim", tau = 0.1, kcfg = K_CFG, trim = f,   lab = sprintf("trim=%.2f", f))))
      for (norm in c(0L, 1L)) for (cfg in grid) {
        # trim is raw-only in the Fortran; mirror that gating here
        if (norm == 1L && cfg$trim > 0) next
        d <- do.call(rbind, parallel::mclapply(seq_len(max(3L, N_SPLITS %/% 2L)), function(s) {
          diag_one(m, norm, n_per, cfg$kcfg, cfg$tau, cfg$trim)
        }, mc.cores = N_CORES))
        if (is.null(d) || !nrow(d)) next
        sweep_rows[[length(sweep_rows) + 1L]] <- data.frame(
          cancer = cname, stage = stage, n_per_group = n_per,
          norm = ifelse(norm == 0L, "raw", "log"), kind = cfg$kind, setting = cfg$lab,
          infl_05 = round(mean(d$p < 0.05) / 0.05, 3),
          infl_01 = round(mean(d$p < 0.01) / 0.01, 3),
          median_p = round(median(d$p), 3),
          n_genes_pool = round(median(d$n_genes_pool_case), 1),
          n_resid_pool = round(median(d$n_resid_pool_case), 1),
          kurt_pool = round(median(d$kurt_pool_case, na.rm = TRUE), 2),
          frac_stop_tau = round(mean(d$stop_case == "tau"), 3),
          stringsAsFactors = FALSE)
      }
    }
  }
}
if (length(sweep_rows)) {
  SW <- as.data.frame(rbindlist(sweep_rows))
  print(SW[order(SW$norm, SW$n_per_group, SW$kind, SW$setting), ], row.names = FALSE)
  write.csv(SW, file.path(OUT_DIR, "part4_parameter_response.csv"), row.names = FALSE)
}

# ==================== Part 5: controlled mechanism experiment ====================

#' Reproduce the failure from ONE controlled cause.
#'
#' Simulates an H0 cohort in which gene means are right-skewed and the per-gene CV
#' is drawn with dispersion `cv_het` -- so `cv_het` is exactly "how much variance
#' heterogeneity survives at a fixed mean", the quantity a mean-neighbourhood
#' cannot condition away. Everything else is held fixed. With `cv_het = 0` the raw
#' exact model should be CONSERVATIVE (the sqrt-scaling result); as `cv_het` rises
#' the pool's kurtosis rises and the alpha = 0.05 inflation follows it, while
#' alpha = 0.01 lags behind -- the real-data signature.
#'
#' Also runs a mean-BOOTSTRAP null on the same pools: that null is a difference of
#' two means-of-n, so it carries the CLT shape the observed statistic has. If the
#' diagnosis is right, it should recover the alpha = 0.05 level; it is not expected
#' to fix alpha = 0.01, because a bootstrap cannot invent tail mass the pool does
#' not contain.
simulate_mechanism <- function(G = 1500L, N = 140L, cv_het = c(0, 0.35, 0.7),
                               n_grid = c(10L, 40L), cv0 = 0.45, B = 4000L, seed = 11L) {
  set.seed(seed)
  make_cohort <- function(G, N, cv_het) {
    mu <- exp(rnorm(G, 2.2, 1.7)); cv <- cv0 * exp(rnorm(G, 0, cv_het))
    m <- matrix(rgamma(N * G, shape = rep(1 / cv^2, each = N),
                       scale = rep(mu * cv^2, each = N)), nrow = N)
    colnames(m) <- paste0("g", seq_len(G)); m
  }
  p_boot <- function(pc, pt, obs, n) {
    a <- colMeans(matrix(sample(pc, n * B, TRUE), nrow = n))
    b <- colMeans(matrix(sample(pt, n * B, TRUE), nrow = n))
    (sum(abs(a - b) >= abs(obs)) + 1) / (B + 1)
  }
  one <- function(m, norm, n, boot) {
    ix <- sample.int(nrow(m), 2 * n); A <- m[ix[1:n], , drop = FALSE]; Bm <- m[ix[n + (1:n)], , drop = FALSE]
    sc <- tox_prepare_sorted(A, norm); st <- tox_prepare_sorted(Bm, norm)
    mc <- colMeans(A); mt <- colMeans(Bm)
    obs <- if (norm == 0L) mc - mt else colMeans(log2(A + 1)) - colMeans(log2(Bm + 1))
    P <- rep(NA_real_, ncol(m)); KU <- P
    for (g in seq_len(ncol(m))) {
      gc_ <- tox_gather(mc[g], sc, 20L, 1L, 50L, 0.1); gt_ <- tox_gather(mt[g], st, 20L, 1L, 50L, 0.1)
      if (length(gc_$pool) < 10 || length(gt_$pool) < 10) next
      P[g] <- if (boot) p_boot(gc_$pool, gt_$pool, obs[g], n)
              else tox_pvalue_exact(gc_$pool / sqrt(n), gt_$pool / sqrt(n), obs[g])
      v <- gc_$pool; KU[g] <- mean(((v - mean(v)) / sd(v))^4) - 3
    }
    k <- !is.na(P)
    data.frame(infl_05 = round(mean(P[k] < 0.05) / 0.05, 3),
               infl_01 = round(mean(P[k] < 0.01) / 0.01, 3),
               median_p = round(median(P[k]), 3),
               kurt_pool = round(median(KU[k], na.rm = TRUE), 2))
  }
  out <- list()
  for (h in cv_het) { m <- make_cohort(G, N, h)
    for (n in n_grid) {
      out[[length(out) + 1L]] <- cbind(cv_het = h, n = n, norm = "raw", null = "exact",     one(m, 0L, n, FALSE))
      out[[length(out) + 1L]] <- cbind(cv_het = h, n = n, norm = "log", null = "exact",     one(m, 1L, n, FALSE))
      out[[length(out) + 1L]] <- cbind(cv_het = h, n = n, norm = "raw", null = "bootstrap", one(m, 0L, n, TRUE))
    } }
  do.call(rbind, out)
}

cat("\n=== Part 5: controlled mechanism experiment ===\n")
cat("cv_het = sd of log(per-gene CV): the variance heterogeneity a mean-neighbourhood cannot remove.\n")
SIM <- simulate_mechanism()
print(SIM, row.names = FALSE)
write.csv(SIM, file.path(OUT_DIR, "part5_mechanism_simulation.csv"), row.names = FALSE)

cat("\nWritten to:", OUT_DIR, "\n")
cat("HOW TO READ THE RESULT\n")
cat("  P1/P2 hold and P3 shows z_own calibrated while z_pool is not\n")
cat("    -> the raw null's SHAPE (a scale-mixture pool) is the cause; tau and k are not.\n")
cat("       Fix direction: give the null the CLT shape (mean-bootstrap), or remove the\n")
cat("       mixture by rescaling each neighbour's residuals to unit sd and re-scaling by\n")
cat("       the target gene's own sd, or keep log as the default for raw-scale data.\n")
cat("  P3 fails (z_own miscalibrated too)\n")
cat("    -> hypothesis (2): observed distances are large relative to within-group noise;\n")
cat("       no neighbourhood change will fix it and the statistic itself must be revisited.\n")
cat("  Part 4 shows a tau or k setting that removes the inflation\n")
cat("    -> hypothesis (1)/(3) after all; adopt that setting and re-run null_calibration.R.\n")
