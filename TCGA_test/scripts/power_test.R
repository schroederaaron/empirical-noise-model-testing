#!/usr/bin/env Rscript
# power_test.R
# =============================================================================
# POWER BENCHMARK for TOX vs edgeR / limma-voom / DESeq2 -- ANALYSIS script.
#
#   Rscript TCGA_test/scripts/power_test.R                 # all parts (R only if POWER_REAL_RDS)
#   POWER_PARTS=Q,C Rscript TCGA_test/scripts/power_test.R # a subset
#   POWER_REAL_RDS=coad_h.rds,coad_s4.rds,... POWER_PARTS=R Rscript TCGA_test/scripts/power_test.R
#   Rscript TCGA_test/scripts/power_report.R               # metrics, tables, plots
#
# Two scripts (docs/power_test_implementation_plan.md section 1):
#   power_test.R    simulate / thin, run every method, SAVE RAW PER-GENE RESULTS. No metrics.
#   power_report.R  reads power_raw/, computes every metric (common/power_metrics.R),
#                   writes power_out/*.csv and every plot. Re-runnable without Fortran.
#
# Companion of calibration_test.R (sourced for draw_counts, common_filter and the
# edgeR / limma / DESeq2 runners; its main block does not run when sourced). Design
# rationale: docs/power_benchmark_design.md. Run from the Tensor-Omics root (rcpp/
# and external/ resolve relative to the working directory); calibration_test.R and
# common/ are found relative to THIS script's location (COMMON_DIR below). Every TOX
# p-value comes from the compiled Fortran.
#
# PARTS
#   R   REAL-DATA SIGNAL INJECTION -- the primary evidence. Each base cohort (an RDS
#       from export_power_cohort.R; plan: COAD / LUAD / KIRC x healthy / Stage IV) is
#       used WHOLE: drop one sample at random if S is odd, split at random into halves
#       A / B (n = floor(S/2) each), inject known log2 FCs by binomial thinning (Gerard
#       2020), balanced so both halves lose the same expected relative abundance
#       (Delta = 0, nulls null in TPM). n_rounds repeats, each a fresh split AND a fresh
#       DE draw. Runs the input x arm grid (as Part N). Checks: gene-length recovery
#       (count/TPM must be constant per gene), one sample per patient (TCGA barcodes),
#       realised vs intended log2 FC per injected gene.
#   P   Main simulated grid: nb / lnpois / tpois x n_rep, balanced effect draw (I3).
#   P+  Dispersion heterogeneity (I2): phi_g = phi_trend(mu_g) exp(sigma_d z - sigma_d^2/2),
#       phi_trend(mu) = phi_inf + c / mu, sigma_d in {0, .5, 1, 2} x sigma_hat. Null arm
#       (pi1 = 0) and power arm (pi1 = 0.1). Parameters from estimate_dispersion_trend.R
#       (POWER_DISP_PARAMS=<rds>); built-in placeholders otherwise (manifest: disp_source).
#   Q   pi1 x effect heterogeneity + oracle (true-null-only) pools + pool diagnostics
#       (tox_diagnose on a gene subset, gated by validate_against_fortran).
#   C   Composition stress, Delta-TARGETED (I4): target Delta in {0, .25, .5, 1}, scored
#       against the biological truth beta. Breakdown point computed by the report.
#   S   STRESS ONLY: `bimodal`. Not for ranking methods.
#   N   Input scale x correction under nuisance (depth, dispersion model).
#   V   Variance-only nulls (case-group dispersion x hv_factor).
#   D   DESeq2 diagnostic (I7): nb, null and power arms, DESeq2 default vs
#       independentFiltering = TRUE vs cooksCutoff = FALSE. No TOX.
#
# RAW OUTPUT (PCFG$raw_dir/<part>/; not committed)
#   results_<job>.csv.gz   one row per (method, gene): p, padj, lfc, tie, floor_th,
#                          nbhd_own_case, nbhd_own_control, lfc_scale
#   truth_<job>.csv.gz     one row per gene: is_de, true_lfc, expr, kept, is_hv
#                          (+ realised_lfc in Part R)
#   diag_<job>.csv.gz      Part Q: pool diagnostics on a gene subset
#   manifest.csv           one row per job: cell keys, round, seed, Delta, effect-draw
#                          diagnostics, Part R cohort fields, git commit
#   pcfg_dump.txt          the full PCFG of the run
# =============================================================================

# ------------------------------------------------------------ dependencies

# ---- locate common/ (config.R, utils.R, ...) from THIS script's own location ----
# Slurm runs the scripts from the Tensor-Omics root, where a plain source("config.R")
# would pick up whatever copy sits in the working directory (e.g. a stale flat one).
# So resolve common/ relative to the script file (<script>/../../common); the other
# entries are fallbacks for interactive use, with a flat copy in "." last.
if (!exists("COMMON_DIR")) COMMON_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand <- c(if (length(f)) file.path(dirname(normalizePath(f[1])), "..", "..", "common"),
            "analysis/common", "experiments/Noise_Model_Test/common", "common", ".")
  hit <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("common/ not found (need config.R); looked in: ", paste(cand, collapse = ", "))
  normalizePath(hit[1])
})

.find_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("None of these files found (run from the calibration_test.R ",
                         "working directory): ", paste(candidates, collapse = ", "))
  hit[1]
}

# Simulator (draw_counts), common_filter, the edgeR / limma / DESeq2 runners,
# load_or_install, LIB_DIR, CFG$run_tox, and the TOX wrappers.
source(.find_file(c(file.path(COMMON_DIR, "..", "Simulated_data", "scripts", "calibration_test.R"),
                     "calibration_test.R")))
source(file.path(COMMON_DIR, "tox_null_reimpl.R"))      # Part Q pool diagnostics
load_or_install("parallel")

#' gzip-compressed CSV via base R (data.table::fwrite compression needs a zlib-enabled build).
write_gz <- function(x, f) { con <- gzfile(f, "w"); on.exit(close(con)); write.csv(x, con, row.names = FALSE) }

# ------------------------------------------------------------ configuration

ALL_PARTS <- c("R", "P", "P+", "Q", "C", "S", "N", "V", "D")
.env_parts <- Sys.getenv("POWER_PARTS", paste(ALL_PARTS, collapse = ","))

# Mean-dependent dispersion for P+ (I2). Estimated from the TCGA healthy cohorts by
# estimate_dispersion_trend.R; sigma_hat is the RAW gene-wise scatter, an UPPER bound
# on true heterogeneity (sampling noise inflates it). The fallback values are
# PLACEHOLDERS, recorded as disp_source = "placeholder" in the manifest.
.disp_params <- local({
  f <- Sys.getenv("POWER_DISP_PARAMS", "")
  if (nzchar(f)) { d <- readRDS(f); d$source <- normalizePath(f); d }
  else list(phi_inf = 0.15, c = 2, sigma_hat = 0.8, source = "placeholder")
})

PCFG <- list(
  parts      = trimws(strsplit(.env_parts, ",")[[1]]),
  # Env overrides for a smoke test: POWER_GENES=1000 POWER_ROUNDS=2 POWER_CORES=4
  n_genes    = as.integer(Sys.getenv("POWER_GENES",  "5000")),
  n_rounds   = as.integer(Sys.getenv("POWER_ROUNDS", "20")),  # >= 20 for anything reported
  n_cores    = as.integer(Sys.getenv("POWER_CORES",  "32")),
  seed       = 20260923L,
  lib_size   = 4e7,
  disp       = 0.2,

  # |log2 FC| ~ log-uniform on this range, sign symmetric, then BALANCED (I3): the
  # side that moves more relative abundance is shrunk by a common factor.
  lfc_range  = c(0.25, 4),
  frac_up    = 0.5,

  # Part grids. A cell = one combination; each cell runs n_rounds datasets.
  R = list(input_rds = Sys.getenv("POWER_REAL_RDS", ""), pi1 = 0.10, min_mean_count = 10,
           len_tol = 0.05,     # max per-gene IQR of log2(count/TPM), after sample scaling
           lfc_tol = 0.5),     # |realised - intended log2 FC| flagged above this
  P = list(dists = c("nb", "lnpois", "tpois"), n_reps = c(3L, 5L, 10L), pi1 = 0.10),
  `P+` = list(dists = c("nb", "lnpois"), n_reps = c(3L, 5L, 10L), pi1 = c(0, 0.10),
              disp_model = "hetero", sigma_mult = c(0, 0.5, 1, 2)),
  Q = list(dists = "lnpois", n_reps = 5L,
           pi1 = c(0.01, 0.05, 0.10, 0.30), het = c(1, 0.5)),
  C = list(dists = c("nb", "lnpois"), n_reps = c(3L, 5L, 10L), pi1 = 0.10,
           delta_target = c(0, 0.25, 0.5, 1.0)),
  S = list(dists = "bimodal", n_reps = c(3L, 5L, 10L), pi1 = 0.10),
  N = list(dists = "lnpois", n_reps = c(3L, 5L, 10L), pi1 = 0.10,
           depth = c("equal", "random", "confounded"),
           disp_model = c("const", "trend", "gene")),
  V = list(dists = "lnpois", n_reps = c(3L, 5L, 10L), pi1 = 0.10, het = c(1, 0.5),
           hv_frac = 0.05),
  D = list(dists = "nb", n_reps = c(3L, 5L, 10L), pi1 = c(0, 0.10)),

  disp_hetero = .disp_params,

  # Nuisance models (Parts N, V). Defaults of the other parts are unchanged:
  # equal depth, constant dispersion, no variance-only nulls.
  depth_sd       = 0.3,           # "random": log library size ~ N(log lib_size, 0.3)
  depth_confound = 1.5,           # "confounded": group B libraries this much deeper
  disp_trend     = c(a = 0.02, b = 5), disp_scatter = 0.3,   # phi = (a + b/mu) * lognormal
  disp_gene_sd   = 0.7,           # "gene": sd of log phi_g around log(disp)
  hv_factor      = 4,             # Part V: case-group dispersion multiplier

  # Parts N and R: input x arm grid. src: "tpm" | "cnt" | "tmm" | "mor".
  inputs = list(
    list(name = "TPM",            src = "tpm", centre = FALSE),
    list(name = "TPM-centred",    src = "tpm", centre = TRUE),
    list(name = "counts",         src = "cnt", centre = FALSE),
    list(name = "counts-centred", src = "cnt", centre = TRUE),
    list(name = "TMM",            src = "tmm", centre = FALSE),
    list(name = "MoR",            src = "mor", centre = FALSE)
  ),
  # TOX-log-boot is left out: it matched TOX-log in every Part P cell.
  grid_arms = c("TOX-log", "TOX-log-blocked", "TOX-raw", "TOX-raw-blocked"),

  run_edger  = TRUE,
  run_limma  = TRUE,
  run_deseq2 = TRUE,           # the slow one

  # TOX: the configuration null_calibration.R fixed on evidence (k20_50, tau 0.1).
  kcfg     = list(k_start = 20L, k_step = 1L, k_max = 50L, name = "k20_50"),
  tau      = 0.1,
  max_pool = 70000L,

  # engine: "exact" = tox_noise_model_exact (sqrt-scaled pairwise, deterministic)
  #         "bootstrap" = tox_noise_model (null = 0 pooled iid; null = 1 gene-blocked,
  #                        enumerated exactly for n_rep <= 5 at k_max = 50, sampled above)
  # input:  "tpm" | "tpm_centred" (median-centred obs) | "tmm" (TMM-CPM, linear)
  tox_arms = list(
    list(label = "TOX-log",         engine = "exact",     norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-raw",         engine = "exact",     norm = 0L, null = 0L, input = "tpm"),
    list(label = "TOX-log-boot",    engine = "bootstrap", norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-log-blocked", engine = "bootstrap", norm = 1L, null = 1L, input = "tpm"),
    list(label = "TOX-raw-blocked", engine = "bootstrap", norm = 0L, null = 1L, input = "tpm")
  ),
  # Part C only: the input-scale comparison (TOX-log, exact, pooled).
  tox_input_arms = list(
    list(label = "TOX-log-TPM",         engine = "exact", norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-log-TPM-centred", engine = "exact", norm = 1L, null = 0L, input = "tpm_centred"),
    list(label = "TOX-log-TMM",         engine = "exact", norm = 1L, null = 0L, input = "tmm")
  ),
  # Part Q only: these tox_arms are ALSO run with a true-null-only pool, as
  # "<label>-oracle". Must use input "tpm" or "tmm" (median centring would be
  # computed over a different gene set in the oracle calls).
  oracle_arms = c("TOX-log", "TOX-log-blocked"),
  # Part Q pool diagnostics (exact engine, pooled, log): genes per job, and genes
  # the R re-implementation is checked on against the Fortran first.
  diag = list(n_genes = 200L, n_check = 200L),

  raw_dir = Sys.getenv("POWER_RAW_DIR", if (exists("TOX_TEST_DIR"))
              file.path(dirname(TOX_TEST_DIR), "power_raw") else "power_raw")
)

RUN_TOX <- isTRUE(CFG$run_tox)
if (!RUN_TOX) message("TOX wrappers not loaded -- running reference methods only.")
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# 1. TRUTH AND SIMULATION
# =============================================================================

#' Draw the DE design shared by every part: which genes, sign, |log2 FC|.
draw_effects <- function(G, n_de, frac_up = PCFG$frac_up, lfc_range = PCFG$lfc_range) {
  beta <- numeric(G)
  de <- if (n_de > 0) sort(sample.int(G, n_de)) else integer(0)
  if (n_de > 0) {
    sgn <- ifelse(runif(n_de) < frac_up, 1, -1)
    beta[de] <- sgn * exp(runif(n_de, log(lfc_range[1]), log(lfc_range[2])))
  }
  list(beta = beta, de = de, is_de = seq_len(G) %in% de)
}

#' Balance a signed log2 FC vector q over abundance weights w (I3 / Part R).
#' The relative abundance the up side moves must equal what the down side moves;
#' the side that moves more has its log2 FCs shrunk by one common factor s in (0, 1]
#' (uniroot; a root exists whenever both sides are non-empty, since that side's mass
#' -> 0 as s -> 0). Signs never change.
#'   composition (thinning = FALSE): up gains sum w (2^q - 1), down loses sum w (1 - 2^q)
#'   thinning    (thinning = TRUE):  group A loses sum_up w (1 - 2^-q), group B sum_dn w (1 - 2^q)
#' ok = FALSE when only one side is non-empty (cannot balance).
balance_effects <- function(q, w, thinning = FALSE) {
  up <- q > 0; dn <- q < 0
  if (!any(up) || !any(dn)) return(list(q = q, ok = !any(up | dn)))
  mu <- function(s) sum(w[up] * (if (thinning) 1 - 2^(-s * q[up]) else 2^(s * q[up]) - 1))
  md <- function(s) sum(w[dn] * (1 - 2^(s * q[dn])))
  a <- mu(1); b <- md(1)
  if (a > b)      q[up] <- q[up] * uniroot(function(s) mu(s) - b, c(0, 1), tol = 1e-12)$root
  else if (b > a) q[dn] <- q[dn] * uniroot(function(s) md(s) - a, c(0, 1), tol = 1e-12)$root
  list(q = q, ok = TRUE)
}

#' Effect-draw diagnostics (manifest; I3). `eff` = realised TPM-space log2 FC of the
#' DE genes, `drawn` = the draw before balancing.
effect_diag <- function(eff, drawn, lfc_range = PCFG$lfc_range) {
  if (!length(eff)) return(list(eff_lfc_median_abs = NA_real_, frac_sign_flip = NA_real_,
                                frac_below_floor = NA_real_))
  list(eff_lfc_median_abs = median(abs(eff)),
       frac_sign_flip     = mean(sign(eff) != sign(drawn)),
       frac_below_floor   = mean(abs(eff) < lfc_range[1]))
}

#' Per-sample composition matrices for both groups.
#'
#' mode = "relative": balanced draw (I3), so sum_DE pi_a 2^beta = sum_DE pi_a and
#'   every null keeps pi_B == pi_a exactly (Delta = 0) with NO shift of the DE effects.
#'   The per-sample renormalisation within the DE set is kept for het < 1 (a random
#'   responder subset is balanced only in expectation); at het = 1 it is a no-op.
#' mode = "absolute": DE genes change absolute output, then the whole sample is
#'   renormalised, so every null moves by -Delta in TPM; true_lfc = beta (biology).
#'   delta_target = NULL keeps the raw draw (Part N); a number (Part C, I4) balances
#'   the draw (Delta = 0) and then scales the UP side's log2 FCs by s >= 1 until
#'   log2 sum(pi_a 2^beta) = delta_target.
#'
#' het < 1: each (DE gene, case sample) responds with probability het; a
#'   non-responder keeps pi_A. This inflates the DE gene's WITHIN-case variance,
#'   which is the only route by which DE genes can contaminate TOX's residual pools.
make_truth_power <- function(G, pi1, n_rep, het = 1, mode = c("relative", "absolute"),
                             frac_up = PCFG$frac_up, delta_target = NULL) {
  mode <- match.arg(mode)
  pi_a <- rlnorm(G, 0, 2); pi_a <- pi_a / sum(pi_a)
  ef <- draw_effects(G, round(pi1 * G), frac_up)
  de <- ef$de; drawn <- ef$beta
  beta <- if (mode == "relative" || !is.null(delta_target)) balance_effects(drawn, pi_a)$q else drawn
  if (mode == "absolute" && !is.null(delta_target) && delta_target > 0) {
    up <- beta > 0; bu <- beta[up]
    base <- sum(pi_a[!up] * 2^beta[!up])
    s <- uniroot(function(s) log2(base + sum(pi_a[up] * 2^(s * bu))) - delta_target,
                 c(1, 2), extendInt = "upX", tol = 1e-12)$root
    beta[up] <- s * bu
  }

  resp <- matrix(runif(G * n_rep) < het, G, n_rep)      # only DE rows are used
  pi_B <- matrix(pi_a, G, n_rep)
  for (j in seq_len(n_rep)) {
    eff <- beta * resp[, j]
    if (mode == "relative") {
      if (length(de)) {
        w <- pi_a[de] * 2^eff[de]
        pi_B[de, j] <- w * sum(pi_a[de]) / sum(w)
      }
    } else {
      w <- pi_a * 2^eff
      pi_B[, j] <- w / sum(w)
    }
  }

  if (mode == "relative") {
    cs <- if (length(de)) sum(pi_a[de]) / sum(pi_a[de] * 2^beta[de]) else 1
    true_lfc <- ifelse(ef$is_de, beta + log2(cs), 0)
    stopifnot(all(abs(sweep(pi_B[!ef$is_de, , drop = FALSE], 1, pi_a[!ef$is_de])) == 0))
    delta <- 0
  } else {
    true_lfc <- beta
    delta <- log2(sum(pi_a * 2^beta))       # TPM shift of every null gene is -delta
  }
  c(list(pi_a = pi_a, pi_B = pi_B, true_lfc = true_lfc, is_de = ef$is_de, delta = delta),
    effect_diag(true_lfc[de], drawn[de]))
}

#' Per-gene dispersion for P+ (I2): phi_trend(mu) times a MEAN-ONE lognormal
#' multiplier, so sigma_d changes the scatter around the trend, not the mean level.
disp_hetero <- function(mu, sigma_d, prm = PCFG$disp_hetero)
  (prm$phi_inf + prm$c / pmax(mu, 1e-3)) * exp(sigma_d * rnorm(length(mu)) - sigma_d^2 / 2)

#' Counts and TPM from per-sample compositions. Mirrors simulate_dataset():
#' fragments are drawn per gene with lambda_gj = N_j * pi_gj * L_g / sum(pi * L).
#'
#' depth:      "equal" (N_j = lib_size), "random" (lognormal per sample),
#'             "confounded" (group B = depth_confound x group A).
#' disp_model: "const" (phi = disp for every gene), "trend" (Part N: phi = a + b/mu,
#'             lognormal scatter), "gene" (log phi_g ~ N(log disp, disp_gene_sd)),
#'             "hetero" (P+: disp_hetero at sigma_d).
#' hv:         logical over genes; those get phi * hv_factor in group B only.
#' The defaults draw NO extra random numbers.
simulate_power <- function(truth, dist, n_rep, lib_size = PCFG$lib_size, disp = PCFG$disp,
                           depth = "equal", disp_model = "const", hv = NULL, sigma_d = 0) {
  G <- length(truth$pi_a)
  L <- sample(round(exp(rnorm(3e4, log(2000), 0.8))), G, replace = TRUE)
  lib <- switch(depth,
    equal      = rep(lib_size, 2 * n_rep),
    random     = lib_size * exp(rnorm(2 * n_rep, 0, PCFG$depth_sd)),
    confounded = lib_size * rep(c(1, PCFG$depth_confound), each = n_rep),
    stop("unknown depth: ", depth))
  ec <- function(pi_mat, N) { p <- pi_mat * L; sweep(sweep(p, 2, colSums(p), "/"), 2, N, "*") }
  lamA <- ec(matrix(truth$pi_a, G, n_rep), lib[seq_len(n_rep)])
  lamB <- ec(truth$pi_B, lib[n_rep + seq_len(n_rep)])

  phi <- switch(disp_model,
    const  = disp,
    trend  = (PCFG$disp_trend[["a"]] + PCFG$disp_trend[["b"]] / pmax(rowMeans(lamA), 1e-3)) *
               exp(rnorm(G, 0, PCFG$disp_scatter)),
    gene   = exp(rnorm(G, log(disp), PCFG$disp_gene_sd)),
    hetero = disp_hetero(rowMeans(lamA), sigma_d),
    stop("unknown disp_model: ", disp_model))
  phiB <- if (is.null(hv)) phi else { p2 <- rep_len(phi, G); p2[hv] <- p2[hv] * PCFG$hv_factor; p2 }

  # draw_counts takes a scalar or a per-gene vector (it recycles down the columns;
  # unit-tested in TCGA_test/tests/test_power_sim.R).
  cA <- draw_counts(lamA, dist, phi)
  cB <- draw_counts(lamB, dist, phiB)
  counts <- cbind(cA, cB)
  colnames(counts) <- c(paste0("A", seq_len(n_rep)), paste0("B", seq_len(n_rep)))
  rownames(counts) <- paste0("gene", seq_len(G))
  rate <- counts / L
  tpm  <- sweep(rate, 2, colSums(rate), "/") * 1e6

  list(counts = counts, tpm = tpm, lengths = L,
       group = rep(c("A", "B"), each = n_rep),
       true_lfc = truth$true_lfc, is_de = truth$is_de,
       is_hv = if (is.null(hv)) rep(FALSE, G) else hv,
       expr = log10(truth$pi_a * 1e6),          # strata by TRUE control abundance
       delta = truth$delta)
}

# ------------------------------------------------------ Part R: real cohort

#' Load one base cohort. RDS = list(counts = genes x samples, tpm = samples x genes
#' or genes x samples, [lengths], [project_id, stage]) from export_power_cohort.R.
#'
#' Checks (plan I10), each reported in the manifest:
#'   one sample per patient -- TCGA barcodes: patient = first 12 characters; later
#'     duplicates are dropped (n_dup_removed). Other ID formats: NOT checked
#'     (patient_check = "not_checked"), with a message.
#'   gene-length recovery   -- lengths = per-gene median of count/TPM after removing
#'     each sample's constant. Valid only if the TPM was built from THESE counts, so
#'     the per-gene IQR of log2(count/TPM) must be ~0; genes above len_tol are dropped
#'     (frac_len_bad).
load_real_cohort <- function(path) {
  d <- readRDS(path)
  cnt <- as.matrix(d$counts)
  # Everything here is genes x samples (edgeR convention; run_tox_arm transposes for the
  # Fortran). Raw TCGA `expression_vectors` are samples x genes and must be t()'d first.
  if (nrow(cnt) < ncol(cnt))
    stop(sprintf("%s: counts are %d x %d -- expected genes x samples", path, nrow(cnt), ncol(cnt)))
  if (is.null(rownames(cnt)) || is.null(colnames(cnt)))
    stop(path, ": counts need gene rownames and sample colnames")
  label <- if (!is.null(d$project_id)) paste0(d$project_id, "/", d$stage)
           else sub("\\.rds$", "", basename(path))

  ids <- colnames(cnt)
  if (all(grepl("^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}", ids))) {
    dup <- duplicated(substr(ids, 1, 12))
    if (any(dup)) message(sprintf("%s: %d sample(s) from an already-used patient dropped: %s",
                                  label, sum(dup), paste(head(ids[dup], 5), collapse = ", ")))
    cnt <- cnt[, !dup, drop = FALSE]
    patient_check <- "tcga_barcode"; n_dup <- sum(dup)
  } else {
    message(sprintf("%s: sample IDs are not TCGA barcodes (e.g. %s) -- one-sample-per-patient ",
                    label, ids[1]), "check NOT done; confirm independence before using this cohort")
    patient_check <- "not_checked"; n_dup <- NA_integer_
  }

  frac_len_bad <- NA_real_
  if (!is.null(d$lengths)) {
    L <- as.numeric(d$lengths[rownames(cnt)])
  } else if (!is.null(d$tpm)) {
    tp <- as.matrix(d$tpm)
    if (setequal(rownames(tp), colnames(d$counts))) tp <- t(tp)    # samples x genes -> genes x samples
    miss <- c(setdiff(rownames(cnt), rownames(tp)), setdiff(colnames(cnt), colnames(tp)))
    if (length(miss))
      stop(sprintf("%s: %d count gene/sample IDs absent from tpm, e.g. %s",
                   path, length(miss), paste(head(miss, 3), collapse = ", ")))
    tp <- tp[rownames(cnt), colnames(cnt)]
    r <- cnt / tp; r[!is.finite(r) | r <= 0] <- NA
    r <- sweep(r, 2, apply(r, 2, median, na.rm = TRUE), "/")
    L <- apply(r, 1, median, na.rm = TRUE)
    spread <- apply(log2(r), 1, IQR, na.rm = TRUE)
    len_ok <- is.finite(spread) & spread <= PCFG$R$len_tol
  } else stop("Real cohort RDS needs $lengths or $tpm.")

  ok <- is.finite(L) & L > 0 & rowMeans(cnt) >= PCFG$R$min_mean_count
  if (exists("len_ok", inherits = FALSE)) {
    frac_len_bad <- mean(!len_ok[ok])
    ok <- ok & len_ok
  }
  message(sprintf("%s: %d genes x %d samples; %d kept (length-recovery failures among expressed: %.2f%%)",
                  label, nrow(cnt), ncol(cnt), sum(ok), 100 * frac_len_bad))
  list(counts = cnt[ok, , drop = FALSE], lengths = L[ok], label = label, path = path,
       patient_check = patient_check, n_dup_removed = n_dup, frac_len_bad = frac_len_bad)
}

#' One signal-injected dataset from a WHOLE real cohort (option A).
simulate_real <- function(cohort, pi1) {
  S <- ncol(cohort$counts); n <- S %/% 2L
  cols <- sample.int(S, 2L * n)                    # random halves; one dropped when S is odd
  cnt <- cohort$counts[, cols, drop = FALSE]; L <- cohort$lengths; G <- nrow(cnt)
  grp <- rep(c("A", "B"), each = n)
  tpm_of <- function(x) { r <- x / L; sweep(r, 2, colSums(r), "/") * 1e6 }
  lr <- function(t) log2(rowMeans(t[, grp == "B", drop = FALSE]) / rowMeans(t[, grp == "A", drop = FALSE]))

  tpm0 <- tpm_of(cnt)
  pi_hat <- rowMeans(tpm0) / 1e6                    # cohort composition
  ef <- draw_effects(G, round(pi1 * G))
  bal <- balance_effects(ef$beta, pi_hat, thinning = TRUE); q <- bal$q
  lr0 <- lr(tpm0)

  # Up-regulation of gene g thins group A by 2^-q, down-regulation thins group B by 2^q.
  thin <- function(y, p) matrix(rbinom(length(y), y, p), nrow(y))
  pa <- ifelse(q > 0, 2^-q, 1); pb <- ifelse(q < 0, 2^q, 1)
  cnt[, grp == "A"] <- thin(cnt[, grp == "A", drop = FALSE], pa)
  cnt[, grp == "B"] <- thin(cnt[, grp == "B", drop = FALSE], pb)
  tpm <- tpm_of(cnt)

  # Realised composition shift: median log2 ratio over true nulls (target 0), and
  # realised injected effect = change of the group log-ratio caused by thinning.
  lr1 <- lr(tpm)
  mA <- rowMeans(tpm[, grp == "A", drop = FALSE]); mB <- rowMeans(tpm[, grp == "B", drop = FALSE])
  okn <- !ef$is_de & mA > 1 & mB > 1
  realised <- ifelse(ef$is_de, lr1 - lr0, 0); realised[!is.finite(realised)] <- NA

  colnames(cnt) <- colnames(tpm) <- c(paste0("A", seq_len(n)), paste0("B", seq_len(n)))
  dg <- effect_diag(q[ef$is_de], ef$beta[ef$is_de])
  dev <- abs(realised - q)[ef$is_de]
  c(list(counts = cnt, tpm = tpm, lengths = L, group = grp,
         true_lfc = q, realised_lfc = realised, is_de = ef$is_de,
         is_hv = rep(FALSE, G), expr = log10(pi_hat * 1e6),
         delta = if (bal$ok) 0 else NA_real_, delta_realised = median(lr1[okn] - lr0[okn]),
         S = S, n_per_group = n, n_dropped = S - 2L * n,
         frac_lfc_mismatch = if (length(dev)) mean(!is.finite(dev) | dev > PCFG$R$lfc_tol) else NA_real_),
    dg)
}

# =============================================================================
# 2. METHODS
# Each returns data.frame(gene, p, padj, lfc, tie, floor_th, nbhd_own_case,
# nbhd_own_control, lfc_scale). p = raw p (NA = untested); tie = |log2 FC of group
# means|, the p_eff tie-break.
# =============================================================================

group_lfc <- function(tpm, group)
  abs(log2((rowMeans(tpm[, group == "B", drop = FALSE]) + 1) /
           (rowMeans(tpm[, group == "A", drop = FALSE]) + 1)))

tmm_cpm <- function(counts)
  edgeR::cpm(normLibSizes(DGEList(counts = counts)), normalized.lib.sizes = TRUE, log = FALSE)

#' DESeq2 median-of-ratios: counts (genes x samples) divided by per-sample size factors.
mor_counts <- function(counts)
  sweep(counts, 2L, DESeq2::estimateSizeFactorsForMatrix(counts), "/")

#' Part D (I7): one DESeq() fit, three results() variants. "DESeq2" is identical to
#' calibration_test.R's run_deseq2 (independentFiltering = FALSE, default Cook's).
run_deseq2_variants <- function(counts, group) {
  dds <- DESeq(DESeqDataSetFromMatrix(round(counts), data.frame(group = factor(group)), ~ group),
               quiet = TRUE)
  v <- list(`DESeq2`        = list(independentFiltering = FALSE),
            `DESeq2-IF`     = list(independentFiltering = TRUE),
            `DESeq2-noCook` = list(independentFiltering = FALSE, cooksCutoff = FALSE))
  lapply(v, function(a) {
    r <- as.data.frame(do.call(results, c(list(dds), a)))
    data.frame(gene = rownames(counts), stat = r$pvalue, padj = r$padj,
               lfc = r$log2FoldChange, row.names = NULL)
  })
}

#' TOX via the production Fortran. `mat` = genes x samples, LINEAR scale.
#' `test_only` (row index): compute the p-value for that gene alone. The pool is
#' still built from ALL rows of `mat` -- only the tested set shrinks.
run_tox_arm <- function(mat, group, arm, tie, test_only = NULL) {
  G <- nrow(mat)
  case_rep <- t(mat[, group == "B", drop = FALSE])
  ctrl_rep <- t(mat[, group == "A", drop = FALSE])
  case_means <- colMeans(case_rep); ctrl_means <- colMeans(ctrl_rep)
  obs <- if (arm$norm == 0L) case_means - ctrl_means
         else colMeans(log2(pmax(case_rep, 0) + 1)) - colMeans(log2(pmax(ctrl_rep, 0) + 1))
  valid <- as.integer(is.finite(obs))
  if (!is.null(test_only)) valid[-test_only] <- 0L
  if (isTRUE(arm$centre) || identical(arm$input, "tpm_centred"))
    obs <- obs - median(obs[is.finite(obs)])
  obs[!is.finite(obs)] <- 0

  fn <- if (arm$engine == "exact") tox_compute_noise_pvalues_pipeline_exact
        else tox_compute_noise_pvalues_pipeline
  res <- fn(case_means = as.numeric(case_means), case_replicates = case_rep,
            control_means = as.numeric(ctrl_means), control_replicates = ctrl_rep,
            obs_own = as.numeric(obs), valid_genes_own = valid,
            norm_method = arm$norm,
            k_start = PCFG$kcfg$k_start, k_step = PCFG$kcfg$k_step, k_max = PCFG$kcfg$k_max,
            tau = PCFG$tau, null_method = arm$null, max_pool_size = PCFG$max_pool)
  if (!is.null(res$ierr) && res$ierr != 0L) stop(sprintf("TOX ierr = %d", res$ierr))

  p <- res$pvalues_own; p[p < 0 | p > 1] <- NA
  padj <- rep(NA_real_, G); padj[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
  nc <- res$neighborhood_size_own_case; nt <- res$neighborhood_size_own_control
  # Theoretical per-gene floor where it is known in closed form.
  floor_th <- if (arm$engine == "exact") ifelse(nc > 0 & nt > 0, 1 / (as.numeric(nc) * nt + 1), NA_real_)
              else if (arm$null == 0L) rep(1 / 10001, G) else rep(NA_real_, G)   # N_BOOTSTRAP_DRAWS
  data.frame(gene = rownames(mat), p = p, padj = padj, lfc = obs, tie = tie,
             floor_th = floor_th, nbhd_own_case = nc, nbhd_own_control = nt,
             lfc_scale = if (arm$norm == 0L) "linear" else "log2", row.names = NULL)
}

#' Oracle pool, computed with the production Fortran.
#'
#' The Fortran builds every gene's pool from ALL genes it is handed
#' (`valid_genes_own` only selects which genes get a p-value, not who is pooled),
#' so "pool from true nulls only" is expressed by what we pass in:
#'   * null genes:  one call with the null genes only -> their pools are all-null.
#'   * each DE gene g: one call with (null genes + g), testing g alone -> g's pool
#'     is its null neighbours plus its OWN residuals, exactly as in production,
#'     where a gene is always its own nearest neighbour.
#' Relative to the production arm the only thing removed is OTHER DE genes from
#' the pool -- which is the contamination being measured.
run_tox_oracle <- function(mat, group, arm, tie, is_null) {
  if (isTRUE(arm$centre) || identical(arm$input, "tpm_centred"))
    stop("oracle arm not defined for median-centred input")
  nul <- which(is_null)
  r0 <- run_tox_arm(mat[nul, , drop = FALSE], group, arm, tie[nul])
  lg <- function(x) if (arm$norm == 0L) x else log2(pmax(x, 0) + 1)
  out <- data.frame(gene = rownames(mat), p = NA_real_, padj = NA_real_,
                    lfc = rowMeans(lg(mat[, group == "B", drop = FALSE])) -
                          rowMeans(lg(mat[, group == "A", drop = FALSE])),
                    tie = tie, floor_th = NA_real_, nbhd_own_case = NA_integer_,
                    nbhd_own_control = NA_integer_, lfc_scale = r0$lfc_scale[1])
  cols <- c("p", "floor_th", "nbhd_own_case", "nbhd_own_control")
  out[nul, cols] <- r0[cols]
  last <- length(nul) + 1L
  for (g in which(!is_null))
    out[g, cols] <- run_tox_arm(mat[c(nul, g), , drop = FALSE], group, arm, tie[c(nul, g)],
                                test_only = last)[last, cols]
  out$padj <- NA_real_; ok <- !is.na(out$p); out$padj[ok] <- p.adjust(out$p[ok], "BH")
  out
}

#' Part Q pool diagnostics (plan section 5): tox_diagnose (R port of the exact,
#' pooled, log engine) on a random gene subset, production pool vs oracle pool.
#' Gated: validate_against_fortran must reproduce the Fortran p-values first; the
#' result is stored with gate_ok and the report uses only gate_ok rows.
run_diag_q <- function(tpm, group, is_de) {
  case <- t(tpm[, group == "B", drop = FALSE]); ctrl <- t(tpm[, group == "A", drop = FALSE])
  k <- PCFG$kcfg
  dg <- function(cm, tm, genes) tox_diagnose(cm, tm, 1L, k$k_start, k$k_step, k$k_max, PCFG$tau,
                                             PCFG$max_pool, genes = genes, verbose = FALSE)
  gate <- validate_against_fortran(case, ctrl, 1L, k$k_start, k$k_step, k$k_max, PCFG$tau,
                                   PCFG$max_pool, n_genes_check = PCFG$diag$n_check)
  sel <- sort(sample.int(ncol(case), min(PCFG$diag$n_genes, ncol(case))))
  nul <- which(!is_de)
  prod <- dg(case, ctrl, sel)
  orc <- list(dg(case[, nul, drop = FALSE], ctrl[, nul, drop = FALSE], match(sel[!is_de[sel]], nul)))
  for (g in sel[is_de[sel]]) {
    idx <- c(nul, g)
    orc[[length(orc) + 1L]] <- dg(case[, idx, drop = FALSE], ctrl[, idx, drop = FALSE], length(idx))
  }
  orc <- do.call(rbind, orc)
  out <- rbind(cbind(pool = "production", prod), cbind(pool = "oracle", orc))
  out$is_de <- is_de[match(out$gene_id, colnames(case))]
  out$gate_ok <- gate$ok; out$gate_max_abs_diff <- gate$max_abs_diff
  out
}

run_methods_power <- function(sim, part) {
  keep <- common_filter(sim)
  counts <- sim$counts[keep, , drop = FALSE]
  tpm    <- sim$tpm[keep, , drop = FALSE]
  grp    <- sim$group
  tie    <- group_lfc(tpm, grp)
  out <- list()
  safe <- function(nm, expr) {
    r <- tryCatch(expr, error = function(e) {
      message(sprintf("    %s failed: %s", nm, conditionMessage(e))); NULL })
    if (!is.null(r)) out[[nm]] <<- r
  }
  ref_tag <- function(df) data.frame(gene = df$gene, p = df$stat, padj = df$padj, lfc = df$lfc,
                                     tie = tie, floor_th = NA_real_, nbhd_own_case = NA_integer_,
                                     nbhd_own_control = NA_integer_, lfc_scale = "log2")

  if (PCFG$run_edger)  safe("edgeR",  ref_tag(run_edger (counts, grp)))
  if (PCFG$run_limma)  safe("limma",  ref_tag(run_limma (counts, grp)))
  if (part == "D") {
    v <- tryCatch(run_deseq2_variants(counts, grp), error = function(e) {
      message("    DESeq2 variants failed: ", conditionMessage(e)); list() })
    for (nm in names(v)) out[[nm]] <- ref_tag(v[[nm]])
  } else if (PCFG$run_deseq2) safe("DESeq2", ref_tag(run_deseq2(counts, grp)))

  diag <- NULL
  if (RUN_TOX && part %in% c("N", "R")) {
    # Input x arm grid. Every input sees the SAME genes (common_filter) and samples.
    mats <- list(tpm = tpm, cnt = counts, tmm = NULL, mor = NULL)
    mats$tmm <- tryCatch(tmm_cpm(counts), error = function(e) NULL)
    mats$mor <- tryCatch(mor_counts(counts), error = function(e) NULL)
    for (a in PCFG$tox_arms) if (a$label %in% PCFG$grid_arms)
      for (inp in PCFG$inputs) {
        if (inp$centre && a$norm == 0L) next      # centring is a log-scale shift
        mat <- mats[[inp$src]]; if (is.null(mat)) next
        ai <- modifyList(a, list(input = inp$src, centre = inp$centre))
        safe(paste0(a$label, "|", inp$name), run_tox_arm(mat, grp, ai, tie))
      }
  } else if (RUN_TOX && part != "D") {
    arms <- if (part == "C") PCFG$tox_input_arms else PCFG$tox_arms
    tmm <- if (any(vapply(arms, function(a) a$input == "tmm", logical(1)))) tmm_cpm(counts)
    for (a in arms) {
      mat <- if (a$input == "tmm") tmm else tpm
      safe(a$label, run_tox_arm(mat, grp, a, tie))
    }
  }

  if (part == "Q" && RUN_TOX) {
    is_null <- !sim$is_de[keep]
    for (a in PCFG$tox_arms) if (a$label %in% PCFG$oracle_arms) {
      mat <- if (a$input == "tmm") tmm_cpm(counts) else tpm
      safe(paste0(a$label, "-oracle"), run_tox_oracle(mat, grp, a, tie, is_null))
    }
    diag <- tryCatch(run_diag_q(tpm, grp, sim$is_de[keep]), error = function(e) {
      message("    pool diagnostics failed: ", conditionMessage(e)); NULL })
  }
  list(results = out, keep = keep, diag = diag)
}

# =============================================================================
# 3. JOBS
# =============================================================================

REAL_COHORTS <- list()

#' One row per dataset. Every CELL key (common/power_metrics.R) is filled for every
#' part; keys a part does not vary hold a constant (NA for delta_target / cohort).
build_jobs <- function() {
  jb <- list()
  add <- function(part, g, cohort = NA_character_) {
    k <- 0L
    for (d in g$dists %||% "real") for (n in g$n_reps) for (pi1 in g$pi1)
      for (h in g$het %||% 1) for (fu in g$frac_up %||% PCFG$frac_up)
        for (dp in g$depth %||% "equal") for (dm in g$disp_model %||% "const")
          for (hvf in g$hv_frac %||% 0) for (sm in g$sigma_mult %||% 0)
            for (dt in g$delta_target %||% NA_real_)
              for (i in seq_len(PCFG$n_rounds)) {
                k <- k + 1L
                jb[[length(jb) + 1L]] <<- data.frame(
                  part = part, dist = d, n_rep = n, pi1 = pi1, het = h, frac_up = fu,
                  depth = dp, disp_model = dm, hv_frac = hvf,
                  sigma_mult = sm, sigma_d = sm * PCFG$disp_hetero$sigma_hat,
                  delta_target = dt, cohort = cohort, round = i, job = k)
              }
  }
  for (pt in setdiff(ALL_PARTS, "R")) if (pt %in% PCFG$parts) add(pt, PCFG[[pt]])
  for (ci in seq_along(REAL_COHORTS)) {        # Part R: one cell per base cohort
    co <- REAL_COHORTS[[ci]]
    before <- length(jb)
    add("R", modifyList(PCFG$R, list(n_reps = ncol(co$counts) %/% 2L)), co$label)
    for (j in (before + 1L):length(jb)) jb[[j]]$job <- (ci - 1L) * PCFG$n_rounds + jb[[j]]$round
  }
  if (!length(jb)) return(NULL)
  j <- do.call(rbind, jb)
  # Seed depends on (part, index within part) only, so running a subset of parts
  # reproduces the same datasets as the full run.
  j$seed <- PCFG$seed + 7919L * j$job + 1000003L * match(j$part, ALL_PARTS)
  j$job_id <- paste0(gsub("[^A-Za-z]", "plus", j$part), "_", j$job)
  j
}

run_job <- function(job) {
  set.seed(job$seed)
  sim <- if (job$part == "R") {
    simulate_real(REAL_COHORTS[[job$cohort]], job$pi1)
  } else {
    mode <- if (job$part %in% c("C", "N")) "absolute" else "relative"
    dt <- if (is.na(job$delta_target)) NULL else job$delta_target
    tr <- make_truth_power(PCFG$n_genes, job$pi1, job$n_rep, job$het, mode, job$frac_up, dt)
    hv <- if (job$hv_frac > 0) {
      nul <- which(!tr$is_de)
      seq_along(tr$is_de) %in% nul[sample.int(length(nul), round(job$hv_frac * length(nul)))]
    }
    c(simulate_power(tr, job$dist, job$n_rep, depth = job$depth,
                     disp_model = job$disp_model, hv = hv, sigma_d = job$sigma_d),
      tr[c("eff_lfc_median_abs", "frac_sign_flip", "frac_below_floor")])
  }

  rr <- run_methods_power(sim, job$part)
  dir <- file.path(PCFG$raw_dir, job$part)
  f <- function(kind) file.path(dir, sprintf("%s_%s.csv.gz", kind, job$job_id))
  truth <- data.frame(gene = rownames(sim$counts), is_de = sim$is_de, true_lfc = sim$true_lfc,
                      expr = sim$expr, kept = rr$keep, is_hv = sim$is_hv)
  if (!is.null(sim$realised_lfc)) truth$realised_lfc <- sim$realised_lfc
  write_gz(truth, f("truth"))
  if (length(rr$results))
    write_gz(do.call(rbind, Map(function(m, r) cbind(method = m, r), names(rr$results), rr$results)),
             f("results"))
  if (!is.null(rr$diag)) write_gz(rr$diag, f("diag"))

  num <- function(x) if (is.null(x)) NA_real_ else x
  cbind(job, delta = num(sim$delta), delta_realised = num(sim$delta_realised),
        n_de_total = sum(sim$is_de), n_kept = sum(rr$keep), n_methods = length(rr$results),
        eff_lfc_median_abs = num(sim$eff_lfc_median_abs), frac_sign_flip = num(sim$frac_sign_flip),
        frac_below_floor = num(sim$frac_below_floor),
        S = num(sim$S), n_per_group = num(sim$n_per_group), n_dropped = num(sim$n_dropped),
        split_seed = if (job$part == "R") job$seed else NA_integer_,
        frac_lfc_mismatch = num(sim$frac_lfc_mismatch),
        cohort_file = if (job$part == "R") REAL_COHORTS[[job$cohort]]$path else NA_character_,
        patient_check = if (job$part == "R") REAL_COHORTS[[job$cohort]]$patient_check else NA_character_,
        n_dup_removed = if (job$part == "R") num(REAL_COHORTS[[job$cohort]]$n_dup_removed) else NA_real_,
        frac_len_bad = if (job$part == "R") num(REAL_COHORTS[[job$cohort]]$frac_len_bad) else NA_real_,
        disp_source = PCFG$disp_hetero$source,
        has_diag = !is.null(rr$diag))
}

#' Commit of the scripts (and of Tensor-Omics, whose Fortran produced the TOX p-values).
#' The server copy is rsynced without .git: set POWER_GIT_COMMIT there.
git_commit <- function() {
  env <- Sys.getenv("POWER_GIT_COMMIT", "")
  if (nzchar(env)) return(env)
  g <- function(d) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(d), "rev-parse", "--short", "HEAD"),
                                                     stdout = TRUE, stderr = FALSE))[1],
                            error = function(e) NA_character_)
  paste0("tests:", g(COMMON_DIR) %||% NA, " tox:", g(".") %||% NA)
}

#' Expected p-value floor per engine, and the BH minimum discovery set m* it implies.
#'   exact:            pool = k_max * n_rep residuals,  p_min = 1/(pool^2 + 1)
#'   bootstrap pooled: p_min = 1/(N_BOOTSTRAP_DRAWS + 1) = 1/10001
#'   blocked, enumerated (n_rep <= 5): W = k_max * n_rep^n_rep, p_min = 1/(W^2 + 1)
#' G is the gene count BEFORE filtering, so m* here is an upper bound on the realised one.
resolution_report_power <- function(n_reps, G = PCFG$n_genes, alpha = 0.05) {
  km <- PCFG$kcfg$k_max
  message("\n--- p-value floor and BH minimum discovery set m* = ceil(G p_min / alpha) ---")
  message(sprintf("%6s %-18s %11s %8s", "n_rep", "engine", "p_min", "m*"))
  for (n in sort(n_reps)) {
    fl <- c(exact = 1 / ((km * n)^2 + 1), `bootstrap pooled` = 1 / 10001,
            `blocked enum` = if (n <= 5) 1 / ((km * n^n)^2 + 1) else 1 / 10001)
    for (e in names(fl))
      message(sprintf("%6d %-18s %11.2e %8d", n, e, fl[[e]], as.integer(ceiling(G * fl[[e]] / alpha))))
  }
  message("  (blocked above n_rep = 5 falls back to sampling at the bootstrap floor)\n")
}

# =============================================================================
if (sys.nframe() == 0L) {
  if ("R" %in% PCFG$parts) {
    rds <- trimws(strsplit(PCFG$R$input_rds, ",")[[1]])
    if (!length(rds)) message("Part R skipped: set POWER_REAL_RDS (comma-separated cohort RDS files).")
    for (p in rds) { co <- load_real_cohort(p); REAL_COHORTS[[co$label]] <- co }
  }
  jobs <- build_jobs()
  if (is.null(jobs) || !nrow(jobs)) stop("No jobs -- check POWER_PARTS / POWER_REAL_RDS.")
  message(sprintf("Parts: %s | datasets: %d | cores: %d | raw -> %s",
                  paste(unique(jobs$part), collapse = ","), nrow(jobs), PCFG$n_cores, PCFG$raw_dir))
  if (PCFG$disp_hetero$source == "placeholder" && "P+" %in% jobs$part)
    message("NOTE: P+ uses PLACEHOLDER dispersion parameters; set POWER_DISP_PARAMS ",
            "(estimate_dispersion_trend.R) for a reportable run.")
  if (RUN_TOX) resolution_report_power(unique(jobs$n_rep))

  # A part's directory holds exactly one run: clear it, so stale jobs never mix in.
  for (pt in unique(jobs$part)) {
    d <- file.path(PCFG$raw_dir, pt); unlink(d, recursive = TRUE)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  commit <- git_commit()

  res <- parallel::mclapply(split(jobs, seq_len(nrow(jobs))), function(j)
           tryCatch(run_job(j), error = function(e) {
             message("job ", j$job_id, " failed: ", conditionMessage(e)); NULL }),
         mc.cores = PCFG$n_cores, mc.preschedule = FALSE)
  res <- Filter(is.data.frame, res)
  if (!length(res)) stop("Every job failed.")
  man <- do.call(rbind, res); man$git_commit <- commit
  for (pt in unique(man$part)) {
    d <- file.path(PCFG$raw_dir, pt)
    write.csv(man[man$part == pt, ], file.path(d, "manifest.csv"), row.names = FALSE)
    writeLines(c(paste("git_commit:", commit), paste("date:", format(Sys.time())),
                 capture.output(dput(PCFG))), file.path(d, "pcfg_dump.txt"))
  }
  message(sprintf("\n%d / %d jobs written to %s/. Next: Rscript power_report.R",
                  nrow(man), nrow(jobs), PCFG$raw_dir))
}
