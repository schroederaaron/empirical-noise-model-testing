#!/usr/bin/env Rscript
# power_test.R
# =============================================================================
# POWER BENCHMARK for TOX vs edgeR / limma-voom / DESeq2
#
#   Rscript power_test.R                 # all simulated parts
#   POWER_PARTS=P,Q Rscript power_test.R # a subset
#   POWER_REAL_RDS=/path/cohort.rds POWER_PARTS=R Rscript power_test.R
#
# Companion of calibration_test.R (sourced for the simulator and the reference
# methods; its main block is guarded by sys.nframe() and does not run). Design
# rationale: claude/power_benchmark_design.md. Run from the Tensor-Omics root,
# next to calibration_test.R and config.R (calibration_test.R resolves
# source("config.R") and source("rcpp/tensoromics_functions.R") relative to the
# working directory). Every TOX p-value comes from the compiled Fortran; there is
# no R re-implementation in the loop.
#
# WHAT THIS TESTS
#   Part P -- MAIN POWER GRID. Relative-abundance truth (make_truth's construction:
#     DE renormalised WITHIN the DE set, so nulls are EXACTLY null in TPM and the
#     composition shift Delta = 0). Observation model x replicates. Answers: how
#     much does each method find at a controlled error rate, and where (effect
#     size, expression level)?
#
#   Part Q -- DE FRACTION x EFFECT HETEROGENEITY (+ oracle pool). TOX builds its
#     null from WITHIN-GROUP residuals of mean-neighbour genes. Within-group
#     centring removes a CONSTANT DE shift, so with homogeneous effects a DE gene
#     contributes ordinary noise to its neighbours' pools and pi1 should NOT
#     matter. With HETEROGENEOUS effects (only a fraction `het` of case samples
#     respond, as in tumours) the DE gene's within-case variance is inflated, it
#     enters its neighbours' case pools, and those nulls widen -> power loss that
#     grows with pi1. A mock null (pi1 = 0) can never show this. The oracle arm
#     re-runs the SAME production Fortran with the pool restricted to true-null
#     genes (see run_tox_oracle), so oracle-minus-production isolates
#     contamination and nothing else.
#
#   Part C -- COMPOSITION STRESS. Absolute truth: DE genes change absolute output,
#     nulls do not, so in TPM every null shifts by -Delta. Scored against the
#     BIOLOGICAL truth beta. Sweeps the up-fraction so Delta grows. TOX-log is run
#     on raw TPM, on median-centred TPM (null_calibration.R's beta_centre = 1) and
#     on TMM-CPM. Expected: uncorrected TPM loses FDR control as |Delta| grows;
#     the question is by how much, and whether centring/TMM repair it.
#
#   Part R -- REAL-DATA SIGNAL INJECTION (off unless POWER_REAL_RDS is set).
#     Binomial thinning (Gerard 2020) on a real homogeneous cohort split into two
#     fake groups: real correlation, real dispersion, real outliers, known truth.
#     Implemented directly with rbinom (one line) rather than via seqgendiff, so
#     the composition balance below is under our control. Up-regulation of gene g
#     thins group A by 2^-q, down-regulation thins group B by 2^q; the side that
#     removes more is then scaled down by a common factor (uniroot) so both groups
#     lose the same expected relative abundance -> Delta = 0 by design, nulls null
#     in TPM. The realised Delta (median null log-ratio) is reported as a check;
#     delta = NA marks a round where balancing was impossible.
#
# METRICS (per method, per simulated dataset; mean +/- SE over rounds)
#   Threshold behaviour
#     fdr_nom_a, tpr_nom_a     observed FDR / TPR of BH calls at nominal a
#     tpr_ach05                TPR at ACHIEVED FDR <= 0.05 -- method-fair power,
#                              independent of whether the method's own alpha holds
#     fpr_null05               fraction of true nulls called at nominal 0.05
#     tpr_all_*                TPR with filtered-out true DE genes counted as missed
#   Ranking (threshold-free); computed twice:
#     *_p   ranked by p-value alone (ties stay tied -> trapezoid = mid-rank AUC)
#     *_pe  p-value, then |log2 FC of group means| as the tie-break
#     The gap between the two is the resolution loss of the discrete p grid.
#     auc, pauc05, pauc10 (McClish-standardised: 0.5 random, 1 perfect), ap
#     (average precision; baseline = pi1), prec_top{k}
#   Direction / effect size
#     sign_err05  wrong sign among true-DE calls at nominal 0.05
#     lfc_rmse    RMSE of the reported log2 FC vs truth, tested DE genes
#                 (NA where the method's lfc is not on a log2 scale, i.e. TOX-raw)
#   Detectable effect
#     lfc80       |log2 FC| at which a logistic fit of P(called @ nominal 0.05)
#                 reaches 0.8 (NA if outside the simulated range)
#   Resolution (the BH floor -- claude/power_benchmark_design.md section 3.1)
#     p_floor     smallest realised p;  n_at_floor  genes tied there
#     m_star      ceil(n_tested * p_floor / 0.05): BH rejects a tied floor block
#                 only if it holds >= m_star genes
#     floor_blocks05  1 if n_at_floor < m_star (floor block unrejectable)
#   Stratified (separate CSV): TPR at nominal 0.05 and at achieved FDR 0.05 by
#     |true log2FC| bin and by expression decile.
#
# OUTPUTS (PCFG$out_dir)
#   power_per_round.csv        one row per (part, cell, round, method)
#   power_summary.csv          mean / SE over rounds
#   power_paired_vs_ref.csv    paired per-round difference vs PCFG$paired_ref
#   power_fdr_tpr_curve.csv    TPR on an achieved-FDR grid, mean / SE
#   power_stratified.csv       stratified TPR, mean / SE
#   power_oracle_gap.csv       paired oracle-minus-production difference (Part Q)
#   *.png                      plots
# =============================================================================

# ------------------------------------------------------------ dependencies

.find_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("None of these files found (run from the calibration_test.R ",
                         "working directory): ", paste(candidates, collapse = ", "))
  hit[1]
}

# Simulator (make_truth / draw_counts / simulate_dataset), common_filter, the
# edgeR / limma / DESeq2 runners, load_or_install, LIB_DIR, CFG$run_tox, and the
# TOX wrappers. Its main block does not run when sourced.
source(.find_file(c("calibration_test.R", "Simulated_data/scripts/calibration_test.R")))
load_or_install("parallel")

# ------------------------------------------------------------ configuration

.env_parts <- Sys.getenv("POWER_PARTS", "P,Q,C")

PCFG <- list(
  parts      = trimws(strsplit(.env_parts, ",")[[1]]),
  # Env overrides for a smoke test: POWER_GENES=1000 POWER_ROUNDS=2 POWER_CORES=4
  n_genes    = as.integer(Sys.getenv("POWER_GENES",  "5000")),
  n_rounds   = as.integer(Sys.getenv("POWER_ROUNDS", "20")),  # >= 20 for anything reported
  n_cores    = as.integer(Sys.getenv("POWER_CORES",  "32")),
  seed       = 20260923L,
  lib_size   = 4e7,
  disp       = 0.2,

  alphas     = c(0.01, 0.05, 0.10),
  fdr_grid   = seq(0, 0.30, by = 0.01),
  topk       = c(100L, 250L, 500L),

  # |true log2 FC| ~ log-uniform on this range. Wider than calibration_test.R's
  # |N(1, 0.6)| on purpose: power needs a low end where the problem is HARD, or
  # every method saturates and the ranking metrics cannot separate them.
  lfc_range  = c(0.25, 4),
  lfc_bins   = c(0.25, 0.5, 1, 2, 4),
  frac_up    = 0.5,

  # Part grids. A cell = one combination; each cell runs n_rounds datasets.
  P = list(dists = c("nb", "lnpois", "tpois", "bimodal"), n_reps = c(3L, 5L, 10L),
           pi1 = 0.10, het = 1),
  Q = list(dists = "lnpois", n_reps = 5L,
           pi1 = c(0.01, 0.05, 0.10, 0.30), het = c(1, 0.5)),
  C = list(dists = "lnpois", n_reps = 5L, pi1 = 0.10,
           frac_up = c(0.5, 0.7, 0.9, 1.0)),
  R = list(input_rds = Sys.getenv("POWER_REAL_RDS", ""), n_reps = c(3L, 5L, 10L),
           pi1 = 0.10, min_mean_count = 10),

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

  paired_ref   = "limma",
  plot_methods = c("limma", "edgeR", "DESeq2", "TOX-log", "TOX-raw",
                   "TOX-log-blocked", "TOX-log-oracle", "TOX-log-blocked-oracle"),
  plot_methods_C = c("limma", "edgeR", "DESeq2", "TOX-log-TPM",
                     "TOX-log-TPM-centred", "TOX-log-TMM"),

  out_dir = if (exists("TOX_TEST_DIR"))
              file.path(dirname(TOX_TEST_DIR), "power_out") else "power_out"
)
dir.create(PCFG$out_dir, showWarnings = FALSE, recursive = TRUE)

RUN_TOX <- isTRUE(CFG$run_tox)
if (!RUN_TOX) message("TOX wrappers not loaded -- running reference methods only.")

# Categorical slots, fixed order (never cycled); a method keeps its colour across plots.
PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")

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

#' Per-sample composition matrices for both groups.
#'
#' mode = "relative": make_truth's construction, per case sample. Renormalise only
#'   WITHIN the DE set, so every null gene has pi_B == pi_A exactly (Delta = 0).
#'   true_lfc = the full-response TPM-space log2 FC (= make_truth's true_lfc).
#' mode = "absolute": DE genes change absolute output, then the whole sample is
#'   renormalised, so every null gene moves by -Delta in TPM. true_lfc = beta, the
#'   biological truth, which is what Part C scores against.
#'
#' het < 1: each (DE gene, case sample) responds with probability het; a
#'   non-responder keeps pi_A. This inflates the DE gene's WITHIN-case variance,
#'   which is the only route by which DE genes can contaminate TOX's residual pools.
make_truth_power <- function(G, pi1, n_rep, het = 1, mode = c("relative", "absolute"),
                             frac_up = PCFG$frac_up) {
  mode <- match.arg(mode)
  pi_a <- rlnorm(G, 0, 2); pi_a <- pi_a / sum(pi_a)
  ef <- draw_effects(G, round(pi1 * G), frac_up)
  de <- ef$de; beta <- ef$beta

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
  list(pi_a = pi_a, pi_B = pi_B, true_lfc = true_lfc, is_de = ef$is_de, delta = delta)
}

#' Counts and TPM from per-sample compositions. Mirrors simulate_dataset():
#' fragments are drawn per gene with lambda_g = N * pi_g * L_g / sum(pi * L).
simulate_power <- function(truth, dist, n_rep, lib_size = PCFG$lib_size, disp = PCFG$disp) {
  G <- length(truth$pi_a)
  L <- sample(round(exp(rnorm(3e4, log(2000), 0.8))), G, replace = TRUE)
  ec <- function(pi_mat) { p <- pi_mat * L; lib_size * sweep(p, 2, colSums(p), "/") }

  cA <- draw_counts(ec(matrix(truth$pi_a, G, n_rep)), dist, disp)
  cB <- draw_counts(ec(truth$pi_B), dist, disp)
  counts <- cbind(cA, cB)
  colnames(counts) <- c(paste0("A", seq_len(n_rep)), paste0("B", seq_len(n_rep)))
  rownames(counts) <- paste0("gene", seq_len(G))
  rate <- counts / L
  tpm  <- sweep(rate, 2, colSums(rate), "/") * 1e6

  list(counts = counts, tpm = tpm, lengths = L,
       group = rep(c("A", "B"), each = n_rep),
       true_lfc = truth$true_lfc, is_de = truth$is_de,
       expr = log10(truth$pi_a * 1e6),          # strata by TRUE control abundance
       delta = truth$delta)
}

# ------------------------------------------------------ Part R: real cohort

#' Load the real cohort once. RDS = list(counts = genes x samples,
#'   lengths = named numeric (optional), tpm = genes x samples (optional)).
#' One of lengths / tpm is required. With TPM only, lengths are recovered as the
#' per-gene median of count/TPM after removing each sample's constant.
#' For a TCGA stage: counts from load_counts_matrix(pid, stage), TPM from the
#' raw-TPM loader in null_calibration.R, same samples and genes.
load_real_cohort <- function(path) {
  d <- readRDS(path)
  cnt <- as.matrix(d$counts)
  # Everything here is genes x samples (edgeR convention; run_tox_arm transposes for the
  # Fortran). load_counts_matrix() already returns that; raw TCGA `expression_vectors`
  # are samples x genes and must be t()'d first.
  if (nrow(cnt) < ncol(cnt))
    stop(sprintf("POWER_REAL_RDS: counts are %d x %d -- expected genes x samples", nrow(cnt), ncol(cnt)))
  if (is.null(rownames(cnt)) || is.null(colnames(cnt)))
    stop("POWER_REAL_RDS: counts need gene rownames and sample colnames")
  if (!is.null(d$lengths)) {
    L <- as.numeric(d$lengths[rownames(cnt)])
  } else if (!is.null(d$tpm)) {
    tp <- as.matrix(d$tpm)
    # TPM from load_stage_data() is samples x genes: orient it by its dimnames.
    if (setequal(rownames(tp), colnames(cnt))) tp <- t(tp)
    miss <- c(setdiff(rownames(cnt), rownames(tp)), setdiff(colnames(cnt), colnames(tp)))
    if (length(miss))
      stop(sprintf("POWER_REAL_RDS: %d count gene/sample IDs absent from tpm, e.g. %s",
                   length(miss), paste(head(miss, 3), collapse = ", ")))
    tp <- tp[rownames(cnt), colnames(cnt)]
    r <- cnt / tp; r[!is.finite(r) | r <= 0] <- NA
    r <- sweep(r, 2, apply(r, 2, median, na.rm = TRUE), "/")
    L <- apply(r, 1, median, na.rm = TRUE)
  } else stop("Real cohort RDS needs $lengths or $tpm.")
  ok <- is.finite(L) & L > 0 & rowMeans(cnt) >= PCFG$R$min_mean_count
  message(sprintf("Real cohort: %d genes x %d samples (%d kept)", nrow(cnt), ncol(cnt), sum(ok)))
  list(counts = cnt[ok, , drop = FALSE], lengths = L[ok])
}

#' One signal-injected dataset from the real cohort.
simulate_real <- function(cohort, n_rep, pi1) {
  S <- ncol(cohort$counts)
  if (S < 2 * n_rep) return(NULL)
  cols <- sample.int(S, 2 * n_rep)
  cnt <- cohort$counts[, cols, drop = FALSE]; L <- cohort$lengths; G <- nrow(cnt)
  grp <- rep(c("A", "B"), each = n_rep)

  rate <- cnt / L
  pi_hat <- rowMeans(sweep(rate, 2, colSums(rate), "/"))       # cohort composition

  ef <- draw_effects(G, round(pi1 * G)); q <- ef$beta
  up <- q > 0; dn <- q < 0
  # Expected relative abundance REMOVED from each group; equal => Delta = 0.
  # The side that removes more has its log2 FCs scaled down by a common factor
  # s in (0, 1] until the two match (a root always exists when both sides are
  # non-empty, since that side's loss -> 0 as s -> 0).
  lossA <- function(s) sum(pi_hat[up] * (1 - 2^-(s * q[up])))
  lossB <- function(s) sum(pi_hat[dn] * (1 - 2^(s * q[dn])))
  balanced <- tryCatch({
    if (lossA(1) > lossB(1)) {
      q[up] <- q[up] * uniroot(function(s) lossA(s) - lossB(1), c(1e-6, 1), tol = 1e-10)$root
    } else if (lossB(1) > lossA(1)) {
      q[dn] <- q[dn] * uniroot(function(s) lossB(s) - lossA(1), c(1e-6, 1), tol = 1e-10)$root
    }
    TRUE
  }, error = function(e) FALSE)

  thin <- function(y, p) matrix(rbinom(length(y), y, p), nrow(y))
  pa <- ifelse(up, 2^-q, 1); pb <- ifelse(dn, 2^q, 1)
  cnt[, grp == "A"] <- thin(cnt[, grp == "A", drop = FALSE], pa)
  cnt[, grp == "B"] <- thin(cnt[, grp == "B", drop = FALSE], pb)

  rate <- cnt / L
  tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6
  # Realised composition shift: median log2 ratio over true nulls (target 0).
  mA <- rowMeans(tpm[, grp == "A", drop = FALSE]); mB <- rowMeans(tpm[, grp == "B", drop = FALSE])
  okn <- !ef$is_de & mA > 1 & mB > 1
  delta_real <- median(log2(mB[okn] / mA[okn]))

  colnames(cnt) <- colnames(tpm) <- c(paste0("A", seq_len(n_rep)), paste0("B", seq_len(n_rep)))
  if (is.null(rownames(cnt))) rownames(cnt) <- rownames(tpm) <- paste0("gene", seq_len(G))
  list(counts = cnt, tpm = tpm, lengths = L, group = grp,
       true_lfc = q, is_de = ef$is_de, expr = log10(pi_hat * 1e6),
       delta = if (balanced) 0 else NA_real_, delta_realised = delta_real)
}

# =============================================================================
# 2. METHODS
# Each returns data.frame(gene, stat, padj, lfc, tie) + attr "lfc_scale".
# stat = raw p (NA = untested); tie = |log2 FC of group means|, the p_eff tie-break.
# =============================================================================

tag <- function(df, scale) { attr(df, "lfc_scale") <- scale; df }

group_lfc <- function(tpm, group)
  abs(log2((rowMeans(tpm[, group == "B", drop = FALSE]) + 1) /
           (rowMeans(tpm[, group == "A", drop = FALSE]) + 1)))

tmm_cpm <- function(counts)
  edgeR::cpm(normLibSizes(DGEList(counts = counts)), normalized.lib.sizes = TRUE, log = FALSE)

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
  if (arm$input == "tpm_centred") obs <- obs - median(obs[is.finite(obs)])
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
  # Theoretical per-gene floor where it is known in closed form.
  floor_th <- if (arm$engine == "exact") {
    nc <- res$neighborhood_size_own_case; nt <- res$neighborhood_size_own_control
    ifelse(nc > 0 & nt > 0, 1 / (as.numeric(nc) * nt + 1), NA_real_)
  } else if (arm$null == 0L) rep(1 / 10001, G) else rep(NA_real_, G)   # N_BOOTSTRAP_DRAWS
  tag(data.frame(gene = rownames(mat), stat = p, padj = padj, lfc = obs, tie = tie,
                 floor_th = floor_th, row.names = NULL),
      if (arm$norm == 0L) "linear" else "log2")
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
#' the pool -- which is the contamination being measured. The per-gene calls
#' cost one prepare + gather pass each (n_de calls per dataset); cheap next to
#' DESeq2, and exact for every engine and null_method.
run_tox_oracle <- function(mat, group, arm, tie, is_null) {
  if (arm$input == "tpm_centred") stop("oracle arm not defined for median-centred input")
  G <- nrow(mat); nul <- which(is_null)
  p <- fl <- rep(NA_real_, G)
  r0 <- run_tox_arm(mat[nul, , drop = FALSE], group, arm, tie[nul])
  p[nul] <- r0$stat; fl[nul] <- r0$floor_th
  last <- length(nul) + 1L
  for (g in which(!is_null)) {
    idx <- c(nul, g)
    r <- run_tox_arm(mat[idx, , drop = FALSE], group, arm, tie[idx], test_only = last)
    p[g] <- r$stat[last]; fl[g] <- r$floor_th[last]
  }
  full_obs <- run_tox_obs(mat, group, arm$norm)
  padj <- rep(NA_real_, G); padj[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
  tag(data.frame(gene = rownames(mat), stat = p, padj = padj, lfc = full_obs, tie = tie,
                 floor_th = fl, row.names = NULL),
      if (arm$norm == 0L) "linear" else "log2")
}

run_tox_obs <- function(mat, group, norm) {
  a <- mat[, group == "B", drop = FALSE]; b <- mat[, group == "A", drop = FALSE]
  if (norm == 0L) rowMeans(a) - rowMeans(b)
  else rowMeans(log2(pmax(a, 0) + 1)) - rowMeans(log2(pmax(b, 0) + 1))
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
  ref_tag <- function(df) { df$tie <- tie; df$floor_th <- NA_real_; tag(df, "log2") }

  if (PCFG$run_edger)  safe("edgeR",  ref_tag(run_edger (counts, grp)))
  if (PCFG$run_limma)  safe("limma",  ref_tag(run_limma (counts, grp)))
  if (PCFG$run_deseq2) safe("DESeq2", ref_tag(run_deseq2(counts, grp)))

  if (RUN_TOX) {
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
  }
  list(results = out, keep = keep)
}

# =============================================================================
# 3. METRICS
# =============================================================================

#' Cumulative TP / FP at every distinct score value (higher = more significant).
#' Tied scores form ONE step, so the trapezoid area equals the mid-rank AUC and no
#' threshold ever splits a tie. Untested genes rank last, tied.
roc_steps <- function(score, y) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE); s <- score[o]; yy <- y[o]
  end <- c(which(s[-1] != s[-length(s)]), length(s))
  list(tp = c(0, cumsum(yy)[end]), fp = c(0, cumsum(!yy)[end]),
       thr = c(Inf, s[end]), P = sum(y), N = sum(!y))
}

#' Area under the ROC up to FPR = fmax (linear interpolation at fmax).
roc_area <- function(r, fmax = 1) {
  x <- r$fp / r$N; yv <- r$tp / r$P
  if (fmax < 1) {
    k <- which(x >= fmax)[1]
    if (!is.na(k)) {
      yk <- yv[k - 1] + (yv[k] - yv[k - 1]) * (fmax - x[k - 1]) / (x[k] - x[k - 1])
      x <- c(x[seq_len(k - 1)], fmax); yv <- c(yv[seq_len(k - 1)], yk)
    }
  }
  sum(diff(x) * (head(yv, -1) + tail(yv, -1)) / 2)
}

pauc_std <- function(r, fmax)                       # McClish: 0.5 random, 1 perfect
  0.5 * (1 + (roc_area(r, fmax) - fmax^2 / 2) / (fmax - fmax^2 / 2))

avg_precision <- function(r) {                      # precision at the end of each step
  prec <- r$tp / pmax(r$tp + r$fp, 1)
  sum(diff(r$tp) * prec[-1]) / r$P
}

#' Index of the last step whose ACHIEVED FDR is <= target (0 = none).
ach_step <- function(r, target) {
  n <- r$tp + r$fp; fdr <- ifelse(n > 0, r$fp / n, 0)
  k <- which(fdr <= target & n > 0)
  if (length(k)) k[which.max(r$tp[k])] else 0L
}

#' Score used for the p-then-effect ranking: strict order, ties broken by `tie`.
score_pe <- function(p, tie) {
  o <- order(p, -tie, na.last = TRUE)
  s <- numeric(length(p)); s[o] <- rev(seq_along(o))
  s[is.na(p)] <- NA
  s
}

#' All scalar metrics for one method on one dataset.
metrics_one <- function(res, is_de, true_lfc, n_de_total) {
  p <- res$stat; padj <- res$padj
  tested <- !is.na(p)
  r_p  <- roc_steps(-p, is_de)
  r_pe <- roc_steps(score_pe(p, res$tie), is_de)
  P <- sum(is_de)
  row <- list(n_genes = length(p), n_tested = sum(tested), n_de = P, n_de_total = n_de_total)

  for (a in PCFG$alphas) {
    called <- !is.na(padj) & padj < a
    tp <- sum(called & is_de); fp <- sum(called & !is_de)
    s <- sprintf("%02d", round(100 * a))
    row[[paste0("n_called_", s)]] <- tp + fp
    row[[paste0("fdr_nom_", s)]]  <- if (tp + fp > 0) fp / (tp + fp) else NA_real_
    row[[paste0("tpr_nom_", s)]]  <- tp / max(1, P)
    row[[paste0("tpr_all_nom_", s)]] <- tp / max(1, n_de_total)
    row[[paste0("zero_disc_", s)]] <- as.integer(tp + fp == 0)
  }
  called05 <- !is.na(padj) & padj < 0.05
  row$fpr_null05 <- sum(called05 & !is_de) / max(1, sum(!is_de))

  k05 <- ach_step(r_p, 0.05)
  row$tpr_ach05     <- if (k05 > 0) r_p$tp[k05] / max(1, P) else 0
  row$tpr_all_ach05 <- if (k05 > 0) r_p$tp[k05] / max(1, n_de_total) else 0

  for (nm in c("p", "pe")) {
    r <- if (nm == "p") r_p else r_pe
    row[[paste0("auc_", nm)]]    <- roc_area(r)
    row[[paste0("pauc05_", nm)]] <- pauc_std(r, 0.05)
    row[[paste0("pauc10_", nm)]] <- pauc_std(r, 0.10)
    row[[paste0("ap_", nm)]]     <- avg_precision(r)
  }
  spe <- score_pe(p, res$tie); o <- order(spe, decreasing = TRUE, na.last = TRUE)
  for (k in PCFG$topk) row[[paste0("prec_top", k)]] <- mean(is_de[o[seq_len(min(k, length(o)))]])

  cd <- called05 & is_de
  row$sign_err05 <- if (any(cd)) mean(sign(res$lfc[cd]) != sign(true_lfc[cd])) else NA_real_
  row$lfc_rmse <- if (identical(attr(res, "lfc_scale"), "log2") && any(tested & is_de))
                    sqrt(mean((res$lfc[tested & is_de] - true_lfc[tested & is_de])^2)) else NA_real_

  # Detectable effect: logistic P(called @ nominal 0.05) on log2|true log2 FC|.
  x <- log2(abs(true_lfc[is_de])); yv <- as.integer(called05[is_de])
  row$lfc80 <- NA_real_
  if (length(unique(yv)) == 2L) {
    fit <- suppressWarnings(glm(yv ~ x, family = binomial()))
    b <- coef(fit)
    if (all(is.finite(b)) && b[2] > 0) {
      x80 <- (qlogis(0.8) - b[1]) / b[2]
      if (x80 >= min(x) && x80 <= max(x)) row$lfc80 <- unname(2^x80)
    }
  }

  # Resolution / BH floor.
  pf <- if (any(tested)) min(p[tested]) else NA_real_
  row$p_floor <- pf
  row$n_at_floor <- if (is.finite(pf)) sum(p[tested] == pf) else NA_integer_
  row$de_at_floor <- if (is.finite(pf)) sum(p[tested] == pf & is_de[tested]) else NA_integer_
  row$m_star <- if (is.finite(pf)) ceiling(sum(tested) * pf / 0.05) else NA_real_
  row$floor_blocks05 <- if (is.finite(pf)) as.integer(row$n_at_floor < row$m_star) else NA_integer_
  row$floor_th_median <- if (any(is.finite(res$floor_th))) median(res$floor_th, na.rm = TRUE) else NA_real_
  as.data.frame(row)
}

#' TPR on the achieved-FDR grid (for averaged FDR-TPR curves).
curve_one <- function(res, is_de) {
  r <- roc_steps(-res$stat, is_de)
  data.frame(fdr = PCFG$fdr_grid, tpr = vapply(PCFG$fdr_grid, function(t) {
    k <- ach_step(r, t); if (k > 0) r$tp[k] / r$P else 0 }, numeric(1)))
}

#' TPR by |true log2 FC| bin and by expression decile (true DE genes only).
strata_one <- function(res, is_de, true_lfc, expr) {
  called_nom <- !is.na(res$padj) & res$padj < 0.05
  r <- roc_steps(-res$stat, is_de); k <- ach_step(r, 0.05)
  s <- -res$stat; s[is.na(s)] <- -Inf
  called_ach <- if (k > 0) s >= r$thr[k] else rep(FALSE, length(s))
  dec <- cut(expr, unique(quantile(expr, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE)
  bins <- cut(abs(true_lfc), PCFG$lfc_bins, include.lowest = TRUE)
  mk <- function(type, f) {
    f <- f[is_de]; if (!length(f)) return(NULL)
    agg <- lapply(split(seq_along(f), f), function(ix) c(
      n = length(ix),
      tpr_nom05 = mean(called_nom[is_de][ix]),
      tpr_ach05 = mean(called_ach[is_de][ix])))
    data.frame(stratum_type = type, stratum = names(agg),
               do.call(rbind, agg), row.names = NULL)
  }
  rbind(mk("abs_lfc_bin", as.character(bins)), mk("expr_decile", dec))
}

# =============================================================================
# 4. JOBS
# =============================================================================

build_jobs <- function() {
  jb <- list()
  add <- function(part, g, extra = list())
    for (d in g$dists %||% "real") for (n in g$n_reps) for (pi1 in g$pi1)
      for (h in g$het %||% 1) for (fu in extra$frac_up %||% PCFG$frac_up)
        for (i in seq_len(PCFG$n_rounds))
          jb[[length(jb) + 1L]] <<- data.frame(part = part, dist = d, n_rep = n, pi1 = pi1,
                                               het = h, frac_up = fu, round = i)
  if ("P" %in% PCFG$parts) add("P", PCFG$P)
  if ("Q" %in% PCFG$parts) add("Q", PCFG$Q)
  if ("C" %in% PCFG$parts) add("C", PCFG$C, list(frac_up = PCFG$C$frac_up))
  if ("R" %in% PCFG$parts && nzchar(PCFG$R$input_rds)) add("R", PCFG$R)
  j <- do.call(rbind, jb); j$job_id <- seq_len(nrow(j)); j
}
`%||%` <- function(a, b) if (is.null(a)) b else a

REAL_COHORT <- NULL

run_job <- function(job) {
  set.seed(PCFG$seed + 7919L * job$job_id)
  sim <- if (job$part == "R") {
    simulate_real(REAL_COHORT, job$n_rep, job$pi1)
  } else {
    mode <- if (job$part == "C") "absolute" else "relative"
    tr <- make_truth_power(PCFG$n_genes, job$pi1, job$n_rep, job$het, mode, job$frac_up)
    simulate_power(tr, job$dist, job$n_rep)
  }
  if (is.null(sim)) return(NULL)

  rr <- run_methods_power(sim, job$part)
  is_de <- sim$is_de[rr$keep]; tl <- sim$true_lfc[rr$keep]; ex <- sim$expr[rr$keep]
  n_de_total <- sum(sim$is_de)
  meta <- cbind(job, delta = sim$delta,
                delta_realised = sim$delta_realised %||% NA_real_)

  m <- cu <- st <- list()
  for (nm in names(rr$results)) {
    res <- rr$results[[nm]]
    m[[nm]]  <- cbind(meta, method = nm, metrics_one(res, is_de, tl, n_de_total))
    cu[[nm]] <- cbind(meta, method = nm, curve_one(res, is_de))
    s <- strata_one(res, is_de, tl, ex)
    if (!is.null(s)) st[[nm]] <- cbind(meta, method = nm, s)
  }
  list(metrics = do.call(rbind, m), curve = do.call(rbind, cu),
       strata = do.call(rbind, st))
}

# =============================================================================
# 5. SUMMARIES, REPORT, PLOTS
# =============================================================================

CELL <- c("part", "dist", "n_rep", "pi1", "het", "frac_up")

mean_se <- function(df, keys, cols) {
  g <- interaction(df[keys], drop = TRUE, lex.order = TRUE)
  do.call(rbind, lapply(split(df, g), function(d) {
    out <- d[1, keys, drop = FALSE]
    out$n_rounds <- nrow(d)
    for (c in cols) {
      x <- d[[c]]; ok <- is.finite(x)
      out[[paste0(c, "_mean")]] <- if (any(ok)) mean(x[ok]) else NA_real_
      out[[paste0(c, "_se")]]   <- if (sum(ok) > 1) sd(x[ok]) / sqrt(sum(ok)) else NA_real_
    }
    out
  }))
}

paired_vs_ref <- function(M, ref, cols = c("tpr_ach05", "pauc05_pe", "ap_pe", "tpr_nom_05")) {
  R <- M[M$method == ref, c(CELL, "round", cols)]
  if (!nrow(R)) return(NULL)
  X <- merge(M[M$method != ref, c(CELL, "round", "method", cols)], R,
             by = c(CELL, "round"), suffixes = c("", ".ref"))
  for (c in cols) X[[paste0("d_", c)]] <- X[[c]] - X[[paste0(c, ".ref")]]
  mean_se(X, c(CELL, "method"), paste0("d_", cols))
}

#' Colour by the method's FIXED slot in `ref` (never by rank among those present),
#' so a method keeps its colour whether or not its neighbours ran.
#' Paired per-round oracle-minus-production difference for each oracle arm.
oracle_gap <- function(M, cols = c("tpr_ach05", "tpr_nom_05", "pauc05_pe", "fpr_null05")) {
  orc <- unique(M$method[grepl("-oracle$", M$method)])
  if (!length(orc)) return(NULL)
  do.call(rbind, lapply(orc, function(o) {
    base <- sub("-oracle$", "", o)
    X <- merge(M[M$method == o, c(CELL, "round", cols)],
               M[M$method == base, c(CELL, "round", cols)],
               by = c(CELL, "round"), suffixes = c("", ".prod"))
    if (!nrow(X)) return(NULL)
    for (c in cols) X[[paste0("d_", c)]] <- X[[c]] - X[[paste0(c, ".prod")]]
    X$method <- o
    mean_se(X, c(CELL, "method"), paste0("d_", cols))
  }))
}

method_colours <- function(methods, ref = PCFG$plot_methods) {
  m <- methods[methods %in% ref[seq_len(min(8L, length(ref)))]]
  setNames(PALETTE[match(m, ref)], m)
}

make_plots <- function(M, CU, ST) {
  dir <- PCFG$out_dir
  th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())

  # (1) FDR-TPR curves, Part P, with nominal-alpha points.
  if (any(M$part == "P")) {
    cols <- method_colours(intersect(PCFG$plot_methods, unique(M$method)))
    cs <- mean_se(CU[CU$part == "P" & CU$method %in% names(cols), ],
                  c(CELL, "method", "fdr"), "tpr")
    ms <- mean_se(M[M$part == "P" & M$method %in% names(cols), ], c(CELL, "method"),
                  c(paste0("fdr_nom_", c("01", "05", "10")), paste0("tpr_nom_", c("01", "05", "10"))))
    pts <- do.call(rbind, lapply(c("01", "05", "10"), function(s)
      data.frame(ms[c(CELL, "method")], alpha = paste0("0.", s),
                 fdr = ms[[paste0("fdr_nom_", s, "_mean")]],
                 tpr = ms[[paste0("tpr_nom_", s, "_mean")]])))
    p <- ggplot(cs, aes(fdr, tpr_mean, colour = method)) +
      geom_vline(xintercept = c(0.01, 0.05, 0.10), colour = "grey85", linewidth = 0.3) +
      geom_line(linewidth = 0.6) +
      geom_point(data = pts[is.finite(pts$fdr), ], aes(fdr, tpr, shape = alpha), size = 2) +
      scale_colour_manual(values = cols) +
      facet_grid(dist ~ n_rep, labeller = label_both) +
      coord_cartesian(xlim = range(PCFG$fdr_grid)) +
      labs(title = "TPR vs achieved FDR (lines); BH calls at nominal alpha (points)",
           subtitle = "A point right of its alpha line = nominal FDR not honoured",
           x = "achieved FDR", y = "TPR (mean over rounds)", shape = "nominal alpha") + th
    ggsave(file.path(dir, "fdr_tpr_curves.png"), p, width = 10, height = 9, dpi = 150)

    # (2) power by |log2FC| and by expression decile, at achieved FDR 0.05.
    for (typ in c("abs_lfc_bin", "expr_decile")) {
      sd_ <- ST[ST$part == "P" & ST$stratum_type == typ & ST$method %in% names(cols), ]
      if (!nrow(sd_)) next
      ss <- mean_se(sd_, c(CELL, "method", "stratum"), "tpr_ach05")
      if (typ == "expr_decile") ss$stratum <- as.integer(ss$stratum)
      else ss$stratum <- factor(ss$stratum, levels = levels(cut(1, PCFG$lfc_bins, include.lowest = TRUE)))
      p <- ggplot(ss, aes(stratum, tpr_ach05_mean, colour = method, group = method)) +
        geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
        scale_colour_manual(values = cols) +
        facet_grid(dist ~ n_rep, labeller = label_both) +
        labs(title = sprintf("Power at achieved FDR 0.05 by %s",
                             if (typ == "expr_decile") "control expression decile (1 = lowest)"
                             else "|true log2 FC|"),
             x = NULL, y = "TPR") + th
      ggsave(file.path(dir, sprintf("power_by_%s.png", typ)), p, width = 10, height = 9, dpi = 150)
    }

    # (3) full AUC vs pAUC(0.05): which one separates the methods?
    a <- mean_se(M[M$part == "P" & M$method %in% names(cols), ], c(CELL, "method"),
                 c("auc_p", "pauc05_p"))
    p <- ggplot(a, aes(auc_p_mean, pauc05_p_mean, colour = method)) +
      geom_point(size = 2) + scale_colour_manual(values = cols) +
      facet_grid(dist ~ n_rep, labeller = label_both) +
      labs(title = "Full ROC-AUC vs standardised pAUC (FPR <= 0.05)",
           x = "AUC", y = "pAUC 0.05 (McClish)") + th
    ggsave(file.path(dir, "auc_vs_pauc.png"), p, width = 10, height = 9, dpi = 150)
  }

  # (4) Part Q: power vs pi1, by heterogeneity; oracle vs all.
  if (any(M$part == "Q")) {
    cols <- method_colours(intersect(PCFG$plot_methods, unique(M$method[M$part == "Q"])))
    q <- mean_se(M[M$part == "Q" & M$method %in% names(cols), ], c(CELL, "method"),
                 c("tpr_ach05", "fpr_null05"))
    p <- ggplot(q, aes(pi1, tpr_ach05_mean, colour = method)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
      geom_errorbar(aes(ymin = tpr_ach05_mean - tpr_ach05_se, ymax = tpr_ach05_mean + tpr_ach05_se),
                    width = 0.05, linewidth = 0.4) +
      scale_x_log10() + scale_colour_manual(values = cols) +
      facet_wrap(~ het, labeller = label_both) +
      labs(title = "Power vs DE fraction: het = 1 homogeneous effects, het < 1 responder fraction",
           subtitle = "*-oracle pools true nulls only (same Fortran); its gap to the production arm is pool contamination",
           x = "pi1 (fraction DE, log scale)", y = "TPR at achieved FDR 0.05 (+/- SE)") + th
    ggsave(file.path(dir, "power_vs_pi1.png"), p, width = 10, height = 5.5, dpi = 150)
  }

  # (5) Part C: composition stress, two small multiples on one x.
  if (any(M$part == "C")) {
    cols <- method_colours(intersect(PCFG$plot_methods_C, unique(M$method[M$part == "C"])),
                           PCFG$plot_methods_C)
    cc <- mean_se(M[M$part == "C" & M$method %in% names(cols), ], c(CELL, "method"),
                  c("fdr_nom_05", "tpr_ach05", "delta"))
    lg <- rbind(data.frame(cc[c("frac_up", "method")], delta = cc$delta_mean,
                           metric = "observed FDR @ nominal 0.05", v = cc$fdr_nom_05_mean),
                data.frame(cc[c("frac_up", "method")], delta = cc$delta_mean,
                           metric = "TPR @ achieved FDR 0.05", v = cc$tpr_ach05_mean))
    p <- ggplot(lg, aes(delta, v, colour = method)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
      scale_colour_manual(values = cols) +
      facet_wrap(~ metric, scales = "free_y") +
      labs(title = "Composition stress: every null gene shifts by -Delta in TPM",
           subtitle = "Scored against the biological truth beta. Delta grows with the up-regulated fraction.",
           x = "Delta = log2 sum(pi_a 2^beta)", y = NULL) + th
    ggsave(file.path(dir, "composition_stress.png"), p, width = 10, height = 5, dpi = 150)
  }
}

report <- function(S, PV, OG) {
  line <- function() message(strrep("=", 78))
  cols <- c("tpr_ach05", "tpr_nom_05", "fdr_nom_05", "fpr_null05", "auc_p", "pauc05_p",
            "pauc05_pe", "ap_pe", "lfc80", "floor_blocks05", "zero_disc_05")
  keep <- c(CELL, "method", "n_rounds", paste0(cols, "_mean"))
  for (pt in unique(S$part)) {
    line(); message("PART ", pt, " -- mean over rounds (SE in power_summary.csv)"); line()
    x <- S[S$part == pt, intersect(keep, names(S))]
    num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], round, 4)
    print(x[order(x$dist, x$n_rep, x$pi1, x$het, x$frac_up, x$method), ], row.names = FALSE)
  }
  message("\n  tpr_ach05   TPR at ACHIEVED FDR 0.05: the method-fair power number.")
  message("  tpr_nom_05 / fdr_nom_05  what a user gets at nominal 0.05.")
  message("  auc_p vs pauc05_p  if AUC spreads < pAUC across methods, the AUC was")
  message("              not measuring the difference (design doc P6).")
  message("  pauc05_pe - pauc05_p  resolution loss from ties at the p floor.")
  message("  floor_blocks05  fraction of rounds where the tied floor block was too")
  message("              small for BH to reject at all (m_star rule).")

  # P6: spread of AUC vs pAUC across methods within each cell.
  sp <- do.call(rbind, lapply(split(S, interaction(S[CELL], drop = TRUE)), function(d)
    data.frame(d[1, CELL], auc_range = diff(range(d$auc_p_mean, na.rm = TRUE)),
               pauc05_range = diff(range(d$pauc05_p_mean, na.rm = TRUE)))))
  line(); message("SPREAD ACROSS METHODS per cell (P6: pauc05_range > auc_range?)"); line()
  sp[c("auc_range", "pauc05_range")] <- lapply(sp[c("auc_range", "pauc05_range")], round, 4)
  print(sp, row.names = FALSE)

  if (!is.null(PV)) {
    line(); message("PAIRED DIFFERENCE vs ", PCFG$paired_ref, " (mean +/- SE over rounds)"); line()
    x <- PV; num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], round, 4)
    print(x, row.names = FALSE)
  }
  if (!is.null(OG)) {
    line(); message("ORACLE GAP (Part Q): oracle minus production, paired per round"); line()
    message("  d_tpr_ach05 > 0 = contaminated pools cost power. Expected ~0 at het = 1.")
    x <- OG; num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], round, 4)
    print(x, row.names = FALSE)
  }
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
  jobs <- build_jobs()
  if (is.null(jobs) || !nrow(jobs)) stop("No jobs -- check POWER_PARTS / POWER_REAL_RDS.")
  if ("R" %in% jobs$part) REAL_COHORT <- load_real_cohort(PCFG$R$input_rds)
  message(sprintf("Parts: %s | cells: %d | datasets: %d | cores: %d",
                  paste(unique(jobs$part), collapse = ","),
                  nrow(unique(jobs[CELL])), nrow(jobs), PCFG$n_cores))
  if (RUN_TOX) resolution_report_power(unique(jobs$n_rep))

  res <- parallel::mclapply(split(jobs, jobs$job_id), function(j)
           tryCatch(run_job(j), error = function(e) {
             message("job ", j$job_id, " failed: ", conditionMessage(e)); NULL }),
         mc.cores = PCFG$n_cores, mc.preschedule = FALSE)
  res <- Filter(function(x) is.list(x) && !is.null(x$metrics), res)
  if (!length(res)) stop("Every job failed.")

  M   <- do.call(rbind, lapply(res, `[[`, "metrics"))
  CU  <- do.call(rbind, lapply(res, `[[`, "curve"))
  ST  <- do.call(rbind, lapply(res, `[[`, "strata"))

  num_cols <- setdiff(names(M)[vapply(M, is.numeric, logical(1))], c(CELL, "round", "job_id"))
  S  <- mean_se(M, c(CELL, "method"), num_cols)
  PV <- paired_vs_ref(M, PCFG$paired_ref)

  write.csv(M, file.path(PCFG$out_dir, "power_per_round.csv"), row.names = FALSE)
  write.csv(S, file.path(PCFG$out_dir, "power_summary.csv"), row.names = FALSE)
  if (!is.null(PV)) write.csv(PV, file.path(PCFG$out_dir, "power_paired_vs_ref.csv"), row.names = FALSE)
  write.csv(mean_se(CU, c(CELL, "method", "fdr"), "tpr"),
            file.path(PCFG$out_dir, "power_fdr_tpr_curve.csv"), row.names = FALSE)
  if (!is.null(ST))
    write.csv(mean_se(ST, c(CELL, "method", "stratum_type", "stratum"), c("n", "tpr_nom05", "tpr_ach05")),
              file.path(PCFG$out_dir, "power_stratified.csv"), row.names = FALSE)
  OG <- oracle_gap(M)
  if (!is.null(OG)) write.csv(OG, file.path(PCFG$out_dir, "power_oracle_gap.csv"), row.names = FALSE)

  report(S, PV, OG)
  tryCatch(make_plots(M, CU, ST), error = function(e) message("plotting failed: ", conditionMessage(e)))
  message("\nWrote CSVs and plots to ", PCFG$out_dir, "/")
}
