#!/usr/bin/env Rscript
# calibration_test.R
# =============================================================================
# NULL AND FDR CALIBRATION TEST for TOX vs edgeR / limma-voom / DESeq2
#
#   Rscript calibration_test.R
#
# WHAT THIS TESTS
#   Part A -- MOCK NULL. No differentially expressed genes at all. Every
#     rejection is a false positive and the p-value histogram must be UNIFORM.
#     Needs no distributional assumption and no effect-size threshold, and it
#     exercises the whole chain at once: residual construction, neighbourhood
#     placement, null scaling, p-value resolution, BH denominator.
#
#   Part B -- FDR CALIBRATION. DE genes present. Sweep the nominal FDR and
#     compare observed FDR against it; AUC measures ranking independently of
#     any threshold.
#
# READING PART A -- THREE NUMBERS, IN ORDER
#   1. tail_ratio = frac_p_lt_05 / frac_p_lt_01.  Target ~5 under a uniform
#      null. This is the SPIKE DETECTOR and it matters more than frac_p_lt_05
#      alone: a method can show frac_p_lt_05 = 0.05 while every one of those
#      p-values sits below 0.01, which is a spike at zero, not calibration.
#      BH reads the extreme tail, so the spike is what produces false
#      positives. Ratio near 1 = broken tail even if the aggregate looks fine.
#   2. frac_p_lt_05.  Target 0.05.  >> means the null is too NARROW
#      (anti-conservative); << means too WIDE (conservative -> the
#      high-precision / low-recall pattern).
#   3. n_rejected.  Target ~0. Under a COMPLETE null, BH's FDR control
#      collapses to FWER control, so any rejection at all is a failure.
#
# READING PART B
#   auc  vs  tpr:
#     high auc + low tpr -> ranking is fine, THRESHOLD is miscalibrated
#     low auc            -> the ranking itself is broken
#   n_zero_disc: rounds that produced NO discoveries at all. A non-zero value
#     means the method is UNSTABLE, not merely conservative -- and it is also
#     the signature of the p-value resolution cliff (see resolution_report()).
#
# WHY SIMULATE IN RELATIVE-ABUNDANCE SPACE
#   TPM estimates the molar fraction pi. Defining the truth directly on pi and
#   renormalising only WITHIN the DE set leaves null genes with pi_B == pi_A
#   exactly, so their true TPM log-fold-change is exactly zero. Simulating in
#   count space and converting to TPM does not have this property: the shared
#   denominator shifts every null gene by a global constant and the "nulls"
#   stop being null.
# =============================================================================

source("config.R")

# Fixed library location for both loading AND installing packages. `.libPaths()`
# silently DROPS non-existent directories, so the folder must be created FIRST or
# the prepend is a no-op (installs would fall back to the default user library).
# The path is anchored absolutely so it does not depend on the working directory.
LIB_DIR <- normalizePath("external/docker_r_libs", mustWork = FALSE)
dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIB_DIR, .libPaths()))

# Name the package that is ACTUALLY missing from a load error. `library()` names
# the package being loaded FIRST, e.g.
#   "package or namespace load failed for 'DESeq2': there is no package called 'locfit'"
# so grabbing the first quoted token yields 'DESeq2' and reinstalls it forever. Match
# the specific phrases that name the missing dependency instead. Returns NA when the
# failure is not a missing-package error (e.g. a shared-object load failure).
.missing_pkg <- function(msg) {
  q <- "[‘’'\"`]([[:alnum:]._]+)[‘’'\"`]"
  for (pat in c(paste0("there is no package called ", q),
                paste0("package ", q, " required by"))) {
    m <- regmatches(msg, regexec(pat, msg, perl = TRUE))[[1]]
    if (length(m) >= 2L) return(m[2L])
  }
  NA_character_
}

load_or_install <- function(package_name) {
  install_pkg <- function(pkg) {
    if (requireNamespace("BiocManager", quietly = TRUE))
      BiocManager::install(pkg, lib = LIB_DIR, update = FALSE, ask = FALSE)
    else   # bootstrap path: install BiocManager itself (or any pure-CRAN pkg) from CRAN
      install.packages(pkg, repos = "https://cloud.r-project.org", lib = LIB_DIR, dependencies = TRUE)
  }
  if (!requireNamespace(package_name, quietly = TRUE)) install_pkg(package_name)
  for (attempt in seq_len(8L)) {
    err <- tryCatch({ suppressPackageStartupMessages(library(package_name, character.only = TRUE))
                      return(invisible()) },
                    error = function(e) e)
    miss <- .missing_pkg(conditionMessage(err))
    # A distinct missing dependency -> install it and retry (resolves a chain). If we
    # can't identify one (or it's the package itself, which install won't fix because
    # it's already current), stop and surface the REAL error rather than spinning 8x.
    if (is.na(miss) || identical(miss, package_name))
      stop("Could not load '", package_name, "': ", conditionMessage(err))
    message("  load_or_install: '", package_name, "' -> installing missing dependency '", miss, "' ...")
    install_pkg(miss)
  }
  stop("Could not load '", package_name, "' -- unresolved dependencies after 8 attempts.")
}

load_or_install("BiocManager")
load_or_install("edgeR")
load_or_install("limma")
load_or_install("DESeq2")
load_or_install("ggplot2")
load_or_install("goftest")   # Anderson-Darling / Cramer-von Mises goodness-of-fit tests

# Headless containers often ship with NO fonts, so the cairo PNG device renders
# every glyph as a "tofu" box. If none are present, install DejaVu + fontconfig so
# plot labels render. (Ephemeral containers re-run this each time; bake it into the
# image -- e.g. `RUN pacman -Sy --noconfirm fontconfig ttf-dejavu && fc-cache -f`.)
ensure_fonts <- function() {
  n_fonts <- suppressWarnings(tryCatch(length(system("fc-list", intern = TRUE, ignore.stderr = TRUE)),
                                       error = function(e) 0L))
  if (isTRUE(n_fonts > 0L)) return(invisible())
  message("No system fonts found -- plot text would render as boxes. Installing a font ...")
  cmd <- if (nzchar(Sys.which("pacman")))       "pacman -Sy --noconfirm fontconfig ttf-dejavu"
         else if (nzchar(Sys.which("apt-get"))) "apt-get update && apt-get install -y fontconfig fonts-dejavu-core"
         else if (nzchar(Sys.which("apk")))     "apk add --no-cache fontconfig ttf-dejavu"
         else NA_character_
  if (is.na(cmd)) { warning("No known package manager; install a font + fontconfig in the image."); return(invisible()) }
  try(system(paste(cmd, "&& fc-cache -f"), ignore.stdout = TRUE, ignore.stderr = TRUE), silent = TRUE)
}
ensure_fonts()

options(width = 250, max.print = 1e6)

# ------------------------------------------------------------ configuration

CFG <- list(
  n_genes     = 5000,          # 20000 for the final run; 5000 keeps DESeq2 sane
  n_rounds    = 5,             # independent datasets per cell
  n_reps      = c(3, 5),
  lib_sizes   = c(4e7),        # add 1e6 to probe the low-depth regime
  dists       = c("nb", "lnpois", "tpois", "bimodal"),
  alpha       = 0.05,
  n_de_partB  = 500,
  seed        = 20240101,

  out_dir     = if (exists("TOX_TEST_DIR"))
                  file.path(dirname(TOX_TEST_DIR), "calibration_out") else "calibration_out",

  # `bimodal` on/off probability. At 0.7 with n_rep = 3 roughly 30% of values
  # are exact zeros and whole groups go all-zero, so the arm tests SPARSITY
  # rather than multimodality and every method fails. 0.9 keeps the mixture
  # structure while leaving enough non-zero observations to be a fair test.
  bimodal_on_prob = 0.9,

  run_tox     = TRUE,          # FALSE to sanity-check the harness alone
  tox_verbose = FALSE,         # per-call TOX internals; very noisy
  tox_variants = c("exact", "bootstrap"),
  tox_norm     = c(0L, 1L),

  # NEIGHBOURHOOD SWEEP. k trades off two things that break in opposite
  # directions:
  #   larger k -> finer p-value resolution, more stable null
  #   larger k -> wider mean range in the neighbourhood -> weaker local
  #               exchangeability -> worse AUC
  # There should be a knee. Plot AUC, tail_ratio and n_rejected against k.
  tox_k_grid = list(
    list(name = "k20_50", k_start = 20L, k_step = 5L, k_max = 50L),
    list(name = "k30_70", k_start = 30L, k_step = 5L, k_max = 70L),
    list(name = "k40_70", k_start = 40L, k_step = 5L, k_max = 70L)
  ),
  tox_tau      = 0.1,
  tox_max_pool = 5000L,

  # Residual-pool tail trim (raw norm only). When > 0, each raw (norm 0) arm is
  # ALSO run with this per-tail trim fraction as a separate "_trim" arm, so the
  # trimmed-vs-untrimmed raw comparison appears in the same sweep. 0 disables it.
  tox_trim_frac = 0.05
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

if (CFG$run_tox) {
  ok <- tryCatch({ source("rcpp/tensoromics_functions.R"); TRUE },
                 error = function(e) { message("TOX wrapper not found: ", e$message); FALSE })
  CFG$run_tox <- ok
}

#' Expected p-value floor and the resulting BH cliff, printed before the run.
#'
#' The exact variant enumerates all pairs, so p_min = 1/(pool^2 + 1) with
#' pool = k_max * n_rep. BH rejects k genes only when p_(k) <= alpha*k/G, so
#' genes pinned at the floor need k >= G*p_min/alpha before ANY of them pass.
#' Below that count the method returns zero discoveries no matter how strong
#' the signal -- which reads as conservatism but is pure arithmetic.
#'
#' NOTE: if the baseline (non-exact) variant falls back to Monte Carlo, its
#' floor is 1/(n_MC + 1) and does NOT improve with k_max. Check that constant
#' in the Fortran; it may bind before the pool size does.
resolution_report <- function() {
  message("\n--- p-value resolution (exact variant) ---")
  message(sprintf("%-10s %6s %8s %11s %14s", "k_max", "n_rep", "pool",
                  "p_min", "genes@floor"))
  for (kg in CFG$tox_k_grid) for (n in CFG$n_reps) {
    pool  <- kg$k_max * n
    p_min <- 1 / (pool^2 + 1)
    need  <- CFG$n_genes * p_min / CFG$alpha
    message(sprintf("%-10s %6d %8d %11.2e %14.1f", kg$name, n, pool, p_min, need))
  }
  message("  genes@floor = how many genes must sit at the floor before BH ",
          "returns anything.\n  Keep it well below the expected number of ",
          "true positives.\n")
}

# =============================================================================
# 1. SIMULATION
# =============================================================================

make_truth <- function(n_genes, n_de, lfc_mean = 1.0, lfc_sd = 0.6,
                       frac_up = 0.5, pi_a = NULL) {
  if (is.null(pi_a)) pi_a <- rlnorm(n_genes, 0, 2)
  pi_a <- pi_a / sum(pi_a)

  beta <- numeric(n_genes); de <- integer(0)
  if (n_de > 0) {
    de   <- sample.int(n_genes, n_de)
    sgn  <- ifelse(runif(n_de) < frac_up, 1, -1)
    beta[de] <- sgn * abs(rnorm(n_de, lfc_mean, lfc_sd))
  }

  pi_b <- pi_a
  if (n_de > 0) {
    cs <- sum(pi_a[de]) / sum(pi_a[de] * 2^beta[de])
    pi_b[de] <- pi_a[de] * 2^beta[de] * cs
  }

  stopifnot(abs(sum(pi_b) - 1) < 1e-10)
  true_lfc <- log2(pi_b / pi_a)
  is_de <- seq_len(n_genes) %in% de
  stopifnot(max(abs(true_lfc[!is_de])) == 0)     # nulls are EXACTLY null

  list(pi_a = pi_a, pi_b = pi_b, true_lfc = true_lfc, is_de = is_de)
}

#' Observation model -- the only thing that differs between arms.
#'
#'   nb       negative binomial -- edgeR / DESeq2's exact model (home turf)
#'   lnpois   Poisson-lognormal -- multiplicative biology, then Poisson sampling
#'   tpois    as lnpois but with STANDARDIZED, winsorized t_3 log-noise -- a
#'            heavier-tailed multiplicative noise at the SAME scale as lnpois
#'   bimodal  gene on/off per sample -- violates unimodality outright
#'   poisson  pure shot noise, no biological variance -- the technical floor
#'
#' `tpois` note: a raw t_3 has variance 3, so `exp(t_3 * sd)` is both ~sqrt(3)x
#' too wide AND has NO finite mean (polynomial tails times exp). The next line
#' divides by `mean(mult)`, which a single extreme draw then dominates, so every
#' gene's lambda collapses toward zero and all methods fail (the old bug at
#' n_rep = 5). We STANDARDIZE the t to unit variance (so its scale matches
#' lnpois, controlled by `disp`) and WINSORIZE it to +/- TPOIS_Z_CLAMP so the
#' multiplier keeps a finite, stable mean while staying genuinely heavy-tailed.
TPOIS_DF      <- 3          # Student-t degrees of freedom (heavy tails)
TPOIS_Z_CLAMP <- 6          # winsorize the unit-variance t at +/- this (finite mean)

draw_counts <- function(lambda, dist, disp = 0.2,
                        on_prob = CFG$bimodal_on_prob) {
  G <- nrow(lambda); S <- ncol(lambda)
  if (dist == "nb")
    return(matrix(rnbinom(G*S, mu = as.vector(lambda), size = 1/disp), G, S))

  mult <- switch(dist,
    poisson = matrix(1, G, S),
    lnpois  = exp(matrix(rnorm(G*S, 0, sqrt(log(1+disp))), G, S)),
    tpois   = {
      z <- rt(G*S, df = TPOIS_DF) / sqrt(TPOIS_DF / (TPOIS_DF - 2))  # -> unit variance
      z <- pmax(pmin(z, TPOIS_Z_CLAMP), -TPOIS_Z_CLAMP)             # finite-mean guard
      exp(matrix(z, G, S) * sqrt(log(1+disp)))
    },
    bimodal = matrix(rbinom(G*S, 1, on_prob), G, S) / on_prob,
    stop("unknown dist: ", dist))
  mult <- mult / mean(mult)
  matrix(rpois(G*S, as.vector(lambda * mult)), G, S)
}

simulate_dataset <- function(truth, dist, n_rep, lib_size, disp = 0.2) {
  G <- length(truth$pi_a)
  L <- sample(round(exp(rnorm(3e4, log(2000), 0.8))), G, replace = TRUE)

  # Reads are drawn per FRAGMENT: lambda_g = N * pi_g * L_g / sum(pi*L)
  ec <- function(pi) { p <- pi * L; lib_size * p / sum(p) }

  cA <- draw_counts(matrix(ec(truth$pi_a), G, n_rep), dist, disp)
  cB <- draw_counts(matrix(ec(truth$pi_b), G, n_rep), dist, disp)

  counts <- cbind(cA, cB)
  colnames(counts) <- c(paste0("A", seq_len(n_rep)), paste0("B", seq_len(n_rep)))
  rownames(counts) <- paste0("gene", seq_len(G))

  rate <- counts / L
  tpm  <- sweep(rate, 2, colSums(rate), "/") * 1e6

  list(counts = counts, tpm = tpm, lengths = L,
       group = rep(c("A", "B"), each = n_rep),   # A = control, B = case
       true_lfc = truth$true_lfc, is_de = truth$is_de)
}

# =============================================================================
# 2. METHODS
# Each returns data.frame(gene, stat, padj, lfc); `stat` = raw p-value.
# =============================================================================

#' One filter for every method, expression-only, never touching the truth.
common_filter <- function(sim, min_tpm = 1, min_count = 10)
  rowMeans(sim$tpm) >= min_tpm & rowSums(sim$counts) >= min_count

# edgeR 4.x QL standard: normLibSizes (TMM), and glmQLFit (legacy = FALSE default)
# estimates the dispersion internally so estimateDisp() is not needed; robust = TRUE
# guards against hypervariable genes. Genes are already filtered upstream by the
# shared common_filter() (one gene set for every method), so no filterByExpr here.
run_edger <- function(counts, group) {
  design <- model.matrix(~ group)
  y <- normLibSizes(DGEList(counts = counts, group = group))
  fit <- glmQLFit(y, design, robust = TRUE)
  tt  <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(counts), stat = tt$PValue, padj = tt$FDR,
             lfc = tt$logFC, row.names = NULL)
}

# limma-voom standard: voomLmFit() (modern voom + lmFit, zero-count df aware) + robust
# eBayes. trend = FALSE is the default and correct for voom output.
run_limma <- function(counts, group) {
  design <- model.matrix(~ group)
  y <- normLibSizes(DGEList(counts = counts, group = group))
  fit <- eBayes(voomLmFit(y, design), robust = TRUE)
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(counts), stat = tt$P.Value, padj = tt$adj.P.Val,
             lfc = tt$logFC, row.names = NULL)
}

run_deseq2 <- function(counts, group) {
  dds <- DESeqDataSetFromMatrix(round(counts),
                                data.frame(group = factor(group)), ~ group)
  res <- as.data.frame(results(DESeq(dds, quiet = TRUE),
                               independentFiltering = FALSE))
  data.frame(gene = rownames(counts), stat = res$pvalue, padj = res$padj,
             lfc = res$log2FoldChange, row.names = NULL)
}

#' TOX. group "B" = case (cancer), group "A" = control (healthy).
#'
#' @param kcfg list(name, k_start, k_step, k_max) from CFG$tox_k_grid
run_tox <- function(tpm, group, norm_method = 1L, variant = "exact",
                    kcfg = CFG$tox_k_grid[[1]], trim_frac = 0.05) {
  n_genes <- nrow(tpm)

  # Replicate matrices are SAMPLES x GENES and PRE-LOG (linear TPM) for BOTH
  # norm methods -- the Fortran applies log2(x + 1) internally when
  # norm_method != 0, so we must not log here.
  case_rep <- t(tpm[, group == "B", drop = FALSE])
  ctrl_rep <- t(tpm[, group == "A", drop = FALSE])

  # Means place the kNN neighbourhood; they stay LINEAR in both modes.
  # log2(x+1) is monotone in x, so the ordering is unaffected either way.
  case_means <- colMeans(case_rep, na.rm = TRUE)
  ctrl_means <- colMeans(ctrl_rep, na.rm = TRUE)

  # The observed statistic MUST sit on the residual scale.
  #   norm 0 -> residuals are x - mean(x)                -> linear difference
  #   norm 1 -> residuals are log2(x+1) - mean(log2(x+1))
  #             -> difference of FRECHET means, NOT log2(mean(x)+1). The two
  #             differ by Jensen and the gap grows with the gene's variance,
  #             i.e. worst exactly where it does most damage.
  if (norm_method == 0L) {
    obs_own <- case_means - ctrl_means
  } else {
    frechet <- function(m) colMeans(log2(pmax(m, 0) + 1.0), na.rm = TRUE)
    obs_own <- frechet(case_rep) - frechet(ctrl_rep)
  }
  obs_own[!is.finite(obs_own)] <- 0

  fn <- if (variant == "exact") tox_compute_noise_pvalues_pipeline_exact
        else                    tox_compute_noise_pvalues_pipeline

  res <- fn(
    case_means         = as.numeric(case_means),
    case_replicates    = case_rep,
    control_means      = as.numeric(ctrl_means),
    control_replicates = ctrl_rep,
    obs_own            = as.numeric(obs_own),
    valid_genes_own    = rep(1L, n_genes),
    norm_method        = as.integer(norm_method),
    k_start = kcfg$k_start, k_step = kcfg$k_step,
    k_max   = kcfg$k_max,   tau    = CFG$tox_tau,
    trim_frac = trim_frac,
    max_pool_size = CFG$tox_max_pool)

  # Fail loudly: a non-zero ierr alongside a silently returned p-vector would
  # look like a badly calibrated method rather than a failed call.
  if (!is.null(res$ierr) && res$ierr != 0L)
    stop(sprintf("TOX returned ierr = %d", res$ierr))

  p <- res$pvalues_own

  if (isTRUE(CFG$tox_verbose)) {
    cat(sprintf("    [%s n%d %s] n_success=%s  real p=%d\n",
                variant, norm_method, kcfg$name,
                as.character(res$n_success), sum(p >= 0)))
    print(summary(res$neighborhood_size_case))
    print(summary(p[p >= 0]))
  }

  # BH over TESTED genes only. -1 is the not-computed sentinel; including those
  # in the denominator makes every threshold too strict by the ratio of total
  # to tested, which reads as a uniformly conservative method.
  tested <- p >= 0
  padj <- rep(NA_real_, n_genes)
  padj[tested] <- p.adjust(p[tested], method = "BH")
  p[!tested] <- NA_real_

  n_success <- if (!is.null(res$n_success)) res$n_success else NA_integer_
  if (!is.na(n_success) && n_success != sum(tested))
    message(sprintf("    note: n_success=%d but %d genes carry a p-value",
                    n_success, sum(tested)))

  # Observed pool size, inferred from the realised floor. Compare against the
  # expected k_max * n_rep: a mismatch means the adaptive schedule is not
  # stopping where the configuration says it should.
  pmin_obs  <- suppressWarnings(min(p, na.rm = TRUE))
  pool_infer <- if (is.finite(pmin_obs) && pmin_obs > 0)
                  sqrt(max(0, 1/pmin_obs - 1)) else NA_real_

  data.frame(gene = rownames(tpm), stat = p, padj = padj, lfc = obs_own,
             row.names = NULL,
             frac_tested   = mean(tested),
             pool_infer    = pool_infer,
             nbhd_median   = suppressWarnings(
                               median(res$neighborhood_size_case[
                                        res$neighborhood_size_case >= 0])))
}

#' All methods on one dataset. Returns a named list of result frames.
run_all_methods <- function(sim) {
  keep <- common_filter(sim)
  counts <- sim$counts[keep, , drop = FALSE]
  tpm    <- sim$tpm[keep, , drop = FALSE]
  grp    <- sim$group

  out <- list()
  safe <- function(nm, expr) {
    r <- tryCatch(expr, error = function(e) {
      message(sprintf("    %s failed: %s", nm, e$message)); NULL })
    if (!is.null(r)) out[[nm]] <<- r
  }

  safe("edgeR",  run_edger (counts, grp))
  safe("limma",  run_limma (counts, grp))
  safe("DESeq2", run_deseq2(counts, grp))

  if (CFG$run_tox)
    for (v in CFG$tox_variants) for (nm in CFG$tox_norm) for (kg in CFG$tox_k_grid) {
      safe(sprintf("TOX_%s_n%d_%s", v, nm, kg$name),
           run_tox(tpm, grp, nm, v, kg))
      # Raw-only tail-trim companion arm (trim is a no-op under log norm).
      if (nm == 0L && CFG$tox_trim_frac > 0)
        safe(sprintf("TOX_%s_n%d_%s_trim", v, nm, kg$name),
             run_tox(tpm, grp, nm, v, kg, CFG$tox_trim_frac))
    }

  list(results = out, keep = keep, is_de = sim$is_de[keep])
}

# =============================================================================
# 3. METRICS
# =============================================================================

null_metrics <- function(res, alpha) {
  p <- res$stat[is.finite(res$stat)]
  if (length(p) < 10) return(NULL)
  f05 <- mean(p < 0.05); f01 <- mean(p < 0.01)
  data.frame(
    n_tested     = length(p),
    frac_p_lt_05 = f05,                       # target 0.05
    frac_p_lt_01 = f01,                       # target 0.01
    # SPIKE DETECTOR: target ~5 under a uniform null. Near 1 means almost every
    # sub-0.05 p-value is also sub-0.01 -- a spike at zero rather than a ramp.
    # BH lives in that spike, so this predicts n_rejected better than f05 does.
    tail_ratio   = if (f01 > 0) f05 / f01 else NA_real_,
    median_p     = median(p),                 # target 0.5
    min_p        = min(p),
    ks_p         = suppressWarnings(ks.test(p, "punif")$p.value),
    n_rejected   = sum(!is.na(res$padj) & res$padj < alpha),
    frac_tested  = if ("frac_tested" %in% names(res)) res$frac_tested[1] else NA,
    pool_infer   = if ("pool_infer" %in% names(res)) res$pool_infer[1] else NA,
    nbhd_median  = if ("nbhd_median" %in% names(res)) res$nbhd_median[1] else NA
  )
}

fdr_curve <- function(res, is_de, alphas = c(0.01, 0.05, 0.10, 0.20)) {
  do.call(rbind, lapply(alphas, function(a) {
    called <- !is.na(res$padj) & res$padj < a
    tp <- sum(called & is_de); fp <- sum(called & !is_de)
    data.frame(nominal_fdr = a,
               observed_fdr = if (tp + fp > 0) fp / (tp + fp) else NA_real_,
               tpr = tp / max(1, sum(is_de)),
               n_called = tp + fp,
               # 1 when the round produced NOTHING. Averaged over rounds this
               # separates INSTABILITY from conservatism -- and it is the
               # signature of the p-value resolution cliff.
               zero_disc = as.integer(tp + fp == 0))
  }))
}

auc_from_stat <- function(stat, is_de) {
  ok <- is.finite(stat); s <- -stat[ok]; y <- is_de[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(s)
  (sum(r[y]) - sum(y) * (sum(y) + 1) / 2) / (sum(y) * sum(!y))
}

# ---- p-value uniformity: distribution shape + tests -------------------------

#' Quantify how close a p-value vector is to Uniform(0,1).
#'
#' KS (one-sample, vs "punif") is the correct test for "is this sample uniform?".
#' NOTE on Wilcoxon: the rank-SUM (Mann-Whitney) test is TWO-sample, so it cannot
#' test a single sample against a distribution. What is reported here is the
#' one-sample signed-RANK test against mu = 0.5 (a uniform's median) -- a pure
#' LOCATION check that complements KS's distribution-shape check.
#'
#' Three one-sample goodness-of-fit tests against Uniform(0,1) -- Kolmogorov-
#' Smirnov (ks_D), Anderson-Darling (ad_A2, tail-sensitive), Cramer-von Mises
#' (cvm_W2) -- plus a Wilcoxon signed-rank LOCATION check (median = 0.5). The
#' *statistics* (ks_D/ad_A2/cvm_W2, 0 = perfectly uniform) are the effect sizes to
#' read; the p-values collapse to ~0 at the n we have here even for trivially small
#' deviations. p is clamped off {0,1} for AD (which diverges at the boundary).
uniformity_stats <- function(p) {
  p <- p[is.finite(p)]
  n <- length(p)
  if (n < 10)
    return(data.frame(n = n, median_p = NA_real_, ks_D = NA_real_, ks_p = NA_real_,
                      ad_A2 = NA_real_, ad_p = NA_real_, cvm_W2 = NA_real_, cvm_p = NA_real_,
                      wilx_p = NA_real_))
  pc  <- pmin(pmax(p, 1e-9), 1 - 1e-9)
  ks  <- suppressWarnings(ks.test(p, "punif"))
  ad  <- suppressWarnings(tryCatch(goftest::ad.test(pc, "punif"),  error = function(e) NULL))
  cvm <- suppressWarnings(tryCatch(goftest::cvm.test(pc, "punif"), error = function(e) NULL))
  wx  <- suppressWarnings(tryCatch(wilcox.test(p, mu = 0.5)$p.value, error = function(e) NA_real_))
  data.frame(n = n, median_p = median(p),
             ks_D  = unname(ks$statistic),  ks_p  = ks$p.value,
             ad_A2 = if (!is.null(ad))  unname(ad$statistic)  else NA_real_,
             ad_p  = if (!is.null(ad))  ad$p.value            else NA_real_,
             cvm_W2 = if (!is.null(cvm)) unname(cvm$statistic) else NA_real_,
             cvm_p  = if (!is.null(cvm)) cvm$p.value           else NA_real_,
             wilx_p = wx)
}

#' 20-bin unicode sparkline of a p-value vector -- a one-line histogram.
spark <- function(p, nbins = 20L) {
  p <- p[is.finite(p)]
  blocks <- c("▁","▂","▃","▄","▅","▆","▇","█")
  if (!length(p)) return(strrep(" ", nbins))
  h  <- hist(pmin(pmax(p, 0), 1), breaks = seq(0, 1, length.out = nbins + 1L),
             plot = FALSE)$counts
  mx <- max(h)
  if (mx == 0) return(strrep(blocks[1], nbins))
  paste(blocks[pmax(1L, ceiling(h / mx * length(blocks)))], collapse = "")
}

#' Multi-line ASCII histogram of a p-value vector (safe in any terminal).
text_hist <- function(p, nbins = 20L, width = 48L) {
  p <- p[is.finite(p)]
  if (!length(p)) return(character(0))
  br <- seq(0, 1, length.out = nbins + 1L)
  h  <- hist(pmin(pmax(p, 0), 1), breaks = br, plot = FALSE)$counts
  mx <- max(h, 1L)
  vapply(seq_len(nbins), function(i)
    sprintf("  %.2f |%-*s %d", br[i], width, strrep("#", round(h[i] / mx * width)), h[i]),
    character(1))
}

# =============================================================================
# 4. PART A -- MOCK NULL
# =============================================================================

part_a_mock_null <- function() {
  message("\n=== PART A: MOCK NULL (no DE genes) ===")
  set.seed(CFG$seed)
  truth <- make_truth(CFG$n_genes, n_de = 0)

  rows <- list(); pvals <- list(); k <- 1L
  for (d in CFG$dists) for (n in CFG$n_reps) for (N in CFG$lib_sizes) {
    for (i in seq_len(CFG$n_rounds)) {
      message(sprintf("  [%s] n_rep=%d lib=%.0e round %d/%d",
                      d, n, N, i, CFG$n_rounds))
      sim <- simulate_dataset(truth, d, n, N)
      rr  <- run_all_methods(sim)
      for (m in names(rr$results)) {
        mt <- null_metrics(rr$results[[m]], CFG$alpha)
        if (is.null(mt)) next
        rows[[k]] <- cbind(dist = d, n_rep = n, lib = N, round = i,
                           method = m, mt)
        if (i == 1L) pvals[[paste(d, n, m)]] <-
          data.frame(dist = d, n_rep = n, method = m,
                     p = rr$results[[m]]$stat[is.finite(rr$results[[m]]$stat)])
        k <- k + 1L
      }
    }
  }
  list(metrics = do.call(rbind, rows), pvals = do.call(rbind, pvals))
}

# =============================================================================
# 5. PART B -- FDR CALIBRATION WITH SIGNAL
# =============================================================================

part_b_fdr <- function() {
  message("\n=== PART B: FDR CALIBRATION (with DE genes) ===")
  set.seed(CFG$seed + 1)
  truth <- make_truth(CFG$n_genes, n_de = CFG$n_de_partB)

  rows <- list(); k <- 1L
  for (d in CFG$dists) for (n in CFG$n_reps) for (N in CFG$lib_sizes) {
    for (i in seq_len(max(1L, CFG$n_rounds %/% 2L))) {
      message(sprintf("  [%s] n_rep=%d lib=%.0e round %d", d, n, N, i))
      sim <- simulate_dataset(truth, d, n, N)
      rr  <- run_all_methods(sim)
      for (m in names(rr$results)) {
        res <- rr$results[[m]]
        fc  <- fdr_curve(res, rr$is_de)
        fc$auc <- auc_from_stat(res$stat, rr$is_de)
        rows[[k]] <- cbind(dist = d, n_rep = n, lib = N, round = i,
                           method = m, fc)
        k <- k + 1L
      }
    }
  }
  do.call(rbind, rows)
}

# =============================================================================
# 6. REPORT
# =============================================================================

#' aggregate() with na.rm, so that ONE round with no discoveries does not
#' blank the whole cell. The count of such rounds is reported separately as
#' n_zero_disc -- conflating the two hides instability behind an NA.
agg_narm <- function(formula, data)
  aggregate(formula, data, function(x) mean(x, na.rm = TRUE), na.action = na.pass)

report <- function(A, B) {
  message("\n", strrep("=", 78))
  message("PART A -- NULL CALIBRATION (mean over rounds)")
  message(strrep("=", 78))
  agg <- agg_narm(cbind(frac_p_lt_05, frac_p_lt_01, tail_ratio, median_p,
                        min_p, n_rejected, frac_tested,
                        pool_infer, nbhd_median) ~ method + dist + n_rep,
                  A$metrics)
  agg <- agg[order(agg$dist, agg$n_rep, agg$method), ]
  num <- sapply(agg, is.numeric); agg[num] <- round(agg[num], 4)
  print(agg, row.names = FALSE)
  message("\n  tail_ratio  target ~5.  Near 1 = SPIKE AT ZERO: the deep tail is")
  message("              broken even when frac_p_lt_05 looks correct, and BH")
  message("              reads the deep tail. Check this BEFORE frac_p_lt_05.")
  message("  frac_p_lt_05 target 0.05.  >> too NARROW (anti-conservative);")
  message("              << too WIDE (conservative -> low recall).")
  message("  n_rejected  target ~0. Under a COMPLETE null BH controls FWER,")
  message("              so any rejection is a failure.")
  message("  pool_infer  pool size implied by the realised floor, = ",
          "sqrt(1/min_p - 1).")
  message("              Compare with k_max * n_rep; a mismatch means the")
  message("              adaptive schedule is not stopping where configured.")

  # ---- p-value uniformity (KS / Wilcoxon). The histograms themselves are
  #      saved to PNG (null_pvalue_histograms.png), not printed here. ---------
  message("\n", strrep("=", 78))
  message("PART A -- P-VALUE UNIFORMITY (mock null, round 1; vs Uniform(0,1))")
  message(strrep("=", 78))
  us <- do.call(rbind, lapply(
    split(A$pvals, list(A$pvals$method, A$pvals$dist, A$pvals$n_rep), drop = TRUE),
    function(df) data.frame(method = df$method[1], dist = df$dist[1], n_rep = df$n_rep[1],
                            uniformity_stats(df$p), stringsAsFactors = FALSE)))
  us <- us[order(us$dist, us$n_rep, us$method), ]
  numc <- c("median_p", "ks_D", "ks_p", "ad_A2", "ad_p", "cvm_W2", "cvm_p", "wilx_p")
  usp <- us; usp[numc] <- lapply(usp[numc], function(x) signif(x, 3))
  print(usp, row.names = FALSE)
  message("\n  ks_D / ad_A2 / cvm_W2  distance of the p-values to Uniform(0,1) -- Kolmogorov-")
  message("         Smirnov, Anderson-Darling (tail-sensitive) and Cramer-von Mises. 0 = uniform.")
  message("         These STATISTICS are the effect sizes -- read them, not the p-values.")
  message("  *_p    the matching GoF p-values; at this n they collapse to ~0 for even tiny")
  message("         deviations, so they are nearly always negligible; judge by the statistics.")
  message("  wilx_p Wilcoxon signed-rank p (H0: median = 0.5) -- a LOCATION check only.")
  message("         (The rank-SUM test is two-sample and cannot test one sample for")
  message("          uniformity; KS/AD/CvM are the right tools. See uniformity_stats() doc.)")
  message("  Histograms + QQ-plots are written to null_pvalue_{histograms,qqplots}.png.")

  message("\n", strrep("=", 78))
  message("PART B -- OBSERVED vs NOMINAL FDR")
  message(strrep("=", 78))
  aggB <- agg_narm(cbind(observed_fdr, tpr, auc, n_called, zero_disc) ~
                     method + dist + nominal_fdr, B)
  names(aggB)[names(aggB) == "zero_disc"] <- "frac_zero_disc"
  aggB <- aggB[order(aggB$dist, aggB$nominal_fdr, aggB$method), ]
  num <- sapply(aggB, is.numeric); aggB[num] <- round(aggB[num], 4)
  print(aggB, row.names = FALSE)
  message("\n  observed_fdr should sit AT OR BELOW nominal_fdr.")
  message("  frac_zero_disc = fraction of rounds returning NOTHING. Non-zero")
  message("              means UNSTABLE, not conservative -- and it is the")
  message("              signature of the p-value resolution cliff.")
  message("  auc vs tpr: high auc + low tpr -> ranking fine, THRESHOLD wrong;")
  message("              low auc -> the ranking itself is broken.")

  # ---- neighbourhood sweep summary ---------------------------------------
  if (CFG$run_tox && any(grepl("^TOX_", aggB$method))) {
    message("\n", strrep("=", 78))
    message("NEIGHBOURHOOD SWEEP -- the k trade-off")
    message(strrep("=", 78))
    sw <- aggB[grepl("^TOX_", aggB$method) & aggB$nominal_fdr == CFG$alpha,
               c("method", "dist", "auc", "tpr", "observed_fdr", "frac_zero_disc")]
    swA <- agg[grepl("^TOX_", agg$method),
               c("method", "dist", "tail_ratio", "n_rejected", "pool_infer")]
    swA <- aggregate(cbind(tail_ratio, n_rejected, pool_infer) ~ method + dist,
                     swA, mean, na.action = na.pass)
    m <- merge(sw, swA, by = c("method", "dist"), all.x = TRUE)
    m <- m[order(m$dist, m$method), ]
    num <- sapply(m, is.numeric); m[num] <- round(m[num], 4)
    print(m, row.names = FALSE)
    message("\n  Larger k -> finer resolution and a more stable null, but a")
    message("  wider mean range in the neighbourhood, so weaker local")
    message("  exchangeability and worse AUC. Look for the knee: the largest k")
    message("  whose AUC has not yet started falling.")
  }

  # ---- plots -------------------------------------------------------------
  n_meth <- length(unique(A$pvals$method))
  ph <- ggplot(A$pvals, aes(p)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey40") +
    facet_grid(method ~ dist + n_rep, scales = "free_y") +
    labs(title = "Mock null: p-value histograms (must be flat)",
         x = "p-value", y = "count") + theme_bw(base_size = 7)
  ggsave(file.path(CFG$out_dir, "null_pvalue_histograms.png"), ph,
         width = 11, height = max(9, 0.6 * n_meth), dpi = 150, limitsize = FALSE)

  # QQ-plots: observed p vs expected uniform quantile (on the diagonal = calibrated)
  qq_points <- function(p, npts = 2000L) {
    p <- sort(p[is.finite(p)]); n <- length(p)
    if (n < 2L) return(NULL)
    idx <- unique(round(seq(1, n, length.out = min(npts, n))))
    data.frame(expected = (idx - 0.5) / n, observed = p[idx])
  }
  qq <- do.call(rbind, lapply(
    split(A$pvals, list(A$pvals$method, A$pvals$dist, A$pvals$n_rep), drop = TRUE),
    function(df) { q <- qq_points(df$p); if (is.null(q)) NULL else
                   cbind(method = df$method[1], dist = df$dist[1], n_rep = df$n_rep[1], q) }))
  if (!is.null(qq) && nrow(qq)) {
    pq <- ggplot(qq, aes(expected, observed)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(size = 0.25, alpha = 0.5, colour = "grey20") +
      facet_grid(method ~ dist + n_rep) + coord_equal() +
      labs(title = "Mock null: p-value QQ-plots (observed vs Uniform(0,1); diagonal = calibrated)",
           subtitle = "bowing ABOVE the diagonal at small expected = anti-conservative (spike at zero)",
           x = "expected p (uniform quantile)", y = "observed p") + theme_bw(base_size = 7)
    ggsave(file.path(CFG$out_dir, "null_pvalue_qqplots.png"), pq,
           width = 11, height = max(9, 0.6 * n_meth), dpi = 150, limitsize = FALSE)
  }

  pf <- ggplot(aggB, aes(nominal_fdr, observed_fdr, colour = method)) +
    geom_abline(linetype = 2, colour = "grey50") +
    geom_line() + geom_point() + facet_wrap(~ dist) +
    labs(title = "Observed vs nominal FDR (below the diagonal = controlled)",
         x = "nominal FDR", y = "observed FDR") + theme_bw(base_size = 9)
  ggsave(file.path(CFG$out_dir, "fdr_calibration.png"), pf,
         width = 11, height = 8, dpi = 150)

  if (CFG$run_tox) {
    ks <- aggB[grepl("^TOX_", aggB$method) & aggB$nominal_fdr == CFG$alpha, ]
    if (nrow(ks)) {
      ks$kcfg    <- sub("^TOX_[a-z]+_n[01]_", "", ks$method)
      ks$variant <- sub("^TOX_([a-z]+)_n([01])_.*$", "\\1 n\\2", ks$method)
      pk <- ggplot(ks, aes(kcfg, auc, colour = variant, group = variant)) +
        geom_line() + geom_point() + facet_wrap(~ dist) +
        labs(title = "AUC vs neighbourhood size (look for the knee)",
             x = "k configuration", y = "AUC") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(CFG$out_dir, "k_sweep_auc.png"), pk,
             width = 10, height = 7, dpi = 150)
    }
  }

  write.csv(A$metrics, file.path(CFG$out_dir, "null_metrics.csv"), row.names = FALSE)
  write.csv(B,        file.path(CFG$out_dir, "fdr_calibration.csv"), row.names = FALSE)
  message("\nWrote CSVs and plots to ", CFG$out_dir, "/")
}

# =============================================================================
if (sys.nframe() == 0L) {
  message("TOX enabled: ", CFG$run_tox)
  n_arms <- 3 + if (CFG$run_tox) {
    base <- length(CFG$tox_variants) * length(CFG$tox_norm) * length(CFG$tox_k_grid)
    # raw-only "_trim" companion arms (one extra per variant x k_grid x raw norm)
    trim <- if (CFG$tox_trim_frac > 0)
      length(CFG$tox_variants) * sum(CFG$tox_norm == 0L) * length(CFG$tox_k_grid) else 0
    base + trim
  } else 0
  message(sprintf("Arms per dataset: %d  |  datasets: %d (A) + %d (B)",
                  n_arms,
                  length(CFG$dists) * length(CFG$n_reps) *
                    length(CFG$lib_sizes) * CFG$n_rounds,
                  length(CFG$dists) * length(CFG$n_reps) *
                    length(CFG$lib_sizes) * max(1L, CFG$n_rounds %/% 2L)))
  message("  Trim CFG$tox_variants / tox_norm / tox_k_grid to cut runtime.")
  resolution_report()

  A <- part_a_mock_null()
  B <- part_b_fdr()
  report(A, B)
}