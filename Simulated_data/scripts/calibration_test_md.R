#!/usr/bin/env Rscript
# calibration_test_md.R
# =============================================================================
# NULL CALIBRATION AND POWER TEST for the MULTI-AXIS noise model (noise_model_md)
#
#   Rscript calibration_test_md.R
#
# The multi-axis counterpart of calibration_test.R, and built on it: that script
# is sourced for the simulator, the metrics and the uniformity tests, so the two
# arms share one observation model and one definition of "calibrated".
#
# DESIGN
#   d = 1, 2, 3 independent axes (tissues). Every tissue has its own samples (the
#   model's scope assumption) and its own baseline abundance profile. A gene is
#   ONE test across all d tissues.
#
#   Part A -- MOCK NULL. No shift in any tissue. p must be Uniform(0,1); read
#     tail_ratio, frac_p_lt_05, n_rejected exactly as in calibration_test.R.
#
#   Part B -- POWER. The same DE genes carry the same log2 shift in n_shift of
#     the d tissues (n_shift = 1..d) and none in the rest. A gene is a true
#     positive if it is shifted in ANY tissue. The question is how power decays
#     when the signal sits in fewer axes than are tested: one shifted tissue out
#     of three carries the same signal as d = 1 but pays for two noise axes.
#
# COMPARATOR
#   edgeR QL F-test of all tissue:group interaction coefficients
#   (~ tissue + tissue:group), i.e. "shift in any tissue" with d df -- the
#   count-model analogue of the D_std test. At d = 1 it reduces to ~ group.
#
# d = 1 runs the enumerated path and reproduces noise_model_exact; d >= 2 is
# sampled (systematic stride, M = n_draws_max), so its p-value floor is
# 1 / (M + 1) whatever the neighbourhood size.
# =============================================================================

# ---- shared machinery from the scalar test (its main block does not run when sourced)
SCRIPT_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand <- c(if (length(f)) dirname(normalizePath(f[1])),
            "experiments/Noise_Model_Test/Simulated_data/scripts", ".")
  hit <- cand[file.exists(file.path(cand, "calibration_test.R"))]
  if (!length(hit)) stop("calibration_test.R not found; looked in: ", paste(cand, collapse = ", "))
  normalizePath(hit[1])
})
source(file.path(SCRIPT_DIR, "calibration_test.R"))

if (!CFG$run_tox) stop("TOX wrapper unavailable -- the MD test has nothing to test.")

MD_CFG <- list(
  n_axes      = 1:3,
  n_reps      = CFG$n_reps,               # 3 and 5
  dists       = c("nb", "lnpois"),
  n_genes     = CFG$n_genes,
  n_de        = CFG$n_de_partB,
  rounds_null = 5,
  rounds_pow  = 3,
  alpha       = CFG$alpha,
  lib_size    = CFG$lib_sizes[1],
  tissue_sd   = 0.5,    # sd of the per-tissue log-abundance perturbation
  k           = CFG$tox_k_grid[[1]],
  n_draws_max = 20000L,
  seed        = 20240202,
  out_dir     = paste0(CFG$out_dir, "_md")
)
dir.create(MD_CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. SIMULATION
# =============================================================================

#' One tissue's truth in the layout simulate_dataset() expects. Same
#' renormalisation as make_truth(): only the DE set is rescaled, so null genes
#' (and every gene of an unshifted tissue) are exactly null.
tissue_truth <- function(pi_a, de, beta) {
  pi_a <- pi_a / sum(pi_a)
  pi_b <- pi_a
  if (length(de)) {
    cs <- sum(pi_a[de]) / sum(pi_a[de] * 2^beta[de])
    pi_b[de] <- pi_a[de] * 2^beta[de] * cs
  }
  list(pi_a = pi_a, pi_b = pi_b, true_lfc = log2(pi_b / pi_a),
       is_de = seq_along(pi_a) %in% de)
}

#' d tissues; the DE genes are shifted in tissues 1..n_shift (the axes are
#' exchangeable, so which ones does not matter). One sim per tissue, with
#' tissue-prefixed sample names so the axes share no sample IDs.
simulate_md <- function(pi_base, de, beta, d, n_shift, dist, n_rep) {
  lapply(seq_len(d), function(t) {
    pi_t <- pi_base * exp(rnorm(length(pi_base), 0, MD_CFG$tissue_sd))
    sim  <- simulate_dataset(tissue_truth(pi_t, if (t <= n_shift) de else integer(0), beta),
                             dist, n_rep, MD_CFG$lib_size)
    colnames(sim$counts) <- colnames(sim$tpm) <- paste0("T", t, "_", colnames(sim$counts))
    sim
  })
}

# =============================================================================
# 2. METHODS  (same return layout as calibration_test.R: gene, stat, padj, lfc)
# =============================================================================

run_tox_md <- function(sims, keep) {
  side <- function(g) lapply(sims, function(s) t(s$tpm[keep, s$group == g, drop = FALSE]))
  res <- tox_compute_noise_pvalues_pipeline_md(
    case_replicates    = side("B"),
    control_replicates = side("A"),
    norm_method        = 1L,
    k_start = MD_CFG$k$k_start, k_step = MD_CFG$k$k_step,
    k_max   = MD_CFG$k$k_max,   tau    = CFG$tox_tau,
    n_draws_max   = MD_CFG$n_draws_max,
    max_pool_size = CFG$tox_max_pool)
  if (!is.null(res$ierr) && res$ierr != 0L)
    stop(sprintf("TOX md returned ierr = %d", res$ierr))

  p <- res$pvalues_own
  tested <- p >= 0                                  # -1 = not computed
  padj <- rep(NA_real_, length(p))
  padj[tested] <- p.adjust(p[tested], method = "BH")
  p[!tested] <- NA_real_
  nb <- res$neighborhood_size_case
  data.frame(gene = rownames(sims[[1]]$tpm)[keep], stat = p, padj = padj,
             lfc = res$d_obs, row.names = NULL,
             # per-axis effect size, one column per tissue: log2fc_T1, log2fc_T2, ...
             setNames(as.data.frame(t(res$log2fc)), paste0("log2fc_T", seq_along(sims))),
             frac_tested = mean(tested),
             nbhd_median = suppressWarnings(median(nb[nb >= 0])))
}

run_edger_md <- function(sims, keep) {
  counts <- do.call(cbind, lapply(sims, function(s) s$counts[keep, , drop = FALSE]))
  group  <- factor(unlist(lapply(sims, `[[`, "group")))
  tissue <- factor(rep(seq_along(sims), vapply(sims, function(s) ncol(s$counts), 1L)))
  design <- if (nlevels(tissue) == 1L) model.matrix(~ group)
            else model.matrix(~ tissue + tissue:group)
  coefs  <- grep("group", colnames(design))
  y   <- normLibSizes(DGEList(counts = counts, group = group))
  fit <- glmQLFit(y, design, robust = TRUE)
  tt  <- topTags(glmQLFTest(fit, coef = coefs), n = Inf, sort.by = "none")$table
  data.frame(gene = rownames(counts), stat = tt$PValue, padj = tt$FDR,
             lfc = NA_real_, row.names = NULL)
}

run_all_md <- function(sims) {
  keep <- Reduce(`&`, lapply(sims, common_filter))   # expressed in EVERY tissue
  out <- list()
  safe <- function(nm, expr) {
    r <- tryCatch(expr, error = function(e) {
      message(sprintf("    %s failed: %s", nm, e$message)); NULL })
    if (!is.null(r)) out[[nm]] <<- r
  }
  safe("edgeR_F", run_edger_md(sims, keep))
  safe("TOX_md",  run_tox_md  (sims, keep))
  list(results = out, keep = keep)
}

# =============================================================================
# 3. PARTS
# =============================================================================

part_a_md <- function() {
  message("\n=== MD PART A: MOCK NULL ===")
  set.seed(MD_CFG$seed)
  pi_base <- rlnorm(MD_CFG$n_genes, 0, 2)
  rows <- list(); pvals <- list()
  for (dist in MD_CFG$dists) for (d in MD_CFG$n_axes) for (n in MD_CFG$n_reps)
    for (i in seq_len(MD_CFG$rounds_null)) {
      message(sprintf("  [%s] d=%d n_rep=%d round %d/%d", dist, d, n, i, MD_CFG$rounds_null))
      rr <- run_all_md(simulate_md(pi_base, integer(0), 0, d, 0L, dist, n))
      for (m in names(rr$results)) {
        res <- rr$results[[m]]
        mt  <- null_metrics(res, MD_CFG$alpha)
        if (is.null(mt)) next
        rows[[length(rows) + 1L]] <- cbind(dist = dist, d = d, n_rep = n, round = i,
                                           method = m, mt)
        if (i == 1L) pvals[[length(pvals) + 1L]] <-
          data.frame(dist = dist, d = d, n_rep = n, method = m,
                     p = res$stat[is.finite(res$stat)])
      }
    }
  list(metrics = do.call(rbind, rows), pvals = do.call(rbind, pvals))
}

part_b_md <- function() {
  message("\n=== MD PART B: POWER (shift in n_shift of d tissues) ===")
  set.seed(MD_CFG$seed + 1)
  G <- MD_CFG$n_genes
  pi_base <- rlnorm(G, 0, 2)
  de   <- sample.int(G, MD_CFG$n_de)
  beta <- numeric(G)                                 # same as make_truth()'s draw
  beta[de] <- ifelse(runif(length(de)) < 0.5, 1, -1) * abs(rnorm(length(de), 1.0, 0.6))
  is_de <- seq_len(G) %in% de

  rows <- list()
  for (dist in MD_CFG$dists) for (d in MD_CFG$n_axes) for (s in seq_len(d))
    for (n in MD_CFG$n_reps) for (i in seq_len(MD_CFG$rounds_pow)) {
      message(sprintf("  [%s] d=%d n_shift=%d n_rep=%d round %d", dist, d, s, n, i))
      rr <- run_all_md(simulate_md(pi_base, de, beta, d, s, dist, n))
      truth <- is_de[rr$keep]
      for (m in names(rr$results)) {
        res <- rr$results[[m]]
        fc  <- fdr_curve(res, truth)
        fc$auc <- auc_from_stat(res$stat, truth)
        rows[[length(rows) + 1L]] <- cbind(dist = dist, d = d, n_shift = s, n_rep = n,
                                           round = i, method = m, fc)
      }
    }
  do.call(rbind, rows)
}

# =============================================================================
# 4. REPORT
# =============================================================================

report_md <- function(A, B) {
  rnd <- function(df) { num <- sapply(df, is.numeric); df[num] <- round(df[num], 4); df }

  message("\n", strrep("=", 78), "\nMD PART A -- NULL CALIBRATION (mean over rounds)\n", strrep("=", 78))
  agg <- agg_narm(cbind(frac_p_lt_05, frac_p_lt_01, tail_ratio, median_p, min_p,
                        n_rejected, frac_tested, nbhd_median) ~ method + dist + d + n_rep,
                  A$metrics)
  print(rnd(agg[order(agg$dist, agg$d, agg$n_rep, agg$method), ]), row.names = FALSE)
  message("  targets: tail_ratio ~5, frac_p_lt_05 0.05, n_rejected 0 (see calibration_test.R).")
  message(sprintf("  d >= 2 is sampled: p floor = 1/(M+1) = %.1e.", 1 / (MD_CFG$n_draws_max + 1)))

  us <- do.call(rbind, lapply(
    split(A$pvals, list(A$pvals$method, A$pvals$dist, A$pvals$d, A$pvals$n_rep), drop = TRUE),
    function(df) data.frame(method = df$method[1], dist = df$dist[1], d = df$d[1],
                            n_rep = df$n_rep[1], uniformity_stats(df$p))))
  message("\nUniformity (round 1): ks_D / ad_A2 / cvm_W2 -- 0 = uniform")
  print(rnd(us[order(us$dist, us$d, us$n_rep, us$method),
               c("method", "dist", "d", "n_rep", "n", "median_p", "ks_D", "ad_A2", "cvm_W2")]),
        row.names = FALSE)

  message("\n", strrep("=", 78), "\nMD PART B -- POWER AND FDR at nominal ", MD_CFG$alpha,
          "\n", strrep("=", 78))
  aggB <- agg_narm(cbind(tpr, observed_fdr, auc, n_called, zero_disc) ~
                     method + dist + d + n_shift + n_rep + nominal_fdr, B)
  names(aggB)[names(aggB) == "zero_disc"] <- "frac_zero_disc"
  at <- aggB[aggB$nominal_fdr == MD_CFG$alpha, ]
  print(rnd(at[order(at$dist, at$n_rep, at$d, at$n_shift, at$method), ]), row.names = FALSE)
  message("  tpr should rise with n_shift at fixed d, and fall with d at fixed n_shift = 1")
  message("  (same signal, more noise axes). observed_fdr must stay <= nominal throughout.")

  # ---- plots
  qq <- do.call(rbind, lapply(
    split(A$pvals, list(A$pvals$method, A$pvals$dist, A$pvals$d, A$pvals$n_rep), drop = TRUE),
    function(df) { p <- sort(df$p); n <- length(p); idx <- unique(round(seq(1, n, length.out = min(2000, n))))
                   data.frame(df[1, c("method", "dist", "d", "n_rep")], row.names = NULL,
                              expected = (idx - 0.5) / n, observed = p[idx]) }))
  pq <- ggplot(qq, aes(expected, observed, colour = factor(d))) +
    geom_abline(linetype = 2, colour = "grey50") + geom_line() +
    facet_grid(method ~ dist + n_rep, labeller = label_both) + coord_equal() +
    labs(title = "MD mock null: p-value QQ (diagonal = calibrated)", colour = "d axes",
         x = "expected p", y = "observed p") + theme_bw(base_size = 8)
  ggsave(file.path(MD_CFG$out_dir, "md_null_qqplots.png"), pq, width = 10, height = 6, dpi = 150)

  pp <- ggplot(at, aes(n_shift, tpr, colour = method, linetype = factor(n_rep))) +
    geom_line() + geom_point() + facet_grid(dist ~ d, labeller = label_both) +
    scale_x_continuous(breaks = 1:3) +
    labs(title = sprintf("MD power at BH %.2f vs number of shifted tissues", MD_CFG$alpha),
         x = "shifted tissues", y = "TPR", linetype = "n_rep") + theme_bw(base_size = 9)
  ggsave(file.path(MD_CFG$out_dir, "md_power.png"), pp, width = 10, height = 6, dpi = 150)

  write.csv(A$metrics, file.path(MD_CFG$out_dir, "md_null_metrics.csv"), row.names = FALSE)
  write.csv(B,         file.path(MD_CFG$out_dir, "md_power_fdr.csv"),    row.names = FALSE)
  message("\nWrote CSVs and plots to ", MD_CFG$out_dir, "/")
}

# =============================================================================
if (sys.nframe() == 0L) {
  A <- part_a_md()
  B <- part_b_md()
  report_md(A, B)
}
