# power_metrics.R
# -----------------------------------------------------------------------------
# Pure metric functions for the power benchmark. Sourced by power_report.R; kept
# in a file of its own so it can be unit-tested (TCGA_test/tests/) without the
# Fortran build or any reference method. No package dependency beyond base R.
#
# Conventions
#   p      raw p-value per gene, NA = untested (ranked last, tied)
#   padj   BH-adjusted p (NA = untested)
#   tie    |log2 FC of group means|: the tie-break of the "_pe" ranking
#   is_de  logical truth over the SAME genes (the common_filter universe)
# -----------------------------------------------------------------------------

PM <- list(
  alphas   = c(0.01, 0.05, 0.10),
  fdr_grid = seq(0, 0.30, by = 0.01),
  pr_grid  = seq(0, 1, by = 0.05),
  topk     = c(100L, 250L, 500L),
  lfc_bins = c(0.25, 0.5, 1, 2, 4)
)

# Cell keys: one combination of simulation parameters. Keys unused by a part
# hold a constant (never cycled), so they group harmlessly.
CELL <- c("part", "dist", "n_rep", "pi1", "het", "frac_up", "depth", "disp_model",
          "hv_frac", "sigma_d", "delta_target", "cohort")

# ------------------------------------------------------------------ ranking

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

#' Step-function AUPRC (average precision): sum_k (R_k - R_{k-1}) * P_k over
#' tie-block ends. No interpolation (not valid in PR space). `denom` = the recall
#' denominator: r$P (tested universe) or n_de_total (filtered-out DE never found).
auprc <- function(r, denom = r$P) {
  prec <- r$tp / pmax(r$tp + r$fp, 1)
  sum(diff(r$tp) * prec[-1]) / denom
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

# ------------------------------------------------------------------ per dataset

#' Kolmogorov-Smirnov distance of p-values to U(0, 1): sup |F_n(x) - x| (exact on
#' the empirical CDF; ties need no correction because the sup is taken at the jumps).
ks_unif <- function(p) {
  p <- sort(p[is.finite(p)]); n <- length(p)
  if (!n) return(NA_real_)
  max(seq_len(n) / n - p, p - (seq_len(n) - 1) / n)
}

#' Null-calibration set (note D WP4) over the true nulls (all genes when pi1 = 0).
#' Untested nulls (NA p) are excluded and counted.
null_calibration <- function(p, is_null) {
  pn <- p[is_null]; ok <- !is.na(pn); x <- pn[ok]
  f <- function(a) if (length(x)) mean(x <= a) else NA_real_
  list(p_le01_null = f(0.01), p_le05_null = f(0.05), p_le10_null = f(0.10),
       ks_null = ks_unif(x), n_p_na_null = sum(!ok))
}

#' All scalar metrics for one method on one dataset.
#'   res: data.frame(p, padj, lfc, tie, floor_th[, nbhd_own_case, nbhd_own_control])
#'   lfc_log2: is res$lfc on a log2 scale (FALSE for TOX-raw)
#' Without DE genes (pi1 = 0 null arms) every metric that needs them (TPR, AUC,
#' AUPRC, precision at k, LFC metrics) is NA; the calibration metrics are computed.
metrics_one <- function(res, is_de, true_lfc, n_de_total, is_hv = NULL, lfc_log2 = TRUE,
                        alphas = PM$alphas, topk = PM$topk) {
  p <- res$p; padj <- res$padj
  tested <- !is.na(p)
  r_p  <- roc_steps(-p, is_de)
  r_pe <- roc_steps(score_pe(p, res$tie), is_de)
  P <- sum(is_de)
  row <- list(n_genes = length(p), n_tested = sum(tested), n_de = P, n_de_total = n_de_total,
              n_padj_na = sum(tested & is.na(padj)))

  for (a in alphas) {
    called <- !is.na(padj) & padj < a
    tp <- sum(called & is_de); fp <- sum(called & !is_de)
    s <- sprintf("%02d", round(100 * a))
    row[[paste0("n_called_", s)]] <- tp + fp
    # I1: zero-filled FDP, so the mean over rounds is E[V / max(R, 1)] (the FDR BH
    # controls). fdr_cond_* is the old conditional E[V / R | R > 0], NA without calls.
    row[[paste0("fdr_nom_", s)]]  <- fp / max(1, tp + fp)
    row[[paste0("fdr_cond_", s)]] <- if (tp + fp > 0) fp / (tp + fp) else NA_real_
    row[[paste0("tpr_nom_", s)]]  <- tp / max(1, P)
    row[[paste0("tpr_all_nom_", s)]] <- tp / max(1, n_de_total)
    row[[paste0("zero_disc_", s)]] <- as.integer(tp + fp == 0)
    if (P == 0) row[[paste0("tpr_nom_", s)]] <- row[[paste0("tpr_all_nom_", s)]] <- NA_real_
  }
  # P(at least one BH call); under pi1 = 0 this is P(any false call) = the FDR.
  row$any_call_05 <- 1 - row$zero_disc_05
  called05 <- !is.na(padj) & padj < 0.05
  row$fpr_null05 <- sum(called05 & !is_de) / max(1, sum(!is_de))
  # Unadjusted calibration of the true nulls (the P+ / D / RS null arms): P(p <= a | null).
  row <- c(row, null_calibration(p, !is_de))
  row$fpr_hv05 <- if (any(is_hv)) mean(called05[is_hv]) else NA_real_
  row$fpr_nonhv05 <- if (any(is_hv)) sum(called05 & !is_de & !is_hv) / max(1, sum(!is_de & !is_hv))
                     else NA_real_

  k05 <- ach_step(r_p, 0.05)
  row$tpr_ach05     <- if (k05 > 0) r_p$tp[k05] / max(1, P) else 0
  row$tpr_all_ach05 <- if (k05 > 0) r_p$tp[k05] / max(1, n_de_total) else 0
  if (P == 0) row$tpr_ach05 <- row$tpr_all_ach05 <- NA_real_

  for (nm in c("p", "pe")) {
    r <- if (nm == "p") r_p else r_pe
    row[[paste0("auc_", nm)]]    <- roc_area(r)
    row[[paste0("pauc05_", nm)]] <- pauc_std(r, 0.05)
    row[[paste0("pauc10_", nm)]] <- pauc_std(r, 0.10)
    row[[paste0("auprc_", nm)]]     <- auprc(r)
    row[[paste0("auprc_all_", nm)]] <- auprc(r, n_de_total)
  }
  # Random-ranking baseline. The ranked universe is every gene in `res` (untested
  # ones included, last and tied), so an all-tied ranking scores exactly P / n_genes.
  row$prevalence <- P / length(p)
  row$auprc_excess_pe <- row$auprc_pe - row$prevalence
  row$auprc_lift_pe   <- row$auprc_pe / row$prevalence
  spe <- score_pe(p, res$tie); o <- order(spe, decreasing = TRUE, na.last = TRUE)
  for (k in topk) row[[paste0("prec_top", k)]] <- if (P > 0) mean(is_de[o[seq_len(min(k, length(o)))]]) else NA_real_
  row$prec_top_nde <- if (P > 0) mean(is_de[o[seq_len(P)]]) else NA_real_   # R-precision

  cd <- called05 & is_de
  row$sign_err05 <- if (any(cd)) mean(sign(res$lfc[cd]) != sign(true_lfc[cd])) else NA_real_
  row$lfc_rmse <- if (lfc_log2 && any(tested & is_de))
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
  # TOX pool size in residuals (section 5 of the plan; NA for the reference tools).
  for (nb in c("nbhd_own_case", "nbhd_own_control"))
    row[[paste0(nb, "_median")]] <- if (!is.null(res[[nb]]) && any(is.finite(res[[nb]])))
                                      median(res[[nb]], na.rm = TRUE) else NA_real_
  # Variance-band arms (implementation notes A/B; vb_* columns, NA for every other
  # method). New columns only -- no existing metric changes.
  #   frac_short_pool   tested genes whose smaller side pool holds < k genes
  #   n_scanned_median  candidate genes examined per (case-side) pool
  #   frac_skipped      kept genes without a p-value (band pool < 10 genes, etc.)
  #   tau_var, d0, d0_inf, w   the band width and note B's prior df / weight
  vb <- !is.null(res$vb_k) && any(is.finite(res$vb_k))
  gmin <- if (vb) pmin(res$vb_genes_case, res$vb_genes_control) else NULL
  row$frac_short_pool  <- if (vb && any(tested)) mean(gmin[tested] < res$vb_k[1], na.rm = TRUE) else NA_real_
  row$n_scanned_median <- if (vb && any(is.finite(res$vb_scan_case))) median(res$vb_scan_case, na.rm = TRUE)
                          else NA_real_
  row$frac_skipped     <- if (vb) mean(!tested) else NA_real_
  row$tau_var          <- if (vb) res$vb_tau[1] else NA_real_
  pos <- function(x) if (vb && is.finite(x) && x >= 0 || (vb && is.infinite(x))) x else NA_real_
  row$d0               <- pos(if (vb) res$vb_d0[1] else NA_real_)  # Inf = no excess scatter; NA in mode A
  row$d0_inf           <- if (vb && !is.na(row$d0)) as.numeric(is.infinite(row$d0)) else NA_real_
  row$w                <- pos(if (vb) res$vb_w[1] else NA_real_)
  # Scan radius of the case-side pool in log2 units (note D WP4; power_test.R
  # scan_radius). NA for runs written before the column existed.
  rad <- if (vb && !is.null(res$vb_scan_radius)) res$vb_scan_radius[tested] else NULL
  rad <- rad[is.finite(rad)]
  row$scan_radius_median <- if (length(rad)) median(rad) else NA_real_
  row$frac_radius_gt2    <- if (length(rad)) mean(rad > 2) else NA_real_
  as.data.frame(row)
}

#' Scan radius by expression decile (note D WP4), tested genes of a variance arm.
radius_strata <- function(res, expr) {
  if (is.null(res$vb_scan_radius) || !any(is.finite(res$vb_scan_radius))) return(NULL)
  dec <- cut(expr, unique(quantile(expr, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE)
  ok <- !is.na(res$p) & is.finite(res$vb_scan_radius) & !is.na(dec)
  if (!any(ok)) return(NULL)
  a <- lapply(split(res$vb_scan_radius[ok], dec[ok]), function(r)
    c(n = length(r), scan_radius_median = median(r), frac_radius_gt2 = mean(r > 2)))
  data.frame(expr_decile = as.integer(names(a)), do.call(rbind, a), row.names = NULL)
}

#' TPR on the achieved-FDR grid (for averaged FDR-TPR curves).
curve_one <- function(res, is_de, grid = PM$fdr_grid) {
  if (!any(is_de)) return(NULL)
  r <- roc_steps(-res$p, is_de)
  data.frame(fdr = grid, tpr = vapply(grid, function(t) {
    k <- ach_step(r, t); if (k > 0) r$tp[k] / r$P else 0 }, numeric(1)))
}

#' PR curve for DISPLAY: interpolated-precision envelope (max precision over steps
#' with recall >= r) on a recall grid, ranking by score_pe. Defined up to recall 1
#' because untested genes form a final tied block. The scalar AUPRC never uses this.
curve_pr_one <- function(res, is_de, grid = PM$pr_grid) {
  if (!any(is_de)) return(NULL)
  r <- roc_steps(score_pe(res$p, res$tie), is_de)
  rec <- r$tp[-1] / r$P; prec <- r$tp[-1] / (r$tp[-1] + r$fp[-1])
  env <- rev(cummax(rev(prec)))                        # max precision at recall >= rec[i]
  data.frame(recall = grid,
             precision = vapply(grid, function(g) env[which(rec >= g - 1e-12)[1]], numeric(1)))
}

#' TPR by |true log2 FC| bin and by expression decile (true DE genes only).
strata_one <- function(res, is_de, true_lfc, expr, lfc_bins = PM$lfc_bins) {
  called_nom <- !is.na(res$padj) & res$padj < 0.05
  r <- roc_steps(-res$p, is_de); k <- ach_step(r, 0.05)
  s <- -res$p; s[is.na(s)] <- -Inf
  called_ach <- if (k > 0) s >= r$thr[k] else rep(FALSE, length(s))
  dec <- cut(expr, unique(quantile(expr, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE)
  bins <- cut(abs(true_lfc), lfc_bins, include.lowest = TRUE)
  mk <- function(type, f) {
    f <- f[is_de]; if (!length(f) || all(is.na(f))) return(NULL)
    agg <- lapply(split(seq_along(f), f), function(ix) c(
      n = length(ix),
      tpr_nom05 = mean(called_nom[is_de][ix]),
      tpr_ach05 = mean(called_ach[is_de][ix])))
    data.frame(stratum_type = type, stratum = names(agg),
               do.call(rbind, agg), row.names = NULL)
  }
  rbind(mk("abs_lfc_bin", as.character(bins)), mk("expr_decile", dec))
}

# ------------------------------------------------------------------ aggregation

#' Mean and SE over rounds per group. Non-finite values are skipped (so every
#' metric that can be undefined must be defined as NA, not 0, and vice versa:
#' fdr_nom_* is zero-filled on purpose, see metrics_one). NA keys group as "NA".
mean_se <- function(df, keys, cols) {
  g <- do.call(paste, c(lapply(df[keys], as.character), sep = "\r"))
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

#' Paired per-round difference of every method against `ref`.
paired_vs_ref <- function(M, ref, cols = c("tpr_ach05", "pauc05_pe", "auprc_pe", "auprc_p",
                                           "tpr_nom_05", "fdr_nom_05")) {
  R <- M[M$method == ref, c(CELL, "round", cols)]
  if (!nrow(R)) return(NULL)
  X <- merge(M[M$method != ref, c(CELL, "round", "method", cols)], R,
             by = c(CELL, "round"), suffixes = c("", ".ref"))
  if (!nrow(X)) return(NULL)
  for (c in cols) X[[paste0("d_", c)]] <- X[[c]] - X[[paste0(c, ".ref")]]
  mean_se(X, c(CELL, "method"), paste0("d_", cols))
}

#' Paired per-round oracle-minus-production difference for each "<arm>-oracle".
#' Calibration columns first: contamination acts on calibration (P3'), and
#' tpr_ach05 is the expected-null control (invariant to monotone p inflation).
oracle_gap <- function(M, cols = c("fdr_nom_05", "fpr_null05", "n_called_05", "tpr_nom_05",
                                   "tpr_ach05", "pauc05_pe", "auprc_pe", "auprc_p")) {
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

#' Part C breakdown point: smallest target Delta whose mean zero-filled FDR at
#' nominal 0.05 exceeds `limit`, with a percentile bootstrap CI over rounds
#' (rounds resampled within each Delta). Inf = no breakdown inside the grid.
#' d: per-round rows of ONE method in ONE cell family, columns delta_target, fdr_nom_05.
breakdown_point <- function(d, limit = 0.10, B = 2000L, level = 0.95) {
  by_d <- split(d$fdr_nom_05, d$delta_target)
  deltas <- as.numeric(names(by_d))
  bp <- function(means) { k <- which(means > limit); if (length(k)) deltas[min(k)] else Inf }
  est <- bp(vapply(by_d, mean, numeric(1)))
  boot <- replicate(B, bp(vapply(by_d, function(x) mean(x[sample.int(length(x), replace = TRUE)]),
                                 numeric(1))))
  q <- quantile(boot, c((1 - level) / 2, (1 + level) / 2), type = 1, names = FALSE)
  data.frame(breakdown_delta = est, ci_lo = q[1], ci_hi = q[2],
             frac_boot_no_breakdown = mean(!is.finite(boot)))
}
