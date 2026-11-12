# variance_diagnostics.R
# -----------------------------------------------------------------------------
# Implementation note G: does gene-wise noise (overall, or its imbalance between the
# two groups) drive TOX's false positives? Base R only (stats::loess, cor).
#
#   vdiag_frame()      per gene, method-independent: group variances and means on
#                      y = log2(x + 1) of the input the TOX arms use (power_test.R
#                      writes it as vdiag_<job>.csv.gz)
#   vdiag_quantities() the trends and the three quantities of note G section 1:
#                        L  = log2(s2_A / T_A(m_A)) - log2(s2_B / T_B(m_B))  signed imbalance
#                        |L|                                                 absolute imbalance
#                        v  = log2(s2 / T(m))                                overall relative noise
#   vdiag_round()      per round and method: A1 rank correlations, A2 FPR by decile
#                      and top-5 % vs rest, A3 flag utility
#   gene_classes()     per round and method: TP / FP / FN / TN per gene with v, |v|
#                      and the expression spread; keeps FP, FN and a random sample
#
# |v| = |log2(s2 / T(m))| is the distance of a gene's variance from the LOESS trend
# (0 = exactly on the trend). Its rank correlation with p (true negatives: do noisy
# or quiet genes give false positives; DE genes: are on-trend genes found more
# easily) and the call precision / TPR by |v| decile answer "are genes whose
# variance is close to the smoothed variance more often true positives?".
#
# Group A = control, B = case (as labelled by the sampler). Each group is normalised
# by its own LOESS trend at its own mean, so variance differences explained by a
# mean difference (DE genes) or by a systematic group difference are removed.
# -----------------------------------------------------------------------------

VD <- list(span = 0.3, degree = 1L, floor = 1e-8, alphas = c(0.001, 0.01, 0.05), top = 0.05)

#' Method-independent per-gene variance file of one dataset (note G section 2).
#' x: genes x samples (linear input), group: "A" / "B" per column.
vdiag_frame <- function(x, group, input = "TPM") {
  y <- log2(pmax(x, 0) + 1)
  a <- group == "A"; b <- group == "B"
  rv <- function(z) { mu <- rowMeans(z); rowSums((z - mu)^2) / (ncol(z) - 1) }
  s2A <- rv(y[, a, drop = FALSE]); s2B <- rv(y[, b, drop = FALSE])
  data.frame(gene = rownames(x), input = input, s2_A = s2A, s2_B = s2B,
             s2 = ((sum(a) - 1) * s2A + (sum(b) - 1) * s2B) / (ncol(x) - 2),
             m_A = log2(rowMeans(x[, a, drop = FALSE]) + 1), m_B = log2(rowMeans(x[, b, drop = FALSE]) + 1),
             m = log2(rowMeans(x) + 1), row.names = NULL)
}

#' LOESS trend of a variance on the mean axis, on the variance scale (non-robust,
#' degree 1, span 0.3), floored.
vd_trend <- function(s2, m) {
  ok <- is.finite(s2) & is.finite(m)
  fit <- stats::loess(s2 ~ m, data = data.frame(s2 = s2[ok], m = m[ok]), span = VD$span,
                      degree = VD$degree, family = "gaussian")
  t <- rep(NA_real_, length(s2)); t[ok] <- stats::fitted(fit)
  pmax(t, VD$floor)
}

#' The three quantities of note G section 1.
vdiag_quantities <- function(vd) {
  fl <- function(s) pmax(s, VD$floor)
  vd$T_A <- vd_trend(vd$s2_A, vd$m_A); vd$T_B <- vd_trend(vd$s2_B, vd$m_B); vd$T <- vd_trend(vd$s2, vd$m)
  vd$L <- log2(fl(vd$s2_A) / vd$T_A) - log2(fl(vd$s2_B) / vd$T_B)
  vd$absL <- abs(vd$L)
  vd$v <- log2(fl(vd$s2) / vd$T)
  vd$absv <- abs(vd$v)
  vd
}

#' Note G A1-A3 for one round. vq = vdiag_quantities() of the round; res = one
#' method's per-gene results (p, padj), aligned to vq$gene; is_de aligned too.
#' thr = null-calibrated |L| threshold (95th percentile in the pi1 = 0 rounds of
#' the same cohort and n), or NA.
vdiag_round <- function(vq, res, is_de, thr = NA_real_) {
  tn <- !is_de & !is.na(res$p) & is.finite(vq$L)
  sp <- function(q) if (sum(tn) > 10) suppressWarnings(stats::cor(q[tn], res$p[tn], method = "spearman")) else NA_real_
  de <- is_de & !is.na(res$p) & is.finite(vq$absv)
  corr <- data.frame(rho_signed = sp(vq$L), rho_abs = sp(vq$absL), rho_overall = sp(vq$v), rho_absv = sp(vq$absv),
                     rho_absv_de = if (sum(de) > 10) suppressWarnings(stats::cor(vq$absv[de], res$p[de], method = "spearman"))
                                   else NA_real_, n_tn = sum(tn))

  called <- !is.na(res$padj) & res$padj < 0.05
  bins <- do.call(rbind, lapply(c("absL", "L", "v", "absv"), function(qn) {
    q <- vq[[qn]]; ok <- !is.na(res$p) & is.finite(q)
    if (sum(ok) < 20) return(NULL)
    br <- unique(stats::quantile(q[ok], seq(0, 1, 0.1)))
    dec <- cut(q, br, include.lowest = TRUE, labels = FALSE)
    do.call(rbind, lapply(sort(unique(dec[ok])), function(d) {
      i <- which(ok & dec == d); t <- i[!is_de[i]]
      nc <- sum(called[i])
      data.frame(quantity = qn, decile = d, n_tn = length(t),
                 fpr_001 = mean(res$p[t] <= 0.001), fpr_01 = mean(res$p[t] <= 0.01), fpr_05 = mean(res$p[t] <= 0.05),
                 n_called_05 = nc, false_share_05 = if (nc > 0) sum(called[i] & !is_de[i]) / nc else NA_real_,
                 tpr_05 = if (any(is_de[i])) mean(called[i][is_de[i]]) else NA_real_)
    }))
  }))

  # top 5 % vs rest (empirical), and above vs below the null-calibrated |L| threshold
  split_counts <- function(top, label) {
    if (!any(tn & top) || !any(tn & !top)) return(NULL)
    do.call(rbind, lapply(VD$alphas, function(a) data.frame(
      split = label, alpha = a,
      fp_top = sum(tn & top & res$p <= a), n_top = sum(tn & top),
      fp_rest = sum(tn & !top & res$p <= a), n_rest = sum(tn & !top))))
  }
  ok <- !is.na(res$p)
  qtop <- function(q) { cut <- stats::quantile(q[ok & is.finite(q)], 1 - VD$top); is.finite(q) & q > cut }
  enrich <- rbind(split_counts(qtop(vq$absL), "absL_top5"), split_counts(qtop(vq$L), "L_top5"),
                  split_counts(qtop(vq$v), "v_top5"),
                  if (is.finite(thr)) split_counts(is.finite(vq$absL) & vq$absL > thr, "absL_nullcal"))

  # A3: flag utility on rounds with DE genes: BH on the unflagged genes only
  flag <- NULL
  if (any(is_de)) {
    fl_rows <- function(flagged, label) {
      keep <- ok & !flagged
      padj <- rep(NA_real_, length(res$p)); padj[keep] <- stats::p.adjust(res$p[keep], "BH")
      cr <- !is.na(padj) & padj < 0.05; ck <- called
      data.frame(flag = label, share_flagged = mean(flagged[ok]), share_de_flagged = mean(flagged[is_de & ok]),
                 fdr_kept = if (any(ck)) sum(ck & !is_de) / sum(ck) else 0,
                 fdr_removed = if (any(cr)) sum(cr & !is_de) / sum(cr) else 0,
                 tpr_kept = mean(ck[is_de]), tpr_removed = mean(cr[is_de]))
    }
    flag <- rbind(fl_rows(qtop(vq$absL), "empirical_top5"),
                  if (is.finite(thr)) fl_rows(is.finite(vq$absL) & vq$absL > thr, "null_calibrated"))
  }
  list(corr = corr, bins = bins, enrich = enrich, flag = flag)
}

#' Per-gene classes of one round and method (BH < 0.05): TP / FP / FN / TN with p,
#' padj, log2fc, expression, v, |v| and the expression spread (NA when the method
#' has none). Keeps every FP and FN and `n_sample` random genes (sample = TRUE, any
#' class) as the "all genes" reference.
gene_classes <- function(vq, res, is_de, expr, n_sample = 200L) {
  called <- !is.na(res$padj) & res$padj < 0.05
  cls <- ifelse(is_de, ifelse(called, "TP", "FN"), ifelse(called, "FP", "TN"))
  smp <- seq_along(cls) %in% sample.int(length(cls), min(n_sample, length(cls)))
  k <- cls %in% c("FP", "FN") | smp
  data.frame(gene = res$gene, class = cls, sample = smp, p = res$p, padj = res$padj, log2fc = res$log2fc,
             expr = expr, v = vq$v, absv = vq$absv,
             spread = if (is.null(res$expr_spread)) NA_real_ else res$expr_spread)[k, ]
}

#' Enrichment of false positives in the flagged set: ratio of FPRs, top / rest,
#' pooled over rounds, with a percentile bootstrap CI over rounds.
enrichment_ci <- function(d, B = 1000L) {
  r <- function(x) (sum(x$fp_top) / sum(x$n_top)) / (sum(x$fp_rest) / sum(x$n_rest))
  est <- r(d)
  bt <- replicate(B, r(d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]))
  q <- stats::quantile(bt[is.finite(bt)], c(0.025, 0.975), names = FALSE)
  data.frame(n_rounds = nrow(d), fpr_top = sum(d$fp_top) / sum(d$n_top), fpr_rest = sum(d$fp_rest) / sum(d$n_rest),
             enrichment = est, ci_lo = q[1], ci_hi = q[2])
}
