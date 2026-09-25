# heldout.R -- K-fold held-out predictive log score (brief 8.2).
#
# Folds are over samples, stratified by the design factor so every level stays in
# training; for n < 16, leave-one-sample-out. The held-out sample's offset comes
# from the full-data offsets (a small leakage; README notes it).

make_folds <- function(n, group, K, seed) {
  if (n < 16L) return(seq_len(n))                        # leave-one-sample-out
  g <- if (is.null(group)) factor(rep("all", n)) else group
  f <- integer(n)
  with_seed(seed, for (lev in levels(g)) {
    idx <- which(g == lev)
    f[idx[sample.int(length(idx))]] <- rep_len(seq_len(K), length(idx))
  })
  for (k in unique(f)) {
    tr <- g[f != k]
    if (length(setdiff(levels(g), unique(as.character(tr)))))
      stop("fold ", k, " removes every sample of a design level; lower K_FOLDS")
  }
  f
}

#' Held-out log p and PIT bounds for every observation of `genes`.
heldout_score <- function(Y, X, O, fam, genes, folds, n_cores = 1L) {
  n <- ncol(Y)
  one <- function(g) {
    y <- as.numeric(Y[g, ]); off <- as.numeric(O[g, ])
    lp <- a <- b <- rep(NA_real_, n); nfail <- 0L
    for (k in sort(unique(folds))) {
      te <- which(folds == k); tr <- which(folds != k)
      fr <- fit_one_gene(fam, y[tr], X[tr, , drop = FALSE], off[tr])
      if (!fr$ok) { nfail <- nfail + 1L; next }
      par <- c(list(mu = fam$mu_fn(fr$beta, fr$extra, X[te, , drop = FALSE], off[te])), fr$extra)
      lp[te] <- fam$logpmf(y[te], par)
      a[te] <- fam$cdf(y[te] - 1, par); b[te] <- fam$cdf(y[te], par)
    }
    list(lp = lp, a = a, b = b, nfail = nfail)
  }
  res <- if (n_cores > 1L) parallel::mclapply(genes, one, mc.cores = n_cores, mc.preschedule = TRUE)
         else lapply(genes, one)
  pick <- function(field) {
    M <- t(vapply(res, function(r) if (is.list(r)) r[[field]] else rep(NA_real_, n), numeric(n)))
    dimnames(M) <- list(genes, colnames(Y)); M
  }
  list(lp = pick("lp"), a = pick("a"), b = pick("b"),
       fold_failures = stats::setNames(vapply(res, function(r) if (is.list(r)) r$nfail else NA_integer_, 0L), genes))
}

#' Score summaries across models on the genes every model scored completely.
#' The interval from resampling genes is labelled approximate: genes are not independent.
heldout_summary <- function(ho, ref = "nb", n_resample = 2000L, seed = 1L) {
  models <- names(ho)
  full <- Reduce(intersect, lapply(ho, function(h) rownames(h$lp)[rowSums(!is.finite(h$lp)) == 0]))
  per_gene <- sapply(ho, function(h) rowMeans(h$lp[full, , drop = FALSE]))
  if (length(full) == 1L) per_gene <- matrix(per_gene, 1, dimnames = list(full, models))
  tab <- data.table::data.table(model = models, n_genes = length(full),
                                total = vapply(models, function(m) sum(ho[[m]]$lp[full, ]), 0),
                                per_obs = colMeans(per_gene))
  if (ref %in% models) {
    d <- per_gene - per_gene[, ref]
    ci <- with_seed(seed, t(replicate(n_resample, colMeans(d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]))))
    tab[, `:=`(delta_vs_nb_mean = colMeans(d), delta_vs_nb_median = apply(d, 2, stats::median),
               frac_genes_better_than_nb = colMeans(d > 0),
               delta_vs_nb_approx_lo = apply(ci, 2, stats::quantile, 0.025),
               delta_vs_nb_approx_hi = apply(ci, 2, stats::quantile, 0.975))]
  }
  tab
}
