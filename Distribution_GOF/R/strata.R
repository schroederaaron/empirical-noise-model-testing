# strata.R -- stratification rules (brief 8.3).
#
# Each gene-level rule is a function (Y, O, X) -> factor over genes. It is applied
# to the observed data AND re-applied to every bootstrap dataset: stratifying by
# empirical dispersion or zero fraction is selection on the same data that
# produce the PIT, and re-applying the rule in each replicate is what calibrates
# that selection. `condition` is per OBSERVATION (design level), not per gene.

STRATUM_RULES <- list(
  mean = function(Y, O, X) rank_bins(log2_mean_cpm(Y, O), 3L, "mean_T"),
  disp = function(Y, O, X) {
    y <- edgeR::DGEList(counts = Y); y$offset <- O
    d <- edgeR::estimateDisp(y, design = X, tagwise = TRUE)$tagwise.dispersion
    rank_bins(d, 3L, "disp_T")
  },
  zero = function(Y, O, X) {
    z <- rowMeans(Y == 0)
    factor(ifelse(z == 0, "zero_0", ifelse(z <= 0.1, "zero_(0,0.1]", "zero_>0.1")),
           levels = c("zero_0", "zero_(0,0.1]", "zero_>0.1"))
  }
)

#' Gene strata for every rule on (Y, O, X), restricted to `genes`.
apply_strata <- function(Y, O, X, genes) {
  Yg <- Y[genes, , drop = FALSE]; Og <- O[genes, , drop = FALSE]
  lapply(STRATUM_RULES, function(f) stats::setNames(f(Yg, Og, X), genes))
}

#' Statistics per stratum. U: genes x samples PIT; gs: output of apply_strata;
#' obs_group: factor over samples (the `condition` rule), NULL for ~ 1 designs.
strata_stats <- function(U, gs, obs_group = NULL) {
  out <- list()
  for (rule in names(gs)) {
    f <- gs[[rule]][rownames(U)]
    for (lev in levels(f)) {
      s <- gof_all(U[which(f == lev), , drop = FALSE])
      out[[length(out) + 1L]] <- data.table::data.table(rule = rule, stratum = lev, N = s$N,
                                                         W2 = s$W2, A2 = s$A2, D = s$D)
    }
  }
  og <- if (is.null(obs_group)) factor(rep("all", ncol(U))) else obs_group
  for (lev in levels(og)) {
    s <- gof_all(U[, og == lev, drop = FALSE])
    out[[length(out) + 1L]] <- data.table::data.table(rule = "condition", stratum = lev, N = s$N,
                                                       W2 = s$W2, A2 = s$A2, D = s$D)
  }
  data.table::rbindlist(out)
}
