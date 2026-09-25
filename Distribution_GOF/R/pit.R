# pit.R -- randomised PIT from stored bounds (brief 7.1).
#
# For every observation a = F(y - 1) and b = F(y) with F(-1) = 0. Both are stored,
# not only u. u = a + (b - a) * U(0,1) under a recorded seed.
#
# G6 asserts: 0 <= a <= b <= 1, and b - a equal to exp(logpmf(y)).
# The second assert CANNOT hold as a pure 1e-10 RELATIVE criterion in double
# precision: once b is close to 1 (the upper tail), b - a is a difference of two
# numbers near 1 and carries an absolute rounding error of order eps * b whatever
# the family (it already fails for R's own pnbinom/dnbinom). It is therefore
# implemented as
#     |(b - a) - p| <= 1e-10 * p + 64 * eps * b + 64 * double.xmin,
# i.e. 1e-10 relative wherever b - a is representable to that precision, and the
# representational floor elsewhere (the double.xmin term covers subnormal
# probabilities, ~1e-308 and below, which carry no relative precision at all). The number of observations where the floor
# was the binding term is reported with every G6 result.

G6_REL <- 1e-10
G6_ABS <- 64 * .Machine$double.eps
G6_TINY <- 64 * .Machine$double.xmin     # below the smallest normal double no relative precision exists

#' PIT bounds for `genes` of a fit. Returns a, b, logp (genes x samples) and G6 facts.
#' Stops (hard gate) on any violation.
pit_bounds <- function(Y, fitres, fam, genes, label = "") {
  n <- ncol(Y)
  A <- B <- P <- matrix(NA_real_, length(genes), n, dimnames = list(genes, colnames(Y)))
  for (g in genes) {
    par <- fit_par(fitres$par[[g]])
    y <- as.numeric(Y[g, ])
    A[g, ] <- fam$cdf(y - 1, par)
    B[g, ] <- fam$cdf(y, par)
    P[g, ] <- fam$logpmf(y, par)
  }
  p <- exp(P)
  bad_order <- !(A >= 0 & A <= B & B <= 1)
  dev <- abs((B - A) - p)
  tol <- G6_REL * p + G6_ABS * B + G6_TINY
  bad_inc <- dev > tol
  if (any(bad_order, na.rm = TRUE) || any(is.na(bad_order)) || any(bad_inc, na.rm = TRUE)) {
    i <- which(bad_order | is.na(bad_order) | bad_inc, arr.ind = TRUE)[1, ]
    g <- genes[i[1]]; j <- i[2]
    stop(sprintf(paste0("G6 FAILED [%s %s]: gene %s obs %d: y=%g a=%.17g b=%.17g p=%.17g ",
                        "|(b-a)-p|=%.3g tol=%.3g (%d order / %d increment violations)"),
                 fam$name, label, g, j, Y[g, j], A[i[1], j], B[i[1], j], p[i[1], j],
                 dev[i[1], j], tol[i[1], j], sum(bad_order | is.na(bad_order)), sum(bad_inc, na.rm = TRUE)))
  }
  rel_ok <- dev <= G6_REL * p
  list(a = A, b = B, logp = P,
       g6 = list(model = fam$name, n_obs = length(A), passed = TRUE,
                 n_floor_binding = sum(!rel_ok),
                 max_rel_dev_where_relative = if (any(rel_ok & p > 0)) max((dev / p)[rel_ok & p > 0]) else NA_real_))
}

randomise_pit <- function(a, b, seed) {
  u <- with_seed(seed, a + (b - a) * stats::runif(length(a)))
  dim(u) <- dim(a); dimnames(u) <- dimnames(a)
  u
}
