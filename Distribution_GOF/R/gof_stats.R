# gof_stats.R -- CvM, AD, KS on PIT values (brief 7.2).
#
# The three formulas are the brief's, verbatim; they agree with goftest::cvm.test,
# goftest::ad.test and stats::ks.test (unit-tested). Textbook p-values of these
# statistics are NOT computed anywhere in the pipeline (brief rule 4): the only
# calibrated reference is the parametric bootstrap.

cvm_u <- function(u) { u <- sort(u); N <- length(u)
  1/(12*N) + sum((u - (2*seq_len(N) - 1)/(2*N))^2) }
ad_u  <- function(u, eps = 1e-12) { u <- pmin(pmax(sort(u), eps), 1 - eps)
  N <- length(u); i <- seq_len(N)
  -N - mean((2*i - 1) * (log(u) + log1p(-rev(u)))) }
ks_u  <- function(u) { u <- sort(u); N <- length(u)
  max(seq_len(N)/N - u, u - (seq_len(N) - 1)/N) }

#' W2, A2, D, their N-normalised forms, and how many u the AD clamp touched.
gof_all <- function(u, eps = 1e-12) {
  u <- as.numeric(u); u <- u[is.finite(u)]
  N <- length(u)
  if (N < 2L) return(list(N = N, W2 = NA_real_, A2 = NA_real_, D = NA_real_,
                          W2_N = NA_real_, A2_N = NA_real_, n_clamped = NA_integer_))
  W2 <- cvm_u(u); A2 <- ad_u(u, eps)
  list(N = N, W2 = W2, A2 = A2, D = ks_u(u), W2_N = W2 / N, A2_N = A2 / N,
       n_clamped = sum(u < eps | u > 1 - eps))
}

#' Redraw u R times from the same (a, b) and report every draw's statistics.
#' The headline uses draw 1 (the draw made with the model's PIT seed).
rerandomise <- function(a, b, R, seed) {
  data.table::rbindlist(lapply(seq_len(R), function(r) {
    s <- gof_all(randomise_pit(a, b, if (r == 1L) seed else seed_for(seed, "rerand", r)))
    data.table::data.table(draw = r, W2 = s$W2, A2 = s$A2, D = s$D, N = s$N)
  }))
}

#' 50-bin PIT histogram counts (for the calibrated bootstrap band).
pit_hist <- function(u, bins = 50L)
  tabulate(pmin(bins, pmax(1L, ceiling(u[is.finite(u)] * bins))), nbins = bins)

#' Per-sample mean of qnorm(clamp(u)) (brief 8.3, sample-level diagnostic).
sample_qnorm_mean <- function(U, eps = 1e-12) colMeans(stats::qnorm(pmin(pmax(U, eps), 1 - eps)), na.rm = TRUE)
