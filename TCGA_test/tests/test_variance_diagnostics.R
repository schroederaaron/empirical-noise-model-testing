# Unit tests for common/variance_diagnostics.R (implementation note G section 4).
library(testthat)
source("../../common/variance_diagnostics.R")

sim_null <- function(G = 3000, n = 5, imb_sd = 0, seed = 1) {
  set.seed(seed)
  mu <- rlnorm(G, 4, 1.5); sd0 <- 0.3 * exp(0.5 * rnorm(G))     # log2-scale noise, trend + scatter
  sdA <- sd0 * exp(imb_sd * rnorm(G))                           # gene-specific imbalance in group A
  y <- cbind(matrix(log2(mu + 1) + sdA * rnorm(G * n), G), matrix(log2(mu + 1) + sd0 * rnorm(G * n), G))
  x <- 2^y - 1; rownames(x) <- paste0("g", 1:G)
  list(x = x, group = rep(c("A", "B"), each = n), sdA = sdA, sdB = sd0)
}
#' z-test with the TRUE SDs: its p-value ignores the estimated variances, so any
#' correlation with L, |L| or v would be manufactured by the diagnostic itself.
ztest_p <- function(s) {
  y <- log2(s$x + 1); a <- s$group == "A"
  z <- (rowMeans(y[, !a]) - rowMeans(y[, a])) / sqrt(s$sdA^2 / sum(a) + s$sdB^2 / sum(!a))
  2 * pnorm(-abs(z))
}

test_that("vdiag_frame: one row per gene and input, variances match var()", {
  s <- sim_null(G = 200)
  vd <- vdiag_frame(s$x, s$group, "TPM")
  expect_identical(nrow(vd), 200L); expect_identical(anyDuplicated(vd[c("gene", "input")]), 0L)
  y <- log2(s$x[7, ] + 1)
  expect_equal(vd$s2_A[7], var(y[s$group == "A"])); expect_equal(vd$s2_B[7], var(y[s$group == "B"]))
  expect_equal(vd$s2[7], (var(y[1:5]) + var(y[6:10])) / 2)
})

test_that("LOESS trends and the three quantities are finite", {
  s <- sim_null()
  vq <- vdiag_quantities(vdiag_frame(s$x, s$group))
  for (c in c("T_A", "T_B", "T", "L", "absL", "v")) expect_true(all(is.finite(vq[[c]])), label = c)
  expect_true(all(vq$T_A >= VD$floor))
})

test_that("equal variances: |rho| < 0.03 for all three quantities (variance-blind test)", {
  r <- sapply(1:4, function(sd) {
    s <- sim_null(seed = sd); vq <- vdiag_quantities(vdiag_frame(s$x, s$group))
    p <- ztest_p(s)
    unlist(vdiag_round(vq, data.frame(p = p, padj = p.adjust(p, "BH")), rep(FALSE, 3000))$corr[1:3])
  })
  expect_true(all(abs(rowMeans(r)) < 0.03))
})

test_that("vdiag_round: bins, enrichment rows, flags on DE rounds, enrichment CI", {
  s <- sim_null(); vq <- vdiag_quantities(vdiag_frame(s$x, s$group))
  p <- runif(3000); is_de <- seq_len(3000) <= 300; p[is_de] <- p[is_de] / 1000
  g <- vdiag_round(vq, data.frame(p = p, padj = p.adjust(p, "BH")), is_de, thr = 2)
  expect_setequal(unique(g$bins$quantity), c("absL", "L", "v"))
  expect_true(all(g$bins$decile %in% 1:10))
  expect_setequal(unique(g$enrich$split), c("absL_top5", "L_top5", "v_top5", "absL_nullcal"))
  expect_setequal(g$flag$flag, c("empirical_top5", "null_calibrated"))
  e <- enrichment_ci(g$enrich[g$enrich$split == "v_top5", ][rep(2, 10), ], B = 200)
  expect_true(e$ci_lo <= e$enrichment && e$enrichment <= e$ci_hi)
})

#| Sign convention: L > 0 = group A (control in power_test.R) noisier than its trend
#| relative to B. Noise injected into A puts noisy genes at large L and (for the
#| current model, whose control pool comes from quieter neighbours) small p, so
#| rho_signed < 0. Injected into B, the sign flips (+0.08 measured 28/10/2026).
test_that("current TOX model: injected gene-specific imbalance gives a negative rho_signed", {
  skip_if_not(Sys.getenv("TOX_TEST_FORTRAN") == "1", "set TOX_TEST_FORTRAN=1 (needs build/libtensor-omics.so)")
  root <- normalizePath("../../../..")
  old <- setwd(root); on.exit(setwd(old))
  .libPaths(c(file.path(root, "external/docker_r_libs"), .libPaths()))
  suppressMessages(source("rcpp/tensoromics_functions.R"))
  r <- sapply(1:3, function(sd) {
    s <- sim_null(n = 5, imb_sd = 0.35, seed = sd)
    a <- s$group == "A"
    res <- tox_compute_noise_pvalues_pipeline_exact(t(s$x[, !a]), t(s$x[, a]), norm_method = 1L,
             k_start = 20L, k_step = 1L, k_max = 50L, tau = 0.1, max_pool_size = 70000L)
    p <- res$pvalues_own; p[p < 0] <- NA
    vq <- vdiag_quantities(vdiag_frame(s$x, s$group))
    vdiag_round(vq, data.frame(p = p, padj = p.adjust(p, "BH")), rep(FALSE, 3000))$corr$rho_signed
  })
  expect_lt(mean(r), -0.02)
})
