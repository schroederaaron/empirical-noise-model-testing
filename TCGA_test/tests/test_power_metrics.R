# Unit tests for common/power_metrics.R (plan section 7).
#   Rscript -e 'testthat::test_dir("TCGA_test/tests")'     (from experiments/Noise_Model_Test)
# test_dir() runs with the working directory set to this folder.
library(testthat)
source("../../common/power_metrics.R")

mk_res <- function(p, tie = rep(0, length(p)), padj = p.adjust(p, "BH"))
  data.frame(p = p, padj = padj, lfc = rep(1, length(p)), tie = tie, floor_th = NA_real_)

# O(n^2) reference: average precision over distinct thresholds, straight from the definition.
ap_brute <- function(score, y) {
  thr <- sort(unique(score), decreasing = TRUE); P <- sum(y); prev_tp <- 0; s <- 0
  for (t in thr) {
    called <- score >= t; tp <- sum(called & y); prec <- tp / sum(called)
    s <- s + (tp - prev_tp) / P * prec; prev_tp <- tp
  }
  s
}

test_that("AUPRC: perfect ranking = 1, all tied = prevalence exactly", {
  y <- c(rep(TRUE, 7), rep(FALSE, 93))
  expect_equal(auprc(roc_steps(rev(seq_along(y)), y)), 1)
  expect_identical(auprc(roc_steps(rep(1, 100), y)), 7 / 100)
  m <- metrics_one(mk_res(rep(0.5, 100)), y, ifelse(y, 1, 0), 7)
  expect_identical(m$auprc_p, m$prevalence)
  expect_equal(m$auprc_excess_pe + m$prevalence, m$auprc_pe)
})

test_that("AUPRC matches scikit-learn average_precision_score (Python fixture)", {
  d <- read.csv("fixtures/auprc_sklearn.csv"); e <- read.csv("fixtures/auprc_sklearn_expected.csv")
  for (k in e$case) {
    x <- d[d$case == k, ]
    expect_equal(auprc(roc_steps(x$score, x$y == 1)), e$ap[e$case == k], tolerance = 1e-12, label = k)
  }
})

test_that("AUPRC matches an O(n^2) brute force, with and without ties", {
  set.seed(1)
  for (i in 1:20) {
    n <- 60; y <- runif(n) < 0.3; if (!any(y)) y[1] <- TRUE
    s <- if (i %% 2) rnorm(n) + y else round(rnorm(n) + y)
    expect_equal(auprc(roc_steps(s, y)), ap_brute(s, y), tolerance = 1e-12)
  }
})

test_that("auprc_* equals the old ap_* estimator; untested genes rank last", {
  old_ap <- function(r) { prec <- r$tp / pmax(r$tp + r$fp, 1); sum(diff(r$tp) * prec[-1]) / r$P }
  set.seed(2); n <- 300; y <- runif(n) < 0.1
  p <- ifelse(y, rbeta(n, 0.3, 3), runif(n)); p[sample(n, 20)] <- NA
  tie <- runif(n)
  m <- metrics_one(mk_res(p, tie), y, ifelse(y, 1, 0), sum(y) + 5)
  expect_equal(m$auprc_p, old_ap(roc_steps(-p, y)))
  expect_equal(m$auprc_pe, old_ap(roc_steps(score_pe(p, tie), y)))
  expect_equal(m$auprc_all_pe, m$auprc_pe * sum(y) / (sum(y) + 5))
})

test_that("random ranking: AUPRC ~ prevalence", {
  set.seed(3); y <- runif(5000) < 0.05
  ex <- replicate(50, metrics_one(mk_res(runif(5000)), y, ifelse(y, 1, 0), sum(y))$auprc_excess_pe)
  expect_lt(abs(mean(ex)), 3 * sd(ex) / sqrt(length(ex)) + 1e-3)
})

test_that("I1: a round with zero calls contributes FDR 0, fdr_cond is NA", {
  y <- c(rep(TRUE, 5), rep(FALSE, 95))
  m0 <- metrics_one(mk_res(rep(0.9, 100)), y, ifelse(y, 1, 0), 5)
  expect_identical(m0$fdr_nom_05, 0); expect_true(is.na(m0$fdr_cond_05)); expect_identical(m0$zero_disc_05, 1L)
  # one false call only -> FDP 1 both ways
  p <- rep(0.9, 100); p[10] <- 1e-9
  m1 <- metrics_one(mk_res(p), y, ifelse(y, 1, 0), 5)
  expect_identical(m1$fdr_nom_05, 1); expect_identical(m1$fdr_cond_05, 1)
  # mean over the two rounds: zero-filled 0.5, the old conditional mean was 1
  M <- rbind(cbind(round = 1, m0), cbind(round = 2, m1)); M$k <- "a"
  s <- mean_se(M, "k", c("fdr_nom_05", "fdr_cond_05"))
  expect_equal(s$fdr_nom_05_mean, 0.5); expect_equal(s$fdr_cond_05_mean, 1)
})

test_that("prec_top_nde is R-precision", {
  y <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, rep(FALSE, 10))
  p <- c(0.001, 0.002, 0.003, 0.004, 0.5, 0.6, rep(0.9, 10))
  expect_equal(metrics_one(mk_res(p), y, ifelse(y, 1, 0), 3)$prec_top_nde, 2 / 3)
})

test_that("PR curve: reaches recall 1 at precision >= prevalence; envelope is non-increasing", {
  set.seed(4); n <- 400; y <- runif(n) < 0.1
  p <- ifelse(y, rbeta(n, 0.5, 2), runif(n)); p[sample(n, 30)] <- NA
  pr <- curve_pr_one(mk_res(p, runif(n)), y)
  expect_identical(pr$recall[nrow(pr)], 1)
  expect_false(anyNA(pr$precision))
  expect_gte(pr$precision[nrow(pr)], mean(y) - 1e-12)
  expect_true(all(diff(pr$precision) <= 1e-12))
  expect_null(curve_pr_one(mk_res(p), rep(FALSE, n)))
})

test_that("mean_se groups NA keys together instead of dropping them", {
  d <- data.frame(a = c(NA, NA, "x"), v = c(1, 3, 5))
  s <- mean_se(d, "a", "v")
  expect_equal(nrow(s), 2L); expect_equal(s$v_mean[is.na(s$a)], 2)
})

test_that("breakdown_point finds the first Delta above the limit", {
  set.seed(5)
  d <- data.frame(delta_target = rep(c(0, 0.25, 0.5, 1), each = 20),
                  fdr_nom_05 = c(rnorm(20, 0.04, 0.01), rnorm(20, 0.06, 0.01),
                                 rnorm(20, 0.30, 0.01), rnorm(20, 0.50, 0.01)))
  b <- breakdown_point(d, 0.10, B = 200)
  expect_identical(b$breakdown_delta, 0.5); expect_identical(b$ci_lo, 0.5); expect_identical(b$ci_hi, 0.5)
  d$fdr_nom_05 <- 0.01
  expect_identical(breakdown_point(d, 0.10, B = 50)$breakdown_delta, Inf)
})
