# Unit tests for the power_test.R simulator (plan section 7): per-gene dispersion
# recycling in draw_counts, the balanced effect draw (I3), Delta targeting (I4),
# and the thinning balance of Part R. Only the named top-level definitions are
# evaluated, so no package, Fortran build or data is needed.
library(testthat)

source_defs <- function(file, names, env) {
  for (e in parse(file))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% names)
      eval(e, env)
}
E <- environment()     # the definitions must see this file's PCFG
source_defs("../../Simulated_data/scripts/calibration_test.R",
            c("draw_counts", "TPOIS_DF", "TPOIS_Z_CLAMP"), E)
source_defs("../scripts/power_test.R",
            c("draw_effects", "balance_effects", "effect_diag", "make_truth_power",
              "disp_hetero", "simulate_power", "simulate_real", "%||%"), E)
PCFG <- list(frac_up = 0.5, lfc_range = c(0.25, 4), lib_size = 4e7, disp = 0.2,
             disp_hetero = list(phi_inf = 0.1, c = 1, sigma_hat = 0.8),
             R = list(lfc_tol = 0.5))

test_that("draw_counts recycles a per-gene dispersion down the columns (nb, lnpois)", {
  set.seed(1)
  G <- 4; S <- 20000; mu <- 200
  phi <- c(0.05, 0.2, 0.8, 1.6)
  lam <- matrix(mu, G, S)
  for (dist in c("nb", "lnpois")) {
    x <- draw_counts(lam, dist, phi)
    # NB / Poisson-lognormal: var = mu + phi mu^2  ->  phi_hat = (var - mean) / mean^2
    phi_hat <- (apply(x, 1, var) - rowMeans(x)) / rowMeans(x)^2
    expect_equal(phi_hat, phi, tolerance = 0.08, label = dist)
  }
})

test_that("I3 balanced draw: nulls exactly null, no sign flips, Delta = 0", {
  set.seed(2)
  flips <- below <- med <- up <- numeric(0)
  for (i in 1:1000) {
    tr <- make_truth_power(500, 0.1, 3)
    de <- tr$is_de
    expect_true(all(tr$true_lfc[!de] == 0))
    expect_equal(sum(tr$pi_a[de] * 2^tr$true_lfc[de]), sum(tr$pi_a[de]), tolerance = 1e-9)
    flips <- c(flips, tr$frac_sign_flip); below <- c(below, tr$frac_below_floor)
    med <- c(med, tr$eff_lfc_median_abs); up <- c(up, mean(tr$true_lfc[de] > 0))
  }
  expect_identical(max(flips), 0)
  expect_equal(mean(up), 0.5, tolerance = 0.01)
  # target: median of log-uniform(0.25, 4) = 1. Balancing shrinks one side only, and
  # it is (almost) always the UP side: an up gene gains pi (2^e - 1), unbounded, a down
  # gene loses at most pi. Measured at G = 5000: median |eff| 0.64 (up 0.41, down 1.01),
  # 16.5% below the 0.25 floor. Known trade-off (plan I3), reported per job in the manifest.
  expect_gt(median(med), 0.6); expect_lt(median(med), 1.05)
  expect_lt(median(below), 0.2)
})

test_that("het < 1: nulls still exactly null in every case sample", {
  set.seed(3)
  tr <- make_truth_power(400, 0.3, 5, het = 0.5)
  expect_true(all(tr$pi_B[!tr$is_de, ] == tr$pi_a[!tr$is_de]))
})

test_that("I4: absolute mode hits the target Delta", {
  set.seed(4)
  for (dt in c(0, 0.25, 0.5, 1)) {
    tr <- make_truth_power(2000, 0.1, 3, mode = "absolute", delta_target = dt)
    expect_equal(tr$delta, dt, tolerance = 1e-9)
    expect_true(all(tr$true_lfc[!tr$is_de] == 0))
  }
  tr <- make_truth_power(2000, 0.1, 3, mode = "absolute")           # Part N: raw draw kept
  expect_false(isTRUE(all.equal(tr$delta, 0)))
})

test_that("balance_effects: thinning losses equal, signs kept, one-sided input flagged", {
  set.seed(5)
  w <- rlnorm(300); w <- w / sum(w); q <- numeric(300); i <- sample(300, 60)
  q[i] <- sample(c(-1, 1), 60, TRUE) * exp(runif(60, log(0.25), log(4)))
  b <- balance_effects(q, w, thinning = TRUE)
  expect_true(b$ok); expect_identical(sign(b$q), sign(q))
  up <- b$q > 0; dn <- b$q < 0
  expect_equal(sum(w[up] * (1 - 2^-b$q[up])), sum(w[dn] * (1 - 2^b$q[dn])), tolerance = 1e-9)
  expect_false(balance_effects(abs(q), w)$ok)
})

test_that("P+ dispersion multiplier is mean-one around the trend", {
  set.seed(6)
  mu <- rep(100, 2e5); tr <- PCFG$disp_hetero$phi_inf + PCFG$disp_hetero$c / 100
  for (s in c(0, 0.4, 0.8, 1.6))
    expect_equal(mean(disp_hetero(mu, s)) / tr, 1, tolerance = 0.02 + 0.02 * s^2)
  expect_identical(disp_hetero(mu[1:3], 0), rep(tr, 3))
})

test_that("Part R option A: whole cohort, odd S drops one, halves balanced", {
  set.seed(7)
  G <- 2000; S <- 11
  L <- round(rlnorm(G, log(2000), 0.5))
  mu <- rlnorm(G, 4, 1.5)
  cnt <- matrix(rnbinom(G * S, mu = mu, size = 5), G, S, dimnames = list(paste0("g", 1:G), paste0("s", 1:S)))
  sim <- simulate_real(list(counts = cnt, lengths = L), 0.1)
  expect_identical(sim$n_per_group, 5L); expect_identical(sim$n_dropped, 1L)
  expect_identical(ncol(sim$counts), 10L)
  expect_identical(sim$delta, 0)
  expect_identical(sim$frac_sign_flip, 0)
  # realised injected effect tracks the intended one for well-expressed genes
  de <- sim$is_de & mu > 200
  expect_gt(cor(sim$realised_lfc[de], sim$true_lfc[de]), 0.95)
})
