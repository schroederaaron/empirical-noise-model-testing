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
              "disp_hetero", "simulate_power", "simulate_real", "%||%", "check_disp_params",
              "build_jobs", "scan_radius", "mat_hash", "ALL_PARTS", "REAL_PARTS", "SEED_PARTS"), E)
PCFG <- list(frac_up = 0.5, lfc_range = c(0.25, 4), lib_size = 4e7, disp = 0.2,
             disp_hetero = list(phi_inf = 0.1, c = 1, sigma_hat = 0.8),
             R = list(lfc_tol = 0.5, pi1 = 0.10), seed = 20260923L, n_rounds = 3L,
             RS = list(n_reps = c(3L, 5L), pi1 = c(0, 0.10), n_rounds = 4L))

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

# ---------------------------------------------------------------- note D (Part RS)

fake_cohort <- function(G = 1500, S = 30, seed = 8) {
  set.seed(seed)
  cnt <- matrix(rnbinom(G * S, mu = rlnorm(G, 4, 1.5), size = 5), G, S,
                dimnames = list(paste0("g", 1:G), sprintf("TCGA-AA-%04d-01A", 1:S)))
  list(counts = cnt, lengths = round(rlnorm(G, log(2000), 0.5)), label = "fake")
}

test_that("RS sampler: sizes n, groups disjoint, samples from the cohort, reproducible, 2n > S refused", {
  co <- fake_cohort()
  for (n in c(3L, 5L)) {
    set.seed(11); a <- simulate_real(co, 0.1, n = n)
    set.seed(11); b <- simulate_real(co, 0.1, n = n)
    expect_identical(ncol(a$counts), 2L * n); expect_identical(a$n_per_group, n)
    expect_identical(sum(a$group == "A"), n); expect_identical(sum(a$group == "B"), n)
    expect_identical(anyDuplicated(a$sample_ids), 0L)            # disjoint groups
    expect_true(all(a$sample_ids %in% colnames(co$counts)))
    expect_identical(a$counts, b$counts); expect_identical(a$sample_ids, b$sample_ids)
  }
  expect_error(simulate_real(co, 0.1, n = 16L), "2n = 32 > S = 30")
})

test_that("RS pi1 = 0: no DE, counts untouched, delta 0, DE diagnostics NA", {
  co <- fake_cohort()
  set.seed(12); sim <- simulate_real(co, 0, n = 3L)
  expect_false(any(sim$is_de)); expect_true(all(sim$true_lfc == 0))
  expect_identical(sim$delta, 0)
  expect_identical(unname(sim$counts), unname(co$counts[, sim$sample_ids]))
  expect_true(is.na(sim$eff_lfc_median_abs)); expect_true(is.na(sim$frac_lfc_mismatch))
})

test_that("RS jobs: seeds unique across all jobs and parts; pi1 = 0 and 0.1 jobs differ; 2n > S cells skipped", {
  PCFG$parts <<- c("RS", "R", "P")
  PCFG$P <<- list(dists = "nb", n_reps = c(3L, 5L), pi1 = 0.1)
  REAL_COHORTS <<- list(big = list(counts = matrix(0, 2, 30), label = "big"),
                        small = list(counts = matrix(0, 2, 8), label = "small"))
  j <- build_jobs()
  run <- j[!j$skip, ]
  expect_identical(anyDuplicated(run$seed), 0L)
  rs <- run[run$part == "RS", ]
  expect_identical(nrow(rs), 2L * 2L * 4L + 1L * 2L * 4L)          # big: n 3, 5; small: n 3 only
  expect_identical(rs$job, seq_len(nrow(rs)))                       # one running counter
  expect_false(any(rs$seed[rs$pi1 == 0] %in% rs$seed[rs$pi1 > 0]))
  sk <- j[j$skip, ]
  expect_identical(nrow(sk), 2L); expect_true(all(sk$cohort == "small" & sk$n_rep == 5L))
  # Part R keeps its numbering (cohort, round) and its seeds when RS is added.
  r <- run[run$part == "R", ]
  expect_identical(r$job, seq_len(2L * PCFG$n_rounds))
  expect_identical(r$seed, PCFG$seed + 7919L * r$job + 1000003L * 1L)
  # Seed part index is frozen: P keeps index 2 although RS now sits before it in ALL_PARTS.
  p <- run[run$part == "P", ]
  expect_identical(p$seed, PCFG$seed + 7919L * p$job + 1000003L * 2L)
})

test_that("WP2 guard: stale or missing mu_depth refused, matching file accepted", {
  ok <- list(phi_inf = 0.12, c = 19, sigma_hat = 1, mu_depth = 4e7)
  expect_identical(check_disp_params(ok, 4e7), ok)
  expect_error(check_disp_params(ok[-4], 4e7), "mu_depth = missing")
  expect_error(check_disp_params(modifyList(ok, list(mu_depth = 1e7)), 4e7), "stale")
  expect_error(check_disp_params(modifyList(ok, list(mu_depth = 4e7 * (1 + 1e-6))), 4e7), "stale")
})

test_that("scan_radius equals a brute-force n-th nearest neighbour in linear mean", {
  set.seed(13)
  m <- c(rlnorm(400, 3, 2), 0, 0, 5, 5)                            # ties included
  n <- sample.int(60, length(m), replace = TRUE); n[3] <- NA
  brute <- vapply(seq_along(m), function(i) {
    if (is.na(n[i])) return(NA_real_)
    d <- abs(m - m[i]); j <- order(d)[n[i]]
    abs(log2(m[j] + 1) - log2(m[i] + 1))
  }, numeric(1))
  r <- scan_radius(m, n)
  # With tied distances the brute force may pick either equidistant gene, so the
  # comparison is restricted to genes whose n-th distance is unique.
  dl <- vapply(seq_along(m), function(i) if (is.na(n[i])) NA_real_ else sort(abs(m - m[i]))[n[i]], numeric(1))
  uniq <- vapply(seq_along(m), function(i) !is.na(n[i]) && sum(abs(abs(m - m[i]) - dl[i]) < 1e-9) == 1L, logical(1))
  expect_true(is.na(r[3]))
  expect_gt(sum(uniq), 300)
  expect_equal(r[uniq], brute[uniq], tolerance = 1e-9)
})

test_that("data_hash: same dataset -> same hash, different round -> different hash", {
  co <- fake_cohort()
  set.seed(14); a <- simulate_real(co, 0.1, n = 3L)
  set.seed(14); b <- simulate_real(co, 0.1, n = 3L)
  set.seed(15); c <- simulate_real(co, 0.1, n = 3L)
  expect_identical(mat_hash(a$counts), mat_hash(b$counts))
  expect_false(identical(mat_hash(a$counts), mat_hash(c$counts)))
})
