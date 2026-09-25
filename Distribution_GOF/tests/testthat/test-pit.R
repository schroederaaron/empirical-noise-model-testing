test_that("PIT bounds are exact at y = 0 and consistent with the pmf", {
  fams <- gof_families(GOF_CONFIG)
  pars <- list(poisson = list(mu = c(0.3, 5)), nb = list(mu = c(0.3, 5), size = 2),
               zinb = list(mu = c(0.3, 5), size = 2, pi = 0.2), genpois = list(mu = c(0.3, 5), phi = 3),
               pln = list(mu = c(0.3, 5), sigma = 0.5), cmp = list(mu = c(0.3, 5), nu = 0.7))
  for (m in names(pars)) {
    f <- fams[[m]]; par <- pars[[m]]; y <- c(0, 0)
    a <- f$cdf(y - 1, par); b <- f$cdf(y, par); p <- exp(f$logpmf(y, par))
    expect_equal(a, c(0, 0), info = m)                          # F(-1) = 0
    expect_equal(b, p, tolerance = 1e-12, info = m)             # F(0) = p(0)
  }
})

test_that("randomised PIT lies in [a, b] and is reproducible from its seed", {
  a <- matrix(c(0, 0.2, 0.5), 1); b <- matrix(c(0.1, 0.2, 0.9), 1)
  u1 <- randomise_pit(a, b, 7L); u2 <- randomise_pit(a, b, 7L)
  expect_identical(u1, u2)
  expect_true(all(u1 >= a & u1 <= b))
})
