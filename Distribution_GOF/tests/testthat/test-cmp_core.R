test_that("cmp_solve gives E[Y] = mu to 1e-10 relative", {
  for (mu in c(0.5, 5, 50, 500, 5e3, 5e4)) for (nu in c(0.1, 0.5, 1, 2, 3)) {
    s <- cmp_solve(mu, nu); y <- s[["lo"]]:s[["hi"]]
    p <- exp(cmp_logpmf(y, list(mu = rep(mu, length(y)), nu = nu,
                                cmp_sol = matrix(s, length(y), length(s), byrow = TRUE, dimnames = list(NULL, names(s))))))
    expect_lt(abs(sum(y * p) - mu) / mu, 1e-10)
    expect_lt(abs(sum(p) - 1), 1e-12)
  }
})

test_that("nu = 1 reduces CMP to Poisson", {
  par <- list(mu = c(0.7, 12, 300), nu = 1)
  y <- c(0, 10, 290)
  expect_equal(cmp_logpmf(y, par), dpois(y, par$mu, log = TRUE), tolerance = 1e-9)
  expect_equal(cmp_cdf(y, par), ppois(y, par$mu), tolerance = 1e-9)
})
