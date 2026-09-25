# Reference: integrate() on the peak-scaled integrand, split at the mode, in log space.
.ref_logpmf <- function(y, eta, s) {
  g <- function(z) dnorm(z, log = TRUE) + dpois(y, exp(eta + s * z), log = TRUE)
  zh <- optimize(function(z) -g(z), c(-60, 60), tol = 1e-13)$minimum
  f <- function(z) exp(g(z) - g(zh))
  g(zh) + log(integrate(f, zh - 40, zh, rel.tol = 1e-13, subdivisions = 5000L)$value +
              integrate(f, zh, zh + 40, rel.tol = 1e-13, subdivisions = 5000L)$value)
}

test_that("PLN pmf matches integrate() on a grid", {
  for (mu in c(0.5, 5, 50, 1000)) for (s in c(0.1, 0.5, 1)) {
    eta <- pln_eta(mu, s)
    ys <- unique(round(c(0, 1, mu * c(0.3, 1, 3))))
    own <- pln_logpmf_eta(ys, rep(eta, length(ys)), s)
    ref <- vapply(ys, .ref_logpmf, 0, eta = eta, s = s)
    expect_lt(max(abs(own - ref)), 1e-9)                 # relative error of p < 1e-9
  }
})

test_that("PLN cdf matches integrate() and the summed pmf", {
  for (mu in c(0.5, 5, 50)) for (s in c(0.1, 0.5, 1)) {
    eta <- pln_eta(mu, s)
    for (q in unique(round(c(0, mu, 3 * mu)))) {
      ref <- integrate(function(z) dnorm(z) * ppois(q, exp(eta + s * z)), -40, 40,
                       rel.tol = 1e-13, subdivisions = 5000L)$value
      expect_lt(abs(pln_cdf_eta(q, eta, s) - ref), 1e-12)
      expect_lt(abs(pln_cdf_eta(q, eta, s) - sum(exp(pln_logpmf_eta(0:q, rep(eta, q + 1), s)))), 1e-12)
    }
  }
})
