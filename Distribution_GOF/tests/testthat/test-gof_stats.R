test_that("CvM / AD / KS match goftest and stats::ks.test", {
  skip_if_not_installed("goftest")
  for (s in 1:5) {
    set.seed(s); u <- rbeta(500, 1.1, 0.9)
    expect_equal(cvm_u(u), unname(goftest::cvm.test(u, "punif")$statistic), tolerance = 1e-10)
    expect_equal(ad_u(u),  unname(goftest::ad.test(u, "punif")$statistic),  tolerance = 1e-10)
    expect_equal(ks_u(u),  unname(suppressWarnings(stats::ks.test(u, "punif"))$statistic), tolerance = 1e-12)
  }
})

test_that("gof_all reports N-normalised forms and the clamp count", {
  u <- c(0, 1e-15, seq(0.1, 0.9, 0.1), 1)
  s <- gof_all(u)
  expect_equal(s$N, length(u)); expect_equal(s$W2_N, s$W2 / s$N); expect_equal(s$A2_N, s$A2 / s$N)
  expect_equal(s$n_clamped, 3L)
})
