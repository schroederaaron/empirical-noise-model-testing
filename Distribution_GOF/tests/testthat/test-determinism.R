# G7: same seed => identical T_obs and T_b[1:5] (exact equality).
test_that("G7: a seeded bootstrap is exactly reproducible", {
  set.seed(11); n <- 12; G <- 8
  X <- model.matrix(~ 1, data.frame(i = 1:n))
  Y <- t(sapply(exp(seq(log(3), log(300), length.out = G)), function(m) rnbinom(n, mu = m, size = 4)))
  dimnames(Y) <- list(paste0("g", 1:G), paste0("s", 1:n))
  cfg <- GOF_CONFIG; cfg$B <- 5L
  run <- function() {
    O <- compute_offsets(Y, "tmm")
    out <- list()
    for (m in c("poisson", "nb")) {
      fam <- gof_families(cfg, m)[[1]]
      f <- fit_dataset(Y, X, O, fam, rownames(Y), 1L)
      pb <- pit_bounds(Y, f, fam, rownames(Y)[f$ok])
      tob <- gof_all(randomise_pit(pb$a, pb$b, seed_for(cfg$BASE_SEED, "pit_gboot", m)))
      ref <- list(fit = f, Y = Y, O = O, X = X, lengths = NULL, genes = rownames(Y), obs_group = NULL)
      bt <- run_bootstrap(m, fam, ref, cfg, B = 5L, n_cores = 1L)
      out[[m]] <- list(T_obs = unlist(tob[c("W2", "A2", "D")]), T_b = as.matrix(bt$T_b[1:5, .(W2, A2, D)]))
    }
    out
  }
  r1 <- run(); r2 <- run()
  expect_identical(r1, r2)
})
