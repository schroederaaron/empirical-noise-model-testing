#!/usr/bin/env Rscript
# null_construction_comparison.R
# -----------------------------------------------------------------------------
# Which null construction should the published tool use?
#
# Five ways to build the `own` null from the same kNN residual pools, compared on
# calibration AND on resolution, under a complete null (no DE anywhere):
#
#   exact pairs      current exact model: all |a - b| over sqrt(n)-scaled
#                    INDIVIDUAL residuals. Deterministic. Floor 1/(n_a*n_b + 1).
#   bootstrap iid    current bootstrap model: resample n_rep residuals iid from the
#                    POOLED neighbourhood, average, difference. Floor 1/(B + 1),
#                    plus RNG dependence.
#   saddlepoint      the B -> infinity limit of `bootstrap iid`, computed
#                    analytically (tox_pvalue_saddlepoint). No floor, no RNG.
#   bootstrap blocked  pick a neighbour GENE, then resample within that gene. Each
#                    null draw then carries ONE coherent noise level, the way a real
#                    gene's mean does.
#   blocked enumerated  the B -> infinity limit of `bootstrap blocked`, computed
#                    exactly by multiset enumeration (tox_pvalue_blocked_exact).
#                    No draws, no RNG, and a far lower floor.
#
# WHY `blocked` EXISTS AT ALL
#   An iid draw of n_rep residuals from the pool can combine one residual from a
#   quiet gene with one from a noisy gene. No real gene's mean is formed that way:
#   its n_rep values all come from a single noise level. That incoherence is a
#   property of the POOLING, not of the averaging -- the averaging (residuals of
#   both signs cancelling) is required, because the observed statistic is itself a
#   mean and cancels identically. Blocking fixes the incoherence and leaves the
#   averaging alone.
#
# WHY `blocked` CAN BE ENUMERATED
#   Resampling n of a gene's n residuals with replacement has only C(2n-1, n)
#   distinct outcomes: 10 at n_rep = 3, 126 at 5, 2.8e6 at 10. So the whole blocked
#   null is writable as k * C(2n-1, n) weighted values per side and scored with the
#   sorted-pool + binary-search tail count the exact model already has. Resolution
#   stops being a function of how long you are willing to sample.
#
# THE RESOLUTION ARITHMETIC (k = 30 neighbour genes, n_rep = 3, G = 12,000 genes)
#   BH at q = 0.05 needs the best p-value to reach q/G = 4.2e-6.
#     exact pairs        90 x 90        = 8,100 pairs      -> floor 1.2e-4
#     bootstrap B=25,000                                   -> floor 4.0e-5
#     bootstrap matching 4.2e-6                            -> needs B ~ 240,000
#     blocked enumerated 810 x 810      = 656,100          -> floor 1.5e-6
#   The combinatorial ceiling is NOT what binds a bootstrap here: distinct multisets
#   of 3 from 90 residuals number 125,580 per side, ~1.6e10 pairs of them, so
#   25,000 draws is nowhere near exhausting the space -- it is simply 25,000 draws,
#   and 1/(B+1) is the floor regardless.
#
# AND THE LIMIT NONE OF THEM REMOVE
#   A floor of 1.5e-6 is not 1.5e-6 worth of information. The effective sample size
#   behind a 90-residual pool at n_rep = 3 is of order 10^2 (see
#   docs/raw_normalisation_diagnosis.md 7b), so p-values far below ~1e-2 are
#   reporting which neighbourhood the gene drew. Removing an ARTIFICIAL floor that
#   binds above the information limit is worth doing; it does not create evidence.
#   That is why this script scores BH false positives under a complete null rather
#   than just reading the floor off a formula.
#
# RUN:  Rscript null_construction_comparison.R          (from TCGA_test/scripts)
# -----------------------------------------------------------------------------

.source_reimpl <- function() {
  cands <- c("tox_null_reimpl.R", "common/tox_null_reimpl.R", "../../common/tox_null_reimpl.R",
             "../common/tox_null_reimpl.R")
  for (p in cands) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Could not locate tox_null_reimpl.R")
}
.source_reimpl()
suppressMessages(library(parallel))
RNGkind("L'Ecuyer-CMRG"); set.seed(77)
options(width = 170)

G_GENES   <- 2500L
N_SAMPLES <- 40L
K_START   <- 30L; K_MAX <- 50L; TAU <- 0.1
B_DRAWS   <- 5000L
N_RUNS    <- 16L
CV_HET    <- c(0.0, 0.4, 0.7)     # variance heterogeneity left inside a mean-neighbourhood
N_REP     <- c(3L, 10L)
N_CORES   <- max(1L, min(8L, detectCores() - 1L))

#' H0 cohort: right-skewed means, per-gene CV dispersed by `cv_het`, gamma counts.
make_cohort <- function(G, N, cv_het, cv0 = 0.45) {
  mu <- exp(rnorm(G, 2.2, 1.7)); cv <- cv0 * exp(rnorm(G, 0, cv_het))
  m <- matrix(rgamma(N * G, shape = rep(1/cv^2, each = N), scale = rep(mu * cv^2, each = N)), nrow = N)
  colnames(m) <- paste0("g", seq_len(G)); m
}

METHODS <- c("exact pairs", "bootstrap iid", "saddlepoint", "bootstrap blocked", "blocked enumerated")

#' One split-half run: every method scored on the SAME pools and the SAME observed
#' statistics, so differences are attributable to the null construction alone.
one_run <- function(cv_het, n) {
  m  <- make_cohort(G_GENES, N_SAMPLES, cv_het)
  ix <- sample.int(nrow(m), 2*n); A <- m[ix[1:n], , drop = FALSE]; Bm <- m[ix[n+(1:n)], , drop = FALSE]
  sc <- tox_prepare_sorted(A, 0L); st <- tox_prepare_sorted(Bm, 0L)
  mc <- colMeans(A); mt <- colMeans(Bm); obs <- mc - mt
  P <- matrix(NA_real_, ncol(m), length(METHODS))
  for (g in seq_len(ncol(m))) {
    ga <- tox_gather(mc[g], sc, K_START, 1L, K_MAX, TAU)
    gb <- tox_gather(mt[g], st, K_START, 1L, K_MAX, TAU)
    if (length(ga$pool) < 10 || length(gb$pool) < 10) next
    va <- ga$pool; vb <- gb$pool
    ma <- sc$resid[, ga$genes, drop = FALSE]; mb <- st$resid[, gb$genes, drop = FALSE]
    P[g, 1] <- tox_pvalue_exact(va/sqrt(n), vb/sqrt(n), obs[g])
    a  <- colMeans(matrix(sample(va, n*B_DRAWS, TRUE), nrow = n))
    b  <- colMeans(matrix(sample(vb, n*B_DRAWS, TRUE), nrow = n))
    P[g, 2] <- (sum(abs(a-b) >= abs(obs[g])) + 1)/(B_DRAWS + 1)
    sp <- tox_pvalue_saddlepoint(va, vb, obs[g], n, n)
    P[g, 3] <- if (is.na(sp)) 1/(length(va)*length(vb) + 1) else sp
    ja <- sample(ncol(ma), B_DRAWS, TRUE); jb <- sample(ncol(mb), B_DRAWS, TRUE)
    a2 <- colMeans(matrix(ma[cbind(sample.int(n, n*B_DRAWS, TRUE), rep(ja, each = n))], nrow = n))
    b2 <- colMeans(matrix(mb[cbind(sample.int(n, n*B_DRAWS, TRUE), rep(jb, each = n))], nrow = n))
    P[g, 4] <- (sum(abs(a2-b2) >= abs(obs[g])) + 1)/(B_DRAWS + 1)
    # Enumeration is only tractable while C(2n-1, n) stays small; above ~n_rep 10
    # fall back to the sampled blocked null.
    P[g, 5] <- if (n <= 10L) tox_pvalue_blocked_exact(ma, mb, obs[g]) else P[g, 4]
  }
  vapply(seq_along(METHODS), function(j) {
    p <- P[, j]; p <- p[!is.na(p)]
    c(f05 = mean(p < 0.05), f01 = mean(p < 0.01), medp = median(p),
      minp = min(p), hits = sum(p.adjust(p, "BH") < 0.05), n = length(p))
  }, numeric(6))
}

cat("Complete-null simulation: no gene is DE, so every BH rejection is a false positive.\n")
cat(sprintf("G = %d genes, %d samples, k_start = %d / k_max = %d, B = %d draws, %d runs per cell.\n",
            G_GENES, N_SAMPLES, K_START, K_MAX, B_DRAWS, N_RUNS))
cat("Targets: inflation 1.00 at both levels, BH hits 0.\n\n")

for (cv_het in CV_HET) for (n in N_REP) {
  rr <- Filter(is.array, mclapply(seq_len(N_RUNS), function(s) one_run(cv_het, n), mc.cores = N_CORES))
  if (!length(rr)) next
  cat(sprintf("--- cv_het = %.1f, n_rep = %d  (floor: pairs %.2e | boot %.2e | blocked-enum %.2e) ---\n",
              cv_het, n, 1/((K_MAX*n)^2 + 1), 1/(B_DRAWS + 1), 1/((K_MAX*n^n)^2 + 1)))
  cat(sprintf("  %-20s %16s %16s %10s %12s %22s\n",
              "null construction", "infl@.05", "infl@.01", "median p", "median min p", "BH hits/run (max)"))
  for (j in seq_along(METHODS)) {
    f5 <- vapply(rr, function(x) x["f05", j], 0); f1 <- vapply(rr, function(x) x["f01", j], 0)
    md <- vapply(rr, function(x) x["medp", j], 0); mp <- vapply(rr, function(x) x["minp", j], 0)
    hh <- vapply(rr, function(x) x["hits", j], 0); NN <- mean(vapply(rr, function(x) x["n", j], 0))
    se5 <- 2*sqrt(mean(f5)*(1-mean(f5))/(NN*length(rr)))/0.05
    se1 <- 2*sqrt(mean(f1)*(1-mean(f1))/(NN*length(rr)))/0.01
    cat(sprintf("  %-20s %8.2f +-%-5.2f %8.2f +-%-5.2f %10.3f %12.2e %10.2f (%d)\n",
                METHODS[j], mean(f5)/0.05, se5, mean(f1)/0.01, se1, mean(md), median(mp),
                mean(hh), max(hh)))
  }
  cat("\n")
}
cat("READING\n")
cat("  infl@.01 is the column that matters most for a published tool: BH reads the tail.\n")
cat("  'BH hits/run' > 0 under a complete null is a false discovery the tool would report.\n")
cat("  A p-value floor is not automatically safe: genes PILE UP at it, and BH rejects a\n")
cat("  tied block of r genes as soon as r >= p_floor * G / q, so a coarse null can\n")
cat("  manufacture rejections rather than prevent them.\n")
