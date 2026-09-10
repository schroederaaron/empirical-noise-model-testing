#!/usr/bin/env Rscript
# pvalue_tail_methods.R
# -----------------------------------------------------------------------------
# Can TOX get a DERIVED tail instead of a counted one, the way edgeR/limma/DESeq2 do?
#
# HOW THE REFERENCE TOOLS AVOID A RESOLUTION FLOOR
#   All three share one structure:
#     1. borrow information across genes to estimate the noise parameter
#        - limma: empirical-Bayes shrinkage of gene-wise residual variances toward
#          a prior fitted across all genes; with trend = TRUE / voom, the prior is
#          expression-dependent (a LOESS mean-variance trend)
#        - edgeR: NB dispersion shrunk toward a mean-dependent trend, then a
#          quasi-likelihood dispersion shrunk the same way
#        - DESeq2: dispersion shrunk toward a mean-dependent trend
#     2. standardise the effect by that borrowed noise estimate
#     3. read the tail off a THEORETICAL distribution whose df is inflated by the
#        borrowing: moderated t on d_g + d_0 df (limma), F(1, d_g + d_0) (edgeR QL),
#        Wald normal (DESeq2)
#   Step 3 is analytic, so p = 1e-30 is arithmetic, not resolution. Nothing is
#   counted, so nothing has a floor.
#
#   TOX ALREADY DOES STEP 1. The kNN mean-neighbourhood is the same device as
#   limma's expression-dependent prior and edgeR's mean-dependent trend -- pooling
#   noise information from genes at comparable expression. What TOX does instead of
#   steps 2-3 is count pairs, which caps the p-value at 1/(n_a*n_b + 1).
#
# THE OPTION THAT KEEPS TOX DISTRIBUTION-FREE
#   A theoretical tail does not have to mean an assumed one. The saddlepoint
#   (Lugannani-Rice) approximation reads the tail off the pool's OWN empirical
#   cumulant generating function, so the reference distribution is still the data's,
#   but the tail is evaluated analytically. And because it is written for the
#   MEAN-difference statistic, it simultaneously supplies the CLT shape the observed
#   statistic has -- the shape mismatch diagnosed in
#   docs/raw_normalisation_diagnosis.md. One change, both problems.
#
# THIS SCRIPT compares four ways of getting the same tail probability:
#   exact pairs   what TOX does now              floor 1/(n_a*n_b + 1)
#   bootstrap     the "bootstrap" model          floor 1/(B + 1), plus RNG noise
#   saddlepoint   tox_pvalue_saddlepoint()       no floor, deterministic
#   GPD           peaks-over-threshold fit       no floor; standard for permutation
#                                                tests (Knijnenburg et al. 2009)
#
# Part A validates against a null whose answer is known in closed form.
# Part B asks the question that actually decides whether this is worth doing:
#   removing the counting floor is only useful if the POOL behind it is solid.
#
# RUN:  Rscript pvalue_tail_methods.R          (from TCGA_test/scripts)
# -----------------------------------------------------------------------------

# tox_null_reimpl.R lives in common/; the repo's scripts are run from more than one
# working directory, so look for it rather than assuming one.
# --- persistent R package library (MUST precede any library()/require()/source() call) ---
# `.libPaths()` silently DROPS non-existent directories, so the folder has to be created
# FIRST or the prepend is a no-op and packages land in the ephemeral container library.
LIB_DIR <- normalizePath("external/docker_r_libs", mustWork = FALSE)
if (!dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE) && !dir.exists(LIB_DIR))
  stop("Could not create package library ", LIB_DIR,
       " -- it must be on a WRITABLE, BIND-MOUNTED path or packages will not persist.")
.libPaths(c(LIB_DIR, .libPaths()))

.source_reimpl <- function() {
  cands <- c("tox_null_reimpl.R", "common/tox_null_reimpl.R", "../../common/tox_null_reimpl.R",
             "../common/tox_null_reimpl.R", "../../../common/tox_null_reimpl.R")
  for (p in cands) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Could not locate tox_null_reimpl.R (looked in: ", paste(cands, collapse = ", "), ")")
}
.source_reimpl()
options(width = 150)
set.seed(8)

#' Generalised-Pareto peaks-over-threshold tail, for comparison.
#' Keeps the empirical null below quantile `q` and extrapolates past it.
gpd_tail <- function(t, nulls, q = 0.90) {
  u  <- unname(quantile(abs(nulls), q)); ex <- abs(nulls)[abs(nulls) > u] - u
  if (length(ex) < 50 || t <= u) return(mean(abs(nulls) >= t))
  nll <- function(par) { xi <- par[1]; sg <- exp(par[2])
    z <- 1 + xi * ex / sg; if (any(z <= 0)) return(1e10)
    length(ex) * log(sg) + (1 + 1/xi) * sum(log(z)) }
  f <- tryCatch(optim(c(0.1, log(mean(ex))), nll), error = function(e) NULL)
  if (is.null(f)) return(mean(abs(nulls) >= t))
  xi <- f$par[1]; sg <- exp(f$par[2])
  (length(ex)/length(nulls)) * (1 + xi * (t - u)/sg)^(-1/xi)
}

#' A residual pool with per-gene scale heterogeneity `het` (0 = single scale).
sim_pool <- function(k, nrep, het) {
  sd_g <- exp(rnorm(k, 0, het))
  as.numeric(sapply(sd_g, function(s) { x <- rnorm(nrep, 0, s); (x - mean(x)) * sqrt(nrep/(nrep - 1)) }))
}

n <- 3L; k <- 50L

# ---------------------------------------------------------------- Part A ----
cat("=== A1. Far-tail accuracy where the answer is known in closed form ===\n")
cat("Pools are large iid N(0,1) samples, so D = mean(n) - mean(n) is exactly N(0, 2/n).\n")
cat("A bounded ratio far into the tail is the property that makes this usable;\n")
cat("an Edgeworth/normal expansion would drift without limit instead.\n\n")
a <- rnorm(4000); b <- rnorm(4000); sdD <- sqrt(2/n)
cat(sprintf("%-10s %11s %14s %14s %8s\n", "true p", "threshold", "saddlepoint", "analytic", "ratio"))
for (tp in c(1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12)) {
  t <- qnorm(1 - tp/2) * sdD
  ps <- tox_pvalue_saddlepoint(a, b, t, n, n); pa <- 2 * pnorm(-t/sdD)
  cat(sprintf("%-10.0e %11.4f %14.3e %14.3e %8.2f\n", tp, t, ps, pa, ps/pa))
}

cat("\n=== A2. The four methods on a realistic heterogeneous pool ===\n")
cat(sprintf("n_rep = %d, k = %d genes -> %d residuals/side, %d pairs, floor = %.2e\n",
            n, k, n*k, (n*k)^2, 1/((n*k)^2 + 1)))
cat("Reference = 2e7 bootstrap mean-differences from the SAME pools.\n")
cat("'exact pairs' is scored on sqrt-scaled INDIVIDUAL residuals, as the model does now,\n")
cat("so its drift away from the reference is the shape mismatch, not an error in the count.\n\n")
pa_ <- sim_pool(k, n, 0.5); pb_ <- sim_pool(k, n, 0.5)
Bref <- 2e7
ref <- abs(colMeans(matrix(sample(pa_, n*Bref, TRUE), nrow = n)) -
           colMeans(matrix(sample(pb_, n*Bref, TRUE), nrow = n)))
boot1e4 <- sample(ref, 1e4)
cat(sprintf("%-10s %10s %13s %13s %13s %13s\n",
            "true p", "threshold", "exact pairs", "boot 1e4", "saddlepoint", "GPD(1e4)"))
for (tp in c(5e-2, 1e-2, 1e-3, 1e-4, 1e-5)) {
  t <- unname(quantile(ref, 1 - tp))
  cat(sprintf("%-10.0e %10.4f %13.3e %13.3e %13.3e %13.3e\n", tp, t,
              tox_pvalue_exact(pa_/sqrt(n), pb_/sqrt(n), t),
              (sum(boot1e4 >= t) + 1)/(1e4 + 1),
              tox_pvalue_saddlepoint(pa_, pb_, t, n, n),
              gpd_tail(t, boot1e4)))
}

# ---------------------------------------------------------------- Part B ----
cat("\n=== B. What an analytic tail does NOT fix ===\n")
cat("The saddlepoint evaluates the tail OF A GIVEN POOL to arbitrary precision. It cannot\n")
cat("reduce uncertainty in the pool. 200 independent neighbourhoods, same true threshold:\n\n")
het <- 0.5
A0 <- sim_pool(20000, n, het); B0 <- sim_pool(20000, n, het)
ref2 <- abs(colMeans(matrix(sample(A0, n*4e6, TRUE), nrow = n)) -
            colMeans(matrix(sample(B0, n*4e6, TRUE), nrow = n)))
cat(sprintf("%-10s %12s %24s %12s %14s\n", "true p", "median p", "IQR", "5-95% span", "within 2x"))
for (tp in c(1e-2, 1e-3, 1e-4)) {
  t <- unname(quantile(ref2, 1 - tp))
  ps <- replicate(200, tox_pvalue_saddlepoint(sim_pool(k, n, het), sim_pool(k, n, het), t, n, n))
  ps <- ps[is.finite(ps) & ps > 0]
  cat(sprintf("%-10.0e %12.2e %10.2e - %-11.2e %11.0fx %13.0f%%\n",
              tp, median(ps), quantile(ps, .25), quantile(ps, .75),
              quantile(ps, .95)/quantile(ps, .05), 100*mean(ps > tp/2 & ps < tp*2)))
}

cat("\nREADING\n")
cat("  A: the counting floor is removable. The saddlepoint tracks a known tail within\n")
cat("     ~1.4x to 1e-12, and on a real pool it is still tracking at 1e-5, where the\n")
cat("     pairwise count and a 1e4 bootstrap have both bottomed out.\n")
cat("  B: at n_rep = 3 removing it is not obviously an improvement. Redrawing the\n")
cat("     neighbourhood moves the answer by orders of magnitude, so an analytic tail\n")
cat("     replaces a visible floor with a number that looks precise and is not.\n")
cat("     The floor is at least honest about where the information stops.\n")
cat("  => adopt the analytic tail where the pool is well estimated (the real\n")
cat("     case-vs-control analysis, tens to hundreds of samples). At n_rep = 3, report\n")
cat("     a bound rather than a small number, whichever tail method is used.\n")
