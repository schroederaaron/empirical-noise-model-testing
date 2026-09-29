#!/usr/bin/env Rscript
# estimate_dispersion_trend.R
# -----------------------------------------------------------------------------
# Parameters of power_test.R Part P+ (plan I2), estimated from real cohorts:
#     phi_trend(mu) = phi_inf + c / mu        (mu = mean normalised count)
#     sigma_hat     = scatter of log gene-wise dispersion around that trend
#
#   Rscript estimate_dispersion_trend.R disp_params.rds coad_healthy.rds luad_healthy.rds kirc_healthy.rds
#   POWER_DISP_PARAMS=disp_params.rds POWER_PARTS=P+ Rscript power_test.R
#
# Inputs are cohort RDS files from export_power_cohort.R (plan: the healthy cohorts).
# Gene-wise dispersions are edgeR's UNSHRUNK estimates (estimateDisp, prior.df = 0);
# the trend is fitted to them with DESeq2's parametric method (Gamma GLM, identity
# link, iterated with outlier exclusion). sigma_hat = mad of log(phi_g / trend).
#
# CAVEAT (record with the estimate): gene-wise estimates carry sampling noise, so
# sigma_hat is an UPPER bound on the true heterogeneity. sigma_hat_corrected removes
# the expected sampling variance of a log dispersion estimate, trigamma((n - 1) / 2)
# (DESeq2's dispersion-prior-variance correction); the truth lies between the two.
# power_test.R uses sigma_hat (the plan brackets it with 0, 0.5x, 1x, 2x).
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: Rscript estimate_dispersion_trend.R <out.rds> <cohort.rds> [...]")
out <- args[1]; files <- args[-1]

.libPaths(c(normalizePath("external/docker_r_libs", mustWork = FALSE), .libPaths()))
suppressPackageStartupMessages(library(edgeR))

MIN_MEAN <- 10          # same expression floor as power_test.R Part R

#' DESeq2's parametricDispersionFit: phi = a0 + a1 / mu, Gamma GLM (identity link),
#' refitted on genes whose residual ratio lies within [1e-4, 15] until stable.
fit_trend <- function(mu, phi) {
  cf <- c(0.1, 1); good <- rep(TRUE, length(mu))
  for (it in 1:20) {
    old <- cf
    fit <- glm(phi[good] ~ I(1 / mu[good]), family = Gamma(link = "identity"), start = cf)
    cf <- coef(fit)
    if (any(cf <= 0)) stop("trend coefficients not positive: ", paste(signif(cf, 3), collapse = ", "))
    r <- phi / (cf[1] + cf[2] / mu); good <- r > 1e-4 & r < 15
    if (sum(log(cf / old)^2) < 1e-6) break
  }
  unname(cf)
}

one <- lapply(files, function(f) {
  d <- readRDS(f)
  cnt <- as.matrix(d$counts)
  if (nrow(cnt) < ncol(cnt)) stop(f, ": counts must be genes x samples")
  cnt <- cnt[rowMeans(cnt) >= MIN_MEAN, , drop = FALSE]
  y <- normLibSizes(DGEList(cnt))
  y <- estimateDisp(y, design = matrix(1, ncol(cnt), 1), prior.df = 0, tagwise = TRUE)
  eff <- y$samples$lib.size * y$samples$norm.factors
  mu <- rowMeans(sweep(cnt, 2, eff / mean(eff), "/"))            # mean normalised count
  phi <- y$tagwise.dispersion
  ok <- is.finite(phi) & phi > 1e-6                              # drop boundary (Poisson) fits
  label <- if (!is.null(d$project_id)) paste0(d$project_id, "/", d$stage) else basename(f)
  message(sprintf("%s: %d samples, %d genes (%d with phi_g > 1e-6)", label, ncol(cnt), nrow(cnt), sum(ok)))
  data.frame(cohort = label, n = ncol(cnt), mu = mu[ok], phi = phi[ok], lib_mean = mean(eff))
})
D <- do.call(rbind, one)

summ <- function(d) {
  cf <- fit_trend(d$mu, d$phi)
  lr <- log(d$phi / (cf[1] + cf[2] / d$mu))
  s_raw <- mad(lr)
  s_cor <- sqrt(max(0, s_raw^2 - mean(trigamma((d$n - 1) / 2))))
  data.frame(phi_inf = cf[1], c = cf[2], sigma_hat = s_raw, sigma_hat_corrected = s_cor,
             n_genes = nrow(d), lib_mean = mean(d$lib_mean))
}
per <- do.call(rbind, lapply(split(D, D$cohort), function(d) cbind(cohort = d$cohort[1], summ(d))))
pooled <- summ(D)
print(per, row.names = FALSE)
message(sprintf("POOLED: phi_inf = %.4f, c = %.3f, sigma_hat = %.3f (upper bound), sigma_hat_corrected = %.3f",
                pooled$phi_inf, pooled$c, pooled$sigma_hat, pooled$sigma_hat_corrected))

saveRDS(list(phi_inf = pooled$phi_inf, c = pooled$c, sigma_hat = pooled$sigma_hat,
             sigma_hat_corrected = pooled$sigma_hat_corrected, lib_mean = pooled$lib_mean,
             per_cohort = per, inputs = normalizePath(files), date = format(Sys.time())), out)
message("wrote ", out)
