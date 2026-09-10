#!/usr/bin/env Rscript
# reanalyse_null_calibration_sweep.R
# -----------------------------------------------------------------------------
# Re-reads an existing null_calibration run and extracts the three facts that
# constrain any explanation of TOX-raw's anti-conservatism. Runs in seconds on
# results already committed -- no data, no Fortran, no re-computation.
#
# Input: either the summary CSV null_calibration.R writes, or the captured stdout
# of a run (results/null_calibration_*.out), which is what is in the repository.
#
# The three readouts:
#   1. inflation at alpha = 0.05 AND at 0.01, per method. If a method is inflated
#      at 0.05 but not at 0.01, its null has the wrong SHAPE, not the wrong width:
#      a too-narrow null inflates both, and inflates 0.01 more.
#   2. inflation vs replicates per group. A fixed scale error is n-invariant
#      (observed and null both shrink as 1/sqrt(n)); an inflation that GROWS with n
#      means the observed statistic is Gaussianising by the CLT while the null
#      keeps the shape of individual residuals.
#   3. inflation vs kNN config, paired within cell. Bounds how much of the excess
#      any kNN retuning could possibly remove.
#
# Plus a two-parameter descriptive decomposition: model the null distance as
# lambda * (standardised t_nu) and the observed as N(0,1), then FPR at 0.05 and
# 0.01 identify lambda (scale: >1 = null too narrow) and nu (shape: small =
# leptokurtic, narrow core with heavy tails). Gaussian reference: lambda = 1,
# quantile ratio 1.314, excess kurtosis 0.
#
# RUN:  Rscript reanalyse_null_calibration_sweep.R [path-to-.out-or-.csv]
# -----------------------------------------------------------------------------

suppressMessages({library(dplyr)})
options(width = 200)

args <- commandArgs(trailingOnly = TRUE)
src  <- if (length(args)) args[1] else "../results/null_calibration_9414.out"
if (!file.exists(src)) stop("Input not found: ", src)

# --- read ---------------------------------------------------------------------
# The .out is the printed summary table: `cancer` and `stage` contain spaces, so
# it cannot be read as whitespace-delimited. The 13 trailing fields are numeric
# and k_config / n_per_group sit just before them, so parse from the RIGHT and
# treat whatever is left in the middle as cancer + stage.
NUMCOLS <- 13L
read_out <- function(path) {
  ln <- readLines(path, warn = FALSE)
  hi <- grep("^\\s*method\\b.*inflation_0\\.05", ln)
  if (!length(hi)) stop("No summary table header found in ", path)
  hdr <- strsplit(trimws(ln[hi[1]]), "\\s+")[[1]]
  keep <- c("TOX-log", "TOX-raw", "TOX-raw-boot", "TOX-log-boot",
            "TOX-raw-trim02", "TOX-raw-trim05", "edgeR", "limma", "DESeq2")
  out <- list()
  for (l in ln) {
    t <- strsplit(trimws(l), "\\s+")[[1]]
    if (length(t) < NUMCOLS + 4L || !(t[1] %in% keep)) next
    v <- suppressWarnings(as.numeric(replace(tail(t, NUMCOLS), tail(t, NUMCOLS) == "NA", NA)))
    if (all(is.na(v))) next
    mid <- t[2:(length(t) - NUMCOLS - 2L)]
    if (length(mid) >= 2L && mid[length(mid) - 1L] == "Stage") {
      stage <- paste(tail(mid, 2), collapse = " "); cancer <- paste(head(mid, -2), collapse = " ")
    } else { stage <- tail(mid, 1); cancer <- paste(head(mid, -1), collapse = " ") }
    out[[length(out) + 1L]] <- c(list(method = t[1], cancer = cancer, stage = stage,
                                      n_per_group = as.integer(t[length(t) - NUMCOLS - 1L]),
                                      k_config = t[length(t) - NUMCOLS]),
                                 setNames(as.list(v), tail(hdr, NUMCOLS)))
  }
  if (!length(out)) stop("No data rows parsed from ", path)
  do.call(rbind, lapply(out, as.data.frame, stringsAsFactors = FALSE))
}
df <- if (grepl("\\.csv$", src)) read.csv(src, stringsAsFactors = FALSE) else read_out(src)
names(df) <- sub("^FPR_0\\.05$", "FPR_0.05", names(df))
df$infl05 <- df$`FPR_0.05` / 0.05
df$infl01 <- df$`FPR_0.01` / 0.01
cat(sprintf("Parsed %d rows | %d methods | %d cancers | %d kNN configs\n",
            nrow(df), dplyr::n_distinct(df$method), dplyr::n_distinct(df$cancer),
            dplyr::n_distinct(df$k_config)))

# --- 1. shape signature: 0.05 vs 0.01 ----------------------------------------
cat("\n=== 1. Inflation at alpha = 0.05 vs 0.01 (median over all cells) ===\n")
cat("Inflated at .05 but NOT at .01  ->  wrong SHAPE (narrow core, heavy tail), not wrong width.\n")
print(df %>% group_by(method) %>%
        summarise(cells = dplyr::n(),
                  `FPR_.05` = round(median(`FPR_0.05`), 4), infl_05 = round(median(infl05), 3),
                  `FPR_.01` = round(median(`FPR_0.01`), 4), infl_01 = round(median(infl01), 3),
                  median_p = round(median(median_p), 3), ad_A2 = round(median(ad_A2, na.rm = TRUE), 1),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)

# --- 2. n-dependence ----------------------------------------------------------
cat("\n=== 2. Inflation vs replicates per group (one kNN config, to hold k fixed) ===\n")
cat("Rising with n  ->  the observed statistic Gaussianises (CLT) while the null does not.\n")
kref <- if ("k20_50" %in% df$k_config) "k20_50" else names(sort(table(df$k_config), TRUE))[1]
tox <- df %>% filter(grepl("^TOX", method), k_config == kref | is.na(k_config))
print(tox %>% mutate(bin = cut(n_per_group, c(0, 12, 25, 50, 100, Inf),
                               labels = c("<=12", "13-25", "26-50", "51-100", ">100"))) %>%
        group_by(method, bin) %>%
        summarise(cells = dplyr::n(), infl_05 = round(median(infl05), 3),
                  infl_01 = round(median(infl01), 3), median_p = round(median(median_p), 3),
                  ad_A2 = round(median(ad_A2, na.rm = TRUE), 1), .groups = "drop") %>%
        as.data.frame(), row.names = FALSE)

# --- 3. how much can k buy? ---------------------------------------------------
if (dplyr::n_distinct(df$k_config) > 1L) {
  cat("\n=== 3. kNN effect, paired within (cancer, stage, n) ===\n")
  cat("Delta = infl@.05 at the LARGEST config minus the SMALLEST. This is the ceiling on\n")
  cat("what any kNN retuning could remove; compare it to the excess (infl - 1).\n")
  # edgeR/limma/DESeq2 carry k_config "NA" (they are kNN-independent); drop anything
  # that is not an actual config label before picking the extremes.
  ks <- sort(unique(df$k_config[grepl("^k[0-9]+_[0-9]+$", df$k_config)]))
  if (length(ks) < 2L) ks <- character(0)
  # order by k_start then k_max, so "smallest" and "largest" are meaningful
  kk <- do.call(rbind, lapply(strsplit(sub("^k", "", ks), "_"), as.integer))
  ks <- ks[order(kk[, 1], kk[, 2])]
  lo <- ks[1]; hi <- ks[length(ks)]
  for (m in sort(unique(df$method[grepl("^TOX", df$method)]))) {
    a <- df %>% filter(method == m, k_config == lo) %>% select(cancer, stage, n_per_group, lo = infl05)
    b <- df %>% filter(method == m, k_config == hi) %>% select(cancer, stage, n_per_group, hi = infl05)
    j <- inner_join(a, b, by = c("cancer", "stage", "n_per_group"))
    if (!nrow(j)) next
    d <- j$hi - j$lo
    cat(sprintf("  %-16s %s -> %s : n=%3d  mean delta = %+0.3f  median = %+0.3f  frac>0 = %.2f  (median excess = %+0.3f)\n",
                m, lo, hi, length(d), mean(d), median(d), mean(d > 0),
                median(df$infl05[df$method == m], na.rm = TRUE) - 1))
  }
}

# --- 4. scale/shape decomposition --------------------------------------------
fit_scale_shape <- function(f05, f01) {
  if (!is.finite(f05) || !is.finite(f01) || f05 <= 0 || f01 <= 0 || f05 >= 1 || f01 >= 1)
    return(c(lambda = NA, nu = NA, excess_kurt = NA, ratio = NA))
  z1 <- qnorm(1 - f05 / 2); z2 <- qnorm(1 - f01 / 2); r <- z2 / z1
  sq <- function(nu, a) if (nu > 200) qnorm(1 - a / 2) else qt(1 - a / 2, nu) * sqrt((nu - 2) / nu)
  g <- function(nu) sq(nu, 0.01) / sq(nu, 0.05) - r
  if (g(2.5) * g(400) > 0) return(c(lambda = z1 / qnorm(0.975), nu = NA, excess_kurt = NA, ratio = r))
  nu <- uniroot(g, c(2.5, 400))$root
  c(lambda = z1 / sq(nu, 0.05), nu = nu, excess_kurt = if (nu > 4) 6 / (nu - 4) else Inf, ratio = r)
}
cat("\n=== 4. Scale / shape decomposition (lambda > 1 = null too NARROW) ===\n")
cat("Gaussian reference: lambda = 1, quantile ratio = 1.314, excess kurtosis = 0.\n")
cat("NA nu = the p-values imply a tail LIGHTER than normal, outside the t family.\n")
dec <- tox %>% mutate(bin = cut(n_per_group, c(0, 12, 25, 50, 100, Inf),
                                labels = c("<=12", "13-25", "26-50", "51-100", ">100"))) %>%
  group_by(method, bin) %>%
  summarise(f05 = median(`FPR_0.05`), f01 = median(`FPR_0.01`), cells = dplyr::n(), .groups = "drop")
res <- t(mapply(fit_scale_shape, dec$f05, dec$f01))
print(cbind(dec[, c("method", "bin", "cells")], round(as.data.frame(res), 3)), row.names = FALSE)
