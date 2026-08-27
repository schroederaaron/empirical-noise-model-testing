#!/usr/bin/env Rscript
# variance_structure.R
# -----------------------------------------------------------------------------
# Characterises the NOISE / VARIANCE STRUCTURE of the filtered TCGA data, per
# cancer x stage, to test the hypothesis that the stages where TOX fails have a
# DISTINCT variance structure -- not simply fewer replicates.
#
# TOX's empirical null assumes that expression-neighbouring genes share a noise
# distribution. It therefore breaks when, within a stage:
#   * the mean-variance relationship is anomalous, and/or
#   * genes of SIMILAR mean expression have very DIFFERENT variances
#     (heteroskedastic / heterogeneous residual variance) -- so pooling them
#     yields a mis-specified null.
#
# For each (cancer, stage) [+ the constant healthy reference] on the SAME
# cross-stage-filtered gene set TOX uses, we compute per-gene mean & SD in raw
# and log2 space and summarise:
#   n_rep, n_genes, median mean/SD, CV, a variance-HETEROGENEITY metric (relative
#   spread of per-gene SD within expression deciles), and the skewness/kurtosis of
#   the pooled log-residuals.
# Plots: mean-SD trend per stage (the money plot), SD densities, and a full
#   mean-SD 2D-density faceted cancer x stage.
#
# Reuses the pipeline loaders/filter by sourcing outlier_significance_analysis.R
# (its main block is `if (sys.nframe()==0)`-guarded, so nothing runs on source).
# -----------------------------------------------------------------------------

source("outlier_significance_analysis.R")   # loaders + compute_gene_keep_mask + globals
suppressMessages({library(dplyr); library(ggplot2); library(data.table)})
options(width = 250, max.print = 1e6)

OUT_DIR  <- file.path(dirname(TOX_TEST_DIR), "variance_structure")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

INCLUDE_HEALTHY <- TRUE   # add the constant healthy reference as an extra "stage"

# ==================== per-gene / per-stat helpers ====================

# Vectorised per-column (per-gene) SD with NA handling and n/(n-1) correction.
col_sd <- function(m) {
  n  <- colSums(!is.na(m))
  mu <- colMeans(m, na.rm = TRUE)
  v  <- (colMeans(m^2, na.rm = TRUE) - mu^2) * n / pmax(n - 1, 1)
  sqrt(pmax(v, 0))
}
skewness <- function(x) { x <- x[is.finite(x)]; n <- length(x)
  if (n < 3) return(NA_real_); s <- sd(x); if (s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^3) }
kurtosis <- function(x) { x <- x[is.finite(x)]; n <- length(x)   # excess kurtosis
  if (n < 4) return(NA_real_); s <- sd(x); if (s == 0) return(NA_real_)
  mean(((x - mean(x)) / s)^4) - 3 }

# Relative spread of per-gene SD WITHIN expression deciles -- the key metric:
# high => genes of similar mean have very different variances, violating TOX's
# expression-neighbourhood homogeneity assumption.
sd_heterogeneity <- function(mean_log, sd_log) {
  ok <- is.finite(mean_log) & is.finite(sd_log)
  mean_log <- mean_log[ok]; sd_log <- sd_log[ok]
  if (length(mean_log) < 50) return(NA_real_)
  br  <- unique(quantile(mean_log, probs = 0:10 / 10, na.rm = TRUE))
  if (length(br) < 3) return(NA_real_)
  bin <- cut(mean_log, br, include.lowest = TRUE)
  rel <- tapply(sd_log, bin, function(s) {
    m <- median(s, na.rm = TRUE); if (is.na(m) || m == 0) NA_real_ else IQR(s, na.rm = TRUE) / m
  })
  median(rel, na.rm = TRUE)
}

# ==================== gather stats over cancer x stage ====================

cancers <- CANCER_TYPES                       # named vector: label -> project id
combo_stages <- STAGES
if (INCLUDE_HEALTHY) combo_stages <- c(STAGES, "Healthy")

gene_rows <- list(); gi <- 1L                 # per-gene (for plots)
summ_rows <- list(); si <- 1L                 # per (cancer, stage) summary

for (cname in names(cancers)) {
  pid <- unname(cancers[cname])
  cat(sprintf("\n>>> %s (%s)\n", cname, pid))

  keep_mask <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL)
  if (is.null(keep_mask) || sum(keep_mask) == 0) { cat("    no keep mask; skip\n"); next }
  kept_ids <- names(keep_mask)[keep_mask]

  for (stg in combo_stages) {
    is_healthy <- identical(stg, "Healthy")
    d <- tryCatch(load_stage_data(pid,
                                  if (is_healthy) STAGES[1] else stg,
                                  data_type = if (is_healthy) "healthy" else "cancer",
                                  use_constant_healthy = is_healthy,
                                  norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || !is.matrix(d$expression_vectors)) { cat(sprintf("    %s: no data\n", stg)); next }

    m <- d$expression_vectors                 # samples x genes (raw)
    if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
    common <- intersect(kept_ids, colnames(m))
    if (length(common) < 50) { cat(sprintf("    %s: <50 kept genes\n", stg)); next }
    m <- m[, common, drop = FALSE]
    n_rep <- nrow(m)

    # raw-space and log2-space per-gene mean / SD
    mean_raw <- colMeans(m, na.rm = TRUE); sd_raw <- col_sd(m)
    lm       <- log2(m + 1)
    mean_log <- colMeans(lm, na.rm = TRUE);  sd_log <- col_sd(lm)
    cv_raw   <- sd_raw / pmax(mean_raw, .Machine$double.eps)

    # pooled log-residuals (mean-removed per gene) for shape stats
    resid <- sweep(lm, 2, mean_log, "-")
    rv <- as.vector(resid)

    gene_rows[[gi]] <- data.frame(
      cancer = cname, stage = stg, n_rep = n_rep,
      mean_raw = mean_raw, sd_raw = sd_raw,
      mean_log = mean_log, sd_log = sd_log, cv_raw = cv_raw,
      stringsAsFactors = FALSE); gi <- gi + 1L

    summ_rows[[si]] <- data.frame(
      cancer = cname, stage = stg, n_rep = n_rep, n_genes = length(common),
      med_mean_log = round(median(mean_log, na.rm = TRUE), 3),
      med_sd_log   = round(median(sd_log, na.rm = TRUE), 3),
      q90_sd_log   = round(quantile(sd_log, 0.9, na.rm = TRUE), 3),
      med_cv_raw   = round(median(cv_raw, na.rm = TRUE), 3),
      sd_heterogeneity = round(sd_heterogeneity(mean_log, sd_log), 3),
      resid_skew   = round(skewness(rv), 3),
      resid_kurt   = round(kurtosis(rv), 3),
      stringsAsFactors = FALSE); si <- si + 1L
    cat(sprintf("    %-10s n_rep=%3d  med_sd_log=%.3f  sd_heterog=%.3f  kurt=%.2f\n",
                stg, n_rep, median(sd_log, na.rm = TRUE),
                sd_heterogeneity(mean_log, sd_log), kurtosis(rv)))
  }
}

genes <- as.data.frame(rbindlist(gene_rows))
summ  <- as.data.frame(rbindlist(summ_rows))
summ$stage <- factor(summ$stage, levels = combo_stages)
genes$stage <- factor(genes$stage, levels = combo_stages)

write.csv(summ, file.path(OUT_DIR, "variance_structure_summary.csv"), row.names = FALSE)
cat("\n", paste(rep("=", 90), collapse = ""), "\n", sep = "")
cat("VARIANCE-STRUCTURE SUMMARY (per cancer x stage)\n")
cat(paste(rep("=", 90), collapse = ""), "\n", sep = "")
print(summ[order(summ$cancer, summ$stage), ], row.names = FALSE)

# ==================== plots ====================
cat("\n>>> writing plots to", PLOT_DIR, "\n")

# (P1) THE MONEY PLOT: mean-SD trend (log space), one loess per stage, faceted by
#      cancer. If a failing stage's trend or spread departs from its siblings, it
#      shows here.
p1 <- ggplot(genes, aes(mean_log, sd_log, colour = stage)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4, linewidth = 0.7) +
  facet_wrap(~ cancer, scales = "free") +
  labs(title = "Mean-SD relationship (log2 space) per stage",
       subtitle = "loess trend of per-gene SD vs mean; diverging stage = distinct variance structure",
       x = "mean log2(expr+1)", y = "per-gene SD (log2)") +
  theme_minimal()
tryCatch(ggsave(file.path(PLOT_DIR, "mean_sd_trend_log.png"), p1, width = 13, height = 9, dpi = 140),
         error = function(e) cat("   P1 failed:", conditionMessage(e), "\n"))

# (P2) same in RAW space (the pp-plots you analysed were raw norm)
p2 <- ggplot(genes, aes(mean_raw, sd_raw, colour = stage)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.4, linewidth = 0.7) +
  facet_wrap(~ cancer, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Mean-SD relationship (raw space, log-log) per stage",
       x = "mean raw expr (log10)", y = "per-gene SD (log10)") +
  theme_minimal()
tryCatch(ggsave(file.path(PLOT_DIR, "mean_sd_trend_raw.png"), p2, width = 13, height = 9, dpi = 140),
         error = function(e) cat("   P2 failed:", conditionMessage(e), "\n"))

# (P3) distribution of per-gene SD (log space) by stage, faceted by cancer
p3 <- ggplot(genes, aes(sd_log, colour = stage)) +
  geom_density(linewidth = 0.6) +
  facet_wrap(~ cancer, scales = "free") +
  labs(title = "Per-gene SD distribution (log2) by stage",
       x = "per-gene SD (log2)", y = "density") +
  theme_minimal()
tryCatch(ggsave(file.path(PLOT_DIR, "sd_density_log.png"), p3, width = 13, height = 9, dpi = 140),
         error = function(e) cat("   P3 failed:", conditionMessage(e), "\n"))

# (P4) full mean-SD 2D density, faceted cancer x stage (the detailed view)
p4 <- ggplot(genes, aes(mean_log, sd_log)) +
  geom_bin2d(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "genes") +
  facet_grid(cancer ~ stage, scales = "free") +
  labs(title = "Mean vs SD (log2) density, per cancer x stage",
       x = "mean log2(expr+1)", y = "per-gene SD (log2)") +
  theme_minimal(base_size = 8) + theme(strip.text.y = element_text(angle = 0))
tryCatch(ggsave(file.path(PLOT_DIR, "mean_sd_density_grid.png"), p4, width = 12, height = 18, dpi = 130),
         error = function(e) cat("   P4 failed:", conditionMessage(e), "\n"))

# (P5) heterogeneity & kurtosis bars per cancer x stage (the two summary metrics)
sm <- summ[summ$stage != "Healthy" | !INCLUDE_HEALTHY, ]
p5 <- ggplot(sm, aes(stage, sd_heterogeneity, fill = stage)) +
  geom_col() + facet_wrap(~ cancer, scales = "free_x") +
  labs(title = "Variance heterogeneity within expression bands (higher = worse for TOX)",
       subtitle = "median relative IQR of per-gene SD within expression deciles",
       x = "", y = "SD heterogeneity") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
tryCatch(ggsave(file.path(PLOT_DIR, "sd_heterogeneity.png"), p5, width = 12, height = 8, dpi = 140),
         error = function(e) cat("   P5 failed:", conditionMessage(e), "\n"))

cat("\n", paste(rep("=", 90), collapse = ""), "\n", sep = "")
cat("DONE. Summary CSV + plots in:", OUT_DIR, "\n")
cat(paste(rep("=", 90), collapse = ""), "\n", sep = "")
