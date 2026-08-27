#!/usr/bin/env Rscript
# variance_space_comparison.R
# -----------------------------------------------------------------------------
# PROVES that the mean-variance relationship (which count-based tools rely on)
# exists in COUNT space but is DESTROYED in TPM space -- the space TOX builds its
# neighbourhoods in.
#
# Mechanism: TPM = counts * (1/gene_length) * scaling, so
#     log2(TPM) = log2(count) - log2(gene_length) + const.
# For a single gene the per-gene SD across samples is UNCHANGED (length is a
# constant), but its MEAN (x-position) is shifted by -log2(length). Gene length
# varies ~100x, so this reshuffles the x-ordering: the tight mean-SD trend in
# count/CPM space becomes a scattered cloud in TPM space, and genes at "equal TPM"
# have very different variances. That directly breaks TOX's assumption that
# variance is comparable at equal (TPM) expression.
#
# We reproduce the limma-voom-style mean-variance plot in BOTH spaces, on the same
# TOX-filtered gene set, at three aggregation levels:
#   (1) per cancer x stage        (2) per cancer (stages + healthy pooled)
#   (3) all cancers (coloured by cancer)
# Each panel: subsampled datapoints (spread) + LOESS trend. Plus a quantitative
# proof: the sd-heterogeneity metric (relative IQR of SD within mean-bins) in each
# space -- LOW in counts (trend holds), HIGH in TPM (trend broken).
#
# COUNT space  = log2( CPM + 1 )   (library-size normalised, like voom's input)
# TPM   space  = log2( TPM + 1 )
# so the ONLY difference is TPM's gene-length normalisation -- isolating the effect.
# -----------------------------------------------------------------------------

source("outlier_significance_analysis.R")   # loaders + compute_gene_keep_mask + globals
suppressMessages({library(dplyr); library(ggplot2); library(data.table)})
options(width = 250)

OUT_DIR  <- file.path(dirname(TOX_TEST_DIR), "variance_space_comparison")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

N_POINTS <- 5000L      # datapoints per panel (subsample)
MAX_TREND_PTS <- 40000L # cap points per group fed to tox_loess (the curve is identical
                        # far below this; keeps the fast netlib LOESS snappy on pooled groups)

# ==================== helpers ====================

col_sd <- function(m) {
  n  <- colSums(!is.na(m)); mu <- colMeans(m, na.rm = TRUE)
  sqrt(pmax((colMeans(m^2, na.rm = TRUE) - mu^2) * n / pmax(n - 1, 1), 0))
}

# median relative IQR of per-gene SD within mean-expression deciles: LOW => genes
# of similar mean share variance (trend holds); HIGH => they don't (trend broken).
sd_heterogeneity <- function(mean_e, sd_e) {
  ok <- is.finite(mean_e) & is.finite(sd_e); mean_e <- mean_e[ok]; sd_e <- sd_e[ok]
  if (length(mean_e) < 50) return(NA_real_)
  br <- unique(quantile(mean_e, 0:10 / 10, na.rm = TRUE)); if (length(br) < 3) return(NA_real_)
  bin <- cut(mean_e, br, include.lowest = TRUE)
  rel <- tapply(sd_e, bin, function(s) { m <- median(s, na.rm = TRUE)
    if (is.na(m) || m == 0) NA_real_ else IQR(s, na.rm = TRUE) / m })
  median(rel, na.rm = TRUE)
}

# per-gene (mean, sd) in a given log-space
gene_stats <- function(logmat) data.frame(mean_expr = colMeans(logmat, na.rm = TRUE),
                                          sd_expr = col_sd(logmat))
count_logcpm <- function(m) { ls <- rowSums(m, na.rm = TRUE); log2(m / ls * 1e6 + 1) }  # samples x genes
tpm_log      <- function(m) log2(m + 1)

# ---- loaders ----
load_counts <- function(project_id, stage = NULL, healthy = FALSE) {
  fn <- if (healthy) paste0("healthy_", project_id, "_counts.rds")
        else paste0(project_id, "-", gsub(" ", "-", stage), "_counts.rds")
  fp <- file.path(BASE_DATA_DIR, project_id, fn)
  if (!file.exists(fp)) return(NULL)
  d <- readRDS(fp); if (!is.matrix(d$expression_vectors)) return(NULL)
  m <- d$expression_vectors                                  # samples x genes
  if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
  m
}
load_tpm <- function(project_id, stage, healthy = FALSE) {
  d <- tryCatch(load_stage_data(project_id, if (healthy) STAGES[1] else stage,
                                data_type = if (healthy) "healthy" else "cancer",
                                use_constant_healthy = healthy,
                                norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
                error = function(e) NULL)
  if (is.null(d) || !is.matrix(d$expression_vectors)) return(NULL)
  m <- d$expression_vectors                                  # samples x genes (TPM)
  if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
  m
}

# ==================== gather per-gene stats in both spaces ====================

cancers   <- CANCER_TYPES
conditions <- c(STAGES, "Healthy")
base_rows <- list(); bi <- 1L
het_rows  <- list(); hi <- 1L

for (cname in names(cancers)) {
  pid <- unname(cancers[cname]); cat(sprintf("\n>>> %s (%s)\n", cname, pid))
  keep_mask <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL)
  if (is.null(keep_mask) || sum(keep_mask) == 0) { cat("    no keep mask; skip\n"); next }
  kept_ids <- names(keep_mask)[keep_mask]

  for (cond in conditions) {
    is_h <- identical(cond, "Healthy")
    cm <- load_counts(pid, if (is_h) NULL else cond, healthy = is_h)
    tm <- load_tpm(pid, cond, healthy = is_h)
    if (is.null(cm) || is.null(tm)) { cat(sprintf("    %-10s missing counts/TPM\n", cond)); next }

    common <- Reduce(intersect, list(kept_ids, colnames(cm), colnames(tm)))
    if (length(common) < 50) { cat(sprintf("    %-10s <50 common genes\n", cond)); next }
    cm <- cm[, common, drop = FALSE]; tm <- tm[, common, drop = FALSE]

    sc <- gene_stats(count_logcpm(cm))          # count/CPM space
    st <- gene_stats(tpm_log(tm))               # TPM space
    base_rows[[bi]] <- data.frame(cancer = cname, stage = cond, space = "counts (log2 CPM)",
                                  gene_id = common, mean_expr = sc$mean_expr, sd_expr = sc$sd_expr,
                                  stringsAsFactors = FALSE); bi <- bi + 1L
    base_rows[[bi]] <- data.frame(cancer = cname, stage = cond, space = "TPM (log2 TPM)",
                                  gene_id = common, mean_expr = st$mean_expr, sd_expr = st$sd_expr,
                                  stringsAsFactors = FALSE); bi <- bi + 1L

    het_c <- sd_heterogeneity(sc$mean_expr, sc$sd_expr)
    het_t <- sd_heterogeneity(st$mean_expr, st$sd_expr)
    het_rows[[hi]] <- data.frame(cancer = cname, stage = cond, n_genes = length(common),
                                 het_counts = round(het_c, 3), het_tpm = round(het_t, 3),
                                 ratio_tpm_over_counts = round(het_t / het_c, 2),
                                 stringsAsFactors = FALSE); hi <- hi + 1L
    cat(sprintf("    %-10s genes=%d  sd_heterog counts=%.3f  TPM=%.3f  (TPM/counts=%.2f)\n",
                cond, length(common), het_c, het_t, het_t / het_c))
  }
}

base <- as.data.frame(rbindlist(base_rows))
het  <- as.data.frame(rbindlist(het_rows))
base$space <- factor(base$space, levels = c("counts (log2 CPM)", "TPM (log2 TPM)"))
base$stage <- factor(base$stage, levels = conditions)
het$stage  <- factor(het$stage,  levels = conditions)

write.csv(het, file.path(OUT_DIR, "sd_heterogeneity_counts_vs_tpm.csv"), row.names = FALSE)
cat("\n", paste(rep("=", 90), collapse = ""), "\n", sep = "")
cat("QUANTITATIVE PROOF: sd-heterogeneity within mean-bins (LOW = trend holds)\n")
cat("counts should be LOW, TPM HIGH; ratio >> 1 confirms the trend is destroyed in TPM.\n")
cat(paste(rep("=", 90), collapse = ""), "\n", sep = "")
print(het[order(het$cancer, het$stage), ], row.names = FALSE)
cat(sprintf("\nMEDIAN sd-heterogeneity  counts=%.3f  TPM=%.3f  (median TPM/counts ratio=%.2f)\n",
            median(het$het_counts, na.rm = TRUE), median(het$het_tpm, na.rm = TRUE),
            median(het$ratio_tpm_over_counts, na.rm = TRUE)))

# ==================== plots ====================
cat("\n>>> writing plots to", PLOT_DIR, "\n")
set.seed(42)
# Subsample up to N_POINTS per panel: shuffle each group (slice_sample prop=1) then
# take the head. slice_head returns the whole group if it has fewer than N_POINTS,
# so this never errors -- unlike slice_sample(n=...), whose n must be a constant.
pts <- base %>% group_by(cancer, stage, space) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = N_POINTS) %>%
  ungroup()

# Precompute the LOESS trend with the project's netlib-backed tox_loess (much faster
# than R's stats::loess used by geom_smooth), default params. tox_loess returns yhat
# at the input points, so we sort by x per group and draw it with geom_line. One
# trend table per aggregation level to match each plot's facets.
fit_trend <- function(df, group_cols) {
  df %>%
    filter(is.finite(mean_expr), is.finite(sd_expr)) %>%
    group_by(across(all_of(group_cols))) %>%
    slice_sample(prop = 1) %>% slice_head(n = MAX_TREND_PTS) %>%
    arrange(mean_expr, .by_group = TRUE) %>%
    mutate(yhat = tox_loess(mean_expr, sd_expr)) %>%   # netlib LOESS, default params
    ungroup()
}
cat(">>> fitting LOESS trends via tox_loess (netlib)\n")
trend_space        <- fit_trend(base, "space")                        # global plot
trend_cancer_space <- fit_trend(base, c("cancer", "space"))           # per-cancer plot
trend_full         <- fit_trend(base, c("cancer", "stage", "space"))  # per-cancer-stage plots

# (3) ALL CANCERS, coloured by cancer, panels = space
p_global <- ggplot(pts, aes(mean_expr, sd_expr)) +
  geom_point(aes(colour = cancer), alpha = 0.18, size = 0.35) +
  geom_line(data = trend_space, aes(mean_expr, yhat), colour = "black", linewidth = 0.8) +
  facet_wrap(~ space, scales = "free_x") +
  labs(title = "Mean-variance relationship: counts vs TPM (all cancers)",
       subtitle = "points coloured by cancer; black = LOESS trend. Tight in counts, scattered in TPM.",
       x = "mean expression (log2)", y = "per-gene SD (log2)") +
  theme_minimal()
tryCatch(ggsave(file.path(PLOT_DIR, "mv_all_cancers.png"), p_global, width = 12, height = 6, dpi = 140),
         error = function(e) cat("   global plot failed:", conditionMessage(e), "\n"))

# (2) PER CANCER, stages+healthy pooled: rows = cancer, cols = space
p_cancer <- ggplot(pts, aes(mean_expr, sd_expr)) +
  geom_point(aes(colour = stage), alpha = 0.2, size = 0.35) +
  geom_line(data = trend_cancer_space, aes(mean_expr, yhat), colour = "black", linewidth = 0.7) +
  facet_grid(cancer ~ space, scales = "free_x") +
  labs(title = "Mean-variance relationship per cancer (stages + healthy pooled)",
       x = "mean expression (log2)", y = "per-gene SD (log2)", colour = "condition") +
  theme_minimal(base_size = 9) + theme(strip.text.y = element_text(angle = 0))
tryCatch(ggsave(file.path(PLOT_DIR, "mv_per_cancer.png"), p_cancer, width = 11, height = 20, dpi = 130),
         error = function(e) cat("   per-cancer plot failed:", conditionMessage(e), "\n"))

# (1) PER CANCER x STAGE: one figure per cancer, rows = stage, cols = space
for (cname in unique(base$cancer)) {
  tsub <- trend_full[trend_full$cancer == cname, ]; psub <- pts[pts$cancer == cname, ]
  p1 <- ggplot(psub, aes(mean_expr, sd_expr)) +
    geom_point(alpha = 0.2, size = 0.35, colour = "steelblue") +
    geom_line(data = tsub, aes(mean_expr, yhat), colour = "black", linewidth = 0.7) +
    facet_grid(stage ~ space, scales = "free_x") +
    labs(title = sprintf("Mean-variance: counts vs TPM  --  %s", cname),
         x = "mean expression (log2)", y = "per-gene SD (log2)") +
    theme_minimal(base_size = 9)
  fn <- sprintf("mv_%s.png", gsub("[^A-Za-z0-9]", "_", cname))
  tryCatch(ggsave(file.path(PLOT_DIR, fn), p1, width = 9, height = 12, dpi = 130),
           error = function(e) cat("   plot failed", cname, "\n"))
}

# quantitative-proof bar plot: sd-heterogeneity counts vs TPM
hl <- rbind(data.frame(cancer = het$cancer, stage = het$stage, space = "counts", het = het$het_counts),
            data.frame(cancer = het$cancer, stage = het$stage, space = "TPM",    het = het$het_tpm))
p_het <- ggplot(hl, aes(stage, het, fill = space)) +
  geom_col(position = "dodge") + facet_wrap(~ cancer, scales = "free_x") +
  labs(title = "SD-heterogeneity within mean-bins: counts vs TPM",
       subtitle = "higher = variance NOT comparable at equal expression (bad for TOX). TPM >> counts = proof.",
       x = "", y = "sd-heterogeneity (relative IQR)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
tryCatch(ggsave(file.path(PLOT_DIR, "sd_heterogeneity_counts_vs_tpm.png"), p_het, width = 12, height = 8, dpi = 140),
         error = function(e) NULL)

cat("\n", paste(rep("=", 90), collapse = ""), "\n", sep = "")
cat("DONE. CSV + plots in:", OUT_DIR, "\n")
cat(paste(rep("=", 90), collapse = ""), "\n", sep = "")
