#!/usr/bin/env Rscript
# variance_space_robustness.R
# -----------------------------------------------------------------------------
# Stress-tests the "TPM is fine" conclusion from variance_space_comparison.R
# against three ways it could have been hiding the effect. Only if TPM still wins
# under ALL THREE is the space concern truly dead.
#
#   (1) FILTER regime. compute_gene_keep_mask strips low-expression genes -- exactly
#       the regime where 1/mu (Poisson noise -> the count mean-variance trend)
#       matters. We recompute the metric on the EXCLUDED (low-expression) genes.
#   (2) PSEUDOCOUNT. Both spaces used log2(x+1); +1 compresses SD at low values, and
#       at fixed mu TPM ~ 1/L so long genes sit lower and take MORE compression --
#       possibly manufacturing TPM's advantage. We rerun with pc = 0.1 and on the
#       TPM > 1 subset (where the pseudocount barely matters).
#   (3) BIN WIDTH. The heterogeneity metric used deciles (thousands of genes/bin);
#       the real kNN neighbourhoods are k ~ 8-15 genes. We recompute the relative
#       IQR at the ACTUAL k -- the number that describes what the model does.
#
# Metric: median over bins of (IQR(SD) / median(SD)) within equal-frequency
#   mean-expression bins. Bin size is set by `resolution`: "decile" -> 10 bins;
#   a number -> that many genes per bin (== a kNN window of that width).
# ratio = het_tpm / het_counts.  ratio < 1 => TPM MORE homogeneous (space fine);
#   ratio > 1 => TPM MORE heterogeneous (concern alive).
# -----------------------------------------------------------------------------

source("outlier_significance_analysis.R")   # loaders + compute_gene_keep_mask + globals
suppressMessages({library(dplyr); library(data.table)})
options(width = 250)

OUT_DIR <- file.path(dirname(TOX_TEST_DIR), "variance_space_comparison")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PCS         <- c(1, 0.1)                 # pseudocounts to test
RESOLUTIONS <- list(decile = "decile", k15 = 15, k8 = 8)   # bin widths (genes/bin)
GENE_SETS   <- c("kept", "excluded")     # filter regimes
TPM_FILTERS <- c("all", "tpm>1")         # remove pseudocount-compression at low end

# ==================== helpers ====================

col_sd <- function(m) {
  n  <- colSums(!is.na(m)); mu <- colMeans(m, na.rm = TRUE)
  sqrt(pmax((colMeans(m^2, na.rm = TRUE) - mu^2) * n / pmax(n - 1, 1), 0))
}

# relative-IQR heterogeneity at a given resolution (bin width). "decile" -> 10 bins;
# numeric -> ~that many genes per bin (a kNN window of that width).
het_metric <- function(mean_e, sd_e, resolution) {
  ok <- is.finite(mean_e) & is.finite(sd_e); mean_e <- mean_e[ok]; sd_e <- sd_e[ok]
  n <- length(mean_e); if (n < 60) return(NA_real_)
  nb <- if (identical(resolution, "decile")) 10L else max(2L, as.integer(round(n / resolution)))
  if (n < 2L * nb) return(NA_real_)
  br <- unique(quantile(mean_e, seq(0, 1, length.out = nb + 1), na.rm = TRUE))
  if (length(br) < 3) return(NA_real_)
  bin <- cut(mean_e, br, include.lowest = TRUE)
  rel <- tapply(sd_e, bin, function(s) { m <- median(s, na.rm = TRUE)
    if (is.na(m) || m == 0) NA_real_ else IQR(s, na.rm = TRUE) / m })
  median(rel, na.rm = TRUE)
}

load_counts <- function(project_id, stage = NULL, healthy = FALSE) {
  fn <- if (healthy) paste0("healthy_", project_id, "_counts.rds")
        else paste0(project_id, "-", gsub(" ", "-", stage), "_counts.rds")
  fp <- file.path(BASE_DATA_DIR, project_id, fn)
  if (!file.exists(fp)) return(NULL)
  d <- readRDS(fp); if (!is.matrix(d$expression_vectors)) return(NULL)
  m <- d$expression_vectors
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
  m <- d$expression_vectors
  if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
  m
}

# ==================== sweep ====================

cancers    <- CANCER_TYPES
conditions <- c(STAGES, "Healthy")
rows <- list(); ri <- 1L

for (cname in names(cancers)) {
  pid <- unname(cancers[cname]); cat(sprintf("\n>>> %s (%s)\n", cname, pid))
  keep_mask <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL)
  if (is.null(keep_mask)) { cat("    no keep mask; skip\n"); next }
  kept_ids     <- names(keep_mask)[keep_mask]      # well-expressed (TOX's set)
  excluded_ids <- names(keep_mask)[!keep_mask]     # low-expression (filtered out)

  for (cond in conditions) {
    is_h <- identical(cond, "Healthy")
    cm <- load_counts(pid, if (is_h) NULL else cond, healthy = is_h)
    tm <- load_tpm(pid, cond, healthy = is_h)
    if (is.null(cm) || is.null(tm)) { cat(sprintf("    %-10s missing\n", cond)); next }
    common <- intersect(colnames(cm), colnames(tm))
    cm <- cm[, common, drop = FALSE]; tm <- tm[, common, drop = FALSE]

    ls       <- rowSums(cm, na.rm = TRUE)          # library size per sample (all genes)
    cpm      <- cm / ls * 1e6
    mean_tpm <- colMeans(tm, na.rm = TRUE)         # raw TPM mean, for the tpm>1 filter

    for (pc in PCS) {
      lc <- log2(cpm + pc); lt <- log2(tm + pc)
      mc <- colMeans(lc, na.rm = TRUE); sc <- col_sd(lc)
      mt <- colMeans(lt, na.rm = TRUE); st <- col_sd(lt)
      for (gs in GENE_SETS) {
        base_sel <- common %in% (if (gs == "kept") kept_ids else excluded_ids)
        for (tf in TPM_FILTERS) {
          sel <- base_sel & (if (tf == "tpm>1") mean_tpm > 1 else TRUE)
          if (sum(sel) < 60) next
          for (rn in names(RESOLUTIONS)) {
            rv <- RESOLUTIONS[[rn]]
            hc <- het_metric(mc[sel], sc[sel], rv)
            ht <- het_metric(mt[sel], st[sel], rv)
            rows[[ri]] <- data.frame(cancer = cname, stage = cond, gene_set = gs,
                                     pseudocount = pc, tpm_filter = tf, resolution = rn,
                                     n = sum(sel), het_counts = round(hc, 4),
                                     het_tpm = round(ht, 4), ratio = round(ht / hc, 3),
                                     stringsAsFactors = FALSE); ri <- ri + 1L
          }
        }
      }
    }
  }
}

res <- as.data.frame(rbindlist(rows))
write.csv(res, file.path(OUT_DIR, "variance_space_robustness.csv"), row.names = FALSE)

# ==================== summary: median over cancer x stage per scenario ====================

summ <- res %>%
  group_by(gene_set, pseudocount, tpm_filter, resolution) %>%
  summarise(n_combos = sum(!is.na(ratio)),
            med_het_counts = round(median(het_counts, na.rm = TRUE), 3),
            med_het_tpm    = round(median(het_tpm, na.rm = TRUE), 3),
            med_ratio      = round(median(ratio, na.rm = TRUE), 3),
            frac_TPM_worse = round(mean(ratio > 1, na.rm = TRUE), 2),   # ratio>1 = TPM more heterogeneous
            .groups = "drop") %>%
  arrange(gene_set, resolution, tpm_filter, pseudocount) %>%
  as.data.frame()

cat("\n", paste(rep("=", 110), collapse = ""), "\n", sep = "")
cat("VARIANCE-SPACE ROBUSTNESS  (median over cancer x stage)\n")
cat("ratio = het_TPM / het_counts.  <1 => TPM more homogeneous (space fine);  >1 => TPM worse (concern alive)\n")
cat("frac_TPM_worse = fraction of cancer x stage combos with ratio > 1\n")
cat(paste(rep("=", 110), collapse = ""), "\n", sep = "")
print(summ, row.names = FALSE)

verdict <- function(lbl, sub) {
  if (nrow(sub) == 0) { cat(sprintf("\n%s: no data\n", lbl)); return(invisible()) }
  cat(sprintf("\n%s -> med ratio %.2f, TPM worse in %.0f%% of combos  [%s]\n",
              lbl, median(sub$med_ratio), 100 * median(sub$frac_TPM_worse),
              if (median(sub$med_ratio) > 1) "CONCERN ALIVE" else "space still fine"))
}
cat("\n---- the three objections ----")
verdict("(1) EXCLUDED (low-expr) genes, pc=1, k15",
        summ[summ$gene_set=="excluded" & summ$pseudocount==1 & summ$tpm_filter=="all" & summ$resolution=="k15", ])
verdict("(2a) KEPT, pc=0.1, k15",
        summ[summ$gene_set=="kept" & summ$pseudocount==0.1 & summ$tpm_filter=="all" & summ$resolution=="k15", ])
verdict("(2b) KEPT, TPM>1, pc=1, k15",
        summ[summ$gene_set=="kept" & summ$pseudocount==1 & summ$tpm_filter=="tpm>1" & summ$resolution=="k15", ])
verdict("(3) KEPT, pc=1, real k (k15 vs decile above)",
        summ[summ$gene_set=="kept" & summ$pseudocount==1 & summ$tpm_filter=="all" & summ$resolution=="k15", ])

cat("\n", paste(rep("=", 110), collapse = ""), "\n", sep = "")
cat("DONE. Full grid: variance_space_robustness.csv in", OUT_DIR, "\n")
cat(paste(rep("=", 110), collapse = ""), "\n", sep = "")
