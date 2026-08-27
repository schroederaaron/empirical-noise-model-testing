#!/usr/bin/env Rscript
# detailed_comparison.R
# -----------------------------------------------------------------------------
# POST-HOC, in-depth comparison from SAVED data -- no TOX / DE re-run.
#
# Loads:
#   * all_models_pvalues.rds   (TOX per-gene p-values, written by compare_noise_models.R)
#   * <method>_results.csv     (edgeR / voom / noiseq, written by deg_comparison.R)
#
# Produces the most detailed picture possible, split by cancer x stage x
# normalization (raw / log):
#
#   (1) P-VALUE CORRELATION DOT-PLOTS. For every model comparison, the reference
#       (or other TOX model) p-value on x vs the TOX p-value on y, as a 2D-density
#       "dot-plot" (geom_bin2d over ~19k genes), faceted cancer x stage, one figure
#       per (comparison, normalization). Plus a Spearman-correlation table per
#       (comparison, cancer, stage, norm).
#         - exact vs bootstrap        (TOX internal)
#         - {exact, bootstrap} vs {edgeR, voom, noiseq}
#
#   (2) STAGE-WISE + CANCER-WISE RECALL / PRECISION for the exact and bootstrap
#       models vs edgeR, voom (limma), noiseq (NOISeqBIO) AND the consensus sets
#       (consensus2 = >=2 methods, consensus3 = all 3), swept over
#         alpha    in {0.01, 0.05}
#         distance in {0.05, 0.10, 0.25}
#       Reported per (cancer, stage) [finest], per cancer, and overall, as CSVs +
#       printed tables + precision-recall plots.
#
# Reference "hit" rule (matches deg_comparison / load_ref_hits): |logFC| > 2 AND
#   edgeR FDR < DEG_FDR | voom adj.P.Val < DEG_FDR | noiseq prob > DEG_NOISEQ_PROB.
# -----------------------------------------------------------------------------

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(data.table)
})
source("config.R")   # standalone: provides TOX_TEST_DIR, STAGES, etc.

# Print wide tables on a single line instead of wrapping into blocks at 80 cols
# (the default). Raise if any table still wraps; view the .out with line-wrap OFF.
options(width = 250, max.print = 1e6)

# ==================== CONFIGURATION ====================

OUT_BASE <- file.path(dirname(TOX_TEST_DIR), "model_comparison")
TOX_RDS  <- file.path(OUT_BASE, "all_models_pvalues.rds")
DE_DIR   <- file.path(dirname(TOX_TEST_DIR), "differential_expression")
OUT_DIR  <- file.path(dirname(TOX_TEST_DIR), "detailed_comparison")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

TOX_MODELS  <- c("bootstrap", "exact")            # evaluated (kept if present in data)
REF_METHODS <- c("edgeR", "voom", "noiseq")       # count-parametric x2 + NOISeqBIO
REF_ALL     <- c(REF_METHODS, "consensus2", "consensus3")
NORMS       <- c("raw", "log")
ALPHAS      <- c(0.01, 0.05)                       # noise BH-adjusted cutoffs
DISTS       <- c(0.05, 0.10, 0.25)                 # effect-size (distance p) cutoffs
COMPARE_COMP<- "own_healthy"

DEG_FDR <- 0.01; DEG_LOGFC <- 2; DEG_NOISEQ_PROB <- 0.95
PLOG_CAP <- 12                                      # cap for -log10(p) axes

# ==================== HELPERS ====================

jaccard <- function(a, b) if (length(a) == 0 && length(b) == 0) NA_real_ else
  length(intersect(a, b)) / length(union(a, b))

nlog10 <- function(p) pmin(-log10(pmax(p, 10^(-PLOG_CAP))), PLOG_CAP)

#' Read one reference CSV -> data.frame(gene_id, ref_p, ref_logFC, ref_hit).
#' ref_p is the method's raw p-value (noiseq: 1 - prob, a confidence proxy).
read_ref <- function(cancer, stage, method) {
  fn <- file.path(DE_DIR, cancer, stage, paste0(method, "_results.csv"))
  if (!file.exists(fn)) return(NULL)
  d <- read.csv(fn, stringsAsFactors = FALSE)
  if (!all(c("gene_id", "logFC") %in% colnames(d))) return(NULL)
  if (method == "edgeR") {
    if (!all(c("PValue", "FDR") %in% colnames(d))) return(NULL)
    pv <- d$PValue; hit <- d$FDR < DEG_FDR
  } else if (method == "voom") {
    if (!all(c("P.Value", "adj.P.Val") %in% colnames(d))) return(NULL)
    pv <- d$P.Value; hit <- d$adj.P.Val < DEG_FDR
  } else {  # noiseq
    if (!all(c("PValue", "prob") %in% colnames(d))) return(NULL)
    pv <- d$PValue; hit <- d$prob > DEG_NOISEQ_PROB
  }
  hit <- hit & (abs(d$logFC) > DEG_LOGFC)
  data.frame(gene_id = as.character(d$gene_id), ref_p = as.numeric(pv),
             ref_logFC = as.numeric(d$logFC), ref_hit = (hit & !is.na(hit)),
             stringsAsFactors = FALSE)
}

# ==================== LOAD ====================

cat(">>> Loading TOX p-values:", TOX_RDS, "\n")
if (!file.exists(TOX_RDS)) stop("Saved TOX rds not found: ", TOX_RDS)
tox <- as.data.frame(readRDS(TOX_RDS), stringsAsFactors = FALSE)
tox <- tox[tox$comparison == COMPARE_COMP &
           !is.na(tox$noise_p_value) & !is.na(tox$distance_p_value), ]

present_models <- intersect(TOX_MODELS, unique(tox$model))
norms   <- intersect(NORMS, unique(tox$norm_method))
cancers <- sort(unique(tox$cancer_id))
stages  <- sort(unique(tox$stage))
combos  <- unique(tox[, c("cancer_id", "stage")])
cat(sprintf("    models=%s  norms=%s  cancers=%d  stages=%d  rows=%d\n",
            paste(present_models, collapse = ","), paste(norms, collapse = ","),
            length(cancers), length(stages), nrow(tox)))

cat(">>> Loading reference DE results + building consensus sets\n")
# ref_pval[[method]] : data.frame(cancer_id, stage, gene_id, ref_p)  -- for scatters
# ref_hits[[i]]      : list(method -> gene set, consensus2, consensus3)  -- for recall
ref_pval <- setNames(vector("list", length(REF_METHODS)), REF_METHODS)
ref_hits <- vector("list", nrow(combos))
for (i in seq_len(nrow(combos))) {
  pid <- combos$cancer_id[i]; stg <- combos$stage[i]
  hitsets <- list()
  for (m in REF_METHODS) {
    rr <- read_ref(pid, stg, m)
    if (is.null(rr)) { hitsets[[m]] <- character(0); next }
    ref_pval[[m]][[length(ref_pval[[m]]) + 1L]] <-
      data.frame(cancer_id = pid, stage = stg, gene_id = rr$gene_id, ref_p = rr$ref_p,
                 stringsAsFactors = FALSE)
    hitsets[[m]] <- unique(rr$gene_id[rr$ref_hit])
  }
  cnt <- table(unlist(hitsets[REF_METHODS], use.names = FALSE))
  hitsets$consensus2 <- names(cnt)[cnt >= 2]
  hitsets$consensus3 <- names(cnt)[cnt >= 3]
  ref_hits[[i]] <- hitsets
}
ref_pval <- lapply(ref_pval, function(x) if (length(x)) rbindlist(x) else NULL)
present_refs <- REF_ALL[vapply(REF_ALL, function(m)
  any(vapply(ref_hits, function(s) length(s[[m]]) > 0, logical(1))), logical(1))]
if (!any(REF_METHODS %in% present_refs)) stop("No usable reference DE result files found.")

# ==================== (1) P-VALUE CORRELATION DOT-PLOTS ====================

cat("\n>>> (1) p-value correlation dot-plots + correlation table\n")

#' Density dot-plot of x_p (reference/other model) vs y_p (TOX model), both as
#' -log10(p), faceted cancer x stage. `pdat` has cols cancer_id, stage, xp, yp.
pp_dotplot <- function(pdat, xlab, ylab, title, fname) {
  if (is.null(pdat) || nrow(pdat) == 0) return(invisible())
  pdat$xl <- nlog10(pdat$xp); pdat$yl <- nlog10(pdat$yp)
  p <- ggplot(pdat, aes(xl, yl)) +
    geom_bin2d(bins = 60) +
    scale_fill_viridis_c(trans = "log10", name = "genes") +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2, linewidth = 0.3) +
    facet_grid(cancer_id ~ stage) +
    labs(title = title, x = sprintf("-log10(%s)", xlab), y = sprintf("-log10(%s)", ylab)) +
    theme_minimal(base_size = 8) +
    theme(strip.text.y = element_text(angle = 0))
  tryCatch(ggsave(file.path(PLOT_DIR, fname), p, width = 11, height = 16, dpi = 130),
           error = function(e) cat("    plot failed:", fname, ":", conditionMessage(e), "\n"))
}

corr_rows <- list(); ci <- 1L
add_corr <- function(pdat, comparison, nm) {
  if (is.null(pdat) || nrow(pdat) == 0) return(invisible())
  s <- pdat %>% group_by(cancer_id, stage) %>%
    summarise(n = sum(is.finite(xp) & is.finite(yp)),
              spearman = suppressWarnings(cor(xp, yp, method = "spearman",
                                              use = "pairwise.complete.obs")),
              pearson_log = suppressWarnings(cor(nlog10(xp), nlog10(yp),
                                              use = "pairwise.complete.obs")),
              .groups = "drop")
  s$comparison <- comparison; s$norm_method <- nm
  corr_rows[[ci]] <<- s; ci <<- ci + 1L
}

for (nm in norms) {
  # ---- TOX vs TOX (bootstrap vs exact) ----
  if (all(c("bootstrap", "exact") %in% present_models)) {
    tb <- tox[tox$norm_method == nm & tox$model == "bootstrap",
              c("cancer_id", "stage", "gene_id", "noise_p_value")]
    te <- tox[tox$norm_method == nm & tox$model == "exact",
              c("cancer_id", "stage", "gene_id", "noise_p_value")]
    names(tb)[4] <- "xp"; names(te)[4] <- "yp"
    pdat <- merge(tb, te, by = c("cancer_id", "stage", "gene_id"))
    pp_dotplot(pdat, "bootstrap p", "exact p",
               sprintf("TOX bootstrap vs exact noise p  (norm=%s)", nm),
               sprintf("pp_TOX_bootstrap_vs_exact_%s.png", nm))
    add_corr(pdat, "exact~bootstrap", nm)
  }
  # ---- TOX vs reference ----
  for (model in present_models) {
    tm <- tox[tox$norm_method == nm & tox$model == model,
              c("cancer_id", "stage", "gene_id", "noise_p_value")]
    names(tm)[4] <- "yp"
    for (method in REF_METHODS) {
      rp <- ref_pval[[method]]
      if (is.null(rp)) next
      rpd <- as.data.frame(rp); names(rpd)[names(rpd) == "ref_p"] <- "xp"
      pdat <- merge(tm, rpd, by = c("cancer_id", "stage", "gene_id"))
      pp_dotplot(pdat, sprintf("%s p", method), sprintf("TOX %s p", model),
                 sprintf("%s p vs TOX %s  (norm=%s)", method, model, nm),
                 sprintf("pp_%s_vs_TOX_%s_%s.png", method, model, nm))
      add_corr(pdat, sprintf("%s~%s", model, method), nm)
    }
  }
}
corr_tbl <- if (length(corr_rows)) as.data.frame(rbindlist(corr_rows)) else data.frame()
if (nrow(corr_tbl)) {
  corr_tbl <- corr_tbl[order(corr_tbl$comparison, corr_tbl$norm_method,
                             corr_tbl$cancer_id, corr_tbl$stage), ]
  write.csv(corr_tbl, file.path(OUT_DIR, "pvalue_correlation_by_cancer_stage.csv"), row.names = FALSE)
  cat("\n-- p-value Spearman correlation, summary (median over cancer x stage) --\n")
  print(as.data.frame(
    corr_tbl %>% group_by(comparison, norm_method) %>%
      summarise(median_spearman = round(median(spearman, na.rm = TRUE), 3),
                min_spearman = round(min(spearman, na.rm = TRUE), 3),
                max_spearman = round(max(spearman, na.rm = TRUE), 3), .groups = "drop")),
    row.names = FALSE)
}

# ==================== (2) RECALL / PRECISION (stage- & cancer-wise) ====================

cat("\n>>> (2) recall / precision vs references + consensus (alpha x distance sweep)\n")

# Per-combo building blocks (additive across combos: keys disjoint per combo).
N <- nrow(combos) * length(norms) * length(present_models) * length(ALPHAS) *
     length(DISTS) * length(present_refs)
v_can <- character(N); v_stg <- character(N); v_norm <- character(N); v_mod <- character(N)
v_a <- numeric(N); v_d <- numeric(N); v_ref <- character(N)
v_nt <- integer(N); v_nr <- integer(N); v_sh <- integer(N)
z <- 0L
for (i in seq_len(nrow(combos))) {
  pid <- combos$cancer_id[i]; stg <- combos$stage[i]
  cb <- tox[tox$cancer_id == pid & tox$stage == stg, ]
  for (nm in norms) {
    cbn <- cb[cb$norm_method == nm, ]
    for (model in present_models) {
      cm <- cbn[cbn$model == model, ]
      for (a in ALPHAS) for (dt in DISTS) {
        # which() drops NA so an NA adjusted-p can never inject an NA "hit".
        q  <- unique(cm$gene_id[which(cm$noise_p_value_adj < a & cm$distance_p_value < dt)])
        nq <- length(q)
        for (refm in present_refs) {
          rk <- ref_hits[[i]][[refm]]
          z <- z + 1L
          v_can[z] <- pid; v_stg[z] <- stg; v_norm[z] <- nm; v_mod[z] <- model
          v_a[z] <- a; v_d[z] <- dt; v_ref[z] <- refm
          v_nt[z] <- nq; v_nr[z] <- length(rk); v_sh[z] <- length(intersect(q, rk))
        }
      }
    }
  }
}
detail <- data.frame(cancer_id = v_can, stage = v_stg, norm_method = v_norm, model = v_mod,
                     alpha = v_a, distance = v_d, reference = v_ref,
                     n_tox = v_nt, n_ref = v_nr, shared = v_sh, stringsAsFactors = FALSE)

metrics <- function(d, group_cols) {
  gv <- c(group_cols, "norm_method", "model", "alpha", "distance", "reference")
  a  <- aggregate(d[c("n_tox", "n_ref", "shared")], by = d[gv], FUN = sum)
  a$recall    <- round(a$shared / pmax(a$n_ref, 1), 4)
  a$precision <- round(a$shared / pmax(a$n_tox, 1), 4)
  a$jaccard   <- round(a$shared / pmax(a$n_tox + a$n_ref - a$shared, 1), 4)
  a$f1        <- round(2 * a$recall * a$precision / pmax(a$recall + a$precision, 1e-9), 4)
  a
}
rp_stage  <- metrics(detail, c("cancer_id", "stage"))   # finest: per cancer x stage
rp_cancer <- metrics(detail, "cancer_id")               # per cancer (summed over stages)
rp_overall<- metrics(detail, character(0))

write.csv(rp_stage,   file.path(OUT_DIR, "recall_precision_per_cancer_stage.csv"), row.names = FALSE)
write.csv(rp_cancer,  file.path(OUT_DIR, "recall_precision_per_cancer.csv"), row.names = FALSE)
write.csv(rp_overall, file.path(OUT_DIR, "recall_precision_overall.csv"), row.names = FALSE)

# ---- printed tables: exhaustive, per reference, split raw/log ----
for (nm in norms) {
  for (refm in present_refs) {
    cat(sprintf("\n============ RECALL/PRECISION  reference=%s  norm=%s ============\n", refm, nm))
    cat("\n-- Per cancer x stage --\n")
    s <- rp_stage[rp_stage$reference == refm & rp_stage$norm_method == nm, ]
    s <- s[order(s$cancer_id, s$stage, s$model, s$alpha, s$distance), ]
    print(s[, c("cancer_id", "stage", "model", "alpha", "distance",
                "n_tox", "n_ref", "shared", "recall", "precision", "f1")], row.names = FALSE)
    cat("\n-- Per cancer (summed over stages) --\n")
    c2 <- rp_cancer[rp_cancer$reference == refm & rp_cancer$norm_method == nm, ]
    c2 <- c2[order(c2$cancer_id, c2$model, c2$alpha, c2$distance), ]
    print(c2[, c("cancer_id", "model", "alpha", "distance",
                 "recall", "precision", "f1")], row.names = FALSE)
  }
}

# ---- precision-recall plots: one figure per (reference, norm), facet cancer x stage ----
cat("\n>>> writing precision-recall plots\n")
for (nm in norms) {
  for (refm in present_refs) {
    pd <- rp_stage[rp_stage$reference == refm & rp_stage$norm_method == nm, ]
    if (nrow(pd) == 0) next
    pd$alpha_f <- factor(pd$alpha); pd$dist_f <- factor(pd$distance)
    p <- ggplot(pd, aes(x = recall, y = precision, colour = model, shape = alpha_f)) +
      geom_path(aes(group = interaction(model, alpha_f)), linewidth = 0.3, alpha = 0.6) +
      geom_point(aes(size = dist_f)) +
      facet_grid(cancer_id ~ stage) +
      scale_size_discrete(name = "distance p<", range = c(1.2, 3)) +
      scale_shape_discrete(name = "alpha") +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = sprintf("Precision-Recall vs %s  (norm=%s)", refm, nm),
           subtitle = "path = distance sweep {0.05, 0.1, 0.25}; shape = alpha {0.01, 0.05}",
           x = "recall", y = "precision") +
      theme_minimal(base_size = 8) + theme(strip.text.y = element_text(angle = 0))
    tryCatch(ggsave(file.path(PLOT_DIR, sprintf("PR_%s_%s.png", refm, nm)),
                    p, width = 11, height = 16, dpi = 130),
             error = function(e) cat("    plot failed PR_", refm, "_", nm, "\n", sep = ""))
  }
}

# ---- one compact overview: recall & precision vs distance, faceted reference x norm ----
ov <- rp_overall
ov$alpha_f <- factor(ov$alpha)
p_rec <- ggplot(ov, aes(x = factor(distance), y = recall, colour = model,
                        linetype = alpha_f, group = interaction(model, alpha_f))) +
  geom_line() + geom_point(size = 2) + ylim(0, 1) +
  facet_grid(reference ~ norm_method) +
  labs(title = "OVERALL recall vs distance (alpha = solid 0.01 / dashed 0.05)",
       x = "distance p-value <", y = "recall", colour = "model", linetype = "alpha") +
  theme_minimal()
tryCatch(ggsave(file.path(PLOT_DIR, "overall_recall_vs_distance.png"), p_rec,
                width = 9, height = 11, dpi = 150), error = function(e) NULL)
p_prec <- p_rec %+% ov + aes(y = precision) +
  labs(title = "OVERALL precision vs distance (alpha = solid 0.01 / dashed 0.05)", y = "precision")
tryCatch(ggsave(file.path(PLOT_DIR, "overall_precision_vs_distance.png"), p_prec,
                width = 9, height = 11, dpi = 150), error = function(e) NULL)

cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("DONE. Detailed comparison written to:", OUT_DIR, "\n")
cat("  CSVs : pvalue_correlation_by_cancer_stage, recall_precision_{per_cancer_stage,per_cancer,overall}\n")
cat("  Plots:", PLOT_DIR, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
