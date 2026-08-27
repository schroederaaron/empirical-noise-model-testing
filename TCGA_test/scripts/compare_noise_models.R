#!/usr/bin/env Rscript
# compare_noise_models.R
# -----------------------------------------------------------------------------
# Runs the Tensor-Omics noise-model variants on the *same* preprocessed inputs and
# compares them:
#
#   1. bootstrap : tox_compute_noise_pvalues_pipeline        (between-neighbourhood
#                  mean-difference bootstrap null)
#   2. exact     : tox_compute_noise_pvalues_pipeline_exact  (exact sorted-pool
#                  pairwise count + 1/sqrt(n) scaling; the production pick)
#
# For each (project, norm, stage) the neighbourhood gathering and observed
# statistics are identical across the calls -- ONLY the null construction differs.
#
# A Tox gene is a "hit" iff BOTH noise_p_value_adj < alpha AND distance_p_value <
# dist_thresh. The whole comparison is swept over alpha (ALPHA_VALUES, the recall
# lever) x distance (DISTANCE_THRESHOLDS, the effect-size cut).
#
# Reports:
#   * significant-gene counts per model x alpha x distance (overall/per-cancer/finest)
#   * Jaccard between Tox models (at PRIMARY_ALPHA) + Spearman of raw p-values
# and -- if reference DE results from deg_comparison.R exist (edgeR / limma-voom /
# NOISeq [semi-parametric NOISeqBIO]) --
#   (A) reference/consensus agreement (edgeR/voom/noiseq + consensus2 [>=2 methods]
#       and consensus3 [all 3]; is the "ground truth" self-consistent?)
#   (B) each Tox model vs each reference AND vs the consensus sets -- recall /
#       precision, swept over alpha x distance, at overall / per-cancer /
#       per-cancer-stage grouping. (Recall-vs-consensus >> recall-vs-single tests
#       whether Tox recovers the reliable core; per-cancer/stage shows whether any
#       recall gap is concentrated rather than global.)
#
# Data loading / preprocessing is reused verbatim from
# outlier_significance_analysis.R.  That script's main block is guarded by
# `if (sys.nframe() == 0)`, so sourcing it only imports helper functions and the
# config globals (CANCER_TYPES, NORM_METHODS, STAGES, COMP_TYPES, ...) WITHOUT
# running the standalone analysis.
# -----------------------------------------------------------------------------

source("outlier_significance_analysis.R")

# ==================== CONFIGURATION ====================

# LFCseq dropped: ~identical results to bootstrap (Jaccard 0.90-0.94) but the
# slowest by far (~17000 s vs exact's ~640 s). Kept only bootstrap (liberal) and
# exact (conservative, sqrt(n), the production pick).
MODELS <- c("bootstrap", "exact")
MODEL_DISPLAY <- c(
  bootstrap = "Bootstrap (mean-diff)",
  exact     = "Exact (+1/sqrt n)"
)

# Scope of the comparison. The full grid (8 cancers x 4 norms x 4 stages x 3
# models) is large; restrict here for a tractable run. Set to the full globals
# (CANCER_TYPES / c("raw","log","std_log","full")) to compare everything.
COMPARE_CANCERS <- CANCER_TYPES                 # e.g. CANCER_TYPES["bladder cancer"]
COMPARE_NORMS   <- c("raw", "log")                 # e.g. c("raw","log","std_log","full")

# Neighbourhood parameters (match outlier_significance_analysis.R).
K_START <- 5L; K_STEP <- 1L; K_MAX <- 20L; TAU <- 0.15; MAX_POOL <- 50000L

# Significance thresholds applied to the BH-adjusted noise p-value. The whole
# comparison is swept over these so we can see how much loosening the noise cut
# (the actual recall lever) buys -- e.g. how many more outliers 0.05 admits over
# 0.001. PRIMARY_ALPHA is the single value used where a scalar is needed (the
# model-vs-model Jaccard, plots).
ALPHA_VALUES  <- c(0.001, 0.01, 0.05)
PRIMARY_ALPHA <- 0.01

# Effect-size (distance) thresholds for the Tox outliers -- the Tensor-Omics
# analogue of edgeR/limma's |logFC| > 2 filter. A gene counts as a Tox "hit" only
# if BOTH  noise_p_value_adj < ALPHA  AND  distance_p_value < <threshold>.
# Smaller distance p-value = larger effect, so these are the "upper X percentile"
# of effect size: 1%, 5%, 10%, 25%, 50%. The whole comparison is swept over them.
DISTANCE_THRESHOLDS <- c(0.01, 0.05, 0.10, 0.25, 0.50)

# Cores for parallel combination processing (mclapply). Combinations run in
# parallel; the three models within a combination stay sequential.
N_CORES <- 64L

# Reference DE methods produced by deg_comparison.R (one CSV per method per
# project/stage). Each defines its own significant-gene set (the "ground truth"
# each is compared against), spanning two method families:
#   edgeR / voom  -- count-parametric (NB GLM / voom-weighted linear model)
#   noiseq        -- semi-parametric (empirical-Bayes NOISeqBIO)
# (SAMseq, a fully non-parametric option, was dropped: it does not scale to TCGA's
#  sample sizes. A scalable non-parametric alternative -- the Wilcoxon rank-sum
#  test -- could be added later if a fully non-parametric reference is wanted.)
# All apply the same effect-size cut |logFC| > DEG_LOGFC; significance uses each
# method's native score (FDR / adj.P.Val / NOISeq prob).
REFERENCE_METHODS <- c("edgeR", "voom", "noiseq")
DEG_FDR         <- 0.01   # edgeR FDR / limma adj.P.Val cutoff
DEG_LOGFC       <- 2      # |logFC| effect-size cutoff (all reference methods)
DEG_NOISEQ_PROB <- 0.95   # NOISeq prob = P(DE) cutoff (its native DE call)

# Only the `own_healthy` comparison is active (family/ortholog are disabled in the
# same way as the reference script), so the model comparison is done on it.
COMPARE_COMP <- "own_healthy"

# Restrict the config globals so the (heavy) loaders only touch the chosen scope.
CANCER_TYPES <- COMPARE_CANCERS
NORM_METHODS <- COMPARE_NORMS

# ==================== MODEL DISPATCH ====================

#' Return the R wrapper function for a given model name.
noise_model_wrapper <- function(model) {
  switch(model,
    bootstrap = tox_compute_noise_pvalues_pipeline,
    exact     = tox_compute_noise_pvalues_pipeline_exact,
    stop(sprintf("Unknown model '%s'", model))
  )
}

#' Build the argument list shared by all three model wrappers for one combination.
#'
#' Mirrors the input marshalling in `compute_noise_pvalues()` (see
#' outlier_significance_analysis.R): pre-log (partially normalized) replicate
#' matrices, the unscaled observed `own` shift, and norm_method collapsed to
#' 0 (linear) / 1 (log2). Family/ortholog have been retired from the model, so
#' `family_stats` / `gene_to_fam` are no longer passed to the pipeline (they are
#' still used elsewhere in the caller for family-size annotation).
build_model_args <- function(cancer_preproc, healthy_preproc, gene_results,
                             family_stats, gene_to_fam, stage, norm_method) {

  n_genes   <- nrow(gene_results)
  stage_col <- gsub(" ", "_", stage)

  obs_own <- gene_results[[paste0("shift_vs_own_healthy_unscaled_", stage_col)]]
  valid_genes_own <- as.integer(!is.na(obs_own))
  obs_own[is.na(obs_own)] <- 0

  norm_method_int <- switch(norm_method,
    "raw" = 0L, "log" = 1L, "std_log" = 1L, "full" = 1L, 0L)

  list(
    case_means         = as.numeric(cancer_preproc$means),
    case_replicates    = cancer_preproc$prelog,
    control_means      = as.numeric(healthy_preproc$means),
    control_replicates = healthy_preproc$prelog,
    obs_own            = as.numeric(obs_own),
    valid_genes_own    = valid_genes_own,
    norm_method        = norm_method_int,
    k_start            = K_START,
    k_step             = K_STEP,
    k_max              = K_MAX,
    tau                = TAU,
    max_pool_size      = MAX_POOL
  )
}

#' Run one model on a pre-built argument list; time ONLY the tox call.
#' Returns the tox result and the elapsed seconds (the caller aggregates timing,
#' since a global accumulator would not survive `mclapply` forks).
run_one_model <- function(model, args, label = "") {
  fn <- noise_model_wrapper(model)
  t0 <- Sys.time()
  res <- do.call(fn, args)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("    [TIMING] %-9s %s : %.3f s\n", model, label, elapsed))
  list(result = res, time = elapsed)
}

# ==================== PER-COMBINATION PROCESSING ====================

#' Run all three models on one (project, norm, stage) combination and return a
#' long data frame with one row per (gene, model) for the active comparison.
process_combination_all_models <- function(cancer_type, project_id, norm_method, stage,
                                            gene_results, family_stats, healthy_preproc,
                                            gene_to_fam, n_families, norm_display,
                                            kept_gene_ids) {

  empty_times <- setNames(numeric(length(MODELS)), MODELS)

  # ---- load + filter + preprocess the cancer stage (as in the reference) ----
  cancer_data_raw <- load_stage_data(
    project_id = project_id, stage = stage, data_type = "cancer",
    use_constant_healthy = FALSE, norm_method = "raw",
    apply_mean = FALSE, normalize = FALSE
  )
  if (is.null(cancer_data_raw) || !is.matrix(cancer_data_raw$expression_vectors))
    return(list(rows = NULL, times = empty_times))

  cancer_raw_mat <- cancer_data_raw$expression_vectors
  if (is.null(colnames(cancer_raw_mat)) && !is.null(cancer_data_raw$gene_ids))
    colnames(cancer_raw_mat) <- cancer_data_raw$gene_ids
  if (length(setdiff(kept_gene_ids, colnames(cancer_raw_mat))) > 0) {
    cat(sprintf("    Warning: kept genes absent from cancer stage %s; skipping\n", stage))
    return(list(rows = NULL, times = empty_times))
  }
  cancer_raw_mat  <- cancer_raw_mat[, kept_gene_ids, drop = FALSE]
  cancer_preproc  <- preprocess_replicates(cancer_raw_mat, norm_method)

  # ---- shared inputs; identical neighbourhoods across models ----
  args <- build_model_args(cancer_preproc, healthy_preproc, gene_results,
                           family_stats, gene_to_fam, stage, norm_method)

  label <- sprintf("%s/%s/%s", project_id, stage, norm_method)

  # Shared per-gene quantities (distance p, shift, family info, neighbourhood).
  distance_p     <- extract_distance_pvalues(gene_results)
  family_n_genes <- tabulate(gene_to_fam, nbins = n_families)
  n_genes        <- nrow(gene_results)
  gene_ids       <- gene_results$gene_id

  rows <- list(); ri <- 1L
  times <- empty_times

  for (model in MODELS) {
    run    <- run_one_model(model, args, label)
    res    <- run$result
    times[[model]] <- run$time

    p_own  <- res$pvalues_own
    p_own[p_own < 0 | p_own > 1] <- NA           # -1 => not computed
    n_cancer  <- res$neighborhood_size_case
    n_own     <- res$neighborhood_size_own_control
    n_cancer[n_cancer < 0] <- NA
    n_own[n_own < 0]       <- NA

    for (g in seq_len(n_genes)) {
      dist_p  <- distance_p[g, stage, COMPARE_COMP]
      noise_p <- p_own[g]
      if (is.na(dist_p) || is.na(noise_p)) next

      shift_col <- paste0("shift_vs_own_healthy_", gsub(" ", "_", stage))
      if (!shift_col %in% colnames(gene_results)) next
      shift <- gene_results[g, shift_col]

      fam_id      <- gene_to_fam[g]
      family_size <- if (!is.na(fam_id) && fam_id > 0) family_n_genes[fam_id] else NA

      rows[[ri]] <- data.frame(
        gene_id                  = gene_ids[g],
        cancer_type              = cancer_type,
        cancer_id                = project_id,
        stage                    = stage,
        normalization            = norm_display,
        norm_method              = norm_method,
        comparison               = COMPARE_COMP,
        model                    = model,
        direction                = ifelse(is.na(shift), NA, ifelse(shift > 0, "up", "down")),
        signed_difference        = shift,
        distance_p_value         = dist_p,
        noise_p_value            = noise_p,
        family_size              = family_size,
        neighborhood_size_cancer = n_cancer[g],
        neighborhood_size_healthy= n_own[g],
        stringsAsFactors = FALSE
      )
      ri <- ri + 1L
    }
  }

  df <- if (length(rows) == 0) NULL else do.call(rbind, rows)
  list(rows = df, times = times)
}

# ==================== MODEL-COMPARISON METRICS ====================
#
# A Tox gene is a "hit" only if it clears BOTH filters:
#   * statistical significance : noise_p_value_adj < alpha
#   * effect size              : distance_p_value  < dist_thresh
# The effect-size (distance) filter is the Tensor-Omics analogue of edgeR/limma's
# |logFC| > 2. Everything below is swept over `DISTANCE_THRESHOLDS`, and reported
# at several grouping levels (overall / per-cancer / per-cancer x stage x norm).

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

#' Composite gene-instance key (unique per gene x cancer x stage x norm).
tox_key <- function(d) paste(d$gene_id, d$cancer_id, d$stage, d$norm_method, sep = "|")

#' The testable base rows (own comparison, finite noise + distance p-values).
tox_base <- function(tbl) {
  tbl %>% filter(comparison == COMPARE_COMP,
                 !is.na(noise_p_value_adj), !is.na(distance_p_value))
}

#' Per-group significant COUNTS for every model, swept over alpha x distance.
#' `group_cols` may be character(0) for the overall summary.
significant_counts <- function(tbl, alphas, dist_thresholds, group_cols) {
  base <- tox_base(tbl)
  grid <- expand.grid(alpha = alphas, t = dist_thresholds)
  purrr::map_dfr(seq_len(nrow(grid)), function(k) {
    a <- grid$alpha[k]; t <- grid$t[k]
    base %>%
      group_by(across(all_of(c(group_cols, "model")))) %>%
      summarise(n_tested = dplyr::n(),
                n_signif = sum(noise_p_value_adj < a & distance_p_value < t),
                .groups = "drop") %>%
      mutate(alpha = a, distance_threshold = t,
             pct_signif = round(100 * n_signif / pmax(n_tested, 1), 2))
  })
}

#' Per-group, per-threshold pairwise model Jaccard + overlap (tidy long df).
pairwise_jaccard <- function(tbl, alpha, dist_thresholds, group_cols) {
  base  <- tox_base(tbl) %>% mutate(.key = paste(gene_id, cancer_id, stage, norm_method, sep = "|"))
  pairs <- utils::combn(MODELS, 2, simplify = FALSE)

  one_group <- function(g_df, g_label) {
    purrr::map_dfr(dist_thresholds, function(t) {
      sg   <- g_df[g_df$noise_p_value_adj < alpha & g_df$distance_p_value < t, ]
      sets <- lapply(MODELS, function(m) unique(sg$.key[sg$model == m])); names(sets) <- MODELS
      purrr::map_dfr(pairs, function(pr) {
        a <- sets[[pr[1]]]; b <- sets[[pr[2]]]
        cbind(g_label, data.frame(
          distance_threshold = t, model_A = pr[1], model_B = pr[2],
          n_A = length(a), n_B = length(b),
          shared = length(intersect(a, b)),
          only_A = length(setdiff(a, b)), only_B = length(setdiff(b, a)),
          jaccard = round(jaccard(a, b), 4), stringsAsFactors = FALSE))
      })
    })
  }

  if (length(group_cols) == 0) {
    one_group(as.data.frame(base), data.frame(scope = "ALL", stringsAsFactors = FALSE))
  } else {
    g    <- base %>% group_by(across(all_of(group_cols)))
    keys <- g %>% group_keys()
    sub  <- g %>% group_split()
    purrr::map_dfr(seq_along(sub),
                   function(i) one_group(as.data.frame(sub[[i]]), keys[i, , drop = FALSE]))
  }
}

#' Full model-vs-model comparison, swept over alpha x distance, at several
#' grouping levels. Counts sweep the full alpha grid (the recall lever); the
#' pairwise Jaccard is reported at PRIMARY_ALPHA only, to bound the output.
compare_models <- function(tbl, alphas, dist_thresholds, output_dir) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("MODEL COMPARISON  (noise p_adj < {%s}  AND  distance p < {%s},  comparison = %s)\n",
              paste(alphas, collapse = ", "), paste(dist_thresholds, collapse = ", "), COMPARE_COMP))
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

  # ---- significant counts (overall / per-cancer / finest), swept over alpha ----
  counts_overall <- significant_counts(tbl, alphas, dist_thresholds, character(0))
  counts_cancer  <- significant_counts(tbl, alphas, dist_thresholds, "cancer_id")
  counts_detail  <- significant_counts(tbl, alphas, dist_thresholds,
                                       c("cancer_id", "stage", "norm_method"))

  cat("\n-- Significant Tox genes per model x alpha x distance threshold (OVERALL) --\n")
  overall_wide <- counts_overall %>%
    select(model, alpha, distance_threshold, n_signif) %>%
    tidyr::pivot_wider(names_from = distance_threshold, values_from = n_signif,
                       names_prefix = "d<") %>%
    arrange(model, alpha)
  print(as.data.frame(overall_wide), row.names = FALSE)

  # ---- pairwise Jaccard at PRIMARY_ALPHA (overall / per-cancer / finest) ----
  jac_overall <- pairwise_jaccard(tbl, PRIMARY_ALPHA, dist_thresholds, character(0))
  jac_cancer  <- pairwise_jaccard(tbl, PRIMARY_ALPHA, dist_thresholds, "cancer_id")
  jac_detail  <- pairwise_jaccard(tbl, PRIMARY_ALPHA, dist_thresholds,
                                  c("cancer_id", "stage", "norm_method"))

  cat(sprintf("\n-- Pairwise model Jaccard by distance threshold (OVERALL, alpha=%.3g) --\n", PRIMARY_ALPHA))
  print(jac_overall %>%
          select(distance_threshold, model_A, model_B, n_A, n_B, shared, jaccard) %>%
          as.data.frame(), row.names = FALSE)

  cat(sprintf("\n-- Per-cancer significant counts (n_signif per model, alpha=%.3g) --\n", PRIMARY_ALPHA))
  print(counts_cancer %>%
          filter(alpha == PRIMARY_ALPHA) %>%
          select(cancer_id, model, distance_threshold, n_signif) %>%
          tidyr::pivot_wider(names_from = model, values_from = n_signif) %>%
          arrange(cancer_id, distance_threshold) %>%
          as.data.frame(), row.names = FALSE)

  # ---- Spearman correlation of raw noise p-values (threshold-independent) ----
  wide <- tox_base(tbl) %>%
    mutate(key = paste(gene_id, cancer_id, stage, norm_method, sep = "|")) %>%
    select(key, model, noise_p_value) %>%
    tidyr::pivot_wider(names_from = model, values_from = noise_p_value)
  present <- intersect(MODELS, colnames(wide))
  corr <- NULL
  if (length(present) >= 2) {
    corr <- cor(as.matrix(wide[, present, drop = FALSE]),
                use = "pairwise.complete.obs", method = "spearman")
    cat("\n-- Spearman correlation of raw noise p-values --\n")
    print(round(corr, 3))
  }

  # ---- save the detailed tables (the exhaustive per-cancer/stage/norm view) ----
  write.csv(counts_overall, file.path(output_dir, "counts_overall_by_alpha_threshold.csv"), row.names = FALSE)
  write.csv(counts_cancer,  file.path(output_dir, "counts_per_cancer_by_alpha_threshold.csv"), row.names = FALSE)
  write.csv(counts_detail,  file.path(output_dir, "counts_per_cancer_stage_norm_by_alpha_threshold.csv"), row.names = FALSE)
  write.csv(jac_overall,    file.path(output_dir, "jaccard_overall_by_threshold.csv"), row.names = FALSE)
  write.csv(jac_cancer,     file.path(output_dir, "jaccard_per_cancer_by_threshold.csv"), row.names = FALSE)
  write.csv(jac_detail,     file.path(output_dir, "jaccard_per_cancer_stage_norm_by_threshold.csv"), row.names = FALSE)

  # ---- plots ----
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # (1) significant counts vs distance threshold, per model, faceted by alpha
  p_cnt <- ggplot(counts_overall,
                  aes(x = factor(distance_threshold), y = n_signif, colour = model, group = model)) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ alpha, labeller = label_both) +
    labs(title = "Significant Tox genes vs effect-size threshold",
         subtitle = "faceted by noise p_adj cutoff (alpha)",
         x = "distance p-value <", y = "n significant", colour = "model") +
    theme_minimal()
  ggsave(file.path(plots_dir, "significant_counts_vs_threshold.png"), p_cnt, width = 10, height = 5, dpi = 300)

  # (2) significant counts per cancer at PRIMARY_ALPHA, faceted
  p_cancer <- ggplot(counts_cancer %>% filter(alpha == PRIMARY_ALPHA),
                     aes(x = factor(distance_threshold), y = n_signif, fill = model)) +
    geom_col(position = "dodge") +
    facet_wrap(~ cancer_id, scales = "free_y") +
    labs(title = sprintf("Significant Tox genes per cancer x threshold (alpha=%.3g)", PRIMARY_ALPHA),
         x = "distance p-value <", y = "n significant") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plots_dir, "significant_counts_per_cancer.png"), p_cancer, width = 11, height = 7, dpi = 300)

  invisible(list(
    counts_overall = counts_overall, counts_cancer = counts_cancer, counts_detail = counts_detail,
    jaccard_overall = jac_overall, jaccard_cancer = jac_cancer, jaccard_detail = jac_detail,
    spearman = corr))
}

#' Symmetric Jaccard matrix over a named list of key sets.
jaccard_matrix <- function(sets) {
  nm <- names(sets)
  M  <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  for (i in seq_along(nm)) for (j in seq_along(nm)) M[i, j] <- jaccard(sets[[i]], sets[[j]])
  M
}

#' Load a reference method's significant gene_ids for one (project, stage), using
#' that method's NATIVE significance score plus the shared |logFC| > DEG_LOGFC cut.
#' edgeR -> FDR, limma-voom -> adj.P.Val, NOISeq -> prob (P(DE)).
load_ref_hits <- function(de_dir, project_id, stage, method) {
  fn <- file.path(de_dir, project_id, stage, paste0(method, "_results.csv"))
  cat(paste("Loading reference hits from: \n", fn, "\n"))
  if (!file.exists(fn)) return(character(0))
  d <- read.csv(fn)
  if (!all(c("gene_id", "logFC") %in% colnames(d))) return(character(0))
  keep <-
    if (method == "edgeR") {
      if (!"FDR" %in% colnames(d)) return(character(0)); d$FDR < DEG_FDR
    } else if (method == "voom") {
      if (!"adj.P.Val" %in% colnames(d)) return(character(0)); d$adj.P.Val < DEG_FDR
    } else {  # noiseq: prob is P(differentially expressed)
      if (!"prob" %in% colnames(d)) return(character(0)); d$prob > DEG_NOISEQ_PROB
    }
  keep <- keep & (abs(d$logFC) > DEG_LOGFC)
  hits <- d$gene_id[which(keep)]
  cat(paste("Keeping ", length(hits), " genes \n"))
  hits                       # <- MUST be the last expression (the function's return value)
}

#' Compare Tox models against the reference DE methods AND their consensus sets,
#' swept over alpha x distance, at three grouping levels (overall / per-cancer /
#' per-cancer-stage). Adds two things over the old version:
#'   * CONSENSUS references -- genes called by >=2 (consensus2) and all 3 (consensus3)
#'     of edgeR/voom/noiseq. Recall-vs-consensus >> recall-vs-single would confirm
#'     that Tox recovers the reliable core and the single-method numbers are pessimistic.
#'   * PER-CANCER / PER-STAGE recall, to see whether the recall gap is concentrated
#'     in specific cancers/stages rather than global.
#' Reference hit sets are per (cancer, stage) gene sets; Tox hits drop the norm
#' dimension to match. recall/precision/jaccard aggregate correctly across combos
#' because each combo's keys are disjoint (cancer/stage are fixed within a combo).
compare_with_references <- function(tbl, alphas, dist_thresholds, output_dir) {
  if (!exists("TOX_TEST_DIR") || !dir.exists(TOX_TEST_DIR)) {
    cat("\n(No DEG output dir found; skipping reference comparison.)\n")
    return(invisible(NULL))
  }
  de_dir <- file.path(dirname(TOX_TEST_DIR), "differential_expression")

  base   <- as.data.frame(tox_base(tbl))
  combos <- unique(base[, c("cancer_id", "stage")])
  norms  <- sort(unique(base$norm_method))
  REF_ALL <- c(REFERENCE_METHODS, "consensus2", "consensus3")

  # ---- per-combo reference + consensus gene sets ----
  ref_sets <- vector("list", nrow(combos))
  for (i in seq_len(nrow(combos))) {
    pid <- combos$cancer_id[i]; stg <- combos$stage[i]
    pm  <- lapply(REFERENCE_METHODS, function(m) unique(load_ref_hits(de_dir, pid, stg, m)))
    names(pm) <- REFERENCE_METHODS
    cnt <- table(unlist(pm, use.names = FALSE))          # #methods calling each gene
    pm$consensus2 <- names(cnt)[cnt >= 2]
    pm$consensus3 <- names(cnt)[cnt >= 3]
    ref_sets[[i]] <- pm
  }
  present_refs <- REF_ALL[vapply(REF_ALL, function(m)
      any(vapply(ref_sets, function(s) length(s[[m]]) > 0, logical(1))), logical(1))]
  if (!any(REFERENCE_METHODS %in% present_refs)) {
    cat("\n(No usable DEG result files found; skipping reference comparison.)\n")
    return(invisible(NULL))
  }
  plots_dir <- file.path(output_dir, "plots"); dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  # ===== (A) reference / consensus agreement (aggregated over combos) =====
  agg_keys <- function(m) {
    ks <- character(0)
    for (i in seq_len(nrow(combos)))
      if (length(ref_sets[[i]][[m]]) > 0)
        ks <- c(ks, paste(ref_sets[[i]][[m]], combos$cancer_id[i], combos$stage[i], sep = "|"))
    unique(ks)
  }
  ref_keys <- setNames(lapply(present_refs, agg_keys), present_refs)
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("REFERENCE vs REFERENCE (+ consensus)  (is the ground truth self-consistent?)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("\nHits per set: %s\n",
      paste(sprintf("%s=%d", present_refs, vapply(ref_keys, length, 0L)), collapse = "  ")))
  Rr <- jaccard_matrix(ref_keys)
  cat("\n-- Jaccard among reference / consensus sets --\n"); print(round(Rr, 3))
  refref_long <- as.data.frame(as.table(Rr), stringsAsFactors = FALSE)
  colnames(refref_long) <- c("set_A", "set_B", "jaccard")
  write.csv(refref_long, file.path(output_dir, "reference_vs_reference_jaccard.csv"), row.names = FALSE)

  # ===== (B) per-combo building blocks: model vs each reference/consensus =====
  n_total <- nrow(combos) * length(norms) * length(MODELS) * length(alphas) *
             length(dist_thresholds) * length(present_refs)
  v_cancer <- character(n_total); v_stage <- character(n_total); v_norm <- character(n_total)
  v_model  <- character(n_total); v_alpha <- numeric(n_total); v_dist  <- numeric(n_total)
  v_ref    <- character(n_total); v_nq <- integer(n_total); v_nr <- integer(n_total); v_sh <- integer(n_total)
  z <- 0L
  for (i in seq_len(nrow(combos))) {
    pid <- combos$cancer_id[i]; stg <- combos$stage[i]
    cb  <- base[base$cancer_id == pid & base$stage == stg, ]
    for (nm in norms) {
      cbn <- cb[cb$norm_method == nm, ]
      for (model in MODELS) {
        cm <- cbn[cbn$model == model, ]
        for (a in alphas) for (dt in dist_thresholds) {
          q  <- unique(cm$gene_id[cm$noise_p_value_adj < a & cm$distance_p_value < dt])
          nq <- length(q)
          for (refm in present_refs) {
            rk <- ref_sets[[i]][[refm]]
            z <- z + 1L
            v_cancer[z] <- pid; v_stage[z] <- stg; v_norm[z] <- nm; v_model[z] <- model
            v_alpha[z]  <- a;   v_dist[z] <- dt;   v_ref[z] <- refm
            v_nq[z] <- nq; v_nr[z] <- length(rk); v_sh[z] <- length(intersect(q, rk))
          }
        }
      }
    }
  }
  detail <- data.frame(cancer_id = v_cancer, stage = v_stage, norm_method = v_norm, model = v_model,
                       alpha = v_alpha, distance_threshold = v_dist, reference = v_ref,
                       n_query = v_nq, n_ref = v_nr, shared = v_sh, stringsAsFactors = FALSE)

  # recall/precision/jaccard aggregate additively across combos (disjoint keys).
  agg_metrics <- function(d, group_cols) {
    gv <- c(group_cols, "norm_method", "model", "alpha", "distance_threshold", "reference")
    a  <- aggregate(d[c("n_query", "n_ref", "shared")], by = d[gv], FUN = sum)
    a$recall    <- round(a$shared / pmax(a$n_ref, 1), 4)     # ref hits recovered
    a$precision <- round(a$shared / pmax(a$n_query, 1), 4)   # model hits in ref
    a$jaccard   <- round(a$shared / pmax(a$n_query + a$n_ref - a$shared, 1), 4)
    a
  }
  mvr_overall <- agg_metrics(detail, character(0))
  mvr_cancer  <- agg_metrics(detail, "cancer_id")
  mvr_detail  <- agg_metrics(detail, c("cancer_id", "stage"))
  write.csv(mvr_overall, file.path(output_dir, "model_vs_reference_overall.csv"), row.names = FALSE)
  write.csv(mvr_cancer,  file.path(output_dir, "model_vs_reference_per_cancer.csv"), row.names = FALSE)
  write.csv(mvr_detail,  file.path(output_dir, "model_vs_reference_per_cancer_stage.csv"), row.names = FALSE)

  rep_norm <- if ("log" %in% norms) "log" else if ("std_log" %in% norms) "std_log" else norms[1]
  mid_dist <- dist_thresholds[ceiling(length(dist_thresholds) / 2)]

  # --- the consensus hypothesis test: does recall jump on the consensus sets? ---
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("TOX vs {edgeR, consensus2, consensus3}  (norm=%s)  -- reliability test\n", rep_norm))
  cat("If recall jumps sharply from edgeR -> consensus3, Tox is finding the reliable\n")
  cat("core and single-method recall understates it.\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  sub <- mvr_overall[mvr_overall$norm_method == rep_norm &
                     mvr_overall$reference %in% c("edgeR", "consensus2", "consensus3"), ]
  sub <- sub[order(sub$model, sub$alpha, sub$distance_threshold), ]
  print(sub[, c("model", "alpha", "distance_threshold", "reference", "n_ref", "recall", "precision")],
        row.names = FALSE)

  # --- alpha effect on Tox hit counts / recall (the recall lever) ---
  cat(sprintf("\n-- ALPHA effect (norm=%s, vs edgeR): how many more hits does looser alpha admit? --\n", rep_norm))
  sub2 <- mvr_overall[mvr_overall$norm_method == rep_norm & mvr_overall$reference == "edgeR", ]
  sub2 <- sub2[order(sub2$model, sub2$distance_threshold, sub2$alpha), ]
  print(sub2[, c("model", "distance_threshold", "alpha", "n_query", "recall", "precision")], row.names = FALSE)

  # --- per-cancer recall (is the gap concentrated?) at PRIMARY_ALPHA, mid distance ---
  cat(sprintf("\n-- Per-cancer recall vs edgeR (norm=%s, alpha=%.3g, distance<%.2g) --\n",
              rep_norm, PRIMARY_ALPHA, mid_dist))
  pc <- mvr_cancer[mvr_cancer$norm_method == rep_norm & mvr_cancer$reference == "edgeR" &
                   mvr_cancer$alpha == PRIMARY_ALPHA & mvr_cancer$distance_threshold == mid_dist, ]
  print(pc[order(pc$cancer_id, pc$model), c("cancer_id", "model", "n_query", "n_ref", "recall", "precision")],
        row.names = FALSE)

  # --- plot: recall vs distance, references + consensus, at PRIMARY_ALPHA ---
  pdat <- mvr_overall[mvr_overall$alpha == PRIMARY_ALPHA, ]
  p_rec <- ggplot(pdat,
                  aes(x = factor(distance_threshold), y = recall, colour = model, group = model)) +
    geom_line() + geom_point(size = 2) + ylim(0, 1) +
    facet_grid(reference ~ norm_method) +
    labs(title = sprintf("Recall of reference/consensus hits vs distance (alpha=%.3g)", PRIMARY_ALPHA),
         x = "distance p-value <", y = "recall", colour = "model") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plots_dir, "recall_vs_reference_by_threshold.png"), p_rec, width = 11, height = 8, dpi = 300)

  invisible(list(reference_jaccard = Rr, model_vs_reference = mvr_overall,
                 model_vs_reference_cancer = mvr_cancer, model_vs_reference_detail = mvr_detail,
                 present_refs = present_refs))
}

# ==================== MAIN DRIVER ====================

compare_all_models <- function(output_dir = NULL, n_cores = 64L) {

  cat("\n", paste(rep("#", 80), collapse = ""), "\n", sep = "")
  cat("### NOISE-MODEL COMPARISON (bootstrap / exact) ###\n")
  cat(sprintf("Scope: cancers=%d  norms={%s}  stages={%s}\n",
              length(CANCER_TYPES), paste(NORM_METHODS, collapse = ","),
              paste(STAGES, collapse = ",")))
  cat(sprintf("Neighbourhood (genes): k_start=%d k_step=%d k_max=%d tau=%g\n",
              K_START, K_STEP, K_MAX, TAU))
  cat(sprintf("Parallel combinations: mc.cores=%d\n", n_cores))
  cat(paste(rep("#", 80), collapse = ""), "\n", sep = "")

  if (is.null(output_dir)) output_dir <- file.path(dirname(TOX_TEST_DIR), "model_comparison")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  all_results <- load_all_gene_results()
  if (length(all_results) == 0) stop("No data found for the selected scope!")

  # One task per (project x norm x stage); each task runs all THREE models on the
  # same preprocessed inputs. Combinations are independent, so they parallelise
  # cleanly over cores (same granularity as create_all_genes_table in the
  # reference script). The three tox calls WITHIN a task stay sequential so the
  # shared input build happens once per combination.
  tasks <- list()
  for (key in names(all_results)) {
    dat <- all_results[[key]]
    for (stage in STAGES) {
      tasks[[length(tasks) + 1L]] <- list(
        cancer_type = dat$cancer_type, project_id = dat$cancer_id,
        norm_method = dat$norm_method, stage = stage,
        gene_results = dat$gene_results, family_stats = dat$family_stats,
        healthy_preproc = dat$healthy_preproc, gene_to_fam = dat$gene_to_fam,
        n_families = dat$n_families, norm_display = dat$norm_display,
        kept_gene_ids = dat$kept_gene_ids)
    }
  }
  cat(sprintf("\n>>> %d combinations x %d models over %d cores\n",
              length(tasks), length(MODELS), n_cores))

  task_results <- mclapply(tasks, function(task) {
    tryCatch(
      process_combination_all_models(
        cancer_type = task$cancer_type, project_id = task$project_id,
        norm_method = task$norm_method, stage = task$stage,
        gene_results = task$gene_results, family_stats = task$family_stats,
        healthy_preproc = task$healthy_preproc, gene_to_fam = task$gene_to_fam,
        n_families = task$n_families, norm_display = task$norm_display,
        kept_gene_ids = task$kept_gene_ids),
      error = function(e) {
        cat(sprintf("    ERROR %s/%s/%s: %s\n",
                    task$project_id, task$stage, task$norm_method, conditionMessage(e)))
        NULL
      })
  }, mc.cores = n_cores, mc.preschedule = FALSE)

  # Drop errored tasks (NULL); the rest are list(rows, times).
  ok <- Filter(is.list, task_results)
  cat(sprintf("    Completed: %d / %d combinations\n", length(ok), length(tasks)))

  row_dfs      <- Filter(Negate(is.null), lapply(ok, `[[`, "rows"))
  total_times  <- Reduce(`+`, lapply(ok, `[[`, "times"),
                         setNames(numeric(length(MODELS)), MODELS))

  if (length(row_dfs) == 0) stop("No rows produced.")
  tbl <- as.data.frame(rbindlist(row_dfs, fill = TRUE), stringsAsFactors = FALSE)

  # BH-adjust the noise p-value WITHIN each (model, cancer, norm, stage, comparison)
  # group -- exactly as the reference script does per model.
  tbl <- tbl %>%
    group_by(model, cancer_id, normalization, stage, comparison) %>%
    mutate(noise_p_value_adj = p.adjust(noise_p_value, method = "BH")) %>%
    ungroup()

  saveRDS(tbl, file.path(output_dir, "all_models_pvalues.rds"))
  write.csv(tbl, file.path(output_dir, "all_models_pvalues.csv"), row.names = FALSE)

  # ---- comparison (swept over alpha x distance; detailed groupings) ----
  cmp     <- compare_models(tbl, ALPHA_VALUES, DISTANCE_THRESHOLDS, output_dir)
  ref_cmp <- compare_with_references(tbl, ALPHA_VALUES, DISTANCE_THRESHOLDS, output_dir)

  # ---- timing summary (summed pure-tox time across all combinations) --
  cat("\n-- Total pure-tox runtime per model (summed over combinations) --\n")
  for (m in MODELS) cat(sprintf("   %-9s : %.3f s\n", m, total_times[[m]]))

  saveRDS(list(counts_overall = cmp$counts_overall, counts_cancer = cmp$counts_cancer,
               counts_detail = cmp$counts_detail,
               jaccard_overall = cmp$jaccard_overall, jaccard_cancer = cmp$jaccard_cancer,
               jaccard_detail = cmp$jaccard_detail,
               spearman = cmp$spearman, references = ref_cmp,
               alphas = ALPHA_VALUES, primary_alpha = PRIMARY_ALPHA,
               distance_thresholds = DISTANCE_THRESHOLDS, time = total_times),
          file.path(output_dir, "model_comparison_summary.rds"))

  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("DONE. Results in:", output_dir, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")

  invisible(list(table = tbl, comparison = cmp, references = ref_cmp))
}

# ==================== EXECUTION ====================

if (sys.nframe() == 0) {
  set.seed(42)
  out_dir <- file.path(dirname(TOX_TEST_DIR), "model_comparison")
  results <- compare_all_models(output_dir = out_dir, n_cores = N_CORES)
}
