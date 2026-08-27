#!/usr/bin/env Rscript
# outlier_significance_analysis.R
# Outlier Significance Analysis mit kNN-basierten Noise-p-Werten
# STAGE‑WEISE IMPLEMENTIERUNG mit unabhängigen Nachbarschaften für Tumor und Gesunde

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(purrr)
# library(patchwork)
# library(data.table)
# library(parallel)

source("config.R")
source("utils.R")
source("rcpp/tensoromics_functions.R")

# Load Fortran shared library
# dyn.load("noise_model.so")

# ==================== KONFIGURATION ====================

CANCER_TYPES <- c(
  "breast cancer" = "TCGA-BRCA",
  "non small cell lung cancer" = "TCGA-LUAD",
  "small cell lung cancer" = "TCGA-LUSC",
  "bladder cancer" = "TCGA-BLCA",
  "colon cancer" = "TCGA-COAD",
  "kidney cancer" = "TCGA-KIRC",
  "stomach cancer" = "TCGA-STAD",
  "thyroid gland cancer" = "TCGA-THCA"
)

NORM_METHODS <- c("raw", "std_log", "full", "log")
NORM_DISPLAY <- c(
  "raw" = "mean(TPM)",
  "log" = "log(mean(TPM))",
  "std_log" = "log(mean(gene-wise scaled TPM))",
  "full" = "log(mean(quantile(gene-wise scaled TPM)))"
)

# Alle drei Vergleichstypen
COMP_TYPES <- c("own_healthy", "family_mean", "ortholog_mean")
COMP_DISPLAY <- c(
  "own_healthy" = "Gene vs own healthy",
  "family_mean" = "Gene vs family mean",
  "ortholog_mean" = "Gene vs ortholog mean"
)

MIN_FAMILY_SIZES <- c(2, 3, 4)

# Für Plot 4: Speichere Perzentile für globale Visualisierung (optional)
STORE_NULL_PERCENTILES <- TRUE

# ==================== DATEN LADEN ====================

load_all_gene_results <- function() {
  all_results <- list()
  gene_keep_masks <- list()   # per-project cross-stage expression filter (cached)

  for (cancer_name in names(CANCER_TYPES)) {
    project_id <- CANCER_TYPES[cancer_name]
    cat(sprintf("\nLade %s (%s)...\n", cancer_name, project_id))

    for (norm_method in NORM_METHODS) {
      norm_suffix <- get_norm_suffix(norm_method)
      
      # Lade Gene Results (normalisiert)
      gene_file <- file.path(
        get_norm_output_dir(TOX_TEST_GENE_OUTPUT_DIR, norm_method),
        project_id,
        paste0(project_id, "_gene_results", norm_suffix, ".rds")
      )
      
      if (!file.exists(gene_file)) {
        cat(sprintf("  Datei nicht gefunden: %s\n", gene_file))
        next
      }
      
      gene_results <- readRDS(gene_file)
      
      # Lade Family LOESS Stats
      family_file <- file.path(
        get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method),
        project_id,
        paste0(project_id, "_family_loess_stats", norm_suffix, ".rds")
      )
      
      if (!file.exists(family_file)) {
        cat(sprintf("  Family Stats nicht gefunden: %s\n", family_file))
        next
      }
      
      family_stats <- readRDS(family_file)
      
      # Lade gesunde Rohdaten (für Noise) – einmal pro Normmethode
      healthy_data_raw <- load_stage_data(
        project_id, 
        STAGES[1], 
        "healthy", 
        use_constant_healthy = TRUE, 
        norm_method = "raw", 
        apply_mean = FALSE,
        normalize = FALSE   # Rohdaten ohne Normalisierung
      )
      
      if (is.null(healthy_data_raw) || !is.matrix(healthy_data_raw$expression_vectors)) {
        cat(sprintf("  Healthy RAW Data nicht gefunden für %s\n", norm_method))
        next
      }

      # ---- Cross-stage expression filter (computed once per project) ----
      if (is.null(gene_keep_masks[[project_id]])) {
        gene_keep_masks[[project_id]] <- compute_gene_keep_mask(project_id)
      }
      keep_mask    <- gene_keep_masks[[project_id]]
      kept_ids_all <- names(keep_mask)[keep_mask]

      # Restrict this norm's gene_results (and everything indexed alongside it) to
      # the kept genes, preserving gene_results row order.
      gene_ids_gr <- gene_results$gene_id
      keep_idx    <- which(gene_ids_gr %in% kept_ids_all)
      if (length(keep_idx) == 0) {
        cat(sprintf("  Keine Gene nach Expression-Filter übrig für %s / %s\n",
                    project_id, norm_method))
        next
      }
      gene_results  <- gene_results[keep_idx, , drop = FALSE]
      gene_to_fam   <- family_stats$gene_to_fam[keep_idx]
      kept_gene_ids <- gene_ids_gr[keep_idx]

      # Filter the healthy raw columns to the kept genes (in gene_results order)
      # BEFORE normalisation, so the gene-wise std-dev / quantile steps run on the
      # filtered gene set.
      healthy_raw_mat <- healthy_data_raw$expression_vectors
      if (is.null(colnames(healthy_raw_mat)) && !is.null(healthy_data_raw$gene_ids))
        colnames(healthy_raw_mat) <- healthy_data_raw$gene_ids
      healthy_raw_mat <- healthy_raw_mat[, kept_gene_ids, drop = FALSE]

      # Preprocess healthy replicates (einmal pro Normmethode)
      cat("    Preprocessing healthy replicates for norm:", norm_method, "\n")
      healthy_preproc <- preprocess_replicates(healthy_raw_mat, norm_method)

      all_results[[paste(project_id, norm_method, sep = "_")]] <- list(
        gene_results = gene_results,
        family_stats = family_stats,
        healthy_replicates_raw = healthy_raw_mat,
        healthy_preproc = healthy_preproc,
        gene_to_fam = gene_to_fam,
        n_families = family_stats$n_families,
        kept_gene_ids = kept_gene_ids,
        cancer_type = cancer_name,
        cancer_id = project_id,
        norm_method = norm_method,
        norm_display = NORM_DISPLAY[norm_method]
      )
    }
  }
  
  cat(sprintf("\n Geladen: %d Datensätze\n", length(all_results)))
  return(all_results)
}

extract_distance_pvalues <- function(gene_results) {
  
  n_genes <- nrow(gene_results)
  n_stages <- length(STAGES)
  n_comparisons <- length(COMP_TYPES)
  
  distance_pvalues <- array(NA, dim = c(n_genes, n_stages, n_comparisons),
                           dimnames = list(
                             gene_results$gene_id,
                             STAGES,
                             COMP_TYPES
                           ))
  
  for (s in 1:n_stages) {
    stage <- STAGES[s]
    stage_col <- gsub(" ", "_", stage)
    
    for (comp_type in COMP_TYPES) {
      p_col <- switch(comp_type,
        "own_healthy" = paste0("p_value_distance_own_", stage_col),
        "family_mean" = paste0("p_value_distance_fam_", stage_col),
        "ortholog_mean" = paste0("p_value_distance_orth_", stage_col)
      )
      
      if (p_col %in% colnames(gene_results)) {
        distance_pvalues[, s, comp_type] <- gene_results[[p_col]]
      }
    }
  }
  
  return(distance_pvalues)
}

# ==================== NEUE HELFER FÜR NOISE MODELL ====================

#' Preprocess replicate matrix into pre‑log space, means, residuals, and sorted order
preprocess_replicates <- function(expr_matrix, norm_method) {
  # expr_matrix: samples × genes (raw TPM values)
  # Returns list with:
  #   prelog: matrix (samples × genes) of normalized values before log transformation
  #   means: vector of per‑gene means in prelog space
  #   residuals: list of vectors (residuals per gene) = prelog - mean
  #   replicate_counts: number of replicates per gene (non‑NA)
  #   sorted_order: order of genes by mean (increasing)
  #   means_sorted: sorted means
  #   residuals_sorted: residuals in same order as sorted_order

  # Apply normalisation up to pre‑log step
  if (norm_method == "raw") {
    prelog <- expr_matrix
  } else if (norm_method == "log") {
    prelog <- expr_matrix   # log will be applied later, residuals are on raw scale
  } else if (norm_method == "std_log") {
    # Scale gene‑wise (no quantile)
    original_names <- colnames(expr_matrix)
    scaled <- tox_normalize_by_std_dev(expr_matrix)
    prelog <- scaled
    colnames(prelog) <- original_names
  } else if (norm_method == "full") {
    original_names <- colnames(expr_matrix)
    scaled <- tox_normalize_by_std_dev(expr_matrix)
    quant <- tox_quantile_normalization(scaled)
    prelog <- quant
    colnames(prelog) <- original_names
  } else {
    stop("Unknown norm_method")
  }

  # Get total number of samples
  n_samples <- nrow(prelog)
  
  # Compute means (across replicates, ignoring NAs)
  means <- colMeans(prelog, na.rm = TRUE)
  names(means) <- colnames(prelog)
  
  # Count valid replicates per gene (non-NA values)
  replicate_counts <- colSums(!is.na(prelog))
  names(replicate_counts) <- colnames(prelog)
  
  # Compute residuals per gene
  residuals <- vector("list", ncol(prelog))
  for (g in seq_len(ncol(prelog))) {
    vec <- prelog[, g]
    vec_clean <- vec[!is.na(vec)]
    if (length(vec_clean) > 0) {
      m <- means[g]
      resid <- vec_clean - m
      residuals[[g]] <- resid
    } else {
      residuals[[g]] <- numeric(0)
    }
  }
  names(residuals) <- colnames(prelog)

  # Keep only genes with at least one valid replicate and non-NA mean
  valid_genes <- !is.na(means) & replicate_counts > 0
  
  # Filter to valid genes only
  means_valid <- means[valid_genes]
  replicate_counts_valid <- replicate_counts[valid_genes]
  residuals_valid <- residuals[valid_genes]
  
  # Sort by mean
  sorted_order <- order(means_valid, na.last = NA)
  means_sorted <- means_valid[sorted_order]
  replicate_counts_sorted <- replicate_counts_valid[sorted_order]
  residuals_sorted <- residuals_valid[sorted_order]
  original_indices <- which(valid_genes)[sorted_order]
  list(
    prelog = prelog,
    means = means,
    residuals = residuals,
    replicate_counts = replicate_counts,  # All genes, with 0 for those with no data
    sorted_order = original_indices,
    means_sorted = means_sorted,
    replicate_counts_sorted = replicate_counts_sorted,
    residuals_sorted = residuals_sorted
  )
}

#' Cross-stage expression filter
#'
#' Pools every replicate of a gene across ALL cancer stages (and, by default, the
#' healthy reference) and keeps the gene only if its expression exceeds
#' `min_expr` in at least `min_frac` of those pooled replicates. The check is run
#' on RAW (un-normalized) values and over all stages jointly, so a gene that is
#' expressed in one stage but silent in another is filtered out consistently
#' rather than being kept for some stages only.
#'
#' @param project_id TCGA project id
#' @param stages Cancer stages to pool over (default: all STAGES)
#' @param min_expr Expression threshold; a replicate "counts" if value > min_expr
#' @param min_frac Minimum fraction of pooled replicates that must exceed min_expr
#' @param include_healthy Also require expression in the healthy reference
#' @return Named logical vector keyed by gene_id (TRUE = keep)
compute_gene_keep_mask <- function(project_id, stages = STAGES,
                                   min_expr = 0.5, min_frac = 0.5,
                                   include_healthy = TRUE) {

  load_raw <- function(stage, data_type, constant = FALSE) {
    d <- load_stage_data(project_id, stage, data_type = data_type,
                         use_constant_healthy = constant,
                         norm_method = "raw", apply_mean = FALSE, normalize = FALSE)
    if (is.null(d) || !is.matrix(d$expression_vectors)) return(NULL)
    m <- d$expression_vectors           # samples × genes (raw TPM)
    if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
    m
  }

  mats <- list()
  for (stage in stages) {
    m <- load_raw(stage, "cancer")
    if (!is.null(m)) mats[[length(mats) + 1]] <- m
  }
  if (include_healthy) {
    m <- load_raw(STAGES[1], "healthy", constant = TRUE)
    if (!is.null(m)) mats[[length(mats) + 1]] <- m
  }

  if (length(mats) == 0) return(setNames(logical(0), character(0)))

  # Only consider genes present in EVERY matrix (all cancer stages, and healthy
  # if included); a gene absent from a stage cannot be "expressed across stages"
  # and would otherwise misalign later. Then accumulate, per such gene, the total
  # pooled replicates and how many exceed min_expr.
  common_ids <- Reduce(intersect, lapply(mats, colnames))
  if (length(common_ids) == 0) return(setNames(logical(0), character(0)))

  total_counts <- setNames(numeric(length(common_ids)), common_ids)
  above_counts <- setNames(numeric(length(common_ids)), common_ids)

  for (m in mats) {
    mm <- m[, common_ids, drop = FALSE]
    total_counts <- total_counts + colSums(!is.na(mm))
    above_counts <- above_counts + colSums(mm > min_expr, na.rm = TRUE)
  }

  frac_above <- ifelse(total_counts > 0, above_counts / total_counts, 0)
  keep <- frac_above >= min_frac

  cat(sprintf("    Expression filter (%s): keeping %d / %d genes (>%g in >=%.0f%% of replicates across %d stages%s)\n",
              project_id, sum(keep), length(keep), min_expr, 100 * min_frac,
              length(stages), if (include_healthy) " + healthy" else ""))
  keep
}

#' R wrapper for the noise model
#'
#' Calls the maintained R/Rcpp wrapper `tox_compute_noise_pvalues_pipeline()`
#' (from tensoromics_functions.R) rather than `.Fortran()` directly, so the
#' argument list stays in sync with the current Fortran ABI.
#'
#' NOTE on scales: `cancer_replicates` / `healthy_replicates` are the PRE-LOG
#' (partially normalized) replicate matrices produced by `preprocess_replicates`.
#' For every norm_method other than "raw" the pipeline applies the log2(x+1)
#' transform and builds the null in log2 space internally, so we must NOT log
#' them here. `norm_method` is passed as 0 (linear) for "raw" and non-zero
#' (log2) for the other three methods — the Fortran only distinguishes those two
#' cases.
compute_noise_pvalues <- function(
    cancer_preproc, healthy_preproc, gene_results, family_stats, gene_to_fam,
    stage, norm_method,
    # k_start / k_step / k_max are now counted in GENES, not residuals, so the
    # neighbourhood spans the same expression window regardless of replicate count.
    k_start = 8, k_step = 1, k_max = 15, tau = 0.1,
    max_pool_size = 50000
) {

  n_genes <- nrow(gene_results)

  # Case / control replicate matrices, samples × genes, in PRE-LOG space.
  cancer_means       <- cancer_preproc$means
  cancer_replicates  <- cancer_preproc$prelog
  healthy_means      <- healthy_preproc$means
  healthy_replicates <- healthy_preproc$prelog

  # Observed shift (already in the same log2 coordinate as the null for the
  # log methods; linear for "raw"). The family/ortholog comparisons have been
  # retired from the noise model, so only the `own` statistic is passed.
  stage_col <- gsub(" ", "_", stage)
  obs_own  <- gene_results[[paste0("shift_vs_own_healthy_unscaled_", stage_col)]]
  valid_genes_own  <- as.integer(!is.na(obs_own))
  obs_own[is.na(obs_own)] <- 0

  # The Fortran only distinguishes linear (0) from log2 (non-zero); "log",
  # "std_log" and "full" all map to non-zero.
  norm_method_int <- switch(norm_method,
    "raw"     = 0L,
    "log"     = 1L,
    "std_log" = 1L,
    "full"    = 1L,
    0L
  )

  # Time ONLY the Tensor-Omics pipeline call (no data loading / normalization /
  # post-processing).
  tox_t0 <- Sys.time()
  result <- tox_compute_noise_pvalues_pipeline_exact(
    case_means         = as.numeric(cancer_means),
    case_replicates    = cancer_replicates,
    control_means      = as.numeric(healthy_means),
    control_replicates = healthy_replicates,
    obs_own            = as.numeric(obs_own),
    valid_genes_own    = valid_genes_own,
    norm_method        = norm_method_int,
    k_start            = k_start,
    k_step             = k_step,
    k_max              = k_max,
    tau                = tau,
    max_pool_size      = max_pool_size
  )
  tox_elapsed <- as.numeric(difftime(Sys.time(), tox_t0, units = "secs"))
  cat(sprintf("    [TIMING] Tensor-Omics pipeline call: %.3f s\n", tox_elapsed))

  pvalues_own  <- result$pvalues_own
  # Family / ortholog comparisons retired from the model; kept as all-NA columns
  # so the downstream table shape is unchanged (they were already disabled).
  pvalues_fam  <- rep(NA_real_, n_genes)
  pvalues_orth <- rep(NA_real_, n_genes)

  # The pipeline reports three neighbourhood sizes; the "own" healthy
  # neighbourhood is the control stratum used for the own comparison (own_control).
  neighborhood_size_cancer <- result$neighborhood_size_case
  neighborhood_size_fam    <- rep(NA_integer_, n_genes)
  neighborhood_size_orth   <- rep(NA_integer_, n_genes)
  neighborhood_size_own    <- result$neighborhood_size_own_control

  neighborhood_size_cancer[neighborhood_size_cancer < 0] <- NA
  neighborhood_size_own[neighborhood_size_own < 0] <- NA

  pvalues_own[pvalues_own < 0 | pvalues_own > 1] <- NA
  # pvalues_fam / pvalues_orth are all-NA (family/ortholog retired) -- nothing to clamp.

  pvalues_matrix <- cbind(
    own_healthy   = pvalues_own,
    family_mean   = pvalues_fam,
    ortholog_mean = pvalues_orth
  )
  rownames(pvalues_matrix) <- gene_results$gene_id

  neighborhood_matrix <- cbind(
    neighborhood_size_cancer = neighborhood_size_cancer,
    neighborhood_size_fam    = neighborhood_size_fam,
    neighborhood_size_orth   = neighborhood_size_orth,
    neighborhood_size_own    = neighborhood_size_own
  )
  rownames(neighborhood_matrix) <- gene_results$gene_id

  list(
    pvalues       = pvalues_matrix,
    n_success     = result$n_success,
    neighborhoods = neighborhood_matrix,
    tox_time_secs = tox_elapsed,  # Pure Tensor-Omics call time (no I/O)
    null_sample   = numeric(0)  # For compatibility with existing code
  )
}

# ==================== GESAMTTABELLE ERSTELLEN (MIT STAGE-WEISER VERARBEITUNG) ====================

process_single_combination <- function(cancer_type, project_id, norm_method, stage,
                                       gene_results, family_stats, healthy_preproc,
                                       gene_to_fam, n_families, norm_display,
                                       kept_gene_ids, output_dir) {

  # Load cancer data for this stage
  cancer_data_raw <- load_stage_data(
    project_id = project_id,
    stage = stage,
    data_type = "cancer",
    use_constant_healthy = FALSE,
    norm_method = "raw",
    apply_mean = FALSE,
    normalize = FALSE
  )
  
  if (is.null(cancer_data_raw) || !is.matrix(cancer_data_raw$expression_vectors)) {
    return(NULL)
  }

  # Apply the cross-stage expression filter to the cancer columns (in the same
  # kept-gene order as gene_results) BEFORE normalisation, so the gene-wise
  # std-dev / quantile steps run on the identical filtered gene set used for the
  # healthy side.
  cancer_raw_mat <- cancer_data_raw$expression_vectors
  if (is.null(colnames(cancer_raw_mat)) && !is.null(cancer_data_raw$gene_ids))
    colnames(cancer_raw_mat) <- cancer_data_raw$gene_ids
  missing_genes <- setdiff(kept_gene_ids, colnames(cancer_raw_mat))
  if (length(missing_genes) > 0) {
    cat(sprintf("    Warning: %d kept genes absent from cancer stage %s; skipping\n",
                length(missing_genes), stage))
    return(NULL)
  }
  cancer_raw_mat <- cancer_raw_mat[, kept_gene_ids, drop = FALSE]

  # Preprocess cancer replicates
  cancer_preproc <- preprocess_replicates(cancer_raw_mat, norm_method)
  
  # Compute noise p-values
  noise_result <- compute_noise_pvalues(
    cancer_preproc = cancer_preproc,
    healthy_preproc = healthy_preproc,
    gene_results = gene_results,
    family_stats = family_stats,
    gene_to_fam = gene_to_fam,
    stage = stage,
    norm_method = norm_method,
    k_start = 8,     # neighbour GENES (not residuals)
    k_step = 1,
    k_max = 15,
    tau = 0.1,
    max_pool_size = 50000
  )

  cat(sprintf("    [TIMING] %s / %s / %s -> Tensor-Omics only: %.3f s\n",
              project_id, stage, norm_method, noise_result$tox_time_secs))

  # Extract distance p-values
  distance_p <- extract_distance_pvalues(gene_results)
  family_n_genes <- tabulate(gene_to_fam, nbins = n_families)
  n_genes <- nrow(gene_results)
  
  # Pre-allocate
  results_list <- list()
  row_idx <- 1
  
  # Get neighborhood size matrix
  neighborhood_mat <- noise_result$neighborhoods
  
  for (g in 1:n_genes) {
    gene_id <- gene_results$gene_id[g]
    for (comp_type in COMP_TYPES) {
      dist_p_val <- distance_p[g, stage, comp_type]
      noise_p_val <- noise_result$pvalues[g, comp_type]
      
      if (is.na(dist_p_val) || is.na(noise_p_val)) next
      
      # Scaled with families sds
      shift_col <- switch(comp_type,
                          "own_healthy" = paste0("shift_vs_own_healthy_", gsub(" ", "_", stage)),
                          "family_mean" = paste0("shift_vs_family_mean_", gsub(" ", "_", stage)),
                          "ortholog_mean" = paste0("shift_vs_ortholog_mean_", gsub(" ", "_", stage)))
      
      if (!shift_col %in% colnames(gene_results)) next
      
      shift <- gene_results[g, shift_col]
      direction <- ifelse(is.na(shift), NA, ifelse(shift > 0, "up", "down"))
      
      fam_id <- gene_to_fam[g]
      family_size <- if (!is.na(fam_id) && fam_id > 0) family_n_genes[fam_id] else NA
      family_mean_val <- if (!is.na(fam_id) && fam_id > 0) 
        family_stats$gene_family_means_all[fam_id] else NA
      
      # Get neighborhood sizes for this gene
      n_cancer <- neighborhood_mat[g, "neighborhood_size_cancer"]
      n_healthy <- switch(comp_type,
                          "own_healthy" = neighborhood_mat[g, "neighborhood_size_own"],
                          "family_mean" = neighborhood_mat[g, "neighborhood_size_fam"],
                          "ortholog_mean" = neighborhood_mat[g, "neighborhood_size_orth"])

      results_list[[row_idx]] <- data.frame(
        gene_id = gene_id,
        cancer_type = cancer_type,
        cancer_id = project_id,
        stage = stage,
        normalization = norm_display,
        norm_method = norm_method,
        comparison = comp_type,
        direction = direction,
        signed_difference = shift,
        distance_p_value = dist_p_val,
        noise_p_value = noise_p_val,
        family_size = family_size,
        family_mean = family_mean_val,
        neighborhood_size_cancer = n_cancer,
        neighborhood_size_healthy = n_healthy,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
  
  # Combine results for this combination
  if (length(results_list) == 0) {
    return(NULL)
  }
  
  result_df <- do.call(rbind, results_list)
  
  # ===== SAVE INTERMEDIATE RESULTS =====
  safe_stage <- gsub(" ", "_", stage)
  filename <- sprintf("%s_%s_%s_results.rds", 
                      project_id, norm_method, safe_stage)
  intermediate_file <- file.path(output_dir, "intermediate", filename)
  
  dir.create(dirname(intermediate_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result_df, intermediate_file)
  
  csv_file <- file.path(output_dir, "intermediate", 
                        gsub(".rds", ".csv", filename))
  write.csv(result_df, csv_file, row.names = FALSE)
  
  done_file <- file.path(output_dir, "intermediate", "done", 
                         sprintf("%s_%s_%s.done", 
                                 project_id, norm_method, safe_stage))
  dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste(
    "Completed:", Sys.time(),
    "\nGenes:", nrow(gene_results),
    "\nSuccessfully processed:", noise_result$n_success,
    "\nRows in output:", nrow(result_df)
  ), done_file)
  
  return(result_df)
}

create_all_genes_table <- function(all_results, n_cores = 64, output_dir) {
  
  cat("\n>>> Erstelle Gesamttabelle (stage-weise) mit paralleler Verarbeitung...\n")
  cat(sprintf("    Verwende %d Kerne für parallele Kombinationen\n", n_cores))
  cat(sprintf("    Intermediate Ergebnisse werden gespeichert in: %s/intermediate/\n", output_dir))
  
  # Create all combination tasks
  tasks <- list()
  task_id <- 1
  
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    
    for (stage in STAGES) {
      tasks[[task_id]] <- list(
        cancer_type = res$cancer_type,
        project_id = res$cancer_id,
        norm_method = res$norm_method,
        norm_display = res$norm_display,
        stage = stage,
        gene_results = res$gene_results,
        family_stats = res$family_stats,
        healthy_preproc = res$healthy_preproc,
        gene_to_fam = res$gene_to_fam,
        n_families = res$n_families,
        kept_gene_ids = res$kept_gene_ids,
        output_dir = output_dir  # Pass output directory
      )
      task_id <- task_id + 1
    }
  }
  
  cat(sprintf("    Gesamt: %d Kombinationen zu verarbeiten\n", length(tasks)))

  combination_results <- mclapply(tasks, function(task) {
    tryCatch({
    process_single_combination(
        cancer_type = task$cancer_type,
        project_id = task$project_id,
        norm_method = task$norm_method,
        norm_display = task$norm_display,
        stage = task$stage,
        gene_results = task$gene_results,
        family_stats = task$family_stats,
        healthy_preproc = task$healthy_preproc,
        gene_to_fam = task$gene_to_fam,
        n_families = task$n_families,
        kept_gene_ids = task$kept_gene_ids,
        output_dir = task$output_dir
    )
    }, error = function(e) {
    cat(sprintf("Error processing %s %s %s: %s\n", 
                task$cancer_type, task$norm_method, task$stage, e$message))
    return(NULL)
    })
  }, mc.cores = n_cores, mc.preschedule = TRUE)
  
  # Combine all results
  cat("\n>>> Kombiniere Ergebnisse...\n")
  valid_results <- combination_results[!sapply(combination_results, is.null)]
  cat(sprintf("    Erfolgreich verarbeitet: %d / %d Kombinationen\n", 
              length(valid_results), length(tasks)))
  
  if (length(valid_results) > 0) {
    final_result <- rbindlist(valid_results, fill = TRUE)
    cat(sprintf("    Finale Tabellengröße: %d Zeilen\n", nrow(final_result)))
    return(final_result)
  } else {
    stop("Keine erfolgreichen Ergebnisse!")
  }
}
# ==================== PLOT FUNKTIONEN ====================

plot_jaccard_heatmap <- function(data, title = "", subtitle = "") {
  
  norm_methods <- c("raw", "log", "std_log", "full")
  norm_sets <- list()
  
  for (norm in norm_methods) {
    norm_data <- data %>% filter(norm_method == norm)
    norm_sets[[norm]] <- unique(paste(norm_data$gene_id, norm_data$cancer_id, norm_data$stage))
  }
  
  jaccard_mat <- matrix(NA, 4, 4)
  rownames(jaccard_mat) <- norm_methods
  colnames(jaccard_mat) <- norm_methods
  
  for (i in 1:4) {
    for (j in i:4) {
      if (length(norm_sets[[i]]) > 0 && length(norm_sets[[j]]) > 0) {
        intersection <- length(intersect(norm_sets[[i]], norm_sets[[j]]))
        union <- length(union(norm_sets[[i]], norm_sets[[j]]))
        jaccard_mat[i, j] <- intersection / union
        jaccard_mat[j, i] <- jaccard_mat[i, j]
      }
    }
  }
  
  plot_data <- as.data.frame(as.table(jaccard_mat))
  colnames(plot_data) <- c("Method1", "Method2", "Jaccard")
  plot_data <- plot_data %>% filter(!is.na(Jaccard))
  
  ggplot(plot_data, aes(x = Method1, y = Method2, fill = Jaccard)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 5) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==================== HAUPTSCHLEIFE ====================

analyze_outlier_significance <- function(output_dir = NULL, n_cores = 64) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("OUTLIER SIGNIFICANCE ANALYSIS (parallel combinations)\n")
  cat(sprintf("Using %d cores for parallel combination processing\n", n_cores))
  cat("Parameters (neighbourhood in genes): k_start=8, k_step=1, k_max=15, tau=0.1\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Load all data (this is still sequential, but only once)
  all_results <- load_all_gene_results()
  if (length(all_results) == 0) stop("Keine Daten gefunden!")
  
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(TOX_TEST_DIR), "outlier_significance")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create results table in parallel
  cat("\n>>> Erstelle Gesamttabelle (parallel)...\n")
  all_genes_table <- create_all_genes_table(all_results, n_cores = n_cores, output_dir = output_dir)
  
  all_genes_table <- all_genes_table %>%
    group_by(cancer_id, normalization, stage, comparison) %>%
    mutate(
      noise_p_value_adj = p.adjust(noise_p_value, method = "BH")
    ) %>%
    ungroup()

  # Save results
  saveRDS(all_genes_table, file.path(output_dir, "all_genes_pvalues.rds"))
  write.csv(all_genes_table, file.path(output_dir, "all_genes_pvalues.csv"), row.names = FALSE)
  
  # ===== 3. BASELINE AUSGABE =====
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("UNKORRIGIERTE WERTE (BASELINE)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  cat("\n--- Summary distance_p_values ---\n")
  print(summary(all_genes_table$distance_p_value))
  cat("\n--- Summary noise_p_values ---\n")
  print(summary(all_genes_table$noise_p_value))
  
  # ===== 4. SELEKTION: NOISE-P-WERT < 0.05 =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SELEKTION: noise_p_adj < 0.01 (sortiert nach distance_p_value)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  ALPHA_NOISE <- 0.01
  
  significant_by_norm <- list()
  
  for (current_norm in unique(all_genes_table$norm_method)) {
    
    cat(sprintf("\n>>> Verarbeite Normalisierung: %s <<<\n", current_norm))
    
    norm_data <- all_genes_table %>% 
      filter(norm_method == current_norm)
    
    # ===== DIAGNOSE PRO VERGLEICHSTYP =====
    cat(sprintf("\n  DIAGNOSE %s:\n", current_norm))
    
    for (comp_type in COMP_TYPES) {
      comp_data <- norm_data %>% filter(comparison == comp_type)
      
      cat(sprintf("\n    --- %s ---\n", comp_type))
      cat(sprintf("      Noise p (unkorr) < 0.01: %d (%.1f%%)\n", 
                  sum(comp_data$noise_p_value < 0.01, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value < 0.01, na.rm = TRUE)))
      cat(sprintf("      Noise p (adj) < 0.01: %d (%.1f%%)\n\n", 
                  sum(comp_data$noise_p_value_adj < 0.01, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value_adj < 0.01, na.rm = TRUE)))
      cat("      Summary of noise_p_value: \n")
      print(summary(comp_data$noise_p_value))
      cat("      Summary of noise_p_value_adj: \n")
      print(summary(comp_data$noise_p_value_adj))
    }
    
    # ===== SELEKTION: Nur noise_p_adj < 0.01 =====
    significant <- norm_data %>%
      filter(noise_p_value_adj < ALPHA_NOISE) %>%
      arrange(distance_p_value)
    
    cat(sprintf("\n  → Signifikante Einträge (noise_p_adj < 0.01): %d\n", nrow(significant)))
    
    if (nrow(significant) > 0) {
      cat("\n  Top 10 nach distance_p_value:\n")
      print(significant %>%
        select(gene_id, stage, comparison, distance_p_value, noise_p_value_adj) %>%
        head(10))
      
      filename <- sprintf("significant_%s.csv", 
                         gsub("[^a-zA-Z0-9]", "_", current_norm))
      write.csv(significant, file.path(output_dir, filename), row.names = FALSE)
      
      significant_by_norm[[current_norm]] <- significant
    }
  }
  
  # ===== 5. GESAMTTABELLE ÜBER ALLE NORMALISIERUNGEN =====
  cat("\n>>> Gesamttabelle über alle Normalisierungen <<<\n")
  
  if (length(significant_by_norm) > 0) {
    all_significant <- rbindlist(significant_by_norm, fill = TRUE)
    all_significant <- all_significant %>% arrange(distance_p_value)
    
    cat(sprintf("\n  Gesamt: %d signifikante Einträge (noise_p_adj < 0.01)\n", nrow(all_significant)))
    
    # Statistik pro Normalisierung
    cat("\n  Statistik pro Normalisierung:\n")
    norm_stats <- all_significant %>%
      group_by(norm_method) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE),
        q95_distance_p = quantile(distance_p_value, 0.95, na.rm = TRUE),
        min_distance_p = min(distance_p_value, na.rm = TRUE)
      ) %>%
      arrange(desc(n))
    print(norm_stats)
    
    # Statistik pro Vergleichstyp
    cat("\n  Statistik pro Vergleichstyp:\n")
    comp_stats <- all_significant %>%
      group_by(comparison) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE),
        q95_distance_p = quantile(distance_p_value, 0.95, na.rm = TRUE)
      ) %>%
      arrange(desc(n))
    print(comp_stats)
    
    # Statistik pro Krebsart
    cat("\n  Statistik pro Krebsart (Top 10):\n")
    cancer_stats <- all_significant %>%
      group_by(cancer_type) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE)
      ) %>%
      arrange(desc(n)) %>%
      head(10)
    print(cancer_stats)
    
    write.csv(all_significant, 
              file.path(output_dir, "significant_all.csv"), 
              row.names = FALSE)
    saveRDS(all_significant, 
            file.path(output_dir, "significant_all.rds"))
  } else {
    cat("\n  Keine signifikanten Einträge gefunden!\n")
    all_significant <- data.frame()
  }
  
  # ===== 6. JACCARD PLOTS =====
  cat("\n>>> Jaccard Plots (basierend auf noise_p_adj < 0.01)\n")
  
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (nrow(all_significant) > 0) {
    for (min_size in MIN_FAMILY_SIZES) {
      for (comp_type in COMP_TYPES) {
        
        plot_data <- all_significant %>%
          filter(
            comparison == comp_type, 
            family_size >= min_size
          )
        
        if (nrow(plot_data) > 0) {
          p_adj <- plot_jaccard_heatmap(
            plot_data,
            title = sprintf("%s (min %d) - KORRIGIERT", comp_type, min_size),
            subtitle = sprintf("noise_p_adj < 0.01, n=%d", nrow(plot_data))
          )
          
          filename_adj <- sprintf("jaccard_ADJ_%s_min%d.png", comp_type, min_size)
          ggsave(file.path(plots_dir, filename_adj), p_adj, 
                 width = 7, height = 6, dpi = 300)
          
          # Rohdaten für Vergleich
          raw_data <- all_genes_table %>%
            filter(
              comparison == comp_type,
              family_size >= min_size,
              noise_p_value < 0.01
            )
          
          if (nrow(raw_data) > 0 && nrow(plot_data) > 0) {
            norm_methods <- c("raw", "log", "std_log", "full")
            
            raw_sets <- list()
            for (norm in norm_methods) {
              norm_data <- raw_data %>% filter(norm_method == norm)
              raw_sets[[norm]] <- unique(paste(norm_data$gene_id, 
                                               norm_data$cancer_id, 
                                               norm_data$stage))
            }
            
            adj_sets <- list()
            for (norm in norm_methods) {
              norm_data <- plot_data %>% filter(norm_method == norm)
              adj_sets[[norm]] <- unique(paste(norm_data$gene_id, 
                                               norm_data$cancer_id, 
                                               norm_data$stage))
            }
            
            comparison_mat <- matrix(NA, 4, 4)
            rownames(comparison_mat) <- norm_methods
            colnames(comparison_mat) <- norm_methods
            
            for (i in 1:4) {
              for (j in 1:4) {
                if (length(raw_sets[[i]]) > 0 && length(adj_sets[[j]]) > 0) {
                  intersection <- length(intersect(raw_sets[[i]], adj_sets[[j]]))
                  union <- length(union(raw_sets[[i]], adj_sets[[j]]))
                  comparison_mat[i, j] <- intersection / union
                }
              }
            }
            
            plot_data_comp <- as.data.frame(as.table(comparison_mat))
            colnames(plot_data_comp) <- c("Raw_Norm", "Adj_Norm", "Jaccard")
            plot_data_comp <- plot_data_comp %>% filter(!is.na(Jaccard))
            
            p_comp <- ggplot(plot_data_comp, aes(x = Adj_Norm, y = Raw_Norm, fill = Jaccard)) +
              geom_tile() +
              geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 4) +
              scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
              labs(title = sprintf("%s (min %d) - Raw data vs. noise_p_adj < 0.01", 
                                   comp_type, min_size),
                   x = "noise_p_adj < 0.01", y = "Raw data (noise_p < 0.01)") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            filename_comp <- sprintf("jaccard_COMP_%s_min%d.png", comp_type, min_size)
            ggsave(file.path(plots_dir, filename_comp), p_comp, 
                   width = 7, height = 6, dpi = 300)
          }
        }
      }
    }
  }
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("ANALYSE ABGESCHLOSSEN\n")
  cat(sprintf("Ergebnisse in: %s\n", output_dir))
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  return(invisible(list(
    data = all_genes_table,
    significant = all_significant
  )))
}

# ==================== AUSFÜHRUNG ====================

if (sys.nframe() == 0) {
  set.seed(42)
  base_dir <- file.path(dirname(TOX_TEST_DIR), "outlier_significance_small_neighborhoods_sqrt")
  
  cat("\n", paste(rep("#", 80), collapse = ""), "\n")
  cat("### ADAPTIVE NOISE MODEL ANALYSIS (STAGE-WISE) ###\n")
  cat("Parameters (neighbourhood in genes): k_start=8, k_step=1, k_max=15, tau=0.1\n")
  cat(paste(rep("#", 80), collapse = ""), "\n")
  
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- analyze_outlier_significance(output_dir = base_dir)
  
  saveRDS(list(
    timestamp = Sys.time(),
    n_significant = nrow(results$significant),
    parameters = list(
      k_start = 5,
      k_step = 1,
      k_max = 13,
      tau = 0.05
    )
  ), file.path(base_dir, "metadata.rds"))
}