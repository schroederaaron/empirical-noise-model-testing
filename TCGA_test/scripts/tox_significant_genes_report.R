#!/usr/bin/env Rscript
# tox_significant_genes_report.R
# -----------------------------------------------------------------------------
# The genes TOX-log still calls after BH, reported with everything we know about
# them and about the null they were scored against.
#
# WHY A SEPARATE SCRIPT
#   The production table (outlier_significance_analysis.R / compare_noise_models.R)
#   carries gene_id, direction, signed_difference, the two p-values and two
#   neighbourhood SIZES. It carries no expression level, no variance, no replicate
#   counts, and nothing about which genes the null was built from -- so a
#   BH survivor cannot currently be judged: there is no way to see whether it is a
#   highly expressed gene with a modest shift, a quiet gene with a large one, or a
#   gene whose neighbourhood was too quiet to score it against.
#
# WHAT IS REPORTED
#   p-values come from the PRODUCTION Fortran pipeline via compute_noise_pvalues(),
#   so the calls here are exactly the project's calls, not a re-derivation. Every
#   added statistic comes from tox_null_reimpl.R, whose agreement with the Fortran
#   is checked at the top of the run and reported before any table is written.
#
#   per survivor, three blocks:
#     EXPRESSION  mean / median / sd / var / CV in both groups, on the TPM scale
#                 and on the log2 scale the log model works in; log2 fold change;
#                 the observed statistic the model actually scored
#     NULL        the neighbourhood the gene was scored against -- pool sd, the
#                 gene's own sd, their ratio rho, the pool's excess kurtosis, the
#                 implied null sd, and the observed distance expressed both as
#                 z_pool (against the pooled null) and z_own (against the gene's
#                 own replicates). rho far from 1 means the null does not describe
#                 this gene's noise; z_own >> z_pool or vice versa says which of
#                 the two disagrees
#     INFERENCE   raw p, BH q, rank, and the RESOLUTION FLOOR
#
#   ...plus a companion long table with one row per (survivor, neighbour gene):
#   the neighbour's id, mean, mean-distance and residual sd. Two survivors that
#   look alike can have very different neighbourhoods, and that is visible only here.
#
# THE RESOLUTION FLOOR -- read this before interpreting any q-value
#   The exact p-value is a pairwise tail count, p = (count + 1) / (n_a * n_b + 1),
#   so the SMALLEST p a gene can receive is p_floor = 1 / (n_a * n_b + 1), set by
#   its own pool sizes. A gene at that floor is CENSORED: its true p may be far
#   smaller, and it is reported as `at_floor = TRUE`.
#
#   The floor does NOT make BH impossible -- it makes rejection QUANTISED. With r
#   genes tied at the floor, BH rejects the whole block as soon as
#   p_floor <= q * r / G, i.e.  r >= p_floor * G / q. At n_rep = 3, k = 50 genes
#   (22,500 pairs, p_floor = 4.44e-5), G = 12,000, q = 0.05 that is r >= 11 --
#   confirmed against p.adjust(): 0 rejections at r = 10, exactly 11 at r = 11.
#   For a SINGLE gene to be rejectable alone at n_rep = 3 you need k >= 164 genes
#   (q = 0.05) or k >= 366 (q = 0.01); at n_rep = 10, k >= 49 / 110.
#
#   Under H0 the expected number of floor genes is G * p_floor ~ 0.53 at k = 50,
#   n_rep = 3, so the r >= 11 block essentially never forms by chance. hits_FDR05
#   = 0 in the split-half calibration runs is therefore the EXPECTED conservative
#   consequence of the floor -- it is not a bug, and it is not evidence of good
#   calibration either.
#
#   A second, softer limit binds long before the floor: the pairs are not
#   independent. 22,500 pairs come from 150 + 150 residuals, and at n_rep = 3 each
#   gene carries only 2 df. Measured by resampling neighbourhoods, the effective
#   sample size is 145-595 depending on how heterogeneous the pool is -- 0.6-2.6%
#   of the nominal pair count, so the printed floor overstates the real resolution
#   by 40-160x. Growing k improves the floor as k^2 but precision only as ~k, so it
#   widens that gap; replicates are what actually buy resolution (at a fixed
#   150-residual pool, n_rep 3 -> 10 raises N_eff from 216 to 3,514).
#   Below roughly 1/N_eff a p-value reports neighbourhood sampling noise rather
#   than evidence about the gene -- around 5e-3 at n_rep = 3, k = 50.
#   The run prints the floor, the BH threshold and the number of genes sitting at
#   the floor per cohort, so an empty survivor list can be told apart from a
#   genuinely null cohort.
#
# RUN:  Rscript tox_significant_genes_report.R          (from TCGA_test/scripts)
# -----------------------------------------------------------------------------

# --- persistent R package library (MUST precede any library()/require()/source() call) ---
# `.libPaths()` silently DROPS non-existent directories, so the folder has to be created
# FIRST or the prepend is a no-op and packages land in the ephemeral container library.
LIB_DIR <- normalizePath("external/docker_r_libs", mustWork = FALSE)
if (!dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE) && !dir.exists(LIB_DIR))
  stop("Could not create package library ", LIB_DIR,
       " -- it must be on a WRITABLE, BIND-MOUNTED path or packages will not persist.")
.libPaths(c(LIB_DIR, .libPaths()))

suppressMessages({library(dplyr); library(data.table)})
source("outlier_significance_analysis.R")
source("config.R")
# tox_null_reimpl.R lives in common/; the repo's scripts are run from more than one
# working directory, so look for it rather than assuming one.
.source_reimpl <- function() {
  cands <- c("tox_null_reimpl.R", "common/tox_null_reimpl.R", "../../common/tox_null_reimpl.R",
             "../common/tox_null_reimpl.R", "../../../common/tox_null_reimpl.R")
  for (p in cands) if (file.exists(p)) { source(p); return(invisible(p)) }
  stop("Could not locate tox_null_reimpl.R (looked in: ", paste(cands, collapse = ", "), ")")
}
.source_reimpl()
options(width = 220)

OUT_DIR <- file.path(dirname(TOX_TEST_DIR), "significant_genes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==================== CONFIG ====================
REPORT_NORM   <- "log"          # the arm that passes calibration in all 12 groups
NORM_INT      <- 1L             # what the Fortran calls "log"
REPORT_STAGES <- STAGES
Q_CUTOFFS     <- c(0.05, 0.01)  # survivors reported at both; TOP-level cut is the first
K_START <- 20L; K_STEP <- 1L; K_MAX <- 50L; TAU <- 0.1; MAX_POOL <- 70000L
N_NEIGHBOURS_REPORTED <- 50L    # cap on rows written per survivor per side
set.seed(42)

# ==================== helpers ====================

#' Descriptive statistics of one gene in one group, on both scales.
#' `sd` is the ordinary unbiased (Bessel) sd, i.e. what the residual pool is
#' built to represent. `cv` is undefined at mean 0 and returned as NA rather than Inf.
gene_stats <- function(v) {
  v <- v[is.finite(v)]; n <- length(v)
  if (n < 2L) return(c(n = n, mean = NA, median = NA, sd = NA, var = NA, cv = NA,
                       mean_log2 = NA, sd_log2 = NA, q25 = NA, q75 = NA, min = NA, max = NA))
  lv <- log2(pmax(v, 0) + 1)
  c(n = n, mean = mean(v), median = median(v), sd = sd(v), var = var(v),
    cv = if (mean(v) > 0) sd(v) / mean(v) else NA_real_,
    mean_log2 = mean(lv), sd_log2 = sd(lv),
    q25 = unname(quantile(v, .25)), q75 = unname(quantile(v, .75)),
    min = min(v), max = max(v))
}

# ==================== load ====================
cat("Loading gene results / healthy reference for norm =", REPORT_NORM, "...\n")
ALL <- load_all_gene_results()
keys <- grep(paste0("_", REPORT_NORM, "$"), names(ALL), value = TRUE)
if (!length(keys)) stop("No datasets loaded for norm_method = ", REPORT_NORM)

# ==================== Part 0: fidelity gate ====================
# The enrichment statistics come from the R port; the p-values come from the
# Fortran. They must describe the same neighbourhoods or the two halves of each
# row would not belong together. Checked once, on the first available cohort.
cat("\n=== Fidelity check: R port vs Fortran ===\n")
{
  d0 <- ALL[[keys[1]]]
  c0 <- tryCatch(load_stage_data(d0$cancer_id, REPORT_STAGES[1], "cancer",
                                 norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
                 error = function(e) NULL)
  if (!is.null(c0) && is.matrix(c0$expression_vectors)) {
    cm <- c0$expression_vectors
    if (is.null(colnames(cm)) && !is.null(c0$gene_ids)) colnames(cm) <- c0$gene_ids
    cm <- cm[, intersect(d0$kept_gene_ids, colnames(cm)), drop = FALSE]
    hm <- d0$healthy_replicates_raw[, colnames(cm), drop = FALSE]
    v <- validate_against_fortran(cm, hm, NORM_INT, K_START, K_STEP, K_MAX, TAU, 0.0, MAX_POOL)
    cat(sprintf("  compared %d genes | max|dp| = %.3g | max|d nbhd| = %d | ok = %s\n",
                v$n_compared, v$max_abs_diff, v$max_abs_diff_nbhd, v$ok))
    if (!isTRUE(v$ok))
      warning("R port does not reproduce the Fortran exactly -- the NULL block of every ",
              "row below is provisional (the p / q columns come from the Fortran and are unaffected).")
  } else cat("  reference cohort unavailable -- gate SKIPPED.\n")
}

# ==================== main loop ====================
sig_rows <- list(); nb_rows <- list(); cohort_rows <- list()

for (k in keys) {
  d <- ALL[[k]]
  cat(sprintf("\n>>> %s (%s)\n", d$cancer_type, d$cancer_id))

  for (stage in REPORT_STAGES) {
    cd <- tryCatch(load_stage_data(d$cancer_id, stage, "cancer", use_constant_healthy = FALSE,
                                   norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
                   error = function(e) NULL)
    if (is.null(cd) || !is.matrix(cd$expression_vectors)) next
    cm <- cd$expression_vectors
    if (is.null(colnames(cm)) && !is.null(cd$gene_ids)) colnames(cm) <- cd$gene_ids
    common <- intersect(d$kept_gene_ids, colnames(cm))
    if (length(common) < 100L) next
    cm <- cm[, common, drop = FALSE]
    hm <- d$healthy_replicates_raw[, common, drop = FALSE]
    gr <- d$gene_results[match(common, d$gene_results$gene_id), , drop = FALSE]

    cancer_preproc <- preprocess_replicates(cm, REPORT_NORM)

    # --- production p-values (Fortran), exactly as the project computes them ---
    res <- tryCatch(compute_noise_pvalues(
             cancer_preproc = cancer_preproc, healthy_preproc = preprocess_replicates(hm, REPORT_NORM),
             gene_results = gr, family_stats = d$family_stats, gene_to_fam = NULL,
             stage = stage, norm_method = REPORT_NORM,
             k_start = K_START, k_step = K_STEP, k_max = K_MAX, tau = TAU,
             max_pool_size = MAX_POOL), error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(res)) next

    p <- as.numeric(res$pvalues[, "own_healthy"])
    tested <- is.finite(p)
    # BH over TESTED genes only. Including not-computed genes in the denominator
    # would silently make the correction more conservative than it should be.
    q <- rep(NA_real_, length(p)); q[tested] <- p.adjust(p[tested], method = "BH")

    # --- BH cliff bookkeeping (see the header note on the resolution floor) ---
    nb_c <- res$neighborhoods[, "neighborhood_size_cancer"]
    nb_h <- res$neighborhoods[, "neighborhood_size_own"]
    floor_p <- 1 / (as.numeric(nb_c) * as.numeric(nb_h) + 1)
    G <- sum(tested)
    bh_thresh <- Q_CUTOFFS[1] / max(G, 1L)                # what p_(1) must beat
    at_floor <- tested & is.finite(floor_p) & (p <= floor_p * (1 + 1e-9))
    cohort_rows[[length(cohort_rows) + 1L]] <- data.frame(
      cancer = d$cancer_type, cancer_id = d$cancer_id, stage = stage, norm = REPORT_NORM,
      n_case = nrow(cm), n_control = nrow(hm), n_genes_tested = G,
      median_pool_case = median(nb_c, na.rm = TRUE), median_pool_control = median(nb_h, na.rm = TRUE),
      median_p_floor = median(floor_p, na.rm = TRUE),
      bh_threshold_p1 = bh_thresh,
      floor_blocks_BH = median(floor_p, na.rm = TRUE) > bh_thresh,
      n_genes_at_floor = sum(at_floor, na.rm = TRUE),
      min_p = suppressWarnings(min(p[tested])),
      n_q05 = sum(q < 0.05, na.rm = TRUE), n_q01 = sum(q < 0.01, na.rm = TRUE),
      stringsAsFactors = FALSE)
    cat(sprintf("   %-10s n_case=%3d n_ctrl=%3d tested=%5d  min_p=%.2e  p_floor(med)=%.2e  BH needs<=%.2e  q<.05: %d\n",
                stage, nrow(cm), nrow(hm), G, min(p[tested]), median(floor_p, na.rm = TRUE),
                bh_thresh, sum(q < 0.05, na.rm = TRUE)))

    # Rank once per cohort, not once per survivor: rank() over ~12k genes inside
    # the per-gene loop would be O(survivors * G log G) for no reason.
    rk <- rep(NA_integer_, length(p)); rk[tested] <- rank(p[tested], ties.method = "min")

    sel <- which(!is.na(q) & q < max(Q_CUTOFFS))
    if (!length(sel)) next

    # --- enrichment: pool / neighbourhood statistics for the survivors only ---
    dg <- tryCatch(tox_diagnose(cm, hm, NORM_INT, K_START, K_STEP, K_MAX, TAU, 0.0,
                                MAX_POOL, genes = common[sel], verbose = FALSE),
                   error = function(e) NULL)
    dgi <- if (!is.null(dg)) match(common[sel], dg$gene_id) else rep(NA_integer_, length(sel))

    for (j in seq_along(sel)) {
      i <- sel[j]; gidx <- common[i]
      sc <- gene_stats(cm[, i]); sh <- gene_stats(hm[, i])
      dd <- if (!is.na(dgi[j])) dg[dgi[j], ] else NULL
      sig_rows[[length(sig_rows) + 1L]] <- data.frame(
        # ---- identity / design
        gene_id = gidx, cancer = d$cancer_type, cancer_id = d$cancer_id, stage = stage,
        norm = REPORT_NORM, n_case = unname(sc["n"]), n_control = unname(sh["n"]),
        # ---- expression, TPM scale
        mean_case = unname(sc["mean"]), mean_control = unname(sh["mean"]),
        median_case = unname(sc["median"]), median_control = unname(sh["median"]),
        sd_case = unname(sc["sd"]), sd_control = unname(sh["sd"]),
        var_case = unname(sc["var"]), var_control = unname(sh["var"]),
        cv_case = unname(sc["cv"]), cv_control = unname(sh["cv"]),
        q25_case = unname(sc["q25"]), q75_case = unname(sc["q75"]),
        q25_control = unname(sh["q25"]), q75_control = unname(sh["q75"]),
        min_case = unname(sc["min"]), max_case = unname(sc["max"]),
        min_control = unname(sh["min"]), max_control = unname(sh["max"]),
        # ---- expression, log2 scale (the scale the log model scores on)
        mean_log2_case = unname(sc["mean_log2"]), mean_log2_control = unname(sh["mean_log2"]),
        sd_log2_case = unname(sc["sd_log2"]), sd_log2_control = unname(sh["sd_log2"]),
        log2FC = unname(sc["mean_log2"] - sh["mean_log2"]),
        obs_statistic = if (!is.null(dd)) dd$obs else NA_real_,
        direction = if (!is.null(dd)) ifelse(dd$obs > 0, "up", "down")
                    else ifelse(sc["mean_log2"] > sh["mean_log2"], "up", "down"),
        # ---- the null this gene was scored against
        sd_pool_case = if (!is.null(dd)) dd$sd_pool_case else NA_real_,
        sd_pool_control = if (!is.null(dd)) dd$sd_pool_control else NA_real_,
        sd_resid_own_case = if (!is.null(dd)) dd$sd_own_case else NA_real_,
        sd_resid_own_control = if (!is.null(dd)) dd$sd_own_control else NA_real_,
        rho_case = if (!is.null(dd)) dd$rho_case else NA_real_,
        rho_control = if (!is.null(dd)) dd$rho_control else NA_real_,
        kurt_pool_case = if (!is.null(dd)) dd$kurt_pool_case else NA_real_,
        kurt_pool_control = if (!is.null(dd)) dd$kurt_pool_control else NA_real_,
        sd_null = if (!is.null(dd)) dd$sd_null else NA_real_,
        z_pool = if (!is.null(dd)) dd$z_pool else NA_real_,
        z_own = if (!is.null(dd)) dd$z_own else NA_real_,
        # ---- neighbourhood
        neighborhood_size_case = unname(nb_c[i]), neighborhood_size_control = unname(nb_h[i]),
        n_genes_pool_case = if (!is.null(dd)) dd$n_genes_pool_case else NA_integer_,
        n_genes_pool_control = if (!is.null(dd)) dd$n_genes_pool_control else NA_integer_,
        mean_span_case = if (!is.null(dd)) dd$mean_span_case else NA_real_,
        mean_span_control = if (!is.null(dd)) dd$mean_span_control else NA_real_,
        frac_below_case = if (!is.null(dd)) dd$frac_below_case else NA_real_,
        frac_below_control = if (!is.null(dd)) dd$frac_below_control else NA_real_,
        stop_case = if (!is.null(dd)) dd$stop_case else NA_character_,
        stop_control = if (!is.null(dd)) dd$stop_control else NA_character_,
        # ---- inference
        p_raw = p[i], q_BH = q[i], p_floor = unname(floor_p[i]), at_floor = unname(at_floor[i]),
        rank_in_cohort = unname(rk[i]),
        n_genes_tested = G,
        stringsAsFactors = FALSE)
    }

    # --- neighbour membership for the top survivors (both sides) ---
    top <- sel[order(p[sel])][seq_len(min(length(sel), 200L))]
    for (i in top) for (side in c("case", "control")) {
      nbt <- tryCatch(tox_neighbours(common[i], cm, hm, NORM_INT, side,
                                     K_START, K_STEP, K_MAX, TAU, MAX_POOL),
                      error = function(e) NULL)
      if (is.null(nbt)) next
      nbt <- head(nbt, N_NEIGHBOURS_REPORTED)
      nbt$cancer <- d$cancer_type; nbt$stage <- stage
      nb_rows[[length(nb_rows) + 1L]] <- nbt
    }
  }
}

# ==================== output ====================
COH <- if (length(cohort_rows)) as.data.frame(rbindlist(cohort_rows)) else NULL
if (!is.null(COH)) {
  cat("\n=== Cohort-level summary (BH feasibility per cohort) ===\n")
  cat("floor_blocks_BH = TRUE means the pairwise p-value floor is ABOVE the p-value BH would\n")
  cat("need from the single best gene, so an empty survivor list carries no evidence either way.\n")
  print(COH, row.names = FALSE)
  write.csv(COH, file.path(OUT_DIR, "cohort_summary.csv"), row.names = FALSE)
}

if (!length(sig_rows)) {
  cat("\nNo gene survived BH at q <", max(Q_CUTOFFS), "in any cohort.\n")
  cat("Check `floor_blocks_BH` above before reading that as an absence of signal.\n")
} else {
  SIG <- as.data.frame(rbindlist(sig_rows, fill = TRUE))
  SIG <- SIG[order(SIG$cancer, SIG$stage, SIG$q_BH, SIG$p_raw), ]
  write.csv(SIG, file.path(OUT_DIR, "significant_genes_full_stats.csv"), row.names = FALSE)
  for (a in Q_CUTOFFS)
    write.csv(SIG[SIG$q_BH < a, ],
              file.path(OUT_DIR, sprintf("significant_genes_q%03d.csv", round(a * 100))),
              row.names = FALSE)

  cat(sprintf("\n=== %d gene x cohort calls at q < %.2f (%d at q < %.2f) ===\n",
              sum(SIG$q_BH < Q_CUTOFFS[1]), Q_CUTOFFS[1],
              sum(SIG$q_BH < Q_CUTOFFS[2]), Q_CUTOFFS[2]))
  cat("\n-- per cohort --\n")
  print(SIG %>% filter(q_BH < Q_CUTOFFS[1]) %>% group_by(cancer, stage) %>%
          summarise(n = dplyr::n(), n_at_floor = sum(at_floor, na.rm = TRUE),
                    up = sum(direction == "up"), down = sum(direction == "down"),
                    median_log2FC = round(median(log2FC), 3),
                    median_mean_case = round(median(mean_case), 2),
                    median_rho_case = round(median(rho_case, na.rm = TRUE), 3),
                    median_z_pool = round(median(z_pool, na.rm = TRUE), 2),
                    median_z_own = round(median(z_own, na.rm = TRUE), 2),
                    .groups = "drop") %>% as.data.frame(), row.names = FALSE)

  cat("\n-- strongest 25 calls --\n")
  show <- c("gene_id", "cancer", "stage", "mean_case", "mean_control", "log2FC",
            "sd_case", "sd_control", "rho_case", "kurt_pool_case",
            "z_pool", "z_own", "neighborhood_size_case", "p_raw", "q_BH", "at_floor")
  print(head(SIG[order(SIG$p_raw), intersect(show, names(SIG))], 25), row.names = FALSE, digits = 4)

  # rho far from 1 = the null does not describe this gene's own noise. rho < 1 means
  # the gene is QUIETER than its neighbourhood, so it was scored against a null that
  # is too wide (a conservative call, safe); rho > 1 means the opposite, and those
  # calls are the ones to treat with caution.
  cat("\n-- calls whose own noise disagrees most with their null (|log rho| largest) --\n")
  SIG$abs_log_rho <- abs(log(SIG$rho_case))
  print(head(SIG[order(-SIG$abs_log_rho), intersect(c(show, "abs_log_rho"), names(SIG))], 15),
        row.names = FALSE, digits = 4)
}

if (length(nb_rows)) {
  NB <- as.data.frame(rbindlist(nb_rows, fill = TRUE))
  write.csv(NB, file.path(OUT_DIR, "significant_genes_neighbourhoods.csv"), row.names = FALSE)
  cat(sprintf("\nNeighbourhood membership: %d rows for %d target genes.\n",
              nrow(NB), length(unique(NB$target_gene))))
}
cat("\nWritten to:", OUT_DIR, "\n")
