#!/usr/bin/env Rscript
# power_report.R
# =============================================================================
# POWER BENCHMARK -- INTERPRETATION script. Reads the raw per-gene output of
# power_test.R (power_raw/<part>/{manifest.csv, truth_*, results_*, diag_*}),
# computes every metric (common/power_metrics.R), writes the tables and every plot.
# Needs no Fortran and no reference-method run; re-run it after any metric change.
#
#   Rscript TCGA_test/scripts/power_report.R
#   POWER_RAW_DIR=... POWER_OUT_DIR=... POWER_CORES=4 Rscript TCGA_test/scripts/power_report.R
#
# METRICS (per method, per dataset; mean +/- SE over rounds)
#   Threshold behaviour
#     fdr_nom_a        observed FDP of BH calls at nominal a, ZERO-FILLED when a
#                      method makes no call (so the mean is the FDR BH controls, I1)
#     fdr_cond_a       the old conditional E[V/R | R > 0] (NA without calls)
#     zero_disc_a      1 = no call in that round; its mean is frac_zero_disc
#     tpr_nom_a        TPR of BH calls at nominal a;  tpr_all_* counts filtered-out DE as missed
#     tpr_ach05        TPR at ACHIEVED FDR <= 0.05, truth-thresholded: the maximum over
#                      all steps with realised FDR <= 0.05 (slightly optimistic on a
#                      non-monotone curve, not user-realisable). tpr_nom_* is the
#                      user-facing number.
#     fpr_null05       fraction of true nulls called at nominal 0.05 (BH)
#     p_le05_null      fraction of true nulls with raw p <= 0.05 (calibration; P+/D null arms)
#     fpr_hv05         (Part V) fraction of variance-only nulls called
#   Ranking (threshold-free); *_p ranked by p alone (ties stay tied), *_pe p then
#     |log2 FC of group means|. The gap is the resolution loss of the discrete p grid.
#     auc, pauc05, pauc10 (McClish: 0.5 random, 1 perfect)
#     auprc            step-function AUPRC (average precision), no interpolation
#     auprc_all        recall denominator n_de_total (filtered-out DE never retrieved)
#     prevalence       n_de / n_genes ranked: the random-ranking AUPRC
#     auprc_excess_pe  auprc_pe - prevalence (comparable across pi1);  auprc_lift_pe ratio
#     prec_top{k}, prec_top_nde (precision at k = n_de, R-precision)
#   Direction / effect size: sign_err05, lfc_rmse (NA for linear-scale TOX-raw), lfc80
#   Resolution: p_floor, n_at_floor, m_star, floor_blocks05, floor_th_median,
#     nbhd_own_*_median (TOX pool size in residuals)
#
# OUTPUTS (POWER_OUT_DIR)
#   power_per_round.csv       one row per (job, method); manifest columns included
#   power_summary.csv         mean / SE over rounds per (cell, method); the variance-band
#                             arms add frac_short_pool, n_scanned_median, frac_skipped,
#                             tau_var, d0 (finite only), d0_inf, w
#   power_paired_vs_tox.csv   variance-band arms minus their plain TOX comparator, paired
#                             per round (P+, V: TOX-log; R, RS: TOX-log|TPM-centred)
#   power_paired_B_vs_A.csv   model B minus model A at equal c (B c1 vs A c1.5, c2, c3, c4),
#                             paired per round, with SE and t
#   power_rs_gates.csv        Part RS decision protocol (note D section 6), steps 1-4
#   power_rs_checks.csv       Part RS per (cohort, n) checks: composition shift, overlap
#   power_scan_radius.csv     variance arms: scan radius by expression decile
#   power_paired_vs_ref.csv   paired per-round difference vs limma
#   power_fdr_tpr_curve.csv   TPR on an achieved-FDR grid
#   power_pr_curve.csv        PR curve (interpolated envelope, display only)
#   power_stratified.csv      TPR by |true log2 FC| bin and by expression decile
#   power_oracle_gap.csv      Part Q oracle - production, paired per round
#   power_pool_diag.csv       Part Q pool diagnostics, oracle - production, paired
#   power_breakdown.csv       Part C breakdown Delta (FDR@0.05 > 0.10), bootstrap CI
#   power_real_checks.csv     Part R per-cohort checks
#   rs_null_calibration.png, sample_size_curve.png, rs_paired_differences.png  (Part RS)
#   fdr_tpr_curves_vmatch_{vs_ref,all}_{Pplus,V,R,RS}.png
#                             variance-band arms' TPR-FDR curves vs limma/edgeR/DESeq2,
#                             and vs every method of the part
#   *.png
# =============================================================================

if (!exists("COMMON_DIR")) COMMON_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand <- c(if (length(f)) file.path(dirname(normalizePath(f[1])), "..", "..", "common"),
            "analysis/common", "experiments/Noise_Model_Test/common", "common", ".")
  hit <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("common/ not found (need config.R); looked in: ", paste(cand, collapse = ", "))
  normalizePath(hit[1])
})
.libPaths(c(normalizePath("external/docker_r_libs", mustWork = FALSE), .libPaths()))
source(file.path(COMMON_DIR, "config.R"))
source(file.path(COMMON_DIR, "power_metrics.R"))
source(file.path(COMMON_DIR, "print_utils.R"))
options(width = 250, max.print = 1e6)          # note C R1 (also set in config.R)
suppressPackageStartupMessages({ library(ggplot2); library(parallel) })
rd <- function(f) read.csv(f, stringsAsFactors = FALSE)      # .csv.gz decompressed by base R

.base <- if (exists("TOX_TEST_DIR")) dirname(TOX_TEST_DIR) else "."
RCFG <- list(
  raw_dir = Sys.getenv("POWER_RAW_DIR", file.path(.base, "power_raw")),
  out_dir = Sys.getenv("POWER_OUT_DIR", file.path(.base, "power_out")),
  n_cores = as.integer(Sys.getenv("POWER_CORES", "4")),
  paired_ref = "limma",
  fdr_limit  = 0.10,             # Part C breakdown: FDR at nominal 0.05 above this
  plot_methods = c("limma", "edgeR", "DESeq2", "TOX-log", "TOX-raw",
                   "TOX-log-blocked", "TOX-log-oracle", "TOX-log-blocked-oracle"),
  plot_methods_C = c("limma", "edgeR", "DESeq2", "TOX-log-TPM",
                     "TOX-log-TPM-centred", "TOX-log-TMM"),
  grid_inputs = c("TPM", "TPM-centred", "counts", "counts-centred", "TMM", "MoR"),
  grid_arms   = c("TOX-log", "TOX-log-blocked", "TOX-raw", "TOX-raw-blocked"),
  # Variance-matched arms (implementation notes A/B) and, per part, the plain TOX arm
  # each is paired against (same input; Part R runs them on TPM-centred).
  # vmatch_c: the band constant of each arm (power_test.R PCFG$vmatch_arms).
  vmatch_c = c(`TOX-log-vband-c1` = 1, `TOX-log-vband` = 1.5, `TOX-log-vband-c2` = 2,
               `TOX-log-vband-c3` = 3, `TOX-log-vband-c4` = 4,
               `TOX-log-vshrunk-c0.5` = 0.5, `TOX-log-vshrunk` = 1, `TOX-log-vshrunk-c2` = 2,
               `TOX-log-vshrunk-c3` = 3, `TOX-log-vshrunk-c4` = 4),
  vmatch_methods = c("TOX-log-vband-c1", "TOX-log-vband", "TOX-log-vband-c2", "TOX-log-vband-c3",
                     "TOX-log-vband-c4", "TOX-log-vshrunk-c0.5", "TOX-log-vshrunk", "TOX-log-vshrunk-c2",
                     "TOX-log-vshrunk-c3", "TOX-log-vshrunk-c4", "TOX-log-vshrunk-classicalvar"),
  vmatch_comparator = c(`P+` = "TOX-log", V = "TOX-log", R = "TOX-log|TPM-centred", RS = "TOX-log|TPM-centred"),
  # Note D WP4: B minus A at equal c (B c1 is paired with A c1.5, the A default).
  b_vs_a = list(c(B = "TOX-log-vshrunk",    A = "TOX-log-vband"),
                c(B = "TOX-log-vshrunk-c2", A = "TOX-log-vband-c2"),
                c(B = "TOX-log-vshrunk-c3", A = "TOX-log-vband-c3"),
                c(B = "TOX-log-vshrunk-c4", A = "TOX-log-vband-c4")),
  # Note D section 6 (confirmed by the project owner 21/10/2026). A cohort-count
  # condition "in at least 5 of 6 cohorts" is applied as n_cohorts - 1.
  gate = list(cal_margin = 0.01,      # step 1: p_le05_null <= comparator + this
              auprc_gain = 0.02,      # steps 3, 4: AUPRC excess gain needed
              auprc_loss = 0.01,      # step 3: no cohort loses more than this
              fdr_margin = 0.02,      # step 3: fdr_nom_05 <= comparator + this, every cohort
              b_cal_margin = 0.01)    # step 4: B's p_le05_null <= A's + this, every cohort
)
dir.create(RCFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# Categorical slots, fixed order (never cycled); a method keeps its colour across plots.
PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")
PART_LABEL <- c(S = "STRESS ARM (bimodal) -- not for ranking methods")
REF_METHODS <- c("limma", "edgeR", "DESeq2")

#' Composite encoding for many methods (plots 11-14): colour = family (fixed slots,
#' so the 8-slot palette is never exceeded), line type = member of the family. Both
#' scales are keyed by method, so ggplot merges them into one legend entry per method.
vm_family <- function(m) ifelse(m %in% REF_METHODS, "reference tools",
                         ifelse(grepl("vband", m), "A: vband",
                         ifelse(grepl("vshrunk", m), "B: vshrunk", sub("\\|.*$", "", m))))
vm_style <- function(m) {
  fam <- vm_family(m)
  lv  <- c("reference tools", "A: vband", "B: vshrunk",
           setdiff(sort(unique(fam)), c("reference tools", "A: vband", "B: vshrunk")))
  o <- order(match(fam, lv), match(m, RCFG$vmatch_methods), m); m <- m[o]; fam <- fam[o]
  slot <- ave(seq_along(m), fam, FUN = seq_along)
  list(methods = m, col = setNames(PALETTE[match(fam, lv)], m),
       lty = setNames(rep_len(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), max(slot))[slot], m))
}

# =============================================================================
# 1. METRICS FROM RAW
# =============================================================================

read_manifests <- function(dir) {
  f <- list.files(dir, "^manifest\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No manifest.csv under ", dir, " -- run power_test.R first.")
  m <- lapply(f, rd); cols <- unique(unlist(lapply(m, names)))
  m <- do.call(rbind, lapply(m, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
  m$part <- as.character(m$part); m$cohort <- as.character(m$cohort)
  m$delta_target <- as.numeric(m$delta_target)
  m
}

#' All metric tables for one job.
job_metrics <- function(job) {
  d <- file.path(RCFG$raw_dir, job$part)
  rf <- file.path(d, sprintf("results_%s.csv.gz", job$job_id))
  if (!file.exists(rf)) return(NULL)
  tr <- rd(file.path(d, sprintf("truth_%s.csv.gz", job$job_id)))
  rs <- rd(rf)
  tk <- tr[tr$kept, ]
  meta <- job[, setdiff(names(job), c("git_commit", "cohort_file")), drop = FALSE]
  n_de_total <- sum(tr$is_de)
  out <- list(metrics = list(), curve = list(), pr = list(), strata = list(), radius = list())
  for (m in unique(rs$method)) {
    r <- rs[rs$method == m, ]
    # Pairing check (note D WP3): every method must have scored exactly the kept genes.
    gene_set_ok <- setequal(r$gene, tk$gene) && !anyDuplicated(r$gene)
    r <- r[match(tk$gene, r$gene), ]
    add <- function(x) if (!is.null(x)) cbind(meta, method = m, x, row.names = NULL)
    out$metrics[[m]] <- add(cbind(metrics_one(r, tk$is_de, tk$true_lfc, n_de_total, tk$is_hv,
                                              lfc_log2 = identical(r$lfc_scale[1], "log2")),
                                  gene_set_ok = gene_set_ok))
    out$radius[[m]] <- add(radius_strata(r, tk$expr))
    out$curve[[m]]  <- add(curve_one(r, tk$is_de))
    out$pr[[m]]     <- add(curve_pr_one(r, tk$is_de))
    out$strata[[m]] <- add(strata_one(r, tk$is_de, tk$true_lfc, tk$expr))
  }
  lapply(out, function(x) do.call(rbind, x))
}

#' Part Q pool diagnostics: per job, pool x (DE / null) means of the width and
#' neighbour-distance measures, then oracle - production paired per job. Only jobs
#' whose validation gate passed (R port == Fortran) are used.
pool_diag <- function(man) {
  q <- man[man$part == "Q" & man$has_diag %in% TRUE, ]
  if (!nrow(q)) return(NULL)
  cols <- c("sd_pool_case", "sd_pool_control", "sd_null", "mean_span_case", "mean_span_control",
            "n_resid_pool_case", "n_resid_pool_control", "kurt_pool_case", "rho_case")
  rows <- lapply(seq_len(nrow(q)), function(i) {
    f <- file.path(RCFG$raw_dir, "Q", sprintf("diag_%s.csv.gz", q$job_id[i]))
    if (!file.exists(f)) return(NULL)
    d <- rd(f)
    if (!isTRUE(all(d$gate_ok))) return(data.frame(q[i, c(CELL, "round")], gate_ok = FALSE))
    a <- aggregate(d[cols], d[c("pool", "is_de")], mean, na.rm = TRUE)
    w <- merge(a[a$pool == "oracle", ], a[a$pool == "production", ], by = "is_de",
               suffixes = c(".orc", ".prod"))
    for (c in cols) w[[paste0("d_", c)]] <- w[[paste0(c, ".orc")]] - w[[paste0(c, ".prod")]]
    data.frame(q[i, c(CELL, "round")], gate_ok = TRUE,
               w[c("is_de", paste0(c(cols), ".prod"), paste0("d_", cols))], row.names = NULL)
  })
  X <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(X)) return(NULL)
  nfail <- sum(!X$gate_ok)
  if (nfail) message(sprintf("POOL DIAGNOSTICS: %d Part Q job(s) failed the Fortran validation gate -- excluded.", nfail))
  X <- X[X$gate_ok, ]
  if (!nrow(X)) return(NULL)
  mean_se(X, c(CELL, "is_de"), c(paste0(cols, ".prod"), paste0("d_", cols)))
}

#' Notes A/B section 7.3: paired per-round differences of every variance-band arm
#' against the plain TOX arm it replaces (RCFG$vmatch_comparator, same input). The
#' differences against limma are already in power_paired_vs_ref.csv.
paired_vs_tox <- function(M) {
  out <- lapply(names(RCFG$vmatch_comparator), function(pt) {
    comp <- RCFG$vmatch_comparator[[pt]]
    d <- M[M$part == pt & M$method %in% c(comp, RCFG$vmatch_methods), ]
    if (!any(d$method %in% RCFG$vmatch_methods)) return(NULL)
    x <- paired_vs_ref(d, comp, cols = c("fdr_nom_05", "p_le05_null", "any_call_05", "tpr_nom_05", "tpr_ach05",
                                         "auprc_excess_pe", "lfc80", "fpr_hv05"))
    if (!is.null(x)) x$comparator <- comp
    x
  })
  out <- Filter(Negate(is.null), out)
  if (length(out)) do.call(rbind, out)
}

#' Note D WP4: model B minus model A at equal c (RCFG$b_vs_a), paired per round,
#' with SE and t = mean / SE. Every part where both arms ran.
paired_b_vs_a <- function(M, cols = c("auprc_excess_pe", "tpr_nom_05", "fdr_nom_05", "p_le05_null", "any_call_05")) {
  out <- lapply(RCFG$b_vs_a, function(pr) {
    if (!all(pr %in% M$method)) return(NULL)
    X <- merge(M[M$method == pr[["B"]], c(CELL, "round", cols)], M[M$method == pr[["A"]], c(CELL, "round", cols)],
               by = c(CELL, "round"), suffixes = c("", ".A"))
    if (!nrow(X)) return(NULL)
    for (c in cols) X[[paste0("d_", c)]] <- X[[c]] - X[[paste0(c, ".A")]]
    X$pair <- sprintf("B c%g - A c%g", RCFG$vmatch_c[[pr[["B"]]]], RCFG$vmatch_c[[pr[["A"]]]])
    X$method <- pr[["B"]]; X$method_A <- pr[["A"]]
    a <- mean_se(X, c(CELL, "pair", "method", "method_A"), paste0("d_", cols))
    for (c in cols) a[[paste0("d_", c, "_t")]] <- a[[paste0("d_", c, "_mean")]] / a[[paste0("d_", c, "_se")]]
    a
  })
  out <- Filter(Negate(is.null), out)
  if (length(out)) do.call(rbind, out)
}

#' Part RS checks per (cohort, n): sampling, composition shift at small n (note D
#' WP3 VERIFY: pi_hat now comes from 2n samples), overlap between rounds.
rs_checks <- function(man) {
  r <- man[man$part == "RS", ]
  if (!nrow(r)) return(NULL)
  do.call(rbind, lapply(split(r, list(r$cohort, r$n_rep), drop = TRUE), function(d) {
    ok <- !d$skip %in% TRUE
    data.frame(cohort = d$cohort[1], n_rep = d$n_rep[1], S = d$S[ok][1], skipped = any(!ok),
               frac_samples_per_round = 2 * d$n_rep[1] / d$S[ok][1],
               rounds_pi0 = sum(ok & d$pi1 == 0), rounds_pi1 = sum(ok & d$pi1 > 0),
               distinct_datasets = length(unique(d$data_hash[ok])),
               n_kept_median = median(d$n_kept[ok]),
               delta_realised_median = median(d$delta_realised[ok & d$pi1 > 0], na.rm = TRUE),
               delta_realised_abs_q90 = unname(quantile(abs(d$delta_realised[ok & d$pi1 > 0]), 0.9, na.rm = TRUE)),
               n_delta_na = sum(ok & d$pi1 > 0 & is.na(d$delta)),
               frac_lfc_mismatch_mean = mean(d$frac_lfc_mismatch[ok], na.rm = TRUE),
               patient_check = d$patient_check[ok][1], n_dup_removed = d$n_dup_removed[ok][1])
  }))
}

#' Note D section 6, the confirmed decision protocol, on Part RS (per base cohort,
#' per n). Uses the per-cell means of power_summary.csv. Returns one table per step:
#'   1 calibration: arm passes in a cohort if p_le05_null <= comparator + cal_margin;
#'     the arm passes at n if it passes in >= n_cohorts - 1 cohorts
#'   2 c per family and n: the smallest c in the grid that passes step 1
#'   3 power at that c (pi1 = 0.1): improvement if d AUPRC excess >= auprc_gain in
#'     >= n_cohorts - 1 cohorts, no cohort below -auprc_loss, and fdr_nom_05 <=
#'     comparator + fdr_margin in every cohort; gaps to limma / edgeR / DESeq2 shown
#'   4 B (its c) vs A (its c): B better if d AUPRC excess >= auprc_gain in
#'     >= n_cohorts - 1 cohorts and B's p_le05_null <= A's + b_cal_margin in every cohort
rs_gates <- function(S) {
  S <- S[S$part == "RS", ]
  comp <- RCFG$vmatch_comparator[["RS"]]
  arms <- intersect(names(RCFG$vmatch_c), unique(S$method))
  if (!nrow(S) || !comp %in% S$method || !length(arms)) return(NULL)
  g <- RCFG$gate; fam <- ifelse(grepl("vband", arms), "A", "B"); names(fam) <- arms
  val <- function(pi0, m, col) {               # cohort-indexed vector of a mean
    x <- S[S$method == m & (if (pi0) S$pi1 == 0 else S$pi1 > 0), ]
    setNames(x[[paste0(col, "_mean")]], paste(x$cohort, x$n_rep))
  }
  key <- function(pi0, n) { x <- S[S$method == comp & (if (pi0) S$pi1 == 0 else S$pi1 > 0) & S$n_rep == n, ]
    paste(x$cohort, x$n_rep) }
  ns <- sort(unique(S$n_rep))
  # Step 1
  s1 <- do.call(rbind, lapply(ns, function(n) do.call(rbind, lapply(arms, function(a) {
    k <- key(TRUE, n); if (!length(k)) return(NULL)
    d <- val(TRUE, a, "p_le05_null")[k] - val(TRUE, comp, "p_le05_null")[k]
    need <- max(1L, length(k) - 1L)
    data.frame(step = 1L, n_rep = n, family = fam[[a]], c = RCFG$vmatch_c[[a]], method = a,
               n_cohorts = length(k), need = need, n_pass = sum(d <= g$cal_margin, na.rm = TRUE),
               d_p_le05_max = max(d, na.rm = TRUE), pass = sum(d <= g$cal_margin, na.rm = TRUE) >= need)
  }))))
  if (is.null(s1)) return(NULL)
  # Step 2
  s2 <- do.call(rbind, lapply(split(s1, list(s1$family, s1$n_rep), drop = TRUE), function(d) {
    ok <- d[d$pass, ]
    data.frame(step = 2L, n_rep = d$n_rep[1], family = d$family[1],
               c = if (nrow(ok)) min(ok$c) else NA_real_,
               method = if (nrow(ok)) ok$method[which.min(ok$c)] else NA_character_,
               pass = nrow(ok) > 0)
  }))
  # Step 3 (per cohort rows + one verdict row per family and n)
  s3 <- do.call(rbind, lapply(seq_len(nrow(s2)), function(i) {
    a <- s2$method[i]; n <- s2$n_rep[i]; if (is.na(a)) return(NULL)
    k <- key(FALSE, n); if (!length(k)) return(NULL)
    dA <- val(FALSE, a, "auprc_excess_pe")[k] - val(FALSE, comp, "auprc_excess_pe")[k]
    dF <- val(FALSE, a, "fdr_nom_05")[k] - val(FALSE, comp, "fdr_nom_05")[k]
    gap <- function(r) val(FALSE, a, "auprc_excess_pe")[k] - val(FALSE, r, "auprc_excess_pe")[k]
    need <- max(1L, length(k) - 1L)
    pass <- sum(dA >= g$auprc_gain, na.rm = TRUE) >= need && all(dA >= -g$auprc_loss, na.rm = TRUE) &&
            all(dF <= g$fdr_margin, na.rm = TRUE)
    data.frame(step = 3L, n_rep = n, family = s2$family[i], c = s2$c[i], method = a,
               cohort = c(sub(" [0-9]+$", "", k), "ALL"), need = need,
               d_auprc_excess = c(dA, NA), n_gain = c(rep(NA, length(k)), sum(dA >= g$auprc_gain, na.rm = TRUE)),
               d_fdr_nom_05 = c(dF, max(dF, na.rm = TRUE)),
               tpr_nom_05 = c(val(FALSE, a, "tpr_nom_05")[k], NA), tpr_ach05 = c(val(FALSE, a, "tpr_ach05")[k], NA),
               lfc80 = c(val(FALSE, a, "lfc80")[k], NA),
               gap_limma = c(gap("limma"), NA), gap_edgeR = c(gap("edgeR"), NA), gap_DESeq2 = c(gap("DESeq2"), NA),
               pass = c(rep(NA, length(k)), pass))
  }))
  # Step 4
  s4 <- do.call(rbind, lapply(ns, function(n) {
    a <- s2$method[s2$family == "A" & s2$n_rep == n]; b <- s2$method[s2$family == "B" & s2$n_rep == n]
    if (!length(a) || !length(b) || is.na(a) || is.na(b)) return(NULL)
    k <- key(FALSE, n); k0 <- key(TRUE, n); need <- max(1L, length(k) - 1L)
    dA <- val(FALSE, b, "auprc_excess_pe")[k] - val(FALSE, a, "auprc_excess_pe")[k]
    dC <- val(TRUE, b, "p_le05_null")[k0] - val(TRUE, a, "p_le05_null")[k0]
    data.frame(step = 4L, n_rep = n, method = b, method_A = a, need = need,
               n_gain = sum(dA >= g$auprc_gain, na.rm = TRUE), d_auprc_excess_mean = mean(dA, na.rm = TRUE),
               d_p_le05_null_max = max(dC, na.rm = TRUE),
               pass = sum(dA >= g$auprc_gain, na.rm = TRUE) >= need && all(dC <= g$b_cal_margin, na.rm = TRUE))
  }))
  list(s1 = s1, s2 = s2, s3 = s3, s4 = s4)
}

# =============================================================================
# 2. PLOTS
# =============================================================================

method_colours <- function(methods, ref = RCFG$plot_methods) {
  m <- methods[methods %in% ref[seq_len(min(8L, length(ref)))]]
  setNames(PALETTE[match(m, ref)], m)
}

zd_note <- function(M, part) {
  z <- M$zero_disc_05[M$part == part]
  sprintf("FDR zero-filled over rounds without calls (frac_zero_disc @0.05 = %.2f over all methods/cells)",
          mean(z, na.rm = TRUE))
}

make_plots <- function(M, CU, PR, ST, PT = NULL, BA = NULL) {
  dir <- RCFG$out_dir
  th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  save <- function(name, p, w = 10, h = 9) ggsave(file.path(dir, name), p, width = w, height = h, dpi = 150)
  sfx <- function(pt) if (pt == "S") "_S_stress" else ""
  ttl <- function(pt, x) if (pt %in% names(PART_LABEL)) paste0(PART_LABEL[[pt]], "\n", x) else x

  # (1)-(3) Parts P and S (S separately, labelled stress).
  for (pt in intersect(c("P", "S"), unique(M$part))) {
    cols <- method_colours(intersect(RCFG$plot_methods, unique(M$method[M$part == pt])))
    cs <- mean_se(CU[CU$part == pt & CU$method %in% names(cols), ], c(CELL, "method", "fdr"), "tpr")
    ms <- mean_se(M[M$part == pt & M$method %in% names(cols), ], c(CELL, "method"),
                  c(paste0("fdr_nom_", c("01", "05", "10")), paste0("tpr_nom_", c("01", "05", "10"))))
    pts <- do.call(rbind, lapply(c("01", "05", "10"), function(s)
      data.frame(ms[c(CELL, "method")], alpha = paste0("0.", s),
                 fdr = ms[[paste0("fdr_nom_", s, "_mean")]],
                 tpr = ms[[paste0("tpr_nom_", s, "_mean")]])))
    p <- ggplot(cs, aes(fdr, tpr_mean, colour = method)) +
      geom_vline(xintercept = c(0.01, 0.05, 0.10), colour = "grey85", linewidth = 0.3) +
      geom_line(linewidth = 0.6) +
      geom_point(data = pts[is.finite(pts$fdr), ], aes(fdr, tpr, shape = alpha), size = 2) +
      scale_colour_manual(values = cols) +
      facet_grid(dist ~ n_rep, labeller = label_both) +
      coord_cartesian(xlim = range(PM$fdr_grid)) +
      labs(title = ttl(pt, "TPR vs achieved FDR (lines); BH calls at nominal alpha (points)"),
           subtitle = paste("A point right of its alpha line = nominal FDR not honoured.", zd_note(M, pt)),
           x = "achieved FDR", y = "TPR (mean over rounds)", shape = "nominal alpha") + th
    save(sprintf("fdr_tpr_curves%s.png", sfx(pt)), p)

    for (typ in c("abs_lfc_bin", "expr_decile")) {
      sd_ <- ST[ST$part == pt & ST$stratum_type == typ & ST$method %in% names(cols), ]
      if (!nrow(sd_)) next
      ss <- mean_se(sd_, c(CELL, "method", "stratum"), "tpr_ach05")
      if (typ == "expr_decile") ss$stratum <- as.integer(ss$stratum)
      else ss$stratum <- factor(ss$stratum, levels = levels(cut(1, PM$lfc_bins, include.lowest = TRUE)))
      p <- ggplot(ss, aes(stratum, tpr_ach05_mean, colour = method, group = method)) +
        geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
        scale_colour_manual(values = cols) +
        facet_grid(dist ~ n_rep, labeller = label_both) +
        labs(title = ttl(pt, sprintf("Power at achieved FDR 0.05 by %s",
                             if (typ == "expr_decile") "control expression decile (1 = lowest)"
                             else "|true log2 FC|")),
             x = NULL, y = "TPR") + th
      save(sprintf("power_by_%s%s.png", typ, sfx(pt)), p)
    }

    # Full AUC vs pAUC(0.05) and vs AUPRC excess: which one separates the methods (P6)?
    a <- mean_se(M[M$part == pt & M$method %in% names(cols), ], c(CELL, "method"),
                 c("auc_p", "pauc05_p", "auprc_excess_pe"))
    la <- rbind(data.frame(a[c(CELL, "method", "auc_p_mean")], y_metric = "pAUC 0.05 (McClish)",
                           y = a$pauc05_p_mean),
                data.frame(a[c(CELL, "method", "auc_p_mean")], y_metric = "AUPRC excess (pe)",
                           y = a$auprc_excess_pe_mean))
    p <- ggplot(la, aes(auc_p_mean, y, colour = method)) +
      geom_point(size = 2) + scale_colour_manual(values = cols) +
      facet_grid(y_metric ~ dist + n_rep, labeller = label_both, scales = "free_y") +
      labs(title = ttl(pt, "Full ROC-AUC vs pAUC (FPR <= 0.05) and vs AUPRC - prevalence"),
           x = "AUC", y = NULL) + th
    save(sprintf("auc_vs_pauc%s.png", sfx(pt)), p, w = 12, h = 6)
  }

  # (4) PR curves: P, S (stress, separate), R (per cohort). Dashed = prevalence.
  for (pt in intersect(c("P", "S", "R"), unique(PR$part))) {
    cols <- method_colours(intersect(c(RCFG$plot_methods, if (pt == "R") paste0(RCFG$grid_arms, "|TPM")),
                                     unique(PR$method[PR$part == pt])),
                           c(RCFG$plot_methods[1:3], if (pt == "R") paste0(RCFG$grid_arms, "|TPM")
                                                     else RCFG$plot_methods[-(1:3)]))
    pr <- mean_se(PR[PR$part == pt & PR$method %in% names(cols), ], c(CELL, "method", "recall"), "precision")
    pv <- mean_se(M[M$part == pt, ], CELL, "prevalence")
    fac <- if (pt == "R") facet_wrap(~ cohort) else facet_grid(dist ~ n_rep, labeller = label_both)
    p <- ggplot(pr, aes(recall, precision_mean, colour = method)) +
      geom_hline(data = pv, aes(yintercept = prevalence_mean), linetype = "dashed", colour = "grey40") +
      geom_line(linewidth = 0.6) + scale_colour_manual(values = cols) + fac +
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = ttl(pt, sprintf("Part %s: precision-recall (interpolated envelope, display only)", pt)),
           subtitle = "Dashed = prevalence (random ranking). Scalar AUPRC uses the raw step estimator.",
           x = "recall", y = "precision (mean over rounds)") + th
    save(sprintf("pr_curves%s%s.png", if (pt == "R") "_R" else "", sfx(pt)), p)
  }

  # (5) Part Q: power and AUPRC excess vs pi1, by heterogeneity; oracle vs all.
  if (any(M$part == "Q")) {
    cols <- method_colours(intersect(RCFG$plot_methods, unique(M$method[M$part == "Q"])))
    q <- mean_se(M[M$part == "Q" & M$method %in% names(cols), ], c(CELL, "method"),
                 c("tpr_ach05", "tpr_nom_05", "auprc_excess_pe"))
    mk <- function(metric, col) data.frame(q[c("pi1", "het", "method")], metric = metric,
                                           v = q[[paste0(col, "_mean")]], se = q[[paste0(col, "_se")]])
    lq <- rbind(mk("TPR @ achieved FDR 0.05 (ranking)", "tpr_ach05"),
                mk("TPR @ nominal 0.05 (threshold)", "tpr_nom_05"),
                mk("AUPRC - prevalence (pe)", "auprc_excess_pe"))
    p <- ggplot(lq, aes(pi1, v, colour = method)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
      geom_errorbar(aes(ymin = v - se, ymax = v + se), width = 0.05, linewidth = 0.4) +
      scale_x_log10() + scale_colour_manual(values = cols) +
      facet_grid(metric ~ het, labeller = labeller(het = label_both), scales = "free_y") +
      labs(title = "Power vs DE fraction: het = 1 homogeneous effects, het < 1 responder fraction",
           subtitle = "*-oracle pools true nulls only (same Fortran); its gap to the production arm is pool contamination",
           x = "pi1 (fraction DE, log scale)", y = "mean +/- SE") + th
    save("power_vs_pi1.png", p, h = 10)
  }

  # (6) Part C: FDR and power vs TARGET Delta.
  if (any(M$part == "C")) {
    cols <- method_colours(intersect(RCFG$plot_methods_C, unique(M$method[M$part == "C"])), RCFG$plot_methods_C)
    cc <- mean_se(M[M$part == "C" & M$method %in% names(cols), ], c(CELL, "method"),
                  c("fdr_nom_05", "tpr_ach05"))
    lg <- rbind(data.frame(cc[c(CELL, "method")], metric = "observed FDR @ nominal 0.05",
                           v = cc$fdr_nom_05_mean, se = cc$fdr_nom_05_se),
                data.frame(cc[c(CELL, "method")], metric = "TPR @ achieved FDR 0.05",
                           v = cc$tpr_ach05_mean, se = cc$tpr_ach05_se))
    hl <- data.frame(metric = "observed FDR @ nominal 0.05", y = c(0.05, RCFG$fdr_limit),
                     lt = c("solid", "dashed"))
    p <- ggplot(lg, aes(delta_target, v, colour = method)) +
      geom_hline(data = hl, aes(yintercept = y, linetype = lt), colour = "grey55") +
      scale_linetype_identity() +
      geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
      geom_errorbar(aes(ymin = v - se, ymax = v + se), width = 0.03, linewidth = 0.4) +
      scale_colour_manual(values = cols) +
      facet_grid(metric ~ dist + n_rep, labeller = label_both, scales = "free_y") +
      labs(title = "Composition stress: every null gene shifts by -Delta in TPM (Delta targeted)",
           subtitle = paste("Scored against the biological truth beta. Dashed = breakdown limit.", zd_note(M, "C")),
           x = "target Delta = log2 sum(pi_a 2^beta)", y = NULL) + th
    save("composition_stress.png", p, w = 12, h = 6)
  }

  # (7) P+: calibration (null arm) and power vs sigma_d.
  if (any(M$part == "P+")) {
    cols <- method_colours(intersect(RCFG$plot_methods, unique(M$method[M$part == "P+"])))
    g <- mean_se(M[M$part == "P+" & M$method %in% names(cols), ], c(CELL, "method"),
                 c("p_le05_null", "fdr_nom_05", "tpr_ach05"))
    mk <- function(sel, metric, col) data.frame(g[sel, c(CELL, "method")], metric = metric,
                                                v = g[sel, paste0(col, "_mean")], se = g[sel, paste0(col, "_se")])
    lg <- rbind(mk(g$pi1 == 0, "null arm: P(p <= 0.05)", "p_le05_null"),
                mk(g$pi1 > 0, "power arm: FDR @ nominal 0.05", "fdr_nom_05"),
                mk(g$pi1 > 0, "power arm: TPR @ achieved FDR 0.05", "tpr_ach05"))
    p <- ggplot(lg, aes(sigma_d, v, colour = method)) +
      geom_hline(data = data.frame(metric = intersect(c("null arm: P(p <= 0.05)", "power arm: FDR @ nominal 0.05"),
                                                    lg$metric), y = 0.05),
                 aes(yintercept = y), colour = "grey55") +
      geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
      geom_errorbar(aes(ymin = v - se, ymax = v + se), width = 0.03, linewidth = 0.4) +
      scale_colour_manual(values = cols) +
      facet_grid(metric ~ dist + n_rep, labeller = label_both, scales = "free_y") +
      labs(title = "P+: gene-wise dispersion scatter around a mean-dependent trend",
           subtitle = sprintf("phi_g = phi_trend(mu) exp(sigma_d z - sigma_d^2/2); dispersion source: %s. %s",
                              paste(unique(M$disp_source[M$part == "P+"]), collapse = ","), zd_note(M, "P+")),
           x = "sigma_d", y = NULL) + th
    save("dispersion_heterogeneity.png", p, w = 12, h = 8)
  }

  # (8) Parts N / R: input x arm grid; reference tools as lines.
  arm_cols <- setNames(PALETTE[seq_along(RCFG$grid_arms)], RCFG$grid_arms)
  refs <- REF_METHODS
  mets <- c(fdr_nom_05 = "observed FDR @ nominal 0.05", tpr_nom_05 = "TPR @ nominal 0.05",
            tpr_ach05 = "TPR @ achieved FDR 0.05")
  for (pt in intersect(c("N", "R"), unique(M$part))) {
    g <- mean_se(M[M$part == pt, ], c(CELL, "method"), names(mets))
    g$arm   <- sub("\\|.*$", "", g$method)
    g$input <- ifelse(grepl("\\|", g$method), sub("^.*\\|", "", g$method), NA)
    panels <- if (pt == "N") split(g, g$n_rep) else list(all = g)
    for (nm in names(panels)) for (mt in names(mets)) {
      gg <- panels[[nm]]
      gi <- gg[!is.na(gg$input), ]; gr <- gg[gg$method %in% refs, ]
      if (!nrow(gi)) next
      gi$input <- factor(gi$input, levels = RCFG$grid_inputs)
      y <- paste0(mt, "_mean")
      fac <- if (pt == "N") facet_grid(disp_model ~ depth, labeller = label_both)
             else facet_wrap(~ cohort + n_rep, labeller = label_both)
      p <- ggplot(gi, aes(input, .data[[y]], colour = arm, group = arm)) +
        geom_hline(data = gr, aes(yintercept = .data[[y]], linetype = method),
                   colour = "grey45", linewidth = 0.4) +
        { if (mt == "fdr_nom_05") geom_hline(yintercept = 0.05, colour = "black", linewidth = 0.3) } +
        geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
        scale_colour_manual(values = arm_cols) + fac +
        labs(title = sprintf("Part %s%s: %s by input", pt,
                             if (pt == "N") paste0(", n_rep = ", nm) else " (whole cohorts, random halves)",
                             mets[[mt]]),
             subtitle = "Grey lines = reference tools on raw counts. Centred inputs are log arms only.",
             x = NULL, y = mets[[mt]], linetype = NULL) +
        th + theme(axis.text.x = element_text(angle = 35, hjust = 1))
      save(if (pt == "N") sprintf("grid_N_%s_n%s.png", mt, nm) else sprintf("grid_R_%s.png", mt), p,
           w = 11, h = if (pt == "N") 9 else 7)
    }
  }

  # (9) Part V: calls on variance-only nulls vs on ordinary nulls.
  if (any(M$part == "V")) {
    v <- mean_se(M[M$part == "V", ], c(CELL, "method"), c("fpr_hv05", "fpr_nonhv05"))
    lv <- rbind(data.frame(v[c("n_rep", "het", "method")], nulls = "variance-only (case dispersion x hv_factor)",
                           f = v$fpr_hv05_mean, se = v$fpr_hv05_se),
                data.frame(v[c("n_rep", "het", "method")], nulls = "ordinary",
                           f = v$fpr_nonhv05_mean, se = v$fpr_nonhv05_se))
    p <- ggplot(lv, aes(method, f, colour = nulls)) +
      geom_pointrange(aes(ymin = f - se, ymax = f + se), position = position_dodge(0.5), size = 0.3) +
      scale_colour_manual(values = PALETTE[1:2]) +
      facet_grid(het ~ n_rep, labeller = label_both) +
      labs(title = "Part V: fraction of null genes called at nominal 0.05",
           subtitle = "Variance-only nulls have no mean change -- every call on them is a false positive",
           x = NULL, y = "fraction called", colour = NULL) +
      th + theme(axis.text.x = element_text(angle = 35, hjust = 1))
    save("variance_only_nulls.png", p, w = 11, h = 6)
  }

  # (10) Variance-matched arms (notes A/B) against limma and their plain comparator.
  # Methods on x, colour = family, so the 8-slot palette is never exhausted.
  vm_fam <- function(m) ifelse(m %in% c("limma", RCFG$vmatch_comparator), "reference",
                        ifelse(grepl("vband", m), "A: vband", "B: vshrunk"))
  vm_mets <- list(`P+` = c(p_le05_null = "null arm: P(p <= 0.05)", fdr_nom_05 = "FDR @ nominal 0.05",
                           tpr_ach05 = "TPR @ achieved FDR 0.05", auprc_excess_pe = "AUPRC excess"),
                  V    = c(fpr_hv05 = "FPR hypervariable nulls", fpr_nonhv05 = "FPR ordinary nulls",
                           fdr_nom_05 = "FDR @ nominal 0.05", tpr_ach05 = "TPR @ achieved FDR 0.05"),
                  R    = c(fdr_nom_05 = "FDR @ nominal 0.05", tpr_ach05 = "TPR @ achieved FDR 0.05",
                           auprc_excess_pe = "AUPRC excess", frac_short_pool = "frac short pool"),
                  RS   = c(p_le05_null = "null arm: P(p <= 0.05)", any_call_05 = "null arm: P(any BH call)",
                           fdr_nom_05 = "FDR @ nominal 0.05", auprc_excess_pe = "AUPRC excess"))
  vm_keys <- list(`P+` = c("dist", "n_rep", "sigma_d"), V = c("n_rep", "het"), R = c("cohort", "n_rep"),
                  RS = c("cohort", "n_rep"))
  null_mets <- c("p_le05_null", "any_call_05")
  for (pt in intersect(names(vm_mets), unique(M$part))) {
    meth <- c("limma", RCFG$vmatch_comparator[[pt]], RCFG$vmatch_methods)
    d <- M[M$part == pt & M$method %in% meth, ]
    if (!any(d$method %in% RCFG$vmatch_methods)) next
    mets <- vm_mets[[pt]]
    lg <- do.call(rbind, lapply(names(mets), function(mt) {
      dd <- if (mt %in% null_mets) d[d$pi1 == 0, ] else if (pt %in% c("P+", "RS")) d[d$pi1 > 0, ] else d
      if (!nrow(dd)) return(NULL)
      a <- mean_se(dd, c(vm_keys[[pt]], "method"), mt)
      data.frame(a[c(vm_keys[[pt]], "method")], metric = mets[[mt]],
                 v = a[[paste0(mt, "_mean")]], se = a[[paste0(mt, "_se")]])
    }))
    if (is.null(lg) || !nrow(lg)) next
    lg$method <- factor(lg$method, levels = meth)
    lg$family <- vm_fam(as.character(lg$method))
    lg$cell <- do.call(paste, c(lapply(vm_keys[[pt]], function(k) paste0(k, "=", lg[[k]])), sep = ", "))
    p <- ggplot(lg, aes(method, v, colour = family)) +
      geom_pointrange(aes(ymin = v - se, ymax = v + se), size = 0.25) +
      { if (any(grepl("FDR|P\\(p|FPR|any BH", lg$metric)))
          geom_hline(data = data.frame(metric = unique(grep("FDR|P\\(p|FPR|any BH", lg$metric, value = TRUE)), y = 0.05),
                     aes(yintercept = y), colour = "grey55") } +
      scale_colour_manual(values = c(reference = PALETTE[1], `A: vband` = PALETTE[2], `B: vshrunk` = PALETTE[3])) +
      facet_grid(metric ~ cell, scales = "free_y") +
      labs(title = sprintf("Part %s: variance-matched neighbourhood (notes A/B) vs %s and limma",
                           pt, RCFG$vmatch_comparator[[pt]]),
           subtitle = "c_tau is chosen on calibration only (null halves, P+ pi1 = 0); paired differences in power_paired_vs_tox.csv",
           x = NULL, y = NULL, colour = NULL) +
      th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save(sprintf("variance_match_%s.png", gsub("[^A-Za-z]", "plus", pt)), p,
         w = max(10, 1.6 * length(unique(lg$cell)) + 4), h = 2.2 * length(mets) + 2)
  }

  # (11) Variance-matched arms: TPR vs achieved FDR curves as in (1), per part. One
  # figure against the reference tools only, one with every method of the part.
  # Composite encoding: vm_style() (colour = family, line type = member).
  vm_fac <- list(`P+` = facet_grid(dist + n_rep ~ sigma_d, labeller = label_both),
                 V    = facet_grid(het ~ n_rep, labeller = label_both),
                 R    = facet_wrap(~ cohort + n_rep, labeller = label_both),
                 RS   = facet_grid(n_rep ~ cohort, labeller = label_both))
  vm_size <- list(`P+` = c(12, 15), V = c(10, 7), R = c(12, 8), RS = c(16, 7))
  for (pt in intersect(names(vm_fac), unique(CU$part))) {
    all_m <- unique(CU$method[CU$part == pt])
    if (!any(all_m %in% RCFG$vmatch_methods)) next
    for (set in c("vs_ref", "all")) {
      st <- vm_style(if (set == "vs_ref") intersect(all_m, c(REF_METHODS, RCFG$vmatch_methods)) else all_m)
      sel <- function(X) X[X$part == pt & X$method %in% st$methods & (is.na(X$pi1) | X$pi1 > 0), ]
      cs <- mean_se(sel(CU), c(CELL, "method", "fdr"), "tpr")
      ms <- mean_se(sel(M), c(CELL, "method"), c("fdr_nom_05", "tpr_nom_05"))
      cs$method <- factor(cs$method, levels = st$methods)
      ms$method <- factor(ms$method, levels = st$methods)
      p <- ggplot(cs, aes(fdr, tpr_mean, colour = method, linetype = method)) +
        geom_vline(xintercept = 0.05, colour = "grey85", linewidth = 0.3) +
        geom_line(linewidth = 0.5) +
        geom_point(data = ms[is.finite(ms$fdr_nom_05_mean), ], aes(fdr_nom_05_mean, tpr_nom_05_mean),
                   size = 1.6, show.legend = FALSE) +
        scale_colour_manual(values = st$col, name = NULL) +
        scale_linetype_manual(values = st$lty, name = NULL) +
        vm_fac[[pt]] + coord_cartesian(xlim = range(PM$fdr_grid)) +
        labs(title = sprintf("Part %s: TPR vs achieved FDR, variance-matched arms vs %s", pt,
                             if (set == "vs_ref") "reference tools" else "every method"),
             subtitle = paste("Points = BH calls at nominal 0.05 (right of the grey line = FDR not honoured).",
                              if (pt %in% c("P+", "RS")) "Null arm (pi1 = 0) has no TPR and is left out." else ""),
             x = "achieved FDR", y = "TPR (mean over rounds)") +
        th + theme(legend.key.width = unit(1.8, "lines"))
      save(sprintf("fdr_tpr_curves_vmatch_%s_%s.png", set, gsub("[^A-Za-z]", "plus", pt)), p,
           w = vm_size[[pt]][1], h = vm_size[[pt]][2] + if (set == "all") 2 else 0)
    }
  }

  if (!any(M$part == "RS")) return(invisible(NULL))
  rs_m <- intersect(c(REF_METHODS, RCFG$vmatch_comparator[["RS"]], RCFG$vmatch_methods), unique(M$method[M$part == "RS"]))

  # (12) Part RS null calibration (note D WP6a): P(p <= 0.05) over the null genes of
  # the pi1 = 0 rounds, mean +/- SE, per method and n, faceted by base cohort.
  nc <- mean_se(M[M$part == "RS" & M$pi1 == 0 & M$method %in% rs_m, ], c("cohort", "n_rep", "method"), "p_le05_null")
  if (!is.null(nc) && nrow(nc)) {
    nc$method <- factor(nc$method, levels = rs_m); nc$n <- factor(nc$n_rep)
    p <- ggplot(nc, aes(method, p_le05_null_mean, colour = n)) +
      geom_hline(yintercept = 0.05, colour = "grey55") +
      geom_pointrange(aes(ymin = p_le05_null_mean - p_le05_null_se, ymax = p_le05_null_mean + p_le05_null_se),
                      position = position_dodge(0.5), size = 0.25) +
      scale_colour_manual(values = PALETTE[1:2]) + facet_wrap(~ cohort) +
      labs(title = "Part RS: null calibration on real data at small n (pi1 = 0)",
           subtitle = "Fraction of null genes with raw p <= 0.05, mean +/- SE over rounds. Line = 0.05.",
           x = NULL, y = "P(p <= 0.05 | null)", colour = "n per group") +
      th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save("rs_null_calibration.png", p, w = 14, h = 8)
  }

  # (13) Sample-size curve (note D WP6b): RS (n = 3, 5) and R (n = S %/% 2) of the
  # same base cohort, power arm (pi1 = 0.1), methods run in both parts.
  both <- intersect(rs_m, unique(M$method[M$part == "R"]))
  sz <- M[M$part %in% c("R", "RS") & M$pi1 > 0 & M$method %in% both, ]
  if (length(both) && length(unique(sz$part)) == 2L) {
    st <- vm_style(both)
    mets <- c(auprc_excess_pe = "AUPRC excess", tpr_nom_05 = "TPR @ nominal 0.05", fdr_nom_05 = "FDR @ nominal 0.05")
    a <- mean_se(sz, c("cohort", "n_rep", "method"), names(mets))
    lg <- do.call(rbind, lapply(names(mets), function(mt)
      data.frame(a[c("cohort", "n_rep", "method")], metric = mets[[mt]], v = a[[paste0(mt, "_mean")]])))
    lg$method <- factor(lg$method, levels = st$methods)
    p <- ggplot(lg, aes(n_rep, v, colour = method, linetype = method)) +
      geom_hline(data = data.frame(metric = mets[["fdr_nom_05"]], y = 0.05), aes(yintercept = y), colour = "grey55") +
      geom_line(linewidth = 0.5) + geom_point(size = 1.4) +
      scale_x_log10() + scale_colour_manual(values = st$col, name = NULL) +
      scale_linetype_manual(values = st$lty, name = NULL) +
      facet_grid(metric ~ cohort, scales = "free_y") +
      labs(title = "Power vs sample size: Part RS (n = 3, 5) and Part R (whole-cohort halves), pi1 = 0.1",
           subtitle = "Same base cohort per column; R points sit at n = floor(S / 2).",
           x = "n per group (log scale)", y = NULL) +
      th + theme(legend.key.width = unit(1.8, "lines"))
    save("sample_size_curve.png", p, w = 16, h = 9)
  }

  # (14) Paired differences with SE (note D WP6c): every variance arm minus
  # TOX-log|TPM-centred, and B minus A at equal c. Calibration metrics on pi1 = 0,
  # power metrics on pi1 = 0.1.
  pd_mets <- c(d_auprc_excess_pe = "AUPRC excess", d_tpr_nom_05 = "TPR @ nominal 0.05",
               d_fdr_nom_05 = "FDR @ nominal 0.05", d_p_le05_null = "P(p <= 0.05 | null)")
  pd <- function(X, lab) {
    if (is.null(X)) return(NULL)
    X <- X[X$part == "RS", ]; if (!nrow(X)) return(NULL)
    do.call(rbind, lapply(names(pd_mets), function(mt) {
      if (!paste0(mt, "_mean") %in% names(X)) return(NULL)
      x <- X[if (mt == "d_p_le05_null") X$pi1 == 0 else X$pi1 > 0, ]
      if (!nrow(x)) return(NULL)
      data.frame(cohort = x$cohort, n = factor(x$n_rep), arm = lab(x), metric = pd_mets[[mt]],
                 v = x[[paste0(mt, "_mean")]], se = x[[paste0(mt, "_se")]])
    }))
  }
  lg <- rbind(pd(PT, function(x) paste(x$method, "- TOX")), pd(BA, function(x) x$pair))
  if (!is.null(lg) && nrow(lg)) {
    lg$arm <- factor(lg$arm, levels = unique(lg$arm))
    p <- ggplot(lg, aes(arm, v, colour = n)) +
      geom_hline(yintercept = 0, colour = "grey55") +
      geom_pointrange(aes(ymin = v - se, ymax = v + se), position = position_dodge(0.5), size = 0.2) +
      scale_colour_manual(values = PALETTE[1:2]) +
      facet_grid(metric ~ cohort, scales = "free_y") +
      labs(title = "Part RS: paired per-round differences, mean +/- SE",
           subtitle = "'<arm> - TOX' = variance arm minus TOX-log|TPM-centred; 'B .. - A ..' = model B minus model A at equal c.",
           x = NULL, y = NULL, colour = "n per group") +
      th + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
    save("rs_paired_differences.png", p, w = 18, h = 10)
  }
}

# =============================================================================
# 3. CONSOLE REPORT
# =============================================================================

spread <- function(x) { x <- x[is.finite(x)]; if (length(x)) diff(range(x)) else NA_real_ }

#' One screen on what this report is built from (note C section 3.3).
report_header <- function(man) {
  ok <- man[!man$skip %in% TRUE, ]
  rr <- aggregate(round ~ part, ok, max)
  cat(strrep("=", 100), "\npower_report.R  ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("raw: %s | out: %s | jobs: %d | parts: %s\n", RCFG$raw_dir, RCFG$out_dir, nrow(ok),
              paste(unique(ok$part), collapse = ",")))
  cat("rounds per cell: ", paste(rr$part, rr$round, sep = "=", collapse = ", "), "\n", sep = "")
  cat("commit: ", paste(unique(man$git_commit), collapse = " / "), "\n", sep = "")
  pp <- ok[ok$part == "P+", ]
  if (nrow(pp)) {
    f <- function(col) if (col %in% names(pp)) paste(unique(format(pp[[col]])), collapse = "/") else "-"
    cat(sprintf("P+ dispersion parameters: %s | date %s | phi_inf %s | c %s | sigma_hat %s | mu_depth %s\n",
                f("disp_source"), f("disp_date"), f("disp_phi_inf"), f("disp_c"), f("disp_sigma_hat"), f("disp_mu_depth")))
  }
  co <- unique(ok$cohort[!is.na(ok$cohort)])
  if (length(co)) cat(sprintf("base cohorts: %d (%s)\n", length(co), paste(co, collapse = ", ")))
  if (any(man$skip %in% TRUE))
    cat("skipped RS cells (2n > S): ", paste(unique(paste0(man$cohort, " n=", man$n_rep)[man$skip %in% TRUE]),
                                             collapse = "; "), "\n", sep = "")
  cat(strrep("=", 100), "\n", sep = "")
}

#' Headline pivots of one part (note C R6): one metric per table, methods side by
#' side. Parts with a null arm (pi1 = 0) get the calibration set on it and the
#' power set on pi1 > 0; the variance arms add their pool diagnostics.
report_part <- function(S, pt) {
  x <- S[S$part == pt, ]
  meth <- vm_style(unique(x$method))$methods
  lab <- if (pt %in% names(PART_LABEL)) paste0(" -- ", PART_LABEL[[pt]]) else ""
  power <- c("tpr_nom_05", "fdr_nom_05", "tpr_ach05", "auprc_excess_pe", if (pt == "V") c("fpr_hv05", "fpr_nonhv05"))
  null  <- c("p_le05_null", "any_call_05", "ks_null")
  has0 <- any(x$pi1 == 0)
  pv <- function(d, mt, arm) print_pivot(d, paste0(mt, "_mean"), methods = meth, csv = "power_summary.csv",
                                         title = sprintf("PART %s%s | %s%s", pt, arm, mt, lab))
  if (has0) for (mt in null) pv(x[x$pi1 == 0, ], mt, " null arm (pi1 = 0)")
  if (any(x$pi1 > 0 | is.na(x$pi1))) for (mt in power) pv(x[!x$pi1 %in% 0, ], mt, if (has0) " power arm" else "")
  vm <- x[x$method %in% RCFG$vmatch_methods, ]
  if (nrow(vm)) for (mt in c("frac_short_pool", "n_scanned_median", "scan_radius_median", "frac_radius_gt2"))
    if (any(is.finite(vm[[paste0(mt, "_mean")]])))
      print_pivot(vm[!vm$pi1 %in% 0 | !has0, ], paste0(mt, "_mean"), methods = meth, csv = "power_summary.csv",
                  title = sprintf("PART %s | variance arms: %s", pt, mt))
}

report <- function(S, PV, OG, PD, BK, RC, PT, BA, RSC, GT, M) {
  for (pt in intersect(ALL_PART_ORDER, unique(S$part))) report_part(S, pt)
  cat("\n  tpr_ach05   TPR at ACHIEVED FDR 0.05 (truth-thresholded, max over steps; not user-realisable).\n",
      "  tpr_nom_05 / fdr_nom_05  what a user gets at nominal 0.05; fdr zero-filled.\n",
      "  any_call_05  P(at least one BH call); under pi1 = 0 it is P(any false call) = the FDR.\n",
      "  auprc_excess_pe  AUPRC minus prevalence: comparable across pi1.\n",
      "  ks_null  Kolmogorov-Smirnov distance of the null p-values to U(0, 1).\n", sep = "")

  bad <- M[!M$gene_set_ok %in% TRUE & "gene_set_ok" %in% names(M), c("part", "job_id", "method")]
  if (nrow(bad)) print_tab(bad, "PAIRING CHECK FAILED: methods that did not score exactly the kept genes",
                           csv = "power_per_round.csv", keys = c("part", "job_id", "method"))

  sp <- do.call(rbind, lapply(split(S, do.call(paste, c(S[CELL], sep = "\r"))), function(d)
    data.frame(d[1, CELL], auc_range = spread(d$auc_p_mean), pauc05_range = spread(d$pauc05_p_mean),
               auprc_range = spread(d$auprc_pe_mean))))
  print_tab(sp, "SPREAD ACROSS METHODS per cell (P6: pauc05_range / auprc_range > auc_range?)")
  if (!is.null(PV)) print_tab(PV, paste0("PAIRED DIFFERENCE vs ", RCFG$paired_ref, " (mean +/- SE over rounds)"),
                              csv = "power_paired_vs_ref.csv")
  if (!is.null(PT)) print_tab(PT[, c(CELL, "method", "comparator", grep("^d_.*_mean$", names(PT), value = TRUE))],
                              "VARIANCE ARMS minus their plain TOX comparator (paired; SE in the CSV)",
                              csv = "power_paired_vs_tox.csv", keys = c(PRINT_KEYS, "comparator"))
  if (!is.null(BA)) print_tab(BA[, c(CELL, "pair", grep("^d_.*_(mean|t)$", names(BA), value = TRUE))],
                              "MODEL B minus MODEL A at equal c (paired; t = mean / SE)",
                              csv = "power_paired_B_vs_A.csv", keys = c(PRINT_KEYS, "pair"))
  if (!is.null(OG)) print_tab(OG, "ORACLE GAP (Part Q): oracle minus production, paired per round", csv = "power_oracle_gap.csv",
    notes = c("P3': contamination is a CALIBRATION effect -> look at d_fdr_nom_05, d_fpr_null05, d_n_called_05,",
              "d_tpr_nom_05. d_tpr_ach05 is the expected-null control (invariant to monotone p inflation)."))
  if (!is.null(PD)) print_tab(PD, "POOL DIAGNOSTICS (Part Q, exact TOX-log): production values and oracle - production",
    csv = "power_pool_diag.csv", keys = c(PRINT_KEYS, "is_de"),
    notes = c("d_sd_pool_* < 0 = contamination widened the production pool.",
              "d_mean_span_* > 0 = the oracle's neighbours lie further away in mean (confound, not contamination)."))
  if (!is.null(BK)) print_tab(BK, sprintf("PART C BREAKDOWN: smallest target Delta with FDR@0.05 > %.2f (bootstrap CI)",
                                          RCFG$fdr_limit), csv = "power_breakdown.csv",
                              notes = "Inf = no breakdown inside the Delta grid.")
  if (!is.null(RC)) print_tab(RC, "PART R CHECKS per base cohort", csv = "power_real_checks.csv",
    notes = c("frac_len_bad: genes dropped because count/TPM was not constant (length recovery invalid).",
              "frac_lfc_mismatch: injected genes with |realised - intended log2 FC| > lfc_tol.",
              "n_delta_na: rounds where the thinning balance failed."))
  if (!is.null(RSC)) print_tab(RSC, "PART RS CHECKS per base cohort and n", csv = "power_rs_checks.csv",
    keys = c("cohort", "n_rep"),
    notes = c("delta_realised: composition shift of the nulls (target 0); pi_hat now comes from 2n samples.",
              "frac_samples_per_round = 2n / S: rounds overlap, so SEs are conditional on the base cohort.",
              "distinct_datasets should equal rounds_pi0 + rounds_pi1 (data_hash)."))
  if (!is.null(GT)) {
    nc <- length(unique(S$cohort[S$part == "RS"]))
    cat("\n", strrep("#", 100), "\nPART RS DECISION PROTOCOL (note D section 6; ", nc,
        " base cohorts -- every statement below is per base cohort, SEs are never pooled)\n", sep = "")
    g <- RCFG$gate
    print_tab(GT$s1, "STEP 1  calibration (pi1 = 0): pass in a cohort if p_le05_null <= comparator + cal_margin",
              csv = "power_rs_gates.csv", keys = c("n_rep", "family", "c", "method"),
              notes = sprintf("cal_margin = %g; an arm passes at n if it passes in >= need cohorts", g$cal_margin))
    print_tab(GT$s2, "STEP 2  c per family and n: smallest c that passes step 1 (NA = fails calibration at n)",
              csv = "power_rs_gates.csv", keys = c("n_rep", "family"))
    if (!is.null(GT$s3)) print_tab(GT$s3, "STEP 3  power at the chosen c (pi1 = 0.1), per cohort; ALL = verdict",
              csv = "power_rs_gates.csv", keys = c("n_rep", "family", "method", "cohort"),
              notes = sprintf("improvement: d_auprc_excess >= %g in >= need cohorts, none < -%g, d_fdr_nom_05 <= %g in every cohort; gap_* = arm minus reference",
                              g$auprc_gain, g$auprc_loss, g$fdr_margin))
    if (!is.null(GT$s4)) print_tab(GT$s4, "STEP 4  does shrinkage help? B (its c) vs A (its c)",
              csv = "power_rs_gates.csv", keys = c("n_rep", "method"),
              notes = sprintf("B better: d_auprc_excess >= %g in >= need cohorts and B's p_le05_null <= A's + %g in every cohort",
                              g$auprc_gain, g$b_cal_margin))
  }
}

ALL_PART_ORDER <- c("R", "RS", "P", "P+", "Q", "C", "S", "N", "V", "D")

# =============================================================================
if (sys.nframe() == 0L) {
  man <- read_manifests(RCFG$raw_dir)
  report_header(man)
  run <- man[!man$skip %in% TRUE, ]
  message(sprintf("%d jobs in %s; cores %d", nrow(run), RCFG$raw_dir, RCFG$n_cores))
  res <- mclapply(split(run, seq_len(nrow(run))), function(j)
           tryCatch(job_metrics(j), error = function(e) {
             message("job ", j$job_id, ": ", conditionMessage(e)); NULL }),
         mc.cores = RCFG$n_cores, mc.preschedule = TRUE)
  res <- Filter(function(x) is.list(x) && !is.null(x$metrics), res)
  if (!length(res)) stop("No job could be scored.")
  M   <- do.call(rbind, lapply(res, `[[`, "metrics"))
  CU  <- do.call(rbind, lapply(res, `[[`, "curve"))
  PR  <- do.call(rbind, lapply(res, `[[`, "pr"))
  ST  <- do.call(rbind, lapply(res, `[[`, "strata"))
  RAD <- do.call(rbind, lapply(res, `[[`, "radius"))

  skip <- c(CELL, "round", "job", "seed", "split_seed", "sigma_mult")
  num_cols <- setdiff(names(M)[vapply(M, is.numeric, logical(1))], skip)
  S  <- mean_se(M, c(CELL, "method"), num_cols)
  PV <- paired_vs_ref(M, RCFG$paired_ref)
  OG <- oracle_gap(M)
  PD <- tryCatch(pool_diag(man), error = function(e) { message("pool diagnostics: ", conditionMessage(e)); NULL })
  BK <- if (any(M$part == "C")) {
    mc <- M[M$part == "C", ]
    fam <- do.call(paste, c(mc[c("dist", "n_rep", "method")], sep = "\r"))
    set.seed(1L)
    do.call(rbind, lapply(split(mc, fam), function(d)
      cbind(d[1, c("dist", "n_rep", "method")], breakdown_point(d, RCFG$fdr_limit))))
  }
  RC <- if (any(run$part == "R")) {
    r <- run[run$part == "R", ]
    do.call(rbind, lapply(split(r, r$cohort), function(d) data.frame(
      cohort = d$cohort[1], S = d$S[1], n_per_group = d$n_per_group[1], n_dropped = d$n_dropped[1],
      patient_check = d$patient_check[1], n_dup_removed = d$n_dup_removed[1],
      frac_len_bad = d$frac_len_bad[1], n_rounds = nrow(d), n_delta_na = sum(is.na(d$delta)),
      delta_realised_median = median(d$delta_realised, na.rm = TRUE),
      frac_lfc_mismatch_mean = mean(d$frac_lfc_mismatch, na.rm = TRUE))))
  }
  PT  <- tryCatch(paired_vs_tox(M), error = function(e) { message("paired vs TOX: ", conditionMessage(e)); NULL })
  BA  <- tryCatch(paired_b_vs_a(M), error = function(e) { message("paired B vs A: ", conditionMessage(e)); NULL })
  RSC <- rs_checks(man)
  GT  <- tryCatch(rs_gates(S), error = function(e) { message("RS gates: ", conditionMessage(e)); NULL })

  w <- function(x, f) if (!is.null(x)) write.csv(x, file.path(RCFG$out_dir, f), row.names = FALSE)
  w(M, "power_per_round.csv"); w(S, "power_summary.csv"); w(PV, "power_paired_vs_ref.csv")
  w(mean_se(CU, c(CELL, "method", "fdr"), "tpr"), "power_fdr_tpr_curve.csv")
  if (!is.null(PR)) w(mean_se(PR, c(CELL, "method", "recall"), "precision"), "power_pr_curve.csv")
  if (!is.null(ST)) w(mean_se(ST, c(CELL, "method", "stratum_type", "stratum"), c("n", "tpr_nom05", "tpr_ach05")),
                      "power_stratified.csv")
  if (!is.null(RAD)) w(mean_se(RAD, c(CELL, "method", "expr_decile"), c("n", "scan_radius_median", "frac_radius_gt2")),
                       "power_scan_radius.csv")
  w(OG, "power_oracle_gap.csv"); w(PD, "power_pool_diag.csv")
  w(PT, "power_paired_vs_tox.csv"); w(BA, "power_paired_B_vs_A.csv")
  w(BK, "power_breakdown.csv"); w(RC, "power_real_checks.csv"); w(RSC, "power_rs_checks.csv")
  if (!is.null(GT)) w(rbind_fill(GT), "power_rs_gates.csv")

  withCallingHandlers(report(S, PV, OG, PD, BK, RC, PT, BA, RSC, GT, M), warning = function(w) {
    message("report warning: ", conditionMessage(w)); invokeRestart("muffleWarning") })
  tryCatch(make_plots(M, CU, PR, ST, PT, BA), error = function(e) message("plotting failed: ", conditionMessage(e)))
  cat("\nWrote CSVs and plots to ", RCFG$out_dir, "/\n", sep = "")
}
