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
#   power_summary.csv         mean / SE over rounds per (cell, method)
#   power_paired_vs_ref.csv   paired per-round difference vs limma
#   power_fdr_tpr_curve.csv   TPR on an achieved-FDR grid
#   power_pr_curve.csv        PR curve (interpolated envelope, display only)
#   power_stratified.csv      TPR by |true log2 FC| bin and by expression decile
#   power_oracle_gap.csv      Part Q oracle - production, paired per round
#   power_pool_diag.csv       Part Q pool diagnostics, oracle - production, paired
#   power_breakdown.csv       Part C breakdown Delta (FDR@0.05 > 0.10), bootstrap CI
#   power_real_checks.csv     Part R per-cohort checks
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
  grid_arms   = c("TOX-log", "TOX-log-blocked", "TOX-raw", "TOX-raw-blocked")
)
dir.create(RCFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# Categorical slots, fixed order (never cycled); a method keeps its colour across plots.
PALETTE <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
             "#e87ba4", "#008300", "#4a3aa7", "#e34948")
PART_LABEL <- c(S = "STRESS ARM (bimodal) -- not for ranking methods")

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
  out <- list(metrics = list(), curve = list(), pr = list(), strata = list())
  for (m in unique(rs$method)) {
    r <- rs[rs$method == m, ]; r <- r[match(tk$gene, r$gene), ]
    add <- function(x) if (!is.null(x)) cbind(meta, method = m, x, row.names = NULL)
    out$metrics[[m]] <- add(metrics_one(r, tk$is_de, tk$true_lfc, n_de_total, tk$is_hv,
                                        lfc_log2 = identical(r$lfc_scale[1], "log2")))
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

make_plots <- function(M, CU, PR, ST) {
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
  refs <- c("limma", "edgeR", "DESeq2")
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
}

# =============================================================================
# 3. CONSOLE REPORT
# =============================================================================

print_tab <- function(x, title, notes = NULL) {
  message(strrep("=", 78)); message(title); message(strrep("=", 78))
  for (n in notes) message("  ", n)
  num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], signif, 4)
  x <- x[, vapply(x, function(v) !all(is.na(v)), logical(1)), drop = FALSE]   # drop all-NA columns
  print(x, row.names = FALSE)
}

spread <- function(x) { x <- x[is.finite(x)]; if (length(x)) diff(range(x)) else NA_real_ }

report <- function(S, PV, OG, PD, BK, RC) {
  cols <- c("tpr_ach05", "tpr_nom_05", "fdr_nom_05", "zero_disc_05", "fpr_null05", "p_le05_null",
            "fpr_hv05", "auc_p", "pauc05_pe", "auprc_pe", "auprc_excess_pe", "prec_top_nde",
            "lfc80", "floor_blocks05", "eff_lfc_median_abs", "frac_below_floor")
  keep <- c(CELL, "method", "n_rounds", paste0(cols, "_mean"))
  for (pt in unique(S$part)) {
    x <- S[S$part == pt, intersect(keep, names(S))]
    x <- x[do.call(order, unname(as.list(x[intersect(c(CELL, "method"), names(x))]))), ]
    print_tab(x, sprintf("PART %s -- mean over rounds (SE in power_summary.csv)%s", pt,
                         if (pt %in% names(PART_LABEL)) paste0("\n", PART_LABEL[[pt]]) else ""))
  }
  message("\n  tpr_ach05   TPR at ACHIEVED FDR 0.05 (truth-thresholded, max over steps; not user-realisable).")
  message("  tpr_nom_05 / fdr_nom_05  what a user gets at nominal 0.05; fdr zero-filled,")
  message("              zero_disc_05_mean = frac_zero_disc (rounds with no call).")
  message("  auprc_excess_pe  AUPRC minus prevalence: comparable across pi1.")
  message("  floor_blocks05  fraction of rounds where the tied floor block was too small for BH.")

  sp <- do.call(rbind, lapply(split(S, do.call(paste, c(S[CELL], sep = "\r"))), function(d)
    data.frame(d[1, CELL], auc_range = spread(d$auc_p_mean), pauc05_range = spread(d$pauc05_p_mean),
               auprc_range = spread(d$auprc_pe_mean))))
  print_tab(sp, "SPREAD ACROSS METHODS per cell (P6: pauc05_range / auprc_range > auc_range?)")
  if (!is.null(PV)) print_tab(PV, paste0("PAIRED DIFFERENCE vs ", RCFG$paired_ref, " (mean +/- SE over rounds)"))
  if (!is.null(OG)) print_tab(OG, "ORACLE GAP (Part Q): oracle minus production, paired per round",
    c("P3': contamination is a CALIBRATION effect -> look at d_fdr_nom_05, d_fpr_null05, d_n_called_05,",
      "d_tpr_nom_05. d_tpr_ach05 is the expected-null control (invariant to monotone p inflation)."))
  if (!is.null(PD)) print_tab(PD, "POOL DIAGNOSTICS (Part Q, exact TOX-log): production values and oracle - production",
    c("d_sd_pool_* < 0 = contamination widened the production pool.",
      "d_mean_span_* > 0 = the oracle's neighbours lie further away in mean (confound, not contamination)."))
  if (!is.null(BK)) print_tab(BK, sprintf("PART C BREAKDOWN: smallest target Delta with FDR@0.05 > %.2f (bootstrap CI)",
                                          RCFG$fdr_limit), "Inf = no breakdown inside the Delta grid.")
  if (!is.null(RC)) print_tab(RC, "PART R CHECKS per base cohort",
    c("frac_len_bad: genes dropped because count/TPM was not constant (length recovery invalid).",
      "frac_lfc_mismatch: injected genes with |realised - intended log2 FC| > lfc_tol.",
      "n_delta_na: rounds where the thinning balance failed."))
}

# =============================================================================
if (sys.nframe() == 0L) {
  man <- read_manifests(RCFG$raw_dir)
  message(sprintf("%d jobs in %s (parts %s); cores %d", nrow(man), RCFG$raw_dir,
                  paste(unique(man$part), collapse = ","), RCFG$n_cores))
  res <- mclapply(split(man, seq_len(nrow(man))), function(j)
           tryCatch(job_metrics(j), error = function(e) {
             message("job ", j$job_id, ": ", conditionMessage(e)); NULL }),
         mc.cores = RCFG$n_cores, mc.preschedule = TRUE)
  res <- Filter(function(x) is.list(x) && !is.null(x$metrics), res)
  if (!length(res)) stop("No job could be scored.")
  M  <- do.call(rbind, lapply(res, `[[`, "metrics"))
  CU <- do.call(rbind, lapply(res, `[[`, "curve"))
  PR <- do.call(rbind, lapply(res, `[[`, "pr"))
  ST <- do.call(rbind, lapply(res, `[[`, "strata"))

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
  RC <- if (any(man$part == "R")) {
    r <- man[man$part == "R", ]
    do.call(rbind, lapply(split(r, r$cohort), function(d) data.frame(
      cohort = d$cohort[1], S = d$S[1], n_per_group = d$n_per_group[1], n_dropped = d$n_dropped[1],
      patient_check = d$patient_check[1], n_dup_removed = d$n_dup_removed[1],
      frac_len_bad = d$frac_len_bad[1], n_rounds = nrow(d), n_delta_na = sum(is.na(d$delta)),
      delta_realised_median = median(d$delta_realised, na.rm = TRUE),
      frac_lfc_mismatch_mean = mean(d$frac_lfc_mismatch, na.rm = TRUE))))
  }

  w <- function(x, f) if (!is.null(x)) write.csv(x, file.path(RCFG$out_dir, f), row.names = FALSE)
  w(M, "power_per_round.csv"); w(S, "power_summary.csv"); w(PV, "power_paired_vs_ref.csv")
  w(mean_se(CU, c(CELL, "method", "fdr"), "tpr"), "power_fdr_tpr_curve.csv")
  if (!is.null(PR)) w(mean_se(PR, c(CELL, "method", "recall"), "precision"), "power_pr_curve.csv")
  if (!is.null(ST)) w(mean_se(ST, c(CELL, "method", "stratum_type", "stratum"), c("n", "tpr_nom05", "tpr_ach05")),
                      "power_stratified.csv")
  w(OG, "power_oracle_gap.csv"); w(PD, "power_pool_diag.csv")
  w(BK, "power_breakdown.csv"); w(RC, "power_real_checks.csv")

  withCallingHandlers(report(S, PV, OG, PD, BK, RC), warning = function(w) {
    message("report warning: ", conditionMessage(w)); invokeRestart("muffleWarning") })
  tryCatch(make_plots(M, CU, PR, ST), error = function(e) message("plotting failed: ", conditionMessage(e)))
  message("\nWrote CSVs and plots to ", RCFG$out_dir, "/")
}
