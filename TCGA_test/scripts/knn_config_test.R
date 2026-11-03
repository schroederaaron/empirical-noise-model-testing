#!/usr/bin/env Rscript
# knn_config_test.R
# =============================================================================
# WHICH NEIGHBOURHOOD, IN WHICH CONFIGURATION, DOES WHAT? (owner request 03/11/2026)
#
# A factorial grid of the four neighbourhood constructions of noise_model_variance,
# each without and with cross-fitted James-Stein scaling (xfit_js), at several sizes:
#
#   mean      the exact model's adaptive kNN in mean (k_start / k_step / k_max, tau stop);
#             larger k_max, looser tau, no tau, fixed k
#   band      fixed k genes inside the variance band (raw match), k x band constant c
#   adaptive  block-wise growth with the variance band, k x growth radius r_max
#   ellipse   the k nearest in the joint (log2 mean, variance) distance, k x semi-axis r_max
#
# on Part RS (real TCGA noise, 2n samples per round from the six base cohorts) and on
# Part P+ (simulated, dispersion scatter sigma_d 0 and 1.05), n = 3, 5, 10, null arm
# (pi1 = 0) and power arm (pi1 = 0.1). Plus limma / edgeR / DESeq2 and the comparators.
# The mean / plain / k20-50 arm is the comparator itself (a built-in parity check).
#
# Reuses power_test.R (simulation, sampling, methods, raw output; its main block is the
# function power_main()) and power_report.R (all metrics), each on this script's own
# folders, then writes the configuration summary:
#   knn_arms.csv        every arm and its settings
#   knn_summary.csv     per arm, part, n (and sigma_d): null tail (p <= 0.001 / 0.01 /
#                       0.05), FDR and TPR at nominal 0.05, TPR at achieved FDR 0.05,
#                       AUPRC excess, pool size, p-value floor, share of short pools; the
#                       cohort (RS) or distribution (P+) mean, and the E v2 gate counts
#   knn_tradeoff_<part>.png   tail calibration (x, fold over nominal) vs TPR at nominal
#                       0.05 (y), one point per arm; per n
#   knn_by_size_<part>.png    FDR and null tail against the neighbourhood size
#
#   Rscript TCGA_test/scripts/knn_config_test.R
# Environment (all optional):
#   KNN_PARTS=RS,P+   KNN_ROUNDS=20   KNN_N=3,5,10   KNN_FAMILIES=mean,band,adaptive,ellipse
#   KNN_SCALINGS=none,xfit_js   KNN_REFS=1 (0 = skip limma/edgeR/DESeq2)
#   KNN_RAW_DIR, KNN_OUT_DIR (default: knn_raw / knn_out next to power_raw)
#   POWER_REAL_RDS (cohort RDS files, as for power_test.R), POWER_CORES, POWER_GENES
# =============================================================================

SCRIPT_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "TCGA_test/scripts"
})
source(file.path(SCRIPT_DIR, "power_test.R"))          # definitions only (main is power_main())
load_or_install("ggplot2")

env_list <- function(name, default) trimws(strsplit(Sys.getenv(name, default), ",")[[1]])
KNN <- list(
  parts    = env_list("KNN_PARTS", "RS,P+"),
  rounds   = as.integer(Sys.getenv("KNN_ROUNDS", "20")),
  n_reps   = as.integer(env_list("KNN_N", "3,5,10")),
  families = env_list("KNN_FAMILIES", "mean,band,adaptive,ellipse"),
  scalings = env_list("KNN_SCALINGS", "none,xfit_js"),
  refs     = Sys.getenv("KNN_REFS", "1") == "1",
  raw_dir  = Sys.getenv("KNN_RAW_DIR", file.path(dirname(PCFG$raw_dir), "knn_raw")),
  out_dir  = Sys.getenv("KNN_OUT_DIR", file.path(dirname(PCFG$raw_dir), "knn_out"))
)

#' The arm grid. Labels "knn-<family>-<plain|xfit>-<config>".
knn_grid <- function(families, scalings) {
  arms <- list()
  add <- function(fam, sc, tag, size, ...) {
    arms[[length(arms) + 1L]] <<- c(
      list(label = sprintf("knn-%s-%s-%s", fam, if (sc == "none") "plain" else "xfit", tag),
           engine = "variance", norm = 1L, null = 0L, input = "tpm", vmode = fam, scaling = sc,
           match = "raw", target = "median", family = fam, config = tag, size = size), list(...))
  }
  kc <- function(s, st, m) list(k_start = as.integer(s), k_step = as.integer(st), k_max = as.integer(m))
  for (sc in scalings) {
    if ("mean" %in% families) {
      add("mean", sc, "k20-50",        50,  kcfg = kc(20, 1, 50),    tau = 0.1)    # the comparator's
      add("mean", sc, "k20-100",       100, kcfg = kc(20, 1, 100),   tau = 0.1)
      add("mean", sc, "k20-200",       200, kcfg = kc(20, 1, 200),   tau = 0.1)
      add("mean", sc, "k50-200-t0.25", 200, kcfg = kc(50, 1, 200),   tau = 0.25)
      add("mean", sc, "k50-200-notau", 200, kcfg = kc(50, 1, 200),   tau = 1e6)    # no stop
      add("mean", sc, "fix100",        100, kcfg = kc(100, 100, 100), tau = 0.1)
      add("mean", sc, "fix200",        200, kcfg = kc(200, 200, 200), tau = 0.1)
      add("mean", sc, "fix400",        400, kcfg = kc(400, 400, 400), tau = 0.1)
    }
    if ("band" %in% families) for (k in c(50L, 100L, 200L)) for (cc in c(3, 5))
      add("band", sc, sprintf("k%d-c%g", k, cc), k, k = k, c_tau = cc)
    if ("adaptive" %in% families) for (k in c(50L, 100L, 200L)) for (r in c(0.5, 1))
      add("adaptive", sc, sprintf("k%d-r%g", k, r), k, k = k, r_max = r, c_tau = 3)
    if ("ellipse" %in% families) for (k in c(50L, 100L, 200L)) for (r in c(0.5, 1))
      add("ellipse", sc, sprintf("k%d-r%g", k, r), k, k = k, r_max = r, c_tau = 3)
  }
  arms
}

#' Configuration summary from the report's power_summary.csv.
knn_summary <- function(out_dir, arms) {
  S <- read.csv(file.path(out_dir, "power_summary.csv"), stringsAsFactors = FALSE)
  A <- do.call(rbind, lapply(arms, function(a) data.frame(
    method = a$label, family = a$family, scaling = a$scaling, config = a$config, size = a$size)))
  comp <- c("TOX-log|TPM-centred", "TOX-log")
  S <- S[S$method %in% c(A$method, comp, REF_KNN), ]
  S$family  <- ifelse(S$method %in% A$method, A$family[match(S$method, A$method)],
                      ifelse(S$method %in% comp, "comparator", "reference"))
  S$scaling <- ifelse(S$method %in% A$method, A$scaling[match(S$method, A$method)], NA)
  S$config  <- ifelse(S$method %in% A$method, A$config[match(S$method, A$method)], S$method)
  S$size    <- ifelse(S$method %in% A$method, A$size[match(S$method, A$method)], NA)
  key <- c("part", "n_rep", "sigma_d", "method", "family", "scaling", "config", "size")
  m <- function(x, col) if (paste0(col, "_mean") %in% names(x)) x[[paste0(col, "_mean")]] else NA_real_
  null <- S[S$pi1 == 0, ]; pow <- S[S$pi1 > 0, ]
  g <- RCFG_GATE
  one <- function(x, cols) {           # mean over cohorts (RS) / distributions (P+), and gate counts
    k <- do.call(paste, c(x[key], sep = "\r"))
    do.call(rbind, lapply(split(x, k), function(d) {
      r <- d[1, key]
      for (c in cols) r[[c]] <- mean(m(d, c), na.rm = TRUE)
      r$n_units <- nrow(d)
      r
    }))
  }
  N0 <- one(null, c("p_le001_null", "p_le01_null", "p_le05_null"))
  N0$n_tail_pass <- vapply(split(null, do.call(paste, c(null[key], sep = "\r")))[do.call(paste, c(N0[key], sep = "\r"))],
                           function(d) sum(m(d, "p_le001_null") <= g$tail_001 & m(d, "p_le01_null") <= g$tail_01, na.rm = TRUE),
                           numeric(1))
  P1 <- one(pow, c("fdr_nom_05", "tpr_nom_05", "tpr_ach05", "auprc_excess_pe", "nbhd_own_case_median",
                   "floor_th_median", "frac_short_pool"))
  P1$n_fdr_pass <- vapply(split(pow, do.call(paste, c(pow[key], sep = "\r")))[do.call(paste, c(P1[key], sep = "\r"))],
                          function(d) sum(m(d, "fdr_nom_05") <= g$fdr_max, na.rm = TRUE), numeric(1))
  out <- merge(N0, P1, by = key, all = TRUE, suffixes = c("_null", "_power"))
  out$tail_fold <- out$p_le001_null / 0.001
  out[order(out$part, out$n_rep, out$sigma_d, out$family, out$scaling, out$size), ]
}

REF_KNN <- c("limma", "edgeR", "DESeq2")
RCFG_GATE <- list(tail_001 = 0.002, tail_01 = 0.015, fdr_max = 0.075)    # note E v2 section 6

knn_plots <- function(K, out_dir) {
  fam_col <- c(mean = "#2a78d6", band = "#eb6834", adaptive = "#1baf7a", ellipse = "#eda100",
               comparator = "#4a3aa7", reference = "#555555")
  th <- theme_bw(base_size = 9) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  for (pt in unique(K$part)) {
    d <- K[K$part == pt & is.finite(K$tail_fold) & is.finite(K$tpr_nom_05), ]
    if (!nrow(d)) next
    d$scaling[is.na(d$scaling)] <- "other"
    fac <- if (pt == "P+") facet_grid(sigma_d ~ n_rep, labeller = label_both) else facet_wrap(~ n_rep, labeller = label_both)
    p <- ggplot(d, aes(tail_fold, tpr_nom_05, colour = family, shape = scaling)) +
      geom_vline(xintercept = RCFG_GATE$tail_001 / 0.001, colour = "grey55", linetype = "dashed") +
      geom_point(size = 2) +
      geom_text(data = d[d$family %in% c("reference", "comparator"), ], aes(label = config), size = 2.4,
                vjust = -0.9, show.legend = FALSE) +
      scale_x_log10() + scale_colour_manual(values = fam_col) +
      scale_shape_manual(values = c(none = 1, xfit_js = 16, other = 17)) + fac +
      labs(title = sprintf("Part %s: null tail vs power, one point per neighbourhood configuration", pt),
           subtitle = "x: null P(p <= 0.001) / 0.001 (dashed = E v2 tail gate, 2x); y: TPR at nominal 0.05. Open = no scaling, filled = cross-fit JS.",
           x = "null tail, fold over nominal (log)", y = "TPR at nominal 0.05", colour = NULL, shape = NULL) + th
    ggsave(file.path(out_dir, sprintf("knn_tradeoff_%s.png", gsub("[^A-Za-z]", "plus", pt))), p,
           width = 13, height = if (pt == "P+") 8 else 5, dpi = 150)

    s <- d[!is.na(d$size), ]
    if (!nrow(s)) next
    lg <- rbind(data.frame(s[c("n_rep", "sigma_d", "family", "scaling", "config", "size")], metric = "FDR @ nominal 0.05", v = s$fdr_nom_05),
                data.frame(s[c("n_rep", "sigma_d", "family", "scaling", "config", "size")], metric = "null P(p <= 0.001) / 0.001", v = s$tail_fold))
    lg$series <- paste(lg$family, lg$scaling, sub("^k[0-9]+-?", "", lg$config))
    p <- ggplot(lg, aes(size, v, colour = family, linetype = scaling, group = series)) +
      geom_line(linewidth = 0.5) + geom_point(size = 1.4) + scale_x_log10() +
      scale_colour_manual(values = fam_col) +
      facet_grid(metric ~ n_rep + sigma_d, scales = "free_y", labeller = label_both) +
      labs(title = sprintf("Part %s: FDR and null tail against neighbourhood size (genes)", pt),
           subtitle = "Lines connect configurations that differ only in k (mean: k_max or fixed k).",
           x = "neighbourhood size k (log)", y = NULL, colour = NULL, linetype = NULL) + th
    ggsave(file.path(out_dir, sprintf("knn_by_size_%s.png", gsub("[^A-Za-z]", "plus", pt))), p,
           width = 14, height = 7, dpi = 150)
  }
}

# =============================================================================
if (sys.nframe() == 0L) {
  arms <- knn_grid(KNN$families, KNN$scalings)
  PCFG$parts      <- KNN$parts
  PCFG$raw_dir    <- KNN$raw_dir
  PCFG$n_rounds   <- KNN$rounds
  PCFG$RS$n_reps  <- KNN$n_reps
  PCFG$RS$n_rounds <- KNN$rounds
  PCFG$`P+`$n_reps <- KNN$n_reps
  PCFG$`P+`$sigma_mult <- c(0, 1)
  PCFG$vmatch$parts <- c("RS", "P+")
  PCFG$vmatch_arms <- arms
  PCFG$tox_arms <- Filter(function(a) a$label == "TOX-log", PCFG$tox_arms)   # P+ comparator only
  if (!KNN$refs) PCFG$run_edger <- PCFG$run_limma <- PCFG$run_deseq2 <- FALSE
  message(sprintf("knn_config_test: %d arms (%s x %s), parts %s, n %s, %d rounds",
                  length(arms), paste(KNN$families, collapse = "/"), paste(KNN$scalings, collapse = "/"),
                  paste(KNN$parts, collapse = ","), paste(KNN$n_reps, collapse = ","), KNN$rounds))

  power_main()

  dir.create(KNN$out_dir, showWarnings = FALSE, recursive = TRUE)
  Sys.setenv(POWER_RAW_DIR = KNN$raw_dir, POWER_OUT_DIR = KNN$out_dir)
  rc <- system2("Rscript", file.path(SCRIPT_DIR, "power_report.R"))
  if (rc != 0) stop("power_report.R failed (exit ", rc, ")")

  write.csv(do.call(rbind, lapply(arms, function(a) data.frame(
    method = a$label, family = a$family, scaling = a$scaling, config = a$config, size = a$size,
    k = a[["k"]] %||% NA, c_tau = a$c_tau %||% NA, r_max = a$r_max %||% NA,
    k_start = a$kcfg$k_start %||% NA, k_step = a$kcfg$k_step %||% NA, k_max = a$kcfg$k_max %||% NA,
    tau = a$tau %||% NA))), file.path(KNN$out_dir, "knn_arms.csv"), row.names = FALSE)
  K <- knn_summary(KNN$out_dir, arms)
  write.csv(K, file.path(KNN$out_dir, "knn_summary.csv"), row.names = FALSE)
  for (pt in unique(K$part)) for (n in sort(unique(K$n_rep[K$part == pt])))
    print_tab(K[K$part == pt & K$n_rep == n, c("sigma_d", "family", "scaling", "config", "tail_fold", "p_le01_null",
                                              "fdr_nom_05", "tpr_nom_05", "auprc_excess_pe", "n_tail_pass", "n_fdr_pass")],
              sprintf("KNN CONFIGURATIONS | part %s | n = %d (means over cohorts / distributions)", pt, n),
              csv = "knn_summary.csv", n = 80L, keys = c("sigma_d", "family", "scaling", "config"))
  tryCatch(knn_plots(K, KNN$out_dir), error = function(e) message("knn plots failed: ", conditionMessage(e)))
  cat("\nknn_config_test: wrote ", KNN$out_dir, "/knn_{arms,summary}.csv and knn_*.png\n", sep = "")
}
