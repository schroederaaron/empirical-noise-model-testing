#!/usr/bin/env Rscript
# benchmark_results_tables.R
# -----------------------------------------------------------------------------
# Builds docs/benchmark_results.md from docs/benchmark_results.template.md: every
# {{NAME}} line of the template is replaced by a markdown table computed from the
# summary CSVs of the power and kNN runs (no number is copied by hand).
#
#   Rscript common/benchmark_results_tables.R \
#       [TCGA_test/results] [Yeast_test/results] \
#       [docs/benchmark_results.template.md] [docs/benchmark_results.md]
# Reads <tcga>/power_out/power_summary.csv, <tcga>/knn_out/knn_summary.csv,
# <yeast>/yeast_power_out/power_summary.csv, <yeast>/yeast_knn_out/knn_summary.csv.
# Fails if a placeholder is missing or left unfilled.
# -----------------------------------------------------------------------------
def <- c("TCGA_test/results", "Yeast_test/results", "docs/benchmark_results.template.md",
         "docs/benchmark_results.md")
a <- def; a[seq_along(commandArgs(TRUE))] <- commandArgs(TRUE); T_DIR <- a[1]; Y_DIR <- a[2]; TPL <- a[3]; OUT <- a[4]
f3 <- function(x) ifelse(is.na(x), "–", formatC(x, format = "f", digits = 3))
f4 <- function(x) ifelse(is.na(x), "–", formatC(x, format = "f", digits = 4))
esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)
md <- function(df) {
  df[] <- lapply(df, function(z) esc(as.character(z)))
  paste(c(paste0("| ", paste(names(df), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
          apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))), collapse = "\n")
}
agg <- function(d, by, k, f = mean) aggregate(d[k], d[by], f)
ord <- function(df, cols) df[do.call(order, unname(df[cols])), , drop = FALSE]

PT <- read.csv(file.path(T_DIR, "power_out/power_summary.csv"))
PY <- read.csv(file.path(Y_DIR, "yeast_power_out/power_summary.csv"))
KT <- read.csv(file.path(T_DIR, "knn_out/knn_summary.csv"))
KY <- read.csv(file.path(Y_DIR, "yeast_knn_out/knn_summary.csv"))

real_tab <- function(d, part, methods, n_vals = NULL) {
  x <- d[d$part == part & d$method %in% methods, ]
  if (!is.null(n_vals)) x <- x[x$n_rep %in% n_vals, ]
  key <- if (part == "R") "method" else c("n_rep", "method")
  nul <- x[x$pi1 == 0, ]; pw <- x[x$pi1 > 0, ]
  P <- agg(pw, key, c("fdr_nom_05_mean", "tpr_nom_05_mean", "tpr_ach05_mean", "auprc_excess_pe_mean"))
  W <- agg(pw, key, "fdr_nom_05_mean", max); names(W)[ncol(W)] <- "fdr_max"
  P <- merge(P, W)
  has_null <- nrow(nul) > 0
  if (has_null) P <- merge(P, agg(nul, key, c("p_le001_null_mean", "p_le01_null_mean")), all.x = TRUE)
  P$mi <- match(P$method, methods); P <- ord(P, c(setdiff(key, "method"), "mi"))
  out <- data.frame(n = if (part == "R") NA else P$n_rep, method = P$method, check.names = FALSE)
  if (has_null) { out$`null p≤0.001` <- f4(P$p_le001_null_mean); out$`null p≤0.01` <- f4(P$p_le01_null_mean) }
  out$`FDR@0.05` <- f3(P$fdr_nom_05_mean); out$`FDR@0.05 worst cohort` <- f3(P$fdr_max)
  out$`TPR@nominal 0.05` <- f3(P$tpr_nom_05_mean); out$`TPR@achieved 0.05` <- f3(P$tpr_ach05_mean)
  out$`AUPRC excess` <- f3(P$auprc_excess_pe_mean)
  if (part == "R") out$n <- NULL
  md(out)
}
coh_tab <- function(d, part, methods, n_vals = NULL) {
  x <- d[d$part == part & d$pi1 > 0 & d$method %in% methods, ]
  if (!is.null(n_vals)) x <- x[x$n_rep %in% n_vals, ]
  x$mi <- match(x$method, methods); x <- ord(x, c("n_rep", "cohort", "mi"))
  out <- data.frame(n = x$n_rep, cohort = x$cohort, method = x$method,
                    `FDR@0.05` = f3(x$fdr_nom_05_mean), `TPR@nominal` = f3(x$tpr_nom_05_mean),
                    `TPR@achieved` = f3(x$tpr_ach05_mean), `AUPRC excess` = f3(x$auprc_excess_pe_mean),
                    check.names = FALSE)
  if (part == "RS") {
    nul <- d[d$part == part & d$pi1 == 0 & d$method %in% methods, ]
    out$`null p≤0.001` <- f4(nul$p_le001_null_mean[match(paste(x$n_rep, x$cohort, x$method),
                                                          paste(nul$n_rep, nul$cohort, nul$method))])
  }
  md(out)
}
knn_tab <- function(k, part, sd, rows) {
  x <- k[k$part == part & abs(k$sigma_d - sd) < 1e-6, ]
  x$id <- paste(x$family, ifelse(is.na(x$scaling), "", x$scaling), x$config)
  x <- x[x$id %in% rows, ]; x$mi <- match(x$id, rows); x <- ord(x, c("n_rep", "mi"))
  md(data.frame(n = x$n_rep, family = x$family, scaling = ifelse(is.na(x$scaling), "", x$scaling),
                config = x$config, `null p≤0.001` = f4(x$p_le001_null), `tail fold` = formatC(x$tail_fold, format = "f", digits = 2),
                `null p≤0.01` = f4(x$p_le01_null), `FDR@0.05` = f3(x$fdr_nom_05), `TPR@nominal` = f3(x$tpr_nom_05),
                `TPR@achieved` = f3(x$tpr_ach05), `AUPRC excess` = f3(x$auprc_excess_pe),
                `p floor` = ifelse(is.na(x$floor_th_median), "–", formatC(x$floor_th_median, format = "e", digits = 1)),
                check.names = FALSE))
}
sim_tab <- function(d, part, by, methods, filt) {
  x <- d[d$part == part & d$method %in% methods, ]; x <- x[filt(x), ]
  nul <- agg(x[x$pi1 == 0, ], c(by, "method"), c("p_le001_null_mean", "p_le01_null_mean"))
  pw  <- agg(x[x$pi1 > 0, ],  c(by, "method"), c("fdr_nom_05_mean", "tpr_ach05_mean", "auprc_excess_pe_mean"))
  m <- merge(nul, pw, all = TRUE); m$mi <- match(m$method, methods); m <- ord(m, c(by, "mi"))
  for (b in by) if (is.numeric(m[[b]])) m[[b]] <- signif(m[[b]], 3)
  md(data.frame(m[by], method = m$method, `null p≤0.001` = f4(m$p_le001_null_mean),
                `null p≤0.01` = f4(m$p_le01_null_mean), `FDR@0.05` = f3(m$fdr_nom_05_mean),
                `TPR@achieved` = f3(m$tpr_ach05_mean), `AUPRC excess` = f3(m$auprc_excess_pe_mean), check.names = FALSE))
}

REFS <- c("limma", "edgeR", "DESeq2")
HEAD <- c(REFS, "TOX-log|TPM-centred", "TOX-log-xfit-js", "TOX-log-xfit-raw")
VAR  <- c("TOX-log|TPM-centred", "TOX-log-xfit-js", "TOX-log-xfit-raw", "TOX-log-xfit-js-loess",
          "TOX-log-xfit-js-nocf", "TOX-log-vband-c3", "TOX-log-vband-js-c0.5", "TOX-log-vband-js-c1",
          "TOX-log-adapt-r0.25", "TOX-log-adapt-r0.5", "TOX-log-adapt-r1", "TOX-log-adapt-xfit-js")
GRID <- c(REFS, "TOX-log|TPM", "TOX-log|TPM-centred", "TOX-log|counts", "TOX-log|counts-centred",
          "TOX-log|TMM", "TOX-log|MoR", "TOX-log-blocked|TPM-centred", "TOX-log-blocked|TMM",
          "TOX-raw|TPM", "TOX-raw|TMM")
KROWS <- c("reference  limma", "reference  edgeR", "reference  DESeq2",
           "mean none k20-50", "mean none fix200",
           "mean xfit_js k20-50", "mean xfit_js fix100", "mean xfit_js fix200", "mean xfit_js fix400",
           "mean xfit_js k50-200-notau", "mean xfit_js k50-200-t0.25",
           "band xfit_js k200-c5", "adaptive xfit_js k200-r0.5", "adaptive none k200-r0.5")
KROWS_P <- c("reference  limma", "reference  edgeR", "reference  DESeq2", "mean none k20-50",
             "mean xfit_js k20-50", "mean xfit_js fix200", "mean xfit_js fix400")
ell <- KY[KY$family == "ellipse", ]
ELL <- md(aggregate(cbind(`null p≤0.01` = ell$p_le01_null, `FDR@0.05` = ell$fdr_nom_05,
                          `TPR@achieved` = ell$tpr_ach05) ~ n_rep + scaling, ell,
                    function(z) sprintf("%.3f–%.3f", min(z), max(z))))
SIMM <- c("limma", "edgeR", "TOX-log", "TOX-log-xfit-js", "TOX-log-xfit-raw", "TOX-log-xfit-js-nocf")

lam <- function(d, ds) {
  x <- d[d$method == "TOX-log-xfit-js" & d$part %in% c("RS", "R") & d$pi1 > 0, ]
  a <- aggregate(cbind(vb_lambda_pool_mean, vb_lambda_case_mean, vb_lambda_control_mean) ~ part + n_rep, x, mean)
  a <- a[order(a$part != "RS", a$n_rep), ]
  data.frame(data = ds, part = a$part, n = a$n_rep, `lambda pooled` = f3(a$vb_lambda_pool_mean),
             `lambda case` = f3(a$vb_lambda_case_mean), `lambda control` = f3(a$vb_lambda_control_mean),
             check.names = FALSE)
}
TABS <- list(
  LAMBDA = md(rbind(lam(PT, "TCGA"), lam(PY, "yeast"))),
  RS_T = real_tab(PT, "RS", HEAD), RS_Y = real_tab(PY, "RS", HEAD),
  RSVAR_T = real_tab(PT, "RS", VAR, c(3, 5)), RSVAR_Y = real_tab(PY, "RS", VAR, c(3, 5)),
  RSCOH_T = coh_tab(PT, "RS", c("edgeR", "TOX-log-xfit-js"), c(3, 5)),
  RSCOH_Y = coh_tab(PY, "RS", c("edgeR", "TOX-log-xfit-js"), c(3, 5)),
  R_T = real_tab(PT, "R", c(REFS, VAR)), R_Y = real_tab(PY, "R", c(REFS, VAR)),
  RGRID_T = real_tab(PT, "R", GRID), RGRID_Y = real_tab(PY, "R", GRID),
  RCOH_T = coh_tab(PT, "R", c("limma", "edgeR", "TOX-log-xfit-js", "TOX-log|TPM-centred")),
  RCOH_Y = coh_tab(PY, "R", c("limma", "edgeR", "DESeq2", "TOX-log-xfit-js", "TOX-log|TPM-centred", "TOX-log|TMM")),
  KNN_T = knn_tab(KT, "RS", 0, KROWS), KNN_Y = knn_tab(KY, "RS", 0, KROWS),
  KNNP0 = knn_tab(KT, "P+", 0, KROWS_P), KNNP1 = knn_tab(KT, "P+", max(KT$sigma_d), KROWS_P),
  ELL_Y = ELL,
  PPLUS = sim_tab(PT, "P+", c("dist", "n_rep", "sigma_d"), SIMM, function(x) x$n_rep %in% c(3, 10) & x$sigma_d %in% range(x$sigma_d)),
  U = sim_tab(PT, "U", c("dist", "n_rep", "rho"), SIMM[1:5], function(x) x$rho %in% c(1, 1.5, 2))
)
tpl <- readLines(TPL)
for (nm in names(TABS)) {
  i <- which(tpl == paste0("{{", nm, "}}"))
  if (length(i) != 1) stop("placeholder ", nm, " found ", length(i), " times")
  tpl[i] <- TABS[[nm]]
}
left <- grep("\\{\\{[A-Z_0-9]+\\}\\}", tpl, value = TRUE)
if (length(left)) stop("unfilled: ", paste(left, collapse = ", "))
writeLines(tpl, OUT)
