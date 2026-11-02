#!/usr/bin/env Rscript
# variance_ratio_real_contrasts.R
# -----------------------------------------------------------------------------
# Implementation note G, analysis A4: variance imbalance in REAL contrasts (healthy
# vs Stage IV per cancer; no truth). For each contrast:
#   * the distribution of L = log2(s2_A / T_A(m_A)) - log2(s2_B / T_B(m_B)) (A =
#     healthy, B = Stage IV; common/variance_diagnostics.R), and the share of genes
#     beyond a null threshold;
#   * the systematic trend ratio T_A(m) / T_B(m) by expression decile (the
#     trend-level group-variance ratio of note E v2 section 2.5).
# Null threshold: 95th percentile of |L| over N_SPLITS random halves of the healthy
# cohort (group sizes = half the healthy cohort, so it is a within-condition
# reference, not matched to the contrast's n -- state this when reading it).
#
#   Rscript TCGA_test/scripts/variance_ratio_real_contrasts.R <out_dir> \
#       coad_healthy.rds coad_stageiv.rds luad_healthy.rds luad_stageiv.rds ...
# Arguments after <out_dir> are (healthy, Stage IV) pairs of cohort RDS files from
# export_power_cohort.R. Writes variance_ratio_real_contrasts.csv and
# variance_ratio_real_contrasts.png.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) %% 2L != 1L)
  stop("usage: Rscript variance_ratio_real_contrasts.R <out_dir> <healthy.rds> <stageiv.rds> [...]")
out_dir <- args[1]; files <- matrix(args[-1], ncol = 2, byrow = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COMMON_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand <- c(if (length(f)) file.path(dirname(normalizePath(f[1])), "..", "..", "common"),
            "analysis/common", "experiments/Noise_Model_Test/common", "common")
  hit <- cand[file.exists(file.path(cand, "variance_diagnostics.R"))]
  if (!length(hit)) stop("common/variance_diagnostics.R not found")
  normalizePath(hit[1])
})
source(file.path(COMMON_DIR, "variance_diagnostics.R"))
.libPaths(c(normalizePath("external/docker_r_libs", mustWork = FALSE), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

N_SPLITS <- 20L
MIN_MEAN_COUNT <- 10      # as power_test.R Part R

#' TPM (genes x samples) of one cohort RDS: its own $tpm, or counts / lengths.
cohort_tpm <- function(path) {
  d <- readRDS(path)
  cnt <- as.matrix(d$counts)
  tp <- if (!is.null(d$tpm)) {
    t0 <- as.matrix(d$tpm)
    if (setequal(rownames(t0), colnames(cnt))) t0 <- t(t0)
    t0[rownames(cnt), colnames(cnt)]
  } else {
    r <- cnt / as.numeric(d$lengths[rownames(cnt)]); sweep(r, 2, colSums(r), "/") * 1e6
  }
  list(tpm = tp, counts = cnt,
       label = if (!is.null(d$project_id)) paste0(d$project_id, "/", d$stage) else basename(path))
}

set.seed(20261028L)
rows <- list()
for (i in seq_len(nrow(files))) {
  h <- cohort_tpm(files[i, 1]); s4 <- cohort_tpm(files[i, 2])
  genes <- intersect(rownames(h$tpm), rownames(s4$tpm))
  genes <- genes[rowMeans(h$counts[genes, ]) >= MIN_MEAN_COUNT & rowMeans(s4$counts[genes, ]) >= MIN_MEAN_COUNT]
  x <- cbind(h$tpm[genes, ], s4$tpm[genes, ])
  grp <- rep(c("A", "B"), c(ncol(h$tpm), ncol(s4$tpm)))
  vq <- vdiag_quantities(vdiag_frame(x, grp))
  # null threshold from random halves of the healthy cohort
  nh <- ncol(h$tpm) %/% 2L
  thr <- unname(quantile(unlist(lapply(seq_len(N_SPLITS), function(b) {
    cols <- sample.int(ncol(h$tpm), 2L * nh)
    q <- vdiag_quantities(vdiag_frame(h$tpm[genes, cols], rep(c("A", "B"), each = nh)))$absL
    q[is.finite(q)]
  })), 1 - VD$top))
  contrast <- sub("/.*$", "", h$label)
  L <- vq$L[is.finite(vq$L)]
  rows[[length(rows) + 1L]] <- data.frame(contrast = contrast, n_A = sum(grp == "A"), n_B = sum(grp == "B"),
    n_genes = length(L), stat = c(paste0("L_q", c(5, 25, 50, 75, 95)), "share_absL_beyond_null", "null_threshold"),
    decile = NA_integer_, value = c(quantile(L, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE),
                                    mean(abs(L) > thr), thr))
  dec <- cut(vq$m, unique(quantile(vq$m, seq(0, 1, 0.1))), include.lowest = TRUE, labels = FALSE)
  # trend ratio at the common axis: both group trends evaluated at the pooled mean
  ta <- vd_trend(vq$s2_A, vq$m); tb <- vd_trend(vq$s2_B, vq$m)
  agg <- aggregate(data.frame(value = log2(ta / tb), share = abs(vq$L) > thr), list(decile = dec),
                   function(z) mean(z, na.rm = TRUE))
  rows[[length(rows) + 1L]] <- data.frame(contrast = contrast, n_A = sum(grp == "A"), n_B = sum(grp == "B"),
    n_genes = length(L), stat = rep(c("log2_trend_ratio_A_over_B", "share_absL_beyond_null"), each = nrow(agg)),
    decile = rep(agg$decile, 2), value = c(agg$value, agg$share))
  message(sprintf("%s: %d genes, healthy n = %d, Stage IV n = %d, |L| null threshold %.2f, share beyond %.3f",
                  contrast, length(L), sum(grp == "A"), sum(grp == "B"), thr, mean(abs(L) > thr)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "variance_ratio_real_contrasts.csv"), row.names = FALSE)

pd <- res[!is.na(res$decile), ]
p <- ggplot(pd, aes(decile, value, colour = contrast)) +
  geom_hline(data = data.frame(stat = "share_absL_beyond_null", y = VD$top), aes(yintercept = y), colour = "grey55") +
  geom_hline(data = data.frame(stat = "log2_trend_ratio_A_over_B", y = 0), aes(yintercept = y), colour = "grey55") +
  geom_line() + geom_point(size = 1.2) + facet_wrap(~ stat, scales = "free_y") +
  labs(title = "Note G A4: healthy (A) vs Stage IV (B) per cancer",
       subtitle = "Top: log2 of the group variance trends' ratio; bottom: share of genes beyond the healthy-halves null threshold (line = 5 %)",
       x = "expression decile (1 = lowest)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "variance_ratio_real_contrasts.png"), p, width = 10, height = 5, dpi = 150)
message("wrote ", file.path(out_dir, "variance_ratio_real_contrasts.csv"))
