#!/usr/bin/env Rscript
# 03_report.R -- summary table + all plots from SAVED outputs of 02_run_gof.R.
# No model is fitted here.
#
#   Rscript Distribution_GOF/scripts/03_report.R --RUN_DIR=<an 02_run_gof.R output dir>

GOF_ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(Sys.getenv("GOF_ROOT", "Distribution_GOF"))
})
source(file.path(GOF_ROOT, "R", "setup.R"))
load_or_install("ggplot2")
source(file.path(GOF_ROOT, "R", "plots.R"))
cfg <- GOF_CONFIG
D <- cfg$RUN_DIR
if (!nzchar(D) || !dir.exists(D)) stop("--RUN_DIR=<an 02_run_gof.R output dir> is required")
ensure_fonts()

tab <- build_summary_table(D)
data.table::fwrite(tab, file.path(D, "summary_table.csv"))
run_cfg <- readRDS(file.path(D, "config_resolved.rds"))
write_summary_md(tab, file.path(D, "summary_table.md"),
                 sprintf("Run: `%s` -- %s, B = %d", basename(D), gof_mode_tag(run_cfg), run_cfg$B))
print(tab)

PL <- file.path(D, "plots"); dir.create(PL, showWarnings = FALSE)
gb <- readRDS(file.path(D, "gboot_data.rds"))
mean_str <- rank_bins(log2_mean_cpm(gb$Y, gb$O), 3L, "mean_T"); names(mean_str) <- rownames(gb$Y)
save_plot <- function(p, name, w = 6, h = 4) ggplot2::ggsave(file.path(PL, name), p, width = w, height = h, dpi = 120)
for (m in tab$model) {
  pf <- file.path(D, sprintf("pit_gboot_%s.rds", m)); if (!file.exists(pf)) next
  U <- readRDS(pf)$u
  bf <- file.path(D, sprintf("boot_%s.rds", m))
  bt <- if (file.exists(bf)) readRDS(bf) else NULL
  u <- as.numeric(U)
  save_plot(plot_pit_hist(u, bt$hist, sprintf("%s: PIT histogram (G_BOOT)", m)), sprintf("%s_pit_hist.png", m))
  save_plot(plot_pit_qq(u, sprintf("%s: uniform QQ", m)), sprintf("%s_pit_qq.png", m), 5, 5)
  save_plot(plot_pit_ecdf_diff(u, bt$hist, sprintf("%s: ECDF(u) - u", m)), sprintf("%s_pit_ecdf_diff.png", m))
  save_plot(plot_pit_by_mean(U, mean_str[rownames(U)], sprintf("%s: by mean tertile", m)), sprintf("%s_pit_by_mean.png", m))
  if (!is.null(bt) && !is.null(bt$sample_mean) && nrow(bt$sample_mean))
    save_plot(plot_sample_means(sample_qnorm_mean(U), bt$sample_mean, sprintf("%s: per-sample mean qnorm(u)", m)),
              sprintf("%s_sample_means.png", m), 8, 4)
}
cat("Wrote", file.path(D, "summary_table.csv"), "and plots in", PL, "\n")
