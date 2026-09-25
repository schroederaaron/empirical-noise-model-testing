# report.R -- the summary table (brief 8.5), assembled ONLY from files that
# 02_run_gof.R wrote. Nothing here fits a model.

.rd <- function(dir, f) { p <- file.path(dir, f); if (file.exists(p)) data.table::fread(p) else NULL }

#' One row per model:
#'   W2, A2, D (on G_BOOT) . W2/N . A2/N . p_GOF (W2, A2, D) . excess(W2) .
#'   held-out score/obs . delta vs NB . rank_heldout . rank_excess_W2 .
#'   failed fits . GP mass flags
#' Two rank columns, not one composite: which criterion "Rank" uses is decision D5.
build_summary_table <- function(dir) {
  obs  <- .rd(dir, "stats_gboot_obs.csv")
  if (is.null(obs)) stop("stats_gboot_obs.csv not found in ", dir, " -- was 02_run_gof.R run?")
  tab <- obs[, .(model, N, W2, A2, D, W2_N, A2_N)]
  bs <- .rd(dir, "boot_summary.csv")
  if (!is.null(bs)) {
    w <- data.table::dcast(bs, model ~ stat, value.var = c("p_GOF", "excess", "z"))
    tab <- merge(tab, w[, .(model, p_GOF_W2, p_GOF_A2, p_GOF_D, excess_W2, z_W2)], by = "model", all.x = TRUE)
    bo <- unique(bs[, .(model, B_ok, B_failed, n_rep_drop_gt_1pct)])
    tab <- merge(tab, bo, by = "model", all.x = TRUE)
  } else tab[, `:=`(p_GOF_W2 = NA_real_, p_GOF_A2 = NA_real_, p_GOF_D = NA_real_,
                    excess_W2 = NA_real_, z_W2 = NA_real_, B_ok = NA_integer_)]
  ho <- .rd(dir, "heldout_summary.csv")
  if (!is.null(ho)) tab <- merge(tab, ho[, .(model, heldout_per_obs = per_obs,
                                             delta_vs_nb = delta_vs_nb_mean,
                                             frac_genes_better_than_nb)], by = "model", all.x = TRUE)
  else tab[, `:=`(heldout_per_obs = NA_real_, delta_vs_nb = NA_real_)]
  fd <- .rd(dir, "fit_failures.csv")
  if (!is.null(fd)) tab <- merge(tab, fd[, .(model, n_failed_obs = n_failed, gp_mass_flags)], by = "model", all.x = TRUE)
  tab[, rank_heldout := if (all(is.na(heldout_per_obs))) NA_integer_ else
        as.integer(rank(-heldout_per_obs, ties.method = "min", na.last = "keep"))]
  tab[, rank_excess_W2 := if (all(is.na(excess_W2))) NA_integer_ else
        as.integer(rank(excess_W2, ties.method = "min", na.last = "keep"))]
  tab[order(match(model, ALL_MODELS))]
}

TABLE_CAPTION <- paste(
  "W2/A2/D and p_GOF are computed on G_BOOT with the parametric bootstrap as the reference.",
  "KNOWN LIMIT: genes are simulated independently given their fitted parameters; gene-gene",
  "correlation and shared sample-level effects are not in the null, so p_GOF is anti-conservative",
  "if they are present. At RNA-seq N, p_GOF = 1/(B+1) for every model is the expected outcome,",
  "not a finding. excess(W2) = T_obs / median(T_b) is a DESCRIPTIVE cross-model comparison that",
  "corrects for each model's own in-sample estimation effect; it is not a model-selection criterion.",
  "Two ranks are shown (decision D5).")

write_summary_md <- function(tab, path, cfg_line = "") {
  fmt <- function(x) if (is.numeric(x)) formatC(x, digits = 4, format = "g") else as.character(x)
  hdr <- paste("|", paste(names(tab), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(tab)), collapse = " | "), "|")
  body <- apply(tab, 1, function(r) paste("|", paste(vapply(r, fmt, ""), collapse = " | "), "|"))
  writeLines(c("# GOF summary", "", cfg_line, "", hdr, sep, body, "", paste0("*", TABLE_CAPTION, "*")), path)
}
