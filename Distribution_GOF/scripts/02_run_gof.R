#!/usr/bin/env Rscript
# 02_run_gof.R -- main GOF analysis for ONE dataset (issue #188).
#
#   Rscript Distribution_GOF/scripts/02_run_gof.R --DATASET=path/to/ds.rds [--NAME=value ...]
#   Rscript Distribution_GOF/scripts/02_run_gof.R --dry-run          # print the resolved config, exit
#
# Run from the directory that holds external/docker_r_libs (repo convention).
# Every knob is in config_gof.R. Output: <OUT_ROOT>/<label>_<round-..._off-...>/
#
# Steps:
#   1  load + validate the dataset, apply EXCLUDE_SAMPLES (D6), integer audit
#   2  ROUNDING (D1), filterByExpr (once, observed data), offsets (D2)
#   3  G_ALL fits for every model; analysis set = genes where ALL models are OK
#   4  PIT (G6 asserts) + observed statistics on G_ALL (descriptive) and on each
#      model's own OK set (sensitivity); re-randomisation spread
#   5  G_BOOT: stratified subset; reference fits with G_BOOT offsets (the procedure
#      the replicates repeat); T_obs on G_BOOT
#   6  parametric bootstrap per model -> p_GOF, excess, z; strata
#   7  held-out predictive score on G_BOOT
# Nothing here writes a number it did not compute; 03_report.R builds the table.

GOF_ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(Sys.getenv("GOF_ROOT", "Distribution_GOF"))
})
source(file.path(GOF_ROOT, "R", "setup.R"))
cfg <- GOF_CONFIG
print_gof_config(cfg)
if (isTRUE(cfg$DRY_RUN)) { cat("Dry run: config resolved, nothing computed.\n"); quit(save = "no", status = 0) }

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")
t_start <- proc.time()[["elapsed"]]

# ---------------------------------------------------------------- 1. input
ds <- apply_exclusions(read_dataset(cfg$DATASET), cfg$EXCLUDE_SAMPLES)
OUT <- file.path(cfg$OUT_ROOT, paste0(ds$label, "_", gof_mode_tag(cfg)))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say("dataset %s (%s): %d genes x %d samples -> %s", ds$label, ds$source, nrow(ds$counts), ncol(ds$counts), OUT)
if (cfg$OFFSET_MODE == "tmm_length" && is.null(ds$lengths))
  stop("OFFSET_MODE = 'tmm_length' needs a dataset with $lengths (tximport input)")

audit <- integer_audit(ds$counts)
data.table::fwrite(audit, file.path(OUT, "input_audit.csv"))
say("integer audit: %.4g%% non-integer entries", 100 * audit[section == "overall" & metric == "frac_non_integer", value])

# ---------------------------------------------------------------- 2. preprocess
rd <- apply_rounding(ds$counts, cfg$ROUNDING, seed_for(cfg$BASE_SEED, "rounding"))
Yraw <- rd$Y
X <- stats::model.matrix(ds$design, data = ds$samples)
if (qr(X)$rank < ncol(X)) stop("design matrix is rank deficient")
grp <- design_group(ds$samples, ds$design)
flt <- filter_genes(Yraw, grp)
Y <- Yraw[flt$keep, , drop = FALSE]
L <- if (is.null(ds$lengths)) NULL else ds$lengths[flt$keep, , drop = FALSE]
O <- compute_offsets(Y, cfg$OFFSET_MODE, L)
data.table::fwrite(data.table::data.table(n_genes_before = flt$n_before, n_genes_after = flt$n_after,
                                          n_entries_changed_by_rounding = rd$n_changed,
                                          rounding = cfg$ROUNDING, offset_mode = cfg$OFFSET_MODE),
                   file.path(OUT, "preprocess.csv"))
say("filterByExpr: %d -> %d genes; rounding changed %d entries", flt$n_before, flt$n_after, rd$n_changed)

G_ALL <- rownames(Y)
if (cfg$MAX_GENES_OBS > 0L && cfg$MAX_GENES_OBS < length(G_ALL))       # dev cap, stratified like G_BOOT
  G_ALL <- select_boot_genes(G_ALL, Y, O, cfg$MAX_GENES_OBS, seed_for(cfg$BASE_SEED, "devcap"))
dec <- rank_bins(log2_mean_cpm(Y[G_ALL, , drop = FALSE], O[G_ALL, , drop = FALSE]), 10L, "D")
names(dec) <- G_ALL

fams <- gof_families(cfg)
models <- names(fams)

# ---------------------------------------------------------------- 3. G_ALL fits
fits <- list()
for (m in models) {
  t0 <- proc.time()[["elapsed"]]
  fits[[m]] <- fit_dataset(Y, X, O, fams[[m]], G_ALL, cfg$N_CORES)
  saveRDS(fits[[m]], file.path(OUT, sprintf("fits_gall_%s.rds", m)))
  data.table::fwrite(fits[[m]]$diag, file.path(OUT, sprintf("fit_diag_gall_%s.csv", m)))
  data.table::fwrite(params_table(fits[[m]]), file.path(OUT, sprintf("params_gall_%s.csv", m)))
  say("fit %-8s on G_ALL: %d/%d OK  (%.1f s, median %.3f s/gene)", m, sum(fits[[m]]$ok), length(G_ALL),
      proc.time()[["elapsed"]] - t0, stats::median(fits[[m]]$diag$time_s, na.rm = TRUE))
}
all_ok <- Reduce(`&`, lapply(fits, function(f) f$ok[G_ALL]))
A_SET <- G_ALL[all_ok]
writeLines(A_SET, file.path(OUT, "analysis_genes.txt"))

fail <- data.table::rbindlist(lapply(models, function(m) {
  bad <- G_ALL[!fits[[m]]$ok[G_ALL]]
  data.table::data.table(model = m, n_failed = length(bad),
                         gp_mass_flags = if (m == "genpois") sum(fits[[m]]$diag$gp_mass_min < 1 - 1e-8, na.rm = TRUE) else NA_integer_,
                         failed_by_mean_decile = paste(names(table(dec[bad])), table(dec[bad]), sep = ":", collapse = " "),
                         top_messages = paste(utils::head(names(sort(table(fits[[m]]$diag[!(ok)]$conv_msg), decreasing = TRUE)), 3), collapse = " || "))
}))
data.table::fwrite(fail, file.path(OUT, "fit_failures.csv"))
say("analysis set (all models OK): %d of %d genes", length(A_SET), length(G_ALL))
if (length(A_SET) < 10L) stop("fewer than 10 genes where every model fitted; see fit_failures.csv")

# ---------------------------------------------------------------- 4. observed PIT + statistics
pit_seed <- function(m, what) seed_for(cfg$BASE_SEED, what, m)
g6 <- list(); stats_rows <- list(); rerand_rows <- list(); sample_means <- list()
for (m in models) {
  pb <- pit_bounds(Y, fits[[m]], fams[[m]], G_ALL[fits[[m]]$ok[G_ALL]], label = "G_ALL")   # G6
  g6[[m]] <- as.data.table(pb$g6)
  U <- randomise_pit(pb$a, pb$b, pit_seed(m, "pit_gall"))
  saveRDS(list(a = pb$a, b = pb$b, u = U), file.path(OUT, sprintf("pit_gall_%s.rds", m)))
  sa <- gof_all(U[A_SET, , drop = FALSE]); so <- gof_all(U)
  stats_rows[[m]] <- rbind(
    data.table::data.table(model = m, gene_set = "G_ALL_analysis_set", n_genes = length(A_SET), as.data.table(sa)),
    data.table::data.table(model = m, gene_set = "G_ALL_own_ok_set", n_genes = nrow(U), as.data.table(so)))
  sample_means[[m]] <- sample_qnorm_mean(U[A_SET, , drop = FALSE])
  if (isTRUE(cfg$RUN_RERAND))
    rerand_rows[[m]] <- cbind(model = m, rerandomise(pb$a[A_SET, , drop = FALSE], pb$b[A_SET, , drop = FALSE],
                                                     cfg$RERAND_R, pit_seed(m, "pit_gall")))
}
data.table::fwrite(data.table::rbindlist(g6), file.path(OUT, "g6_pit_asserts.csv"))
data.table::fwrite(data.table::rbindlist(stats_rows), file.path(OUT, "stats_gall_obs.csv"))
if (length(rerand_rows)) data.table::fwrite(data.table::rbindlist(rerand_rows), file.path(OUT, "rerandomisation_gall.csv"))
say("G6 PIT asserts passed for all models on G_ALL")

# ---------------------------------------------------------------- 5. G_BOOT reference
G_BOOT <- select_boot_genes(A_SET, Y, O, cfg$N_GENES_BOOT, seed_for(cfg$BASE_SEED, "gboot"))
writeLines(G_BOOT, file.path(OUT, "gboot_genes.txt"))
Yb0 <- Y[G_BOOT, , drop = FALSE]
Lb0 <- if (is.null(L)) NULL else L[G_BOOT, , drop = FALSE]
# The replicates re-estimate offsets on their own G_BOOT counts (D3 = reestimate),
# so the observed reference is processed the same way: offsets from the G_BOOT
# counts, not the G_ALL offsets. With OFFSETS_BOOT = fixed both use these.
Ob0 <- compute_offsets(Yb0, cfg$OFFSET_MODE, Lb0)
saveRDS(list(Y = Yb0, O = Ob0, X = X, lengths = Lb0, samples = ds$samples, design = ds$design),
        file.path(OUT, "gboot_data.rds"))
gs_obs <- apply_strata(Yb0, Ob0, X, G_BOOT)
ref <- list(); tobs <- list(); strata_obs <- list()
for (m in models) {
  f <- fit_dataset(Yb0, X, Ob0, fams[[m]], G_BOOT, cfg$N_CORES)
  saveRDS(f, file.path(OUT, sprintf("fits_gboot_%s.rds", m)))
  okg <- G_BOOT[f$ok]
  pb <- pit_bounds(Yb0, f, fams[[m]], okg, label = "G_BOOT")
  U <- randomise_pit(pb$a, pb$b, pit_seed(m, "pit_gboot"))
  saveRDS(list(a = pb$a, b = pb$b, u = U), file.path(OUT, sprintf("pit_gboot_%s.rds", m)))
  s <- gof_all(U)
  tobs[[m]] <- data.table::data.table(model = m, n_genes = length(okg), n_failed_ref = length(G_BOOT) - length(okg),
                                      as.data.table(s))
  strata_obs[[m]] <- cbind(model = m, strata_stats(U, gs_obs, grp))
  ref[[m]] <- list(fit = f, Y = Yb0[okg, , drop = FALSE], O = Ob0[okg, , drop = FALSE], X = X,
                   lengths = if (is.null(Lb0)) NULL else Lb0[okg, , drop = FALSE], genes = okg, obs_group = grp)
  say("G_BOOT %-8s: %d/%d OK, W2 = %.4g", m, length(okg), length(G_BOOT), s$W2)
}
data.table::fwrite(data.table::rbindlist(tobs), file.path(OUT, "stats_gboot_obs.csv"))
data.table::fwrite(data.table::rbindlist(strata_obs), file.path(OUT, "strata_gboot_obs.csv"))
saveRDS(sample_means, file.path(OUT, "sample_qnorm_means_obs.rds"))

# ---------------------------------------------------------------- 6. bootstrap
if (isTRUE(cfg$RUN_BOOTSTRAP)) {
  bsum <- list(); bstr <- list(); timing <- list()
  for (m in models) {
    t0 <- proc.time()[["elapsed"]]
    bt <- run_bootstrap(m, fams[[m]], ref[[m]], cfg, cfg$B, cfg$N_CORES)
    saveRDS(bt, file.path(OUT, sprintf("boot_%s.rds", m)))
    bsum[[m]] <- bootstrap_summary(bt, tobs[[m]])
    so <- data.table::copy(strata_obs[[m]])[, model := NULL]
    bstr[[m]] <- cbind(model = m, bootstrap_strata_summary(bt, so))
    timing[[m]] <- data.table::data.table(model = m, B = cfg$B, B_ok = nrow(bt$T_b), B_failed = nrow(bt$failed),
                                          wall_s = proc.time()[["elapsed"]] - t0,
                                          median_s_per_replicate = stats::median(bt$T_b$time_s),
                                          n_genes = length(ref[[m]]$genes), n_cores = cfg$N_CORES)
    if (nrow(bt$failed)) data.table::fwrite(bt$failed, file.path(OUT, sprintf("boot_failed_%s.csv", m)))
    say("bootstrap %-8s: %d/%d replicates, p_GOF(W2) = %.4g, %.1f s/replicate", m, nrow(bt$T_b), cfg$B,
        bsum[[m]][stat == "W2", p_GOF], timing[[m]]$median_s_per_replicate)
  }
  data.table::fwrite(data.table::rbindlist(bsum), file.path(OUT, "boot_summary.csv"))
  data.table::fwrite(data.table::rbindlist(bstr, fill = TRUE), file.path(OUT, "boot_strata_summary.csv"))
  data.table::fwrite(data.table::rbindlist(timing), file.path(OUT, "boot_timing.csv"))
}

# ---------------------------------------------------------------- 7. held-out score
if (isTRUE(cfg$RUN_HELDOUT)) {
  folds <- make_folds(ncol(Y), grp, cfg$K_FOLDS, seed_for(cfg$BASE_SEED, "folds"))
  ho <- list()
  for (m in models) {
    t0 <- proc.time()[["elapsed"]]
    ho[[m]] <- heldout_score(Yb0, X, Ob0, fams[[m]], G_BOOT, folds, cfg$N_CORES)
    saveRDS(ho[[m]], file.path(OUT, sprintf("heldout_%s.rds", m)))
    say("held-out %-8s: %d fold-fit failures (%.1f s)", m, sum(ho[[m]]$fold_failures, na.rm = TRUE),
        proc.time()[["elapsed"]] - t0)
  }
  hs <- heldout_summary(ho, "nb", seed = seed_for(cfg$BASE_SEED, "heldout_ci"))
  data.table::fwrite(hs, file.path(OUT, "heldout_summary.csv"))
  hpit <- data.table::rbindlist(lapply(models, function(m) {
    ok <- is.finite(ho[[m]]$a) & is.finite(ho[[m]]$b)
    u <- randomise_pit(ho[[m]]$a, ho[[m]]$b, pit_seed(m, "pit_heldout"))[ok]
    data.table::data.table(model = m, as.data.table(gof_all(u)))
  }))
  data.table::fwrite(hpit, file.path(OUT, "heldout_pit_stats.csv"))
}

write_provenance(OUT, cfg, extra = list(dataset = cfg$DATASET, label = ds$label, source = ds$source,
                                        n_samples = ncol(Y), n_genes_filtered = nrow(Y),
                                        n_G_ALL = length(G_ALL), n_analysis_set = length(A_SET),
                                        n_G_BOOT = length(G_BOOT),
                                        wall_seconds = round(proc.time()[["elapsed"]] - t_start)))
say("done: %s", OUT)
