#!/usr/bin/env Rscript
# 01_calibration_sim.R -- gate G5: does the whole test work when the truth is known?
#
#   Rscript Distribution_GOF/scripts/01_calibration_sim.R --G5_REF_DIR=<an 02_run_gof.R output dir>
#
# For truth = NB, then truth = PLN: G5_M datasets (default 40), each of G5_G genes
# (default 300) whose parameters are drawn (with replacement) from the REFERENCE
# run's fits of that model, at the reference run's n, design and library sizes.
# Each dataset then goes through the full procedure -- offsets re-estimated,
# fit, PIT, T_obs, parametric bootstrap with G5_B replicates (default 99) -- for
# the TRUE model and for Poisson.
#
# Pass criteria (brief section 9), operationalised as:
#   * true model: fraction of p_GOF(W2) < 0.05 whose binomial 95% CI includes or
#     lies below 0.05 (i.e. not demonstrably concentrated near 0); the ECDF of
#     p_GOF is written in full;
#   * Poisson: rejected (p_GOF(W2) < 0.05) in >= 95% of datasets.
# A run with M, G or B below the defaults is written as a SMOKE run and is never
# reported as the gate.
#
# Library term for simulation: per sample, the median over genes of the reference
# run's G_BOOT offsets (a pure library-size term; for OFFSET_MODE = tmm every gene
# has exactly this offset). Simulated datasets are integer, and offsets are
# re-estimated with OFFSET_MODE = "tmm" (no transcript lengths exist for them).

GOF_ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(Sys.getenv("GOF_ROOT", "Distribution_GOF"))
})
source(file.path(GOF_ROOT, "R", "setup.R"))
cfg <- GOF_CONFIG
print_gof_config(cfg)
if (isTRUE(cfg$DRY_RUN)) quit(save = "no", status = 0)
if (!nzchar(cfg$G5_REF_DIR) || !dir.exists(cfg$G5_REF_DIR))
  stop("G5 needs --G5_REF_DIR=<an 02_run_gof.R output dir> (the source of the NB/PLN parameters)")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")

full_gate <- cfg$G5_M >= 40L && cfg$G5_G >= 300L && cfg$G5_B >= 99L
OUT <- file.path(cfg$OUT_ROOT, if (full_gate) "G5_calibration" else "G5_calibration_SMOKE")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

gb <- readRDS(file.path(cfg$G5_REF_DIR, "gboot_data.rds"))
X <- gb$X; n <- nrow(X)
lib_off <- apply(gb$O, 2, stats::median)
grp <- design_group(gb$samples, gb$design)
fams <- gof_families(cfg, c("poisson", "nb", "pln"))
cfg_sim <- cfg; cfg_sim$OFFSET_MODE <- "tmm"

rows <- list()
for (truth in c("nb", "pln")) {
  ref <- readRDS(file.path(cfg$G5_REF_DIR, sprintf("fits_gall_%s.rds", truth)))
  pool <- names(ref$ok)[ref$ok]
  if (length(pool) < 10L) stop("reference run has fewer than 10 OK ", truth, " fits")
  for (k in seq_len(cfg$G5_M)) {
    t0 <- proc.time()[["elapsed"]]
    pick <- with_seed(seed_for(cfg$BASE_SEED, "G5pick", truth, k), sample(pool, cfg$G5_G, replace = TRUE))
    genes <- sprintf("sim%04d", seq_len(cfg$G5_G))
    Y <- with_seed(seed_for(cfg$BASE_SEED, "G5sim", truth, k), t(vapply(pick, function(g) {
      fr <- ref$par[[g]]
      par <- c(list(mu = fams[[truth]]$mu_fn(fr$beta, fr$extra, X, lib_off)), fr$extra)
      as.numeric(fams[[truth]]$rgen(par))
    }, numeric(n))))
    dimnames(Y) <- list(genes, rownames(X) %||% paste0("s", seq_len(n)))
    O <- compute_offsets(Y, "tmm")
    for (m in unique(c(truth, "poisson"))) {
      f <- fit_dataset(Y, X, O, fams[[m]], genes, cfg$N_CORES)
      okg <- genes[f$ok]
      pb <- pit_bounds(Y, f, fams[[m]], okg, label = sprintf("G5 %s k=%d", truth, k))
      U <- randomise_pit(pb$a, pb$b, seed_for(cfg$BASE_SEED, "G5pit", truth, k, m))
      tob <- gof_all(U)
      cfg_k <- cfg_sim; cfg_k$BASE_SEED <- seed_for(cfg$BASE_SEED, "G5boot", truth, k)
      bt <- run_bootstrap(m, fams[[m]], list(fit = f, Y = Y[okg, , drop = FALSE], O = O[okg, , drop = FALSE],
                                             X = X, lengths = NULL, genes = okg, obs_group = grp),
                          cfg_k, cfg$G5_B, cfg$N_CORES)
      bs <- bootstrap_summary(bt, tob)
      rows[[length(rows) + 1L]] <- data.table::data.table(truth = truth, dataset = k, model = m,
                                                           n_genes_ok = length(okg), bs)
    }
    say("G5 truth=%s dataset %d/%d (%.0f s)", truth, k, cfg$G5_M, proc.time()[["elapsed"]] - t0)
    data.table::fwrite(data.table::rbindlist(rows), file.path(OUT, "G5_per_dataset.csv"))   # checkpoint
  }
}
res <- data.table::rbindlist(rows)
w2 <- res[stat == "W2"]
summ <- w2[, {
  rej <- sum(p_GOF < 0.05); ci <- stats::binom.test(rej, .N)$conf.int
  .(M = .N, frac_p_below_0.05 = rej / .N, ci_lo = ci[1], ci_hi = ci[2],
    p_GOF_q10 = stats::quantile(p_GOF, 0.1), p_GOF_median = stats::median(p_GOF))
}, by = .(truth, model)]
summ[, role := ifelse(model == truth, "true model", "Poisson")]
summ[, pass := ifelse(role == "true model", ci_lo <= 0.05, frac_p_below_0.05 >= 0.95)]
data.table::fwrite(summ, file.path(OUT, "G5_summary.csv"))
ecdf_tab <- w2[, .(p = sort(p_GOF), ecdf = seq_len(.N) / .N), by = .(truth, model)]
data.table::fwrite(ecdf_tab, file.path(OUT, "G5_pGOF_ecdf.csv"))
write_provenance(OUT, cfg, extra = list(reference_dir = cfg$G5_REF_DIR, full_gate = full_gate))
print(summ)
if (!full_gate) cat("\nSMOKE RUN (M, G or B below the brief's defaults): this is NOT gate G5.\n")
if (full_gate && !all(summ$pass)) { cat("\nG5 FAILED.\n"); quit(save = "no", status = 1) }
