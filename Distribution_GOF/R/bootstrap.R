# bootstrap.R -- parametric bootstrap (brief 8.1).
#
# Parallel over replicates b; genes run serially inside each worker (the choice
# the brief asks to document: one replicate is the natural unit of work, and it
# keeps every replicate's RNG stream a pure function of (model, b)).
#
# Stream of replicate b of model m: set.seed(seed_for(BASE_SEED, "boot", m, b))
# under L'Ecuyer-CMRG; its PIT draw uses seed_for(BASE_SEED, "bootpit", m, b). So a
# single replicate can be rerun on its own with run_one_replicate().
#
# KNOWN LIMIT (also in the README and the table caption): genes are simulated
# independently given their fitted parameters. Gene-gene correlation and shared
# sample-level effects are not in the null; if they are present in the real data,
# p_GOF is anti-conservative. At RNA-seq N, p_GOF = 1/(B+1) for every model is the
# expected outcome, not a finding.

#' One bootstrap replicate. `ref` bundles the reference fit and data (see run_bootstrap).
run_one_replicate <- function(b, model, fam, ref, cfg) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(seed_for(cfg$BASE_SEED, "boot", model, b))
  genes <- ref$genes
  Yb <- t(vapply(genes, function(g) as.numeric(fam$rgen(fit_par(ref$fit$par[[g]]))),
                 numeric(ncol(ref$Y))))
  dimnames(Yb) <- list(genes, colnames(ref$Y))
  Ob <- if (cfg$OFFSETS_BOOT == "fixed") ref$O else
    compute_offsets(Yb, cfg$OFFSET_MODE, ref$lengths)
  fb <- fit_dataset(Yb, ref$X, Ob, fam, genes, n_cores = 1L)
  okg <- genes[fb$ok]
  pb <- pit_bounds(Yb, fb, fam, okg, label = paste0("boot", b))
  U <- randomise_pit(pb$a, pb$b, seed_for(cfg$BASE_SEED, "bootpit", model, b))
  s <- gof_all(U)
  gs <- apply_strata(Yb, Ob, ref$X, okg)
  list(b = b, W2 = s$W2, A2 = s$A2, D = s$D, N = s$N,
       drop_frac = 1 - length(okg) / length(genes),
       strata = strata_stats(U, gs, ref$obs_group),
       hist = pit_hist(U),
       sample_mean = sample_qnorm_mean(U),
       time_s = proc.time()[["elapsed"]] - t0)
}

#' Full bootstrap for one model. ref = list(fit, Y, O, X, lengths, genes, obs_group),
#' where `fit` is the model's fit on `genes` with offsets O computed by the same
#' procedure the replicates use.
run_bootstrap <- function(model, fam, ref, cfg, B = cfg$B, n_cores = cfg$N_CORES) {
  one <- function(b) tryCatch(run_one_replicate(b, model, fam, ref, cfg),
                              error = function(e) list(b = b, error = conditionMessage(e)))
  res <- parallel::mclapply(seq_len(B), one, mc.cores = n_cores, mc.preschedule = FALSE)
  failed <- vapply(res, function(r) !is.list(r) || !is.null(r$error), NA)
  ok <- res[!failed]
  list(
    model = model,
    T_b = data.table::rbindlist(lapply(ok, function(r)
      data.table::data.table(b = r$b, W2 = r$W2, A2 = r$A2, D = r$D, N = r$N,
                             drop_frac = r$drop_frac, time_s = r$time_s))),
    strata = data.table::rbindlist(lapply(ok, function(r) cbind(b = r$b, r$strata))),
    hist = do.call(rbind, lapply(ok, `[[`, "hist")),
    sample_mean = do.call(rbind, lapply(ok, `[[`, "sample_mean")),
    failed = data.table::data.table(
      b = which(failed),
      error = vapply(res[failed], function(r) if (is.list(r)) r$error else as.character(r), "")))
}

#' p_GOF = (1 + #{T_b >= T_obs}) / (B + 1), over the replicates that completed.
#' excess = T_obs / median(T_b); z = (T_obs - mean T_b) / sd(T_b). Descriptive
#' cross-model comparisons, NOT a model-selection criterion.
bootstrap_summary <- function(boot, T_obs) {
  Tb <- boot$T_b
  Bok <- nrow(Tb)
  data.table::rbindlist(lapply(c("W2", "A2", "D"), function(st) {
    tb <- Tb[[st]]; to <- T_obs[[st]]
    data.table::data.table(model = boot$model, stat = st, T_obs = to, B_ok = Bok,
                           B_failed = nrow(boot$failed),
                           p_GOF = (1 + sum(tb >= to)) / (Bok + 1),
                           excess = to / stats::median(tb),
                           z = (to - mean(tb)) / stats::sd(tb),
                           max_drop_frac = max(Tb$drop_frac),
                           n_rep_drop_gt_1pct = sum(Tb$drop_frac > 0.01))
  }))
}

#' Per-stratum p_GOF against the replicates' per-stratum statistics.
bootstrap_strata_summary <- function(boot, strata_obs) {
  if (!nrow(boot$strata)) return(NULL)
  m <- merge(strata_obs, boot$strata, by = c("rule", "stratum"), suffixes = c("_obs", "_b"))
  m[, .(N_obs = N_obs[1], W2_obs = W2_obs[1],
        p_GOF_W2 = (1 + sum(W2_b >= W2_obs[1])) / (.N + 1), excess_W2 = W2_obs[1] / stats::median(W2_b),
        p_GOF_A2 = (1 + sum(A2_b >= A2_obs[1])) / (.N + 1),
        p_GOF_D  = (1 + sum(D_b  >= D_obs[1]))  / (.N + 1), B_ok = .N),
    by = .(rule, stratum)]
}
