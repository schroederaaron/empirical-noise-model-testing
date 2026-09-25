#!/usr/bin/env Rscript
# 00_validate_families.R -- gates G1-G4 (brief section 9). HARD STOPS: exits with
# status 1 if G1, G2 or G3 fails; G4 has no threshold and is printed prominently.
#
#   Rscript Distribution_GOF/scripts/00_validate_families.R [--DATASET=ds.rds] [--N_TARGET=96]
#
# G1  pmf/cdf, on mu in {0.5, 5, 50, 500, 5e3, 5e4} x three dispersion levels:
#     |sum pmf - 1| < 1e-8 over the support (GP under-dispersed: mass reported);
#     F(q) - F(q-1) = p(q) to 1e-10 relative -- with the representational floor
#     64 eps F(q) that also governs G6 (see R/pit.R); checked on q = 0..1000 plus
#     2000 points spread over the support.
# G2  rgen: 1e5 draws per grid point vs the own cdf: max |ECDF - F| < 3 sqrt(1/1e5);
#     mean and variance within 4 standard errors.
# G3  own logpmf summed at the fitted parameters vs the engine's logLik
#     (glmmTMB families |d| < 1e-6); Poisson glm.fit vs glmmTMB; own CMP fitter vs
#     glmmTMB::compois on 200 genes with mu <= 50 (|dlogLik| < 1e-4, relative
#     d beta, nu < 1e-3); genes where the glmmTMB REFERENCE does not converge or errors
#     are excluded and listed (G3_excluded_reference_failures.csv); own PLN vs
#     glmmTMB (1|obs) Laplace -- REPORTED only.
# G4  parameter recovery: 200 genes per family at n = target n; median relative
#     bias of the dispersion parameter (no threshold).

GOF_ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(Sys.getenv("GOF_ROOT", "Distribution_GOF"))
})
source(file.path(GOF_ROOT, "R", "setup.R"))
cfg <- GOF_CONFIG
print_gof_config(cfg)
if (isTRUE(cfg$DRY_RUN)) quit(save = "no", status = 0)
OUT <- file.path(cfg$OUT_ROOT, "validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")
fams <- gof_families(cfg, ALL_MODELS)
EPS_F <- 64 * .Machine$double.eps
TINY  <- 64 * .Machine$double.xmin   # subnormal floor, as in R/pit.R

MUS <- c(0.5, 5, 50, 500, 5e3, 5e4)
# Three dispersion levels per family: over- to under-dispersion where the family has it.
DISP <- list(poisson = list(list()),
             nb      = lapply(c(0.5, 5, 50), function(s) list(size = s)),
             zinb    = lapply(c(0.5, 5, 50), function(s) list(size = s, pi = 0.1)),
             genpois = lapply(c(0.5, 2, 20), function(p) list(phi = p)),     # 0.5 = under-dispersed
             pln     = lapply(c(0.1, 0.5, 1), function(s) list(sigma = s)),
             cmp     = lapply(c(0.1, 1, 3), function(v) list(nu = v)))      # 3 = under-dispersed
disp_label <- function(d) if (!length(d)) "-" else paste(names(d), unlist(d), sep = "=", collapse = ",")

#' Upper end of the support to sum over: mean + 60 sd, the lognormal 7-sigma point for PLN,
#' and for GP at least the family's own (tail-widened) window.
support_max <- function(m, par) {
  mo <- fams[[m]]$moments(par)
  hi <- ceiling(mo$mean + 60 * sqrt(mo$var)) + 50
  if (m == "pln") hi <- max(hi, ceiling(par$mu * exp(7 * par$sigma)) + 50)
  if (m == "genpois") hi <- max(hi, .gp_window(par$mu, par$phi)[["hi"]])     # GP's right tail can be heavy
  if (m == "genpois" && par$phi < 1) hi <- min(hi, ceiling(-par$mu / sqrt(par$phi) / (1 - 1 / sqrt(par$phi))))
  hi
}

# ---------------------------------------------------------------- G1
say("G1: pmf / cdf")
GRID <- do.call(c, lapply(ALL_MODELS, function(m) do.call(c, lapply(DISP[[m]], function(d)
  lapply(MUS, function(mu) list(m = m, d = d, mu = mu))))))
g1 <- parallel::mclapply(GRID, function(gp) {
  m <- gp$m; d <- gp$d; mu <- gp$mu
  par <- c(list(mu = mu), d)
  f <- fams[[m]]
  Q <- support_max(m, par)
  mass <- 0; chunk <- 2e5
  for (s in seq(0, Q, by = chunk)) {
    q <- s:min(Q, s + chunk - 1)
    mass <- mass + sum(exp(f$logpmf(q, c(list(mu = rep(mu, length(q))), d))))
  }
  qs <- unique(c(0:min(Q, 1000), round(seq(0, Q, length.out = 2000))))
  pq <- c(list(mu = rep(mu, length(qs))), d)
  Fq <- f$cdf(qs, pq); Fm <- f$cdf(qs - 1, pq); pp <- exp(f$logpmf(qs, pq))
  dev <- abs(Fq - Fm - pp); tol <- 1e-10 * pp + EPS_F * Fq + TINY
  ratio <- ifelse(dev == 0, 0, dev / tol)
  under <- m == "genpois" && d$phi < 1
  # Under-dispersed GP whose total mass exceeds 1 is not a distribution (no cdf
  # exists); such a point is reported as invalid, with pass = NA, never as a pass.
  invalid <- under && mass > 1 + 1e-8
  data.table::data.table(
    family = m, mu = mu, dispersion = disp_label(d), support_max = Q, mass = mass,
    mass_err = abs(mass - 1), n_q_checked = length(qs),
    max_inc_dev_over_tol = max(ratio), n_floor_binding = sum(dev > 1e-10 * pp),
    gp_underdispersed = under, invalid_distribution = invalid,
    pass = if (invalid) NA else (under || abs(mass - 1) < 1e-8) && max(ratio) <= 1)
}, mc.cores = cfg$N_CORES, mc.preschedule = FALSE)
.chk <- function(l, what) { bad <- !vapply(l, data.table::is.data.table, NA)
  if (any(bad)) stop(what, " crashed on ", sum(bad), " grid points: ", as.character(l[[which(bad)[1]]])); l }
g1 <- data.table::rbindlist(.chk(g1, "G1"))
data.table::fwrite(g1, file.path(OUT, "G1_pmf_cdf.csv"))
say("G1: %d pass, %d fail, %d invalid (GP mass > 1)", sum(g1$pass, na.rm = TRUE), sum(!g1$pass, na.rm = TRUE), sum(g1$invalid_distribution))

# ---------------------------------------------------------------- G2
say("G2: rgen (1e5 draws per grid point)")
NDRAW <- 1e5
g2 <- parallel::mclapply(GRID, function(gp) {
  m <- gp$m; d <- gp$d; mu <- gp$mu
  f <- fams[[m]]
  par <- c(list(mu = rep(mu, NDRAW)), d)
  x <- with_seed(seed_for(cfg$BASE_SEED, "G2", m, disp_label(d), mu), f$rgen(par))
  ux <- sort(unique(x))
  p1 <- c(list(mu = rep(mu, length(ux))), d)
  Fx <- f$cdf(ux, p1); Fxm <- f$cdf(ux - 1, p1)
  en <- stats::ecdf(x)
  Dmax <- max(abs(en(ux) - Fx), abs(c(0, en(ux)[-length(ux)]) - Fxm))
  mo <- f$moments(c(list(mu = mu), d))
  se_m <- sqrt(mo$var / NDRAW)
  m4 <- mean((x - mean(x))^4); s2 <- stats::var(x)
  se_v <- sqrt(max(m4 - s2^2, 0) / NDRAW)
  under <- m == "genpois" && d$phi < 1
  invalid <- under && sum(exp(gp_logpmf(0:support_max(m, c(list(mu = mu), d)), mu, d$phi))) > 1 + 1e-8
  data.table::data.table(
    family = m, mu = mu, dispersion = disp_label(d), max_ecdf_dev = Dmax, invalid_distribution = invalid,
    mean_z = (mean(x) - mo$mean) / se_m, var_z = (s2 - mo$var) / se_v,
    gp_underdispersed = under,
    pass = if (invalid) NA else Dmax < 3 * sqrt(1 / NDRAW) && abs((mean(x) - mo$mean) / se_m) < 4 &&
           abs((s2 - mo$var) / se_v) < 4)
}, mc.cores = cfg$N_CORES, mc.preschedule = FALSE)
g2 <- data.table::rbindlist(.chk(g2, "G2"))
data.table::fwrite(g2, file.path(OUT, "G2_rgen.csv"))
say("G2: %d pass, %d fail, %d invalid (GP mass > 1)", sum(g2$pass, na.rm = TRUE), sum(!g2$pass, na.rm = TRUE), sum(g2$invalid_distribution))

# ---------------------------------------------------------------- design for G3/G4
if (nzchar(cfg$DATASET)) {
  ds <- read_dataset(cfg$DATASET)
  Xd <- stats::model.matrix(ds$design, ds$samples); n_t <- nrow(Xd)
  lib <- colSums(ds$counts); design_src <- paste("dataset", ds$label)
} else {
  n_t <- cfg$N_TARGET
  Xd <- stats::model.matrix(~ condition, data.frame(condition = factor(rep(c("A", "B"), length.out = n_t))))
  lib <- with_seed(seed_for(cfg$BASE_SEED, "G34lib"), stats::runif(n_t, 0.5, 2) * 1e7)
  design_src <- sprintf("synthetic 2-condition design, n = %d (no DATASET given)", n_t)
}
offv <- log(lib)
say("G3/G4 design: %s", design_src)

# ---------------------------------------------------------------- G3
say("G3: logLik consistency")
tmb_ll <- function(y, X, off, family, zi = FALSE, re = FALSE) {
  d <- data.frame(y = y, off = off, X); names(d)[-(1:2)] <- xs <- paste0("x", seq_len(ncol(X)))
  d$obs <- factor(seq_along(y))
  f <- stats::as.formula(paste("y ~ 0 +", paste(xs, collapse = "+"), "+ offset(off)", if (re) "+ (1|obs)" else ""))
  suppressWarnings(glmmTMB::glmmTMB(f, data = d, family = family, ziformula = if (zi) ~1 else ~0,
                                    control = glmmTMB::glmmTMBControl(parallel = 1L)))
}
g3 <- list()
mu_levels <- c(3, 30, 300, 3000)
for (m in c("nb", "zinb", "genpois")) for (lv in mu_levels) for (r in 1:5) {
  y <- with_seed(seed_for(cfg$BASE_SEED, "G3", m, lv, r),
                 stats::rnbinom(n_t, mu = lv * lib / mean(lib), size = 4))
  fr <- fams[[m]]$fit(y, Xd, offv)
  if (!isTRUE(fr$conv_ok)) { g3[[length(g3) + 1L]] <- data.table::data.table(check = m, mu = lv, rep = r, delta = NA_real_, pass = FALSE, note = fr$conv_msg); next }
  own <- sum(fams[[m]]$logpmf(y, fit_par(fr)))
  g3[[length(g3) + 1L]] <- data.table::data.table(check = paste0(m, ": own logpmf vs logLik"), mu = lv, rep = r,
                                                   delta = abs(own - fr$loglik), pass = abs(own - fr$loglik) < 1e-6, note = "")
}
for (lv in mu_levels) for (r in 1:5) {
  y <- with_seed(seed_for(cfg$BASE_SEED, "G3pois", lv, r), stats::rpois(n_t, lv * lib / mean(lib)))
  a <- fams$poisson$fit(y, Xd, offv); b <- tmb_ll(y, Xd, offv, stats::poisson())
  d <- abs(a$loglik - as.numeric(stats::logLik(b)))
  g3[[length(g3) + 1L]] <- data.table::data.table(check = "poisson: glm.fit vs glmmTMB logLik", mu = lv, rep = r,
                                                   delta = d, pass = d < 1e-6, note = "")
}
# CMP: own fitter vs glmmTMB::compois, 200 genes with mu <= 50, intercept-only
# design (a near-zero contrast coefficient would make "relative delta beta" undefined).
X1 <- matrix(1, n_t, 1, dimnames = list(NULL, "(Intercept)"))
cmp_rows <- parallel::mclapply(1:200, function(r) {
  mu <- with_seed(seed_for(cfg$BASE_SEED, "G3cmpmu", r), stats::runif(1, 1, 50))
  nu <- with_seed(seed_for(cfg$BASE_SEED, "G3cmpnu", r), exp(stats::runif(1, log(0.3), log(3))))
  mui <- mu * lib / mean(lib)
  if (max(mui) > 50) mui <- mui * 50 / max(mui)
  y <- with_seed(seed_for(cfg$BASE_SEED, "G3cmpy", r), fams$cmp$rgen(list(mu = mui, nu = nu)))
  own <- tryCatch(fams$cmp$fit(y, X1, log(mui / mean(mui) * mean(lib))), error = function(e) NULL)
  tm <- tryCatch(tmb_ll(y, X1, log(mui / mean(mui) * mean(lib)), glmmTMB::compois()), error = function(e) NULL)
  # Reference-side failures are EXCLUDED from the gate (decision 24.09.2026, Aaron):
  # glmmTMB::compois non-convergence (fit$convergence != 0) or an error is a failure
  # of the cross-check's reference, not of the own fitter. Such rows get
  # excluded = TRUE and pass = NA, and are listed in the verdict. A failure of the
  # OWN fitter is never excluded.
  tmb_status <- if (is.null(tm)) "error" else if (tm$fit$convergence != 0)
    paste0("not converged (code ", tm$fit$convergence, ")") else "converged"
  ll_t <- if (!is.null(tm) && tm$fit$convergence == 0) as.numeric(stats::logLik(tm)) else NA_real_
  own_ok <- !is.null(own) && isTRUE(own$conv_ok)
  base <- data.table::data.table(check = "cmp: own vs glmmTMB::compois", mu = mu, rep = r,
                                 own_status = if (own_ok) "converged" else "failed",
                                 own_loglik = if (is.null(own)) NA_real_ else own$loglik,
                                 tmb_status = tmb_status, tmb_loglik = ll_t)
  if (!own_ok)
    return(cbind(base, delta = NA_real_, excluded = FALSE, pass = FALSE, note = "own fitter failed"))
  if (tmb_status != "converged")
    return(cbind(base, delta = NA_real_, excluded = TRUE, pass = NA,
                 note = sprintf("excluded: glmmTMB %s (glmmTMB nu = %s)", tmb_status,
                                if (is.null(tm)) "NA" else format(1 / stats::sigma(tm), digits = 5))))
  b_t <- unname(glmmTMB::fixef(tm)$cond); nu_t <- 1 / stats::sigma(tm)
  d_ll <- abs(own$loglik - ll_t)
  rb <- max(abs(own$beta - b_t) / abs(b_t)); rn <- abs(own$extra$nu - nu_t) / nu_t
  cbind(base, delta = d_ll, excluded = FALSE, pass = d_ll < 1e-4 && rb < 1e-3 && rn < 1e-3,
        note = sprintf("rel dbeta=%.2e rel dnu=%.2e", rb, rn))
}, mc.cores = cfg$N_CORES)
g3 <- data.table::rbindlist(c(g3, cmp_rows), fill = TRUE)
if (!"excluded" %in% names(g3)) g3[, excluded := FALSE]
g3[is.na(excluded), excluded := FALSE]
# PLN: own vs Laplace (reported only, not a gate: Laplace differs from quadrature).
pln_rows <- parallel::mclapply(1:200, function(r) {
  sg <- with_seed(seed_for(cfg$BASE_SEED, "G3plns", r), stats::runif(1, 0.1, 1))
  mu <- with_seed(seed_for(cfg$BASE_SEED, "G3plnm", r), exp(stats::runif(1, log(2), log(50))))
  mui <- mu * lib / mean(lib)
  y <- with_seed(seed_for(cfg$BASE_SEED, "G3plny", r), fams$pln$rgen(list(mu = mui, sigma = sg)))
  own <- fams$pln$fit(y, Xd, offv)
  tm <- tryCatch(tmb_ll(y, Xd, offv, stats::poisson(), re = TRUE), error = function(e) NULL)
  s_t <- if (is.null(tm)) NA_real_ else sqrt(unlist(glmmTMB::VarCorr(tm)$cond))
  data.table::data.table(sigma_true = sg, sigma_own = own$extra$sigma, sigma_laplace = unname(s_t),
                         own_ok = fit_ok(own))
}, mc.cores = cfg$N_CORES)
pln_cmp <- data.table::rbindlist(pln_rows)
data.table::fwrite(g3, file.path(OUT, "G3_loglik.csv"))
data.table::fwrite(pln_cmp, file.path(OUT, "G3_pln_vs_laplace_REPORTED_ONLY.csv"))
g3_pass <- all(g3$pass, na.rm = TRUE)
g3_excl <- g3[excluded == TRUE]
say("G3: %d pass, %d fail, %d excluded (glmmTMB::compois reference did not converge)",
    sum(g3$pass, na.rm = TRUE), sum(!g3$pass, na.rm = TRUE), nrow(g3_excl))
if (nrow(g3_excl)) {
  data.table::fwrite(g3_excl, file.path(OUT, "G3_excluded_reference_failures.csv"))
  cat("G3 excluded genes (glmmTMB reference failures):\n"); print(g3_excl[, .(rep, mu, tmb_status, own_loglik, note)])
}

# ---------------------------------------------------------------- G4
say("G4: parameter recovery (200 genes per family, n = %d)", n_t)
truth <- list(nb = list(size = 5), zinb = list(size = 5, pi = 0.1), genpois = list(phi = 4),
              pln = list(sigma = 0.4), cmp = list(nu = 0.3))
disp_key <- c(nb = "size", zinb = "size", genpois = "phi", pln = "sigma", cmp = "nu")
g4 <- list()
for (m in names(truth)) {
  est <- parallel::mclapply(1:200, function(r) {
    mu <- with_seed(seed_for(cfg$BASE_SEED, "G4mu", m, r), exp(stats::runif(1, log(5), log(500))))
    par <- c(list(mu = mu * lib / mean(lib)), truth[[m]])
    y <- with_seed(seed_for(cfg$BASE_SEED, "G4y", m, r), fams[[m]]$rgen(par))
    fr <- fit_one_gene(fams[[m]], y, Xd, offv)
    if (fr$ok) c(fr$extra[[disp_key[[m]]]], if (m == "zinb") fr$extra$pi else NA) else c(NA, NA)
  }, mc.cores = cfg$N_CORES)
  e <- do.call(rbind, est)
  tv <- truth[[m]][[disp_key[[m]]]]
  g4[[m]] <- data.table::data.table(family = m, parameter = disp_key[[m]], true = tv,
                                    n_ok = sum(is.finite(e[, 1])),
                                    median_estimate = stats::median(e[, 1], na.rm = TRUE),
                                    median_relative_bias = stats::median((e[, 1] - tv) / tv, na.rm = TRUE))
  if (m == "zinb") g4[["zinb_pi"]] <- data.table::data.table(family = m, parameter = "pi", true = 0.1,
                                    n_ok = sum(is.finite(e[, 2])), median_estimate = stats::median(e[, 2], na.rm = TRUE),
                                    median_relative_bias = stats::median((e[, 2] - 0.1) / 0.1, na.rm = TRUE))
}
g4 <- data.table::rbindlist(g4)
data.table::fwrite(g4, file.path(OUT, "G4_recovery.csv"))
cat("\n================ G4: PARAMETER RECOVERY (no threshold) ================\n")
print(g4, digits = 4)
cat("=======================================================================\n\n")

# ---------------------------------------------------------------- verdict
verdict <- data.table::data.table(gate = c("G1", "G2", "G3"),
                                  passed = c(all(g1$pass, na.rm = TRUE), all(g2$pass, na.rm = TRUE), g3_pass),
                                  n_checks = c(nrow(g1), nrow(g2), nrow(g3)),
                                  n_failed = c(sum(!g1$pass, na.rm = TRUE), sum(!g2$pass, na.rm = TRUE), sum(!g3$pass, na.rm = TRUE)),
                                  n_invalid_distribution = c(sum(g1$invalid_distribution), sum(g2$invalid_distribution), 0L),
                                  n_excluded_reference_failures = c(0L, 0L, nrow(g3_excl)))
data.table::fwrite(verdict, file.path(OUT, "verdict.csv"))
write_provenance(OUT, cfg, extra = list(g34_design = design_src))
print(verdict)
if (!all(verdict$passed)) {
  cat("\nGATE FAILURE -- see the failing rows:\n")
  if (!all(g1$pass, na.rm = TRUE)) print(g1[pass == FALSE])
  if (!all(g2$pass, na.rm = TRUE)) print(g2[pass == FALSE])
  if (!g3_pass) print(g3[pass == FALSE])
  quit(save = "no", status = 1)
}
say("G1-G3 passed. Report: %s", OUT)
