# fit_dataset.R -- per-gene loop, error capture, parallel (brief section 6).

#' A fit counts as OK when the optimiser converged (convergence == 0), loglik is
#' finite, and every parameter is finite. pdHess is not a criterion.
fit_ok <- function(fr)
  isTRUE(fr$conv_ok) && is.finite(fr$loglik) && all(is.finite(fr$beta)) &&
    all(is.finite(unlist(fr$extra))) && all(is.finite(fr$mu)) && all(fr$mu > 0)

fit_one_gene <- function(fam, y, X, off) {
  t0 <- proc.time()[["elapsed"]]
  fr <- tryCatch(fam$fit(y, X, off),
                 error = function(e) list(beta = NA_real_, extra = list(), mu = NA_real_,
                                          loglik = NA_real_, conv_ok = FALSE,
                                          conv_msg = paste("ERROR:", conditionMessage(e)), engine = NA_character_))
  fr$time_s <- proc.time()[["elapsed"]] - t0
  fr$ok <- fit_ok(fr)
  fr
}

#' Fit `fam` to every gene in `genes`. n_cores = 1 runs serially (used inside
#' bootstrap workers, which are themselves the parallel unit).
fit_dataset <- function(Y, X, O, fam, genes, n_cores = 1L) {
  one <- function(g) fit_one_gene(fam, as.numeric(Y[g, ]), X, as.numeric(O[g, ]))
  res <- if (n_cores > 1L && length(genes) > 1L)
    parallel::mclapply(genes, one, mc.cores = n_cores, mc.preschedule = TRUE)
  else lapply(genes, one)
  bad <- !vapply(res, is.list, NA)                       # a worker that died outright
  if (any(bad)) res[bad] <- lapply(which(bad), function(i)
    list(beta = NA_real_, extra = list(), mu = NA_real_, loglik = NA_real_, conv_ok = FALSE,
         conv_msg = paste("ERROR: worker failed:", as.character(res[[i]])), engine = NA_character_,
         time_s = NA_real_, ok = FALSE))
  names(res) <- genes
  diag <- data.table::data.table(
    gene     = genes,
    loglik   = vapply(res, function(r) as.numeric(r$loglik)[1], 0),
    conv_ok  = vapply(res, function(r) isTRUE(r$conv_ok), NA),
    ok       = vapply(res, function(r) isTRUE(r$ok), NA),
    conv_msg = vapply(res, function(r) substr(paste(r$conv_msg, collapse = " "), 1, 300), ""),
    engine   = vapply(res, function(r) as.character(r$engine)[1], ""),
    gp_mass_min = vapply(res, function(r) if (is.null(r$gp_mass_min)) NA_real_ else r$gp_mass_min, 0),
    time_s   = vapply(res, function(r) as.numeric(r$time_s)[1], 0))
  list(par = res, diag = diag, ok = stats::setNames(diag$ok, genes))
}

#' Flat per-gene parameter table (for saving / G4 / G5 parameter draws).
params_table <- function(fitres) {
  data.table::rbindlist(lapply(names(fitres$par), function(g) {
    r <- fitres$par[[g]]
    if (!isTRUE(r$ok)) return(NULL)
    b <- as.list(stats::setNames(r$beta, paste0("beta", seq_along(r$beta))))
    data.table::as.data.table(c(list(gene = g), b, r$extra, list(loglik = r$loglik)))
  }), fill = TRUE)
}
