# preprocess.R -- identical for all six models (brief section 4).

#' Fractional counts (decision D1). Integer input passes through unchanged under
#' every mode; the audit, not an assumption, decides whether anything happens.
#'   error      : stop if any entry is non-integer (the default)
#'   stochastic : floor(y) + Bernoulli(y - floor(y)), ONCE, under its own seed.
#'                Unbiased in the mean; adds <= 0.25 variance per entry.
#'   floor      : floor(y). Sensitivity arm only; biased downward.
apply_rounding <- function(Y, mode, seed) {
  nonint <- abs(Y - round(Y)) > 1e-8
  if (!any(nonint)) return(list(Y = round(Y), n_changed = 0L))
  if (mode == "error")
    stop(sum(nonint), " non-integer counts (", signif(100 * mean(nonint), 3), "% of entries); ",
         "ROUNDING = 'error'. Choose 'stochastic' or 'floor' explicitly (decision D1).")
  fl <- floor(Y)
  out <- if (mode == "floor") fl else
    with_seed(seed, fl + (stats::runif(length(Y)) < (Y - fl)))
  out[!nonint] <- round(Y[!nonint])          # exact integers stay exactly what they were
  dim(out) <- dim(Y); dimnames(out) <- dimnames(Y)
  list(Y = out, n_changed = sum(out != Y))
}

#' The factor filterByExpr groups by: the interaction of the design's
#' factor/character variables, or NULL for a single-cohort design (~ 1).
design_group <- function(samples, design) {
  vars <- all.vars(design)
  fac <- vars[vapply(vars, function(v) is.factor(samples[[v]]) || is.character(samples[[v]]), NA)]
  if (!length(fac)) return(NULL)
  droplevels(interaction(samples[fac], drop = TRUE))
}

#' edgeR::filterByExpr on the OBSERVED data only (brief 4.1). The resulting gene
#' set is then fixed for every model, bootstrap replicate and the held-out score.
filter_genes <- function(Y, group) {
  keep <- edgeR::filterByExpr(edgeR::DGEList(counts = Y), group = group)
  list(keep = keep, n_before = nrow(Y), n_after = sum(keep))
}

#' Offsets o_gi on the log scale, as a genes x samples matrix (decision D2).
#'   tmm        : log(lib.size * norm.factors) from normLibSizes on these counts.
#'   tmm_length : the tximport-vignette recipe, which in tximport 1.40 / edgeR 4.10
#'                is edgeR::DGEListFromTximport() (offset.prior = centred log
#'                transcript length) followed by normLibSizes(); getOffset() then
#'                returns offset.prior + log(lib.size * norm.factors).
compute_offsets <- function(Y, mode, lengths = NULL) {
  y <- if (mode == "tmm_length") {
    if (is.null(lengths)) stop("OFFSET_MODE = 'tmm_length' needs a dataset with $lengths (tximport)")
    edgeR::DGEListFromTximport(list(counts = Y, length = lengths[rownames(Y), colnames(Y), drop = FALSE],
                                    countsFromAbundance = "no"))
  } else edgeR::DGEList(counts = Y)
  y <- edgeR::normLibSizes(y)
  o <- edgeR::getOffset(y)
  if (is.null(dim(o))) o <- matrix(o, nrow(Y), ncol(Y), byrow = TRUE)
  dimnames(o) <- dimnames(Y)
  o
}

#' log2 mean CPM per gene, using the offsets as the library-size measure.
log2_mean_cpm <- function(Y, O) log2(rowMeans(Y / exp(O) * 1e6) + 1e-8)

#' Rank-based quantile bins (ties broken by position, so bins are equal-sized).
rank_bins <- function(x, k, prefix) {
  b <- ceiling(k * rank(x, ties.method = "first") / length(x))
  factor(sprintf("%s%d", prefix, b), levels = sprintf("%s%d", prefix, seq_len(k)))
}

#' G_BOOT: a fixed random subset of `genes`, stratified by mean-expression decile
#' (brief 4.4), proportional allocation, under its own seed.
select_boot_genes <- function(genes, Y, O, n, seed) {
  if (n >= length(genes)) return(genes)
  dec <- rank_bins(log2_mean_cpm(Y[genes, , drop = FALSE], O[genes, , drop = FALSE]), 10L, "D")
  alloc <- round(n * table(dec) / length(genes))
  alloc[which.max(alloc)] <- alloc[which.max(alloc)] + (n - sum(alloc))      # exact total
  pick <- with_seed(seed, unlist(lapply(names(alloc), function(d) {
    pool <- genes[dec == d]
    pool[sample.int(length(pool), min(length(pool), alloc[[d]]))]
  })))
  genes[genes %in% pick]                                                      # keep input order
}
