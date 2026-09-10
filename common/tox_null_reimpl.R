# tox_null_reimpl.R
# -----------------------------------------------------------------------------
# A faithful R re-implementation of the `own`-path of tox_noise_model_exact.F90,
# instrumented with every per-gene intermediate the Fortran does NOT return.
#
# WHY THIS EXISTS
# ---------------
# The Fortran pipeline returns only `pvalues_own` and three neighbourhood SIZES.
# Every quantity needed to explain a calibration failure -- the width, the SHAPE
# and the composition of the residual pool that forms the null, the gene's own
# variance, which neighbour genes were pooled, where the expansion stopped and
# why -- is computed inside `gather_residuals_helper` and discarded. Adding those
# outputs to the Fortran means an ABI change and a rebuild in every consumer.
# Re-implementing the same ~150 lines in R costs one file, needs no rebuild, and
# lets every diagnostic below be written as ordinary R.
#
# The price is that a re-implementation can drift from the thing it models, so
# `validate_against_fortran()` is not optional: it runs both implementations on
# the same input and reports the max absolute p-value difference. Treat any
# result from this file as void until that check passes.
#
# WHAT IS REPRODUCED (tox_noise_model_exact.F90, branch 145-empirical-noise-model)
#   prepare_sorted_data_helper  -- mean-sort, Bessel-corrected centred residuals,
#                                  Frechet-mean centring in log2 space
#   find_closest_helper / choose_index / gather_residuals_helper
#                               -- adaptive kNN pool growth in GENE units
#   trim_pool_tails_helper      -- symmetric quantile trim (raw only)
#   1/sqrt(n_rep) scaling       -- individual-residual -> mean-difference scale
#   compute_pvalue_helper       -- exact pairwise tail count, add-one corrected
#
# NOT reproduced: the bootstrap variant (tox_noise_model.F90), the family /
# ortholog comparisons (retired), max_pool_size truncation semantics beyond the
# plain cap.
# -----------------------------------------------------------------------------

NOISE_LOG_OFFSET <- 1.0   # matches the Fortran parameter of the same name

# ---------------------------------------------------------------- sorted data

#' Mean-sort genes and build the Bessel-corrected, centred residual matrix.
#'
#' Mirrors `prepare_sorted_data_helper`. Residuals are
#'   raw : sqrt(n/(n-1)) * (x_ig - mu_g)
#'   log : sqrt(n/(n-1)) * (log2(x_ig + c) - ghat_g),  ghat_g = mean_i log2(x_ig + c)
#' i.e. centred on the FRECHET mean of log space, not on log2(mu_g + c), so the
#' residuals sum to exactly zero (Jensen).
#'
#' @param replicates samples x genes matrix in PRE-LOG space (raw TPM)
#' @param norm_method 0 = linear residuals, non-zero = log2-space residuals
#' @return list(order, means_sorted, resid, n_samples, gene_ids) where `resid` is
#'   n_samples x n_genes with columns in mean-sorted order, and `order` maps a
#'   sorted slot back to the original column index.
tox_prepare_sorted <- function(replicates, norm_method) {
  stopifnot(is.matrix(replicates))
  n <- nrow(replicates)
  if (n < 2L) stop("prepare_sorted_data: need >= 2 replicates (Bessel correction).")
  means <- colMeans(replicates)
  ord   <- order(means)                      # ascending, ties in first-seen order
  x     <- replicates[, ord, drop = FALSE]
  bessel <- sqrt(n / (n - 1))
  if (norm_method == 0L) {
    centre <- colMeans(x)
  } else {
    x      <- log2(pmax(x, 0) + NOISE_LOG_OFFSET)
    centre <- colMeans(x)                    # Frechet mean in log2 space
  }
  resid <- bessel * sweep(x, 2L, centre, "-")
  list(order = ord, means_sorted = means[ord], resid = resid,
       n_samples = n, gene_ids = colnames(replicates)[ord],
       centre_sorted = centre, norm_method = norm_method)
}

# ---------------------------------------------------------------- pool growth

#' Index of the mean-sorted gene closest to `target` (mirrors find_closest_helper).
.tox_find_closest <- function(target, means_sorted) {
  n <- length(means_sorted)
  if (n == 0L) return(0L)
  if (target <= means_sorted[1L]) return(1L)
  if (target >= means_sorted[n]) return(n)
  j <- findInterval(target, means_sorted)    # means_sorted[j] <= target < means_sorted[j+1]
  if (abs(means_sorted[j] - target) <= abs(means_sorted[j + 1L] - target)) j else j + 1L
}

#' Adaptive kNN residual pool, in GENE units (mirrors gather_residuals_helper).
#'
#' Phase 1 grows to `k_start` genes unconditionally. Phase 2 adds `k_step` genes
#' per round and commits the round only if the relative change in the pool's MEAN
#' ABSOLUTE RESIDUAL, (S_new - S_old)/S_old, does not exceed `tau`; a rejected
#' round is discarded whole and expansion stops. Growth also stops at `k_max`
#' genes, at `max_pool_size` residuals, or when the gene list is exhausted.
#'
#' NOTE the stopping rule is ONE-SIDED: only an INCREASE beyond tau stops growth.
#' A round that lowers the mean absolute residual is always accepted, so the pool
#' can only ever be stopped by neighbours that are noisier than the current pool.
#' That asymmetry is one of the things the diagnostics below measure.
#'
#' @return list with the pooled residuals, the sorted-slot indices of the genes
#'   contributing to them, and the reason expansion ended.
tox_gather <- function(target_mean, sorted, k_start, k_step, k_max, tau,
                       max_pool_size = 70000L) {
  ms <- sorted$means_sorted
  ng <- length(ms)
  nr <- sorted$n_samples
  pos <- .tox_find_closest(target_mean, ms)
  if (pos == 0L) return(list(pool = numeric(0), genes = integer(0), stop_reason = "empty"))

  genes <- pos
  left  <- pos - 1L
  right <- pos + 1L

  # Next neighbour gene, alternating left/right by |mean - target|; ties go left,
  # exactly as `choose_index` does.
  take_next <- function() {
    if (left >= 1L && right <= ng) {
      if (abs(ms[left] - target_mean) <= abs(ms[right] - target_mean)) {
        idx <- left;  left  <<- left - 1L
      } else {
        idx <- right; right <<- right + 1L
      }
    } else if (left >= 1L) {
      idx <- left;  left  <<- left - 1L
    } else if (right <= ng) {
      idx <- right; right <<- right + 1L
    } else return(NA_integer_)
    idx
  }

  cap_genes <- max_pool_size %/% nr           # residual cap expressed in genes
  stop_reason <- "k_max"

  # ---- Phase 1: unconditional growth to k_start genes
  while (length(genes) < k_start && length(genes) < cap_genes &&
         (left >= 1L || right <= ng)) {
    idx <- take_next(); if (is.na(idx)) break
    genes <- c(genes, idx)
  }
  n_pool <- length(genes) * nr
  if (n_pool < 10L)
    return(list(pool = numeric(0), genes = genes, stop_reason = "under_10"))

  pool_abs_sum <- sum(abs(sorted$resid[, genes, drop = FALSE]))
  S_old <- pool_abs_sum / n_pool
  if (S_old == 0)
    return(list(pool = as.numeric(sorted$resid[, genes, drop = FALSE]),
                genes = genes, stop_reason = "zero_spread"))

  # ---- Phase 2: adaptive rounds of k_step genes
  n_rounds <- 0L
  repeat {
    if (length(genes) >= k_max)     { stop_reason <- "k_max";     break }
    if (length(genes) >= cap_genes) { stop_reason <- "pool_cap";  break }
    if (left < 1L && right > ng)    { stop_reason <- "exhausted"; break }
    cand <- integer(0)
    while (length(cand) < k_step && length(genes) + length(cand) < k_max &&
           length(genes) + length(cand) < cap_genes && (left >= 1L || right <= ng)) {
      idx <- take_next(); if (is.na(idx)) break
      cand <- c(cand, idx)
    }
    if (!length(cand)) { stop_reason <- "exhausted"; break }
    trial_abs_sum <- pool_abs_sum + sum(abs(sorted$resid[, cand, drop = FALSE]))
    trial_size    <- (length(genes) + length(cand)) * nr
    S_new         <- trial_abs_sum / trial_size
    if ((S_new - S_old) / S_old > tau) { stop_reason <- "tau"; break }
    genes <- c(genes, cand); pool_abs_sum <- trial_abs_sum; S_old <- S_new
    n_rounds <- n_rounds + 1L
  }

  list(pool = as.numeric(sorted$resid[, genes, drop = FALSE]),
       genes = genes, stop_reason = stop_reason, n_rounds = n_rounds)
}

#' Symmetric quantile trim of a residual pool (mirrors trim_pool_tails_helper).
tox_trim <- function(pool, trim_frac) {
  if (!length(pool) || trim_frac <= 0) return(pool)
  k <- floor(length(pool) * trim_frac)
  if (k < 1L || (length(pool) - 2L * k) < 1L) return(pool)
  s <- sort(pool)
  s[(k + 1L):(length(pool) - k)]
}

# ------------------------------------------------------------------ p-value

#' Exact add-one-corrected pairwise tail p-value (mirrors compute_pvalue_helper).
#'
#'   p = ( #{(a,b) : |a - b| >= |obs|} + 1 ) / ( n_a * n_b + 1 )
#'
#' Counted with a sorted control pool and two binary searches per case residual,
#' as the Fortran does, so it is O(n_a log n_b) rather than O(n_a n_b).
#' `#{b < a+t}` is strict and `#{b <= a-t}` inclusive, matching the Fortran's
#' `count_below_helper` calls, so the boundary |a-b| == t counts as EXCEEDING.
tox_pvalue_exact <- function(pool_case, pool_ctrl, obs) {
  na <- length(pool_case); nb <- length(pool_ctrl)
  if (!na || !nb) return(NA_real_)
  t <- abs(obs)
  b <- sort(pool_ctrl)
  hi <- findInterval(pool_case + t, b, left.open = TRUE)   # #{b <  a+t}
  lo <- findInterval(pool_case - t, b, left.open = FALSE)  # #{b <= a-t}
  inside <- pmax(0L, hi - lo)                              # #{a-t < b < a+t}
  count_ge <- sum(as.numeric(nb - inside))
  (count_ge + 1) / (as.numeric(na) * as.numeric(nb) + 1)
}

# --------------------------------------------------------------- diagnostics

#' Per-gene run of the exact `own` pipeline WITH full instrumentation.
#'
#' Returns one row per gene carrying, besides the p-value the Fortran would give:
#'
#'   own-scale        sd_case / sd_ctrl        the gene's OWN replicate sd (Bessel)
#'   pool-scale       sd_pool_case/_ctrl       sd of the neighbourhood pool it is
#'                                             scored against
#'   scale ratio      rho_case / rho_ctrl      sd_own / sd_pool. rho = 1 means the
#'                                             pool describes this gene's noise;
#'                                             rho spread across genes is variance
#'                                             heterogeneity WITHIN a mean-neighbourhood
#'   pool SHAPE       kurt_pool_*              excess kurtosis of the pool. A pool
#'                                             mixing genes of different sd is a SCALE
#'                                             MIXTURE: leptokurtic, i.e. a NARROW CORE
#'                                             with heavy tails, even at correct variance
#'   null scale       sd_null                  sd of the null distance distribution,
#'                                             sqrt(sd_pool_case^2 + sd_pool_ctrl^2)/sqrt(n)
#'   standardised     z_pool = |obs| / sd_null       scored against the pooled null
#'   observed         z_own  = |obs| / sd_own_null   scored against the gene's OWN sd
#'                                             (the two together separate "the pool is
#'                                             wrong" from "the observed distance is
#'                                             large for this gene")
#'   neighbourhood    n_genes_pool_*, mean span, asymmetry (share of neighbours BELOW
#'                    the target in mean), stop_reason
#'   resolution       p_floor = 1/(n_a n_b + 1), the smallest p this gene COULD get
#'
#' @param case_mat,ctrl_mat samples x genes, same gene columns, PRE-LOG space
#' @param obs optional observed statistic per gene; default = the difference of
#'   group means on the residual scale the model works in (linear for raw,
#'   log2(x+1) means for log) -- i.e. what null_calibration.R passes.
#' @param genes optional subset (indices or names) to instrument, for speed.
tox_diagnose <- function(case_mat, ctrl_mat, norm_method,
                         k_start = 20L, k_step = 1L, k_max = 50L, tau = 0.1,
                         trim_frac = 0.0, max_pool_size = 70000L,
                         obs = NULL, genes = NULL, verbose = TRUE) {
  stopifnot(ncol(case_mat) == ncol(ctrl_mat))
  ng <- ncol(case_mat)
  gid <- colnames(case_mat); if (is.null(gid)) gid <- paste0("g", seq_len(ng))

  sc <- tox_prepare_sorted(case_mat, norm_method)
  st <- tox_prepare_sorted(ctrl_mat, norm_method)
  mc <- colMeans(case_mat); mt <- colMeans(ctrl_mat)

  if (is.null(obs)) {
    obs <- if (norm_method == 0L) mc - mt
           else colMeans(log2(case_mat + 1)) - colMeans(log2(ctrl_mat + 1))
  }
  # sorted-slot lookup so a gene's own residuals can be found in O(1)
  slot_c <- integer(ng); slot_c[sc$order] <- seq_len(ng)
  slot_t <- integer(ng); slot_t[st$order] <- seq_len(ng)

  idx <- if (is.null(genes)) seq_len(ng)
         else if (is.character(genes)) match(genes, gid) else as.integer(genes)
  idx <- idx[!is.na(idx)]

  na_c <- nrow(case_mat); na_t <- nrow(ctrl_mat)
  sqc <- 1 / sqrt(na_c); sqt <- 1 / sqrt(na_t)
  exkurt <- function(v) { n <- length(v); if (n < 4L) return(NA_real_)
                          m <- mean(v); s <- sd(v); if (!is.finite(s) || s == 0) return(NA_real_)
                          sum(((v - m) / s)^4) / n - 3 }

  out <- vector("list", length(idx))
  for (ii in seq_along(idx)) {
    g <- idx[ii]
    if (verbose && ii %% 2000L == 0L) message("  ... gene ", ii, "/", length(idx))
    gc_ <- tox_gather(mc[g], sc, k_start, k_step, k_max, tau, max_pool_size)
    gt_ <- tox_gather(mt[g], st, k_start, k_step, k_max, tau, max_pool_size)
    pc <- tox_trim(gc_$pool, trim_frac); pt <- tox_trim(gt_$pool, trim_frac)
    if (length(pc) < 10L || length(pt) < 10L) next

    own_c <- sc$resid[, slot_c[g]]      # already Bessel-corrected
    own_t <- st$resid[, slot_t[g]]
    sd_own_c <- sqrt(sum(own_c^2) / na_c)   # = Bessel-corrected sd of replicates
    sd_own_t <- sqrt(sum(own_t^2) / na_t)

    spc <- sd(pc); spt <- sd(pt)
    sd_null     <- sqrt(spc^2 / na_c + spt^2 / na_t)
    sd_own_null <- sqrt(sd_own_c^2 / na_c + sd_own_t^2 / na_t)

    p <- tox_pvalue_exact(pc * sqc, pt * sqt, obs[g])

    # per-neighbour-gene sd, to measure the scale MIXTURE the pool is made of
    sd_nb_c <- sqrt(colMeans(sc$resid[, gc_$genes, drop = FALSE]^2))
    sd_nb_t <- sqrt(colMeans(st$resid[, gt_$genes, drop = FALSE]^2))

    out[[ii]] <- data.frame(
      gene_id = gid[g], mean_case = mc[g], mean_control = mt[g], obs = obs[g],
      n_case = na_c, n_control = na_t, p = p,
      p_floor = 1 / (length(pc) * length(pt) + 1),
      sd_own_case = sd_own_c, sd_own_control = sd_own_t,
      sd_pool_case = spc, sd_pool_control = spt,
      rho_case = sd_own_c / spc, rho_control = sd_own_t / spt,
      kurt_pool_case = exkurt(pc), kurt_pool_control = exkurt(pt),
      # spread of the per-neighbour sd inside one pool: the mixture width that
      # makes the pool leptokurtic. sd of log sd is scale-free and comparable
      # between raw and log normalisation.
      sd_log_nb_sd_case    = if (length(sd_nb_c) > 2) sd(log(pmax(sd_nb_c, 1e-12))) else NA_real_,
      sd_log_nb_sd_control = if (length(sd_nb_t) > 2) sd(log(pmax(sd_nb_t, 1e-12))) else NA_real_,
      sd_null = sd_null, sd_own_null = sd_own_null,
      z_pool = abs(obs[g]) / sd_null, z_own = abs(obs[g]) / sd_own_null,
      n_resid_pool_case = length(pc), n_resid_pool_control = length(pt),
      n_genes_pool_case = length(gc_$genes), n_genes_pool_control = length(gt_$genes),
      mean_span_case    = diff(range(sc$means_sorted[gc_$genes])),
      mean_span_control = diff(range(st$means_sorted[gt_$genes])),
      # 0.5 = neighbours balanced around the target in mean; deviation = the
      # neighbourhood is one-sided, which in a right-skewed mean distribution
      # systematically pairs a gene with LOWER-variance neighbours
      frac_below_case    = mean(sc$means_sorted[gc_$genes] < mc[g]),
      frac_below_control = mean(st$means_sorted[gt_$genes] < mt[g]),
      stop_case = gc_$stop_reason, stop_control = gt_$stop_reason,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# --------------------------------------------------------------- validation

#' Compare this re-implementation's p-values against the Fortran pipeline's.
#'
#' MUST be run before any conclusion is drawn from `tox_diagnose`. Requires
#' `tox_compute_noise_pvalues_pipeline_exact` (rcpp/tensoromics_functions.R) to be
#' loaded. Returns the max/median absolute difference over genes with a p-value
#' from both, and the count of genes where only one implementation produced one.
#'
#' Expect max_abs_diff to be 0 exactly, or at most a few ULP: both compute the
#' same integer pair count, so any real discrepancy is a logic difference, not
#' floating-point noise.
validate_against_fortran <- function(case_mat, ctrl_mat, norm_method,
                                     k_start = 20L, k_step = 1L, k_max = 50L,
                                     tau = 0.1, trim_frac = 0.0,
                                     max_pool_size = 70000L, n_genes_check = 500L) {
  if (!exists("tox_compute_noise_pvalues_pipeline_exact"))
    stop("Load rcpp/tensoromics_functions.R first.")
  obs <- if (norm_method == 0L) colMeans(case_mat) - colMeans(ctrl_mat)
         else colMeans(log2(case_mat + 1)) - colMeans(log2(ctrl_mat + 1))
  obs <- as.numeric(obs); valid <- as.integer(is.finite(obs)); obs[!is.finite(obs)] <- 0
  fres <- tox_compute_noise_pvalues_pipeline_exact(
    case_means = as.numeric(colMeans(case_mat)), case_replicates = case_mat,
    control_means = as.numeric(colMeans(ctrl_mat)), control_replicates = ctrl_mat,
    obs_own = obs, valid_genes_own = valid, norm_method = as.integer(norm_method),
    k_start = k_start, k_step = k_step, k_max = k_max, tau = tau,
    trim_frac = trim_frac, max_pool_size = max_pool_size)
  pf <- fres$pvalues_own; pf[pf < 0 | pf > 1] <- NA

  sel <- sort(sample.int(ncol(case_mat), min(n_genes_check, ncol(case_mat))))
  d <- tox_diagnose(case_mat, ctrl_mat, norm_method, k_start, k_step, k_max, tau,
                    trim_frac, max_pool_size, obs = obs, genes = sel, verbose = FALSE)
  pr <- setNames(d$p, d$gene_id)
  gid <- colnames(case_mat); if (is.null(gid)) gid <- paste0("g", seq_len(ncol(case_mat)))
  common <- intersect(names(pr)[!is.na(pr)], gid[sel][!is.na(pf[sel])])
  diff <- abs(pr[common] - pf[match(common, gid)])
  # Neighbourhood sizes are a second, independent check on the gather logic.
  nb <- fres$neighborhood_size_own_case[match(common, gid)]
  list(n_compared = length(common),
       max_abs_diff = if (length(diff)) max(diff) else NA_real_,
       median_abs_diff = if (length(diff)) median(diff) else NA_real_,
       n_only_r = sum(!is.na(pr) & is.na(pf[match(names(pr), gid)])),
       n_only_fortran = sum(is.na(pr[gid[sel]]) & !is.na(pf[sel])),
       max_abs_diff_nbhd = max(abs(d$n_resid_pool_case[match(common, d$gene_id)] - nb)),
       ok = isTRUE(length(diff) > 0 && max(diff) < 1e-12))
}

# ------------------------------------------------------- neighbour expansion

#' The neighbour genes that made up one gene's null, with their own statistics.
#'
#' `tox_diagnose` summarises a pool; this returns the pool's MEMBERSHIP -- which
#' genes were pooled, at what mean distance, with what sd -- so a significant
#' gene can be read together with the genes its null was actually built from.
#' That matters because the null is other genes' noise: a call is only as
#' trustworthy as the neighbourhood behind it, and a neighbourhood dominated by
#' one much quieter gene, or spread over a wide mean range, is a caveat on the
#' call, not a detail.
#'
#' @param side "case" or "control"
#' @return data.frame, one row per neighbour gene, ordered by |mean - target|.
tox_neighbours <- function(target_gene, case_mat, ctrl_mat, norm_method, side = "case",
                           k_start = 20L, k_step = 1L, k_max = 50L, tau = 0.1,
                           max_pool_size = 70000L) {
  m <- if (side == "case") case_mat else ctrl_mat
  gid <- colnames(m); g <- if (is.character(target_gene)) match(target_gene, gid) else target_gene
  if (is.na(g)) return(NULL)
  s  <- tox_prepare_sorted(m, norm_method)
  tm <- colMeans(m)[g]
  gg <- tox_gather(tm, s, k_start, k_step, k_max, tau, max_pool_size)
  if (!length(gg$genes)) return(NULL)
  sd_nb <- sqrt(colMeans(s$resid[, gg$genes, drop = FALSE]^2))   # Bessel-corrected
  o <- order(abs(s$means_sorted[gg$genes] - tm))
  data.frame(target_gene = gid[g], side = side,
             neighbour_gene = s$gene_ids[gg$genes][o],
             neighbour_mean = s$means_sorted[gg$genes][o],
             mean_distance  = abs(s$means_sorted[gg$genes] - tm)[o],
             neighbour_sd_resid = unname(sd_nb[o]),
             is_target = s$gene_ids[gg$genes][o] == gid[g],
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------- analytic (derived) tails

#' Saddlepoint tail probability for the mean-difference null -- an ANALYTIC
#' alternative to counting pairs or bootstrap draws.
#'
#' WHY
#' ---
#' The exact pairwise p-value cannot go below `1/(n_a*n_b + 1)`, and a bootstrap
#' cannot go below `1/(B + 1)`. Both floors are properties of the COUNTING, not of
#' the data. edgeR, limma and DESeq2 have no such floor because they never count:
#' they borrow variance information across genes to estimate a noise parameter,
#' then read the tail off a THEORETICAL reference distribution (t, F, normal)
#' analytically, so a p-value of 1e-30 is arithmetic rather than resolution.
#'
#' TOX already does the borrowing -- the kNN mean-neighbourhood does exactly the
#' job limma's expression-dependent prior and edgeR's dispersion trend do. What it
#' lacks is the analytic step. This function supplies it WITHOUT importing a
#' distributional assumption: the reference distribution is the pool's own, entered
#' through its empirical cumulant generating function.
#'
#' WHAT IT COMPUTES
#' Let D = mean(n_case draws from pool_case) - mean(n_control draws from
#' pool_control), i.e. the mean-bootstrap null, evaluated exactly instead of
#' resampled. Its CGF under iid resampling is
#'
#'   K(s) = n_a * log E[exp(s/n_a * A)] + n_b * log E[exp(-s/n_b * B)]
#'
#' with the expectations taken over the empirical pools. Lugannani-Rice then gives
#' P(D >= t). Two-sided by summing both tails. Unlike an Edgeworth expansion the
#' RELATIVE error stays bounded as t grows, which is the whole point here.
#'
#' TWO PROBLEMS, ONE CHANGE. Because D is a difference of MEANS, this null carries
#' the CLT shape the observed statistic has -- so it also removes the shape
#' mismatch that makes the raw arm anti-conservative at alpha = 0.05 (see
#' docs/raw_normalisation_diagnosis.md). It is deterministic (no RNG, no seed) and
#' costs one 1-D root find per gene.
#'
#' VALIDATED (scripts/pvalue_tail_methods.R):
#'   * against an analytically known null (large N(0,1) pools, so D ~ N(0, 2/n)):
#'     ratio to truth 1.15 / 1.30 / 1.41 / 1.41 / 1.21 / 0.81 at true p = 1e-2
#'     down to 1e-12. Bounded relative error, no floor.
#'   * against 2e7 bootstrap draws from a heterogeneous n_rep = 3 pool: within a
#'     factor 1.6 at true p = 1e-5, where the pairwise count has already hit its
#'     floor and a 1e4 bootstrap has hit its own.
#'
#' THE LIMIT IT DOES NOT REMOVE -- read before using it
#' This computes the tail of the pool you give it to arbitrary precision. It does
#' NOT reduce uncertainty in the POOL. At n_rep = 3, k = 50 the effective sample
#' size behind a pool is ~200, and redrawing the neighbourhood moves the answer by
#' 21x at true p = 1e-2, ~500x at 1e-3 and ~2e4x at 1e-4 (only 18% of
#' neighbourhoods land within 2x of truth there). So at low replicate counts an
#' analytic tail replaces an honest, visible floor with a number that LOOKS precise
#' and is not. Use it where the pool is well estimated -- the real case-vs-control
#' analysis, with tens to hundreds of samples -- and at n_rep = 3 prefer reporting
#' a bound ("p < 5e-3") over a small number nothing supports.
#'
#' @param pool_case,pool_control residual pools, UNSCALED (do not pre-divide by
#'   sqrt(n)): the averaging is in the CGF, so passing sqrt-scaled pools would
#'   apply the correction twice.
#' @param obs observed statistic (sign ignored)
#' @param n_case,n_control replicates per side
#' @param two_sided sum both tails (default TRUE, matching the pairwise p-value)
#' @return tail probability, or NA if the saddlepoint equation has no root (obs
#'   beyond the pool's attainable mean range) -- callers should fall back to the
#'   pairwise floor in that case.
tox_pvalue_saddlepoint <- function(pool_case, pool_control, obs,
                                   n_case, n_control, two_sided = TRUE) {
  # log E[exp(u * x)], computed by max-subtraction so large u cannot overflow
  lmgf  <- function(u, x) { z <- u * x; m <- max(z); m + log(mean(exp(z - m))) }
  tilt  <- function(u, x) { z <- u * x; z <- z - max(z); w <- exp(z); w / sum(w) }
  Kf  <- function(s, a, b, na, nb) na * lmgf(s/na, a) + nb * lmgf(-s/nb, b)
  K1f <- function(s, a, b, na, nb) sum(tilt(s/na, a) * a) - sum(tilt(-s/nb, b) * b)
  K2f <- function(s, a, b, na, nb) {
    wa <- tilt(s/na, a); wb <- tilt(-s/nb, b)
    va <- sum(wa * a^2) - sum(wa * a)^2; vb <- sum(wb * b^2) - sum(wb * b)^2
    va/na + vb/nb
  }
  one_tail <- function(t, a, b, na, nb) {
    if (!is.finite(t)) return(NA_real_)
    if (t <= K1f(0, a, b, na, nb)) return(NA_real_)   # not in the upper tail
    # K1 is strictly increasing in s; bracket the root, giving up if the target
    # exceeds what the pools can produce (t >= max(a) - min(b)).
    hi <- 1e-8
    for (i in 1:200) { if (K1f(hi, a, b, na, nb) > t) break; hi <- hi * 2 }
    if (K1f(hi, a, b, na, nb) <= t) return(NA_real_)
    s <- tryCatch(uniroot(function(s) K1f(s, a, b, na, nb) - t, c(0, hi),
                          tol = .Machine$double.eps^0.6)$root,
                  error = function(e) NA_real_)
    if (!is.finite(s) || s <= 0) return(NA_real_)
    k <- Kf(s, a, b, na, nb); k2 <- K2f(s, a, b, na, nb)
    if (!is.finite(k2) || k2 <= 0) return(NA_real_)
    w <- sqrt(max(2 * (s * t - k), 0)); u <- s * sqrt(k2)
    if (w < 1e-8 || u < 1e-12) return(0.5)            # t at the mean: LR degenerates
    p <- pnorm(w, lower.tail = FALSE) + dnorm(w) * (1/w - 1/u)
    min(max(p, 0), 1)
  }
  t <- abs(obs)
  # An observed distance at or below the null's centre is not in either tail; the
  # Lugannani-Rice formula degenerates there and the honest answer is p = 1.
  mu0 <- mean(pool_case) - mean(pool_control)
  if (!is.finite(t) || t <= abs(mu0)) return(1)
  up <- one_tail(t, pool_case, pool_control, n_case, n_control)
  if (!two_sided) return(up)
  # lower tail of D = upper tail of -D, and -D is the same construction with the
  # pools negated and swapped
  lo <- one_tail(t, -pool_control, -pool_case, n_control, n_case)
  s  <- sum(c(up, lo), na.rm = TRUE)
  if (all(is.na(c(up, lo)))) return(NA_real_)
  min(max(s, 0), 1)
}

# ------------------------------------------ gene-blocked mean null (enumerated)

#' All distinct bootstrap-mean values of one gene's residuals, with multiplicities.
#'
#' Resampling `n` of a gene's `n` residuals with replacement gives `n^n` ordered
#' tuples but only `C(2n-1, n)` distinct MULTISETS, each with a multinomial weight.
#' Enumerating the multisets is exact and far cheaper: 10 instead of 27 at n = 3,
#' 462 instead of 46,656 at n = 6, 2.8e6 instead of 3e11 at n = 10.
.tox_gene_mean_multiset <- function(r) {
  n <- length(r)
  comps <- function(rem, slot) {                    # compositions of n into n parts
    if (slot == n) return(matrix(rem, 1, 1))
    do.call(rbind, lapply(0:rem, function(c) cbind(c, comps(rem - c, slot + 1))))
  }
  cnt <- comps(n, 1L)
  list(v = as.numeric(cnt %*% r) / n,
       w = apply(cnt, 1, function(c) exp(lfactorial(n) - sum(lfactorial(c)))))
}

#' Exact p-value from the GENE-BLOCKED mean-difference null.
#'
#' WHAT IT FIXES
#' The bootstrap null (`tox_noise_model.F90`) draws `n_rep` residuals iid from the
#' POOLED neighbourhood, so one null "mean" can combine a residual from a quiet
#' gene with one from a noisy gene. A real gene's mean is never formed that way:
#' its `n_rep` values all come from ONE noise level. Blocking the draw -- pick a
#' neighbour gene, then resample within that gene -- makes each null draw as
#' coherent as the observed statistic it is scored against.
#'
#' WHY IT CAN BE ENUMERATED
#' Because the within-gene resample has only `C(2n-1, n)` distinct outcomes, the
#' ENTIRE blocked null can be written down: `k * C(2n-1, n)` weighted values per
#' side. There is no draw count and no RNG. The p-value is then the same weighted
#' pairwise tail count the exact model already performs, so this reuses the sorted
#' pool + binary-search machinery rather than adding new numerics.
#'
#' RESOLUTION
#' The floor is `1 / (W_case * W_control + 1)` with `W = k * n^n` (the total
#' weight, not the multiset count). At k = 30 neighbour genes and n_rep = 3 that is
#' 810 x 810 = 656,100 -> floor 1.5e-6, BELOW the 4.2e-6 that BH at q = 0.05 needs
#' from 12,000 genes. The individual-residual pairwise null at the same k and n_rep
#' gives 90 x 90 -> 1.2e-4, and a 25,000-draw bootstrap gives 4.0e-5; matching
#' 1.5e-6 by sampling would need ~250,000 draws per gene.
#'
#' MEASURED BEHAVIOUR (H0 split-half simulation, raw normalisation, k_start = 30,
#' k_max = 50, n_rep = 3, 5,000 gene-observations per cell, +- 2 MC SE):
#'
#'   cv_het   exact pairs (current)          enumerated blocked
#'   0.0      infl05 0.90+-0.12 / infl01 0.54    0.62+-0.10 / 0.08
#'   0.4      infl05 0.88+-0.12 / infl01 0.64    0.70+-0.10 / 0.26
#'   0.7      infl05 1.16+-0.13 / infl01 0.72    0.91+-0.12 / 0.28
#'
#' So at n_rep = 3 it is CONSERVATIVE -- the safe direction, but it costs power. At
#' n_rep = 10 the same blocked null (sampled) was the best-calibrated of four
#' constructions tested, at both levels (infl05 1.16, infl01 1.07). Verified to
#' reproduce the blocked bootstrap as B -> infinity: 0.1002 vs 0.1001, 0.01130 vs
#' 0.01120, 9.25e-4 vs 9.18e-4 at B = 2e6.
#'
#' COST
#' `C(2n-1, n) * k` values per side: 300 at n_rep = 3 / k = 30, 3,780 at n_rep = 5,
#' 2.8e6 (22 MB) at n_rep = 10. Enumerate up to about n_rep = 10 and fall back to
#' sampling above it -- by which point replicate counts are no longer the binding
#' problem.
#'
#' @param resid_case,resid_control n_rep x k_genes residual matrices for the
#'   neighbourhood (columns = neighbour genes), as produced by `tox_gather`'s
#'   `genes` index into `sorted$resid`. NOT sqrt-scaled: the averaging is explicit.
#' @param obs observed statistic (sign ignored)
tox_pvalue_blocked_exact <- function(resid_case, resid_control, obs) {
  if (!is.matrix(resid_case))    resid_case    <- as.matrix(resid_case)
  if (!is.matrix(resid_control)) resid_control <- as.matrix(resid_control)
  build <- function(m) {
    parts <- lapply(seq_len(ncol(m)), function(j) .tox_gene_mean_multiset(m[, j]))
    list(v = unlist(lapply(parts, `[[`, "v")), w = unlist(lapply(parts, `[[`, "w")))
  }
  A <- build(resid_case); B <- build(resid_control)
  t <- abs(obs)
  o <- order(B$v); bv <- B$v[o]; bw <- B$w[o]
  cw <- c(0, cumsum(bw))                       # cw[i+1] = weight of the i smallest
  Wb <- cw[length(cw)]; Wa <- sum(A$w)
  hi <- findInterval(A$v + t, bv, left.open = TRUE)   # weight of b <  a+t
  lo <- findInterval(A$v - t, bv, left.open = FALSE)  # weight of b <= a-t
  inside <- cw[hi + 1L] - cw[lo + 1L]
  count_ge <- sum(A$w * (Wb - pmax(inside, 0)))
  (count_ge + 1) / (Wa * Wb + 1)
}
