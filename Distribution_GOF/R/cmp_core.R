# cmp_core.R -- mean-parameterised Conway-Maxwell-Poisson (brief 5.6).
#
#   p(y) = lambda^y / (y!)^nu / Z(lambda, nu),   lambda solved so that E[Y] = mu.
#
# This is the glmmTMB::compois parameterisation (Huang 2017), NOT the
# COMPoissonReg / Sellers-Shmueli lambda-parameterisation (there an offset enters
# log lambda and E[Y] != lambda unless nu = 1). We work with l = log(lambda).
#
# Normalising sums use log-sum-exp over a WINDOW [max(0, mu - K s), mu + K s],
# s = sqrt(mu / nu) (the asymptotic sd), so the cost is O(s) rather than O(mu).
# The window depends only on (mu, nu), not on l, so every Newton step sums the
# same terms. After convergence the boundary terms are required to be < 1e-15 of
# the maximum term; otherwise K is doubled and the solve repeated.

# lgamma(y + 1) table, grown on demand (each fork keeps its own copy).
.cmp_lg <- new.env()
.cmp_lg$tab <- lgamma(seq_len(1e5))            # tab[y + 1] = lgamma(y + 1), y = 0..99999
cmp_lgfact <- function(y) {
  top <- max(y) + 1
  if (top > length(.cmp_lg$tab)) .cmp_lg$tab <- lgamma(seq_len(max(top, 2 * length(.cmp_lg$tab))))
  .cmp_lg$tab[y + 1]
}

cmp_window <- function(mu, nu, K) {
  s <- sqrt(max(mu / nu, 1e-8))
  c(lo = max(0, floor(mu - K * s)), hi = ceiling(mu + K * s) + 10)
}

# Moments of the windowed CMP at log-lambda l. Returns logZ, E, V, and the
# boundary ratios needed for the window assertion.
# log-terms up to a constant: y l - nu lgamma(y+1) = nu * log dpois(y, e^(l/nu)) + nu e^(l/nu).
# The dpois form (Loader's algorithm) has no catastrophic cancellation; the plain
# form subtracts numbers of size ~ y log y and, via a shared logZ, biased every
# probability by the same ~1e-12 factor at mu ~ 5e4 (the cdf then hit 1 early).
# The constant cancels in every normalised quantity. The plain form is used only
# if e^(l/nu) is not finite (transient Newton iterates).
.cmp_terms <- function(l, nu, y, lg) {
  lam <- exp(l / nu)
  if (is.finite(lam) && lam > 0) nu * stats::dpois(y, lam, log = TRUE) else y * l - nu * lg
}

.cmp_moments <- function(l, nu, y, lg, full = FALSE) {
  t <- .cmp_terms(l, nu, y, lg)
  mx <- max(t)
  w <- exp(t - mx); S <- sum(w)
  E <- sum(y * w) / S
  V <- sum((y - E)^2 * w) / S
  out <- c(logZ = mx + log(S), E = E, V = V,
           edge_lo = if (y[1] > 0) w[1] else 0, edge_hi = w[length(w)])
  if (full) {                                     # for the analytic score (cmp_fit)
    Elg <- sum(lg * w) / S
    out <- c(out, Elg = Elg, C = sum((y - E) * (lg - Elg) * w) / S)
  }
  out
}

#' Solve l = log(lambda) with E[Y] = mu for one observation.
#' Newton on l (dE/dl = Var(Y)), started from the asymptotic approximation
#' log lambda ~ nu * log(mu + (nu - 1) / (2 nu)), or from `l0` (warm start).
#' Falls back to uniroot on a bracket if Newton does not converge.
cmp_solve <- function(mu, nu, K = 40, l0 = NA_real_, max_window = 5e6) {
  if (!(is.finite(mu) && mu > 0 && is.finite(nu) && nu > 0)) stop("cmp_solve: need mu > 0 and nu > 0")
  repeat {
    wn <- cmp_window(mu, nu, K)
    if (wn[["hi"]] - wn[["lo"]] > max_window)
      stop(sprintf("cmp_solve: window %.3g wider than CMP_MAX_WINDOW (mu=%.4g, nu=%.4g)",
                   wn[["hi"]] - wn[["lo"]], mu, nu))
    y <- wn[["lo"]]:wn[["hi"]]
    lg <- cmp_lgfact(y)
    a <- mu + (nu - 1) / (2 * nu)
    l <- if (is.finite(l0)) l0 else if (a > 0) nu * log(a) else nu * log(mu)
    ok <- FALSE
    for (it in seq_len(60L)) {
      m <- .cmp_moments(l, nu, y, lg, full = TRUE)
      d <- m[["E"]] - mu
      if (abs(d) <= 1e-12 * mu) { ok <- TRUE; break }   # the term arithmetic carries ~1e-13 noise at mu ~ 1e4
      if (!(is.finite(m[["V"]]) && m[["V"]] > 0)) break
      step <- d / m[["V"]]
      l <- l - max(-2, min(2, step))
    }
    if (!ok) {
      f <- function(l) .cmp_moments(l, nu, y, lg)[["E"]] - mu
      lo <- l - 1; hi <- l + 1
      for (k in seq_len(60L)) { if (f(lo) < 0) break; lo <- lo - 2^k }
      for (k in seq_len(60L)) { if (f(hi) > 0) break; hi <- hi + 2^k }
      l <- stats::uniroot(f, c(lo, hi), tol = 1e-14, maxiter = 1000L)$root
      m <- .cmp_moments(l, nu, y, lg, full = TRUE)
      if (abs(m[["E"]] - mu) > 1e-10 * mu) stop(sprintf("cmp_solve: no root (mu=%.4g, nu=%.4g)", mu, nu))
    }
    if (m[["edge_lo"]] < 1e-15 && m[["edge_hi"]] < 1e-15)
      return(c(l = l, logZ = m[["logZ"]], lo = wn[["lo"]], hi = wn[["hi"]], K = K,
               V = m[["V"]], Elg = m[["Elg"]], C = m[["C"]]))
    K <- 2 * K                                   # boundary not negligible: widen
  }
}

#' Solve for a vector of observations sharing one nu. `l0` (optional) warm-starts.
cmp_solve_vec <- function(mu, nu, K = 40, l0 = NULL, max_window = 5e6) {
  out <- matrix(NA_real_, length(mu), 8,
                dimnames = list(NULL, c("l", "logZ", "lo", "hi", "K", "V", "Elg", "C")))
  for (i in seq_along(mu))
    out[i, ] <- cmp_solve(mu[i], nu, K, if (is.null(l0)) NA_real_ else l0[i], max_window)
  out
}

#' The per-observation solution is cached on `par` (fit stores it at the optimum),
#' so repeated logpmf / cdf / rgen calls on the same parameters do not re-solve.
cmp_prepare <- function(par, K = 40, max_window = 5e6) {
  if (is.null(par$cmp_sol) || nrow(par$cmp_sol) != length(par$mu))
    par$cmp_sol <- cmp_solve_vec(par$mu, par$nu, K, NULL, max_window)
  par
}

cmp_logpmf <- function(y, par, K = 40, max_window = 5e6) {
  par <- cmp_prepare(par, K, max_window)
  s <- par$cmp_sol
  yy <- pmax(y, 0)
  out <- vapply(seq_along(y), function(i) .cmp_terms(s[i, "l"], par$nu, yy[i], cmp_lgfact(yy[i])), 0) - s[, "logZ"]
  ifelse(y < 0 | y != floor(y), -Inf, out)
}

#' F(q) per observation. Below the window the (tiny) partial sum is computed
#' exactly; inside it, a cumulative sum; above it, 1 (the mass beyond the window
#' is below 1e-15 of the maximum term by construction).
cmp_cdf <- function(q, par, K = 40, max_window = 5e6) {
  par <- cmp_prepare(par, K, max_window)
  s <- par$cmp_sol
  q <- floor(q)
  out <- numeric(length(q))
  for (i in seq_along(q)) {
    if (q[i] < 0) { out[i] <- 0; next }
    lo <- if (q[i] < s[i, "lo"]) 0 else s[i, "lo"]
    hi <- min(q[i], s[i, "hi"])
    y <- lo:hi
    out[i] <- sum(exp(.cmp_terms(s[i, "l"], par$nu, y, cmp_lgfact(y)) - s[i, "logZ"]))
  }
  pmin(out, 1)
}

#' Inversion over the window. A uniform draw above the window's total mass (which
#' is 1 - O(1e-15)) maps to the window's top.
cmp_rgen <- function(par, K = 40, max_window = 5e6) {
  par <- cmp_prepare(par, K, max_window)
  s <- par$cmp_sol
  u <- stats::runif(length(par$mu))
  out <- integer(length(u))
  for (i in seq_along(u)) {
    y <- s[i, "lo"]:s[i, "hi"]
    cs <- cumsum(exp(.cmp_terms(s[i, "l"], par$nu, y, cmp_lgfact(y)) - s[i, "logZ"]))
    out[i] <- as.integer(y[min(length(y), findInterval(u[i], cs) + 1L)])
  }
  out
}

#' Own ML fit: nlminb on (beta, log nu) with the analytic score, started from the NB fit with nu0 chosen
#' so the CMP variance ~ mu/nu matches the NB variance mu + mu^2/size at the mean
#' fitted mu. Newton solutions are warm-started across likelihood evaluations.
#' log nu is box-constrained to [log 1e-6, log 100] purely for numerical safety;
#' a fit that ends on that box is flagged in conv_msg.
cmp_fit <- function(y, X, off, beta0, size0, K = 40, max_window = 5e6) {
  p <- ncol(X)
  mu0 <- mean(exp(X %*% beta0 + off))
  nu0 <- if (is.finite(size0) && size0 > 0) 1 / (1 + mu0 / size0) else 1
  nu0 <- min(max(nu0, 1e-5), 50)
  lgy <- cmp_lgfact(y)
  warm <- new.env(); warm$l <- NULL
  # Objective and analytic score from one solve + one moment pass per observation:
  #   d log p / d log mu = (y - mu) mu / V            (dl/dmu = 1/V from E[Y] = mu)
  #   d log p / d nu     = (y - mu) C / V - lg(y) + E[lg(Y)],  C = Cov(Y, lg(Y)), lg = lgamma(. + 1)
  cache <- new.env()
  eval_at <- function(th) {
    if (!identical(cache$th, th)) {
      cache$th <- th; cache$f <- Inf; cache$g <- rep(NA_real_, p + 1L)
      beta <- th[seq_len(p)]; nu <- exp(th[p + 1L])
      mu <- as.vector(exp(X %*% beta + off))
      if (all(is.finite(mu)) && all(mu > 0)) {
        sol <- tryCatch(cmp_solve_vec(mu, nu, K, warm$l, max_window), error = function(e) NULL)
        if (!is.null(sol)) {
          warm$l <- sol[, "l"]
          mom <- sol                                   # V, Elg, C come from the solver's last pass
          cache$f <- -sum(vapply(seq_along(y), function(i) .cmp_terms(sol[i, "l"], nu, y[i], lgy[i]), 0) - sol[, "logZ"])
          dlm <- (y - mu) * mu / mom[, "V"]
          dnu <- (y - mu) * mom[, "C"] / mom[, "V"] - lgy + mom[, "Elg"]
          cache$g <- -c(crossprod(X, dlm), nu * sum(dnu))
        }
      }
    }
    cache
  }
  negll <- function(th) eval_at(th)$f
  grad  <- function(th) eval_at(th)$g
  lb <- c(rep(-Inf, p), log(1e-6)); ub <- c(rep(Inf, p), log(100))
  opt <- stats::nlminb(c(beta0, log(nu0)), negll, grad, lower = lb, upper = ub,
                       control = list(eval.max = 600L, iter.max = 400L))
  beta <- opt$par[seq_len(p)]; nu <- exp(opt$par[p + 1L])
  mu <- as.vector(exp(X %*% beta + off))
  sol <- cmp_solve_vec(mu, nu, K, warm$l, max_window)
  at_bound <- opt$par[p + 1L] <= lb[p + 1L] + 1e-8 || opt$par[p + 1L] >= ub[p + 1L] - 1e-8
  ll <- sum(vapply(seq_along(y), function(i) .cmp_terms(sol[i, "l"], nu, y[i], lgy[i]), 0) - sol[, "logZ"])
  list(beta = beta, nu = nu, mu = mu, sol = sol, loglik = ll,
       convergence = opt$convergence,
       msg = paste0(opt$message, if (at_bound) " | log(nu) on its numerical box" else ""))
}
