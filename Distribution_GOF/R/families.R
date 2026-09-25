# families.R -- six count families behind one interface (brief section 5).
#
#   fam$fit(y, X, off)  -> list(beta, extra = list(...), mu, loglik, conv_ok, conv_msg, engine)
#   fam$logpmf(y, par)  -> log p(y)          par = c(list(mu = <per obs>), extra)
#   fam$cdf(q, par)     -> F(q), 0 for q < 0, exact at integers, clamped to [0, 1]
#   fam$rgen(par)       -> one integer draw per observation
#   fam$mu_fn(beta, extra, X, off) -> mu for new observations (held-out score)
#   fam$moments(par)    -> list(mean, var) per observation (validation only)
#
# `mu_fn` and `moments` are the two additions to the brief's interface: the
# held-out score needs means for samples the fit did not see, and G2 needs the true
# mean and variance. `mu` is always the family's mean parameter: the conditional
# NB mean for ZINB, and E[Y] for PLN.
#
# glmmTMB is called with se = FALSE (no sdreport): pdHess is never a failure
# criterion (brief section 6), and skipping the Hessian roughly halves the fit time.

.tmb_fit <- function(y, X, off, family, zi = FALSE) {
  p <- ncol(X)
  d <- data.frame(y = y, off = off, X)
  names(d)[-(1:2)] <- xs <- paste0("x", seq_len(p))
  f <- stats::as.formula(paste("y ~ 0 +", paste(xs, collapse = " + "), "+ offset(off)"))
  msgs <- character(0)
  fit <- withCallingHandlers(
    glmmTMB::glmmTMB(f, data = d, family = family, ziformula = if (zi) ~1 else ~0,
                     control = glmmTMB::glmmTMBControl(parallel = 1L), se = FALSE),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(fit = fit, beta = unname(glmmTMB::fixef(fit)$cond),
       conv = fit$fit$convergence, msg = paste(c(fit$fit$message, msgs), collapse = " | "))
}

.loglink <- function(beta, X, off) as.vector(exp(X %*% beta + off))

.result <- function(beta, extra, mu, loglik, conv, msg, engine)
  list(beta = beta, extra = extra, mu = mu, loglik = loglik,
       conv_ok = identical(as.integer(conv), 0L), conv_msg = msg, engine = engine)

# ---------------------------------------------------------------- Poisson
fam_poisson <- function(cfg) list(
  name = "poisson", n_extra = 0L,
  fit = function(y, X, off) {
    f <- stats::glm.fit(X, y, offset = off, family = stats::poisson())
    mu <- as.vector(f$fitted.values)
    .result(unname(f$coefficients), list(), mu, sum(stats::dpois(y, mu, log = TRUE)),
            if (isTRUE(f$converged)) 0L else 1L, if (isTRUE(f$converged)) "" else "glm.fit did not converge",
            "glm.fit")
  },
  logpmf  = function(y, par) stats::dpois(y, par$mu, log = TRUE),
  cdf     = function(q, par) stats::ppois(q, par$mu),
  rgen    = function(par) stats::rpois(length(par$mu), par$mu),
  mu_fn   = function(beta, extra, X, off) .loglink(beta, X, off),
  moments = function(par) list(mean = par$mu, var = par$mu)
)

# ---------------------------------------------------------------- NB
# size = sigma(fit); V = mu + mu^2/size. The issue's phi is 1/size; both are reported.
fam_nb <- function(cfg) list(
  name = "nb", n_extra = 1L,
  fit = function(y, X, off) {
    r <- .tmb_fit(y, X, off, glmmTMB::nbinom2())
    size <- stats::sigma(r$fit)
    mu <- .loglink(r$beta, X, off)
    .result(r$beta, list(size = size, nb_phi_issue = 1 / size), mu, as.numeric(stats::logLik(r$fit)),
            r$conv, r$msg, "glmmTMB::nbinom2")
  },
  logpmf  = function(y, par) stats::dnbinom(y, mu = par$mu, size = par$size, log = TRUE),
  cdf     = function(q, par) stats::pnbinom(q, mu = par$mu, size = par$size),
  rgen    = function(par) stats::rnbinom(length(par$mu), mu = par$mu, size = par$size),
  mu_fn   = function(beta, extra, X, off) .loglink(beta, X, off),
  moments = function(par) list(mean = par$mu, var = par$mu + par$mu^2 / par$size)
)

# ---------------------------------------------------------------- ZINB
# pi is constant per gene and does NOT scale with library size -- a modelling choice.
fam_zinb <- function(cfg) list(
  name = "zinb", n_extra = 2L,
  fit = function(y, X, off) {
    r <- .tmb_fit(y, X, off, glmmTMB::nbinom2(), zi = TRUE)
    size <- stats::sigma(r$fit)
    pi0 <- unname(stats::plogis(glmmTMB::fixef(r$fit)$zi))
    mu <- .loglink(r$beta, X, off)
    .result(r$beta, list(size = size, pi = pi0), mu, as.numeric(stats::logLik(r$fit)),
            r$conv, r$msg, "glmmTMB::nbinom2+zi")
  },
  logpmf = function(y, par) {
    nb <- stats::dnbinom(y, mu = par$mu, size = par$size, log = TRUE)
    ifelse(y == 0, log(par$pi + (1 - par$pi) * exp(nb)), log1p(-par$pi) + nb)
  },
  cdf  = function(q, par) ifelse(q < 0, 0, pmin(1, par$pi + (1 - par$pi) * stats::pnbinom(q, mu = par$mu, size = par$size))),
  rgen = function(par) {
    n <- length(par$mu)
    ifelse(stats::runif(n) < par$pi, 0L, stats::rnbinom(n, mu = par$mu, size = par$size))
  },
  mu_fn   = function(beta, extra, X, off) .loglink(beta, X, off),
  moments = function(par) {
    m <- par$mu; v <- m + m^2 / par$size
    list(mean = (1 - par$pi) * m, var = (1 - par$pi) * (v + m^2) - ((1 - par$pi) * m)^2)
  }
)

# ---------------------------------------------------------------- generalised Poisson
# Joe-Zhu parameterisation with the mapping verified in brief section 1:
#   phi = predict(fit, type = "disp"), alpha = 1 - 1/sqrt(phi),
#   lambda1 = mu (1 - alpha), lambda2 = alpha,
#   log p(y) = log lambda1 + (y-1) log(lambda1 + y lambda2) - (lambda1 + y lambda2) - lgamma(y+1).
# Var = mu * phi. Under-dispersion (phi < 1 => lambda2 < 0): terms with
# lambda1 + y lambda2 <= 0 are ZERO, so the pmf need not sum to 1. The mass S is
# computed and stored, never renormalised; S < 1 - 1e-8 flags the gene.
# Written as p(y) = (lambda1 / t) * dpois(y, t), t = lambda1 + y lambda2 -- the same
# function, but R's dpois (Loader's algorithm) avoids the cancellation of the
# textbook form, whose terms reach ~y log y and made the pmf sum to 1 + 8e-12.
gp_logpmf <- function(y, mu, phi) {
  a <- 1 - 1 / sqrt(phi); l1 <- mu * (1 - a); t <- l1 + y * a
  out <- log(l1) - log(pmax(t, 1e-300)) + stats::dpois(pmax(y, 0), pmax(t, 1e-300), log = TRUE)
  out[t <= 0 | y < 0] <- -Inf
  out
}
# Support window: mean +- K sd, then WIDENED until the edge terms are < 1e-17 of
# the term at the mean (strongly over-dispersed GP has a heavy right tail -- at
# phi = 20 the plain window ended where p was still 6.5e-7). The upper end also
# stops where lambda1 + y lambda2 hits 0 (under-dispersion).
.gp_window <- function(mu, phi, K = 40) {
  sd <- sqrt(mu * phi)
  a <- 1 - 1 / sqrt(phi)
  cap <- if (a < 0) ceiling(-mu * (1 - a) / a) else Inf
  lo <- max(0, floor(mu - K * sd)); hi <- min(cap, ceiling(mu + K * sd) + 10)
  ref <- max(gp_logpmf(unique(c(floor(mu), ceiling(mu))), mu, phi))
  while (hi < cap && gp_logpmf(hi, mu, phi) > ref + log(1e-17)) hi <- min(cap, 2 * hi + 10)
  while (lo > 0 && gp_logpmf(lo, mu, phi) > ref + log(1e-17)) lo <- max(0, lo - ceiling(K * sd))
  c(lo = lo, hi = hi)
}
gp_mass <- function(mu, phi) vapply(seq_along(mu), function(i) {
  w <- .gp_window(mu[i], phi); y <- 0:w[["hi"]]
  sum(exp(gp_logpmf(y, mu[i], phi)))
}, 0)
gp_cdf <- function(q, mu, phi) {
  q <- floor(q)
  vapply(seq_along(q), function(i) {
    if (q[i] < 0) return(0)
    w <- .gp_window(mu[i], phi)
    lo <- if (q[i] < w[["lo"]]) 0 else w[["lo"]]
    y <- lo:min(q[i], w[["hi"]])
    min(1, sum(exp(gp_logpmf(y, mu[i], phi))))
  }, 0)
}
gp_rgen <- function(mu, phi) vapply(seq_along(mu), function(i) {
  w <- .gp_window(mu[i], phi); y <- w[["lo"]]:w[["hi"]]
  cs <- cumsum(exp(gp_logpmf(y, mu[i], phi)))
  as.integer(y[min(length(y), findInterval(stats::runif(1), cs) + 1L)])
}, 0L)

fam_genpois <- function(cfg) list(
  name = "genpois", n_extra = 1L,
  fit = function(y, X, off) {
    r <- .tmb_fit(y, X, off, glmmTMB::genpois())
    phi <- unname(stats::predict(r$fit, type = "disp")[1])
    mu <- .loglink(r$beta, X, off)
    res <- .result(r$beta, list(phi = phi), mu, as.numeric(stats::logLik(r$fit)),
                   r$conv, r$msg, "glmmTMB::genpois")
    mass <- gp_mass(mu, phi)
    res$gp_mass_min <- min(mass)
    # Under-dispersed GP (Consul's truncation) can also put MORE than 1 in total
    # mass (e.g. mu = 0.5, phi = 0.5: 1.021). That is not a distribution: no valid
    # PIT exists, so the fit is not OK (reported with the reason). Mass < 1 is
    # kept and flagged, as the brief specifies.
    if (max(mass) > 1 + 1e-8) {
      res$conv_ok <- FALSE
      res$conv_msg <- paste(res$conv_msg, sprintf("| GP total mass %.6g > 1: not a distribution", max(mass)))
    }
    res
  },
  logpmf  = function(y, par) gp_logpmf(y, par$mu, par$phi),
  cdf     = function(q, par) gp_cdf(q, par$mu, par$phi),
  rgen    = function(par) gp_rgen(par$mu, par$phi),
  mu_fn   = function(beta, extra, X, off) .loglink(beta, X, off),
  moments = function(par) list(mean = par$mu, var = par$mu * par$phi)
)

# ---------------------------------------------------------------- PLN (own fitter)
fam_pln <- function(cfg) {
  nn <- cfg$PLN_NODES; P <- cfg$PLN_CDF_PANELS; gl <- cfg$PLN_CDF_GL
  list(
    name = "pln", n_extra = 1L,
    fit = function(y, X, off) {
      pf <- stats::glm.fit(X, y, offset = off, family = stats::poisson())
      r <- pln_fit(y, X, off, unname(pf$coefficients), nn)
      .result(r$beta, list(sigma = r$sigma), r$mu, r$loglik, r$convergence, r$msg, "own:pln_core")
    },
    logpmf  = function(y, par) pln_logpmf_eta(y, pln_eta(par$mu, par$sigma), par$sigma, nn),
    cdf     = function(q, par) pln_cdf_eta(q, pln_eta(par$mu, par$sigma), par$sigma, P, gl, nn),
    rgen    = function(par) pln_rgen_eta(pln_eta(par$mu, par$sigma), par$sigma),
    mu_fn   = function(beta, extra, X, off) as.vector(exp(X %*% beta + off + extra$sigma^2 / 2)),
    moments = function(par) list(mean = par$mu, var = par$mu + par$mu^2 * (exp(par$sigma^2) - 1))
  )
}

# ---------------------------------------------------------------- CMP (own fitter)
fam_cmp <- function(cfg) {
  K <- cfg$CMP_K; mw <- cfg$CMP_MAX_WINDOW
  nb <- fam_nb(cfg)
  list(
    name = "cmp", n_extra = 1L,
    fit = function(y, X, off) {
      st <- nb$fit(y, X, off)                              # start values (brief 5.6)
      if (!all(is.finite(st$beta))) stop("CMP start: NB start fit has non-finite coefficients")
      r <- cmp_fit(y, X, off, st$beta, st$extra$size, K, mw)
      res <- .result(r$beta, list(nu = r$nu), r$mu, r$loglik, r$convergence, r$msg, "own:cmp_core")
      res$cmp_sol <- r$sol                                  # cached lambda solutions (see cmp_prepare)
      res
    },
    logpmf  = function(y, par) cmp_logpmf(y, par, K, mw),
    cdf     = function(q, par) cmp_cdf(q, par, K, mw),
    rgen    = function(par) cmp_rgen(par, K, mw),
    mu_fn   = function(beta, extra, X, off) .loglink(beta, X, off),
    moments = function(par) {
      par <- cmp_prepare(par, K, mw)
      v <- vapply(seq_along(par$mu), function(i) {
        s <- par$cmp_sol[i, ]; y <- s[["lo"]]:s[["hi"]]
        p <- exp(.cmp_terms(s[["l"]], par$nu, y, cmp_lgfact(y)) - s[["logZ"]])
        sum((y - par$mu[i])^2 * p)
      }, 0)
      list(mean = par$mu, var = v)
    }
  )
}

#' All requested families, by name.
gof_families <- function(cfg = GOF_CONFIG, models = cfg$MODELS) {
  all <- list(poisson = fam_poisson, nb = fam_nb, zinb = fam_zinb,
              genpois = fam_genpois, pln = fam_pln, cmp = fam_cmp)
  lapply(stats::setNames(models, models), function(m) all[[m]](cfg))
}

#' The par list a family's logpmf/cdf/rgen take, from one gene's fit result.
fit_par <- function(fr) {
  p <- c(list(mu = fr$mu), fr$extra)
  if (!is.null(fr$cmp_sol)) p$cmp_sol <- fr$cmp_sol
  p
}
