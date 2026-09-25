# pln_core.R -- Poisson-lognormal (brief 5.5).
#
#   Y | Z ~ Pois(exp(eta + sigma Z)),  Z ~ N(0, 1),  eta = o + x beta.
#
# beta is then on the MEDIAN scale; E[Y] = exp(eta + sigma^2 / 2). Throughout this
# code `mu` means E[Y] (as for every other family), and eta = log(mu) - sigma^2/2.
#
# pmf: ADAPTIVE Gauss-Hermite, as specified (60 nodes by default, not 30 -- see PLN_NODES). The integrand in z,
#   g(z) = -z^2/2 + y (eta + sigma z) - exp(eta + sigma z),
# is strictly concave (g'' = -1 - sigma^2 e^(eta + sigma z) <= -1), so its mode is
# found by Newton from any start, and the nodes are centred on the mode with the
# Laplace scale tau = (-g'')^(-1/2).
#
# cdf: NOT the brief's "sum_k w_k ppois(q, .) on the same adapted nodes". That was
# checked and is numerically wrong where it matters: nodes adapted to the pmf peak
# span only +-~8 tau in z, and when tau is small (large y, moderate sigma) the
# integrand phi(z) ppois(q, .) has most of its mass OUTSIDE that span -- e.g.
# y = 1000, sigma = 0.5 gives tau = 0.063, nodes within |z| < 0.61, while
# F(999) ~ 0.5 is carried by phi over the whole real line. Instead:
#   F(q) = Phi(z*) - int_{-inf}^{z*} phi(z) (1 - S(z)) dz + int_{z*}^{inf} phi(z) S(z) dz,
# with S(z) = ppois(q, e^(eta + sigma z)) and z* the point where e^(eta + sigma z*) = q + 1.
# The Phi(z*) term is exact; the two corrections are smooth one-sided integrands
# that decay away from z*, integrated by composite Gauss-Legendre on panels graded
# towards z* (uniform panels when z* lies outside [-20, 20] and was clipped off).
# Where F < 1e-25 (deep lower tail) it is instead summed exactly from the pmf. That keeps F accurate to ~1e-16 absolute across the whole support,
# which is what the G1 increment check (F(q) - F(q-1) = pmf(q)) needs.

.pln_gh <- new.env()
pln_gh_nodes <- function(n) {
  key <- paste0("gh", n)
  if (is.null(.pln_gh[[key]])) {
    g <- statmod::gauss.quad(n, kind = "hermite")       # int e^{-x^2} f(x) dx
    .pln_gh[[key]] <- list(x = g$nodes, logw = log(g$weights) + g$nodes^2)
  }
  .pln_gh[[key]]
}
pln_gl_nodes <- function(n) {
  key <- paste0("gl", n)
  if (is.null(.pln_gh[[key]])) {
    g <- statmod::gauss.quad(n, kind = "legendre")      # on [-1, 1]
    .pln_gh[[key]] <- list(x = g$nodes, w = g$weights)
  }
  .pln_gh[[key]]
}

pln_eta <- function(mu, sigma) log(mu) - sigma^2 / 2

#' Mode of g(z) for each observation (vectorised Newton; g' is concave and
#' decreasing, so after at most one overshoot the iterates converge monotonically).
.pln_mode <- function(y, eta, sigma) {
  z <- pmin(pmax((log(y + 0.5) - eta) / sigma, -40), 40)
  for (it in seq_len(200L)) {
    e <- exp(eta + sigma * z)
    g1 <- -z + sigma * y - sigma * e
    g2 <- -1 - sigma^2 * e
    step <- g1 / g2
    z <- z - step
    if (all(abs(step) < 1e-12 * pmax(1, abs(z)))) break
  }
  list(z = z, g2 = -1 - sigma^2 * exp(eta + sigma * z))
}

#' log p(y) by adaptive Gauss-Hermite. y, eta vectors (same length); sigma scalar.
# Evaluated in row blocks so the (obs x nodes) matrices stay small on modest machines.
PLN_BLOCK <- 5000L
.pln_blocks <- function(n) split(seq_len(n), ceiling(seq_len(n) / PLN_BLOCK))

pln_logpmf_eta <- function(y, eta, sigma, n_nodes = 60L, grad = FALSE) {
  if (length(y) > PLN_BLOCK) {
    parts <- lapply(.pln_blocks(length(y)), function(i) pln_logpmf_eta(y[i], eta[i], sigma, n_nodes, grad))
    if (!grad) return(unlist(parts, use.names = FALSE))
    return(list(lp = unlist(lapply(parts, `[[`, "lp")), d_eta = unlist(lapply(parts, `[[`, "d_eta")),
                d_sigma = unlist(lapply(parts, `[[`, "d_sigma"))))
  }
  gh <- pln_gh_nodes(n_nodes)
  md <- .pln_mode(y, eta, sigma)
  tau <- 1 / sqrt(-md$g2)
  Z <- md$z + sqrt(2) * tau %o% gh$x                     # obs x nodes
  Lam <- exp(eta + sigma * Z)
  L <- -Z^2 / 2 - 0.5 * log(2 * pi) + y * (eta + sigma * Z) - Lam - lgamma(y + 1)
  L <- sweep(L, 2, gh$logw, "+")
  mx <- apply(L, 1, max)
  E <- exp(L - mx); S <- rowSums(E)
  lp <- log(sqrt(2) * tau) + mx + log(S)
  if (!grad) return(lp)
  # Posterior moments on the same nodes give the exact score:
  #   d log p / d eta = y - E[lambda | y],  d log p / d sigma = E[z (y - lambda) | y].
  Wn <- E / S
  list(lp = lp, d_eta = y - rowSums(Wn * Lam), d_sigma = rowSums(Wn * Z * (y - Lam)))
}

#' Composite Gauss-Legendre nodes/weights on [A, B] per row, panels graded
#' quadratically towards the `fine` end (A if fine_at_A, else B).
.gl_panels <- function(A, B, fine_at_A, P, n, graded) {
  gl <- pln_gl_nodes(n)
  len <- B - A
  eg <- ((0:P) / P)^2                                     # graded edges (fine at the z* end)
  eu <- (0:P) / P                                         # uniform edges (z* clipped off the interval)
  Zs <- vector("list", P); Ws <- vector("list", P)
  for (k in seq_len(P)) {
    a <- ifelse(graded, eg[k], eu[k]); b <- ifelse(graded, eg[k + 1], eu[k + 1])   # per row
    t <- (a + b) / 2 + ((b - a) / 2) %o% gl$x             # rows x nodes, in [0, 1]
    wt <- ((b - a) / 2) %o% gl$w
    Zs[[k]] <- if (fine_at_A) A + len * t else B - len * t
    Ws[[k]] <- len * wt
  }
  list(Z = do.call(cbind, Zs), W = do.call(cbind, Ws))
}

#' F(q) for PLN. q, eta vectors; sigma scalar.
PLN_ZLIM <- 20            # phi(20) ~ 1e-88: the z-range the GL corrections cover
PLN_SUM_BELOW <- 1e-25   # below this F is summed from the pmf (see pln_cdf_eta)

pln_cdf_eta <- function(q, eta, sigma, P = 12L, n = 16L, n_nodes = 60L) {
  if (length(q) > PLN_BLOCK)
    return(unlist(lapply(.pln_blocks(length(q)), function(i) pln_cdf_eta(q[i], eta[i], sigma, P, n, n_nodes)),
                  use.names = FALSE))
  q <- floor(q)
  out <- numeric(length(q))
  pos <- q >= 0
  if (!any(pos)) return(out)
  qq <- q[pos]; ee <- eta[pos]
  zs <- (log(qq + 1) - ee) / sigma
  w <- 1 / (sigma * sqrt(qq + 1))                         # transition width in z (large q)
  # Left correction (1 - S) decays like a Gaussian of sd w for large q, but only
  # exponentially, at rate (q + 1) sigma, for small q -- take the wider reach.
  reachL <- pmax(12 * w, 80 / ((qq + 1) * sigma))
  reachR <- 12 * w                                       # S -> 0 at least as fast (double-exponentially for small q)
  lam <- function(Z, e) exp(e + sigma * Z)
  # left: [max(-ZLIM, z* - reachL), min(z*, ZLIM)]
  A <- pmax(-PLN_ZLIM, zs - reachL); B <- pmin(zs, PLN_ZLIM)
  CL <- numeric(length(qq)); iL <- which(B > A)
  if (length(iL)) {
    g <- .gl_panels(A[iL], B[iL], fine_at_A = FALSE, P, n, graded = B[iL] == zs[iL])
    S1 <- stats::ppois(qq[iL], lam(g$Z, ee[iL]), lower.tail = FALSE)
    CL[iL] <- rowSums(g$W * stats::dnorm(g$Z) * S1)
  }
  # right: [max(z*, -ZLIM), min(ZLIM, z* + reachR)]
  A <- pmax(zs, -PLN_ZLIM); B <- pmin(PLN_ZLIM, zs + reachR)
  CR <- numeric(length(qq)); iR <- which(B > A)
  if (length(iR)) {
    g <- .gl_panels(A[iR], B[iR], fine_at_A = TRUE, P, n, graded = A[iR] == zs[iR])
    S <- stats::ppois(qq[iR], lam(g$Z, ee[iR]))
    CR[iR] <- rowSums(g$W * stats::dnorm(g$Z) * S)
  }
  Fq <- pmin(pmax(stats::pnorm(zs) - CL + CR, 0), 1)
  # Deep lower tail: the integrand's mass can sit beyond z = -ZLIM (e.g. y = 0 for a
  # gene with mu ~ 1e3), where the window above truncates it. There F is tiny and q
  # is below the bulk, so F(q) = sum_{j <= q} p(j) is both exact and cheap.
  deep <- which(Fq < PLN_SUM_BELOW)
  for (i in deep) {
    j <- 0:qq[i]
    Fq[i] <- sum(exp(pln_logpmf_eta(j, rep(ee[i], length(j)), sigma, n_nodes)))
  }
  out[pos] <- Fq
  out
}

pln_rgen_eta <- function(eta, sigma)
  stats::rpois(length(eta), exp(eta + sigma * stats::rnorm(length(eta))))

#' Own ML fit: nlminb on (beta, log sigma) with the analytic gradient, started from
#' the Poisson fit and sigma = 0.3 (brief 5.5). log sigma is box-constrained to [log 1e-4, log 10]
#' for numerical safety; a fit ending on that box is flagged in conv_msg.
#' Note the Poisson start is on the MEAN scale; the PLN intercept is on the median
#' scale, so the start is shifted by -sigma0^2/2 (the offset enters eta unchanged).
pln_fit <- function(y, X, off, beta_pois, n_nodes = 60L) {
  p <- ncol(X)
  s0 <- 0.3
  b0 <- beta_pois
  ic <- which(apply(X, 2, function(v) all(v == 1)))
  if (length(ic)) b0[ic[1]] <- b0[ic[1]] - s0^2 / 2
  # Objective and analytic gradient share one quadrature pass (cached by theta).
  cache <- new.env()
  eval_at <- function(th) {
    if (!identical(cache$th, th)) {
      eta <- as.vector(X %*% th[seq_len(p)] + off); sg <- exp(th[p + 1L])
      r <- pln_logpmf_eta(y, eta, sg, n_nodes, grad = TRUE)
      cache$th <- th
      cache$f <- -sum(r$lp)
      cache$g <- -c(crossprod(X, r$d_eta), sg * sum(r$d_sigma))
    }
    cache
  }
  negll <- function(th) { v <- eval_at(th)$f; if (is.finite(v)) v else Inf }
  grad  <- function(th) eval_at(th)$g
  lb <- c(rep(-Inf, p), log(1e-4)); ub <- c(rep(Inf, p), log(10))
  opt <- stats::nlminb(c(b0, log(s0)), negll, grad, lower = lb, upper = ub,
                       control = list(eval.max = 600L, iter.max = 400L))
  beta <- opt$par[seq_len(p)]; sigma <- exp(opt$par[p + 1L])
  eta <- as.vector(X %*% beta + off)
  at_bound <- opt$par[p + 1L] <= lb[p + 1L] + 1e-8 || opt$par[p + 1L] >= ub[p + 1L] - 1e-8
  list(beta = beta, sigma = sigma, mu = exp(eta + sigma^2 / 2), loglik = -opt$objective,
       convergence = opt$convergence,
       msg = paste0(opt$message, if (at_bound) " | log(sigma) on its numerical box" else ""))
}
