# plots.R -- figures from saved outputs only (called from 03_report.R).

ensure_fonts <- function() {                      # copied from TCGA_test/scripts/null_calibration.R
  n_fonts <- suppressWarnings(tryCatch(length(system("fc-list", intern = TRUE, ignore.stderr = TRUE)),
                                       error = function(e) 0L))
  if (isTRUE(n_fonts > 0L)) return(invisible())
  message("No system fonts found -- plot text would render as boxes. Installing a font ...")
  cmd <- if (nzchar(Sys.which("pacman")))       "pacman -Sy --noconfirm fontconfig ttf-dejavu"
         else if (nzchar(Sys.which("apt-get"))) "apt-get update && apt-get install -y fontconfig fonts-dejavu-core"
         else if (nzchar(Sys.which("apk")))     "apk add --no-cache fontconfig ttf-dejavu"
         else NA_character_
  if (is.na(cmd)) { warning("No known package manager; install a font + fontconfig in the image."); return(invisible()) }
  try(system(paste(cmd, "&& fc-cache -f"), ignore.stdout = TRUE, ignore.stderr = TRUE), silent = TRUE)
}

#' PIT histogram (50 bins) with a per-bin bootstrap 2.5-97.5% band -- calibrated,
#' unlike a binomial band -- plus the mean-tertile panels.
plot_pit_hist <- function(u, boot_hist, title, strata_u = NULL) {
  bins <- 50L
  h <- pit_hist(u, bins) / length(u) * bins
  d <- data.frame(mid = (seq_len(bins) - 0.5) / bins, dens = h)
  p <- ggplot2::ggplot(d, ggplot2::aes(mid, dens))
  if (!is.null(boot_hist) && nrow(boot_hist)) {
    bh <- boot_hist / rowSums(boot_hist) * bins
    d$lo <- apply(bh, 2, stats::quantile, 0.025); d$hi <- apply(bh, 2, stats::quantile, 0.975)
    p <- p + ggplot2::geom_ribbon(data = d, ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80")
  }
  p + ggplot2::geom_col(width = 1 / bins, fill = NA, colour = "black", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::labs(x = "PIT u", y = "density", title = title,
                  subtitle = "grey = bootstrap 2.5-97.5% per bin") + ggplot2::theme_bw()
}

#' Uniform QQ, thinned to 1e4 quantiles.
plot_pit_qq <- function(u, title) {
  pr <- (seq_len(1e4) - 0.5) / 1e4
  d <- data.frame(theo = pr, emp = stats::quantile(u, pr, names = FALSE))
  ggplot2::ggplot(d, ggplot2::aes(theo, emp)) + ggplot2::geom_line() +
    ggplot2::geom_abline(linetype = 2) + ggplot2::coord_equal() +
    ggplot2::labs(x = "Uniform quantile", y = "PIT quantile", title = title) + ggplot2::theme_bw()
}

#' ECDF - u, with a bootstrap envelope at the 50 bin edges.
plot_pit_ecdf_diff <- function(u, boot_hist, title) {
  e <- (1:50) / 50
  d <- data.frame(u = e, diff = stats::ecdf(u)(e) - e)
  p <- ggplot2::ggplot(d, ggplot2::aes(u, diff))
  if (!is.null(boot_hist) && nrow(boot_hist)) {
    cb <- t(apply(boot_hist, 1, function(h) cumsum(h) / sum(h))) - matrix(e, nrow(boot_hist), 50, byrow = TRUE)
    d$lo <- apply(cb, 2, stats::quantile, 0.025); d$hi <- apply(cb, 2, stats::quantile, 0.975)
    p <- p + ggplot2::geom_ribbon(data = d, ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80")
  }
  p + ggplot2::geom_line() + ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(x = "u", y = "ECDF(u) - u", title = title) + ggplot2::theme_bw()
}

#' ECDF by mean tertile (the mean-stratum panel).
plot_pit_by_mean <- function(U, mean_stratum, title) {
  d <- do.call(rbind, lapply(levels(mean_stratum), function(l) {
    uu <- as.numeric(U[names(mean_stratum)[mean_stratum == l], ])
    e <- seq(0, 1, length.out = 201)
    data.frame(stratum = l, u = e, diff = stats::ecdf(uu)(e) - e)
  }))
  ggplot2::ggplot(d, ggplot2::aes(u, diff, colour = stratum)) + ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(x = "u", y = "ECDF(u) - u", title = title) + ggplot2::theme_bw()
}

#' Per-sample mean of qnorm(u) against its bootstrap envelope.
plot_sample_means <- function(sm_obs, sm_boot, title) {
  d <- data.frame(sample = factor(names(sm_obs), levels = names(sm_obs)), obs = sm_obs,
                  lo = apply(sm_boot, 2, stats::quantile, 0.025), hi = apply(sm_boot, 2, stats::quantile, 0.975))
  d$out <- d$obs < d$lo | d$obs > d$hi
  ggplot2::ggplot(d, ggplot2::aes(sample, obs)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = out)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
    ggplot2::labs(x = NULL, y = "mean qnorm(u)", title = title, colour = "outside") +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_blank())
}
