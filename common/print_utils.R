# print_utils.R
# -----------------------------------------------------------------------------
# Console tables for the benchmark scripts (implementation note C). Every table a
# script prints goes through print_tab() or print_pivot(), so that
#   * titles and tables share one stream (stdout, cat), and stay together in a log;
#   * a table wider than getOption("width") is split into column chunks that ALL
#     repeat the key columns (R3), so no row is ever separated from its keys;
#   * the console is a capped summary: all-NA columns are dropped, constant columns
#     move into one "constant:" line under the title, rows are capped at `n` with a
#     pointer to the CSV that holds the full table (R5).
# No package dependency beyond base R.
# -----------------------------------------------------------------------------

PRINT_KEYS <- c("part", "dist", "n_rep", "pi1", "het", "frac_up", "depth", "disp_model", "hv_frac",
                "sigma_d", "delta_target", "cohort", "method")

#' Columns with a single distinct non-NA value, and their values as "name=value".
.const_cols <- function(x, keep = "method") {
  cc <- names(x)[vapply(x, function(v) length(unique(v[!is.na(v)])) <= 1L, logical(1))]
  cc <- setdiff(cc, keep)
  list(cols = cc, text = if (length(cc))
    paste(cc, vapply(x[cc], function(v) as.character(v[!is.na(v)][1]), ""), sep = "=", collapse = ", "))
}

#' Print `x` in column chunks: every chunk = keys + as many value columns as fit in
#' `width` (and at most `max_cols`). R's own wrapping is switched off inside, so a
#' chunk is never split a second time.
.print_chunked <- function(x, keys, width = getOption("width"), max_cols = Inf) {
  fmt_w <- function(v, nm) max(nchar(format(v)), nchar(nm), 0L) + 1L
  vals <- setdiff(names(x), keys)
  kw <- sum(vapply(keys, function(k) fmt_w(x[[k]], k), 0L))
  chunks <- list(); cur <- character(0); used <- kw
  for (v in vals) {
    w <- fmt_w(x[[v]], v)
    if (length(cur) && (used + w > width || length(cur) >= max_cols)) {
      chunks[[length(chunks) + 1L]] <- cur; cur <- character(0); used <- kw
    }
    cur <- c(cur, v); used <- used + w
  }
  if (length(cur) || !length(chunks)) chunks[[length(chunks) + 1L]] <- cur
  old <- options(width = 10000L); on.exit(options(old))
  for (i in seq_along(chunks)) {
    if (i > 1L) cat(strrep("-", 20), "\n", sep = "")
    print(x[, c(keys, chunks[[i]]), drop = FALSE], row.names = FALSE)
  }
}

.title <- function(title, width) cat("\n", strrep("=", min(width, 100L)), "\n", title, "\n", sep = "")

#' One console table (note C section 3.1).
#'   csv   file name of the full table (shown when rows are cut)
#'   n     row cap
#'   keys  candidate key columns, in order; those present and non-constant come
#'         first, are sorted on, and are repeated in every column chunk
print_tab <- function(x, title, csv = NULL, n = 60L, notes = NULL, keys = PRINT_KEYS,
                      width = getOption("width")) {
  .title(title, width)
  if (is.null(x) || !nrow(x)) { cat("(empty)\n"); return(invisible(x)) }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- x[, vapply(x, function(v) !all(is.na(v)), logical(1)), drop = FALSE]
  cst <- if (nrow(x) > 1L) .const_cols(x) else list(cols = character(0), text = NULL)
  if (length(cst$cols)) cat("constant: ", cst$text, "\n", sep = "")
  for (nt in notes) cat("  ", nt, "\n", sep = "")
  x <- x[, setdiff(names(x), cst$cols), drop = FALSE]
  keys <- intersect(keys, names(x))
  num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], signif, 4)
  if (length(keys)) x <- x[do.call(order, unname(as.list(x[keys]))), c(keys, setdiff(names(x), keys)), drop = FALSE]
  total <- nrow(x)
  .print_chunked(head(x, n), keys, width)
  if (total > n) cat(sprintf("(showing %d of %d rows%s)\n", n, total,
                             if (length(csv)) paste0("; full table: ", csv) else ""))
  invisible(x)
}

#' Methods side by side (note C section 3.2): rows = cell, columns = methods, one
#' metric. `S` holds ONE part. Row keys default to the non-constant cell keys;
#' constant ones go into the "constant:" line. Methods keep the order of `methods`
#' (default: order of appearance) and are chunked at `max_methods` per block.
print_pivot <- function(S, metric, title, rows = NULL, methods = NULL, digits = 3, csv = NULL,
                        max_methods = 12L, notes = NULL, width = getOption("width")) {
  .title(title, width)
  if (is.null(S) || !nrow(S) || !metric %in% names(S)) { cat("(empty)\n"); return(invisible(NULL)) }
  S <- as.data.frame(S, stringsAsFactors = FALSE)
  methods <- intersect(methods %||% unique(S$method), unique(S$method))
  S <- S[S$method %in% methods, , drop = FALSE]
  if (is.null(rows)) rows <- setdiff(intersect(PRINT_KEYS, names(S)), "method")
  rows <- intersect(rows, names(S))
  cst <- .const_cols(S[rows], keep = character(0))
  if (length(cst$cols)) cat("constant: ", cst$text, "\n", sep = "")
  for (nt in notes) cat("  ", nt, "\n", sep = "")
  rows <- setdiff(rows, cst$cols)
  rk <- if (length(rows)) do.call(paste, c(lapply(S[rows], as.character), sep = "\r")) else rep("", nrow(S))
  R <- S[!duplicated(rk), rows, drop = FALSE]; ukey <- rk[!duplicated(rk)]
  if (length(rows)) { o <- do.call(order, unname(as.list(R))); R <- R[o, , drop = FALSE]; ukey <- ukey[o] }
  out <- R
  for (m in methods) {
    s <- S[S$method == m, , drop = FALSE]
    out[[m]] <- signif(s[[metric]][match(ukey, rk[S$method == m])], digits)
  }
  if (!length(rows)) out <- out[, methods, drop = FALSE]
  .print_chunked(out, rows, width, max_cols = max_methods)
  if (length(csv)) cat("(", metric, "; SE and all metrics: ", csv, ")\n", sep = "")
  invisible(out)
}

#' Row-bind data frames whose columns differ (missing ones filled with NA).
rbind_fill <- function(xs) {
  xs <- Filter(Negate(is.null), xs)
  cols <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
