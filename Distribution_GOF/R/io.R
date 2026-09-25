# io.R -- input contract (brief section 3), builders, integer audit, provenance.

DATASET_SOURCES <- c("salmon_tximport", "featureCounts", "tcga_counts_rds", "simulated")

#' Assemble a dataset object and check it against the input contract.
new_dataset <- function(counts, samples, design, lengths = NULL, source, label) {
  ds <- list(counts = counts, samples = samples, design = design,
             lengths = lengths, source = source, label = label)
  validate_dataset(ds)
}

validate_dataset <- function(ds) {
  need <- c("counts", "samples", "design", "lengths", "source", "label")
  miss <- setdiff(need, names(ds))
  if (length(miss)) stop("dataset is missing: ", paste(miss, collapse = ", "))
  Y <- ds$counts
  if (!is.matrix(Y) || !is.numeric(Y)) stop("counts must be a numeric matrix (genes x samples)")
  if (is.null(rownames(Y)) || anyDuplicated(rownames(Y))) stop("counts needs unique rownames (gene ids)")
  if (is.null(colnames(Y)) || anyDuplicated(colnames(Y))) stop("counts needs unique colnames (sample ids)")
  if (any(!is.finite(Y))) stop("counts contains non-finite values")
  if (any(Y < 0)) stop("counts contains negative values")
  if (!is.data.frame(ds$samples) || nrow(ds$samples) != ncol(Y))
    stop("samples must be a data.frame with one row per column of counts")
  if (!is.null(ds$samples$sample_id) && !identical(as.character(ds$samples$sample_id), colnames(Y)))
    stop("samples$sample_id does not match colnames(counts) (same order required)")
  if (!inherits(ds$design, "formula")) stop("design must be a formula, e.g. ~ condition")
  dv <- all.vars(ds$design)
  if (length(setdiff(dv, names(ds$samples))))
    stop("design uses columns not in samples: ", paste(setdiff(dv, names(ds$samples)), collapse = ", "))
  if (!is.null(ds$lengths)) {
    L <- ds$lengths
    if (!is.matrix(L) || !identical(dim(L), dim(Y)) || !identical(dimnames(L), dimnames(Y)))
      stop("lengths must be a matrix with the same dim and dimnames as counts")
    if (any(!is.finite(L)) || any(L <= 0)) stop("lengths must be finite and > 0")
  }
  if (!(ds$source %in% DATASET_SOURCES))
    stop("source must be one of ", paste(DATASET_SOURCES, collapse = " / "))
  if (!is.character(ds$label) || !grepl("^[A-Za-z0-9._-]+$", ds$label))
    stop("label must be a short id of [A-Za-z0-9._-]")
  ds
}

read_dataset <- function(path) {
  if (!nzchar(path) || !file.exists(path)) stop("dataset file not found: '", path, "'")
  validate_dataset(readRDS(path))
}

#' Drop samples by id (decision D6). Unknown ids are an error, not a no-op.
apply_exclusions <- function(ds, ids) {
  if (!length(ids)) return(ds)
  unk <- setdiff(ids, colnames(ds$counts))
  if (length(unk)) stop("EXCLUDE_SAMPLES contains unknown sample ids: ", paste(unk, collapse = ", "))
  keep <- !(colnames(ds$counts) %in% ids)
  ds$counts  <- ds$counts[, keep, drop = FALSE]
  ds$samples <- ds$samples[keep, , drop = FALSE]
  if (!is.null(ds$lengths)) ds$lengths <- ds$lengths[, keep, drop = FALSE]
  validate_dataset(ds)
}

# ---------------------------------------------------------------- builders

#' Salmon quantifications via tximport. The offset recipe is NOT reimplemented
#' here: the installed tximport vignette (1.40.0) delegates it to
#' edgeR::DGEListFromTximport(txi), which compute_offsets() calls for
#' OFFSET_MODE = "tmm_length". So this builder only keeps $counts and $length.
build_from_tximport <- function(quant_files, tx2gene, samples, design, label) {
  if (!requireNamespace("tximport", quietly = TRUE)) load_or_install("tximport")
  if (is.null(names(quant_files))) stop("quant_files must be named by sample id")
  txi <- tximport::tximport(quant_files, type = "salmon", tx2gene = tx2gene,
                            countsFromAbundance = "no")
  new_dataset(counts = txi$counts, samples = samples, design = design,
              lengths = txi$length, source = "salmon_tximport", label = label)
}

#' featureCounts output table. Sample columns are matched to samples$sample_id
#' after stripping the directory and a trailing ".bam"; every column must match.
build_from_featurecounts <- function(file, samples, design, label) {
  if (is.null(samples$sample_id)) stop("samples needs a sample_id column for featureCounts input")
  tab <- data.table::fread(file, skip = "Geneid")
  fixed <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  if (!all(fixed %in% names(tab))) stop("not a featureCounts table: missing ", setdiff(fixed, names(tab)))
  sc <- setdiff(names(tab), fixed)
  ids <- sub("\\.bam$", "", basename(sc))
  if (!setequal(ids, samples$sample_id))
    stop("featureCounts columns do not match samples$sample_id: ",
         paste(setdiff(union(ids, samples$sample_id), intersect(ids, samples$sample_id)), collapse = ", "))
  Y <- as.matrix(tab[, sc, with = FALSE])
  dimnames(Y) <- list(tab$Geneid, ids)
  storage.mode(Y) <- "double"
  Y <- Y[, as.character(samples$sample_id), drop = FALSE]
  new_dataset(counts = Y, samples = samples, design = design, lengths = NULL,
              source = "featureCounts", label = label)
}

#' Raw COUNTS (genes x samples) for a TCGA cancer stage OR the matched-normal cohort.
#' Copied from TCGA_test/scripts/count_distribution.R (that script cannot be
#' sourced: it runs its analysis on load). BASE_DATA_DIR comes from the repo's
#' config.R, looked up in the working directory first, then in common/.
load_counts_matrix <- function(project_id, stage) {
  fn <- if (identical(stage, "healthy"))
          file.path(BASE_DATA_DIR, project_id, paste0("healthy_", project_id, "_counts.rds"))
        else
          file.path(BASE_DATA_DIR, project_id,
                    paste0(project_id, "-", gsub(" ", "-", stage), "_counts.rds"))
  if (!file.exists(fn)) return(NULL)
  obj <- readRDS(fn)
  if (!is.matrix(obj$expression_vectors)) return(NULL)
  m <- t(obj$expression_vectors)               # samples x genes -> genes x samples
  colnames(m) <- rownames(obj$expression_vectors)
  rownames(m) <- if (!is.null(obj$gene_ids)) obj$gene_ids else colnames(obj$expression_vectors)
  m
}

#' A TCGA cohort as a single group (design ~ 1), as in count_distribution.R.
#' Unmodelled between-patient heterogeneity then counts as lack of fit for every
#' model (README states this).
build_from_tcga <- function(project_id, stage, label) {
  if (!exists("BASE_DATA_DIR")) {
    cand <- c("config.R", file.path(GOF_ROOT, "..", "common", "config.R"))
    cand <- cand[file.exists(cand)]
    if (!length(cand)) stop("config.R (BASE_DATA_DIR) not found in the working dir or common/")
    source(cand[1])
  }
  Y <- load_counts_matrix(project_id, stage)
  if (is.null(Y)) stop("no counts for ", project_id, " / ", stage)
  storage.mode(Y) <- "double"
  new_dataset(counts = Y, samples = data.frame(sample_id = colnames(Y)), design = ~ 1,
              lengths = NULL, source = "tcga_counts_rds", label = label)
}

# ---------------------------------------------------------------- integer audit

#' Integer audit (brief section 3), always run. Long-format table:
#' section / stratum / metric / value. Nothing is assumed about the input.
integer_audit <- function(Y) {
  frac <- Y - floor(Y)
  nonint <- abs(Y - round(Y)) > 1e-8
  rows <- list()
  add <- function(section, stratum, metric, value)
    rows[[length(rows) + 1L]] <<- data.table::data.table(section = section, stratum = stratum,
                                                         metric = metric, value = as.numeric(value))
  add("overall", "all", "n_entries", length(Y))
  add("overall", "all", "frac_non_integer", mean(nonint))
  add("overall", "all", "n_non_integer", sum(nonint))
  add("overall", "all", "min_count", min(Y))
  add("overall", "all", "max_count", max(Y))
  ls <- colSums(Y)
  add("overall", "all", "lib_size_min", min(ls))
  add("overall", "all", "lib_size_median", stats::median(ls))
  add("overall", "all", "lib_size_max", max(ls))
  fq <- if (any(nonint)) stats::quantile(frac[nonint], c(0, .1, .25, .5, .75, .9, 1)) else rep(NA, 7)
  for (i in seq_along(fq)) add("fractional_part", "non_integer_entries",
                               paste0("q", c(0, 10, 25, 50, 75, 90, 100)[i]), fq[i])
  h <- if (any(nonint)) graphics::hist(frac[nonint], breaks = seq(0, 1, 0.1), plot = FALSE)$counts else rep(0, 10)
  for (i in seq_along(h)) add("fractional_part_hist", sprintf("[%.1f,%.1f)", (i - 1) / 10, i / 10), "count", h[i])
  dec <- ceiling(10 * rank(rowMeans(Y), ties.method = "first") / nrow(Y))
  for (d in sort(unique(dec))) {
    s <- dec == d
    add("by_mean_decile", sprintf("D%02d", d), "frac_non_integer", mean(nonint[s, ]))
    add("by_mean_decile", sprintf("D%02d", d), "median_frac_part_non_integer",
        if (any(nonint[s, ])) stats::median(frac[s, ][nonint[s, ]]) else NA)
    add("by_mean_decile", sprintf("D%02d", d), "mean_count", mean(Y[s, ]))
  }
  data.table::rbindlist(rows)
}

audit_has_fractional <- function(audit)
  audit[section == "overall" & metric == "n_non_integer", value] > 0

# ---------------------------------------------------------------- seeds & provenance

#' Deterministic integer seed from BASE_SEED and a tag, e.g.
#' seed_for(188, "boot", "nb", 17). Lets any single stochastic step (one model's
#' one bootstrap replicate, one PIT draw) be rerun on its own.
seed_for <- function(base, ...) {
  s <- paste(c(as.character(base), vapply(list(...), as.character, "")), collapse = "|")
  h <- 0
  for (b in utf8ToInt(s)) h <- (h * 131 + b) %% 2147483647
  as.integer(h)
}

#' Run `expr` under set.seed(seed) and restore the caller's RNG state afterwards.
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
          else rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  expr
}

#' Provenance for an output directory (brief rule 8): sessionInfo, package
#' versions, the full resolved config, RNG kind and base seed, git HEAD.
write_provenance <- function(dir, cfg = GOF_CONFIG, extra = list()) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
  head <- Sys.getenv("GOF_GIT_HEAD", "")
  if (!nzchar(head))
    head <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(GOF_ROOT), "rev-parse", "HEAD"),
                                              stdout = TRUE, stderr = FALSE)),
                     error = function(e) character(0))
  if (!length(head) || !nzchar(head[1])) head <- "unavailable (git not reachable; set GOF_GIT_HEAD)"
  lines <- c(
    sprintf("written:      %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("git HEAD:     %s", head[1]),
    sprintf("R:            %s", R.version.string),
    sprintf("RNG kind:     %s", paste(RNGkind(), collapse = " / ")),
    sprintf("BASE_SEED:    %s", cfg$BASE_SEED),
    "packages:",
    sprintf("  %-10s %s", c("glmmTMB", "TMB", "edgeR", "tximport", "statmod", "data.table"),
            vapply(c("glmmTMB", "TMB", "edgeR", "tximport", "statmod", "data.table"), pv, "")),
    "resolved config:",
    sprintf("  %-17s = %s", names(cfg), vapply(cfg, function(v)
      if (length(v)) paste(as.character(v), collapse = ", ") else "<none>", "")))
  if (length(extra)) lines <- c(lines, "run facts:",
                                sprintf("  %-17s = %s", names(extra), vapply(extra, function(v)
                                  paste(as.character(v), collapse = ", "), "")))
  writeLines(lines, file.path(dir, "provenance.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(dir, "sessionInfo.txt"))
  saveRDS(cfg, file.path(dir, "config_resolved.rds"))
  invisible(dir)
}

#' A simulated NB dataset in the contract's format (source = "simulated"), for
#' development runs and for checking the pipeline end to end without real data.
#' Two conditions, per-gene NB with size ~ lognormal around 5, means log-uniform
#' over [1, 1e4] counts per 1e7 reads, 10% of genes changed 1.5-fold in B.
build_simulated <- function(n_genes, n_per_group, seed, label = "simNB") {
  with_seed(seed, {
    n <- 2L * n_per_group
    cond <- factor(rep(c("A", "B"), each = n_per_group))
    lib <- stats::runif(n, 0.5, 2) * 1e7
    base <- exp(stats::runif(n_genes, log(1), log(1e4)))
    size <- exp(stats::rnorm(n_genes, log(5), 0.5))
    fc <- ifelse(stats::runif(n_genes) < 0.1, 1.5, 1)
    Y <- t(vapply(seq_len(n_genes), function(g)
      as.numeric(stats::rnbinom(n, mu = base[g] * lib / 1e7 * ifelse(cond == "B", fc[g], 1), size = size[g])),
      numeric(n)))
  })
  dimnames(Y) <- list(sprintf("gene%05d", seq_len(n_genes)), sprintf("s%03d", seq_len(2L * n_per_group)))
  new_dataset(counts = Y, samples = data.frame(sample_id = colnames(Y), condition = cond),
              design = ~ condition, lengths = NULL, source = "simulated", label = label)
}
