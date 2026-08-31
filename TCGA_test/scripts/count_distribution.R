#!/usr/bin/env Rscript
# count_distribution.R
# -----------------------------------------------------------------------------
# Is the raw TCGA count data Poisson or Negative Binomial?
#
# This is a property of the DATA, independent of the null-calibration experiment
# (no A/B splits here) -- so it lives in its own script. For each adequately
# expressed gene, in each cohort (the 4 cancer stages + healthy, across all 8
# cancers), we fit an intercept-only Poisson vs NB GLM with a log(library-size)
# offset and ask which fits better. RNA-seq is expected to be overdispersed => NB,
# which is exactly why the reference tools (edgeR / DESeq2) use NB, not Poisson.
#
# Pure R + edgeR only -- deliberately does NOT depend on the TOX Fortran/rcpp build.
# Output: count_distribution_poisson_vs_nb.csv + a printed per-cohort summary.
# -----------------------------------------------------------------------------

LIB_DIR <- normalizePath("external/docker_r_libs", mustWork = FALSE)
if (!dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE) && !dir.exists(LIB_DIR))
  stop("Could not create package library ", LIB_DIR,
       " -- it must be on a WRITABLE, BIND-MOUNTED path or packages will not persist.")
.libPaths(c(LIB_DIR, .libPaths()))

# Name the package that is ACTUALLY missing from a load error (see null_calibration.R
# for the full rationale) so we install the missing dependency, not the top package.
.missing_pkg <- function(msg) {
  q <- "[‘’'\"`]([[:alnum:]._]+)[‘’'\"`]"
  for (pat in c(paste0("there is no package called ", q),
                paste0("package ", q, " required by"))) {
    m <- regmatches(msg, regexec(pat, msg, perl = TRUE))[[1]]
    if (length(m) >= 2L) return(m[2L])
  }
  NA_character_
}
load_or_install <- function(package_name) {
  install_pkg <- function(pkg) {
    if (requireNamespace("BiocManager", quietly = TRUE))
      BiocManager::install(pkg, lib = LIB_DIR, update = FALSE, ask = FALSE)
    else
      install.packages(pkg, repos = "https://cloud.r-project.org", lib = LIB_DIR, dependencies = TRUE)
  }
  if (!requireNamespace(package_name, quietly = TRUE)) install_pkg(package_name)
  for (attempt in seq_len(8L)) {
    err <- tryCatch({ suppressPackageStartupMessages(library(package_name, character.only = TRUE))
                      return(invisible()) },
                    error = function(e) e)
    miss <- .missing_pkg(conditionMessage(err))
    if (is.na(miss) || identical(miss, package_name))
      stop("Could not load '", package_name, "': ", conditionMessage(err))
    message("  load_or_install: '", package_name, "' -> installing missing dependency '", miss, "' ...")
    install_pkg(miss)
  }
  stop("Could not load '", package_name, "' -- unresolved dependencies after 8 attempts.")
}

load_or_install("BiocManager")
load_or_install("edgeR")                 # DGEList / filterByExpr / normLibSizes / estimateDisp
load_or_install("data.table")            # rbindlist
load_or_install("parallel")              # mclapply -- per-gene fits in parallel

# MASS provides glm.nb. Ensure it is installed but do NOT attach it (keeps the namespace clean);
# we call MASS::glm.nb fully qualified.
if (!requireNamespace("MASS", quietly = TRUE))
  install.packages("MASS", repos = "https://cloud.r-project.org", lib = LIB_DIR)

# Each fork must run single-threaded so N_CORES forks don't oversubscribe.
Sys.setenv(OMP_NUM_THREADS = "1")

source("config.R")                       # BASE_DATA_DIR, STAGES, TOX_TEST_DIR (no rcpp/TOX pulled in)
options(width = 200)

# Friendly-name -> TCGA project map. Kept INLINE (mirrors CANCER_TYPES in
# common/outlier_significance_analysis.R) so this pure-R diagnostic does not have to
# source outlier_significance_analysis.R, which would drag in the TOX rcpp build.
CANCER_TYPES <- c(
  "breast cancer" = "TCGA-BRCA",
  "non small cell lung cancer" = "TCGA-LUAD",
  "small cell lung cancer" = "TCGA-LUSC",
  "bladder cancer" = "TCGA-BLCA",
  "colon cancer" = "TCGA-COAD",
  "kidney cancer" = "TCGA-KIRC",
  "stomach cancer" = "TCGA-STAD",
  "thyroid gland cancer" = "TCGA-THCA"
)

# ==================== CONFIG ====================
DIST_CANCERS   <- CANCER_TYPES           # all 8; subset for a quick run
DIST_STAGES    <- STAGES                 # the 4 cancer stages
DIST_HEALTHY   <- TRUE                   # also test the matched-normal (healthy) cohort per cancer
DIST_MIN_MEAN  <- 5                      # min mean raw count for a gene to be testable
DIST_MAX_GENES <- 1500L                  # cap tested genes per cohort (glm.nb is iterative)
N_CORES        <- 32L                    # cores for the parallel per-gene fits

OUT_DIR <- file.path(dirname(TOX_TEST_DIR), "count_distribution")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==================== data loader ====================

#' Load raw COUNTS (genes x samples) for a cancer stage OR the matched-normal cohort.
#' `stage == "healthy"` reads healthy_<pid>_counts.rds; any other value reads
#' <pid>-<stage>_counts.rds. (Same loader as null_calibration.R.)
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

# ==================== the test ====================
#
# For each adequately-expressed gene we fit an intercept-only Poisson and NB GLM, EACH
# with a log(library-size) offset so that sequencing-depth differences between samples
# do NOT masquerade as overdispersion, then compare them three ways:
#   * LRT  = 2*(logLik_NB - logLik_Poisson). Poisson is NB with phi = 0, which lies on
#     the boundary of the parameter space, so under H0 (Poisson) the statistic follows a
#     50:50 mixture of chi^2_0 and chi^2_1 -> p = 0.5 * P(chi^2_1 > LRT). Small p => NB.
#   * AIC  : NB is preferred when AIC_NB < AIC_Poisson (this penalises NB's extra param).
#   * phi  = 1/theta: the fitted NB dispersion itself (0 = Poisson); BCV = sqrt(phi) is
#     the biological coefficient of variation, the standard RNA-seq overdispersion scale.
# edgeR's common dispersion for the whole cohort is reported as an independent cross-check.
test_count_distribution <- function(counts, cancer, stage,
                                    min_mean = DIST_MIN_MEAN, max_genes = DIST_MAX_GENES,
                                    n_cores = N_CORES) {
  if (is.null(counts) || ncol(counts) < 3L) return(NULL)
  libsize <- colSums(counts)
  if (any(libsize <= 0)) return(NULL)
  logls <- log(libsize)

  testable <- which(rowMeans(counts) >= min_mean)
  if (length(testable) < 10L) return(NULL)
  # Deterministic, evenly-spaced thinning (no RNG) to cap the per-gene glm.nb cost.
  if (length(testable) > max_genes)
    testable <- testable[round(seq(1, length(testable), length.out = max_genes))]

  per_gene <- function(g) {
    y <- as.numeric(counts[g, ])
    if (stats::var(y) == 0) return(NULL)
    fitP  <- tryCatch(glm(y ~ offset(logls), family = poisson()), error = function(e) NULL)
    if (is.null(fitP)) return(NULL)
    fitNB <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ offset(logls))), error = function(e) NULL)
    if (is.null(fitNB)) return(NULL)
    llP <- as.numeric(logLik(fitP)); llNB <- as.numeric(logLik(fitNB))
    lrt <- 2 * (llNB - llP)
    p   <- if (lrt <= 0) 1 else 0.5 * pchisq(lrt, df = 1, lower.tail = FALSE)
    c(p = p, nb_aic = as.numeric(AIC(fitNB) < AIC(fitP)), phi = 1 / fitNB$theta)
  }
  fits <- parallel::mclapply(testable, per_gene, mc.cores = n_cores, mc.preschedule = TRUE)
  fits <- do.call(rbind, Filter(function(x) is.numeric(x) && length(x) == 3L, fits))
  if (is.null(fits) || nrow(fits) < 10L) return(NULL)

  # edgeR common dispersion cross-check (one robust cohort-level phi via the NB GLM /
  # qCML framework, library-size-normalised by TMM). NA if it fails.
  edger_phi <- tryCatch(suppressWarnings(suppressMessages({
    y <- DGEList(counts = counts)
    y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
    y <- normLibSizes(y)
    estimateDisp(y)$common.dispersion
  })), error = function(e) NA_real_)

  frac_fdr <- mean(p.adjust(fits[, "p"], method = "BH") < 0.05)
  data.frame(
    cancer = cancer, stage = stage, n_samples = ncol(counts),
    n_genes_tested   = nrow(fits),
    frac_reject_pois = round(mean(fits[, "p"] < 0.05), 4),  # raw p<0.05: ~1 => NB, ~0 => Poisson
    frac_reject_fdr  = round(frac_fdr, 4),                  # BH-FDR<0.05 (stricter)
    frac_nb_by_aic   = round(mean(fits[, "nb_aic"]), 4),    # NB preferred by AIC
    median_phi       = round(median(fits[, "phi"]), 4),     # NB dispersion (0 = Poisson)
    median_bcv       = round(median(sqrt(pmax(fits[, "phi"], 0))), 4),
    edger_common_phi = round(edger_phi, 4),                 # edgeR cohort dispersion (cross-check)
    verdict = ifelse(frac_fdr >= 0.5, "Negative Binomial (overdispersed)", "Poisson-like"),
    stringsAsFactors = FALSE)
}

# ==================== run: every cohort (stages + healthy) x every cancer ====================
dist_rows <- list()
for (cname in names(DIST_CANCERS)) {
  pid <- unname(DIST_CANCERS[cname]); cat(sprintf("\n>>> %s (%s)\n", cname, pid))
  for (stage in c(DIST_STAGES, if (DIST_HEALTHY) "healthy")) {
    m_cnt <- load_counts_matrix(pid, stage)
    if (is.null(m_cnt) || ncol(m_cnt) < 3L) { cat(sprintf("    %-10s no/too-few counts; skip\n", stage)); next }
    dt <- tryCatch(test_count_distribution(m_cnt, cname, stage), error = function(e) NULL)
    if (!is.null(dt)) {
      dist_rows <- c(dist_rows, list(dt))
      cat(sprintf("    %-10s  %-33s reject-Poisson(FDR<.05)=%.0f%% of %d genes | median phi=%.3f | edgeR phi=%.3f\n",
                  stage, dt$verdict, 100 * dt$frac_reject_fdr, dt$n_genes_tested,
                  dt$median_phi, dt$edger_common_phi))
    } else {
      cat(sprintf("    %-10s  test not computable (too few testable genes)\n", stage))
    }
  }
}

# ==================== summary ====================
if (!length(dist_rows)) stop("No cohorts produced a result -- check data / config.")
dist_summary <- as.data.frame(data.table::rbindlist(dist_rows))
write.csv(dist_summary, file.path(OUT_DIR, "count_distribution_poisson_vs_nb.csv"), row.names = FALSE)

cat("\n", paste(rep("=", 116), collapse = ""), "\n", sep = "")
cat("COUNT DISTRIBUTION -- is the raw TCGA count data Poisson or Negative Binomial?\n")
cat("Per gene: intercept-only Poisson vs NB GLM with a log(library-size) offset; compared by LRT, AIC and dispersion.\n")
cat("frac_reject_pois / _fdr = fraction of genes significantly OVERDISPERSED (=> NB) at raw p<0.05 / BH-FDR<0.05.\n")
cat("frac_nb_by_aic = fraction where NB beats Poisson by AIC.  median_phi = NB dispersion (0 = Poisson; Var = mu + phi*mu^2).\n")
cat("median_bcv = sqrt(phi) = biological coeff. of variation.  edger_common_phi = edgeR cohort dispersion (independent check).\n")
cat("RNA-seq expectation: strong overdispersion => Negative Binomial (this is why edgeR/DESeq2 use NB, not Poisson).\n")
cat(paste(rep("=", 116), collapse = ""), "\n", sep = "")
print(dist_summary[order(dist_summary$cancer, dist_summary$stage), ], row.names = FALSE)
cat("\n--- overall: cohorts called NB vs Poisson ---\n")
cat(sprintf("  NB (overdispersed): %d / %d cohorts   |   median BCV across cohorts: %.3f\n",
            sum(grepl("^Negative", dist_summary$verdict)), nrow(dist_summary),
            median(dist_summary$median_bcv, na.rm = TRUE)))
cat("\nWrote:", file.path(OUT_DIR, "count_distribution_poisson_vs_nb.csv"), "\n")
