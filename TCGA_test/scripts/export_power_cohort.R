#!/usr/bin/env Rscript
# export_power_cohort.R
# -----------------------------------------------------------------------------
# Writes the input for power_test.R Part R (real-data signal injection): ONE
# homogeneous TCGA cohort as
#     list(counts = genes x samples (raw counts), tpm = samples x genes (raw TPM))
# over the samples and genes the two loaders have in common. power_test.R orients
# the TPM itself and recovers gene lengths from count/TPM.
#
#   Rscript export_power_cohort.R TCGA-KIRC "Stage IV" kirc_stage4.rds
#   Rscript export_power_cohort.R TCGA-KIRC healthy    kirc_healthy.rds
#   Rscript export_power_cohort.R all <out_dir>   # the six base cohorts of the plan:
#                                                 # COAD / LUAD / KIRC x healthy / Stage IV
#   POWER_REAL_RDS=a.rds,b.rds,... POWER_PARTS=R Rscript power_test.R
#
# Run from the same place as null_calibration.R (resolves outlier_significance_analysis.R,
# config.R, utils.R). Gene filtering is left to power_test.R (mean count >= 10).
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) == 3L || (length(args) == 2L && args[1] == "all")))
  stop('usage: Rscript export_power_cohort.R <project_id> <stage|healthy> <out.rds>\n',
       '       Rscript export_power_cohort.R all <out_dir>')

# ---- locate common/ (config.R, utils.R, ...) from THIS script's own location ----
# Slurm runs the scripts from the Tensor-Omics root, where a plain source("config.R")
# would pick up whatever copy sits in the working directory (e.g. a stale flat one).
# So resolve common/ relative to the script file (<script>/../../common); the other
# entries are fallbacks for interactive use, with a flat copy in "." last.
if (!exists("COMMON_DIR")) COMMON_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cand <- c(if (length(f)) file.path(dirname(normalizePath(f[1])), "..", "..", "common"),
            "analysis/common", "experiments/Noise_Model_Test/common", "common", ".")
  hit <- cand[file.exists(file.path(cand, "config.R"))]
  if (!length(hit)) stop("common/ not found (need config.R); looked in: ", paste(cand, collapse = ", "))
  normalizePath(hit[1])
})

source(file.path(COMMON_DIR, "outlier_significance_analysis.R"))   # load_stage_data, config, BASE_DATA_DIR

export_one <- function(pid, stage, out) {
  # Same file convention as load_counts_matrix() in null_calibration.R. Kept as a copy
  # because null_calibration.R runs its sweep when sourced.
  cnt_file <- if (identical(stage, "healthy"))
    file.path(BASE_DATA_DIR, pid, paste0("healthy_", pid, "_counts.rds")) else
    file.path(BASE_DATA_DIR, pid, paste0(pid, "-", gsub(" ", "-", stage), "_counts.rds"))
  if (!file.exists(cnt_file)) stop("counts file not found: ", cnt_file)
  co  <- readRDS(cnt_file)
  cnt <- t(co$expression_vectors)                          # samples x genes -> genes x samples
  colnames(cnt) <- rownames(co$expression_vectors)
  rownames(cnt) <- if (!is.null(co$gene_ids)) co$gene_ids else colnames(co$expression_vectors)

  d <- if (identical(stage, "healthy"))
    load_stage_data(pid, STAGES[1], "healthy", use_constant_healthy = TRUE,
                    norm_method = "raw", apply_mean = FALSE, normalize = FALSE) else
    load_stage_data(pid, stage, "cancer", use_constant_healthy = FALSE,
                    norm_method = "raw", apply_mean = FALSE, normalize = FALSE)
  if (is.null(d)) stop("TPM not found for ", pid, " / ", stage)
  tpm <- d$expression_vectors                              # samples x genes
  if (is.null(colnames(tpm)) && !is.null(d$gene_ids)) colnames(tpm) <- d$gene_ids

  samp <- intersect(rownames(tpm), colnames(cnt))
  gene <- intersect(colnames(tpm), rownames(cnt))
  cat(sprintf("%s / %s: %d matched samples (TPM %d, counts %d), %d matched genes (TPM %d, counts %d)\n",
              pid, stage, length(samp), nrow(tpm), ncol(cnt), length(gene), ncol(tpm), nrow(cnt)))
  if (length(samp) < 6L || length(gene) < 1000L)
    stop("too few matched samples/genes -- check that both loaders use the same ID conventions")
  stopifnot(length(gene) > length(samp))                   # genes x samples sanity

  saveRDS(list(counts = cnt[gene, samp, drop = FALSE], tpm = tpm[samp, gene, drop = FALSE],
               project_id = pid, stage = stage), out)
  cat("wrote ", out, "\n", sep = "")
}

if (args[1] == "all") {
  dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
  for (pid in c("TCGA-COAD", "TCGA-LUAD", "TCGA-KIRC")) for (stage in c("healthy", "Stage IV"))
    tryCatch(export_one(pid, stage, file.path(args[2], sprintf("%s_%s.rds", tolower(sub("TCGA-", "", pid)),
                                                          tolower(gsub(" ", "", stage))))),
             error = function(e) message("FAILED ", pid, " / ", stage, ": ", conditionMessage(e)))
} else export_one(args[1], args[2], args[3])
