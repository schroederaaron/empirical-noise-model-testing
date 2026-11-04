#!/usr/bin/env Rscript
# make_yeast_cohorts.R
# -----------------------------------------------------------------------------
# Yeast 48x48 (Schurch et al. 2016, bartongroup/profDGE48) as base cohorts for the
# power benchmark. Reads the dataset Distribution_GOF/scripts/make_dataset.R built
# on the server (list(counts = genes x samples, samples = data.frame(sample_id,
# condition = WT / Snf2, ...), lengths = NULL, ...)) and writes ONE cohort RDS per
# condition in the format power_test.R's load_real_cohort() reads (as
# TCGA_test/scripts/export_power_cohort.R does for TCGA):
#
#   list(counts = genes x samples (integer), lengths = named per-gene vector,
#        project_id = "YEAST48", stage = "WT" | "Snf2", length_source, source_dataset)
#
# Parts R / RS then work on yeast unchanged: random halves (R) or 2n random samples
# (RS) of ONE condition, with known log2 FCs injected by binomial thinning.
#
# Gene lengths (TPM = counts / length, rescaled): the htseq counts carry none.
#   --GTF=<gtf>          union exon length per gene_id (the GTF the counts were made with)
#   --LENGTHS=<tsv>      two columns, gene and length, no header
#   neither              every length = 1, i.e. TPM = CPM (length_source = "uniform")
# Genes without a length in the GTF / table are dropped (reported).
#
#   Rscript Yeast_test/scripts/make_yeast_cohorts.R \
#     --DATASET=/media/BioNAS2/TCGA_TOX_TEST/data/yeast48/yeast48_clean.rds \
#     --OUT_DIR=/media/BioNAS2/TCGA_TOX_TEST/data/yeast48/cohorts [--GTF=...]
# -----------------------------------------------------------------------------

arg <- function(name, default = NULL) {
  a <- grep(paste0("^--", name, "="), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else default
}
DATASET <- arg("DATASET", "/media/BioNAS2/TCGA_TOX_TEST/data/yeast48/yeast48_clean.rds")
OUT_DIR <- arg("OUT_DIR", file.path(dirname(DATASET), "cohorts"))
GTF     <- arg("GTF"); LENGTHS <- arg("LENGTHS")

#' Union exon length per gene_id from a GTF (base R).
gtf_gene_lengths <- function(gtf) {
  x <- utils::read.delim(gtf, header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE)
  x <- x[x[[3]] == "exon", ]
  gid <- sub('.*gene_id "([^"]+)".*', "\\1", x[[9]])
  sp <- split(data.frame(s = x[[4]], e = x[[5]]), gid)
  vapply(sp, function(d) {
    d <- d[order(d$s), ]; tot <- 0; cs <- d$s[1]; ce <- d$e[1]
    for (i in seq_len(nrow(d))[-1]) {
      if (d$s[i] <= ce + 1) ce <- max(ce, d$e[i]) else { tot <- tot + ce - cs + 1; cs <- d$s[i]; ce <- d$e[i] }
    }
    tot + ce - cs + 1
  }, numeric(1))
}

ds <- readRDS(DATASET)
cnt <- as.matrix(ds$counts)
stopifnot(nrow(cnt) > ncol(cnt), !is.null(rownames(cnt)), identical(colnames(cnt), ds$samples$sample_id),
          all(cnt == round(cnt)), all(cnt >= 0))
storage.mode(cnt) <- "integer"

if (!is.null(GTF)) {
  L <- gtf_gene_lengths(GTF); src <- paste0("gtf:", normalizePath(GTF))
} else if (!is.null(LENGTHS)) {
  t <- utils::read.delim(LENGTHS, header = FALSE, stringsAsFactors = FALSE)
  L <- setNames(as.numeric(t[[2]]), t[[1]]); src <- paste0("table:", normalizePath(LENGTHS))
} else {
  L <- setNames(rep(1, nrow(cnt)), rownames(cnt)); src <- "uniform (TPM = CPM)"
}
miss <- setdiff(rownames(cnt), names(L)[is.finite(L) & L > 0])
if (length(miss)) message(sprintf("%d of %d genes have no length and are dropped (e.g. %s)",
                                  length(miss), nrow(cnt), paste(head(miss, 3), collapse = ", ")))
cnt <- cnt[!rownames(cnt) %in% miss, , drop = FALSE]

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
for (cond in levels(factor(ds$samples$condition))) {
  keep <- ds$samples$condition == cond
  co <- list(counts = cnt[, keep, drop = FALSE], lengths = L[rownames(cnt)],
             project_id = "YEAST48", stage = cond, length_source = src,
             source_dataset = normalizePath(DATASET), date = format(Sys.time()))
  f <- file.path(OUT_DIR, sprintf("yeast48_%s.rds", tolower(cond)))
  saveRDS(co, f)
  message(sprintf("%s: %d genes x %d samples (%s), lengths: %s -> %s",
                  cond, nrow(co$counts), ncol(co$counts), if (is.null(ds$label)) basename(DATASET) else ds$label, src, f))
}
