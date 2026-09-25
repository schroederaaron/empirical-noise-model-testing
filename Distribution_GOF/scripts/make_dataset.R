#!/usr/bin/env Rscript
# make_dataset.R -- build one dataset object (brief section 3) and save it as .rds.
#
#   TCGA cohort (a stage, or the matched normals):
#     Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=tcga --PROJECT=TCGA-KIRC --STAGE="Stage I"
#     Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=tcga --PROJECT=TCGA-KIRC --STAGE=healthy
#   Simulated NB (pipeline checks without real data):
#     Rscript Distribution_GOF/scripts/make_dataset.R --SOURCE=simulated --N_GENES=2000 --N_PER_GROUP=24
#
# Optional: --LABEL=<id> (default derived from the inputs), --OUT=<path.rds>
# (default <OUT_ROOT>/datasets/<label>.rds, i.e. under the git-ignored results/).
# The integer audit is printed, so a ROUNDING decision (D1) can be made before
# 02_run_gof.R is started. tximport / featureCounts input needs a sample sheet;
# use build_from_tximport() / build_from_featurecounts() in R/io.R directly.
# TCGA needs the repo's config.R (BASE_DATA_DIR) in the working dir, and the NAS mount.

GOF_ROOT <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(Sys.getenv("GOF_ROOT", "Distribution_GOF"))
})
source(file.path(GOF_ROOT, "R", "setup.R"))
arg <- function(name, default = "") { v <- Sys.getenv(paste0("GOF_", name), ""); if (nzchar(v)) v else default }

src <- arg("SOURCE")
ds <- switch(src,
  tcga = {
    pid <- arg("PROJECT"); stage <- arg("STAGE")
    if (!nzchar(pid) || !nzchar(stage)) stop("--SOURCE=tcga needs --PROJECT=TCGA-XXXX and --STAGE=\"Stage I\" (or healthy)")
    build_from_tcga(pid, stage, arg("LABEL", paste0(sub("^TCGA-", "", pid), "_", gsub("[^A-Za-z0-9]", "", stage))))
  },
  simulated = build_simulated(as.integer(arg("N_GENES", "2000")), as.integer(arg("N_PER_GROUP", "24")),
                              seed = GOF_CONFIG$BASE_SEED, label = arg("LABEL", "simNB")),
  stop("--SOURCE must be 'tcga' or 'simulated' (got '", src, "')"))

out <- arg("OUT", file.path(GOF_CONFIG$OUT_ROOT, "datasets", paste0(ds$label, ".rds")))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(ds, out)

au <- integer_audit(ds$counts)
cat(sprintf("dataset '%s' (%s): %d genes x %d samples, design %s\n", ds$label, ds$source,
            nrow(ds$counts), ncol(ds$counts), deparse(ds$design)))
print(au[section == "overall"], row.names = FALSE)
if (audit_has_fractional(au))
  cat("NOTE: non-integer counts present -- 02_run_gof.R stops under ROUNDING = 'error';",
      "choose --ROUNDING=stochastic or floor (decision D1).\n")
cat("saved:", out, "\n")
