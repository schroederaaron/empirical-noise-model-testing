# setup.R -- shared bootstrap for every script in Distribution_GOF/ (issue #188).
#
# The LIB_DIR / .missing_pkg() / load_or_install() block below is copied verbatim
# from TCGA_test/scripts/count_distribution.R (brief section 2). Two additions,
# both marked:
#   * GOF_LIB_DIR overrides the library path (default unchanged:
#     external/docker_r_libs, relative to the working directory -- i.e. run the
#     scripts from the directory that holds external/, as on the cluster);
#   * USE_BUNDLED_LIBUV=1, because `fs` (a testthat dependency) otherwise needs the
#     system libuv headers, which the arch-gfortran image does not have. With it,
#     fs builds a static libuv and has no runtime system dependency.
# Separately: glmmTMB -> lme4 -> nloptr needs CMake AT BUILD TIME to compile its
# bundled static nlopt. The image has no cmake; install it once in the container
# that performs the first install (see README, "Packages").

LIB_DIR <- normalizePath(Sys.getenv("GOF_LIB_DIR", "external/docker_r_libs"), mustWork = FALSE)  # (addition: env override)
if (!dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE) && !dir.exists(LIB_DIR))
  stop("Could not create package library ", LIB_DIR,
       " -- it must be on a WRITABLE, BIND-MOUNTED path or packages will not persist.")
.libPaths(c(LIB_DIR, .libPaths()))
Sys.setenv(USE_BUNDLED_LIBUV = "1")                                                         # (addition)

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

# Each fork runs single-threaded, so N_CORES forks do not oversubscribe (repo convention).
Sys.setenv(OMP_NUM_THREADS = "1")

load_or_install("BiocManager")
load_or_install("edgeR")        # DGEList / filterByExpr / normLibSizes / getOffset / estimateDisp
load_or_install("data.table")
load_or_install("parallel")
load_or_install("glmmTMB")      # nbinom2 / zinb / genpois engines; compois + Laplace-PLN only as G3 cross-checks
load_or_install("statmod")      # gauss.quad(): Gauss-Hermite / Gauss-Legendre nodes (an edgeR dependency)

# Bootstrap streams are derived per (model, replicate) with set.seed() under this kind.
RNGkind("L'Ecuyer-CMRG")

source(file.path(GOF_ROOT, "config_gof.R"))
for (.f in c("io", "preprocess", "cmp_core", "pln_core", "families", "fit_dataset", "pit",
             "gof_stats", "strata", "bootstrap", "heldout", "report"))
  source(file.path(GOF_ROOT, "R", paste0(.f, ".R")))
rm(.f)
