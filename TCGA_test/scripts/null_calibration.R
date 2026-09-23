#!/usr/bin/env Rscript
# null_calibration.R
# -----------------------------------------------------------------------------
# THE decisive calibration experiment, now for ALL methods. Take ONE homogeneous
# group of samples (all tumours of a cancer stage), randomly split it into two
# fake groups ("A" vs "B"). There is no true differential expression between two
# random halves of the same population, so this is the null hypothesis H0.
#
# We run the core TOX arms (raw, log), edgeR, limma-voom and DESeq2 on the SAME H0 splits,
# and sweep TOX's kNN neighbourhood size (K_GRID) to see how it affects calibration
# (expected: little for raw, more for log). Optional TOX variants --
# edgeR-gene-set (degfilt) arms -- are commented out in TOX_ARMS for a cheaper sweep;
# re-enable them at the single best kNN config once this identifies one.
#
# We ALSO vary the number of replicates per fake group: besides the full half
# split, we draw fixed-size subsamples (e.g. 10 / 20 / 40 per group) so the
# influence of replicate count on calibration is visible for every method. A
# subsample arm only runs when the cohort has enough samples (>= 2 * size); a
# small stage (e.g. small-cell lung Stage IV, ~13 samples) simply skips the
# larger arms. Each result row is tagged with its per-group replicate count.
#
# NOTE on the p-value floor: it now differs by arm, so compare `min_p` against the
# arm's own floor, not a single number. Per gene:
#   exact                        1 / (n_pool_case * n_pool_control + 1)
#   bootstrap, null = 0 or       1 / (N_BOOTSTRAP_DRAWS + 1)  (= 1e-4 at 10,000 draws)
#     null = 1 above the cap
#   bootstrap, null = 1 enumerated  1 / (W_case * W_control + 1),
#                                W = (n_pool / n_rep) * n_rep^n_rep
# e.g. k = 30 genes at n_rep = 3: 1.2e-4 (exact), 1.0e-4 (sampled), 1.5e-6 (enumerated).
# BH at q needs the smallest p to reach q/G, so the floor decides whether an empty
# hit list means "no signal" or "not resolvable".
#
#   Under H0 a well-calibrated method must give:
#     * raw p-values ~UNIFORM  => FPR@alpha ~= alpha  (edgeR / limma / DESeq2 / TOX)
#     * ~0 calls at its own significance threshold      (hits@FDR<0.05 ~= 0)
#
# What the comparison decides:
#   * ALL methods ~calibrated  -> the DATASET is clean; the real-analysis
#     precision gap is power / effect-size philosophy, not false positives.
#   * only TOX inflated        -> a TOX null bug.
#   * edgeR/limma/DESeq2 inflated (esp. at large n, Stage I) -> the DEG tools
#     over-call under H0 -> the "reference" set contains false positives and is
#     NOT ground truth (confirms the Stage-I hypothesis directly).
#
# TOX faithfulness: observed statistic on the residual scale (linear for raw;
# log2(x+1) Frechet-mean difference for log/std_log/full). DEG tools use raw
# counts + their native filters, exactly as in deg_comparison.R.
# -----------------------------------------------------------------------------

LIB_DIR <- normalizePath("external/docker_r_libs", mustWork = FALSE)
if (!dir.create(LIB_DIR, recursive = TRUE, showWarnings = FALSE) && !dir.exists(LIB_DIR))
  stop("Could not create package library ", LIB_DIR,
       " -- it must be on a WRITABLE, BIND-MOUNTED path or packages will not persist.")
.libPaths(c(LIB_DIR, .libPaths()))

cat(sprintf("R package library: %s  (exists=%s, libPaths[1]=%s)\n",
            LIB_DIR, dir.exists(LIB_DIR), .libPaths()[1]))

# Name the package that is ACTUALLY missing from a load error. `library()` names
# the package being loaded FIRST, e.g.
#   "package or namespace load failed for 'DESeq2': there is no package called 'locfit'"
# so grabbing the first quoted token yields 'DESeq2' and reinstalls it forever. Match
# the specific phrases that name the missing dependency instead. Returns NA when the
# failure is not a missing-package error (e.g. a shared-object load failure).
.missing_pkg <- function(msg) {
  q <- "[‘’'\"`]([[:alnum:]._]+)[‘’'\"`]"
  for (pat in c(paste0("there is no package called ", q),
                paste0("package ", q, " required by"))) {
    m <- regmatches(msg, regexec(pat, msg, perl = TRUE))[[1]]
    if (length(m) >= 2L) return(m[2L])
  }
  NA_character_
}

# Load a package, installing it AND any MISSING DEPENDENCIES. Even when the package
# itself is present, library() can fail because a dependency is absent (e.g.
# 'generics' under 'BiocGenerics'); we install whatever package the error names and
# retry, resolving a whole chain. BiocManager::install resolves BOTH CRAN and
# Bioconductor packages + their deps, so no per-package "bioconductor" flag is needed.
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
    # A distinct missing dependency -> install it and retry (resolves a chain). If we
    # can't identify one (or it's the package itself, which install won't fix because
    # it's already current), stop and surface the REAL error rather than spinning 8x.
    if (is.na(miss) || identical(miss, package_name))
      stop("Could not load '", package_name, "': ", conditionMessage(err))
    message("  load_or_install: '", package_name, "' -> installing missing dependency '", miss, "' ...")
    install_pkg(miss)
  }
  stop("Could not load '", package_name, "' -- unresolved dependencies after 8 attempts.")
}

load_or_install("BiocManager")
load_or_install("edgeR")
load_or_install("limma")
load_or_install("DESeq2")
load_or_install("ggplot2")
load_or_install("data.table")
load_or_install("dplyr")
load_or_install("goftest")               # Anderson-Darling / Cramer-von Mises GoF tests
load_or_install("parallel")              # mclapply -- fork-based parallel run sweep

# Each fork must run single-threaded: the TOX Fortran uses OpenMP `do concurrent`,
# so N_CORES forks x N_CORES OpenMP threads would oversubscribe (e.g. 32x32 = 1024
# threads) and run SLOWER than serial. The mclapply forks provide the parallelism.
Sys.setenv(OMP_NUM_THREADS = "1")

# Headless containers often ship with NO fonts, so the cairo PNG device renders
# every glyph as a "tofu" box. If none are present, install DejaVu + fontconfig so
# plot labels render. (Containers are ephemeral, so this runs each fresh container;
# bake it into the image -- e.g. `RUN pacman -Sy --noconfirm fontconfig ttf-dejavu
# && fc-cache -f` -- to avoid the per-run cost.)
ensure_fonts <- function() {
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
ensure_fonts()


source("outlier_significance_analysis.R")   # loaders, preprocess_replicates, tox wrappers
source("config.R")
suppressMessages({library(dplyr); library(ggplot2); library(data.table)})
options(width = 200)

OUT_DIR  <- file.path(dirname(TOX_TEST_DIR), "null_calibration")
PLOT_DIR <- file.path(OUT_DIR, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==================== CONFIG ====================

# Back to the focused 3-cancer subset: a mixture cancer (KIRC), a high-signal one (LUAD)
# and a heavy-tailed-control one (STAD). Swap in `CANCER_TYPES` for the full 8-cancer run.
CALIB_CANCERS <- CANCER_TYPES[c("kidney cancer",
                                "non small cell lung cancer",
                                "stomach cancer")]
CALIB_STAGES  <- STAGES                 # the 4 cancer stages; restrict for speed, e.g. c("Stage I")
RUN_HEALTHY   <- TRUE                   # ALSO run the matched-normal (healthy) cohort per cancer
                                        # (loaded via healthy_<pid>_counts.rds / the "healthy" TPM ref)

# TOX arms: each is a (norm, shared, label).
# For the kNN sweep we run only the CORE log + raw arms (no degfilt): the k effect is
# a norm question, so if k doesn't move log/raw calibration here it won't move the
# degfilt variants either -- and if it DOES, we re-enable those variants at the
# single best kNN config identified here (much cheaper than sweeping all arms x all k).
# The commented arms below (results already collected) can be switched back on then.
#   shared = TRUE runs on the edgeR/limma filterByExpr gene set (apples-to-apples).
# `model` selects the null construction: "exact" (sqrt-scaled individual-residual
# pairs) or "bootstrap" (differences of resampled MEANS). The bootstrap arms are
# the direct test of the shape diagnosis in raw_anticonservative_diagnosis.R: the
# raw model's anti-conservatism at alpha = 0.05 is attributed there to the exact
# null keeping the INDIVIDUAL-residual shape while the observed statistic is a
# difference of means, so a null built at the mean level should recover the 0.05
# level. It is NOT expected to fix alpha = 0.01, since resampling cannot add tail
# mass the pool does not contain -- if TOX-raw-boot fixes both, the diagnosis is
# incomplete; if it fixes neither, it is wrong.
# (The raw tail-trim arms that used to sit here are gone: trimming narrows an
# already-too-narrow core, it measurably made raw calibration worse, and the knob
# has been removed from the model entirely.)
# `null` selects the bootstrap module's null construction (ignored by "exact"):
#   0 = POOLED  -- resample n_rep residuals iid from the whole neighbourhood pool
#   1 = BLOCKED -- pick one neighbour GENE, resample within it, so each null mean
#                  carries a single coherent noise level; exactly ENUMERATED (no RNG,
#                  no 1/(n_boot+1) floor) wherever C(2n-1,n)*n_genes_pool fits the
#                  Fortran cap, which covers n_rep <= 5 at k_max = 50.
# The *-blocked arms are the ones to watch: the simulation work says the pooled
# bootstrap is anti-conservative at alpha = 0.01 (1.4-2.0x) while blocking restores
# it, and that the blocked null produces far fewer BH false positives under a
# complete null. Both claims are about REAL data here for the first time.
# THIS RUN: only raw + log, each in BOTH engines (exact and bootstrap), pooled and
# blocked nulls. No degfilt, no kNN sweep -- a clean engine x normalization x null
# comparison on the 3-cancer subset. The degfilt arms below are kept ready to
# re-enable (just uncomment).
TOX_ARMS <- list(
  list(norm = "log", shared = FALSE, model = "exact",     null = 0L, label = "TOX-log"),
  list(norm = "raw", shared = FALSE, model = "exact",     null = 0L, label = "TOX-raw"),
  list(norm = "raw", shared = FALSE, model = "bootstrap", null = 0L, label = "TOX-raw-boot"),
  list(norm = "log", shared = FALSE, model = "bootstrap", null = 0L, label = "TOX-log-boot"),
  list(norm = "raw", shared = FALSE, model = "bootstrap", null = 1L, label = "TOX-raw-blocked"),
  list(norm = "log", shared = FALSE, model = "bootstrap", null = 1L, label = "TOX-log-blocked")
  # , list(norm = "log", shared = TRUE,  model = "exact",     null = 0L, label = "TOX-log-degfilt")
  # , list(norm = "raw", shared = TRUE,  model = "exact",     null = 0L, label = "TOX-raw-degfilt")
)
# Only build the (costly) shared edgeR/limma gene set per cohort if a *-degfilt arm needs it.
ANY_SHARED_ARM <- any(vapply(TOX_ARMS, function(a) isTRUE(a$shared), logical(1)))

# ---------------- normalisation comparison: counts vs TPM ----------------
# Plan section 6 / validation item 12 of md_noise_model_plan_independent_axes.md, run
# here at d = 1 where it is cheap and directly comparable to everything else in this
# sweep. The question it settles is which input scale TOX should take, decided
# EMPIRICALLY rather than argued from the literature.
#
# Why it matters. TPM is compositional: log2FC = beta - Delta with Delta identical
# for every gene, and TOX's null -- built from WITHIN-group residuals -- contains no
# information about a BETWEEN-group offset. Running uncorrected is therefore not the
# assumption-free option, it is the assumption Delta = 0. Two ways to do better:
#   * correct on the TPM itself  -> `beta_centre = 1`, i.e. subtract the per-axis
#     MEDIAN of beta across genes (median centring of log-ratios). Assumes only that
#     the median gene is unchanged, which is strictly weaker than Delta = 0.
#   * normalise the COUNTS instead -> TMM (edgeR) or median-of-ratios (DESeq2), the
#     same estimator family but with weights that come from count magnitudes, which
#     do not exist on TPM.
# Naming, so the comparison is labelled honestly: what the `-centred` arm does is
# median centring of log-ratios, NOT TMM.
#
# A second, separate reason to expect counts to win shows up in POWER rather than
# calibration (so this null run will not see it): TPM removes the count magnitude,
# which is what predicts technical precision, and its length correction works against
# the kNN matching -- at equal TPM a longer gene carries more counts and is more
# precise, roughly a 4-5x CV spread across 500 bp - 10 kb among genes the matcher
# treats as equivalent.
#
# Kept deliberately small: 4 arms, ONE kNN config, the full half split only, and a
# 3-cancer x 3-stage subset. Everything except the input scale is held fixed --
# same samples, same genes, same split -- or the comparison confounds normalisation
# with cohort differences.
RUN_NORM_TEST <- TRUE
NORM_STAGES   <- c("healthy", "Stage I", "Stage IV")   # a low and a high stage, plus the normals
NORM_SIZES    <- c(0L)                                 # full half split only (0L); n-dependence
                                                       # is already answered by the main sweep
NORM_ARMS <- list(
  list(input = "tpm", centre = 0L, label = "NORM-TPM"),
  list(input = "tpm", centre = 1L, label = "NORM-TPM-centred"),
  list(input = "tmm", centre = 0L, label = "NORM-TMM"),
  list(input = "mor", centre = 0L, label = "NORM-MoR")
)
# All four are TOX-log in the `exact` engine with the pooled null -- only the input
# scale and `centre` vary. NORM-TPM is therefore the same computation as the TOX-log
# arm above (on the matched gene/sample subset), which makes it a built-in
# consistency check on this whole block.

# Which methods to run.
RUN_TOX    <- TRUE
RUN_EDGER  <- TRUE
RUN_LIMMA  <- TRUE
RUN_DESEQ2 <- TRUE                       # SLOW reference (median-of-ratios NB GLM + Wald)

# Independent random partitions ("runs"). The reported result is the MEAN over
# these runs, with the across-run SD so a single lucky/unlucky split is visible
# and averaged away; every run's per-cell metrics are ALSO saved individually to
# null_calibration_per_run.csv (not printed). WARNING: runtime scales linearly --
# 100 runs over all cells and methods is a long job.
N_SPLITS         <- 100L                 # runs for TOX / edgeR / limma
N_CORES          <- 32L                  # cores for the parallel (mclapply) run sweep
MIN_SAMPLES      <- 16L                  # need enough samples for the full half split
ALPHAS           <- c(0.05, 0.01)

# NOTE: the Poisson-vs-NB count-distribution diagnostic now lives in its own script,
# count_distribution.R (it's a property of the data, not the null-split experiment).

# Replicate-count sweep. Besides the full half split (recorded at its own
# per-group size), draw these FIXED per-group replicate counts so we can see how
# calibration depends on n. An arm needs 2 * size samples (size per fake group);
# arms that don't fit the cohort are skipped. MIN_PER_GROUP guards the full arm
# on small cohorts. Each output row carries `n_per_group` = the per-group count.
SUBSAMPLE_SIZES  <- c(10L, 20L, 40L)     # per-group replicate counts to probe
MIN_PER_GROUP    <- 3L                   # smallest usable group (full arm floor)

# TOX kNN neighbourhood -- FIXED at the best config the 9-config sweep identified.
# Evidence (docs/raw_normalisation_diagnosis.md §1.1, reproduce with
# reanalyse_null_calibration_sweep.R): paired within (cancer, stage, replicate arm),
# going from the smallest config to k20_50 improves inflation in 128/128 comparisons
# (TOX-raw -0.093, TOX-log -0.109), i.e. LARGER k is consistently better -- so k20_50,
# the largest tested, is the best of the grid. Note this also SETTLES the question:
# k is NOT the cause of the raw arm's anti-conservatism (it closes only ~1/7 of the
# 0.64 excess; k_max 30->50 moves it by 0.005) -- that is a null-SHAPE problem, which
# is what the blocked null addresses. So no more sweeping: one config, evidence-backed.
# (Uncomment the sweep block only if the grid itself ever needs revisiting.)
K_GRID <- list(list(k_start = 20L, k_step = 1L, k_max = 50L, name = "k20_50"))
# --- full sweep (uncomment to re-run; slow): ---
# K_GRID <- do.call(c, lapply(c(8L, 15L, 20L), function(ks)
#   lapply(c(30L, 40L, 50L), function(km)
#     list(k_start = ks, k_step = 1L, k_max = km, name = sprintf("k%d_%d", ks, km)))))
TAU <- 0.1; MAX_POOL <- 70000L
# Both null constructions, selected per arm by `model` (see TOX_ARMS above).
TOX_MODEL_FNS <- list(exact     = tox_compute_noise_pvalues_pipeline_exact,
                      bootstrap = tox_compute_noise_pvalues_pipeline)
TOX_MODEL_FN  <- TOX_MODEL_FNS$exact          # default for callers that pass no model

set.seed(42)

if (RUN_EDGER || RUN_LIMMA || RUN_NORM_TEST) suppressMessages({library(edgeR); library(limma)})
if (RUN_DESEQ2 || RUN_NORM_TEST)             suppressMessages(library(DESeq2))

# ==================== data loaders ====================

#' Load raw COUNTS (genes x samples) for a cancer stage OR the matched-normal cohort,
#' mirroring deg_comparison.R. `stage == "healthy"` reads healthy_<pid>_counts.rds (the
#' matched normals, one group per cancer); any other value reads <pid>-<stage>_counts.rds.
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

# ==================== metric helpers ====================

#' Per-run metrics for a raw p-value vector: false-positive rates plus three
#' goodness-of-fit DISTANCES to Uniform(0,1) -- ks_D (Kolmogorov-Smirnov), ad_A2
#' (Anderson-Darling), cvm_W2 (Cramer-von Mises). These statistics (0 = perfectly
#' uniform) are the effect sizes to read; their p-values collapse to ~0 at this n
#' and are not aggregated. p is clamped off {0,1} for AD/CvM (AD -> Inf at the
#' boundary); the fpr/median use the raw p.
pval_metrics <- function(p) {
  p <- p[is.finite(p)]
  n <- length(p)
  if (n < 10)
    return(list(n = n, fpr_0.05 = NA_real_, fpr_0.01 = NA_real_, median_p = NA_real_,
                hits_fdr05 = NA_real_, ks_D = NA_real_, ad_A2 = NA_real_, cvm_W2 = NA_real_))
  pc  <- pmin(pmax(p, 1e-9), 1 - 1e-9)
  ks  <- suppressWarnings(ks.test(p, "punif"))
  ad  <- suppressWarnings(tryCatch(goftest::ad.test(pc, "punif")$statistic,  error = function(e) NA_real_))
  cvm <- suppressWarnings(tryCatch(goftest::cvm.test(pc, "punif")$statistic, error = function(e) NA_real_))
  list(n = n, fpr_0.05 = mean(p < 0.05), fpr_0.01 = mean(p < 0.01),
       median_p = median(p), hits_fdr05 = mean(p.adjust(p, "BH") < 0.05),
       ks_D = unname(ks$statistic), ad_A2 = unname(ad), cvm_W2 = unname(cvm))
}
# Build one A/B partition of `n` samples at a target per-group replicate count.
# `size == 0L` means the FULL half split (g = n %/% 2); any other size draws
# exactly `size` samples per group (2 * size total). Returns NULL when the
# cohort is too small for this arm, so the caller simply skips it. `g` is the
# realised per-group count, recorded as n_per_group in the output.
pick_split <- function(n, size) {
  g <- if (size == 0L) n %/% 2L else size
  if (g < MIN_PER_GROUP || n < 2L * g) return(NULL)
  ix <- sample.int(n, 2L * g)
  list(a = ix[seq_len(g)], b = ix[g + seq_len(g)], g = g)
}

# ==================== TOX ====================

#' One TOX run on one A/B split.
#'
#' `beta_centre = 1L` subtracts the per-gene MEDIAN of the observed statistic across
#' genes before scoring -- the composition correction of plan section 6. At d = 1
#' that is exactly what `noise_model_md`'s `beta_centre` does (its `delta_hat` is the
#' median of beta over genes), so it is done here in R rather than routed through the
#' multidimensional entry point: the scalar `exact` module computes the same p-value
#' in O((n+m) log m) via its sorted-pool binary search, while the d-dimensional module
#' would have to enumerate or sample the n_pool_case * n_pool_ctrl pairs explicitly.
#' Same answer, and much faster at the pool sizes this sweep uses.
run_tox_once <- function(raw_mat, norm_method, sp, kcfg = K_GRID[[1]],
                         model = "exact", null_method = 0L, beta_centre = 0L) {
  ng <- ncol(raw_mat)
  pa <- preprocess_replicates(raw_mat[sp$a, , drop = FALSE], norm_method)
  pb <- preprocess_replicates(raw_mat[sp$b, , drop = FALSE], norm_method)
  norm_int <- if (norm_method == "raw") 0L else 1L
  if (norm_int == 0L) {
    obs_own <- colMeans(pa$prelog, na.rm = TRUE) - colMeans(pb$prelog, na.rm = TRUE)
  } else {
    obs_own <- colMeans(log2(pa$prelog + 1), na.rm = TRUE) -
               colMeans(log2(pb$prelog + 1), na.rm = TRUE)
  }
  obs_own <- as.numeric(obs_own); valid <- as.integer(is.finite(obs_own))
  # Composition centring, before the non-finite genes are zeroed (they must not
  # contribute to the median, and `valid` already excludes them from the test).
  if (isTRUE(beta_centre == 1L)) {
    delta_hat <- median(obs_own[is.finite(obs_own)])
    if (is.finite(delta_hat)) obs_own <- obs_own - delta_hat
  }
  obs_own[!is.finite(obs_own)] <- 0
  model_fn <- TOX_MODEL_FNS[[model]]
  if (is.null(model_fn)) stop("run_tox_once: unknown model '", model, "'")
  res <- model_fn(
    case_means = as.numeric(pa$means), case_replicates = pa$prelog,
    control_means = as.numeric(pb$means), control_replicates = pb$prelog,
    obs_own = obs_own, valid_genes_own = valid,
    norm_method = norm_int, k_start = kcfg$k_start, k_step = kcfg$k_step, k_max = kcfg$k_max,
    tau = TAU, null_method = null_method, max_pool_size = MAX_POOL)
  p <- res$pvalues_own; p[p < 0 | p > 1] <- NA; p[!is.na(p)]
}

#' Build the three input scales for the normalisation comparison, matched.
#'
#' Returns `list(tpm =, tmm =, mor =)`, each samples x genes over the SAME samples
#' and the SAME genes, or NULL when the cohort cannot supply them. Matching is the
#' whole point: `m_tpm` and `m_cnt` arrive with different sample sets (different
#' loaders) and different gene sets (TOX's keep mask vs the raw count rows), and the
#' main sweep gives them independent splits. Comparing normalisations across those
#' would confound the input scale with the cohort and the gene set.
#'
#' Normalisation factors are estimated on the FULL cohort, before splitting -- that
#' is how they would be used in practice, and since both fake groups are drawn from
#' one cohort it cannot manufacture a difference between them.
#'
#'   tmm: edgeR TMM factors -> CPM (`normLibSizes` + `cpm`), linear scale.
#'   mor: DESeq2 median-of-ratios size factors -> counts divided per sample.
#'
#' Both are left on a LINEAR scale: TOX applies log2(x + 1) internally for
#' `norm_method = "log"`.
build_norm_inputs <- function(m_tpm, m_cnt) {
  if (is.null(m_tpm) || is.null(m_cnt)) return(NULL)
  samp <- intersect(rownames(m_tpm), colnames(m_cnt))
  gene <- intersect(colnames(m_tpm), rownames(m_cnt))
  if (length(samp) < 2L * MIN_PER_GROUP || length(gene) < 100L) {
    # Say WHY. The TPM and count matrices come from different loaders, so an empty
    # intersection almost always means the two use different sample- or gene-ID
    # conventions, not that the cohort is too small -- and those look identical
    # from a bare "skipped" line.
    message(sprintf("  build_norm_inputs: only %d matched samples (of %d TPM / %d count) and %d matched genes (of %d / %d)",
                    length(samp), nrow(m_tpm), ncol(m_cnt),
                    length(gene), ncol(m_tpm), nrow(m_cnt)))
    return(NULL)
  }

  cnt <- m_cnt[gene, samp, drop = FALSE]
  tpm <- m_tpm[samp, gene, drop = FALSE]

  # Each scale fails independently: a broken MoR must not silently take the TPM and
  # TMM arms down with it, or the block would look like it simply did not run.
  tmm <- tryCatch({
    dge <- edgeR::normLibSizes(edgeR::DGEList(counts = cnt))
    t(edgeR::cpm(dge, normalized.lib.sizes = TRUE))
  }, error = function(e) { message("  build_norm_inputs: TMM failed -- ", conditionMessage(e)); NULL })

  mor <- tryCatch({
    sf <- DESeq2::estimateSizeFactorsForMatrix(cnt)
    if (!all(is.finite(sf)) || any(sf <= 0)) stop("non-finite or non-positive size factors")
    t(sweep(cnt, 2L, sf, "/"))          # cnt is genes x samples -> divide per SAMPLE
  }, error = function(e) { message("  build_norm_inputs: MoR failed -- ", conditionMessage(e)); NULL })

  list(tpm = tpm, tmm = tmm, mor = mor)
}

# ==================== reference methods (native counts) ====================

# edgeR 4.x quasi-likelihood pipeline (current recommended standard). With the
# default legacy = FALSE, glmQLFit() estimates the NB/QL dispersion internally, so
# estimateDisp() is no longer required (results are identical with or without it in
# edgeR >= 4.0). normLibSizes() is the current name for calcNormFactors() (TMM), and
# robust = TRUE guards the dispersion fit against hypervariable genes. coef = 2 of
# ~group is the group effect, i.e. the treatment-vs-control contrast.
ref_edgeR <- function(counts, group) {
  design <- model.matrix(~ group)
  y <- DGEList(counts = counts, group = group)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  fit <- glmQLFit(y, design, robust = TRUE)
  glmQLFTest(fit, coef = 2)$table$PValue
}
# limma-voom (current recommended standard): voomLmFit() is the modern combined
# voom + lmFit entry point (it also accounts for residual df lost to zero counts, as
# glmQLFit does); robust = TRUE moderates against variance outliers. trend = FALSE is
# the eBayes default and correct here (voom already models the mean-variance trend).
ref_limma <- function(counts, group) {
  design <- model.matrix(~ group)
  dge <- DGEList(counts = counts, group = group)
  dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
  dge <- normLibSizes(dge)
  fit <- eBayes(voomLmFit(dge, design), robust = TRUE)
  topTable(fit, coef = 2, number = Inf, sort.by = "none")$P.Value
}

# DESeq2 (current standard): median-of-ratios size factors + NB GLM with empirical-Bayes
# dispersion shrinkage; Wald test on the group coefficient. `results()` defaults to the
# last coefficient (group B-vs-A here) and applies BH + independent filtering -- the latter
# only NAs out `padj` for low-count genes; the raw `pvalue` we score for uniformity is still
# returned (NAs are dropped downstream by pval_metrics). `quiet = TRUE` + suppressMessages
# keep DESeq2's progress chatter out of the log. Note: DESeq2 is the SLOW reference.
ref_deseq2 <- function(counts, group) {
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData   = data.frame(group = group),
                                design    = ~ group)
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  as.numeric(results(dds)$pvalue)
}

# ==================== sweep ====================

# Pure builders (return data.frames) so results can be collected OUT of parallel
# workers -- mclapply children cannot write to a parent-global accumulator.
make_row <- function(method, cancer, stage, run, n_samples, n_per_group, k_config, mt) {
  data.frame(method = method, cancer = cancer, stage = stage, run = run,
             n_samples = n_samples, n_per_group = n_per_group, k_config = k_config, n_genes = mt$n,
             fpr_0.05 = mt$fpr_0.05, fpr_0.01 = mt$fpr_0.01,
             median_p = mt$median_p, hits_fdr05 = mt$hits_fdr05,
             ks_D = mt$ks_D, ad_A2 = mt$ad_A2, cvm_W2 = mt$cvm_W2,
             stringsAsFactors = FALSE)
}
# A capped random sample of raw p-values, kept for the histograms and QQ-plots.
make_pval <- function(method, cancer, stage, n_per_group, k_config, pv) {
  pv <- pv[is.finite(pv)]
  if (!length(pv)) return(NULL)
  data.frame(method = method, cancer = cancer, stage = stage, n_per_group = n_per_group,
             k_config = k_config, p = sample(pv, min(500L, length(pv))), stringsAsFactors = FALSE)
}

# One run = one set of random partitions, one per replicate-count arm. For each
# arm (full half + every feasible fixed size) it draws an independent A/B split
# and runs every method at that per-group replicate count, tagging each row with
# n_per_group. Pure (no globals) so it is safe inside an mclapply worker.
run_one_split <- function(s, m_tpm, m_cnt, cname, stage, n_tox, n_cnt, m_tpm_shared = NULL,
                          m_norm = NULL) {
  rws <- list(); pvs <- list()
  arms <- c(0L, SUBSAMPLE_SIZES)      # 0L = full half split; rest = fixed per-group sizes

  # ---- normalisation comparison: counts vs TPM, everything else held fixed ----
  # One split per size, SHARED by all four arms, so the only thing that differs
  # between them is the input scale and `centre`.
  if (!is.null(m_norm)) for (size in NORM_SIZES) {
    spn <- pick_split(nrow(m_norm$tpm), size)
    if (is.null(spn)) next
    kcfg <- K_GRID[[1]]
    for (arm in NORM_ARMS) {
      mat <- m_norm[[arm$input]]
      if (is.null(mat)) next
      pv <- tryCatch(run_tox_once(mat, "log", spn, kcfg, "exact", 0L, arm$centre),
                     error = function(e) numeric(0))
      if (length(pv)) {
        rws <- c(rws, list(make_row(arm$label, cname, stage, s, nrow(mat), spn$g,
                                    kcfg$name, pval_metrics(pv))))
        pvs <- c(pvs, list(make_pval(arm$label, cname, stage, spn$g, kcfg$name, pv)))
      }
    }
  }

  for (size in arms) {
    if (!is.null(m_tpm)) {
      sp <- pick_split(n_tox, size)   # split is over SAMPLES; the same sp serves every TOX arm
      # Each TOX arm runs at EVERY kNN config (K_GRID); the config name is tagged as k_config.
      if (!is.null(sp)) for (arm in TOX_ARMS) {
        # `shared` arms run on the edgeR/limma gene set (`m_tpm_shared`: same samples/rows as
        # m_tpm, columns = the edgeR keep genes). Skip them when there is no shared set.
        if (isTRUE(arm$shared) && (is.null(m_tpm_shared) || !ncol(m_tpm_shared))) next
        mt_arm <- if (isTRUE(arm$shared)) m_tpm_shared else m_tpm
        for (kcfg in K_GRID) {
          pv <- tryCatch(run_tox_once(mt_arm, arm$norm, sp, kcfg,
                                      if (is.null(arm$model)) "exact" else arm$model,
                                      if (is.null(arm$null)) 0L else arm$null),
                         error = function(e) numeric(0))
          if (length(pv)) {
            rws <- c(rws, list(make_row(arm$label, cname, stage, s, n_tox, sp$g, kcfg$name, pval_metrics(pv))))
            pvs <- c(pvs, list(make_pval(arm$label, cname, stage, sp$g, kcfg$name, pv)))
          }
        }
      }
    }
    if (!is.null(m_cnt)) {
      sp <- pick_split(n_cnt, size)
      if (!is.null(sp)) {
        grp <- factor(c(rep("A", length(sp$a)), rep("B", length(sp$b))))
        cnt <- m_cnt[, c(sp$a, sp$b), drop = FALSE]
        # edgeR/limma/DESeq2 are kNN-independent -> k_config = NA.
        if (RUN_EDGER) { pv <- tryCatch(ref_edgeR(cnt, grp), error = function(e) numeric(0))
          if (length(pv)) { rws <- c(rws, list(make_row("edgeR", cname, stage, s, n_cnt, sp$g, NA_character_, pval_metrics(pv))))
                            pvs <- c(pvs, list(make_pval("edgeR", cname, stage, sp$g, NA_character_, pv))) } }
        if (RUN_LIMMA) { pv <- tryCatch(ref_limma(cnt, grp), error = function(e) numeric(0))
          if (length(pv)) { rws <- c(rws, list(make_row("limma", cname, stage, s, n_cnt, sp$g, NA_character_, pval_metrics(pv))))
                            pvs <- c(pvs, list(make_pval("limma", cname, stage, sp$g, NA_character_, pv))) } }
        if (RUN_DESEQ2) { pv <- tryCatch(ref_deseq2(cnt, grp), error = function(e) numeric(0))
          if (length(pv)) { rws <- c(rws, list(make_row("DESeq2", cname, stage, s, n_cnt, sp$g, NA_character_, pval_metrics(pv))))
                            pvs <- c(pvs, list(make_pval("DESeq2", cname, stage, sp$g, NA_character_, pv))) } }
      }
    }
  }
  list(rows = rws, pvals = pvs)
}

RNGkind("L'Ecuyer-CMRG")   # parallel-safe RNG streams so the mclapply runs are reproducible
set.seed(42)

rows <- list(); pval_pool <- list()

for (cname in names(CALIB_CANCERS)) {
  pid <- unname(CALIB_CANCERS[cname]); cat(sprintf("\n>>> %s (%s)\n", cname, pid))
  keep_mask <- tryCatch(compute_gene_keep_mask(pid), error = function(e) NULL)
  kept_ids  <- if (!is.null(keep_mask)) names(keep_mask)[keep_mask] else NULL

  # Cohorts: the 4 cancer stages plus (optionally) the matched-normal "healthy" group.
  for (stage in c(CALIB_STAGES, if (RUN_HEALTHY) "healthy")) {
    # ---- TOX input: raw TPM (samples x genes). m_tpm_full keeps ALL genes; m_tpm is
    #      restricted to TOX's own compute_gene_keep_mask filter. "healthy" loads the
    #      constant matched-normal reference; a stage loads that stage's cancer TPM. ----
    m_tpm <- NULL; m_tpm_full <- NULL
    if (RUN_TOX && !is.null(kept_ids)) {
      d <- tryCatch(
        if (identical(stage, "healthy"))
          load_stage_data(pid, STAGES[1], "healthy", use_constant_healthy = TRUE,
                          norm_method = "raw", apply_mean = FALSE, normalize = FALSE)
        else
          load_stage_data(pid, stage, "cancer", use_constant_healthy = FALSE,
                          norm_method = "raw", apply_mean = FALSE, normalize = FALSE),
        error = function(e) NULL)
      if (!is.null(d) && is.matrix(d$expression_vectors)) {
        m <- d$expression_vectors
        if (is.null(colnames(m)) && !is.null(d$gene_ids)) colnames(m) <- d$gene_ids
        m_tpm_full <- m
        m_tpm <- m[, intersect(kept_ids, colnames(m)), drop = FALSE]
      }
    }
    # ---- reference input: raw counts (genes x samples), native filters ----
    m_cnt <- if (RUN_EDGER || RUN_LIMMA || RUN_DESEQ2 || RUN_NORM_TEST)
               load_counts_matrix(pid, stage) else NULL

    n_tox <- if (!is.null(m_tpm)) nrow(m_tpm) else 0L
    n_cnt <- if (!is.null(m_cnt)) ncol(m_cnt) else 0L
    if (max(n_tox, n_cnt) < MIN_SAMPLES) { cat(sprintf("    %-10s too few samples; skip\n", stage)); next }

    # ---- shared gene set: the edgeR/limma keep mask (cohort-level filterByExpr on the
    #      raw counts) intersected with the TPM genes. TOX's `*-degfilt` arms run on this,
    #      so TOX and the DE tools operate on the EXACT same genes. Same samples/rows as
    #      m_tpm, so a sample split applies to it unchanged. ----
    m_tpm_shared <- NULL
    if (ANY_SHARED_ARM && !is.null(m_tpm_full) && !is.null(m_cnt)) {   # skip the cost if no *-degfilt arm is active
      deg_keep <- tryCatch(suppressWarnings(rownames(m_cnt)[filterByExpr(DGEList(counts = m_cnt))]),
                           error = function(e) NULL)
      shared_ids <- if (!is.null(deg_keep)) intersect(colnames(m_tpm_full), deg_keep) else character(0)
      if (length(shared_ids) >= 10L) m_tpm_shared <- m_tpm_full[, shared_ids, drop = FALSE]
    }

    # ---- matched input scales for the normalisation comparison (TPM / TMM / MoR).
    #      Built once per cohort, on the full sample set, and only for the stages the
    #      comparison actually runs on. ----
    m_norm <- NULL
    if (RUN_NORM_TEST && stage %in% NORM_STAGES) {
      m_norm <- tryCatch(build_norm_inputs(m_tpm, m_cnt), error = function(e) NULL)
      if (is.null(m_norm))
        cat(sprintf("    %-10s  norm-test inputs unavailable (no matched samples/genes); skipped\n",
                    stage))
    }

    # ---- N_SPLITS runs IN PARALLEL over N_CORES. mc.silent = TRUE swallows ALL
    #      child stdout, which also silences the noisy per-call debug prints that
    #      used to flood the log. ----
    res_list <- mclapply(seq_len(N_SPLITS), run_one_split,
                         m_tpm = m_tpm, m_cnt = m_cnt, cname = cname, stage = stage,
                         n_tox = n_tox, n_cnt = n_cnt, m_tpm_shared = m_tpm_shared,
                         m_norm = m_norm,
                         mc.cores = N_CORES, mc.preschedule = TRUE, mc.silent = TRUE)
    ok <- vapply(res_list, is.list, logical(1))        # drop any worker that errored out
    for (r in res_list[ok]) { rows <- c(rows, r$rows); pval_pool <- c(pval_pool, r$pvals) }

    cat(sprintf("    %-10s  n(tox)=%d n(cnt)=%d%s  [%d/%d runs OK]\n",
                stage, n_tox, n_cnt,
                if (is.null(m_norm)) "" else sprintf(" norm(n=%d g=%d)",
                                                     nrow(m_norm$tpm), ncol(m_norm$tpm)),
                sum(ok), N_SPLITS))
  }
}

if (length(rows) == 0) stop("No null runs produced -- check data / config.")
detail <- as.data.frame(rbindlist(rows))
# Individual per-run results (all N_SPLITS partitions per cell) -- saved, not printed.
write.csv(detail, file.path(OUT_DIR, "null_calibration_per_run.csv"), row.names = FALSE)

# ==================== summary + verdict ====================

# Aggregate per (method, cancer, stage, n_per_group, k_config): each replicate-count
# arm and each kNN config is its own row so the influence of n AND of the neighbourhood
# size on calibration is directly readable. (edgeR/limma have k_config = NA.)
summ <- detail %>%
  group_by(method, cancer, stage, n_per_group, k_config) %>%
  summarise(n_samples = first(n_samples), runs = dplyr::n(),
            n_genes = round(mean(n_genes)),
            FPR_0.05    = round(mean(fpr_0.05), 4),
            FPR_0.05_sd = round(sd(fpr_0.05),   4),   # across-run spread (chance variability)
            FPR_0.01    = round(mean(fpr_0.01), 4),
            median_p    = round(mean(median_p), 3),
            hits_FDR05    = round(mean(hits_fdr05), 4),
            hits_FDR05_sd = round(sd(hits_fdr05),   4),
            # goodness-of-fit distance to Uniform(0,1), averaged over runs (0 = uniform)
            ks_D   = round(mean(ks_D,   na.rm = TRUE), 4),
            ad_A2  = round(mean(ad_A2,  na.rm = TRUE), 4),
            cvm_W2 = round(mean(cvm_W2, na.rm = TRUE), 4),
            .groups = "drop") %>%
  mutate(inflation_0.05 = ifelse(is.na(FPR_0.05), NA, round(FPR_0.05 / 0.05, 2)),
         verdict = dplyr::case_when(
           # all remaining methods carry raw p-values -> judged on FPR inflation.
           !is.na(inflation_0.05) & inflation_0.05 > 1.5  ~ "ANTI-CONSERVATIVE",
           !is.na(inflation_0.05) & inflation_0.05 < 0.67 ~ "conservative",
           !is.na(inflation_0.05)                         ~ "calibrated",
           hits_FDR05 > 0.02                              ~ "ANTI-CONSERVATIVE (calls under H0)",
           TRUE                                           ~ "calibrated")) %>%
  as.data.frame()
write.csv(summ, file.path(OUT_DIR, "null_calibration_summary.csv"), row.names = FALSE)

cat("\n", paste(rep("=", 116), collapse = ""), "\n", sep = "")
cat("NULL CALIBRATION -- all methods on random halves of a homogeneous group\n")
cat("Result = MEAN over `runs` independent random partitions; *_sd is the across-run SD (chance spread).\n")
cat("Individual per-run metrics are saved to null_calibration_per_run.csv (not printed here).\n")
cat("Under H0: FPR@alpha ~= alpha  AND  hits_FDR05 ~= 0.  inflation = FPR@.05/0.05 (~1 good, >1 over-calls).\n")
cat("ks_D / ad_A2 / cvm_W2 = distance of the p-values to Uniform(0,1) (KS / Anderson-Darling /\n")
cat("  Cramer-von Mises), averaged over runs; 0 = perfectly uniform. See QQ-plots in plots/.\n")
cat("n_per_group = replicates per fake group for that arm (full half + fixed subsamples).\n")
cat("k_config = TOX kNN neighbourhood config (k_start _ k_max, step 1); NA for edgeR/limma.\n")
cat("Cohorts include the 4 cancer stages AND healthy (matched normals), across CALIB_CANCERS.\n")
cat(paste(rep("=", 116), collapse = ""), "\n", sep = "")
print(summ[order(summ$cancer, summ$stage, summ$n_per_group, summ$k_config, summ$method), ], row.names = FALSE)

# The key new readout: how the kNN neighbourhood size affects TOX calibration, per arm
# (medians across cancers/stages/replicate-counts). Expect little effect for raw, more for log.
cat("\n--- kNN SWEEP: per (method x k_config) medians across cancers/stages/replicate-counts (TOX only) ---\n")
ksw <- summ %>% filter(!is.na(k_config)) %>% group_by(method, k_config) %>%
  summarise(FPR_0.05 = round(median(FPR_0.05, na.rm = TRUE), 4),
            ks_D  = round(median(ks_D,  na.rm = TRUE), 4),
            ad_A2 = round(median(ad_A2, na.rm = TRUE), 4),
            hits_FDR05     = round(median(hits_FDR05), 4),
            worst_FPR_0.05 = round(max(FPR_0.05, na.rm = TRUE), 4), .groups = "drop") %>% as.data.frame()
print(ksw[order(ksw$method, ksw$k_config), ], row.names = FALSE)

cat("\n--- per-method overall (median across ALL cells: cancers/stages/replicate-counts/k) ---\n")
ov <- summ %>% group_by(method) %>%
  summarise(FPR_0.05 = round(median(FPR_0.05, na.rm = TRUE), 4),
            hits_FDR05 = round(median(hits_FDR05), 4),
            worst_FPR_0.05 = round(max(FPR_0.05, na.rm = TRUE), 4),
            worst_hits_FDR05 = round(max(hits_FDR05), 4), .groups = "drop") %>% as.data.frame()
print(ov[order(ov$method), ], row.names = FALSE)

# --- NORMALISATION COMPARISON: counts vs TPM, everything else held fixed ---------
# The four NORM-* arms differ ONLY in the input scale and whether the composition
# shift was centred out; same samples, same genes, same split, same kNN config. Read
# them against each other, not against the main table.
norm_summ <- summ[grepl("^NORM-", summ$method), ]
if (nrow(norm_summ)) {
  cat("\n", paste(rep("=", 116), collapse = ""), "\n", sep = "")
  cat("NORMALISATION COMPARISON -- TOX-log on TPM vs TMM-CPM vs median-of-ratios counts\n")
  cat("  NORM-TPM          raw TPM, no composition correction (the assumption Delta = 0)\n")
  cat("  NORM-TPM-centred  raw TPM, per-gene MEDIAN of beta subtracted (median centring of log-ratios)\n")
  cat("  NORM-TMM          counts normalised with edgeR TMM factors -> CPM\n")
  cat("  NORM-MoR          counts normalised with DESeq2 median-of-ratios size factors\n")
  cat("Under H0 all four should be calibrated; the question is which is CLOSEST to uniform.\n")
  cat("A composition shift inflates FPR via Delta/sigma, so if TPM is materially worse than\n")
  cat("the count arms, that is the shift showing up -- and NORM-TPM-centred says whether\n")
  cat("median centring recovers it without leaving TPM.\n")
  cat(paste(rep("=", 116), collapse = ""), "\n", sep = "")
  print(norm_summ[order(norm_summ$cancer, norm_summ$stage, norm_summ$method),
                  c("method", "cancer", "stage", "n_samples", "n_per_group", "n_genes",
                    "FPR_0.05", "FPR_0.05_sd", "FPR_0.01", "median_p", "hits_FDR05",
                    "ks_D", "ad_A2", "inflation_0.05", "verdict")],
        row.names = FALSE)

  cat("\n--- normalisation arms: median across the (cancer x stage) cells ---\n")
  nov <- norm_summ %>% group_by(method) %>%
    summarise(cells = dplyr::n(),
              FPR_0.05 = round(median(FPR_0.05, na.rm = TRUE), 4),
              FPR_0.01 = round(median(FPR_0.01, na.rm = TRUE), 4),
              ks_D     = round(median(ks_D,  na.rm = TRUE), 4),
              ad_A2    = round(median(ad_A2, na.rm = TRUE), 4),
              hits_FDR05 = round(median(hits_FDR05), 4),
              worst_FPR_0.05 = round(max(FPR_0.05, na.rm = TRUE), 4), .groups = "drop") %>%
    as.data.frame()
  print(nov[order(nov$ks_D), ], row.names = FALSE)
  cat("(sorted by ks_D: smallest = closest to Uniform(0,1) = best calibrated)\n")
}

# ==================== plots ====================
pval_pool <- Filter(Negate(is.null), pval_pool)          # make_pval returns NULL for empty runs
pv_all <- if (length(pval_pool)) as.data.frame(rbindlist(pval_pool)) else NULL

# The per-cohort money/hist/QQ plots can't also facet the kNN grid without exploding, so
# they show a single REFERENCE k config (plus the k-independent refs, k_config = NA). The
# full k grid lives in the CSV, the kNN-sweep table above, and the dedicated kNN plot below.
REF_K   <- if ("k20_50" %in% summ$k_config) "k20_50" else stats::na.omit(unique(summ$k_config))[1]
summ_ref <- summ[is.na(summ$k_config) | summ$k_config == REF_K, ]
pv_ref   <- if (!is.null(pv_all)) pv_all[is.na(pv_all$k_config) | pv_all$k_config == REF_K, ] else NULL

# KNN-SWEEP PLOT: how neighbourhood size moves TOX calibration, per arm (TOX only) -----
ksw_long <- summ[!is.na(summ$k_config), ]
if (nrow(ksw_long)) {
  p_ksw <- ggplot(ksw_long, aes(k_config, ks_D, colour = factor(n_per_group), group = factor(n_per_group))) +
    stat_summary(fun = median, geom = "line") + stat_summary(fun = median, geom = "point") +
    facet_wrap(~ method) +
    labs(title = "kNN sweep: p-value non-uniformity vs neighbourhood size (median over cancers/stages)",
         subtitle = "y = KS distance to Uniform(0,1) (0 = perfect); expect ~flat for raw, more movement for log",
         x = "kNN config (k_start _ k_max)", y = "ks_D (median)", colour = "reps/group") +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                         legend.position = "bottom")
  tryCatch(ggsave(file.path(PLOT_DIR, "null_knn_sweep.png"), p_ksw, width = 13, height = 9, dpi = 140),
           error = function(e) NULL)
}

# MONEY PLOT 1: false calls under H0 vs replicate count (must be ~0 for every method) -
p_hits <- ggplot(summ_ref, aes(factor(n_per_group), hits_FDR05, fill = method)) +
  geom_col(position = position_dodge(preserve = "single")) +
  geom_hline(yintercept = 0.05, linetype = 3, colour = "grey40") +
  facet_grid(cancer ~ stage) +
  labs(title = "False calls under H0 (fraction significant at each method's FDR<0.05)",
       subtitle = sprintf("must be ~0. TOX shown at k_config=%s (see CSV for the full kNN grid).", REF_K),
       x = "replicates per group", y = "hits @ FDR<0.05 (H0)", fill = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
tryCatch(ggsave(file.path(PLOT_DIR, "null_false_calls.png"), p_hits, width = 12, height = 8, dpi = 140),
         error = function(e) NULL)

# MONEY PLOT 2: raw-p FPR inflation vs replicate count (p-value methods) ------------
sp2 <- summ_ref[!is.na(summ_ref$inflation_0.05), ]
p_inf <- ggplot(sp2, aes(factor(n_per_group), inflation_0.05, fill = method)) +
  geom_col(position = position_dodge(preserve = "single")) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid(cancer ~ stage) +
  labs(title = "Type-I calibration under H0: FPR@0.05 / 0.05  (1.0 = perfect)",
       subtitle = ">1 = anti-conservative (raw p-values skewed toward 0)",
       x = "replicates per group", y = "inflation", fill = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
tryCatch(ggsave(file.path(PLOT_DIR, "null_fpr_inflation.png"), p_inf, width = 12, height = 8, dpi = 140),
         error = function(e) NULL)

# p-value histograms, faceted method x replicate count (pooled over cohorts, since
# each panel must hold a single replicate count for the shape to be meaningful) ----
if (!is.null(pv_ref)) {
  p_hist <- ggplot(pv_ref, aes(p)) +
    geom_histogram(boundary = 0, bins = 20, fill = "steelblue", colour = "white") +
    facet_grid(method ~ n_per_group, scales = "free_y") +
    labs(title = "Null p-value distributions (should be UNIFORM under H0)",
         subtitle = sprintf("columns = replicates per group; TOX at k_config=%s; left spike = anti-conservative", REF_K),
         x = "raw p-value", y = "count") +
    theme_minimal(base_size = 8)
  tryCatch(ggsave(file.path(PLOT_DIR, "null_pvalue_histograms.png"), p_hist,
                  width = 16, height = 10, dpi = 130), error = function(e) NULL)
}

# QQ-plots: observed p vs expected uniform quantile (should hug the diagonal) ------
if (!is.null(pv_ref)) {
  qq_points <- function(p, npts = 2000L) {
    p <- sort(p[is.finite(p)]); n <- length(p)
    if (n < 2L) return(NULL)
    idx <- unique(round(seq(1, n, length.out = min(npts, n))))
    data.frame(expected = (idx - 0.5) / n, observed = p[idx])
  }
  qq <- do.call(rbind, lapply(
    split(pv_ref, list(pv_ref$method, pv_ref$n_per_group), drop = TRUE),
    function(df) { q <- qq_points(df$p)
                  if (is.null(q)) NULL else
                    cbind(method = df$method[1], n_per_group = df$n_per_group[1], q) }))
  if (!is.null(qq) && nrow(qq)) {
    p_qq <- ggplot(qq, aes(expected, observed)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
      geom_point(size = 0.3, alpha = 0.5, colour = "steelblue") +
      facet_grid(method ~ n_per_group) + coord_equal() +
      labs(title = "Null p-value QQ-plots (observed vs Uniform(0,1); on the diagonal = calibrated)",
           subtitle = "columns = replicates per group; bowing ABOVE the diagonal at small expected = anti-conservative",
           x = "expected p (uniform quantile)", y = "observed p") +
      theme_minimal(base_size = 8)
    tryCatch(ggsave(file.path(PLOT_DIR, "null_pvalue_qqplots.png"), p_qq,
                    width = 12, height = 10, dpi = 130), error = function(e) NULL)
  }
}

cat("\n", paste(rep("=", 116), collapse = ""), "\n", sep = "")
cat("DONE. CSVs + plots in:", OUT_DIR, "\n")
cat("  * all methods ~calibrated        -> dataset is clean; precision gap is power/philosophy.\n")
cat("  * only TOX inflated              -> TOX null bug.\n")
cat("  * edgeR/limma/DESeq2 inflated    -> the DEG references over-call under H0 (contain false positives).\n")
cat(paste(rep("=", 116), collapse = ""), "\n", sep = "")
