#!/usr/bin/env Rscript
# power_test.R
# =============================================================================
# POWER BENCHMARK for TOX vs edgeR / limma-voom / DESeq2 -- ANALYSIS script.
#
#   Rscript TCGA_test/scripts/power_test.R                 # all parts (R only if POWER_REAL_RDS)
#   POWER_PARTS=Q,C Rscript TCGA_test/scripts/power_test.R # a subset
#   POWER_REAL_RDS=coad_h.rds,coad_s4.rds,... POWER_PARTS=R Rscript TCGA_test/scripts/power_test.R
#   Rscript TCGA_test/scripts/power_report.R               # metrics, tables, plots
#
# Two scripts (docs/power_test_implementation_plan.md section 1):
#   power_test.R    simulate / thin, run every method, SAVE RAW PER-GENE RESULTS. No metrics.
#   power_report.R  reads power_raw/, computes every metric (common/power_metrics.R),
#                   writes power_out/*.csv and every plot. Re-runnable without Fortran.
#
# Companion of calibration_test.R (sourced for draw_counts, common_filter and the
# edgeR / limma / DESeq2 runners; its main block does not run when sourced). Design
# rationale: docs/power_benchmark_design.md. Run from the Tensor-Omics root (rcpp/
# and external/ resolve relative to the working directory); calibration_test.R and
# common/ are found relative to THIS script's location (COMMON_DIR below). Every TOX
# p-value comes from the compiled Fortran.
#
# PARTS
#   R   REAL-DATA SIGNAL INJECTION -- the primary evidence. Each base cohort (an RDS
#       from export_power_cohort.R; plan: COAD / LUAD / KIRC x healthy / Stage IV) is
#       used WHOLE: drop one sample at random if S is odd, split at random into halves
#       A / B (n = floor(S/2) each), inject known log2 FCs by binomial thinning (Gerard
#       2020), balanced so both halves lose the same expected relative abundance
#       (Delta = 0, nulls null in TPM). n_rounds repeats, each a fresh split AND a fresh
#       DE draw. Runs the input x arm grid (as Part N). Checks: gene-length recovery
#       (count/TPM must be constant per gene), one sample per patient (TCGA barcodes),
#       realised vs intended log2 FC per injected gene.
#   RS  REAL DATA AT SMALL n (note D WP3): the Part R cohorts, but each round draws
#       2n samples at random (n = 3, 5 per group; disjoint groups), null arm (pi1 = 0,
#       calibration) and power arm (pi1 = 0.1), 50 rounds (POWER_ROUNDS_RS). Short
#       arm list: limma / edgeR / DESeq2, TOX-log|TPM-centred and the variance arms.
#   P   Main simulated grid: nb / lnpois / tpois x n_rep, balanced effect draw (I3).
#   P+  Dispersion heterogeneity (I2): phi_g = phi_trend(mu_g) exp(sigma_d z - sigma_d^2/2),
#       phi_trend(mu) = phi_inf + c / mu, sigma_d in {0, .5, 1, 2} x sigma_hat. Null arm
#       (pi1 = 0) and power arm (pi1 = 0.1). Parameters from estimate_dispersion_trend.R
#       (POWER_DISP_PARAMS=<rds>); built-in placeholders otherwise (manifest: disp_source).
#   Q   pi1 x effect heterogeneity + oracle (true-null-only) pools + pool diagnostics
#       (tox_diagnose on a gene subset, gated by validate_against_fortran).
#   C   Composition stress, Delta-TARGETED (I4): target Delta in {0, .25, .5, 1}, scored
#       against the biological truth beta. Breakdown point computed by the report.
#   S   STRESS ONLY: `bimodal`. Not for ranking methods.
#   N   Input scale x correction under nuisance (depth, dispersion model).
#   V   Variance-only nulls (case-group dispersion x hv_factor).
#   U   Unequal group variances (note E v2 section 5, note F section 7): P+ at
#       sigma_d = sigma_hat with the case-group dispersion x rho^2 for EVERY gene
#       (case extra-Poisson SD = rho x control), rho in {1, 1.25, 1.5, 2}; tests the
#       equal-variance assumption of the cross-fitted scaling. Cross-fit arms only.
#   D   DESeq2 diagnostic (I7): nb, null and power arms, DESeq2 default vs
#       independentFiltering = TRUE vs cooksCutoff = FALSE. No TOX.
#
# RAW OUTPUT (PCFG$raw_dir/<part>/; not committed)
#   results_<job>.csv.gz   one row per (method, gene): p, padj, lfc, log2fc, tie,
#                          floor_th, nbhd_own_case, nbhd_own_control, lfc_scale; the
#                          variance-aware arms add vb_u, vb_scan_case/_control,
#                          vb_genes_case/_control, vb_k, vb_tau, vb_scatter,
#                          vb_n_floored, vb_scan_radius, vb_lambda_pool/_case/_control,
#                          vb_p_hat, vb_participation_ratio, vb_stop_case/_control,
#                          vb_radius_case (NA for every other method)
#   truth_<job>.csv.gz     one row per gene: is_de, true_lfc, expr, kept, is_hv
#                          (+ realised_lfc in Part R)
#   diag_<job>.csv.gz      Part Q: pool diagnostics on a gene subset
#   vdiag_<job>.csv.gz     parts with variance-aware arms: per gene s2_A, s2_B, s2,
#                          m_A, m_B, m on log2(TPM + 1) (note G section 2)
#   manifest.csv           one row per job: cell keys, round, seed, Delta, effect-draw
#                          diagnostics, Part R cohort fields, variance-band scalars
#                          (vb_*: band width, lambdas, p_hat), dispersion
#                          parameter file (disp_*), data_hash, elapsed_sec, raw_bytes,
#                          skipped (cells with 2n > S), git commit
#   pcfg_dump.txt          the full PCFG of the run
# =============================================================================

# ------------------------------------------------------------ dependencies

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

.find_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("None of these files found (run from the calibration_test.R ",
                         "working directory): ", paste(candidates, collapse = ", "))
  hit[1]
}

# Simulator (draw_counts), common_filter, the edgeR / limma / DESeq2 runners,
# load_or_install, LIB_DIR, CFG$run_tox, and the TOX wrappers.
source(.find_file(c(file.path(COMMON_DIR, "..", "Simulated_data", "scripts", "calibration_test.R"),
                     "calibration_test.R")))
source(file.path(COMMON_DIR, "tox_null_reimpl.R"))      # Part Q pool diagnostics
source(file.path(COMMON_DIR, "power_constants.R"))      # POWER_LIB_SIZE
source(file.path(COMMON_DIR, "print_utils.R"))          # print_tab (note C)
source(file.path(COMMON_DIR, "variance_diagnostics.R")) # vdiag_frame (note G)
load_or_install("parallel")

#' gzip-compressed CSV via base R (data.table::fwrite compression needs a zlib-enabled build).
write_gz <- function(x, f) { con <- gzfile(f, "w"); on.exit(close(con)); write.csv(x, con, row.names = FALSE) }

# ------------------------------------------------------------ configuration

ALL_PARTS <- c("R", "RS", "P", "P+", "U", "Q", "C", "S", "N", "V", "D")
REAL_PARTS <- c("R", "RS")                   # parts drawn from the real cohorts
# The seed's part index. FROZEN in the order parts were added (RS last), so adding a
# part never changes the datasets of the existing ones; ALL_PARTS only sets the order.
SEED_PARTS <- c("R", "P", "P+", "Q", "C", "S", "N", "V", "D", "RS", "U")
.env_parts <- Sys.getenv("POWER_PARTS", paste(ALL_PARTS, collapse = ","))
#' Comma-separated environment list, `default` when unset.
env_list <- function(var, default) {
  v <- Sys.getenv(var, ""); if (nzchar(v)) trimws(strsplit(v, ",")[[1]]) else default
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Mean-dependent dispersion for P+ (I2). Estimated from the TCGA healthy cohorts by
# estimate_dispersion_trend.R; sigma_hat is the RAW gene-wise scatter, an UPPER bound
# on true heterogeneity (sampling noise inflates it). The fallback values are
# PLACEHOLDERS, recorded as disp_source = "placeholder" in the manifest.

#' Guard against a stale parameter file (note D WP2): mu is defined at the depth the
#' estimator saved as mu_depth, which must equal the simulator's lib_size. A file
#' without mu_depth predates the guard and counts as stale.
check_disp_params <- function(d, lib_size, path = "<rds>") {
  md <- d$mu_depth
  if (!is.numeric(md) || length(md) != 1L || !is.finite(md) || abs(md / lib_size - 1) >= 1e-9)
    stop(sprintf("%s: mu_depth = %s but lib_size = %g -- stale dispersion parameters; re-run estimate_dispersion_trend.R",
                 path, if (is.null(md)) "missing" else format(md), lib_size), call. = FALSE)
  d
}
.disp_params <- local({
  f <- Sys.getenv("POWER_DISP_PARAMS", "")
  if (nzchar(f)) { d <- check_disp_params(readRDS(f), POWER_LIB_SIZE, f); d$source <- normalizePath(f); d }
  else list(phi_inf = 0.15, c = 2, sigma_hat = 0.8, source = "placeholder")
})

PCFG <- list(
  parts      = trimws(strsplit(.env_parts, ",")[[1]]),
  # Env overrides for a smoke test: POWER_GENES=1000 POWER_ROUNDS=2 POWER_CORES=4
  n_genes    = as.integer(Sys.getenv("POWER_GENES",  "5000")),
  n_rounds   = as.integer(Sys.getenv("POWER_ROUNDS", "20")),  # >= 20 for anything reported
  n_cores    = as.integer(Sys.getenv("POWER_CORES",  "32")),
  seed       = 20260923L,
  lib_size   = POWER_LIB_SIZE,
  disp       = 0.2,

  # |log2 FC| ~ log-uniform on this range, sign symmetric, then BALANCED (I3): the
  # side that moves more relative abundance is shrunk by a common factor.
  lfc_range  = c(0.25, 4),
  frac_up    = 0.5,

  # Part grids. A cell = one combination; each cell runs n_rounds datasets.
  R = list(input_rds = Sys.getenv("POWER_REAL_RDS", ""), pi1 = 0.10, min_mean_count = 10,
           len_tol = 0.05,     # max per-gene IQR of log2(count/TPM), after sample scaling
           lfc_tol = 0.5,      # |realised - intended log2 FC| flagged above this
           # E1 (note D WP7, run only on request): keep only genes with no zero count in
           # any sample of the cohort, for ALL methods of Parts R and RS.
           zero_free = Sys.getenv("POWER_R_ZERO_FREE", "0") == "1"),
  # Part RS (note D WP3): n per group, rounds per cell (own setting; the global
  # n_rounds stays the default of every other part), TOX comparator = Part R's.
  RS = list(n_reps = c(3L, 5L, 7L, 10L), pi1 = c(0, 0.10),     # n = 7, 10: note E v2 section 5
            n_rounds = as.integer(Sys.getenv("POWER_ROUNDS_RS", "50")),
            # POWER_RS_GRID=1: Part RS runs Part R's input grid (plain grid arms and
            # vmatch$grid_arms on every input) instead of the single comparator below.
            grid = Sys.getenv("POWER_RS_GRID", "0") == "1",
            tox_arms = list(list(label = "TOX-log|TPM-centred", engine = "exact", norm = 1L, null = 0L,
                                 input = "tpm", centre = TRUE))),
  P = list(dists = c("nb", "lnpois", "tpois"), n_reps = c(3L, 5L, 10L), pi1 = 0.10),
  `P+` = list(dists = c("nb", "lnpois"), n_reps = c(3L, 5L, 10L), pi1 = c(0, 0.10),
              disp_model = "hetero", sigma_mult = c(0, 0.5, 1)),     # sigma_d 0 / 0.53 / 1.05 (note F)
  U = list(dists = c("nb", "lnpois"), n_reps = c(3L, 5L, 10L), pi1 = c(0, 0.10),
           disp_model = "hetero", sigma_mult = 1, rho = c(1, 1.25, 1.5, 2)),
  Q = list(dists = "lnpois", n_reps = 5L,
           pi1 = c(0.01, 0.05, 0.10, 0.30), het = c(1, 0.5)),
  C = list(dists = c("nb", "lnpois"), n_reps = c(3L, 5L, 10L), pi1 = 0.10,
           delta_target = c(0, 0.25, 0.5, 1.0)),
  S = list(dists = "bimodal", n_reps = c(3L, 5L, 10L), pi1 = 0.10),
  N = list(dists = "lnpois", n_reps = c(3L, 5L, 10L), pi1 = 0.10,
           depth = c("equal", "random", "confounded"),
           disp_model = c("const", "trend", "gene")),
  V = list(dists = "lnpois", n_reps = c(3L, 5L, 10L), pi1 = 0.10, het = c(1, 0.5),
           hv_frac = 0.05),
  D = list(dists = "nb", n_reps = c(3L, 5L, 10L), pi1 = c(0, 0.10)),

  disp_hetero = .disp_params,

  # Nuisance models (Parts N, V). Defaults of the other parts are unchanged:
  # equal depth, constant dispersion, no variance-only nulls.
  depth_sd       = 0.3,           # "random": log library size ~ N(log lib_size, 0.3)
  depth_confound = 1.5,           # "confounded": group B libraries this much deeper
  disp_trend     = c(a = 0.02, b = 5), disp_scatter = 0.3,   # phi = (a + b/mu) * lognormal
  disp_gene_sd   = 0.7,           # "gene": sd of log phi_g around log(disp)
  hv_factor      = 4,             # Part V: case-group dispersion multiplier

  # Parts N and R (and RS with POWER_RS_GRID=1): input x arm grid. src: "tpm" |
  # "cnt" | "tmm" | "mor". POWER_GRID_INPUTS / POWER_GRID_ARMS restrict both lists
  # (e.g. the lean yeast run); the code of every arm stays.
  inputs = list(
    list(name = "TPM",            src = "tpm", centre = FALSE),
    list(name = "TPM-centred",    src = "tpm", centre = TRUE),
    list(name = "counts",         src = "cnt", centre = FALSE),
    list(name = "counts-centred", src = "cnt", centre = TRUE),
    list(name = "TMM",            src = "tmm", centre = FALSE),
    list(name = "MoR",            src = "mor", centre = FALSE)
  ),
  # TOX-log-boot is left out: it matched TOX-log in every Part P cell.
  grid_arms = env_list("POWER_GRID_ARMS", c("TOX-log", "TOX-log-blocked", "TOX-raw", "TOX-raw-blocked")),
  n_boot = as.integer(Sys.getenv("POWER_N_BOOT", "10000")),   # variance arms with nullc = "bootstrap"

  run_edger  = TRUE,
  run_limma  = TRUE,
  run_deseq2 = TRUE,           # the slow one

  # TOX: the configuration null_calibration.R fixed on evidence (k20_50, tau 0.1).
  kcfg     = list(k_start = 20L, k_step = 1L, k_max = 50L, name = "k20_50"),
  tau      = 0.1,
  max_pool = 70000L,

  # engine: "exact" = tox_noise_model_exact (sqrt-scaled pairwise, deterministic)
  #         "bootstrap" = tox_noise_model (null = 0 pooled iid; null = 1 gene-blocked,
  #                        enumerated exactly for n_rep <= 5 at k_max = 50, sampled above)
  # input:  "tpm" | "tpm_centred" (median-centred obs) | "tmm" (TMM-CPM, linear)
  tox_arms = list(
    list(label = "TOX-log",         engine = "exact",     norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-raw",         engine = "exact",     norm = 0L, null = 0L, input = "tpm"),
    list(label = "TOX-log-boot",    engine = "bootstrap", norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-log-blocked", engine = "bootstrap", norm = 1L, null = 1L, input = "tpm"),
    list(label = "TOX-raw-blocked", engine = "bootstrap", norm = 0L, null = 1L, input = "tpm")
  ),
  # Part C only: the input-scale comparison (TOX-log, exact, pooled).
  tox_input_arms = list(
    list(label = "TOX-log-TPM",         engine = "exact", norm = 1L, null = 0L, input = "tpm"),
    list(label = "TOX-log-TPM-centred", engine = "exact", norm = 1L, null = 0L, input = "tpm_centred"),
    list(label = "TOX-log-TMM",         engine = "exact", norm = 1L, null = 0L, input = "tmm")
  ),
  # Part Q only: these tox_arms are ALSO run with a true-null-only pool, as
  # "<label>-oracle". Must use input "tpm" or "tmm" (median centring would be
  # computed over a different gene set in the oracle calls).
  oracle_arms = c("TOX-log", "TOX-log-blocked"),

  # Variance-aware arms (implementation notes A and F; tox_compute_noise_pvalues_pipeline_variance,
  # eval doc §21/§22). Note F section 7, first round:
  #   (a) TOX-log-vband-c3       the band as it is now (raw match, c = 3)
  #       TOX-log-vband-js-c*    band on the James-Stein-shrunk variance (replaces note B's vshrunk)
  #   (b) TOX-log-xfit-raw       cross-fitted scaling, raw variances (lambda = 0)
  #   (c) TOX-log-xfit-js        cross-fitted scaling, James-Stein, median target   PRIMARY
  #       TOX-log-xfit-js-loess  same, LOESS target
  #       TOX-log-xfit-js-nocf   negative control: same-group scaling (RS n = 3, 5 and P+ only)
  #       TOX-log-adapt-r*       adaptive variance-aware growth (raw match, c = 3)
  #       TOX-log-adapt-xfit-js  adaptive growth (JS match, c_js = 1) + cross-fit JS
  # "mean" arms use the comparator's neighbourhood (kcfg, tau). In Parts R and RS
  # every arm takes the comparator's input, `TOX-log|TPM-centred` (centre = TRUE);
  # elsewhere the comparator is `TOX-log` (TPM). Free parameters (c_js, r_max) are
  # chosen on the pi1 = 0 rounds only (power_report.R, note E v2 section 6).
  # `parts` restricts an arm; `rs_n` restricts it to these n in Part RS.
  vmatch = list(k = 50L, span = 0.3, degree = 1L, parts = c("P+", "V", "R", "RS", "U"),
                # Parts R (and RS with the grid): these variance arms ALSO run on every other
                # grid input as "<label>|<input>" (their own input, TPM-centred, is the plain
                # label). They also export the neighbour genes of their false positives and
                # false negatives (nbhd_<job>.csv.gz).
                grid_arms = env_list("POWER_VGRID_ARMS", "TOX-log-xfit-js")),
  adapt  = list(block_size = 5L, k_min = 20L, tau_stab = 0.10, w_blocks = 2L, a_min = 0.2),
  vmatch_arms = local({
    four <- c("P+", "V", "R", "RS")
    v <- function(label, ...) c(list(label = label, engine = "variance", norm = 1L, null = 0L, input = "tpm"),
                                list(...))
    list(
      v("TOX-log-vband-c3",      vmode = "band", match = "raw", scaling = "none", c_tau = 3, parts = four),
      v("TOX-log-vband-js-c0.5", vmode = "band", match = "js",  scaling = "none", c_js = 0.5, parts = four),
      v("TOX-log-vband-js-c1",   vmode = "band", match = "js",  scaling = "none", c_js = 1.0, parts = four),
      v("TOX-log-xfit-raw",      vmode = "mean", scaling = "xfit_raw"),
      v("TOX-log-xfit-js",       vmode = "mean", scaling = "xfit_js", target = "median"),
      v("TOX-log-xfit-js-loess", vmode = "mean", scaling = "xfit_js", target = "loess"),
      # cross-fit JS with the bootstrap model's Monte Carlo difference-of-means null
      # (PCFG$n_boot draws, floor 1/(n_boot + 1)); real data only (runtime)
      v("TOX-log-xfit-js-boot",  vmode = "mean", scaling = "xfit_js", target = "median", nullc = "bootstrap",
        parts = c("R", "RS")),
      v("TOX-log-xfit-js-nocf",  vmode = "mean", scaling = "same_group", target = "median",
        parts = c("RS", "P+"), rs_n = c(3L, 5L)),
      v("TOX-log-adapt-r0.25",   vmode = "adaptive", match = "raw", scaling = "none", c_tau = 3, r_max = 0.25, parts = four),
      v("TOX-log-adapt-r0.5",    vmode = "adaptive", match = "raw", scaling = "none", c_tau = 3, r_max = 0.5,  parts = four),
      v("TOX-log-adapt-r1",      vmode = "adaptive", match = "raw", scaling = "none", c_tau = 3, r_max = 1.0,  parts = four),
      v("TOX-log-adapt-xfit-js", vmode = "adaptive", match = "js",  scaling = "xfit_js", c_js = 1.0, r_max = 0.5,
        parts = four)
    )
  }),
  # Part Q pool diagnostics (exact engine, pooled, log): genes per job, and genes
  # the R re-implementation is checked on against the Fortran first.
  diag = list(n_genes = 200L, n_check = 200L),

  raw_dir = Sys.getenv("POWER_RAW_DIR", if (exists("TOX_TEST_DIR"))
              file.path(dirname(TOX_TEST_DIR), "power_raw") else "power_raw")
)

PCFG$inputs <- Filter(function(i) i$name %in% env_list("POWER_GRID_INPUTS", i$name), PCFG$inputs)

RUN_TOX <- isTRUE(CFG$run_tox)
if (!RUN_TOX) message("TOX wrappers not loaded -- running reference methods only.")

# =============================================================================
# 1. TRUTH AND SIMULATION
# =============================================================================

#' Draw the DE design shared by every part: which genes, sign, |log2 FC|.
draw_effects <- function(G, n_de, frac_up = PCFG$frac_up, lfc_range = PCFG$lfc_range) {
  beta <- numeric(G)
  de <- if (n_de > 0) sort(sample.int(G, n_de)) else integer(0)
  if (n_de > 0) {
    sgn <- ifelse(runif(n_de) < frac_up, 1, -1)
    beta[de] <- sgn * exp(runif(n_de, log(lfc_range[1]), log(lfc_range[2])))
  }
  list(beta = beta, de = de, is_de = seq_len(G) %in% de)
}

#' Balance a signed log2 FC vector q over abundance weights w (I3 / Part R).
#' The relative abundance the up side moves must equal what the down side moves;
#' the side that moves more has its log2 FCs shrunk by one common factor s in (0, 1]
#' (uniroot; a root exists whenever both sides are non-empty, since that side's mass
#' -> 0 as s -> 0). Signs never change.
#'   composition (thinning = FALSE): up gains sum w (2^q - 1), down loses sum w (1 - 2^q)
#'   thinning    (thinning = TRUE):  group A loses sum_up w (1 - 2^-q), group B sum_dn w (1 - 2^q)
#' ok = FALSE when only one side is non-empty (cannot balance).
balance_effects <- function(q, w, thinning = FALSE) {
  up <- q > 0; dn <- q < 0
  if (!any(up) || !any(dn)) return(list(q = q, ok = !any(up | dn)))
  mu <- function(s) sum(w[up] * (if (thinning) 1 - 2^(-s * q[up]) else 2^(s * q[up]) - 1))
  md <- function(s) sum(w[dn] * (1 - 2^(s * q[dn])))
  a <- mu(1); b <- md(1)
  if (a > b)      q[up] <- q[up] * uniroot(function(s) mu(s) - b, c(0, 1), tol = 1e-12)$root
  else if (b > a) q[dn] <- q[dn] * uniroot(function(s) md(s) - a, c(0, 1), tol = 1e-12)$root
  list(q = q, ok = TRUE)
}

#' Effect-draw diagnostics (manifest; I3). `eff` = realised TPM-space log2 FC of the
#' DE genes, `drawn` = the draw before balancing.
effect_diag <- function(eff, drawn, lfc_range = PCFG$lfc_range) {
  if (!length(eff)) return(list(eff_lfc_median_abs = NA_real_, frac_sign_flip = NA_real_,
                                frac_below_floor = NA_real_))
  list(eff_lfc_median_abs = median(abs(eff)),
       frac_sign_flip     = mean(sign(eff) != sign(drawn)),
       frac_below_floor   = mean(abs(eff) < lfc_range[1]))
}

#' Per-sample composition matrices for both groups.
#'
#' mode = "relative": balanced draw (I3), so sum_DE pi_a 2^beta = sum_DE pi_a and
#'   every null keeps pi_B == pi_a exactly (Delta = 0) with NO shift of the DE effects.
#'   The per-sample renormalisation within the DE set is kept for het < 1 (a random
#'   responder subset is balanced only in expectation); at het = 1 it is a no-op.
#' mode = "absolute": DE genes change absolute output, then the whole sample is
#'   renormalised, so every null moves by -Delta in TPM; true_lfc = beta (biology).
#'   delta_target = NULL keeps the raw draw (Part N); a number (Part C, I4) balances
#'   the draw (Delta = 0) and then scales the UP side's log2 FCs by s >= 1 until
#'   log2 sum(pi_a 2^beta) = delta_target.
#'
#' het < 1: each (DE gene, case sample) responds with probability het; a
#'   non-responder keeps pi_A. This inflates the DE gene's WITHIN-case variance,
#'   which is the only route by which DE genes can contaminate TOX's residual pools.
make_truth_power <- function(G, pi1, n_rep, het = 1, mode = c("relative", "absolute"),
                             frac_up = PCFG$frac_up, delta_target = NULL) {
  mode <- match.arg(mode)
  pi_a <- rlnorm(G, 0, 2); pi_a <- pi_a / sum(pi_a)
  ef <- draw_effects(G, round(pi1 * G), frac_up)
  de <- ef$de; drawn <- ef$beta
  beta <- if (mode == "relative" || !is.null(delta_target)) balance_effects(drawn, pi_a)$q else drawn
  if (mode == "absolute" && !is.null(delta_target) && delta_target > 0) {
    up <- beta > 0; bu <- beta[up]
    base <- sum(pi_a[!up] * 2^beta[!up])
    s <- uniroot(function(s) log2(base + sum(pi_a[up] * 2^(s * bu))) - delta_target,
                 c(1, 2), extendInt = "upX", tol = 1e-12)$root
    beta[up] <- s * bu
  }

  resp <- matrix(runif(G * n_rep) < het, G, n_rep)      # only DE rows are used
  pi_B <- matrix(pi_a, G, n_rep)
  for (j in seq_len(n_rep)) {
    eff <- beta * resp[, j]
    if (mode == "relative") {
      if (length(de)) {
        w <- pi_a[de] * 2^eff[de]
        pi_B[de, j] <- w * sum(pi_a[de]) / sum(w)
      }
    } else {
      w <- pi_a * 2^eff
      pi_B[, j] <- w / sum(w)
    }
  }

  if (mode == "relative") {
    cs <- if (length(de)) sum(pi_a[de]) / sum(pi_a[de] * 2^beta[de]) else 1
    true_lfc <- ifelse(ef$is_de, beta + log2(cs), 0)
    stopifnot(all(abs(sweep(pi_B[!ef$is_de, , drop = FALSE], 1, pi_a[!ef$is_de])) == 0))
    delta <- 0
  } else {
    true_lfc <- beta
    delta <- log2(sum(pi_a * 2^beta))       # TPM shift of every null gene is -delta
  }
  c(list(pi_a = pi_a, pi_B = pi_B, true_lfc = true_lfc, is_de = ef$is_de, delta = delta),
    effect_diag(true_lfc[de], drawn[de]))
}

#' Per-gene dispersion for P+ (I2): phi_trend(mu) times a MEAN-ONE lognormal
#' multiplier, so sigma_d changes the scatter around the trend, not the mean level.
disp_hetero <- function(mu, sigma_d, prm = PCFG$disp_hetero)
  (prm$phi_inf + prm$c / pmax(mu, 1e-3)) * exp(sigma_d * rnorm(length(mu)) - sigma_d^2 / 2)

#' Counts and TPM from per-sample compositions. Mirrors simulate_dataset():
#' fragments are drawn per gene with lambda_gj = N_j * pi_gj * L_g / sum(pi * L).
#'
#' depth:      "equal" (N_j = lib_size), "random" (lognormal per sample),
#'             "confounded" (group B = depth_confound x group A).
#' disp_model: "const" (phi = disp for every gene), "trend" (Part N: phi = a + b/mu,
#'             lognormal scatter), "gene" (log phi_g ~ N(log disp, disp_gene_sd)),
#'             "hetero" (P+: disp_hetero at sigma_d).
#' hv:         logical over genes; those get phi * hv_factor in group B only.
#' rho:        Part U: every gene's group-B (case) dispersion x rho^2, so the case
#'             extra-Poisson SD is rho x the control one.
#' The defaults draw NO extra random numbers.
simulate_power <- function(truth, dist, n_rep, lib_size = PCFG$lib_size, disp = PCFG$disp,
                           depth = "equal", disp_model = "const", hv = NULL, sigma_d = 0, rho = 1) {
  G <- length(truth$pi_a)
  L <- sample(round(exp(rnorm(3e4, log(2000), 0.8))), G, replace = TRUE)
  lib <- switch(depth,
    equal      = rep(lib_size, 2 * n_rep),
    random     = lib_size * exp(rnorm(2 * n_rep, 0, PCFG$depth_sd)),
    confounded = lib_size * rep(c(1, PCFG$depth_confound), each = n_rep),
    stop("unknown depth: ", depth))
  ec <- function(pi_mat, N) { p <- pi_mat * L; sweep(sweep(p, 2, colSums(p), "/"), 2, N, "*") }
  lamA <- ec(matrix(truth$pi_a, G, n_rep), lib[seq_len(n_rep)])
  lamB <- ec(truth$pi_B, lib[n_rep + seq_len(n_rep)])

  phi <- switch(disp_model,
    const  = disp,
    trend  = (PCFG$disp_trend[["a"]] + PCFG$disp_trend[["b"]] / pmax(rowMeans(lamA), 1e-3)) *
               exp(rnorm(G, 0, PCFG$disp_scatter)),
    gene   = exp(rnorm(G, log(disp), PCFG$disp_gene_sd)),
    hetero = disp_hetero(rowMeans(lamA), sigma_d),
    stop("unknown disp_model: ", disp_model))
  phiB <- if (is.null(hv)) phi else { p2 <- rep_len(phi, G); p2[hv] <- p2[hv] * PCFG$hv_factor; p2 }
  phiB <- phiB * rho^2

  # draw_counts takes a scalar or a per-gene vector (it recycles down the columns;
  # unit-tested in TCGA_test/tests/test_power_sim.R).
  cA <- draw_counts(lamA, dist, phi)
  cB <- draw_counts(lamB, dist, phiB)
  counts <- cbind(cA, cB)
  colnames(counts) <- c(paste0("A", seq_len(n_rep)), paste0("B", seq_len(n_rep)))
  rownames(counts) <- paste0("gene", seq_len(G))
  rate <- counts / L
  tpm  <- sweep(rate, 2, colSums(rate), "/") * 1e6

  list(counts = counts, tpm = tpm, lengths = L,
       group = rep(c("A", "B"), each = n_rep),
       true_lfc = truth$true_lfc, is_de = truth$is_de,
       is_hv = if (is.null(hv)) rep(FALSE, G) else hv,
       expr = log10(truth$pi_a * 1e6),          # strata by TRUE control abundance
       delta = truth$delta)
}

# ------------------------------------------------------ Part R: real cohort

#' Load one base cohort. RDS = list(counts = genes x samples, tpm = samples x genes
#' or genes x samples, [lengths], [project_id, stage]) from export_power_cohort.R.
#'
#' Checks (plan I10), each reported in the manifest:
#'   one sample per patient -- TCGA barcodes: patient = first 12 characters; later
#'     duplicates are dropped (n_dup_removed). Other ID formats: NOT checked
#'     (patient_check = "not_checked"), with a message.
#'   gene-length recovery   -- lengths = per-gene median of count/TPM after removing
#'     each sample's constant. Valid only if the TPM was built from THESE counts, so
#'     the per-gene IQR of log2(count/TPM) must be ~0; genes above len_tol are dropped
#'     (frac_len_bad).
load_real_cohort <- function(path) {
  d <- readRDS(path)
  cnt <- as.matrix(d$counts)
  # Everything here is genes x samples (edgeR convention; run_tox_arm transposes for the
  # Fortran). Raw TCGA `expression_vectors` are samples x genes and must be t()'d first.
  if (nrow(cnt) < ncol(cnt))
    stop(sprintf("%s: counts are %d x %d -- expected genes x samples", path, nrow(cnt), ncol(cnt)))
  if (is.null(rownames(cnt)) || is.null(colnames(cnt)))
    stop(path, ": counts need gene rownames and sample colnames")
  label <- if (!is.null(d$project_id)) paste0(d$project_id, "/", d$stage)
           else sub("\\.rds$", "", basename(path))

  ids <- colnames(cnt)
  if (all(grepl("^TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}", ids))) {
    dup <- duplicated(substr(ids, 1, 12))
    if (any(dup)) message(sprintf("%s: %d sample(s) from an already-used patient dropped: %s",
                                  label, sum(dup), paste(head(ids[dup], 5), collapse = ", ")))
    cnt <- cnt[, !dup, drop = FALSE]
    patient_check <- "tcga_barcode"; n_dup <- sum(dup)
  } else {
    message(sprintf("%s: sample IDs are not TCGA barcodes (e.g. %s) -- one-sample-per-patient ",
                    label, ids[1]), "check NOT done; confirm independence before using this cohort")
    patient_check <- "not_checked"; n_dup <- NA_integer_
  }

  frac_len_bad <- NA_real_
  if (!is.null(d$lengths)) {
    L <- as.numeric(d$lengths[rownames(cnt)])
  } else if (!is.null(d$tpm)) {
    tp <- as.matrix(d$tpm)
    if (setequal(rownames(tp), colnames(d$counts))) tp <- t(tp)    # samples x genes -> genes x samples
    miss <- c(setdiff(rownames(cnt), rownames(tp)), setdiff(colnames(cnt), colnames(tp)))
    if (length(miss))
      stop(sprintf("%s: %d count gene/sample IDs absent from tpm, e.g. %s",
                   path, length(miss), paste(head(miss, 3), collapse = ", ")))
    tp <- tp[rownames(cnt), colnames(cnt)]
    r <- cnt / tp; r[!is.finite(r) | r <= 0] <- NA
    r <- sweep(r, 2, apply(r, 2, median, na.rm = TRUE), "/")
    L <- apply(r, 1, median, na.rm = TRUE)
    spread <- apply(log2(r), 1, IQR, na.rm = TRUE)
    len_ok <- is.finite(spread) & spread <= PCFG$R$len_tol
  } else stop("Real cohort RDS needs $lengths or $tpm.")

  ok <- is.finite(L) & L > 0 & rowMeans(cnt) >= PCFG$R$min_mean_count
  if (exists("len_ok", inherits = FALSE)) {
    frac_len_bad <- mean(!len_ok[ok])
    ok <- ok & len_ok
  }
  n_zero_dropped <- NA_integer_
  if (isTRUE(PCFG$R$zero_free)) {                  # E1: zero-free universe
    zf <- rowSums(cnt == 0) == 0
    n_zero_dropped <- sum(ok & !zf); ok <- ok & zf
  }
  message(sprintf("%s: %d genes x %d samples; %d kept (length-recovery failures among expressed: %.2f%%)%s",
                  label, nrow(cnt), ncol(cnt), sum(ok), 100 * frac_len_bad,
                  if (is.na(n_zero_dropped)) "" else sprintf("; zero-free filter dropped %d", n_zero_dropped)))
  list(counts = cnt[ok, , drop = FALSE], lengths = L[ok], label = label, path = path,
       patient_check = patient_check, n_dup_removed = n_dup, frac_len_bad = frac_len_bad,
       n_zero_dropped = n_zero_dropped)
}

#' One signal-injected dataset from a real cohort. Part R (n = NULL): the WHOLE
#' cohort split into random halves (option A). Part RS: 2n distinct samples drawn at
#' random, the first n to group A, the next n to group B (disjoint by construction).
#' pi1 = 0 needs no special case: no DE gene, nothing balanced, every thinning
#' probability is 1 (counts unchanged), true_lfc = 0, delta = 0.
simulate_real <- function(cohort, pi1, n = NULL) {
  S <- ncol(cohort$counts); n <- n %||% (S %/% 2L)
  if (2L * n > S) stop(sprintf("2n = %d > S = %d samples", 2L * n, S))
  cols <- sample.int(S, 2L * n)                    # random halves; one dropped when S is odd
  cnt <- cohort$counts[, cols, drop = FALSE]; L <- cohort$lengths; G <- nrow(cnt)
  grp <- rep(c("A", "B"), each = n)
  tpm_of <- function(x) { r <- x / L; sweep(r, 2, colSums(r), "/") * 1e6 }
  lr <- function(t) log2(rowMeans(t[, grp == "B", drop = FALSE]) / rowMeans(t[, grp == "A", drop = FALSE]))

  tpm0 <- tpm_of(cnt)
  pi_hat <- rowMeans(tpm0) / 1e6                    # cohort composition
  ef <- draw_effects(G, round(pi1 * G))
  bal <- balance_effects(ef$beta, pi_hat, thinning = TRUE); q <- bal$q
  lr0 <- lr(tpm0)

  # Up-regulation of gene g thins group A by 2^-q, down-regulation thins group B by 2^q.
  thin <- function(y, p) matrix(rbinom(length(y), y, p), nrow(y))
  pa <- ifelse(q > 0, 2^-q, 1); pb <- ifelse(q < 0, 2^q, 1)
  cnt[, grp == "A"] <- thin(cnt[, grp == "A", drop = FALSE], pa)
  cnt[, grp == "B"] <- thin(cnt[, grp == "B", drop = FALSE], pb)
  tpm <- tpm_of(cnt)

  # Realised composition shift: median log2 ratio over true nulls (target 0), and
  # realised injected effect = change of the group log-ratio caused by thinning.
  lr1 <- lr(tpm)
  mA <- rowMeans(tpm[, grp == "A", drop = FALSE]); mB <- rowMeans(tpm[, grp == "B", drop = FALSE])
  okn <- !ef$is_de & mA > 1 & mB > 1
  realised <- ifelse(ef$is_de, lr1 - lr0, 0); realised[!is.finite(realised)] <- NA

  colnames(cnt) <- colnames(tpm) <- c(paste0("A", seq_len(n)), paste0("B", seq_len(n)))
  dg <- effect_diag(q[ef$is_de], ef$beta[ef$is_de])
  dev <- abs(realised - q)[ef$is_de]
  c(list(counts = cnt, tpm = tpm, lengths = L, group = grp,
         true_lfc = q, realised_lfc = realised, is_de = ef$is_de,
         is_hv = rep(FALSE, G), expr = log10(pi_hat * 1e6),
         delta = if (bal$ok) 0 else NA_real_, delta_realised = median(lr1[okn] - lr0[okn]),
         S = S, n_per_group = n, n_dropped = S - 2L * n, sample_ids = colnames(cohort$counts)[cols],
         frac_lfc_mismatch = if (length(dev)) mean(!is.finite(dev) | dev > PCFG$R$lfc_tol) else NA_real_),
    dg)
}

# =============================================================================
# 2. METHODS
# Each returns data.frame(gene, p, padj, lfc, log2fc, tie, floor_th, nbhd_own_case,
# nbhd_own_control, lfc_scale) -- the SAME columns for every method, because all
# methods' frames are rbind-ed into one results file. p = raw p (NA = untested);
# lfc = the method's own statistic; log2fc = its effect size (TOX: the model's
# log2((mean_B + 1)/(mean_A + 1)); edgeR/limma/DESeq2: their logFC); tie = |log2 FC
# of group means|, the p_eff tie-break.
# =============================================================================

group_lfc <- function(tpm, group)
  abs(log2((rowMeans(tpm[, group == "B", drop = FALSE]) + 1) /
           (rowMeans(tpm[, group == "A", drop = FALSE]) + 1)))

tmm_cpm <- function(counts)
  edgeR::cpm(normLibSizes(DGEList(counts = counts)), normalized.lib.sizes = TRUE, log = FALSE)

#' DESeq2 median-of-ratios: counts (genes x samples) divided by per-sample size factors.
mor_counts <- function(counts)
  sweep(counts, 2L, DESeq2::estimateSizeFactorsForMatrix(counts), "/")

#' Part D (I7): one DESeq() fit, three results() variants. "DESeq2" is identical to
#' calibration_test.R's run_deseq2 (independentFiltering = FALSE, default Cook's).
run_deseq2_variants <- function(counts, group) {
  dds <- DESeq(DESeqDataSetFromMatrix(round(counts), data.frame(group = factor(group)), ~ group),
               quiet = TRUE)
  v <- list(`DESeq2`        = list(independentFiltering = FALSE),
            `DESeq2-IF`     = list(independentFiltering = TRUE),
            `DESeq2-noCook` = list(independentFiltering = FALSE, cooksCutoff = FALSE))
  lapply(v, function(a) {
    r <- as.data.frame(do.call(results, c(list(dds), a)))
    data.frame(gene = rownames(counts), stat = r$pvalue, padj = r$padj,
               lfc = r$log2FoldChange, row.names = NULL)
  })
}

#' TOX via the production Fortran. `mat` = genes x samples, LINEAR scale.
#' `test_only` (row index): compute the p-value for that gene alone. The pool is
#' still built from ALL rows of `mat` -- only the tested set shrinks.
run_tox_arm <- function(mat, group, arm, tie, test_only = NULL) {
  G <- nrow(mat)
  case_rep <- t(mat[, group == "B", drop = FALSE])
  ctrl_rep <- t(mat[, group == "A", drop = FALSE])
  # The model computes the means and the observed statistic itself (mean
  # difference, or log2 Frechet-mean difference) and centres it on request.
  valid <- rep(1L, G)
  if (!is.null(test_only)) valid[-test_only] <- 0L
  centre <- as.integer(isTRUE(arm$centre) || identical(arm$input, "tpm_centred"))

  # Per-arm neighbourhood overrides (knn_config_test.R): arm$kcfg (k_start, k_step,
  # k_max), arm$tau, arm$k (band / adaptive / ellipse size), arm$adapt (adaptive settings).
  kc <- modifyList(PCFG$kcfg, arm$kcfg %||% list())
  tau <- arm$tau %||% PCFG$tau
  if (arm$engine == "variance") {
    ad <- modifyList(PCFG$adapt, arm$adapt %||% list())
    res <- tox_compute_noise_pvalues_pipeline_variance(
             case_replicates = case_rep, control_replicates = ctrl_rep,
             variance_mode = arm$vmode, match_variable = arm$match %||% "raw",
             scaling = arm$scaling %||% "none", shrink_target = arm$target %||% "median",
             valid_genes_own = valid, beta_centre = centre, norm_method = arm$norm,
             k_start = kc$k_start, k_step = kc$k_step, k_max = kc$k_max, tau = tau,
             k = arm[["k"]] %||% PCFG$vmatch$k, c_tau = arm$c_tau %||% 3, c_js = arm$c_js %||% 1,
             block_size = ad$block_size, k_min = min(ad$k_min, arm[["k"]] %||% PCFG$vmatch$k), r_max = arm$r_max %||% 0.5,
             tau_stab = ad$tau_stab, w_blocks = ad$w_blocks, a_min = ad$a_min,
             loess_span = PCFG$vmatch$span, loess_degree = PCFG$vmatch$degree, max_pool_size = PCFG$max_pool,
             null_construction = arm$nullc %||% "exact", n_boot = PCFG$n_boot,
             neighbour_cap = if (isTRUE(arm$nbr)) kc$k_max else 0L)
  } else {
    fn <- if (arm$engine == "exact") tox_compute_noise_pvalues_pipeline_exact
          else tox_compute_noise_pvalues_pipeline
    res <- fn(case_replicates = case_rep, control_replicates = ctrl_rep,
              valid_genes_own = valid, beta_centre = centre,
              norm_method = arm$norm,
              k_start = kc$k_start, k_step = kc$k_step, k_max = kc$k_max,
              tau = tau, null_method = arm$null, max_pool_size = PCFG$max_pool)
  }
  if (!is.null(res$ierr) && res$ierr != 0L) stop(sprintf("TOX ierr = %d", res$ierr))

  p <- res$pvalues_own; p[p < 0 | p > 1] <- NA
  padj <- rep(NA_real_, G); padj[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
  nc <- res$neighborhood_size_own_case; nt <- res$neighborhood_size_own_control
  # Theoretical per-gene floor where it is known in closed form.
  floor_th <- if (identical(arm$nullc, "bootstrap")) rep(1 / (PCFG$n_boot + 1), G)
              else if (arm$engine %in% c("exact", "variance"))
                ifelse(nc > 0 & nt > 0, 1 / (as.numeric(nc) * nt + 1), NA_real_)
              else if (arm$null == 0L) rep(1 / 10001, G) else rep(NA_real_, G)   # N_BOOTSTRAP_DRAWS
  out <- data.frame(gene = rownames(mat), p = p, padj = padj, lfc = res$obs_own,
                    log2fc = res$log2fc, tie = tie,
                    floor_th = floor_th, nbhd_own_case = nc, nbhd_own_control = nt,
                    lfc_scale = if (arm$norm == 0L) "linear" else "log2", row.names = NULL)
  if (arm$engine == "variance") {
    # Per-gene band diagnostics (vb_*); the scalars are repeated per row so the
    # report needs no second file. Other methods get NA (bind_methods()).
    na_neg <- function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }
    out$vb_u            <- res$match_u
    out$vb_scan_case    <- na_neg(res$n_scanned_own_case)
    out$vb_scan_radius  <- scan_radius(colMeans(case_rep), out$vb_scan_case)
    out$vb_scan_control <- na_neg(res$n_scanned_own_control)
    out$vb_genes_case    <- na_neg(nc) / nrow(case_rep)     # pool size in GENES
    out$vb_genes_control <- na_neg(nt) / nrow(ctrl_rep)
    out$vb_k   <- arm[["k"]] %||% PCFG$vmatch$k
    out$vb_tau <- res$tau_var
    out$vb_scatter   <- res$scatter
    out$vb_n_floored <- res$n_floored
    for (nm in c("lambda_pool", "lambda_case", "lambda_control", "p_hat", "participation_ratio"))
      out[[paste0("vb_", nm)]] <- res[[nm]]
    out$vb_stop_case    <- na_neg(res$stop_own_case)
    out$vb_stop_control <- na_neg(res$stop_own_control)
    out$vb_radius_case  <- na_neg(res$radius_own_case)
    out$expr_spread <- na_neg(res$expr_spread_own)
    # neighbour genes (column indices of mat), for nbhd_rows(); not written per gene
    if (!is.null(res$neighbours_own_case))
      attr(out, "nbr") <- list(case = res$neighbours_own_case, control = res$neighbours_own_control)
  }
  out
}

#' Neighbour genes of one arm's false positives and false negatives (BH < 0.05):
#' one row per gene and side, the neighbours as ";"-joined gene IDs (own gene first).
nbhd_rows <- function(r, is_de) {
  nb <- attr(r, "nbr"); if (is.null(nb)) return(NULL)
  called <- !is.na(r$padj) & r$padj < 0.05
  cls <- ifelse(called & !is_de, "FP", ifelse(!called & is_de, "FN", NA))
  i <- which(!is.na(cls)); if (!length(i)) return(NULL)
  ids <- function(m) apply(m[, i, drop = FALSE], 2, function(x) paste(r$gene[x[x > 0]], collapse = ";"))
  data.frame(gene = rep(r$gene[i], 2), class = rep(cls[i], 2), side = rep(c("case", "control"), each = length(i)),
             p = rep(r$p[i], 2), neighbours = c(ids(nb$case), ids(nb$control)))
}

#' Scan radius of the case-side band pool (note D WP4). The Fortran visits candidates
#' in increasing |linear mean - target mean| (find_closest_helper, then choose_index's
#' two-pointer expansion), so its n_scanned slots are the n_scanned genes nearest in
#' LINEAR mean. Returned: |log2(m_j + 1) - log2(m_0 + 1)| of the n_scanned-th of them
#' (the target itself is the first). NA where n_scanned is NA.
scan_radius <- function(m, n_scan) {
  s <- sort(m); out <- rep(NA_real_, length(m))
  ok <- which(is.finite(m) & is.finite(n_scan) & n_scan >= 1)
  if (!length(ok)) return(out)
  m0 <- m[ok]; n <- n_scan[ok]
  # Smallest d whose closed ball [m0 - d, m0 + d] holds >= n genes, by bisection.
  inside <- function(d) findInterval(m0 + d, s) - findInterval(m0 - d, s, left.open = TRUE)
  lo <- rep(0, length(ok)); hi <- rep(diff(range(s)), length(ok))
  for (it in 1:60) { mid <- (lo + hi) / 2; big <- inside(mid) >= n; hi[big] <- mid[big]; lo[!big] <- mid[!big] }
  # The n-th nearest is the farther (in linear distance) of the ball's two edge genes.
  r_gene <- s[pmax(findInterval(m0 + hi, s), 1L)]
  l_gene <- s[pmin(findInterval(m0 - hi, s, left.open = TRUE) + 1L, length(s))]
  mj <- ifelse(r_gene - m0 >= m0 - l_gene, r_gene, l_gene)
  out[ok] <- abs(log2(pmax(mj, 0) + 1) - log2(pmax(m0, 0) + 1))
  out
}

#' MD5 of an object (manifest data_hash: identifies the dataset every method saw).
mat_hash <- function(x) { f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, compress = FALSE)
  unname(tools::md5sum(f)) }

#' Row-bind every method's result frame into ONE results table. Methods may carry
#' extra columns (e.g. the vb_* band diagnostics); missing ones are filled with NA,
#' so adding a column to one method can never break the write again.
bind_methods <- function(results) {
  cols <- unique(unlist(lapply(results, names)))
  do.call(rbind, Map(function(m, r) {
    r[setdiff(cols, names(r))] <- NA
    cbind(method = m, r[cols])
  }, names(results), results))
}

#' Oracle pool, computed with the production Fortran.
#'
#' The Fortran builds every gene's pool from ALL genes it is handed
#' (`valid_genes_own` only selects which genes get a p-value, not who is pooled),
#' so "pool from true nulls only" is expressed by what we pass in:
#'   * null genes:  one call with the null genes only -> their pools are all-null.
#'   * each DE gene g: one call with (null genes + g), testing g alone -> g's pool
#'     is its null neighbours plus its OWN residuals, exactly as in production,
#'     where a gene is always its own nearest neighbour.
#' Relative to the production arm the only thing removed is OTHER DE genes from
#' the pool -- which is the contamination being measured.
run_tox_oracle <- function(mat, group, arm, tie, is_null) {
  if (isTRUE(arm$centre) || identical(arm$input, "tpm_centred"))
    stop("oracle arm not defined for median-centred input")
  nul <- which(is_null)
  r0 <- run_tox_arm(mat[nul, , drop = FALSE], group, arm, tie[nul])
  lg <- function(x) if (arm$norm == 0L) x else log2(pmax(x, 0) + 1)
  out <- data.frame(gene = rownames(mat), p = NA_real_, padj = NA_real_,
                    lfc = rowMeans(lg(mat[, group == "B", drop = FALSE])) -
                          rowMeans(lg(mat[, group == "A", drop = FALSE])),
                    log2fc = NA_real_, tie = tie, floor_th = NA_real_, nbhd_own_case = NA_integer_,
                    nbhd_own_control = NA_integer_, lfc_scale = r0$lfc_scale[1])
  cols <- c("p", "log2fc", "floor_th", "nbhd_own_case", "nbhd_own_control")
  out[nul, cols] <- r0[cols]
  last <- length(nul) + 1L
  for (g in which(!is_null))
    out[g, cols] <- run_tox_arm(mat[c(nul, g), , drop = FALSE], group, arm, tie[c(nul, g)],
                                test_only = last)[last, cols]
  out$padj <- NA_real_; ok <- !is.na(out$p); out$padj[ok] <- p.adjust(out$p[ok], "BH")
  out
}

#' Part Q pool diagnostics (plan section 5): tox_diagnose (R port of the exact,
#' pooled, log engine) on a random gene subset, production pool vs oracle pool.
#' Gated: validate_against_fortran must reproduce the Fortran p-values first; the
#' result is stored with gate_ok and the report uses only gate_ok rows.
run_diag_q <- function(tpm, group, is_de) {
  case <- t(tpm[, group == "B", drop = FALSE]); ctrl <- t(tpm[, group == "A", drop = FALSE])
  k <- PCFG$kcfg
  dg <- function(cm, tm, genes) tox_diagnose(cm, tm, 1L, k$k_start, k$k_step, k$k_max, PCFG$tau,
                                             PCFG$max_pool, genes = genes, verbose = FALSE)
  gate <- validate_against_fortran(case, ctrl, 1L, k$k_start, k$k_step, k$k_max, PCFG$tau,
                                   PCFG$max_pool, n_genes_check = PCFG$diag$n_check)
  sel <- sort(sample.int(ncol(case), min(PCFG$diag$n_genes, ncol(case))))
  nul <- which(!is_de)
  prod <- dg(case, ctrl, sel)
  orc <- list(dg(case[, nul, drop = FALSE], ctrl[, nul, drop = FALSE], match(sel[!is_de[sel]], nul)))
  for (g in sel[is_de[sel]]) {
    idx <- c(nul, g)
    orc[[length(orc) + 1L]] <- dg(case[, idx, drop = FALSE], ctrl[, idx, drop = FALSE], length(idx))
  }
  orc <- do.call(rbind, orc)
  out <- rbind(cbind(pool = "production", prod), cbind(pool = "oracle", orc))
  out$is_de <- is_de[match(out$gene_id, colnames(case))]
  out$gate_ok <- gate$ok; out$gate_max_abs_diff <- gate$max_abs_diff
  out
}

run_methods_power <- function(sim, part) {
  keep <- common_filter(sim)
  counts <- sim$counts[keep, , drop = FALSE]
  tpm    <- sim$tpm[keep, , drop = FALSE]
  grp    <- sim$group
  tie    <- group_lfc(tpm, grp)
  out <- list()
  safe <- function(nm, expr) {
    r <- tryCatch(expr, error = function(e) {
      message(sprintf("    %s failed: %s", nm, conditionMessage(e))); NULL })
    if (!is.null(r)) out[[nm]] <<- r
  }
  ref_tag <- function(df) data.frame(gene = df$gene, p = df$stat, padj = df$padj, lfc = df$lfc,
                                     log2fc = df$lfc, tie = tie, floor_th = NA_real_, nbhd_own_case = NA_integer_,
                                     nbhd_own_control = NA_integer_, lfc_scale = "log2")

  if (PCFG$run_edger)  safe("edgeR",  ref_tag(run_edger (counts, grp)))
  if (PCFG$run_limma)  safe("limma",  ref_tag(run_limma (counts, grp)))
  if (part == "D") {
    v <- tryCatch(run_deseq2_variants(counts, grp), error = function(e) {
      message("    DESeq2 variants failed: ", conditionMessage(e)); list() })
    for (nm in names(v)) out[[nm]] <- ref_tag(v[[nm]])
  } else if (PCFG$run_deseq2) safe("DESeq2", ref_tag(run_deseq2(counts, grp)))

  diag <- NULL
  grid <- part %in% c("N", "R") || (part == "RS" && isTRUE(PCFG$RS$grid))
  vgrid <- function(a) a$label %in% PCFG$vmatch$grid_arms
  if (RUN_TOX && grid) {
    # Input x arm grid. Every input sees the SAME genes (common_filter) and samples.
    mats <- list(tpm = tpm, cnt = counts, tmm = NULL, mor = NULL)
    mats$tmm <- tryCatch(tmm_cpm(counts), error = function(e) NULL)
    mats$mor <- tryCatch(mor_counts(counts), error = function(e) NULL)
    for (a in PCFG$tox_arms) if (a$label %in% PCFG$grid_arms)
      for (inp in PCFG$inputs) {
        if (inp$centre && a$norm == 0L) next      # centring is a log-scale shift
        mat <- mats[[inp$src]]; if (is.null(mat)) next
        ai <- modifyList(a, list(input = inp$src, centre = inp$centre))
        safe(paste0(a$label, "|", inp$name), run_tox_arm(mat, grp, ai, tie))
      }
    # variance grid arms on the other inputs (their own, TPM-centred, runs below)
    if (part != "N") for (a in Filter(vgrid, PCFG$vmatch_arms))
      for (inp in PCFG$inputs) {
        if (inp$src == "tpm" && inp$centre) next
        mat <- mats[[inp$src]]; if (is.null(mat)) next
        safe(paste0(a$label, "|", inp$name),
             run_tox_arm(mat, grp, modifyList(a, list(centre = inp$centre, nbr = TRUE)), tie))
      }
  } else if (RUN_TOX && part == "RS") {
    # Short arm list (note D WP3): Part R's comparator only, no input grid.
    for (a in PCFG$RS$tox_arms) safe(a$label, run_tox_arm(tpm, grp, a, tie))
  } else if (RUN_TOX && part != "D") {
    arms <- if (part == "C") PCFG$tox_input_arms
            else if (part == "U") Filter(function(a) a$label == "TOX-log", PCFG$tox_arms)
            else PCFG$tox_arms
    tmm <- if (any(vapply(arms, function(a) a$input == "tmm", logical(1)))) tmm_cpm(counts)
    for (a in arms) {
      mat <- if (a$input == "tmm") tmm else tpm
      safe(a$label, run_tox_arm(mat, grp, a, tie))
    }
  }

  # Variance-aware arms (notes A, F). Parts R, RS: the comparator's input, TPM-centred.
  vdiag <- NULL
  if (RUN_TOX && part %in% PCFG$vmatch$parts) {
    n_rep <- sum(grp == "B")
    for (a in PCFG$vmatch_arms) {
      if (!is.null(a$parts) && !part %in% a$parts) next
      if (part == "RS" && !is.null(a$rs_n) && !n_rep %in% a$rs_n) next
      safe(a$label, run_tox_arm(tpm, grp, modifyList(a, list(centre = part %in% REAL_PARTS, nbr = vgrid(a))), tie))
    }
    vdiag <- vdiag_frame(tpm, grp, "TPM")             # note G: method-independent
  }
  is_de <- sim$is_de[keep]
  nbhd <- do.call(rbind, lapply(names(out), function(m) {
    x <- nbhd_rows(out[[m]], is_de); if (!is.null(x)) cbind(method = m, x) }))

  if (part == "Q" && RUN_TOX) {
    is_null <- !sim$is_de[keep]
    for (a in PCFG$tox_arms) if (a$label %in% PCFG$oracle_arms) {
      mat <- if (a$input == "tmm") tmm_cpm(counts) else tpm
      safe(paste0(a$label, "-oracle"), run_tox_oracle(mat, grp, a, tie, is_null))
    }
    diag <- tryCatch(run_diag_q(tpm, grp, sim$is_de[keep]), error = function(e) {
      message("    pool diagnostics failed: ", conditionMessage(e)); NULL })
  }
  list(results = out, keep = keep, diag = diag, vdiag = vdiag, nbhd = nbhd, data_hash = mat_hash(counts))
}

# =============================================================================
# 3. JOBS
# =============================================================================

REAL_COHORTS <- list()

#' One row per dataset. Every CELL key (common/power_metrics.R) is filled for every
#' part; keys a part does not vary hold a constant (NA for delta_target / cohort).
#' Rounds per cell: the part's own `n_rounds` (Part RS), else PCFG$n_rounds. RS cells
#' with 2n > S get ONE row with skip = TRUE (recorded in the manifest, never run).
build_jobs <- function() {
  jb <- list()
  add <- function(part, g, cohort = NA_character_, skip = FALSE) {
    for (d in g$dists %||% "real") for (n in g$n_reps) for (pi1 in g$pi1)
      for (h in g$het %||% 1) for (fu in g$frac_up %||% PCFG$frac_up)
        for (dp in g$depth %||% "equal") for (dm in g$disp_model %||% "const")
          for (hvf in g$hv_frac %||% 0) for (sm in g$sigma_mult %||% 0) for (rh in g$rho %||% 1)
            for (dt in g$delta_target %||% NA_real_)
              for (i in if (skip) 1L else seq_len(g$n_rounds %||% PCFG$n_rounds))
                jb[[length(jb) + 1L]] <<- data.frame(
                  part = part, dist = d, n_rep = n, pi1 = pi1, het = h, frac_up = fu,
                  depth = dp, disp_model = dm, hv_frac = hvf,
                  sigma_mult = sm, sigma_d = sm * PCFG$disp_hetero$sigma_hat, rho = rh,
                  delta_target = dt, cohort = cohort, round = i, skip = skip)
  }
  for (pt in setdiff(ALL_PARTS, REAL_PARTS)) if (pt %in% PCFG$parts) add(pt, PCFG[[pt]])
  for (co in REAL_COHORTS) {
    S <- ncol(co$counts)
    if ("R" %in% PCFG$parts)                   # Part R: one cell per base cohort
      add("R", modifyList(PCFG$R, list(n_reps = S %/% 2L)), co$label)
    if ("RS" %in% PCFG$parts) {                # Part RS: cohort x n x pi1
      fits <- 2L * PCFG$RS$n_reps <= S
      if (any(fits))  add("RS", modifyList(PCFG$RS, list(n_reps = PCFG$RS$n_reps[fits])), co$label)
      if (any(!fits)) add("RS", modifyList(PCFG$RS, list(n_reps = PCFG$RS$n_reps[!fits])), co$label, skip = TRUE)
    }
  }
  if (!length(jb)) return(NULL)
  j <- do.call(rbind, jb)
  # One running job counter per part (for R and RS across cohorts, n, pi1 and rounds),
  # so no two jobs of a part share a seed. Same numbering as before for every part
  # that existed (each simulated part is one add() call; Part R was (cohort, round)
  # in this order). Seed depends on (part, index within part) only, so running a
  # subset of parts reproduces the same datasets as the full run; across parts the
  # seeds cannot collide (7919 and 1000003 are coprime, job < 1000003).
  j$job <- ave(seq_len(nrow(j)), j$part, FUN = seq_along)
  j$seed <- PCFG$seed + 7919L * j$job + 1000003L * match(j$part, SEED_PARTS)
  j$job_id <- paste0(gsub("[^A-Za-z]", "plus", j$part), "_", j$job)
  j
}

#' The dispersion-parameter fields recorded in the manifest and the log header.
disp_fields <- function(d = PCFG$disp_hetero)
  data.frame(disp_source = d$source, disp_date = d$date %||% NA_character_,
             disp_phi_inf = d$phi_inf, disp_c = d$c, disp_sigma_hat = d$sigma_hat,
             disp_mu_depth = d$mu_depth %||% NA_real_)

run_job <- function(job) {
  if (isTRUE(job$skip)) return(cbind(job, disp_fields()))
  t0 <- proc.time()[["elapsed"]]
  set.seed(job$seed)
  sim <- if (job$part %in% REAL_PARTS) {
    simulate_real(REAL_COHORTS[[job$cohort]], job$pi1, n = if (job$part == "RS") job$n_rep)
  } else {
    mode <- if (job$part %in% c("C", "N")) "absolute" else "relative"
    dt <- if (is.na(job$delta_target)) NULL else job$delta_target
    tr <- make_truth_power(PCFG$n_genes, job$pi1, job$n_rep, job$het, mode, job$frac_up, dt)
    hv <- if (job$hv_frac > 0) {
      nul <- which(!tr$is_de)
      seq_along(tr$is_de) %in% nul[sample.int(length(nul), round(job$hv_frac * length(nul)))]
    }
    c(simulate_power(tr, job$dist, job$n_rep, depth = job$depth,
                     disp_model = job$disp_model, hv = hv, sigma_d = job$sigma_d, rho = job$rho),
      tr[c("eff_lfc_median_abs", "frac_sign_flip", "frac_below_floor")])
  }

  rr <- run_methods_power(sim, job$part)
  dir <- file.path(PCFG$raw_dir, job$part)
  f <- function(kind) file.path(dir, sprintf("%s_%s.csv.gz", kind, job$job_id))
  truth <- data.frame(gene = rownames(sim$counts), is_de = sim$is_de, true_lfc = sim$true_lfc,
                      expr = sim$expr, kept = rr$keep, is_hv = sim$is_hv)
  if (!is.null(sim$realised_lfc)) truth$realised_lfc <- sim$realised_lfc
  write_gz(truth, f("truth"))
  if (length(rr$results))
    write_gz(bind_methods(rr$results), f("results"))
  if (!is.null(rr$diag)) write_gz(rr$diag, f("diag"))
  if (!is.null(rr$vdiag)) write_gz(rr$vdiag, f("vdiag"))
  if (!is.null(rr$nbhd)) write_gz(rr$nbhd, f("nbhd"))
  raw_bytes <- sum(file.size(Filter(file.exists, f(c("truth", "results", "diag", "vdiag", "nbhd")))))

  num <- function(x) if (is.null(x)) NA_real_ else x
  # Scalars of the reference variance arms (note F): the band width of the band as
  # it is now, and the James-Stein shrinkage of the primary cross-fit arm (the
  # lambdas depend on the data and the target only, not on the arm).
  vb <- function(lab, col) { r <- rr$results[[lab]]
    if (is.null(r) || is.null(r[[col]])) NA_real_ else r[[col]][1] }
  real <- job$part %in% REAL_PARTS; co <- if (real) REAL_COHORTS[[job$cohort]]
  cbind(job,
        vb_tau_band_c3 = vb("TOX-log-vband-c3", "vb_tau"),
        vb_lambda_pool = vb("TOX-log-xfit-js", "vb_lambda_pool"),
        vb_lambda_case = vb("TOX-log-xfit-js", "vb_lambda_case"),
        vb_lambda_control = vb("TOX-log-xfit-js", "vb_lambda_control"),
        vb_lambda_pool_loess = vb("TOX-log-xfit-js-loess", "vb_lambda_pool"),
        vb_p_hat = vb("TOX-log-xfit-js", "vb_p_hat"),
        vb_participation_ratio = vb("TOX-log-xfit-js", "vb_participation_ratio"),
        vb_n_floored = vb("TOX-log-xfit-js", "vb_n_floored"),
        vb_span = PCFG$vmatch$span, vb_k = PCFG$vmatch$k,
        delta = num(sim$delta), delta_realised = num(sim$delta_realised),
        n_de_total = sum(sim$is_de), n_kept = sum(rr$keep), n_methods = length(rr$results),
        eff_lfc_median_abs = num(sim$eff_lfc_median_abs), frac_sign_flip = num(sim$frac_sign_flip),
        frac_below_floor = num(sim$frac_below_floor),
        S = num(sim$S), n_per_group = num(sim$n_per_group), n_dropped = num(sim$n_dropped),
        split_seed = if (real) job$seed else NA_integer_,
        frac_lfc_mismatch = num(sim$frac_lfc_mismatch),
        cohort_file = if (real) co$path else NA_character_,
        patient_check = if (real) co$patient_check else NA_character_,
        n_dup_removed = if (real) num(co$n_dup_removed) else NA_real_,
        frac_len_bad = if (real) num(co$frac_len_bad) else NA_real_,
        zero_free = if (real) isTRUE(PCFG$R$zero_free) else NA,
        n_zero_dropped = if (real) num(co$n_zero_dropped) else NA_real_,
        disp_fields(),
        data_hash = rr$data_hash,
        has_diag = !is.null(rr$diag),
        elapsed_sec = proc.time()[["elapsed"]] - t0, raw_bytes = raw_bytes)
}

#' Commit of the scripts (and of Tensor-Omics, whose Fortran produced the TOX p-values).
#' The server copy is rsynced without .git: set POWER_GIT_COMMIT there.
git_commit <- function() {
  env <- Sys.getenv("POWER_GIT_COMMIT", "")
  if (nzchar(env)) return(env)
  g <- function(d) tryCatch(suppressWarnings(system2("git", c("-C", shQuote(d), "rev-parse", "--short", "HEAD"),
                                                     stdout = TRUE, stderr = FALSE))[1],
                            error = function(e) NA_character_)
  paste0("tests:", g(COMMON_DIR) %||% NA, " tox:", g(".") %||% NA)
}

#' Expected p-value floor per engine, and the BH minimum discovery set m* it implies.
#'   exact:            pool = k_max * n_rep residuals,  p_min = 1/(pool^2 + 1)
#'   bootstrap pooled: p_min = 1/(N_BOOTSTRAP_DRAWS + 1) = 1/10001
#'   blocked, enumerated (n_rep <= 5): W = k_max * n_rep^n_rep, p_min = 1/(W^2 + 1)
#' G is the gene count BEFORE filtering, so m* here is an upper bound on the realised one.
resolution_report_power <- function(n_reps, G = PCFG$n_genes, alpha = 0.05) {
  km <- PCFG$kcfg$k_max
  cat("\n--- p-value floor and BH minimum discovery set m* = ceil(G p_min / alpha) ---\n")
  cat(sprintf("%6s %-18s %11s %8s\n", "n_rep", "engine", "p_min", "m*"))
  for (n in sort(n_reps)) {
    fl <- c(exact = 1 / ((km * n)^2 + 1), `bootstrap pooled` = 1 / 10001,
            `blocked enum` = if (n <= 5) 1 / ((km * n^n)^2 + 1) else 1 / 10001)
    for (e in names(fl))
      cat(sprintf("%6d %-18s %11.2e %8d\n", n, e, fl[[e]], as.integer(ceiling(G * fl[[e]] / alpha))))
  }
  cat("  (blocked above n_rep = 5 falls back to sampling at the bootstrap floor)\n\n")
}

#' Log header (note C section 3.3, note D WP2.3): one screen on what this run is.
run_header <- function(jobs, commit) {
  d <- PCFG$disp_hetero
  cat(strrep("=", 100), "\npower_test.R  ", format(Sys.time()), "  commit ", commit, "\n", sep = "")
  cat(sprintf("parts: %s | datasets: %d (%d skipped cells) | cores: %d | raw -> %s\n",
              paste(unique(jobs$part), collapse = ","), sum(!jobs$skip), sum(jobs$skip), PCFG$n_cores, PCFG$raw_dir))
  rr <- aggregate(round ~ part, jobs[!jobs$skip, ], max)
  cat("rounds per cell: ", paste(rr$part, rr$round, sep = "=", collapse = ", "), "\n", sep = "")
  cat(sprintf("dispersion parameters (P+): %s | date %s | phi_inf %s | c %s | sigma_hat %s | mu_depth %s\n",
              d$source, d$date %||% "-", format(d$phi_inf), format(d$c), format(d$sigma_hat), format(d$mu_depth %||% NA)))
  if (length(REAL_COHORTS))
    cat(sprintf("real cohorts: %d (%s)%s\n", length(REAL_COHORTS), paste(names(REAL_COHORTS), collapse = ", "),
                if (isTRUE(PCFG$R$zero_free)) " | E1 zero-free universe ON" else ""))
  cat(strrep("=", 100), "\n", sep = "")
}

#' Runtime and storage per part, extrapolated to a full run (note D WP3: the smoke
#' run must print them). Full run = POWER_EXTRAP_JOBS datasets (default: this run's
#' cells at the part's configured rounds) on POWER_EXTRAP_CORES cores (default 32).
cost_report <- function(man) {
  man <- man[!man$skip %in% TRUE, ]
  cores <- as.numeric(Sys.getenv("POWER_EXTRAP_CORES", "32"))
  x <- do.call(rbind, lapply(split(man, man$part), function(d) {
    cells <- nrow(unique(d[c("dist", "n_rep", "pi1", "het", "frac_up", "depth", "disp_model", "hv_frac",
                             "sigma_d", "delta_target", "cohort")]))
    full <- as.numeric(Sys.getenv("POWER_EXTRAP_JOBS",
                                  cells * (PCFG[[d$part[1]]]$n_rounds %||% PCFG$n_rounds)))
    data.frame(part = d$part[1], jobs = nrow(d), methods_mean = mean(d$n_methods),
               sec_per_job = mean(d$elapsed_sec), sec_per_job_method = mean(d$elapsed_sec / pmax(d$n_methods, 1)),
               mb_per_job = mean(d$raw_bytes) / 2^20, full_jobs = full,
               full_hours = full * mean(d$elapsed_sec) / cores / 3600,
               full_gb = full * mean(d$raw_bytes) / 2^30)
  }))
  print_tab(x, sprintf("RUNTIME AND STORAGE per part, extrapolated to full_jobs datasets on %g cores", cores),
            notes = "set POWER_EXTRAP_JOBS / POWER_EXTRAP_CORES to extrapolate to another run size", keys = "part")
}

# =============================================================================
#' The whole run: load cohorts, build and run the jobs, write the manifests. A
#' function so knn_config_test.R can call it with its own PCFG.
power_main <- function() {
  if (any(REAL_PARTS %in% PCFG$parts)) {
    rds <- trimws(strsplit(PCFG$R$input_rds, ",")[[1]])
    if (!length(rds)) message("Parts R / RS skipped: set POWER_REAL_RDS (comma-separated cohort RDS files).")
    for (p in rds) { co <- load_real_cohort(p); REAL_COHORTS[[co$label]] <<- co }
  }
  jobs <- build_jobs()
  if (is.null(jobs) || !nrow(jobs)) stop("No jobs -- check POWER_PARTS / POWER_REAL_RDS.")
  commit <- git_commit()
  run_header(jobs, commit)
  if (any(jobs$skip))
    message("RS cells skipped (2n > S): ", paste(unique(paste0(jobs$cohort, " n=", jobs$n_rep)[jobs$skip]), collapse = "; "))
  if (PCFG$disp_hetero$source == "placeholder" && "P+" %in% jobs$part)
    message("NOTE: P+ uses PLACEHOLDER dispersion parameters; set POWER_DISP_PARAMS ",
            "(estimate_dispersion_trend.R) for a reportable run.")
  if (RUN_TOX) resolution_report_power(unique(jobs$n_rep))

  # A part's directory holds exactly one run: clear it, so stale jobs never mix in.
  for (pt in unique(jobs$part)) {
    d <- file.path(PCFG$raw_dir, pt); unlink(d, recursive = TRUE)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  res <- parallel::mclapply(split(jobs, seq_len(nrow(jobs))), function(j)
           tryCatch(run_job(j), error = function(e) {
             message("job ", j$job_id, " failed: ", conditionMessage(e)); NULL }),
         mc.cores = PCFG$n_cores, mc.preschedule = FALSE)
  res <- Filter(is.data.frame, res)
  if (!length(res)) stop("Every job failed.")
  man <- rbind_fill(res); man$git_commit <- commit
  for (pt in unique(man$part)) {
    d <- file.path(PCFG$raw_dir, pt)
    write.csv(man[man$part == pt, ], file.path(d, "manifest.csv"), row.names = FALSE)
    writeLines(c(paste("git_commit:", commit), paste("date:", format(Sys.time())),
                 capture.output(dput(PCFG))), file.path(d, "pcfg_dump.txt"))
  }
  cost_report(man)
  message(sprintf("\n%d / %d jobs written to %s/. Next: Rscript power_report.R",
                  sum(!man$skip), sum(!jobs$skip), PCFG$raw_dir))
}

if (sys.nframe() == 0L) power_main()
