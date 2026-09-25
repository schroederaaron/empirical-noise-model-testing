# config_gof.R -- every knob of the GOF analysis (issue #188), in one place.
#
# Each value can be overridden by the environment variable GOF_<NAME>, and every
# command-line argument `--NAME=value` is mapped onto GOF_<NAME> before the config
# is read (`--dry-run` sets DRY_RUN). The resolved config is printed by every script
# and written to each output directory's provenance.
#
# The decisions of brief section 12 (D1-D6) are exposed here, each with the
# brief's proposed value as the default. None is applied silently: all of them are
# printed with the config and recorded in provenance, and the modes that change the
# data (ROUNDING, OFFSET_MODE) are written into every output file name.

.gof_apply_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  for (x in args) {
    if (identical(x, "--dry-run")) { Sys.setenv(GOF_DRY_RUN = "TRUE"); next }
    m <- regmatches(x, regexec("^--([A-Za-z0-9_]+)=(.*)$", x))[[1]]
    if (length(m) != 3L) stop("Unrecognised argument '", x, "'. Use --NAME=value or --dry-run.")
    do.call(Sys.setenv, stats::setNames(list(m[3]), paste0("GOF_", toupper(m[2]))))
  }
}

.cfg <- function(name, default, choices = NULL, vector = FALSE) {
  raw <- Sys.getenv(paste0("GOF_", name), NA_character_)
  v <- if (is.na(raw)) default
       else if (vector) { s <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]]); s[nzchar(s)] }
       else if (is.logical(default)) as.logical(raw)
       else if (is.integer(default)) suppressWarnings(as.integer(raw))
       else if (is.numeric(default)) suppressWarnings(as.numeric(raw))
       else raw
  if (!vector && length(v) == 1L && is.na(v) && !is.na(raw))
    stop("GOF_", name, "='", raw, "' is not a valid ", class(default)[1], ".")
  if (!is.null(choices) && !all(v %in% choices))
    stop("GOF_", name, " must be ", if (vector) "a subset of " else "one of ",
         paste(choices, collapse = " / "), "; got '", paste(v, collapse = ","), "'.")
  v
}

.gof_apply_cli()

ALL_MODELS <- c("poisson", "nb", "zinb", "genpois", "pln", "cmp")

GOF_CONFIG <- list(
  # ---- input / output ----
  DATASET          = .cfg("DATASET", ""),              # path to a dataset .rds (brief section 3)
  OUT_ROOT         = .cfg("OUT_ROOT", file.path(GOF_ROOT, "results")),
  DRY_RUN          = .cfg("DRY_RUN", FALSE),
  RUN_DIR          = .cfg("RUN_DIR", ""),              # an 02_run_gof.R output dir (03_report.R input)

  # ---- decisions (brief section 12) ----
  ROUNDING         = .cfg("ROUNDING", "error", c("error", "stochastic", "floor")),         # D1
  OFFSET_MODE      = .cfg("OFFSET_MODE", "tmm", c("tmm", "tmm_length")),                   # D2
  OFFSETS_BOOT     = .cfg("OFFSETS_BOOT", "reestimate", c("reestimate", "fixed")),         # D3
  N_GENES_BOOT     = .cfg("N_GENES_BOOT", 2000L),                                          # D4
  RANK_BY          = .cfg("RANK_BY", "heldout", c("heldout", "excess_W2")),                # D5
  EXCLUDE_SAMPLES  = .cfg("EXCLUDE_SAMPLES", character(0), vector = TRUE),                 # D6

  # ---- models ----
  MODELS           = .cfg("MODELS", ALL_MODELS, ALL_MODELS, vector = TRUE),

  # ---- seeds (every stochastic step derives its own stream from BASE_SEED) ----
  BASE_SEED        = .cfg("BASE_SEED", 188L),

  # ---- sizes ----
  MAX_GENES_OBS    = .cfg("MAX_GENES_OBS", 0L),        # 0 = all filtered genes; >0 = dev cap (M2: 500)
  B                = .cfg("B", 100L),                  # bootstrap replicates: B_DEV = 100, B_FINAL = 500
  RERAND_R         = .cfg("RERAND_R", 20L),            # PIT re-randomisations
  K_FOLDS          = .cfg("K_FOLDS", 8L),              # held-out folds (leave-one-out below n = 16)
  N_CORES          = .cfg("N_CORES", max(1L, parallel::detectCores() - 1L)),

  # ---- numerics ----
  PLN_NODES        = .cfg("PLN_NODES", 60L),           # adaptive Gauss-Hermite nodes for the PLN pmf. The brief's
                                                       # default is 30; at sigma = 1, mu = 0.5 that leaves a 2.5e-10
                                                       # relative pmf error (vs a trapezoid reference) and fails the
                                                       # G1 increment check; 60 passes the whole G1 grid.
  PLN_CDF_PANELS   = .cfg("PLN_CDF_PANELS", 12L),       # graded Gauss-Legendre panels per side (PLN cdf)
  PLN_CDF_GL       = .cfg("PLN_CDF_GL", 16L),          # Gauss-Legendre nodes per panel
  CMP_K            = .cfg("CMP_K", 40),                # CMP normalising window half-width, in sd units
  CMP_MAX_WINDOW   = .cfg("CMP_MAX_WINDOW", 5e6),      # refuse a CMP window wider than this (fit fails, recorded)

  # ---- stages (for development milestones) ----
  RUN_BOOTSTRAP    = .cfg("RUN_BOOTSTRAP", TRUE),
  RUN_HELDOUT      = .cfg("RUN_HELDOUT", TRUE),
  RUN_RERAND       = .cfg("RUN_RERAND", TRUE),

  # ---- validation (00_validate_families.R, 01_calibration_sim.R) ----
  N_TARGET         = .cfg("N_TARGET", 96L),            # n for G4/G5 when no DATASET is given
  G5_REF_DIR       = .cfg("G5_REF_DIR", ""),           # an 02_run_gof.R output dir: source of NB/PLN parameters
  G5_M             = .cfg("G5_M", 40L),
  G5_G             = .cfg("G5_G", 300L),
  G5_B             = .cfg("G5_B", 99L)
)

if (GOF_CONFIG$N_GENES_BOOT < 1L || GOF_CONFIG$B < 1L || GOF_CONFIG$N_CORES < 1L)
  stop("N_GENES_BOOT, B and N_CORES must be >= 1.")

print_gof_config <- function(cfg = GOF_CONFIG) {
  cat("---- resolved config (GOF_CONFIG) ----\n")
  for (k in names(cfg)) {
    v <- cfg[[k]]
    cat(sprintf("  %-17s = %s\n", k,
                if (length(v) == 0L) "<none>" else paste(as.character(v), collapse = ", ")))
  }
  cat("--------------------------------------\n")
  invisible(cfg)
}

# The data-changing modes go into every output file name (brief section 4.2).
gof_mode_tag <- function(cfg = GOF_CONFIG)
  sprintf("round-%s_off-%s", cfg$ROUNDING, cfg$OFFSET_MODE)
