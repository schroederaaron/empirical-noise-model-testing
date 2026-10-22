# Implementation note C — readable console output (power_report.R, power_test.R, Slurm logs)

**For:** Claude Code. Independent of notes A and B; can be done first. No statistics change; output formatting only.
**Evidence:** Slurm log `power_full_11874.out` (20,783 lines, all of it unreadable as a result report). Written without executing R: the helper below is a sketch, harden it with the tests in section 5.

---

## 1. What is wrong (measured on the log)

1. **Every wide table is split into 3–5 column blocks, and each block repeats ALL rows.** Blocks of one table sit hundreds to thousands of lines apart
   (e.g. a 621-row table: block headers at lines 2361 / 2983 / 3605 / 4227 / 4849; a 1,770-row table: 13289 / 15060 / 16831 / 18602).
   R printed with `print(x, row.names = FALSE)` (`power_report.R`, `print_tab()` and the VARIANCE-MATCHED ARMS block), so **nothing links a row in a later block to its keys**.
   The last lines of the log are single columns of numbers with no identifier at all.
2. **Key columns are split across blocks.** `method` (and `cohort`, `delta_target`) often appear only in the SECOND block, so even the identity of a row is separated from its cell.
3. **Titles and notes are missing from the file.** `print_tab()` writes titles with `message()` (stderr) and the table with `print()` (stdout); the sbatch scripts send them to two files
   (`#SBATCH --output=...out`, `#SBATCH --error=...err`). The `.out` therefore has bare tables with no heading telling which table is which.
4. **Noise:** all-NA columns printed hundreds of times (`d0_mean`, `d0_inf_mean`, `w_mean` for the A arms), constant columns (`depth`, `disp_model`, `hv_frac`, `frac_up`) in every row,
   and tables with 1,770 rows dumped to the console.
5. The current width is already about 144 characters (longest line in the log), so raising it alone cannot work: these tables need 300-450 characters.

---

## 2. Rules (apply to ALL console output of the benchmark scripts)

* **R1 — width:** `options(width = 250, max.print = 1e6)` set once in `common/config.R` (sourced by every script), with a comment saying why. Also set it at the top of `power_report.R`.
* **R2 — one print helper for tables:** every table printed to the console goes through ONE function `print_tab()`; no direct `print(df, row.names = FALSE)` anywhere (add a grep check, section 5).
* **R3 — keys are never separated from their values:** if a table is wider than `width`, split it into column chunks and **repeat the key columns in every chunk**.
* **R4 — titles and tables go to the same stream** (stdout, via `cat`), so that they stay together in one file.
* **R5 — the console is a summary, the CSVs are the data:** console tables are capped (default 60 rows, argument `n`) and say where the full table is; constant columns are moved into the title; all-NA columns are dropped.
* **R6 — methods side by side:** headline comparisons are printed as a pivot (rows = cell, columns = methods, one metric per table), not one row per method.
* **R7 — one Slurm log file:** remove the `#SBATCH --error` line from `power_smoke.sbatch` and `power_full.sbatch`; without it Slurm writes stderr into the same file as stdout. (Keep `--output`.)

---

## 3. Implementation spec

### 3.1 `print_tab()` (replace the existing one in `power_report.R`; move to `common/power_metrics.R` or a new `common/print_utils.R` if it is also needed by `power_test.R`)

Behaviour:

1. Drop all-NA columns. Detect columns with a single distinct value; remove them from the table and list them in one line under the title (`constant: depth=equal, hv_frac=0, ...`). `method` is never removed.
2. Key columns = those of `c("part","dist","n_rep","pi1","het","frac_up","sigma_d","delta_target","cohort","method")` that are present and not constant, in that order, always first.
3. Sort rows by the keys (methods of one cell adjacent). `signif(., 4)` on numeric columns.
4. Print at most `n` rows (default 60); if cut, print `(showing 60 of N rows; full table: <csv name>)`. The csv name is an argument.
5. Chunk the non-key columns so that `keys + chunk` fits in `getOption("width")` and print each chunk with the keys; separate chunks by a short dashed line.
6. Title and notes via `cat()` to stdout.

Sketch (untested):

```r
print_tab <- function(x, title, csv = NULL, n = 60L, notes = NULL,
                      keys = c("part","dist","n_rep","pi1","het","frac_up","sigma_d","delta_target","cohort","method"),
                      width = getOption("width")) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- x[, vapply(x, function(v) !all(is.na(v)), logical(1)), drop = FALSE]
  const <- setdiff(names(x)[vapply(x, function(v) length(unique(v[!is.na(v)])) <= 1L, logical(1))], "method")
  keys  <- setdiff(intersect(keys, names(x)), const)
  num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], signif, 4)
  if (length(keys)) x <- x[do.call(order, unname(x[keys])), , drop = FALSE]
  total <- nrow(x); shown <- head(x, n); vals <- setdiff(names(shown), c(keys, const))
  cat(strrep("=", min(width, 100)), "\n", title, "\n", sep = "")
  if (length(const)) cat("constant: ", paste(const, vapply(x[const], function(v) as.character(v[!is.na(v)][1]), ""),
                                              sep = "=", collapse = ", "), "\n", sep = "")
  for (nt in notes) cat("  ", nt, "\n", sep = "")
  w   <- function(df, nm) pmax(vapply(df, function(v) max(nchar(format(v)), 0L), 0L), nchar(nm)) + 1L
  kw  <- sum(w(shown[keys], keys)); vw <- w(shown[vals], vals)
  cid <- (cumsum(vw) - 1L) %/% max(width - kw - max(vw), 20L)             # chunk id per value column
  for (k in unique(cid)) { print(shown[, c(keys, vals[cid == k]), drop = FALSE], row.names = FALSE)
                           if (k != max(cid)) cat(strrep("-", 20), "\n") }
  if (total > n) cat(sprintf("(showing %d of %d rows%s)\n", n, total, if (length(csv)) paste0("; full table: ", csv) else ""))
  invisible(x)
}
```

### 3.2 `print_pivot()`

`print_pivot(S, metric, rows = c("dist","n_rep","sigma_d"), methods = NULL, title, digits = 3)`: rows = cell, columns = methods (order fixed by the report's palette order; at most about 12 methods per table),
values = `metric` (`signif`/`round` to `digits`). Use it for the headline summary of each part: `tpr_nom_05_mean`, `fdr_nom_05_mean`, `tpr_ach05_mean`, `auprc_excess_pe_mean`
(and for the new arms the pool diagnostics `frac_short_pool_mean`, `n_scanned_median_mean`). Rows that belong to different `part` values are printed as separate tables.

### 3.3 Replace the raw prints in `power_report.R`

* `print(vs, row.names = FALSE)` (VARIANCE-MATCHED ARMS block) -> `print_tab(vs, "VARIANCE-MATCHED ARMS: pool diagnostics", csv = "power_summary.csv")`; the message text becomes a `notes` argument.
* The existing `report()` section prints: route all through `print_tab()` / `print_pivot()`, and drop the full per-round/per-cell dumps (they stay in the CSVs).
* Add a one-screen header at the top of the log: parts run, `POWER_ROUNDS`, `disp_source` and its parameters (phi_inf, c, sigma_hat), commit hash, number of cohorts.
* Add `message()` -> `cat()` for the section banners that are meant to be read with the tables.

### 3.4 Slurm

Remove `#SBATCH --error=logs/%x_%j.err` from `power_smoke.sbatch` and `power_full.sbatch`. Keep `--output`.

---

## 4. Out of scope

Plots, CSV contents and column names, metric definitions. (CSV files stay the data source for analysis.)

## 5. Acceptance

1. Run `power_report.R` on an existing `power_raw/`; the single log file has **every table preceded by its title**, no table line longer than `width`, and **no table chunk without its key columns** (test: capture the output of `print_tab()` on a 40-column data frame with `width = 100`; every chunk must contain all key columns).
2. The log is at most about 1,500 lines (it is 20,783 today), and each CSV named in a "full table:" hint exists in `power_out/`.
3. `grep -n "print(.*row.names *= *FALSE" TCGA_test/scripts/*.R common/*.R` returns nothing except inside `print_tab()`.
4. A pivot table for one part shows all its methods side by side on one screen at width 250.
5. No all-NA column and no constant column is printed (the constants appear once in the title line).
