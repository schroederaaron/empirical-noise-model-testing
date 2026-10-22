# Unit tests for common/print_utils.R (implementation note C section 5).
library(testthat)
source("../../common/print_utils.R")

wide <- function() {
  set.seed(31)
  x <- data.frame(part = "P+", dist = rep(c("nb", "lnpois"), each = 6), n_rep = rep(c(3L, 5L, 10L), 4),
                  method = rep(c("limma", "TOX-log"), 6), depth = "equal", all_na = NA_real_)
  for (i in 1:40) x[[sprintf("metric_number_%02d_mean", i)]] <- runif(12)
  x
}

test_that("every column chunk carries all key columns, no line exceeds the width", {
  out <- capture.output(print_tab(wide(), "T", width = 100))
  hdr <- grep("metric_number", out, value = TRUE)
  expect_gt(length(hdr), 3)                                         # really split into chunks
  for (h in hdr) expect_true(all(c("dist", "n_rep", "method") %in% strsplit(trimws(h), " +")[[1]]))
  body <- out[-(1:3)]
  expect_true(all(nchar(body) <= 100))
})

test_that("title on stdout, constants in one line, all-NA columns dropped, rows capped with CSV hint", {
  out <- capture.output(print_tab(wide(), "MY TITLE", csv = "x.csv", n = 5, width = 250))
  expect_true("MY TITLE" %in% out)
  expect_true("constant: part=P+, depth=equal" %in% out)
  expect_false(any(grepl("all_na|depth", out[!grepl("^constant", out)])))
  expect_true(any(grepl("showing 5 of 12 rows; full table: x.csv", out, fixed = TRUE)))
})

test_that("print_pivot: methods side by side, values in the right cells", {
  S <- data.frame(part = "RS", cohort = rep(c("A", "B"), each = 3), n_rep = 3L,
                  method = rep(c("limma", "edgeR", "TOX"), 2), v_mean = c(1, 2, 3, 4, 5, 6))
  out <- capture.output(pv <- print_pivot(S, "v_mean", "P", methods = c("TOX", "limma", "edgeR")))
  expect_identical(names(pv), c("cohort", "TOX", "limma", "edgeR"))
  expect_identical(pv$TOX, c(3, 6)); expect_identical(pv$limma, c(1, 4))
  expect_true(any(grepl("constant: part=RS, n_rep=3", out, fixed = TRUE)))
  out <- capture.output(print_pivot(S, "v_mean", "P", max_methods = 2L))
  expect_identical(sum(grepl("^ *cohort", out)), 2L)                 # 3 methods -> 2 blocks, keys in both
})
