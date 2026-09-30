write_gbgout <- function(dir, name, feat, counts) {
  summ <- data.frame(V1 = HTSEQ_SUMMARY_ROWS, V2 = c(10, 5, 3, 2, 1))
  utils::write.table(rbind(data.frame(V1 = feat, V2 = counts), summ), file.path(dir, paste0(name, ".gbgout")),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

htseq_fixture <- function() {
  dir <- tempfile("htseq"); dir.create(dir)
  feat <- c("YAL001C", "YAL002W", "YAL003W", "YAL004W")
  write_gbgout(dir, "WT_rep01_MID39_allLanes_tophat2.0.5.bam", feat, c(0, 3, 100, 7))
  write_gbgout(dir, "WT_rep02_MID40_allLanes_tophat2.0.5.bam", feat, c(1, 4, 120, 9))
  write_gbgout(dir, "Snf2_rep06_MID80_allLanes_tophat2.0.5.bam", feat, c(2, 0, 90, 5))
  excl <- file.path(dir, "exclude.lst")
  writeLines("Snf2_rep06_MID80_allLanes_tophat2.0.5.bam", excl)
  list(dir = dir, excl = excl, feat = feat)
}

test_that("htseq builder drops summary rows, names samples and flags bad replicates", {
  fx <- htseq_fixture()
  ds <- build_from_htseq_dir(fx$dir, fx$excl, "fx_all", drop_bad = FALSE)
  expect_identical(dim(ds$counts), c(4L, 3L))
  expect_identical(rownames(ds$counts), fx$feat)
  expect_false(any(HTSEQ_SUMMARY_ROWS %in% rownames(ds$counts)))
  expect_setequal(colnames(ds$counts), c("WT_rep01", "WT_rep02", "Snf2_rep06"))
  expect_identical(ds$source, "htseq_counts")
  expect_identical(levels(ds$samples$condition), c("WT", "Snf2"))
  expect_identical(ds$samples$bad_replicate[ds$samples$sample_id == "Snf2_rep06"], TRUE)
  expect_equal(ds$samples$lib_size, unname(colSums(ds$counts)))
  expect_identical(ds$samples$replicate[ds$samples$sample_id == "Snf2_rep06"], 6L)

  clean <- build_from_htseq_dir(fx$dir, fx$excl, "fx_clean", drop_bad = TRUE)
  expect_identical(colnames(clean$counts), c("WT_rep01", "WT_rep02"))
  expect_false(any(clean$samples$bad_replicate))
})

test_that("htseq builder rejects differing feature lists and unmatched exclude entries", {
  fx <- htseq_fixture()
  write_gbgout(fx$dir, "WT_rep03_MID41_allLanes_tophat2.0.5.bam", rev(fx$feat), c(1, 1, 1, 1))
  expect_error(build_from_htseq_dir(fx$dir, fx$excl, "fx", FALSE), "feature list")

  fx2 <- htseq_fixture()
  writeLines("WT_rep99_MID1_allLanes_tophat2.0.5.bam", fx2$excl)
  expect_error(build_from_htseq_dir(fx2$dir, fx2$excl, "fx", FALSE), "exactly one")
})
