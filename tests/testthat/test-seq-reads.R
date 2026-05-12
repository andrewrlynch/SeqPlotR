library(testthat)
library(GenomicRanges)

# ---- Constructor guardrails ----

test_that("seq_reads() errors when window exceeds max_width", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  wide_win <- GRanges("chr1", IRanges(1, 200001))
  expect_error(seq_reads("fake.bam", wide_win, max_width = 200000),
               "max_width")
})

test_that("seq_reads() errors when region is not GRanges", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  expect_error(seq_reads("fake.bam", "chr1:1-1000"), "GRanges")
})

test_that("seq_reads() errors when BAM not found", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  win <- GRanges("chr1", IRanges(1, 1000))
  expect_error(seq_reads("/no/such/file.bam", win), "not found")
})

# ---- BAM-dependent tests (skip when no test BAM is configured) ----

bam_path <- Sys.getenv("SEQPLOTR_TEST_BAM", unset = "")
bam_seq  <- Sys.getenv("SEQPLOTR_TEST_BAM_SEQ", unset = "chr20")
bam_pos  <- as.integer(Sys.getenv("SEQPLOTR_TEST_BAM_POS",
                                  unset = "30000000"))

test_that("sort_by defaults to 'insert_length'", {
  skip_if(nchar(bam_path) == 0, "No test BAM configured")
  win <- GRanges(bam_seq, IRanges(bam_pos, bam_pos + 50000L))
  r <- seq_reads(bam_path, win, max_reads = 500L)
  expect_equal(r$sort_by, "insert_length")
})

test_that("sort_by = 'start' accepted", {
  skip_if(nchar(bam_path) == 0, "No test BAM configured")
  win <- GRanges(bam_seq, IRanges(bam_pos, bam_pos + 50000L))
  r <- seq_reads(bam_path, win, sort_by = "start", max_reads = 500L)
  expect_equal(r$sort_by, "start")
})

test_that("link_mates defaults to TRUE", {
  skip_if(nchar(bam_path) == 0, "No test BAM configured")
  win <- GRanges(bam_seq, IRanges(bam_pos, bam_pos + 50000L))
  r <- seq_reads(bam_path, win, max_reads = 500L)
  expect_true(r$link_mates)
})

test_that("seq_reads integrates into seq_plot$plot()", {
  skip_if(nchar(bam_path) == 0, "No test BAM configured")
  skip_if_not(capabilities("png"))
  win <- GRanges(bam_seq, IRanges(bam_pos, bam_pos + 50000L))
  reads <- seq_reads(bam_path, win, max_reads = 500L)
  sp <- seq_plot() %+%
    (seq_track(windows = win, height = 4) %+% reads)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})
