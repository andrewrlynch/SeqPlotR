library(testthat)
library(GenomicRanges)

# ---- open_bigwig ----

test_that("open_bigwig() errors on missing file", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  expect_error(open_bigwig("/no/such.bw"), "not found")
})

test_that("open_bigwig() returns SeqBigWig object", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  bw <- Sys.getenv("SEQPLOTR_TEST_BW", unset = "")
  skip_if(nchar(bw) == 0, "No test bigWig configured")
  obj <- open_bigwig(bw)
  expect_s3_class(obj, "SeqBigWig")
  expect_true(is_seq_file_conn(obj))
})

test_that("open_bigwig() $fetch() errors when region exceeds max_fetch_bp", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  bw <- Sys.getenv("SEQPLOTR_TEST_BW", unset = "")
  skip_if(nchar(bw) == 0, "No test bigWig configured")
  obj  <- open_bigwig(bw, max_fetch_bp = 1000L)
  wide <- GRanges("chr1", IRanges(1, 1e6))
  expect_error(obj$fetch(wide), "max_fetch_bp")
})

# ---- open_bam ----

test_that("open_bam() errors on missing file", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  expect_error(open_bam("/no/such.bam"), "not found")
})

test_that("open_bam() errors when index is missing", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  tf <- tempfile(fileext = ".bam"); file.create(tf)
  on.exit(unlink(tf))
  expect_error(open_bam(tf), "index")
})

test_that("open_bam() returns SeqBam object", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  bam <- Sys.getenv("SEQPLOTR_TEST_BAM", unset = "")
  skip_if(nchar(bam) == 0, "No test BAM configured")
  obj <- open_bam(bam)
  expect_s3_class(obj, "SeqBam")
  expect_true(is_seq_file_conn(obj))
})

test_that("open_bam() $fetch() errors when region exceeds max_fetch_bp", {
  skip_if_not(requireNamespace("Rsamtools", quietly = TRUE))
  bam <- Sys.getenv("SEQPLOTR_TEST_BAM", unset = "")
  skip_if(nchar(bam) == 0, "No test BAM configured")
  obj  <- open_bam(bam, max_fetch_bp = 1000L)
  wide <- GRanges("chr1", IRanges(1, 1e6))
  expect_error(obj$fetch(wide), "max_fetch_bp")
})

# ---- open_hic ----

test_that("open_hic() errors on missing file", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  expect_error(open_hic("/no/such.hic"), "not found")
})

test_that("open_hic() returns SeqHic object", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  hic <- Sys.getenv("SEQPLOTR_TEST_HIC", unset = "")
  skip_if(nchar(hic) == 0, "No test .hic configured")
  obj <- open_hic(hic)
  expect_s3_class(obj, "SeqHic")
  expect_true(is_seq_file_conn(obj))
})

test_that("open_hic() $fetch() errors without resolution", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  hic <- Sys.getenv("SEQPLOTR_TEST_HIC", unset = "")
  skip_if(nchar(hic) == 0, "No test .hic configured")
  obj <- open_hic(hic)
  win <- GRanges("chr1", IRanges(1, 1e6))
  expect_error(obj$fetch(win), "resolution")
})

test_that("open_hic() $fetch() errors when region exceeds max_fetch_bp", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  hic <- Sys.getenv("SEQPLOTR_TEST_HIC", unset = "")
  skip_if(nchar(hic) == 0, "No test .hic configured")
  obj  <- open_hic(hic, max_fetch_bp = 1000L)
  wide <- GRanges("chr1", IRanges(1, 1e6))
  expect_error(obj$fetch(wide, resolution = 5000L), "max_fetch_bp")
})

# ---- is_seq_file_conn ----

test_that("is_seq_file_conn returns FALSE for non-connection objects", {
  expect_false(is_seq_file_conn(list()))
  expect_false(is_seq_file_conn("path/to/file.bw"))
})
