library(testthat)
library(GenomicRanges)

# ── .resolve_bin_fun ──────────────────────────────────────────────────────────

test_that(".resolve_bin_fun returns mean function", {
  f <- SeqPlotR:::.resolve_bin_fun("mean")
  expect_equal(f(c(1, 3, 5)), 3)
})

test_that(".resolve_bin_fun returns median function", {
  f <- SeqPlotR:::.resolve_bin_fun("median")
  expect_equal(f(c(1, 3, 100)), 3)
})

test_that(".resolve_bin_fun returns sum function", {
  f <- SeqPlotR:::.resolve_bin_fun("sum")
  expect_equal(f(c(1, 2, 3)), 6)
})

test_that(".resolve_bin_fun accepts a custom function", {
  f <- SeqPlotR:::.resolve_bin_fun(function(x) max(x))
  expect_equal(f(c(1, 5, 3)), 5)
})

test_that(".resolve_bin_fun errors on unknown string", {
  expect_error(SeqPlotR:::.resolve_bin_fun("geometric"), "must be")
})

# ── .bin_signal_gr ────────────────────────────────────────────────────────────

make_signal_gr <- function() {
  # Four ranges covering bp 1-400 with score = 10, 20, 30, 40
  GRanges("chr1",
    IRanges(start = c(1, 101, 201, 301), end = c(100, 200, 300, 400)),
    score = c(10, 20, 30, 40)
  )
}

test_that(".bin_signal_gr returns one row per bin", {
  gr  <- make_signal_gr()
  win <- GRanges("chr1", IRanges(1, 400))
  df  <- SeqPlotR:::.bin_signal_gr(gr, win, bin_size = 100L,
                                   agg_fun = mean)
  expect_equal(nrow(df), 4L)
})

test_that(".bin_signal_gr mean aggregation is correct for non-overlapping ranges", {
  gr  <- make_signal_gr()
  win <- GRanges("chr1", IRanges(1, 400))
  df  <- SeqPlotR:::.bin_signal_gr(gr, win, bin_size = 100L,
                                   agg_fun = SeqPlotR:::.resolve_bin_fun("mean"))
  # Each bin maps exactly to one range — score equals that range's score
  expect_equal(df$score, c(10, 20, 30, 40))
})

test_that(".bin_signal_gr sum aggregation across overlapping bins is correct", {
  # 200 bp bins: bins [1,200] and [201,400] each capture two 100-bp ranges
  gr  <- make_signal_gr()
  win <- GRanges("chr1", IRanges(1, 400))
  df  <- SeqPlotR:::.bin_signal_gr(gr, win, bin_size = 200L,
                                   agg_fun = SeqPlotR:::.resolve_bin_fun("sum"))
  expect_equal(nrow(df), 2L)
  expect_true(all(is.finite(df$score)))
})

test_that(".bin_signal_gr returns NA for empty bins", {
  # Signal sits entirely in bin 4 (bp 301-400). Bins 1-3 should be NA.
  gr  <- GRanges("chr1", IRanges(301, 400), score = 5)
  win <- GRanges("chr1", IRanges(1, 400))
  df  <- SeqPlotR:::.bin_signal_gr(gr, win, bin_size = 100L,
                                   agg_fun = SeqPlotR:::.resolve_bin_fun("mean"))
  expect_equal(nrow(df), 4L)
  expect_true(is.na(df$score[1]))
  expect_true(is.na(df$score[2]))
  expect_true(is.na(df$score[3]))
  expect_false(is.na(df$score[4]))
})

test_that(".bin_signal_gr output columns are correct", {
  gr  <- make_signal_gr()
  win <- GRanges("chr1", IRanges(1, 400))
  df  <- SeqPlotR:::.bin_signal_gr(gr, win, bin_size = 100L,
                                   agg_fun = SeqPlotR:::.resolve_bin_fun("mean"))
  expect_named(df, c("seqnames", "start", "end", "score"))
})

# ── .rebin_contacts ───────────────────────────────────────────────────────────

make_contacts <- function() {
  data.frame(
    seqnames1 = "chr1", start1 = c(0, 5000, 10000, 15000),
    end1      = c(4999, 9999, 14999, 19999),
    seqnames2 = "chr1", start2 = c(0, 5000, 10000, 15000),
    end2      = c(4999, 9999, 14999, 19999),
    score     = c(10, 20, 30, 40),
    stringsAsFactors = FALSE
  )
}

test_that(".rebin_contacts returns empty df for empty input", {
  df  <- make_contacts()[0, ]
  out <- SeqPlotR:::.rebin_contacts(df, 10000L,
                                    SeqPlotR:::.resolve_bin_fun("sum"))
  expect_equal(nrow(out), 0L)
})

test_that(".rebin_contacts aggregates 5 kb bins into 10 kb bins (sum)", {
  df  <- make_contacts()    # four 5 kb bins
  out <- SeqPlotR:::.rebin_contacts(df, 10000L,
                                    SeqPlotR:::.resolve_bin_fun("sum"))
  # (0,0) and (5000,5000) snap to bin (0,0); (10000,10000) and (15000,15000)
  # snap to bin (10000,10000) — two output bins
  expect_equal(nrow(out), 2L)
  expect_true(all(out$score %in% c(10 + 20, 30 + 40)))
})

test_that(".rebin_contacts mean aggregation is correct", {
  df  <- make_contacts()
  out <- SeqPlotR:::.rebin_contacts(df, 10000L,
                                    SeqPlotR:::.resolve_bin_fun("mean"))
  expect_equal(nrow(out), 2L)
  expect_true(all(out$score %in% c((10 + 20) / 2, (30 + 40) / 2)))
})

test_that(".rebin_contacts output has correct columns", {
  df  <- make_contacts()
  out <- SeqPlotR:::.rebin_contacts(df, 10000L,
                                    SeqPlotR:::.resolve_bin_fun("sum"))
  expect_named(out, c("seqnames1", "start1", "end1",
                      "seqnames2", "start2", "end2", "score"))
})

# ── open_bigwig() $fetch_binned() ─────────────────────────────────────────────

bw_path <- Sys.getenv("SEQPLOTR_TEST_BW", unset = "")

test_that("open_bigwig() $fetch_binned() errors on non-positive bin_size", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  skip_if(nchar(bw_path) == 0, "No test bigWig configured")
  bw  <- open_bigwig(bw_path)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  expect_error(bw$fetch_binned(win, bin_size = 0), "positive")
  expect_error(bw$fetch_binned(win, bin_size = -1), "positive")
})

test_that("open_bigwig() $fetch_binned() returns data.frame with correct cols", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  skip_if(nchar(bw_path) == 0, "No test bigWig configured")
  bw  <- open_bigwig(bw_path)
  win <- GRanges("chr1", IRanges(1e6, 1.1e6))
  df  <- bw$fetch_binned(win, bin_size = 1000L, fun = "mean")
  expect_named(df, c("seqnames", "start", "end", "score"))
})

test_that("open_bigwig() $fetch_binned() returns correct number of bins", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  skip_if(nchar(bw_path) == 0, "No test bigWig configured")
  bw  <- open_bigwig(bw_path)
  win <- GRanges("chr1", IRanges(1e6, 1.1e6))   # 100 kb window
  df  <- bw$fetch_binned(win, bin_size = 10000L)
  expect_equal(nrow(df), 10L)
})

test_that("open_bigwig() $fetch_binned() accepts custom function", {
  skip_if_not(requireNamespace("rtracklayer", quietly = TRUE))
  skip_if(nchar(bw_path) == 0, "No test bigWig configured")
  bw  <- open_bigwig(bw_path)
  win <- GRanges("chr1", IRanges(1e6, 1.1e6))
  df  <- bw$fetch_binned(win, bin_size = 10000L, fun = function(x) max(x))
  expect_named(df, c("seqnames", "start", "end", "score"))
})

# ── open_hic() $fetch_binned() ────────────────────────────────────────────────

hic_path <- Sys.getenv("SEQPLOTR_TEST_HIC", unset = "")

test_that("open_hic() $fetch_binned() errors when bin_size <= resolution", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  skip_if(nchar(hic_path) == 0, "No test .hic configured")
  hic <- open_hic(hic_path)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  expect_error(
    hic$fetch_binned(win, resolution = 5000L, bin_size = 5000L),
    "must be larger"
  )
})

test_that("open_hic() $fetch_binned() returns BEDPE data.frame", {
  skip_if_not(requireNamespace("strawr", quietly = TRUE))
  skip_if(nchar(hic_path) == 0, "No test .hic configured")
  hic <- open_hic(hic_path)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  df  <- hic$fetch_binned(win, resolution = 5000L, bin_size = 25000L)
  expect_named(df, c("seqnames1", "start1", "end1",
                     "seqnames2", "start2", "end2", "score"))
})

# ── open_h5() ─────────────────────────────────────────────────────────────────

h5_path <- Sys.getenv("SEQPLOTR_TEST_H5", unset = "")

test_that("open_h5() errors on missing file", {
  expect_error(open_h5("/no/such.cool"), "not found")
})

test_that("open_h5() returns SeqH5 object", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path)
  expect_s3_class(obj, "SeqH5")
  expect_true(is_seq_file_conn(obj))
})

test_that("open_h5() populates chromosome table", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path)
  expect_true(nrow(obj$chromosomes) > 0L)
  expect_named(obj$chromosomes, c("name", "length"))
})

test_that("open_h5() $fetch() errors without resolution for mcool", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path, resolution = NULL)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  expect_error(obj$fetch(win, resolution = NULL), "resolution")
})

test_that("open_h5() $fetch() returns BEDPE data.frame", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  df  <- obj$fetch(win)
  expect_named(df, c("seqnames1", "start1", "end1",
                     "seqnames2", "start2", "end2", "score"))
})

test_that("open_h5() $fetch() respects max_fetch_bp", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj  <- open_h5(h5_path, max_fetch_bp = 1000L)
  wide <- GRanges("chr1", IRanges(1, 1e6))
  expect_error(obj$fetch(wide), "max_fetch_bp")
})

test_that("open_h5() $fetch_binned() errors when bin_size <= resolution", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path)
  win <- GRanges("chr1", IRanges(1e6, 2e6))
  res <- obj$resolution
  expect_error(obj$fetch_binned(win, bin_size = res), "must be larger")
})

test_that("open_h5() $fetch_binned() returns coarser bins", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj     <- open_h5(h5_path)
  win     <- GRanges("chr1", IRanges(1e6, 2e6))
  native  <- obj$fetch(win)
  coarser <- obj$fetch_binned(win, bin_size = obj$resolution * 5L)
  if (nrow(native) > 0 && nrow(coarser) > 0)
    expect_lte(nrow(coarser), nrow(native))
})

test_that("is_seq_file_conn is TRUE for SeqH5", {
  skip_if_not(requireNamespace("rhdf5", quietly = TRUE))
  skip_if(nchar(h5_path) == 0, "No test H5 configured")
  obj <- open_h5(h5_path)
  expect_true(is_seq_file_conn(obj))
})
