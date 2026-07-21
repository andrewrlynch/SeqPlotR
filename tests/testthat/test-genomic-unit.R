library(testthat)
library(GenomicRanges)

# ── seq_scale_genomic(unit=, digits=) ─────────────────────────────────────────

test_that("unit forces the scale factor (divisor)", {
  win <- GRanges("chr1", IRanges(1e6, 1.3e6))   # 300 kb window -> auto = kb
  s_mb <- seq_scale_genomic(windows = win, unit = "Mb")
  s_bp <- seq_scale_genomic(windows = win, unit = "bp")
  expect_equal(s_mb$scale_factor[1], 1e-6)
  expect_equal(s_bp$scale_factor[1], 1)
  expect_identical(s_mb$unit, "Mb")
})

test_that("unit is case-insensitive", {
  win <- GRanges("chr1", IRanges(1, 1e6))
  expect_identical(seq_scale_genomic(windows = win, unit = "mb")$unit, "Mb")
  expect_identical(seq_scale_genomic(windows = win, unit = "KB")$unit, "Kb")
  expect_identical(seq_scale_genomic(windows = win, unit = "BP")$unit, "bp")
  expect_error(seq_scale_genomic(windows = win, unit = "furlong"), "Mb")
})

test_that("digits is stored and validated", {
  win <- GRanges("chr1", IRanges(1, 1e6))
  expect_equal(seq_scale_genomic(windows = win, digits = 2)$digits, 2L)
  expect_null(seq_scale_genomic(windows = win)$digits)
  expect_error(seq_scale_genomic(windows = win, digits = -1), "digits")
})

test_that(".axis_fmt honours digits", {
  expect_equal(SeqPlotR:::.axis_fmt(0.3, 2), "0.30")
  expect_equal(SeqPlotR:::.axis_fmt(0.3, NULL), "0.3")
  expect_equal(SeqPlotR:::.axis_fmt(1234.5, 0), "1,234")
})

test_that("auto-detect is still the default when unit is unset", {
  # A 300 kb window infers kb (1e-3); a >=1 Mb window infers Mb (1e-6).
  narrow <- SeqPlotR:::.infer_track_scale_factor(3e5)
  wide   <- SeqPlotR:::.infer_track_scale_factor(2e6)
  expect_equal(narrow, 1e-3)
  expect_equal(wide,   1e-6)
})

test_that("unit + digits flow through the plot pipeline", {
  win <- GRanges("chr1", IRanges(1e6, 1.3e6))
  d   <- GRanges("chr1", IRanges(seq(1e6, 1.3e6, 5e4), width = 1),
                 score = rnorm(7))
  p <- seq_plot() %+%
    seq_track(data = d, windows = win,
              scale_x = seq_scale_genomic(windows = win, unit = "Mb",
                                          digits = 2)) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$scale_x$unit, "Mb")
  expect_equal(trk$scale_x$digits, 2L)
  expect_no_error(p$plot())
})

test_that("manual unit is honoured in combine_windows mode", {
  win <- GRanges("chr1", IRanges(c(1e6, 3e6), width = 3e5))
  d   <- GRanges("chr1", IRanges(c(1.1e6, 3.1e6), width = 1), score = c(1, 2))
  p <- seq_plot() %+%
    seq_track(data = d, windows = win, combine_windows = TRUE,
              scale_x = seq_scale_genomic(windows = win, unit = "bp")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  panels <- p$layout$panelBounds[[1]]
  expanded <- SeqPlotR:::.expand_panels_for_combined_axis(panels)
  # unit = "bp" -> divisor 1 in the expanded per-window panels (previously the
  # helper hardcoded 1e-6 and disagreed with the "bp" suffix).
  expect_true(all(vapply(expanded, function(w) w$xScaleFactor, numeric(1)) == 1))
})
