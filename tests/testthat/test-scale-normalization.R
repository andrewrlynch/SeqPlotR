library(testthat)
library(GenomicRanges)

# ---- .infer_track_scale_factor ----

test_that(".infer_track_scale_factor returns 1e-6 for wide windows (>= 100 kb)", {
  expect_equal(SeqPlotR:::.infer_track_scale_factor(c(1e6, 5e6)), 1e-6)
})

test_that(".infer_track_scale_factor returns 1e-3 for medium windows (>= 100 bp, < 100 kb)", {
  expect_equal(SeqPlotR:::.infer_track_scale_factor(c(5e3, 2e4)), 1e-3)
})

test_that(".infer_track_scale_factor returns 1 for narrow windows (< 100 bp)", {
  expect_equal(SeqPlotR:::.infer_track_scale_factor(c(50, 80)), 1)
})

test_that(".infer_track_scale_factor uses narrowest window", {
  # Mix of Mb-range and kb-range - narrowest drives kb result
  expect_equal(SeqPlotR:::.infer_track_scale_factor(c(5e6, 5e3)), 1e-3)
})

# ---- .window_relative_widths: inferred scale ----

make_win <- function(starts, ends) {
  GRanges("chr1", IRanges(start = starts, end = ends))
}

test_that(".window_relative_widths infers Mb for large windows without track mcols", {
  win <- make_win(c(1, 1e7 + 1), c(1e7, 2e7))
  # Provide a minimal mock track with no window_scale
  trk <- seq_track(windows = win)
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(unique(ww$scale), 1e-6)
})

test_that(".window_relative_widths infers kb when narrowest window is in kb range", {
  win <- make_win(c(1, 5e4 + 1), c(5e4, 1e5))
  trk <- seq_track(windows = win)
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(unique(ww$scale), 1e-3)
})

test_that(".window_relative_widths infers bp for narrow windows", {
  win <- make_win(c(1, 51), c(50, 100))
  trk <- seq_track(windows = win)
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(unique(ww$scale), 1)
})

test_that(".window_relative_widths uses mcols$scale when present (priority 1)", {
  win <- make_win(c(1, 1e7 + 1), c(1e7, 2e7))
  S4Vectors::mcols(win)$scale <- c(1e-3, 1e-3)  # manual kb override
  trk <- seq_track(windows = win)
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(ww$scale, c(1e-3, 1e-3))
})

# ---- window_scale on seq_track ----

test_that("seq_track stores window_scale", {
  win <- make_win(c(1, 5e6 + 1), c(5e6, 1e7))
  trk <- seq_track(windows = win, window_scale = 1e-3)
  expect_equal(trk$window_scale, 1e-3)
})

test_that("window_scale length 1 applied to all windows", {
  win <- make_win(c(1, 5e6 + 1, 1e7 + 1), c(5e6, 1e7, 1.5e7))
  trk <- seq_track(windows = win, window_scale = 1e-3)
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(ww$scale, rep(1e-3, 3))
})

test_that("window_scale positional vector applied correctly", {
  win <- make_win(c(1, 5e6 + 1), c(5e6, 1e7))
  trk <- seq_track(windows = win, window_scale = c(1e-6, 1e-3))
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(ww$scale, c(1e-6, 1e-3))
})

test_that("window_scale wrong length warns and recycles", {
  win <- make_win(c(1, 5e6 + 1, 1e7 + 1), c(5e6, 1e7, 1.5e7))
  trk <- seq_track(windows = win, window_scale = c(1e-6, 1e-3))
  expect_warning(
    SeqPlotR:::.window_relative_widths(win, track = trk),
    "Recycling"
  )
})

test_that("window_scale takes priority over inferred scale", {
  win <- make_win(c(1, 5e6 + 1), c(5e6, 1e7))   # would infer Mb
  trk <- seq_track(windows = win, window_scale = 1)  # force bp
  ww  <- SeqPlotR:::.window_relative_widths(win, track = trk)
  expect_equal(unique(ww$scale), 1)
})

# ---- Layout integration: xScaleFactor stored on panels ----

test_that("layoutGrid stores correct xScaleFactor on panels for large windows", {
  win <- GRanges("chr1", IRanges(c(1, 5e7 + 1), width = c(5e7, 5e7)))
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  panels <- sp$layout$panelBounds[[1]]
  # Both panels should have xScaleFactor = 1e-6 (inferred from 50 Mb windows)
  sf1 <- panels[[1]]$xScaleFactor
  sf2 <- panels[[2]]$xScaleFactor
  expect_equal(sf1, 1e-6)
  expect_equal(sf2, 1e-6)
})

test_that("layoutGrid: mixed-size windows share track-level scale from narrowest", {
  # One 50 Mb window, one 500 kb window - narrowest drives kb (1e-3)
  win <- GRanges("chr1", IRanges(c(1, 5e7 + 1), width = c(5e7, 5e5)))
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  panels <- sp$layout$panelBounds[[1]]
  expect_equal(panels[[1]]$xScaleFactor, 1e-3)
  expect_equal(panels[[2]]$xScaleFactor, 1e-3)
})

test_that("layoutGrid: window_scale overrides inferred scale in panels", {
  win <- GRanges("chr1", IRanges(c(1, 5e7 + 1), width = c(5e7, 5e7)))
  trk <- seq_track(windows = win, window_scale = 1e-3)  # force kb
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  panels <- sp$layout$panelBounds[[1]]
  expect_equal(panels[[1]]$xScaleFactor, 1e-3)
  expect_equal(panels[[2]]$xScaleFactor, 1e-3)
})

test_that("layoutGrid: positional window_scale gives different factors per panel", {
  win <- GRanges("chr1", IRanges(c(1, 5e7 + 1), width = c(5e7, 5e7)))
  trk <- seq_track(windows = win, window_scale = c(1e-6, 1e-3))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  panels <- sp$layout$panelBounds[[1]]
  expect_equal(panels[[1]]$xScaleFactor, 1e-6)
  expect_equal(panels[[2]]$xScaleFactor, 1e-3)
})

# ---- Legacy: no track argument preserves 1e-6 default ----

test_that(".window_relative_widths without track defaults to 1e-6 (legacy)", {
  win <- make_win(c(1, 5e6), c(5e6, 1e7))
  ww  <- SeqPlotR:::.window_relative_widths(win)
  expect_equal(unique(ww$scale), 1e-6)
})
