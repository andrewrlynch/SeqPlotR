library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(1, 200, 500), width = 50))

# ---- trackMarginBounds in layout ----

test_that("layoutGrid produces trackMarginBounds", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  expect_true(!is.null(sp$layout$trackMarginBounds))
  expect_length(sp$layout$trackMarginBounds, 1)
})

test_that("trackMarginBounds has four sides per track", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  sides <- names(sp$layout$trackMarginBounds[[1]])
  expect_setequal(sides, c("top", "bottom", "left", "right"))
})

test_that("trackMarginBounds top: y0 equals trackBounds y1", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(
    tracks = list(trk),
    aesthetics = aes(margins = list(top = 0.05, bottom = 0.05,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  tb <- sp$layout$trackBounds[[1]]
  mb <- sp$layout$trackMarginBounds[[1]]
  expect_equal(mb$top$y0, tb$y1)
})

test_that("trackMarginBounds bottom: y1 equals trackBounds y0", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  tb <- sp$layout$trackBounds[[1]]
  mb <- sp$layout$trackMarginBounds[[1]]
  expect_equal(mb$bottom$y1, tb$y0)
})

test_that("trackMarginBounds left: x1 equals trackBounds x0", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  tb <- sp$layout$trackBounds[[1]]
  mb <- sp$layout$trackMarginBounds[[1]]
  expect_equal(mb$left$x1, tb$x0)
})

test_that("trackMarginBounds right: x0 equals trackBounds x1", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  tb <- sp$layout$trackBounds[[1]]
  mb <- sp$layout$trackMarginBounds[[1]]
  expect_equal(mb$right$x0, tb$x1)
})

test_that("trackMarginBounds is consistent for multi-track plots", {
  trk1 <- seq_track(windows = win)
  trk2 <- seq_track(windows = win)
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk1, trk2))
  sp$layoutGrid()
  expect_length(sp$layout$trackMarginBounds, 2)
})

# ---- .draw_legend_track_margin geometry ----

make_margin_rect <- function() list(x0 = 0.1, x1 = 0.9, y0 = 0.82, y1 = 0.95)

test_that(".draw_legend_track_margin returns invisible for zero-size rect", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "track_margin", side = "top")
  r    <- list(x0 = 0.5, x1 = 0.5, y0 = 0.8, y1 = 0.9)  # zero width
  expect_invisible(SeqPlotR:::.draw_legend_track_margin(spec, r))
})

test_that(".draw_legend_track_margin returns invisible for wrong position", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "inside")
  r    <- make_margin_rect()
  expect_invisible(SeqPlotR:::.draw_legend_track_margin(spec, r))
})

# ---- drawLegends dispatch for track_margin ----

test_that("drawLegends dispatches track_margin without error", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "Signal")
  spec <- seq_legend(k, position = "track_margin", side = "top")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends track_margin no longer emits batch 7E warning", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "track_margin", side = "bottom")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_no_warning(sp$drawLegends())
})

test_that("drawLegends left-margin spec uses vertical orientation by default", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "track_margin", side = "left")
  expect_equal(spec$orientation, "vertical")
})
