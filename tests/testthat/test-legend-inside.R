library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(1, 200, 500), width = 50))

# Minimal panel_meta stub matching the shape SeqPlot produces
make_panel <- function(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.8) {
  list(full  = list(x0 = x0, x1 = x1, y0 = y0, y1 = y1),
       inner = list(x0 = x0, x1 = x1, y0 = y0, y1 = y1),
       xscale = c(1, 1000), yscale = c(0, 1))
}

# ── .legend_layout_cells ──────────────────────────────────────────────────────

test_that(".legend_layout_cells returns NULL for zero keys", {
  spec       <- seq_legend(LegendKey(label = "x"))
  spec$keys  <- list()
  expect_null(SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  ))
})

test_that(".legend_layout_cells returns correct number of cells", {
  k1   <- LegendKey(label = "A")
  k2   <- LegendKey(label = "B")
  spec <- seq_legend(list(k1, k2))
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9)
  )
  expect_length(lay$cells, 2)
})

test_that(".legend_layout_cells horizontal default: 1 row, 2 cols for 2 keys", {
  k1   <- LegendKey(label = "A"); k2 <- LegendKey(label = "B")
  spec <- seq_legend(list(k1, k2), orientation = "horizontal")
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  )
  # Second cell must be to the right of first
  expect_gt(lay$cells[[2]]$key_x0, lay$cells[[1]]$key_x0)
})

test_that(".legend_layout_cells vertical default: 2 rows, 1 col for 2 keys", {
  k1   <- LegendKey(label = "A"); k2 <- LegendKey(label = "B")
  spec <- seq_legend(list(k1, k2), orientation = "vertical")
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  )
  # Second cell must be below first (lower y)
  expect_lt(lay$cells[[2]]$y, lay$cells[[1]]$y)
})

test_that(".legend_layout_cells title coords are non-NULL when title set", {
  spec <- seq_legend(LegendKey(label = "A"), title = "Group")
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  )
  expect_false(is.null(lay$title_x))
  expect_false(is.null(lay$title_y))
})

test_that(".legend_layout_cells title coords are NULL when no title", {
  spec <- seq_legend(LegendKey(label = "A"))
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  )
  expect_null(lay$title_x)
})

test_that(".legend_layout_cells nrow/ncol override respected", {
  keys <- list(LegendKey("A"), LegendKey("B"), LegendKey("C"), LegendKey("D"))
  spec <- seq_legend(keys, nrow = 2, ncol = 2)
  lay  <- SeqPlotR:::.legend_layout_cells(
    spec, list(x0 = 0, x1 = 1, y0 = 0, y1 = 1)
  )
  # Item 3 is in row 2 — lower y than item 1
  expect_lt(lay$cells[[3]]$y, lay$cells[[1]]$y)
  # Items 1 and 2 are in the same row — same y
  expect_equal(lay$cells[[2]]$y, lay$cells[[1]]$y)
})

# ── drawLegends() dispatch ────────────────────────────────────────────────────

test_that("SeqPlot$show_legend defaults to TRUE", {
  sp <- SeqPlotR:::SeqPlotR6$new(tracks = list(seq_track(windows = win)))
  expect_true(sp$show_legend)
})

test_that("SeqPlot$show_legend = FALSE suppresses drawLegends", {
  k   <- LegendKey(label = "x")
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = k)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk), show_legend = FALSE)
  sp$layoutGrid()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends runs without error for 'inside' spec", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "Signal", color = "firebrick")
  spec <- seq_legend(k, position = "inside", x = 0.05, y = 0.9)
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends wraps bare LegendKey in default inside spec", {
  k   <- LegendKey(label = "Bare")
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = k)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})


