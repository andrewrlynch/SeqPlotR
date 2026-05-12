library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(1, 200, 500), width = 50))

# ---- canvasMarginBounds in layout ----

test_that("layoutGrid produces canvasMarginBounds", {
  trk <- seq_track(windows = win)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  expect_false(is.null(sp$layout$canvasMarginBounds))
  expect_setequal(names(sp$layout$canvasMarginBounds),
                  c("top", "bottom", "left", "right"))
})

test_that("canvasMarginBounds top: y0 equals 1 - margins$top", {
  sp <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(seq_track(windows = win)),
    aesthetics = aes(margins = list(top = 0.08, bottom = 0.05,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  expect_equal(sp$layout$canvasMarginBounds$top$y0, 1 - 0.08)
  expect_equal(sp$layout$canvasMarginBounds$top$y1, 1)
})

test_that("canvasMarginBounds bottom: y1 equals margins$bottom", {
  sp <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(seq_track(windows = win)),
    aesthetics = aes(margins = list(top = 0.05, bottom = 0.06,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  expect_equal(sp$layout$canvasMarginBounds$bottom$y0, 0)
  expect_equal(sp$layout$canvasMarginBounds$bottom$y1, 0.06)
})

# ---- .collect_canvas_legend_specs ----

test_that(".collect_canvas_legend_specs returns empty list when no canvas specs", {
  k   <- LegendKey(label = "x")
  s   <- seq_legend(k, position = "inside")
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = s)
  trk <- seq_track(windows = win, elements = list(el))
  expect_length(SeqPlotR:::.collect_canvas_legend_specs(list(trk)), 0)
})

test_that(".collect_canvas_legend_specs collects canvas_margin specs", {
  k   <- LegendKey(label = "x")
  s   <- seq_legend(k, position = "canvas_margin", side = "top")
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = s)
  trk <- seq_track(windows = win, elements = list(el))
  res <- SeqPlotR:::.collect_canvas_legend_specs(list(trk))
  expect_length(res, 1)
  expect_equal(res[[1]]$track_idx, 1)
})

test_that(".collect_canvas_legend_specs respects track show_legend = FALSE", {
  k   <- LegendKey(label = "x")
  s   <- seq_legend(k, position = "canvas_margin")
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = s)
  trk <- seq_track(windows = win, elements = list(el), show_legend = FALSE)
  expect_length(SeqPlotR:::.collect_canvas_legend_specs(list(trk)), 0)
})

# ---- .merge_canvas_specs ----

test_that(".merge_canvas_specs returns NULL for empty entries", {
  expect_null(SeqPlotR:::.merge_canvas_specs(list(), "top"))
})

test_that(".merge_canvas_specs merges two same-side specs into one", {
  k1 <- LegendKey(label = "A", title = "Marks")
  k2 <- LegendKey(label = "B", title = "Marks")
  s1 <- seq_legend(k1, position = "canvas_margin", side = "top")
  s2 <- seq_legend(k2, position = "canvas_margin", side = "top")
  entries <- list(list(spec = s1, track_idx = 1L),
                  list(spec = s2, track_idx = 2L))
  merged  <- SeqPlotR:::.merge_canvas_specs(entries, "top")
  expect_s3_class(merged, "SeqLegendSpec")
  expect_gte(length(merged$keys), 2)
})

test_that(".merge_canvas_specs inserts separator for distinct titles", {
  k1 <- LegendKey(label = "A")
  k2 <- LegendKey(label = "B")
  s1 <- seq_legend(k1, title = "Group1", position = "canvas_margin", side = "top")
  s2 <- seq_legend(k2, title = "Group2", position = "canvas_margin", side = "top")
  entries <- list(list(spec = s1, track_idx = 1L),
                  list(spec = s2, track_idx = 2L))
  merged  <- SeqPlotR:::.merge_canvas_specs(entries, "top")
  # 2 separators + 2 real keys = 4
  expect_length(merged$keys, 4)
  shapes <- vapply(merged$keys, function(k) as.character(k$shape), character(1))
  expect_equal(sum(shapes == "none"), 2)
})

# ---- full drawLegends dispatch ----

test_that("drawLegends dispatches canvas_margin without error", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "Signal", color = "firebrick")
  spec <- seq_legend(k, position = "canvas_margin", side = "top",
                     title = "ChIP")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(trk),
    aesthetics = aes(margins = list(top = 0.1, bottom = 0.05,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends canvas_margin no longer emits warning", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "canvas_margin", side = "bottom")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_no_warning(sp$drawLegends())
})

test_that("multi-track canvas aggregation draws without error", {
  skip_if_not(capabilities("png"))
  k1   <- LegendKey(label = "H3K27ac", color = "firebrick")
  k2   <- LegendKey(label = "H3K4me3", color = "steelblue")
  s1   <- seq_legend(k1, position = "canvas_margin", side = "top",
                     title = "Marks")
  s2   <- seq_legend(k2, position = "canvas_margin", side = "top",
                     title = "Marks")
  el1  <- SeqPlotR:::SeqElementR6$new(gr, legend = s1)
  el2  <- SeqPlotR:::SeqElementR6$new(gr, legend = s2)
  trk1 <- seq_track(windows = win, elements = list(el1))
  trk2 <- seq_track(windows = win, elements = list(el2))
  sp   <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(trk1, trk2),
    aesthetics = aes(margins = list(top = 0.1, bottom = 0.05,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})
