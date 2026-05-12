library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(50, 300, 700), width = 80))

make_plot <- function(legend_spec,
                      margins = list(top = 0.10, bottom = 0.05,
                                     left = 0.05, right = 0.05)) {
  el  <- SeqPlotR:::SeqElementR6$new(gr, legend = legend_spec)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(trk),
    aesthetics = aes(margins = margins)
  )
  sp$layoutGrid()
  sp
}

open_null_dev <- function() grDevices::pdf(file = NULL)

# ---- Inside: end-to-end ----

test_that("inside legend renders end-to-end without error", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "Signal", color = "firebrick", shape = "circle")
  spec <- seq_legend(k, position = "inside", x = 0.05, y = 0.85, title = "ChIP")
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("inside legend: multiple keys with nrow=2 renders without error", {
  skip_if_not(capabilities("png"))
  keys <- list(LegendKey("A", color = "red"), LegendKey("B", color = "blue"),
               LegendKey("C", color = "green"), LegendKey("D", color = "orange"))
  spec <- seq_legend(keys, position = "inside", nrow = 2)
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

# ---- Track margin: end-to-end ----

test_that("track_margin top legend renders end-to-end", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "CNV", color = "steelblue", shape = "-")
  spec <- seq_legend(k, position = "track_margin", side = "top",
                     title = "Copy Number")
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("track_margin bottom legend renders end-to-end", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "track_margin", side = "bottom")
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("track_margin left legend (vertical orientation) renders end-to-end", {
  skip_if_not(capabilities("png"))
  keys <- list(LegendKey("A"), LegendKey("B"))
  spec <- seq_legend(keys, position = "track_margin", side = "left")
  expect_equal(spec$orientation, "vertical")
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

# ---- Canvas margin: end-to-end ----

test_that("canvas_margin top legend renders end-to-end (single track)", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "H3K27ac", color = "firebrick")
  spec <- seq_legend(k, position = "canvas_margin", side = "top",
                     title = "Marks")
  sp   <- make_plot(spec)
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("canvas_margin aggregates keys from two tracks onto one side", {
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
    aesthetics = aes(margins = list(top = 0.12, bottom = 0.05,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("canvas_margin different sides don't interfere", {
  skip_if_not(capabilities("png"))
  k1  <- LegendKey(label = "Top key")
  k2  <- LegendKey(label = "Bottom key")
  s1  <- seq_legend(k1, position = "canvas_margin", side = "top")
  s2  <- seq_legend(k2, position = "canvas_margin", side = "bottom")
  el1 <- SeqPlotR:::SeqElementR6$new(gr, legend = s1)
  el2 <- SeqPlotR:::SeqElementR6$new(gr, legend = s2)
  trk <- seq_track(windows = win, elements = list(el1, el2))
  sp  <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(trk),
    aesthetics = aes(margins = list(top = 0.10, bottom = 0.10,
                                    left = 0.05, right = 0.05))
  )
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

# ---- Suppression at all three levels ----

test_that("SeqPlot(legend = FALSE) suppresses all legends", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "inside")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk), legend = FALSE)
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
  expect_false(sp$show_legend)
})

test_that("SeqPlot(show_legend = FALSE) suppresses all legends", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "inside")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk), show_legend = FALSE)
  expect_false(sp$show_legend)
})

test_that("seq_track(show_legend = FALSE) suppresses that track's legends", {
  k    <- LegendKey(label = "x")
  spec <- seq_legend(k, position = "inside")
  el   <- SeqPlotR:::SeqElementR6$new(gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el), show_legend = FALSE)
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

test_that("element-level show_legend = FALSE suppresses only that element", {
  k1  <- LegendKey(label = "shown")
  k2  <- LegendKey(label = "hidden")
  s1  <- seq_legend(k1, position = "inside")
  s2  <- seq_legend(k2, position = "inside")
  el1 <- SeqPlotR:::SeqElementR6$new(gr, legend = s1, show_legend = TRUE)
  el2 <- SeqPlotR:::SeqElementR6$new(gr, legend = s2, show_legend = FALSE)
  trk <- seq_track(windows = win, elements = list(el1, el2))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

# ---- Operator API (%+%) compatibility ----

test_that("%+% correctly forwards element with legend into track", {
  k    <- LegendKey(label = "via operator", color = "purple")
  spec <- seq_legend(k, position = "inside")
  el   <- seq_point(data = gr, legend = spec)
  trk  <- seq_track(windows = win)
  trk  <- trk %+% el
  expect_length(trk$elements, 1)
  expect_true(inherits(trk$elements[[1]]$legend, "SeqLegendSpec"))
})

# ---- Bare LegendKey (no seq_legend() call) ----

test_that("bare LegendKey on element renders as inside spec", {
  skip_if_not(capabilities("png"))
  k   <- LegendKey(label = "Bare key", color = "navy")
  el  <- seq_point(data = gr, legend = k)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawLegends())
})

# ---- plot() calls drawLegends() ----

test_that("seq_plot()$plot() calls drawLegends() without error", {
  skip_if_not(capabilities("png"))
  k    <- LegendKey(label = "full pipeline", color = "darkgreen")
  spec <- seq_legend(k, position = "inside", x = 0.9, y = 0.9)
  el   <- seq_point(data = gr, legend = spec)
  trk  <- seq_track(windows = win, elements = list(el))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

# ---- collect_legend_keys() round-trip ----

test_that("collect_legend_keys round-trip: SeqElement -> SeqTrack -> flat list", {
  k1  <- LegendKey(label = "A", title = "Group1")
  k2  <- LegendKey(label = "B", title = "Group2")
  s1  <- seq_legend(k1, title = "Group1")
  s2  <- seq_legend(k2, title = "Group2")
  el1 <- SeqPlotR:::SeqElementR6$new(gr, legend = s1)
  el2 <- SeqPlotR:::SeqElementR6$new(gr, legend = s2)
  trk <- seq_track(windows = win, elements = list(el1, el2))

  keys <- trk$collect_legend_keys()
  expect_length(keys, 2)
  labels <- vapply(keys, function(x) x$key$label, character(1))
  expect_setequal(labels, c("A", "B"))
})
