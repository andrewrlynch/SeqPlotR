library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 2000))

gr_disc <- GRanges("chr1",
  IRanges(c(100, 400, 800, 1200), width = 80),
  group = c("A", "B", "A", "C"),
  score = c(1.2, 3.4, 2.1, 4.5)
)

gr_cont <- GRanges("chr1",
  IRanges(c(100, 400, 800, 1200), width = 80),
  score = c(0.1, 0.5, 0.8, 0.3)
)

open_null_dev <- function() grDevices::pdf(file = NULL)

# Helper: build a plot pipeline and call drawElements() so prep() runs
run_pipeline <- function(el) {
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev()
  on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  invisible(sp)
}

# ---- .element_shape_for ----

test_that(".element_shape_for returns 'point' for SeqPoint", {
  el <- seq_point(data = gr_disc, mapping = map(x = start, y = score))
  expect_equal(.element_shape_for(el), "point")
})

test_that(".element_shape_for returns 'rect' for SeqBar", {
  el <- seq_bar(data = gr_disc, mapping = map(x = start, y = score))
  expect_equal(.element_shape_for(el), "rect")
})

test_that(".element_shape_for returns 'line' for SeqLine", {
  el <- seq_line(data = gr_disc, mapping = map(x = start, y = score))
  expect_equal(.element_shape_for(el), "line")
})

test_that(".element_shape_for returns 'rect' for unknown class", {
  fake <- structure(list(), class = "UnknownElement")
  expect_equal(.element_shape_for(fake), "rect")
})

# ---- Auto-legend: SeqBar ----

test_that("SeqBar discrete fill map produces auto_legend after prep()", {
  el <- seq_bar(data = gr_disc,
                mapping = map(x = start, y = score, fill = group))
  run_pipeline(el)
  expect_false(is.null(el$auto_legend))
})

test_that("SeqBar auto_legend is a SeqLegendSpec for discrete fill", {
  el <- seq_bar(data = gr_disc,
                mapping = map(x = start, y = score, fill = group))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "SeqLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "SeqLegendSpec")
  expect_equal(length(spec$keys), 3L)  # A, B, C
})

test_that("SeqBar continuous fill map produces GradientLegendSpec", {
  el <- seq_bar(data = gr_cont,
                mapping = map(x = start, y = score, fill = score))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "GradientLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "GradientLegendSpec")
})

test_that("SeqBar fill auto-scale writes hex colors to coords", {
  el <- seq_bar(data = gr_disc,
                mapping = map(x = start, y = score, fill = group))
  run_pipeline(el)
  # All non-NULL coordCanvas fill vectors should be valid hex colors
  for (coords in el$coordCanvas) {
    if (!is.null(coords) && nrow(coords) > 0L)
      expect_true(.looks_like_color(coords$fill))
  }
})

# ---- Auto-legend: SeqLine ----

test_that("SeqLine discrete color map produces auto_legend after prep()", {
  el <- seq_line(data = gr_disc,
                 mapping = map(x = start, y = score, color = group))
  run_pipeline(el)
  expect_false(is.null(el$auto_legend))
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "SeqLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "SeqLegendSpec")
})

test_that("SeqLine continuous color map produces GradientLegendSpec", {
  el <- seq_line(data = gr_cont,
                 mapping = map(x = start, y = score, color = score))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "GradientLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "GradientLegendSpec")
})

# ---- Auto-legend: SeqArea ----

test_that("SeqArea continuous fill map produces GradientLegendSpec", {
  el <- seq_area(data = gr_cont,
                 mapping = map(x = start, y = score, fill = score))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "GradientLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "GradientLegendSpec")
})

test_that("SeqArea discrete color map produces SeqLegendSpec", {
  el <- seq_area(data = gr_disc,
                 mapping = map(x = start, y = score, color = group))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "SeqLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_s3_class(spec, "SeqLegendSpec")
})

# ---- No auto-legend without mapped aesthetic ----

test_that("SeqBar without fill/color mapping produces no auto_legend", {
  el <- seq_bar(data = gr_disc, mapping = map(x = start, y = score))
  run_pipeline(el)
  expect_null(el$auto_legend)
})

test_that("SeqLine without color mapping produces no auto_legend", {
  el <- seq_line(data = gr_disc, mapping = map(x = start, y = score))
  run_pipeline(el)
  expect_null(el$auto_legend)
})

# ---- Manual legend overrides auto-legend ----

test_that("manual legend field suppresses auto_legend for SeqBar", {
  k  <- LegendKey(label = "Manual")
  el <- seq_bar(data = gr_disc,
                mapping = map(x = start, y = score, fill = group),
                legend  = k)
  run_pipeline(el)
  expect_null(el$auto_legend)
  expect_true(inherits(el$legend, "LegendKey"))
})

# ---- show_legend = FALSE suppresses auto-legend ----

test_that("show_legend = FALSE prevents auto_legend for SeqBar", {
  el <- seq_bar(data = gr_disc,
                mapping     = map(x = start, y = score, fill = group),
                show_legend = FALSE)
  run_pipeline(el)
  expect_null(el$auto_legend)
})

test_that("show_legend = FALSE prevents auto_legend for SeqLine", {
  el <- seq_line(data = gr_disc,
                 mapping     = map(x = start, y = score, color = group),
                 show_legend = FALSE)
  run_pipeline(el)
  expect_null(el$auto_legend)
})

# ---- Auto-legend position default is track_margin right ----

test_that("auto-legend SeqLegendSpec has position track_margin by default", {
  el <- seq_bar(data = gr_disc,
                mapping = map(x = start, y = score, fill = group))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "SeqLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_equal(spec$position, "track_margin")
  expect_equal(spec$side, "right")
})

test_that("auto-legend GradientLegendSpec has position track_margin by default", {
  el <- seq_bar(data = gr_cont,
                mapping = map(x = start, y = score, fill = score))
  run_pipeline(el)
  spec <- if (is.list(el$auto_legend) &&
              !inherits(el$auto_legend, "GradientLegendSpec"))
    el$auto_legend[[1L]] else el$auto_legend
  expect_equal(spec$position, "track_margin")
  expect_equal(spec$side, "right")
})

# ---- Colorbar axis ticks ----

test_that("GradientLegendSpec with breaks = NULL uses pretty ticks (no error)", {
  spec <- seq_gradient_legend(palette = "viridis", limits = c(0, 100),
                               title = "Score", breaks = NULL)
  expect_null(spec$breaks)
  skip_if_not(capabilities("png"))
  panel <- list(full  = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9),
                inner = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.draw_gradient_legend(spec, panel))
})

test_that("GradientLegendSpec with breaks = 3L renders tick axis without error", {
  spec <- seq_gradient_legend(palette = "viridis", limits = c(0, 1),
                               title = "Score", breaks = 3L)
  expect_equal(spec$breaks, 3L)
  skip_if_not(capabilities("png"))
  panel <- list(full  = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9),
                inner = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.draw_gradient_legend(spec, panel))
})

test_that("GradientLegendSpec with explicit breaks vector renders without error", {
  spec <- seq_gradient_legend(palette = "plasma", limits = c(0, 10),
                               title = "Depth",
                               breaks = c(0, 2.5, 5, 7.5, 10))
  skip_if_not(capabilities("png"))
  panel <- list(full  = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9),
                inner = list(x0 = 0.05, x1 = 0.45, y0 = 0.1, y1 = 0.9))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.draw_gradient_legend(spec, panel))
})

test_that("colorbar horizontal orientation renders without error", {
  spec <- seq_gradient_legend(palette = "blues", limits = c(0, 1),
                               title = "Score", orientation = "horizontal",
                               x = 0.5, y = 0.9)
  skip_if_not(capabilities("png"))
  panel <- list(full  = list(x0 = 0.05, x1 = 0.95, y0 = 0.1, y1 = 0.9),
                inner = list(x0 = 0.05, x1 = 0.95, y0 = 0.1, y1 = 0.9))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.draw_gradient_legend(spec, panel))
})

# ---- drawLegends() auto_legend dispatch ----

test_that("drawLegends() renders SeqBar fill auto_legend without error", {
  skip_if_not(capabilities("png"))
  el  <- seq_bar(data = gr_disc,
                 mapping = map(x = start, y = score, fill = group))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends() renders SeqLine GradientLegendSpec without error", {
  skip_if_not(capabilities("png"))
  el  <- seq_line(data = gr_cont,
                  mapping = map(x = start, y = score, color = score))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

test_that("seq_plot(legend = FALSE) suppresses SeqBar auto_legend rendering", {
  el  <- seq_bar(data = gr_disc,
                 mapping = map(x = start, y = score, fill = group))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot(legend = FALSE) %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
  expect_false(sp$show_legend)
})

# ---- .SeqLegend_drawKey rect shape ----

test_that(".SeqLegend_drawKey renders rect shape without error", {
  skip_if_not(capabilities("png"))
  k <- LegendKey(label = "Group A", color = "#3D3D3A", fill = "#4385BE",
                 shape = "rect")
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.SeqLegend_drawKey(k, x0 = 0.1, x1 = 0.15, y = 0.5,
                                      height = 0.05))
})

test_that(".SeqLegend_drawKey renders point shape without error", {
  skip_if_not(capabilities("png"))
  k <- LegendKey(label = "Obs", color = "#AF3029", fill = "#AF3029",
                 shape = "point")
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(.SeqLegend_drawKey(k, x0 = 0.1, x1 = 0.15, y = 0.5,
                                      height = 0.05))
})
