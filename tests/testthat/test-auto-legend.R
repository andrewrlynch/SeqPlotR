library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(50, 300, 700, 900), width = 50))
mcols(gr)$type  <- c("A", "B", "A", "B")
mcols(gr)$score <- c(0.1, 0.5, 0.8, 0.2)

open_null_dev <- function() grDevices::pdf(file = NULL)

# ── .looks_like_color ──────────────────────────────────────────────────────────

test_that(".looks_like_color identifies hex colors", {
  expect_true(.looks_like_color(c("#FF0000", "#00FF00", "#0000FF")))
})

test_that(".looks_like_color identifies R color names", {
  expect_true(.looks_like_color(c("red", "blue", "steelblue")))
})

test_that(".looks_like_color rejects non-color character values", {
  expect_false(.looks_like_color(c("typeA", "typeB")))
})

test_that(".looks_like_color handles NA (treated as valid for pass-through)", {
  expect_true(.looks_like_color(NA_character_))
})

test_that(".looks_like_color returns FALSE for numeric", {
  expect_false(.looks_like_color(c(0.1, 0.5, 0.9)))
})

# ── .auto_scale_colors ─────────────────────────────────────────────────────────

test_that(".auto_scale_colors discrete returns one color per observation", {
  res <- .auto_scale_colors(c("A", "B", "A"), col_name = "type")
  expect_length(res$colors, 3L)
})

test_that(".auto_scale_colors discrete creates SeqLegendSpec", {
  res <- .auto_scale_colors(c("A", "B"), col_name = "type")
  expect_s3_class(res$legend, "SeqLegendSpec")
  expect_equal(res$legend$title, "type")
  expect_length(res$legend$keys, 2L)
})

test_that(".auto_scale_colors discrete: same level gets same color", {
  res <- .auto_scale_colors(c("A", "B", "A"), col_name = "type")
  expect_equal(res$colors[1], res$colors[3])
})

test_that(".auto_scale_colors continuous returns GradientLegendSpec", {
  res <- .auto_scale_colors(c(0.1, 0.5, 0.9), col_name = "score")
  expect_s3_class(res$legend, "GradientLegendSpec")
  expect_equal(res$legend$title, "score")
})

test_that(".auto_scale_colors continuous colors are hex", {
  res <- .auto_scale_colors(c(0, 0.5, 1), col_name = "x")
  expect_true(.looks_like_color(res$colors))
})

test_that(".auto_scale_colors already-hex returns NULL legend", {
  res <- .auto_scale_colors(c("#FF0000", "#00FF00"), col_name = "color")
  expect_null(res$legend)
  expect_equal(res$colors, c("#FF0000", "#00FF00"))
})

# ── .auto_scale_shapes ─────────────────────────────────────────────────────────

test_that(".auto_scale_shapes creates SeqLegendSpec", {
  res <- .auto_scale_shapes(c("A", "B", "C"), col_name = "group")
  expect_s3_class(res$legend, "SeqLegendSpec")
  expect_length(res$legend$keys, 3L)
})

test_that(".auto_scale_shapes assigns distinct shapes per level", {
  res <- .auto_scale_shapes(c("A", "B"), col_name = "group")
  shapes <- vapply(res$legend$keys, `[[`, character(1), "shape")
  expect_false(shapes[1] == shapes[2])
})

test_that(".auto_scale_shapes returns NULL for numeric input", {
  res <- .auto_scale_shapes(c(1, 2, 3), col_name = "x")
  expect_null(res$shapes)
  expect_null(res$legend)
})

# ── seq_gradient_legend ────────────────────────────────────────────────────────

test_that("seq_gradient_legend creates GradientLegendSpec", {
  g <- seq_gradient_legend("viridis", limits = c(0, 100))
  expect_s3_class(g, "GradientLegendSpec")
  expect_equal(g$palette, "viridis")
  expect_equal(g$limits, c(0, 100))
  expect_null(g$breaks)
})

test_that("seq_gradient_legend rejects invalid palette", {
  expect_error(seq_gradient_legend("magenta"))
})

test_that("seq_gradient_legend rejects non-length-2 limits", {
  expect_error(seq_gradient_legend(limits = c(0, 1, 2)))
})

test_that("seq_gradient_legend stores breaks", {
  g <- seq_gradient_legend(breaks = 5)
  expect_equal(g$breaks, 5)
})

test_that("seq_gradient_legend infers orientation from side", {
  g_left <- seq_gradient_legend(side = "left")
  expect_equal(g_left$orientation, "vertical")
  g_top <- seq_gradient_legend(side = "top")
  expect_equal(g_top$orientation, "horizontal")
})

# ── .gradient_spec_to_keys ────────────────────────────────────────────────────

test_that(".gradient_spec_to_keys returns NULL when breaks is NULL", {
  g <- seq_gradient_legend("viridis", limits = c(0, 1))
  expect_null(.gradient_spec_to_keys(g))
})

test_that(".gradient_spec_to_keys returns n keys for integer breaks", {
  g    <- seq_gradient_legend("viridis", limits = c(0, 1), breaks = 4)
  keys <- .gradient_spec_to_keys(g)
  expect_length(keys, 4L)
})

test_that(".gradient_spec_to_keys returns one key per explicit break value", {
  g    <- seq_gradient_legend("viridis", limits = c(0, 10),
                               breaks = c(2, 5, 8))
  keys <- .gradient_spec_to_keys(g)
  expect_length(keys, 3L)
})

test_that(".gradient_spec_to_keys keys are LegendKey objects", {
  g    <- seq_gradient_legend("plasma", limits = c(0, 1), breaks = 3)
  keys <- .gradient_spec_to_keys(g)
  expect_true(all(vapply(keys, inherits, logical(1), "LegendKey")))
})

# ── print methods ─────────────────────────────────────────────────────────────

test_that("print.GradientLegendSpec runs without error", {
  g <- seq_gradient_legend("viridis", limits = c(0, 1), title = "Score")
  expect_output(print(g), "GradientLegendSpec")
})

# ── auto-legend from seq_point ─────────────────────────────────────────────────

test_that("seq_point sets auto_legend for discrete color mapping", {
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "SeqLegendSpec")
  expect_equal(el$auto_legend$title, "type")
  expect_length(el$auto_legend$keys, 2L)
})

test_that("seq_point sets auto_legend as GradientLegendSpec for numeric color", {
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = score))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "GradientLegendSpec")
})

test_that("seq_point does NOT generate auto_legend when legend is explicit", {
  k   <- LegendKey(label = "manual")
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = type),
                   legend = k)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_null(el$auto_legend)
  expect_true(inherits(el$legend, "LegendKey"))
})

test_that("seq_point does NOT generate auto_legend when show_legend = FALSE", {
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = type),
                   show_legend = FALSE)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_null(el$auto_legend)
})

test_that("seq_point does NOT generate auto_legend for concrete hex colors", {
  mcols(gr)$hex_color <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF")
  el  <- seq_point(data = gr, mapping = map(x = start, y = score,
                                             color = hex_color))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_null(el$auto_legend)
})

test_that("seq_point auto-scales shape mapping", {
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, shape = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "SeqLegendSpec")
})

# ── auto-legend from seq_tile ──────────────────────────────────────────────────

test_that("seq_tile sets auto_legend as GradientLegendSpec for numeric fill", {
  el  <- seq_tile(data = gr, mapping = map(x = start, fill = score))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "GradientLegendSpec")
})

test_that("seq_tile sets auto_legend as SeqLegendSpec for discrete fill", {
  el  <- seq_tile(data = gr, mapping = map(x = start, fill = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "SeqLegendSpec")
})

# ── auto-legend from seq_segment ──────────────────────────────────────────────

test_that("seq_segment sets auto_legend for discrete color mapping", {
  el  <- seq_segment(data = gr, mapping = map(x = start, x_end = end,
                                               y = score, y_end = score,
                                               color = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_s3_class(el$auto_legend, "SeqLegendSpec")
})

# ── drawLegends() dispatches auto_legend ──────────────────────────────────────

test_that("drawLegends() renders discrete auto-legend without error", {
  skip_if_not(capabilities("png"))
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends() renders gradient auto-legend without error", {
  skip_if_not(capabilities("png"))
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = score))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends() renders GradientLegendSpec with breaks without error", {
  skip_if_not(capabilities("png"))
  g   <- seq_gradient_legend("viridis", limits = c(0, 1), breaks = 3,
                              x = 0.02, y = 0.95)
  el  <- seq_tile(data = gr, mapping = map(x = start, fill = score),
                  legend = g)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

test_that("drawLegends() renders GradientLegendSpec as color bar without error", {
  skip_if_not(capabilities("png"))
  g   <- seq_gradient_legend("plasma", limits = c(0, 1), title = "Intensity",
                              x = 0.02, y = 0.95)
  el  <- seq_tile(data = gr, mapping = map(x = start, fill = score),
                  legend = g)
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  sp$drawElements()
  expect_invisible(sp$drawLegends())
})

# ── End-to-end plot() call ────────────────────────────────────────────────────

test_that("full plot() pipeline with auto-legend completes without error", {
  skip_if_not(capabilities("png"))
  el <- seq_point(data = gr, mapping = map(x = start, y = score, color = type))
  sp <- seq_plot() %+% seq_track(windows = win, elements = list(el))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

test_that("full plot() pipeline with gradient legend (color bar) completes", {
  skip_if_not(capabilities("png"))
  el <- seq_tile(data = gr, mapping = map(x = start, fill = score))
  sp <- seq_plot() %+% seq_track(windows = win, elements = list(el))
  open_null_dev(); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

# ── Global suppress still works with auto_legend ──────────────────────────────

test_that("seq_plot(legend=FALSE) suppresses auto-legend rendering", {
  el  <- seq_point(data = gr, mapping = map(x = start, y = score, color = type))
  trk <- seq_track(windows = win, elements = list(el))
  sp  <- seq_plot(legend = FALSE) %+% trk
  sp$layoutGrid()
  open_null_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid(); sp$drawElements()
  expect_invisible(sp$drawLegends())  # renders nothing (show_legend=FALSE)
  expect_false(sp$show_legend)
})
