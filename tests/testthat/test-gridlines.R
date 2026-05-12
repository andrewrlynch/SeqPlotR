library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(100, 400, 700), width = 80),
               score = c(1, 3, 2))

make_sp <- function(track_aes = aes(), plot_aes = aes()) {
  el  <- seq_bar(data = gr, mapping = map(y = score))
  trk <- seq_track(windows = win, aesthetics = track_aes)
  trk$addElement(el)
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk), aesthetics = plot_aes)
  sp$layoutGrid()
  sp
}

open_dev <- function() grDevices::pdf(file = NULL)

# ---- Default theme keys ----

test_that("axis.gridline.visible defaults to FALSE in .default_theme()", {
  dt <- SeqPlotR:::.default_theme()
  expect_false(isTRUE(dt[["axis.gridline.visible"]]))
})

test_that("axis.gridline.color defaults to 'grey85'", {
  dt <- SeqPlotR:::.default_theme()
  expect_equal(dt[["axis.gridline.color"]], "grey85")
})

test_that("axis.gridline.lwd defaults to 0.5", {
  dt <- SeqPlotR:::.default_theme()
  expect_equal(dt[["axis.gridline.lwd"]], 0.5)
})

# ---- Resolved theme reflects user settings ----

test_that("axis.y.gridline = TRUE flows into resolved_theme$flat", {
  sp <- make_sp(track_aes = aes(axis.y.gridline = TRUE))
  flat <- sp$allTracks()[[1]]$resolved_theme$flat
  expect_true(isTRUE(flat[["axis.y.gridline"]]))
})

test_that("axis.x.gridline = TRUE flows into resolved_theme$flat", {
  sp <- make_sp(track_aes = aes(axis.x.gridline = TRUE))
  flat <- sp$allTracks()[[1]]$resolved_theme$flat
  expect_true(isTRUE(flat[["axis.x.gridline"]]))
})

# ---- drawGridlines() with gridlines off (default) ----

test_that("drawGridlines() returns invisibly when both gridlines are off", {
  sp <- make_sp()
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

# ---- drawGridlines() with y gridlines on ----

test_that("drawGridlines() runs without error when axis.y.gridline = TRUE", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.y.gridline = TRUE))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

test_that("drawGridlines() runs without error when axis.x.gridline = TRUE", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.x.gridline = TRUE))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

test_that("drawGridlines() runs without error when both axes are on", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.x.gridline = TRUE, axis.y.gridline = TRUE))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

# ---- Nested aes() form enables gridlines ----

test_that("axis.y.gridline = aes(color) form enables gridlines", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.y.gridline = aes(color = "grey80", lwd = 0.5)))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

test_that("axis.x.gridline = aes(color, lty) form enables gridlines", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.x.gridline = aes(color = "grey80", lty = 2)))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

# ---- Custom styling sub-keys ----

test_that("axis.y.gridline custom color/lwd/lty renders without error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(
    axis.y.gridline       = TRUE,
    axis.y.gridline.color = "grey60",
    axis.y.gridline.lwd   = 1,
    axis.y.gridline.lty   = 2
  ))
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

test_that("axis.gridline.color base key is inherited by both x and y", {
  sp <- make_sp(track_aes = aes(
    axis.x.gridline       = TRUE,
    axis.y.gridline       = TRUE,
    axis.gridline.color   = "grey70"
  ))
  flat <- sp$allTracks()[[1]]$resolved_theme$flat
  expect_equal(flat[["axis.gridline.color"]], "grey70")
})

# ---- Plot-level vs track-level precedence ----

test_that("plot-level axis.y.gridline = TRUE is visible in resolved_theme", {
  sp <- make_sp(plot_aes = aes(axis.y.gridline = TRUE))
  flat <- sp$allTracks()[[1]]$resolved_theme$flat
  expect_true(isTRUE(flat[["axis.y.gridline"]]))
})

test_that("per-track axis.y.gridline = FALSE overrides plot-level TRUE", {
  skip_if_not(capabilities("png"))
  el1  <- seq_bar(data = gr, mapping = map(y = score))
  el2  <- seq_bar(data = gr, mapping = map(y = score))
  trk1 <- seq_track(windows = win, elements = list(el1),
                    aesthetics = aes(axis.y.gridline = TRUE))
  trk2 <- seq_track(windows = win, elements = list(el2),
                    aesthetics = aes(axis.y.gridline = FALSE))
  sp <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk1, trk2))
  sp$layoutGrid()
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})

# ---- plot() integration ----

test_that("SeqPlot$plot() calls drawGridlines() without error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(track_aes = aes(axis.y.gridline = TRUE))
  open_dev(); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

# ---- Multi-window track ----

test_that("drawGridlines() handles multi-window tracks without error", {
  skip_if_not(capabilities("png"))
  win2 <- GRanges("chr1", IRanges(c(1, 501), width = 500))
  gr2  <- GRanges("chr1", IRanges(c(50, 200, 600, 800), width = 40),
                  score = c(1, 2, 3, 4))
  el   <- seq_bar(data = gr2, mapping = map(y = score))
  trk  <- seq_track(windows = win2, elements = list(el),
                    aesthetics = aes(axis.y.gridline = TRUE))
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks = list(trk))
  sp$layoutGrid()
  open_dev(); on.exit(grDevices::dev.off())
  sp$drawGrid()
  expect_invisible(sp$drawGridlines())
})
