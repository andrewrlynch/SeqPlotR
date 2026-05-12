library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000000))
make_layout <- function() {
  trk <- seq_track(windows = win)
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  sp$layout$panelBounds[[1]]
}

test_that("seq_scalebar() errors on non-positive length_bp", {
  expect_error(seq_scalebar(-100), "positive")
  expect_error(seq_scalebar(0),    "positive")
})

test_that("auto label formats Mb correctly", {
  sb <- seq_scalebar(1e6)
  expect_equal(sb$label, "1 Mb")
})

test_that("auto label formats kb correctly", {
  sb <- seq_scalebar(50000)
  expect_equal(sb$label, "50 kb")
})

test_that("auto label formats bp correctly", {
  sb <- seq_scalebar(200)
  expect_equal(sb$label, "200 bp")
})

test_that("custom label overrides auto", {
  sb <- seq_scalebar(50000, label = "50 kilobases")
  expect_equal(sb$label, "50 kilobases")
})

test_that("prep() produces coordCanvas with correct length", {
  lt <- make_layout()
  sb <- seq_scalebar(50000)
  sb$prep(lt, win)
  expect_length(sb$coordCanvas, 1L)
})

test_that("prep() bar x0 < x1", {
  lt <- make_layout()
  sb <- seq_scalebar(50000)
  sb$prep(lt, win)
  cc <- sb$coordCanvas[[1]]
  expect_lt(cc$x0, cc$x1)
})

test_that("draw() runs without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  sb <- seq_scalebar(50000)
  sb$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sb$draw())
})

test_that("draw() with ticks = FALSE runs without error", {
  skip_if_not(capabilities("png"))
  lt <- make_layout()
  sb <- seq_scalebar(50000, ticks = FALSE)
  sb$prep(lt, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sb$draw())
})

test_that("seq_scalebar integrates into seq_plot$plot()", {
  skip_if_not(capabilities("png"))
  trk <- seq_track(windows = win) %+% seq_scalebar(50000)
  sp  <- seq_plot() %+% trk
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})
