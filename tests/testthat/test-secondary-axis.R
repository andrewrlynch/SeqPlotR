library(testthat)
library(GenomicRanges)

.make_data <- function() {
  set.seed(1)
  GRanges("chr1", IRanges(start = seq(1, 1000, 50), width = 1),
          score = rnorm(20, mean = 0, sd = 1),
          depth = rpois(20, lambda = 100))
}

test_that("element resolve captures axis_x / axis_y as scalar integers", {
  d <- .make_data()
  e <- seq_point(data = d, mapping = map(x = start, y = score))
  e$resolve()
  expect_identical(e$resolved$axis_x, 1L)
  expect_identical(e$resolved$axis_y, 1L)

  e2 <- seq_line(data = d, mapping = map(x = start, y = depth, axis.y = 2))
  e2$resolve()
  expect_identical(e2$resolved$axis_y, 2L)
  expect_identical(e2$resolved$axis_x, 1L)
})

test_that("invalid axis selectors error", {
  d <- .make_data()
  e <- seq_point(data = d, mapping = map(x = start, y = score, axis.y = 3))
  expect_error(e$resolve(), "axis.y")
})

test_that("secondary y scale is inferred independently from primary", {
  d <- .make_data()
  g <- GRanges("chr1", IRanges(1, 1000))
  p <- seq_plot() %+%
    seq_track(data = d, windows = g) %+%
    seq_point(map(x = start, y = score)) %+%
    seq_line(map(x = start, y = depth, axis.y = 2))
  # layoutGrid without drawing
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_false(is.null(trk$scale_y))
  expect_false(is.null(trk$scale_y2))
  # Score range (centered on 0) vs depth range (around 100) differ.
  expect_true(trk$has_axis_y2)
  expect_gt(mean(trk$scale_y2$limits),
            mean(trk$scale_y$limits))
})

test_that("has_axis_y2 is FALSE when no secondary elements exist", {
  d <- .make_data()
  g <- GRanges("chr1", IRanges(1, 1000))
  p <- seq_plot() %+%
    seq_track(data = d, windows = g) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_false(trk$has_axis_y2)
})

test_that("has_axis_y2 is TRUE when theme toggles visibility", {
  d <- .make_data()
  g <- GRanges("chr1", IRanges(1, 1000))
  p <- seq_plot() %+%
    seq_track(data = d, windows = g,
              scale_y2 = seq_scale_continuous(limits = c(0, 10)),
              aesthetics = aes(axis.y2.visible = TRUE)) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_true(trk$has_axis_y2)
})

test_that(".panels_for_element swaps xscale/yscale when selectors are 2", {
  panel <- list(inner = list(x0=0, x1=1, y0=0, y1=1),
                xscale  = c(0, 100),
                yscale  = c(0, 1),
                xscale2 = c(0, 1000),
                yscale2 = c(-10, 10))
  e_pri <- list(axis_x = 1L, axis_y = 1L)
  e_sec <- list(axis_x = 2L, axis_y = 2L)
  out_pri <- .panels_for_element(list(panel), e_pri)[[1]]
  out_sec <- .panels_for_element(list(panel), e_sec)[[1]]
  expect_equal(out_pri$xscale, c(0, 100))
  expect_equal(out_pri$yscale, c(0, 1))
  expect_equal(out_sec$xscale, c(0, 1000))
  expect_equal(out_sec$yscale, c(-10, 10))
})
