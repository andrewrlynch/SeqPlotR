library(testthat)
library(GenomicRanges)

.win <- function() GRanges("chr1", IRanges(1, 1000))
.data <- function() {
  set.seed(1)
  GRanges("chr1", IRanges(seq(1, 1000, 50), width = 1),
          score = rnorm(20), depth = rpois(20, 100))
}

test_that("plot pipeline renders with just a primary axis", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win()) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(p$plot())
})

test_that("plot pipeline renders with a secondary y axis", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win()) %+%
    seq_point(map(x = start, y = score)) %+%
    seq_line(map(x = start, y = depth, axis.y = 2))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(p$plot())
})

test_that("plot pipeline renders with stacked x1/x2 on same side", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.x1.position = "top",
                               axis.x2.position = "top",
                               axis.x2.visible = TRUE),
              scale_x2 = seq_scale_continuous(limits = c(0, 1))) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(p$plot())
})

test_that("track chrome fill flows into resolved theme", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(track.background.fill = "lavender")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$resolved_theme$chrome$background$fill, "lavender")
})

test_that("axis line color propagates via inheritance", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.line.col = "grey30")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$resolved_theme$axes$x1$line$col, "grey30")
  expect_identical(trk$resolved_theme$axes$y1$line$col, "grey30")
})

test_that("more specific axis override wins", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.line.col = "black",
                               axis.y1.line.col = "red")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$resolved_theme$axes$y1$line$col, "red")
  expect_identical(trk$resolved_theme$axes$x1$line$col, "black")
})

test_that("axis.x1.scale.cap affects axis_range", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.x1.scale.cap = "ticks")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$resolved_theme$axes$x1$scale$cap, "ticks")
})
