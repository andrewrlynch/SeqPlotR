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

# ── Tick / label direction (Part B) ───────────────────────────────────────────

test_that(".axis_tick_dir defaults reproduce the historical hardcoded signs", {
  # These are the exact signs the old identifier-keyed code produced.
  expect_equal(SeqPlotR:::.axis_tick_dir("bottom", "out"), -1)  # was x1
  expect_equal(SeqPlotR:::.axis_tick_dir("top",    "out"),  1)  # was x2
  expect_equal(SeqPlotR:::.axis_tick_dir("left",   "out"), -1)  # was y1
  expect_equal(SeqPlotR:::.axis_tick_dir("right",  "out"),  1)  # was y2
})

test_that(".axis_tick_dir flips with direction = 'in'", {
  for (pos in c("bottom", "top", "left", "right")) {
    expect_equal(SeqPlotR:::.axis_tick_dir(pos, "in"),
                 -SeqPlotR:::.axis_tick_dir(pos, "out"))
  }
})

test_that(".axis_label_just tracks direction and matches historical defaults", {
  # x axes carry the vertical justification.
  expect_equal(SeqPlotR:::.axis_label_just("bottom", "out")$v, "top")
  expect_equal(SeqPlotR:::.axis_label_just("top",    "out")$v, "bottom")
  # y axes carry the horizontal justification.
  expect_equal(SeqPlotR:::.axis_label_just("left",  "out")$h, "right")
  expect_equal(SeqPlotR:::.axis_label_just("right", "out")$h, "left")
  # "in" flips the outward/inward side.
  expect_equal(SeqPlotR:::.axis_label_just("bottom", "in")$v, "bottom")
  expect_equal(SeqPlotR:::.axis_label_just("left",   "in")$h, "left")
})

test_that("axis.ticks.direction resolves through the axis inheritance chain", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.x1.ticks.direction = "in")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  expect_identical(trk$resolved_theme$axes$x1$ticks$direction, "in")
  # Untouched axes keep the "out" default.
  expect_identical(trk$resolved_theme$axes$y1$ticks$direction, "out")
})

test_that("default theme carries axis.ticks.direction = 'out'", {
  p <- seq_plot() %+%
    seq_track(data = .data(), windows = .win()) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  p$layoutGrid()
  trk <- p$allTracks()[[1]]
  for (side in c("x1", "x2", "y1", "y2"))
    expect_identical(trk$resolved_theme$axes[[side]]$ticks$direction, "out")
})

test_that("position and direction overrides render without error", {
  p1 <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.x1.position = "top")) %+%
    seq_point(map(x = start, y = score))
  p2 <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.y1.position = "right")) %+%
    seq_point(map(x = start, y = score))
  p3 <- seq_plot() %+%
    seq_track(data = .data(), windows = .win(),
              aesthetics = aes(axis.x1.ticks.direction = "in")) %+%
    seq_point(map(x = start, y = score))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(p1$plot())
  expect_no_error(p2$plot())
  expect_no_error(p3$plot())
})
