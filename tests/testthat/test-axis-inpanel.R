library(testthat)
library(GenomicRanges)

win <- GRanges("chr1", IRanges(1, 1000))
gr  <- GRanges("chr1", IRanges(c(100, 400, 700), width = 80),
               score = c(1, 3, 2))

make_sp <- function(track_aes = list()) {
  trk <- seq_track(data = gr,
                   windows = win,
                   aesthetics = if (length(track_aes)) do.call(aes, track_aes)
                                else aes()) %+%
         seq_bar(map(x = start, y = score))
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  sp
}

test_that("default title position is 'axis' — no error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp()
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

test_that("in-panel y title renders without error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(list("axis.y.title.position" = c(0.02, 0.95),
                     "axis.y.title.text"     = "Score"))
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

test_that("range-style y labels render without error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(list("axis.y.labels.style" = "range"))
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})

test_that("range-style y labels with in-panel position render without error", {
  skip_if_not(capabilities("png"))
  sp <- make_sp(list("axis.y.labels.style"    = "range",
                     "axis.y.labels.position" = c(0.02, 0.02)))
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})
