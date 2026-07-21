library(testthat)
library(GenomicRanges)

win1 <- GRanges("chr1", IRanges(1, 1e6))

# Two tracks in one row (default direction) sit side by side; `track.gap.width`
# controls the npc gap between them, measured on the track bounds.
stacked <- function(gap_aes = list()) {
  t1 <- seq_track(windows = win1)
  t2 <- seq_track(windows = win1)
  sp <- SeqPlotR:::SeqPlotR6$new(tracks = list(t1, t2),
                                 aesthetics = do.call(aes, gap_aes))
  sp$layoutGrid()
  sp
}

track_gap <- function(sp) {
  tb <- sp$layout$trackBounds
  tb[[2]]$x0 - tb[[1]]$x1
}

test_that("default track gap is ~0.01 NPC", {
  expect_equal(track_gap(stacked()), 0.01, tolerance = 1e-6)
})

test_that("aes('track.gap.width' = 0) removes the gap", {
  expect_equal(track_gap(stacked(list("track.gap.width" = 0))),
               0, tolerance = 1e-6)
})

test_that("aes('track.gap.width' = 0.05) sets a 0.05 NPC gap", {
  expect_equal(track_gap(stacked(list("track.gap.width" = 0.05))),
               0.05, tolerance = 1e-6)
})

test_that("track_gaps alias still controls the gap", {
  expect_equal(track_gap(stacked(list(track_gaps = 0.04))),
               0.04, tolerance = 1e-6)
})

test_that("track.gap.width wins over the track_gaps alias", {
  expect_equal(track_gap(stacked(list("track.gap.width" = 0.03,
                                       track_gaps        = 0.08))),
               0.03, tolerance = 1e-6)
})

test_that("larger track gap shrinks the tracks proportionally", {
  width <- function(gap) {
    sp <- stacked(list("track.gap.width" = gap))
    tb <- sp$layout$trackBounds[[1]]
    tb$x1 - tb$x0
  }
  expect_lt(width(0.2), width(0))
})
