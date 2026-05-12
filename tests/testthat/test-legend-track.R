library(testthat)
library(GenomicRanges)

gr  <- GRanges("chr1", IRanges(1, 100))
gr2 <- GRanges("chr1", IRanges(200, 300))
win <- GRanges("chr1", IRanges(1, 1000))

make_el <- function(label, title = NULL, show = TRUE) {
  k  <- LegendKey(label = label, title = title)
  SeqPlotR:::SeqElementR6$new(gr, legend = k, show_legend = show)
}

# ── show_legend field on SeqTrack ─────────────────────────────────────────────

test_that("SeqTrack has show_legend = TRUE by default", {
  trk <- seq_track(windows = win)
  expect_true(trk$show_legend)
})

test_that("SeqTrack show_legend = FALSE can be set at construction", {
  trk <- seq_track(windows = win, show_legend = FALSE)
  expect_false(trk$show_legend)
})

# ── collect_legend_keys ────────────────────────────────────────────────────────

test_that("collect_legend_keys returns NULL for empty track", {
  trk <- seq_track(windows = win)
  expect_null(trk$collect_legend_keys())
})

test_that("collect_legend_keys returns NULL when track show_legend = FALSE", {
  el  <- make_el("Signal")
  trk <- seq_track(windows = win, elements = list(el), show_legend = FALSE)
  expect_null(trk$collect_legend_keys())
})

test_that("collect_legend_keys returns NULL when all elements have show_legend = FALSE", {
  el  <- make_el("Signal", show = FALSE)
  trk <- seq_track(windows = win, elements = list(el))
  expect_null(trk$collect_legend_keys())
})

test_that("collect_legend_keys returns NULL when no element has a legend", {
  el  <- SeqPlotR:::SeqElementR6$new(gr)
  trk <- seq_track(windows = win, elements = list(el))
  expect_null(trk$collect_legend_keys())
})

test_that("collect_legend_keys returns keys from a single element", {
  el  <- make_el("Signal", title = "ChIP")
  trk <- seq_track(windows = win, elements = list(el))
  res <- trk$collect_legend_keys()
  expect_length(res, 1)
  expect_equal(res[[1]]$key$label, "Signal")
  expect_equal(res[[1]]$title, "ChIP")
})

test_that("collect_legend_keys aggregates keys from multiple elements", {
  el1 <- make_el("H3K27ac", title = "Marks")
  el2 <- make_el("H3K4me3", title = "Marks")
  trk <- seq_track(windows = win, elements = list(el1, el2))
  res <- trk$collect_legend_keys()
  expect_length(res, 2)
  labels <- vapply(res, function(x) x$key$label, character(1))
  expect_setequal(labels, c("H3K27ac", "H3K4me3"))
})

test_that("collect_legend_keys respects per-element show_legend = FALSE", {
  el1 <- make_el("Shown")
  el2 <- make_el("Hidden", show = FALSE)
  trk <- seq_track(windows = win, elements = list(el1, el2))
  res <- trk$collect_legend_keys()
  expect_length(res, 1)
  expect_equal(res[[1]]$key$label, "Shown")
})

test_that("collect_legend_keys works with SeqLink elements", {
  k   <- LegendKey(label = "HiC", title = "Contacts")
  lnk <- SeqPlotR:::SeqLinkR6$new(legend = k)
  trk <- seq_track(windows = win, elements = list(lnk))
  res <- trk$collect_legend_keys()
  expect_length(res, 1)
  expect_equal(res[[1]]$key$label, "HiC")
})
