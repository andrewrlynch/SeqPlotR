library(testthat)
library(GenomicRanges)

gr  <- GRanges("chr1", IRanges(1, 100))
gr2 <- GRanges("chr1", IRanges(200, 300))

# ── LegendKey title field ──────────────────────────────────────────────────────

test_that("LegendKey stores title field", {
  k <- LegendKey(label = "Signal", title = "ChIP")
  expect_equal(k$title, "ChIP")
  expect_equal(k$label, "Signal")
})

test_that("LegendKey title defaults to NULL", {
  k <- LegendKey(label = "Signal")
  expect_null(k$title)
})

# ── SeqElement legend fields ───────────────────────────────────────────────────

test_that("SeqElement has legend = NULL and show_legend = TRUE by default", {
  el <- SeqPlotR:::SeqElementR6$new(gr)
  expect_null(el$legend)
  expect_true(el$show_legend)
})

test_that("collect_legend_keys returns NULL when legend is NULL", {
  el <- SeqPlotR:::SeqElementR6$new(gr)
  expect_null(el$collect_legend_keys())
})

test_that("collect_legend_keys returns NULL when show_legend = FALSE", {
  k  <- LegendKey(label = "A")
  el <- SeqPlotR:::SeqElementR6$new(gr, legend = k, show_legend = FALSE)
  expect_null(el$collect_legend_keys())
})

test_that("collect_legend_keys returns single entry for a LegendKey", {
  k  <- LegendKey(label = "Signal", title = "ChIP")
  el <- SeqPlotR:::SeqElementR6$new(gr, legend = k)
  res <- el$collect_legend_keys()
  expect_length(res, 1)
  expect_equal(res[[1]]$title, "ChIP")
  expect_equal(res[[1]]$key$label, "Signal")
  expect_true(inherits(res[[1]]$key, "LegendKey"))
})

test_that("collect_legend_keys returns multiple entries for a list of LegendKeys", {
  k1 <- LegendKey(label = "H3K27ac", title = "Marks")
  k2 <- LegendKey(label = "H3K4me3", title = "Marks")
  el <- SeqPlotR:::SeqElementR6$new(gr, legend = list(k1, k2))
  res <- el$collect_legend_keys()
  expect_length(res, 2)
  expect_equal(res[[1]]$key$label, "H3K27ac")
  expect_equal(res[[2]]$key$label, "H3K4me3")
})

test_that("collect_legend_keys errors on non-LegendKey list elements", {
  el <- SeqPlotR:::SeqElementR6$new(gr, legend = list(LegendKey(label = "ok"), "bad"))
  expect_error(el$collect_legend_keys(), "LegendKey")
})

test_that("element_class field in collect_legend_keys output is correct", {
  k  <- LegendKey(label = "x")
  el <- SeqPlotR:::SeqElementR6$new(gr, legend = k)
  res <- el$collect_legend_keys()
  expect_equal(res[[1]]$element_class, "SeqElement")
})

# ── SeqLink legend fields ──────────────────────────────────────────────────────

test_that("SeqLink inherits show_legend and legend fields", {
  lnk <- SeqPlotR:::SeqLinkR6$new(x1 = gr, x2 = gr2)
  expect_null(lnk$legend)
  expect_true(lnk$show_legend)
})

test_that("SeqLink collect_legend_keys works via inheritance", {
  k   <- LegendKey(label = "HiC", title = "Contacts")
  lnk <- SeqPlotR:::SeqLinkR6$new(x1 = gr, x2 = gr2, legend = k)
  res <- lnk$collect_legend_keys()
  expect_length(res, 1)
  expect_equal(res[[1]]$key$label, "HiC")
})

test_that("SeqLink collect_legend_keys returns NULL when show_legend = FALSE", {
  k   <- LegendKey(label = "HiC")
  lnk <- SeqPlotR:::SeqLinkR6$new(x1 = gr, x2 = gr2, legend = k, show_legend = FALSE)
  expect_null(lnk$collect_legend_keys())
})
