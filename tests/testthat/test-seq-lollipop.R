library(testthat)
library(GenomicRanges)

.win <- GRanges("chr1", IRanges(1, 600))
.pm  <- function(gr) list(list(
  xscale = c(1, 600), yscale = c(-5, 5), yplot_range = c(-5, 5),
  inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
  track_data = gr, track_mapping = map(x = start, y = score)
))

baseline_npc <- 0.1 + (0 - (-5)) / 10 * 0.8   # baseline 0 mapped into [0.1,0.9]

# ── D1: negative values place the head at the value, not the baseline ─────────

test_that("negative value draws its head below the baseline", {
  gr <- GRanges("chr1", IRanges(c(100, 300), width = 1), score = c(3, -3))
  el <- seq_lollipop(gr, map(x = start, y = score))
  el$prep(.pm(gr), .win)
  cc <- el$coordCanvas[[1]]
  # order of idx matches data order
  expect_gt(cc$y_head[1], baseline_npc)   # +3 above
  expect_lt(cc$y_head[2], baseline_npc)   # -3 below (the bug drew this at baseline)
})

test_that("+x and -x place heads symmetrically about the baseline", {
  gr <- GRanges("chr1", IRanges(c(100, 300), width = 1), score = c(2, -2))
  el <- seq_lollipop(gr, map(x = start, y = score))
  el$prep(.pm(gr), .win)
  h <- el$coordCanvas[[1]]$y_head
  expect_equal(abs(h[1] - baseline_npc), abs(h[2] - baseline_npc), tolerance = 1e-8)
})

test_that("stem still spans baseline to value", {
  gr <- GRanges("chr1", IRanges(100, width = 1), score = -3)
  el <- seq_lollipop(gr, map(x = start, y = score))
  el$prep(.pm(gr), .win)
  cc <- el$coordCanvas[[1]]
  expect_equal(cc$y1, baseline_npc, tolerance = 1e-8)  # top of stem = baseline
  expect_lt(cc$y0, baseline_npc)                        # bottom of stem = value
  expect_equal(cc$y_head, cc$y0, tolerance = 1e-8)      # head at the value
})

# ── D2: aesthetic mapping ─────────────────────────────────────────────────────

test_that("map(color = grp) yields per-point colors and a legend", {
  gr <- GRanges("chr1", IRanges(c(100, 300, 500), width = 1),
                score = c(1, 2, 3), grp = c("a", "b", "a"))
  el <- seq_lollipop(gr, map(x = start, y = score, color = grp))
  el$prep(.pm(gr), .win)
  cc <- el$coordCanvas[[1]]
  expect_length(cc$color, 3L)
  expect_equal(length(unique(cc$color)), 2L)  # two groups -> two colors
  expect_false(is.null(el$auto_legend))
})

test_that("draw() runs with and without a color mapping", {
  gr <- GRanges("chr1", IRanges(c(100, 300), width = 1),
                score = c(2, -2), grp = c("a", "b"))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  el1 <- seq_lollipop(gr, map(x = start, y = score))
  el1$prep(.pm(gr), .win); expect_invisible(el1$draw())
  el2 <- seq_lollipop(gr, map(x = start, y = score, color = grp))
  el2$prep(.pm(gr), .win); expect_invisible(el2$draw())
})
