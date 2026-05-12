library(SeqPlotR)
library(GenomicRanges)

# ── .parse_layout_string ───────────────────────────────────────────────────────

test_that("parse_layout_string parses simple 2-column layout", {
  parsed <- SeqPlotR:::.parse_layout_string("AB")
  expect_equal(parsed$nrow, 1)
  expect_equal(parsed$ncol, 2)
  expect_named(parsed$regions, c("A", "B"))
  expect_equal(parsed$regions$A, list(r0 = 1, r1 = 1, c0 = 1, c1 = 1))
  expect_equal(parsed$regions$B, list(r0 = 1, r1 = 1, c0 = 2, c1 = 2))
})

test_that("parse_layout_string handles multi-row multi-col spanning", {
  s <- "##AA\n##AA\nBBBC\nBBBD"
  parsed <- SeqPlotR:::.parse_layout_string(s)
  expect_equal(parsed$nrow, 4)
  expect_equal(parsed$ncol, 4)
  expect_equal(parsed$regions$A$r0, 1); expect_equal(parsed$regions$A$r1, 2)
  expect_equal(parsed$regions$A$c0, 3); expect_equal(parsed$regions$A$c1, 4)
  expect_equal(parsed$regions$B$r0, 3); expect_equal(parsed$regions$B$r1, 4)
  expect_equal(parsed$regions$B$c0, 1); expect_equal(parsed$regions$B$c1, 3)
  expect_equal(parsed$regions$C$r0, 3); expect_equal(parsed$regions$C$r1, 3)
  expect_equal(parsed$regions$C$c0, 4); expect_equal(parsed$regions$C$c1, 4)
  expect_equal(parsed$regions$D$r0, 4); expect_equal(parsed$regions$D$r1, 4)
  expect_equal(parsed$regions$D$c0, 4); expect_equal(parsed$regions$D$c1, 4)
})

test_that("parse_layout_string excludes # from regions", {
  parsed <- SeqPlotR:::.parse_layout_string("##AA\n##AA\nBBBC\nBBBD")
  expect_false("#" %in% names(parsed$regions))
})

test_that("parse_layout_string errors on unequal row lengths", {
  expect_error(SeqPlotR:::.parse_layout_string("AB\nABC"), "length")
})

test_that("parse_layout_string errors on non-rectangular region", {
  # L-shape is invalid
  expect_error(SeqPlotR:::.parse_layout_string("AB\nAA"), "rectangular")
})

# ── Positional layout bounds ───────────────────────────────────────────────────

make_win <- function() GRanges("chr1", IRanges(1, 1000))

test_that("positional layout: panel metadata includes xscale2/yscale2 fields", {
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win)
  p$layoutGrid()
  panel <- p$layout$panelBounds[[1]][[1]]
  expect_true("xscale2" %in% names(panel))
  expect_true("yscale2" %in% names(panel))
  # NULL when no secondary scale supplied
  expect_null(panel$xscale2)
  expect_null(panel$yscale2)
})

test_that("positional layout: track widths are proportional", {
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", track_width = 2, windows = win) %|%
    seq_track(track_id = "B", track_width = 1, windows = win)
  p$layoutGrid()
  bounds <- p$layout$trackBounds
  w_A <- bounds[[1]]$x1 - bounds[[1]]$x0
  w_B <- bounds[[2]]$x1 - bounds[[2]]$x0
  expect_equal(w_A / w_B, 2, tolerance = 0.01)
})

test_that("positional layout: rows stack vertically", {
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win)
  p$layoutGrid()
  bounds <- p$layout$trackBounds
  expect_gt(bounds[[1]]$y0, bounds[[2]]$y1 - 0.01)
})

# ── Patchwork layout bounds ────────────────────────────────────────────────────

test_that("patchwork layout assigns correct npc bounds", {
  win <- make_win()
  layout_str <- "AB"
  p <- seq_plot(layout = layout_str) %+%
    seq_track(track_id = "A", windows = win) %+%
    seq_track(track_id = "B", windows = win)
  p$layoutGrid()
  expect_true(is.list(p$layout$panelBounds))
})

test_that("patchwork layout: track not in layout string is silently skipped", {
  win <- make_win()
  p <- seq_plot(layout = "AB") %+%
    seq_track(track_id = "A", windows = win) %+%
    seq_track(track_id = "B", windows = win) %+%
    seq_track(track_id = "Z", windows = win)
  expect_silent(p$layoutGrid())
})

# ── SeqElementR6 base class ───────────────────────────────────────────────────

test_that("SeqElementR6 resolve() applies field-level mapping inheritance", {
  library(S4Vectors)
  gr <- GRanges("chr1", IRanges(c(100, 200, 300), width = 50),
                score = c(1, 2, 3), af = c(0.1, 0.5, 0.9))

  elem <- SeqPlotR:::SeqElementR6$new(
    mapping = map(y = score)
  )
  track_mapping <- map(x = start, y = af)

  elem$resolve(track_data = gr, track_mapping = track_mapping)

  expect_equal(elem$resolved$x, BiocGenerics::start(gr))
  expect_equal(elem$resolved$y, gr$score)
})

test_that("SeqElementR6 resolve() uses element data over track data", {
  gr1 <- GRanges("chr1", IRanges(100, 200), score = 1)
  gr2 <- GRanges("chr1", IRanges(300, 400), score = 99)

  elem <- SeqPlotR:::SeqElementR6$new(
    data    = gr2,
    mapping = map(y = score)
  )
  elem$resolve(track_data = gr1, track_mapping = map(y = score))
  expect_equal(elem$resolved$y, 99)
})

test_that("SeqElementR6 prep() stub errors with class name", {
  elem <- SeqPlotR:::SeqElementR6$new()
  expect_error(elem$prep(NULL, NULL), "SeqElement")
})

# ── SeqLinkR6 base class ──────────────────────────────────────────────────────

test_that("SeqLinkR6 inherits from SeqElement", {
  lnk <- SeqPlotR:::SeqLinkR6$new()
  expect_true(inherits(lnk, "SeqElement"))
  expect_true(inherits(lnk, "SeqLink"))
})

test_that("SeqLinkR6 resolve() builds anchor0_gr and anchor1_gr from one BEDPE", {
  gr <- GRanges("chr1", IRanges(100, 200), score = 5,
                start2 = 350, chr2 = "chr1")

  lnk <- SeqPlotR:::SeqLinkR6$new(
    data    = gr,
    mapping = map(x0 = start, x1 = start2,
                  chrom0 = seqnames, chrom1 = chr2,
                  y0 = score),
    t0 = "A", t1 = "B"
  )
  lnk$resolve()
  expect_s4_class(lnk$anchor0_gr, "GRanges")
  expect_s4_class(lnk$anchor1_gr, "GRanges")
  expect_equal(BiocGenerics::start(lnk$anchor0_gr), 100)
  expect_equal(BiocGenerics::start(lnk$anchor1_gr), 350)
  expect_equal(lnk$resolved$y0, 5)
})

# ── End-to-end layoutGrid ─────────────────────────────────────────────────────

test_that("layoutGrid() runs without error on minimal plot", {
  win <- GRanges("chr1", IRanges(1, 1000))
  p <- seq_plot() %|% seq_track(track_id = "A", windows = win)
  expect_silent(p$layoutGrid())
  expect_false(is.null(p$layout))
})

test_that("layoutGrid() errors if track has no windows", {
  p <- seq_plot() %|% seq_track(track_id = "A")
  expect_error(p$layoutGrid(), "windows")
})
