library(SeqPlotR)
library(GenomicRanges)

# ── Local helpers (tests are siloed; can't reuse those in test-links.R) ───────

hl_make_win <- function(start = 1L, end = 1000L) {
  GRanges("chr1", IRanges(start, end))
}

hl_make_two_track_layout <- function(win_a = hl_make_win(),
                                     win_b = hl_make_win()) {
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win_a) %__%
    seq_track(track_id = "B", windows = win_b)
  p$layoutGrid()
  list(plot = p, layout = p$layout)
}

hl_make_three_track_layout <- function(win = hl_make_win()) {
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win) %__%
    seq_track(track_id = "C", windows = win)
  p$layoutGrid()
  list(plot = p, layout = p$layout)
}

hl_named_layout <- function(tl) {
  setNames(tl$layout$panelBounds,
           vapply(tl$plot$allTracks(),
                  function(t) t$track_id, character(1)))
}

hl_named_windows <- function(tl) {
  setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
           vapply(tl$plot$allTracks(),
                  function(t) t$track_id, character(1)))
}

hl_make_region <- function() {
  GRanges("chr1", IRanges(c(200, 600), width = 100))
}

# ── Construction & inheritance ────────────────────────────────────────────────

test_that("seq_highlight instantiates", {
  gr <- hl_make_region()
  expect_no_error(
    seq_highlight(data = gr, map(x0 = start, x0_end = end),
                  t0 = "A", t1 = "B")
  )
})

test_that("seq_highlight inherits from SeqLink", {
  expect_true(inherits(seq_highlight(), "SeqLink"))
  expect_true(inherits(seq_highlight(), "SeqElement"))
})

# ── prep() ────────────────────────────────────────────────────────────────────

test_that("seq_highlight prep() runs without error on a 2-track layout", {
  tl <- hl_make_two_track_layout()
  gr <- hl_make_region()
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "A", t1 = "B")
  expect_no_error(
    lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  )
})

test_that("seq_highlight prep() populates coordCanvas with the right columns", {
  tl <- hl_make_two_track_layout()
  gr <- hl_make_region()
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "A", t1 = "B")
  lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  expect_false(is.null(lnk$coordCanvas))
  expect_true(all(c("region_id", "track_pos", "window_idx",
                    "xL", "xR", "y0_npc", "y1_npc") %in%
                  names(lnk$coordCanvas)))
  # Two regions, two tracks, one window each — expect 4 rows.
  expect_equal(nrow(lnk$coordCanvas), 4L)
  expect_setequal(unique(lnk$coordCanvas$track_pos), c(1L, 2L))
})

# ── draw() ────────────────────────────────────────────────────────────────────

test_that("seq_highlight draw() runs without error", {
  tl <- hl_make_two_track_layout()
  gr <- hl_make_region()
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "A", t1 = "B")
  lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(lnk$draw())
  dev.off()
})

# ── End-to-end via seq_plot ───────────────────────────────────────────────────

test_that("seq_highlight end-to-end across 3 tracks", {
  win <- hl_make_win()
  gr  <- hl_make_region()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win) %__%
    seq_track(track_id = "C", windows = win) %+%
    seq_highlight(data = gr, map(x0 = start, x0_end = end),
                  t0 = "A", t1 = "C")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── Single-track mode ─────────────────────────────────────────────────────────

test_that("seq_highlight with t0 == t1 highlights a single track only", {
  tl <- hl_make_two_track_layout()
  gr <- hl_make_region()
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "A", t1 = "A")
  lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  expect_false(is.null(lnk$coordCanvas))
  expect_setequal(unique(lnk$coordCanvas$track_pos), 1L)
})

test_that("seq_highlight with t1 = NULL defaults to t0 (single track)", {
  tl <- hl_make_two_track_layout()
  gr <- hl_make_region()
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "B")
  lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  expect_false(is.null(lnk$coordCanvas))
  expect_setequal(unique(lnk$coordCanvas$track_pos), 2L)
})

# ── Different scales: highlight expands / compresses ──────────────────────────

test_that("seq_highlight projects to different widths when track scales differ", {
  # Track A has window 1..1000, track B has window 1..200 (5x zoom).
  win_a <- hl_make_win(1L, 1000L)
  win_b <- hl_make_win(1L,  200L)
  tl <- hl_make_two_track_layout(win_a = win_a, win_b = win_b)

  # Region 100..200 — fully inside both windows.
  gr <- GRanges("chr1", IRanges(100, 200))
  lnk <- seq_highlight(data = gr, map(x0 = start, x0_end = end),
                       t0 = "A", t1 = "B")
  lnk$prep(hl_named_layout(tl), hl_named_windows(tl))
  df <- lnk$coordCanvas
  expect_equal(nrow(df), 2L)
  width_a <- df$xR[df$track_pos == 1L] - df$xL[df$track_pos == 1L]
  width_b <- df$xR[df$track_pos == 2L] - df$xL[df$track_pos == 2L]
  # Same panel inner width, but B covers a 5x narrower genomic span,
  # so the same 100bp region should be ~5x wider in B than in A.
  expect_gt(width_b, width_a * 3)
})
