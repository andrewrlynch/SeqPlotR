library(testthat)
library(GenomicRanges)

# Two equal-width 500 kb windows, abutting
win2 <- GRanges("chr1", IRanges(c(1, 5e5 + 1), width = 5e5))

make_sp <- function(gap_aes = list(), track_aes = list()) {
  trk <- seq_track(windows = win2, aesthetics = do.call(aes, track_aes))
  SeqPlotR:::SeqPlotR6$new(tracks  = list(trk),
                           aesthetics = do.call(aes, gap_aes))
}

# ── Default value ─────────────────────────────────────────────────────────────

test_that("default window.gap.width resolves to 0.01", {
  sp <- make_sp()
  sp$layoutGrid()
  # The canonical key isn't a default — only the alias is. The effective
  # value comes from the same chain the layout builder uses.
  effective <- sp$flat_theme[["window.gap.width"]] %||%
               sp$flat_theme$window_gaps           %||% 0.01
  expect_equal(effective, 0.01)
})

test_that("default gap between two windows is ~0.01 NPC", {
  sp <- make_sp()
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.01, tolerance = 1e-6)
})

# ── aes("window.gap.width" = x) at plot level ─────────────────────────────────

test_that("aes('window.gap.width' = 0) produces zero gap", {
  sp <- make_sp(gap_aes = list("window.gap.width" = 0))
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0, tolerance = 1e-6)
})

test_that("aes('window.gap.width' = 0.05) produces 0.05 NPC gap", {
  sp <- make_sp(gap_aes = list("window.gap.width" = 0.05))
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.05, tolerance = 1e-6)
})

# ── aes("window.gap.width" = x) at track level ────────────────────────────────

test_that("track-level aes overrides plot-level gap", {
  trk <- seq_track(windows = win2,
                   aesthetics = aes("window.gap.width" = 0.02))
  sp  <- SeqPlotR:::SeqPlotR6$new(tracks     = list(trk),
                                  aesthetics = aes("window.gap.width" = 0.05))
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.02, tolerance = 1e-6)
})

# ── Backward compatibility: window_gaps alias ─────────────────────────────────

test_that("window_gaps alias still controls gap", {
  sp <- make_sp(gap_aes = list(window_gaps = 0.02))
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.02, tolerance = 1e-6)
})

test_that("window.gap.width takes priority over window_gaps alias", {
  sp <- make_sp(gap_aes = list("window.gap.width" = 0.03,
                                window_gaps        = 0.07))
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[[1]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.03, tolerance = 1e-6)
})

# ── window_margin deprecation warning ────────────────────────────────────────

test_that("window_margin constructor arg emits deprecation warning", {
  expect_warning(
    seq_track(windows = win2, window_margin = 0.02),
    "deprecated"
  )
})

test_that("window_margin = NULL does not warn", {
  expect_no_warning(seq_track(windows = win2, window_margin = NULL))
})

# ── Single-window track: gap irrelevant ───────────────────────────────────────

test_that("single-window track unaffected by window.gap.width", {
  win1 <- GRanges("chr1", IRanges(1, 1e6))
  trk  <- seq_track(windows = win1)
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks     = list(trk),
                                   aesthetics = aes("window.gap.width" = 0))
  sp$layoutGrid()
  p <- sp$layout$panelBounds[[1]]
  expect_length(p, 1L)
})

# ── Multi-window: gap applied (nWin - 1) times ───────────────────────────────

test_that("three windows: two gaps of correct width", {
  win3 <- GRanges("chr1", IRanges(c(1, 4e5 + 1, 8e5 + 1), width = 4e5))
  trk  <- seq_track(windows = win3)
  sp   <- SeqPlotR:::SeqPlotR6$new(tracks     = list(trk),
                                   aesthetics = aes("window.gap.width" = 0.02))
  sp$layoutGrid()
  p    <- sp$layout$panelBounds[[1]]
  gap1 <- p[[2]]$full$x0 - p[[1]]$full$x1
  gap2 <- p[[3]]$full$x0 - p[[2]]$full$x1
  expect_equal(gap1, 0.02, tolerance = 1e-6)
  expect_equal(gap2, 0.02, tolerance = 1e-6)
})

# ── Patchwork layout also respects the key ────────────────────────────────────

test_that("patchwork layout respects window.gap.width", {
  win2b <- GRanges("chr1", IRanges(c(1, 5e5 + 1), width = 5e5))
  trk   <- seq_track(windows = win2b, track_id = "A")
  sp    <- SeqPlotR:::SeqPlotR6$new(
    tracks     = list(trk),
    layout     = "A",
    aesthetics = aes("window.gap.width" = 0.03)
  )
  sp$layoutGrid()
  p   <- sp$layout$panelBounds[["A"]]
  gap <- p[[2]]$full$x0 - p[[1]]$full$x1
  expect_equal(gap, 0.03, tolerance = 1e-6)
})
