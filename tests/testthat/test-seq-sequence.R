library(testthat)
library(SeqPlotR)
library(GenomicRanges)

win_short <- GRanges("chr1", IRanges(1, 100))   # 100 bp <= 200 — renders
win_long  <- GRanges("chr1", IRanges(1, 500))   # 500 bp > 200 — suppressed
win_mid   <- GRanges("chr1", IRanges(1, 60))    # 60 bp <= 80 — letters ok

# Build a layout_track list (list of per-window panel metadata) for a window.
# This mirrors what SeqPlotR6$layoutGrid() puts in layout$panelBounds[[i]].
make_layout <- function(win) {
  trk <- seq_track(windows = win)
  sp  <- seq_plot() %|% trk
  sp$layoutGrid()
  # panelBounds is keyed by track; first track = [[1]]
  # That value is the layout_track list passed to prep()
  sp$layout$panelBounds[[1]]
}

seq_str_100 <- paste(rep(c("A", "T", "C", "G"), 25), collapse = "")
seq_str_60  <- paste(rep(c("A", "T", "C", "G"), 15), collapse = "")

# ---- Construction ----

test_that("seq_sequence() errors without genome or sequence", {
  expect_error(seq_sequence(), "genome|sequence")
})

test_that("seq_sequence() accepts sequence string", {
  el <- seq_sequence(sequence = seq_str_100)
  expect_s3_class(el, "SeqSequence")
})

test_that("seq_sequence() accepts genome argument", {
  el <- seq_sequence(genome = "BSgenome.Hsapiens.UCSC.hg38")
  expect_equal(el$genome, "BSgenome.Hsapiens.UCSC.hg38")
})

test_that("colors defaults to UCSC standard", {
  el <- seq_sequence(sequence = "ATCG")
  expect_equal(el$colors["A"], c(A = "#00AA00"))
  expect_equal(el$colors["T"], c(T = "#FF0000"))
})

test_that("custom colors are stored", {
  custom <- c(A = "#111111", T = "#222222", C = "#333333", G = "#444444")
  el <- seq_sequence(sequence = "ATCG", colors = custom)
  expect_equal(el$colors["A"], c(A = "#111111"))
})

# ---- prep(): long window suppression ----

test_that("prep() emits message and produces empty coordCanvas for > 200 bp window", {
  el <- seq_sequence(sequence = paste(rep("A", 500), collapse = ""))
  layout_track <- make_layout(win_long)
  expect_message(el$prep(layout_track, win_long), "suppressed")
  expect_equal(nrow(el$coordCanvas), 0)
})

# ---- prep(): short window renders ----

test_that("prep() produces n rows equal to window width for short window", {
  el <- seq_sequence(sequence = seq_str_100)
  layout_track <- make_layout(win_short)
  el$prep(layout_track, win_short)
  expect_equal(nrow(el$coordCanvas), 100L)
})

test_that("prep() assigns correct colors to bases", {
  el <- seq_sequence(sequence = "ATCG")
  win4 <- GRanges("chr1", IRanges(1, 4))
  layout_track <- make_layout(win4)
  el$prep(layout_track, win4)
  df <- el$coordCanvas
  expect_equal(df$color[df$base == "A"][1], "#00AA00")
  expect_equal(df$color[df$base == "T"][1], "#FF0000")
  expect_equal(df$color[df$base == "C"][1], "#0000FF")
  expect_equal(df$color[df$base == "G"][1], "#FFB300")
})

test_that("prep() pads short sequence with N and warns", {
  el <- seq_sequence(sequence = "ATCG")   # 4 bp for a 10 bp window
  win10 <- GRanges("chr1", IRanges(1, 10))
  layout_track <- make_layout(win10)
  expect_warning(el$prep(layout_track, win10), "shorter")
  expect_equal(nrow(el$coordCanvas), 10L)
})

test_that("prep() stores win_width correctly", {
  el <- seq_sequence(sequence = seq_str_100)
  layout_track <- make_layout(win_short)
  el$prep(layout_track, win_short)
  expect_equal(unique(el$coordCanvas$win_width), 100L)
})

# ---- rect_height ----

test_that("rect_height = NULL uses full track height", {
  el <- seq_sequence(sequence = seq_str_100)
  layout_track <- make_layout(win_short)
  el$prep(layout_track, win_short)
  df <- el$coordCanvas
  p  <- layout_track[[1]]$inner
  expect_equal(df$y0[1], p$y0, tolerance = 1e-6)
  expect_equal(df$y1[1], p$y1, tolerance = 1e-6)
})

test_that("rect_height = 0.5 occupies half track height", {
  el <- seq_sequence(sequence = seq_str_100, rect_height = 0.5)
  layout_track <- make_layout(win_short)
  el$prep(layout_track, win_short)
  df  <- el$coordCanvas
  p   <- layout_track[[1]]$inner
  full_h  <- p$y1 - p$y0
  rect_h  <- df$y1[1] - df$y0[1]
  expect_equal(rect_h, full_h * 0.5, tolerance = 0.01)
})

# ---- show_letters ----

test_that("show_letters defaults to FALSE", {
  el <- seq_sequence(sequence = seq_str_100)
  expect_false(el$show_letters)
})

test_that("show_letters = TRUE stored", {
  el <- seq_sequence(sequence = seq_str_60, show_letters = TRUE)
  expect_true(el$show_letters)
})

# ---- draw() ----

test_that("draw() produces no output for empty coordCanvas", {
  el <- seq_sequence(sequence = paste(rep("A", 500), collapse = ""))
  layout_track <- make_layout(win_long)
  suppressMessages(el$prep(layout_track, win_long))
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() renders rectangles without error for <= 200 bp", {
  skip_if_not(capabilities("png"))
  el <- seq_sequence(sequence = seq_str_100)
  layout_track <- make_layout(win_short)
  el$prep(layout_track, win_short)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() renders letters without error for <= 80 bp and show_letters = TRUE", {
  skip_if_not(capabilities("png"))
  el <- seq_sequence(sequence = seq_str_60, show_letters = TRUE)
  layout_track <- make_layout(win_mid)
  el$prep(layout_track, win_mid)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with background aes renders without error", {
  skip_if_not(capabilities("png"))
  el <- seq_sequence(sequence = seq_str_60,
                     show_letters = TRUE,
                     aesthetics   = aes(background = "white", color = "black"))
  layout_track <- make_layout(win_mid)
  el$prep(layout_track, win_mid)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

# ---- SeqPlot integration ----

test_that("seq_sequence integrates into seq_plot$plot() without error", {
  skip_if_not(capabilities("png"))
  el  <- seq_sequence(sequence = seq_str_100)
  trk <- seq_track(windows = win_short, elements = list(el))
  sp  <- seq_plot() %|% trk
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(sp$plot())
})
