library(testthat)
library(GenomicRanges)

# ── Fixtures ──────────────────────────────────────────────────────────────────
# A two-window track: window 1 holds a single gene; window 2 holds four
# mutually-overlapping genes. Used to exercise the shared-row-height fix and the
# packing modes.

gr_pack <- GRanges(
  "chr1",
  IRanges(
    start = c(100,   5100, 5150, 5200, 5250),
    end   = c(300,   5400, 5450, 5500, 5550)
  ),
  gene_id    = c("A", "B", "C", "D", "E"),
  gene_name  = c("A", "B", "C", "D", "E"),
  strand_col = rep("+", 5),
  feat_type  = rep("exon", 5)
)

.pm <- function(x0, x1) list(
  xscale        = c(x0, x1),
  yscale        = c(0, 4),
  inner         = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
  track_data    = gr_pack,
  track_mapping = map(group = gene_id, strand = strand_col,
                      label = gene_name, type = feat_type)
)

two_win_layout <- list(.pm(1, 1000), .pm(5000, 6000))
two_win        <- GRanges("chr1", IRanges(c(1, 5000), c(1000, 6000)))

win2_layout    <- list(.pm(5000, 6000))
win2           <- GRanges("chr1", IRanges(5000, 6000))

gene_map <- function()
  map(group = gene_id, strand = strand_col, label = gene_name, type = feat_type)

# ── Constructor / match.arg ───────────────────────────────────────────────────

test_that("new packing args default correctly", {
  el <- seq_gene(gr_pack, gene_map())
  expect_equal(el$label_position, "flush")
  expect_equal(el$pack_direction, "bottom_to_top")
  expect_false(el$expanded)
  expect_equal(el$label_row_frac, 0.35)
})

test_that("invalid label_position / pack_direction error via match.arg", {
  expect_error(seq_gene(gr_pack, gene_map(), label_position = "sideways"), "arg")
  expect_error(seq_gene(gr_pack, gene_map(), pack_direction = "left_to_right"), "arg")
})

test_that("label_position / pack_direction / expanded stored", {
  el <- seq_gene(gr_pack, gene_map(),
                 label_position = "above",
                 pack_direction = "top_to_bottom",
                 expanded       = TRUE)
  expect_equal(el$label_position, "above")
  expect_equal(el$pack_direction, "top_to_bottom")
  expect_true(el$expanded)
})

# ── Shared row height across windows (core structural fix) ─────────────────────

test_that("exon/row height is identical across a sparse and a dense window", {
  el <- seq_gene(gr_pack, gene_map())
  el$prep(two_win_layout, two_win)
  df <- el$coordCanvas
  h  <- df$exon_y1 - df$exon_y0
  # A single glyph thickness across both windows — no stretch-to-fill.
  expect_equal(length(unique(round(h, 8))), 1L)
})

# ── expanded = TRUE ───────────────────────────────────────────────────────────

test_that("expanded = TRUE gives one tier per gene regardless of overlap", {
  # Three mutually overlapping genes in their own window.
  gr3 <- GRanges("chr1", IRanges(c(100, 150, 200), c(400, 450, 500)),
                 gene_id = c("g1", "g2", "g3"),
                 gene_name = c("g1", "g2", "g3"),
                 strand_col = rep("+", 3), feat_type = rep("exon", 3))
  pm3 <- list(list(xscale = c(1, 1000), yscale = c(0, 3),
                   inner = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
                   track_data = gr3,
                   track_mapping = map(group = gene_id, strand = strand_col,
                                       label = gene_name, type = feat_type)))
  el <- seq_gene(gr3, map(group = gene_id, strand = strand_col,
                          label = gene_name, type = feat_type),
                 expanded = TRUE)
  el$prep(pm3, GRanges("chr1", IRanges(1, 1000)))
  ymids <- tapply(el$coordCanvas$ymid, el$coordCanvas$gene, function(x) x[1])
  expect_equal(length(unique(round(ymids, 8))), 3L)
})

test_that("expanded packs per window but keeps constant glyph height", {
  # A2: gene A is only in window 1; window 2 (B/C/D/E) must start its own
  # tiering at 1 rather than reserving A's row. Glyph height stays constant.
  el <- seq_gene(gr_pack, gene_map(), expanded = TRUE)
  el$prep(two_win_layout, two_win)
  df <- el$coordCanvas
  h <- df$exon_y1 - df$exon_y0
  expect_equal(length(unique(round(h, 8))), 1L)   # constant height preserved

  # Window 1 has one gene (A) at tier 1; window 2's lowest gene is also tier 1,
  # so both share the same bottom-row ymid — no empty reserved row.
  ya <- df$ymid[df$gene == "A"][1]
  w2_min <- min(df$ymid[df$gene %in% c("B", "C", "D", "E")])
  expect_equal(w2_min, ya, tolerance = 1e-8)
})

# ── pack_direction ────────────────────────────────────────────────────────────

test_that("pack_direction = top_to_bottom draws tier 1 nearest the top", {
  el <- seq_gene(gr_pack, gene_map(), pack_direction = "top_to_bottom")
  el$prep(win2_layout, win2)
  df <- el$coordCanvas
  ym <- tapply(df$ymid, df$gene, function(x) x[1])
  # B starts leftmost -> tier 1; must sit closest to inner$y1 (0.9).
  expect_gt(ym[["B"]], ym[["C"]])
  expect_gt(ym[["C"]], ym[["D"]])
})

test_that("bottom_to_top (default) draws tier 1 nearest the bottom", {
  el <- seq_gene(gr_pack, gene_map())
  el$prep(win2_layout, win2)
  df <- el$coordCanvas
  ym <- tapply(df$ymid, df$gene, function(x) x[1])
  expect_lt(ym[["B"]], ym[["C"]])
})

# ── label_position above/below ────────────────────────────────────────────────

test_that("above/below shrink the glyph height by label_row_frac", {
  elf <- seq_gene(gr_pack, gene_map())
  elf$prep(win2_layout, win2)
  ela <- seq_gene(gr_pack, gene_map(), label_position = "above")
  ela$prep(win2_layout, win2)
  hf <- (elf$coordCanvas$exon_y1 - elf$coordCanvas$exon_y0)[1]
  ha <- (ela$coordCanvas$exon_y1 - ela$coordCanvas$exon_y0)[1]
  expect_equal(ha / hf, 1 - 0.35, tolerance = 1e-6)
})

test_that("above/below populate a finite label_y; flush leaves it NA", {
  ela <- seq_gene(gr_pack, gene_map(), label_position = "below")
  ela$prep(win2_layout, win2)
  expect_true(all(is.finite(ela$coordCanvas$label_y)))

  elf <- seq_gene(gr_pack, gene_map())
  elf$prep(win2_layout, win2)
  expect_true(all(is.na(elf$coordCanvas$label_y)))
})

test_that("label_position modes prep() and draw() without error", {
  for (lp in c("flush", "above", "below")) {
    el <- seq_gene(gr_pack, gene_map(), label_position = lp)
    expect_silent(el$prep(win2_layout, win2))
    grDevices::pdf(file = NULL)
    expect_invisible(el$draw())
    grDevices::dev.off()
  }
})

test_that("expanded + pack_direction + above compose in draw()", {
  el <- seq_gene(gr_pack, gene_map(),
                 expanded = TRUE, pack_direction = "top_to_bottom",
                 label_position = "above")
  el$prep(two_win_layout, two_win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})
