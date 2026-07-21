library(testthat)
library(GenomicRanges)

gene_map2 <- function()
  map(group = gene_id, strand = strand_col, label = gene_name, type = feat_type)

.mk_pm <- function(gr, x0, x1) list(
  xscale = c(x0, x1), yscale = c(0, 4),
  inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
  track_data = gr,
  track_mapping = gene_map2()
)

# ── A1: label-width packing (flush) ───────────────────────────────────────────

test_that("flush packing interval scales with label length", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  gr <- GRanges("chr1", IRanges(c(100, 100), c(200, 200)),
                gene_id = c("s", "l"),
                gene_name = c("X", paste(rep("Q", 20), collapse = "")),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  el <- seq_gene(gr, gene_map2(), label_pad = 0)
  el$prep(list(.mk_pm(gr, 1, 2000)), GRanges("chr1", IRanges(1, 2000)))
  pad <- tapply(el$coordCanvas$pad_bp, el$coordCanvas$gene, function(x) x[1])
  # A 20-char label must receive a larger pad than a 1-char label.
  expect_gt(pad[["l"]], pad[["s"]])
})

test_that("label-colliding, body-disjoint genes land on different tiers (flush)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  # Bodies are far apart, but a very long label on the first gene overlaps
  # the second gene's start once padding is applied.
  gr <- GRanges("chr1", IRanges(c(900, 1000), c(950, 1050)),
                gene_id = c("g1", "g2"),
                gene_name = c("A", "A"),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  el <- seq_gene(gr, gene_map2(), label_pad = 500)  # large pad forces overlap
  el$prep(list(.mk_pm(gr, 1, 2000)), GRanges("chr1", IRanges(1, 2000)))
  ym <- tapply(el$coordCanvas$ymid, el$coordCanvas$gene, function(x) x[1])
  expect_false(isTRUE(all.equal(ym[["g1"]], ym[["g2"]])))
})

test_that("above/below labels apply no horizontal packing padding", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  gr <- GRanges("chr1", IRanges(c(900, 1000), c(950, 1050)),
                gene_id = c("g1", "g2"), gene_name = c("A", "A"),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  el <- seq_gene(gr, gene_map2(), label_position = "above", label_pad = 500)
  el$prep(list(.mk_pm(gr, 1, 2000)), GRanges("chr1", IRanges(1, 2000)))
  # No horizontal label pad -> the two body-disjoint genes share a tier.
  expect_true(all(el$coordCanvas$pad_bp == 0))
  ym <- tapply(el$coordCanvas$ymid, el$coordCanvas$gene, function(x) x[1])
  expect_equal(ym[["g1"]], ym[["g2"]], tolerance = 1e-8)
})

# ── A2: expanded packs per window ─────────────────────────────────────────────

test_that("expanded tiering starts at 1 in each window independently", {
  gr <- GRanges("chr1",
                IRanges(c(100, 150, 200, 5100, 5150),
                        c(400, 450, 500, 5400, 5450)),
                gene_id = c("a", "b", "c", "d", "e"),
                gene_name = c("a", "b", "c", "d", "e"),
                strand_col = rep("+", 5), feat_type = rep("exon", 5))
  lt   <- list(.mk_pm(gr, 1, 1000), .mk_pm(gr, 5000, 6000))
  wins <- GRanges("chr1", IRanges(c(1, 5000), c(1000, 6000)))
  el <- seq_gene(gr, gene_map2(), expanded = TRUE)
  el$prep(lt, wins)
  df <- el$coordCanvas
  # Window 2 (genes d, e) must use two rows starting at the track floor, not
  # rows reserved above window 1's three genes.
  w2 <- df[df$gene %in% c("d", "e"), ]
  expect_equal(length(unique(round(w2$ymid, 6))), 2L)
  # Its lowest row equals window 1's lowest row (both tier 1).
  w1_min <- min(df$ymid[df$gene %in% c("a", "b", "c")])
  expect_equal(min(w2$ymid), w1_min, tolerance = 1e-8)
  # Glyph height still constant across windows.
  h <- df$exon_y1 - df$exon_y0
  expect_equal(length(unique(round(h, 8))), 1L)
})

# ── A4: glyph borders / linewidths ────────────────────────────────────────────

test_that("exon_lwd = NA keeps borderless boxes; a value strokes them", {
  gr <- GRanges("chr1", IRanges(c(100, 400), c(300, 600)),
                gene_id = c("g", "g"), gene_name = c("G", "G"),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  wins <- GRanges("chr1", IRanges(1, 1000))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  el_na <- seq_gene(gr, gene_map2(), exon_lwd = NA)
  el_na$prep(list(.mk_pm(gr, 1, 1000)), wins)
  expect_invisible(el_na$draw())
  el_b <- seq_gene(gr, gene_map2(), exon_lwd = 1.5, exon_border_col = "black",
                   backbone_lwd = 2, arrow_size = 2)
  el_b$prep(list(.mk_pm(gr, 1, 1000)), wins)
  expect_invisible(el_b$draw())
})

# ── A5: TSS marker is a filled triangle (grid.polygon), not a grid arrow ──────

test_that("show_start draws a filled triangle polygon", {
  gr <- GRanges("chr1", IRanges(c(100, 400), c(300, 600)),
                gene_id = c("g", "g"), gene_name = c("G", "G"),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  wins <- GRanges("chr1", IRanges(1, 1000))
  el <- seq_gene(gr, gene_map2(), show_start = TRUE,
                 tss_arrow_size = 4, tss_stem_lwd = 1.5)
  el$prep(list(.mk_pm(gr, 1, 1000)), wins)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  grid::grid.newpage(); grid::pushViewport(grid::viewport()); el$draw()
  g <- grid::grid.ls(print = FALSE)
  expect_gte(sum(grepl("polygon", g$name, ignore.case = TRUE)), 1L)
})

# ── A6: oob handling ──────────────────────────────────────────────────────────

test_that("oob 'remove' drops straddling genes; 'truncate' clips; 'keep' keeps", {
  gr <- GRanges("chr1", IRanges(c(400, 900), c(600, 1400)),
                gene_id = c("in", "straddle"),
                gene_name = c("in", "straddle"),
                strand_col = c("+", "+"), feat_type = c("exon", "exon"))
  win <- GRanges("chr1", IRanges(1, 1000))
  pm  <- list(.mk_pm(gr, 1, 1000))

  el_keep <- seq_gene(gr, gene_map2(), oob = "keep")
  el_keep$prep(pm, win)
  st_keep <- el_keep$coordCanvas[el_keep$coordCanvas$gene == "straddle", ]
  expect_gt(max(st_keep$x1b), 0.9)   # extends past the window's right edge

  el_tr <- seq_gene(gr, gene_map2(), oob = "truncate")
  el_tr$prep(pm, win)
  st_tr <- el_tr$coordCanvas[el_tr$coordCanvas$gene == "straddle", ]
  expect_equal(max(st_tr$x1b), 0.9, tolerance = 1e-6)  # clipped to window

  el_rm <- seq_gene(gr, gene_map2(), oob = "remove")
  el_rm$prep(pm, win)
  expect_false("straddle" %in% el_rm$coordCanvas$gene)  # dropped
  expect_true("in" %in% el_rm$coordCanvas$gene)
})

test_that("oob default is 'truncate'", {
  expect_equal(seq_gene()$oob, "truncate")
})
