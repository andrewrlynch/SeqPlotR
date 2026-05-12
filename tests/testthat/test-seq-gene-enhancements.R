library(testthat)
library(GenomicRanges)

# Two-gene, two-exon dataset with strand info
gr_genes <- GRanges(
  "chr1",
  IRanges(
    start = c(100, 400, 1000, 1500),
    end   = c(300, 600, 1200, 1700)
  ),
  gene_id    = c("geneA", "geneA", "geneB", "geneB"),
  gene_name  = c("GeneA", "GeneA", "GeneB", "GeneB"),
  strand_col = c("+", "+", "-", "-"),
  feat_type  = c("exon", "exon", "exon", "exon")
)

win <- GRanges("chr1", IRanges(1, 2000))

# Construct a layout_track in the same form as make_gene_pm() used elsewhere
make_layout <- function(gr = gr_genes) {
  list(list(
    xscale        = c(1, 2000),
    yscale        = c(0, 4),
    inner         = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(group = gene_id, strand = strand_col,
                        label = gene_name, type = feat_type)
  ))
}

layout_track <- make_layout()

# ---- backbone_type ----

test_that("backbone_type defaults to 'arrow'", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col))
  expect_equal(el$backbone_type, "arrow")
})

test_that("backbone_type = 'solid' accepted", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 backbone_type = "solid")
  expect_equal(el$backbone_type, "solid")
})

test_that("backbone_type = 'dashed' accepted", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 backbone_type = "dashed")
  expect_equal(el$backbone_type, "dashed")
})

test_that("invalid backbone_type errors", {
  expect_error(
    seq_gene(gr_genes, backbone_type = "wavy"),
    "arg"
  )
})

test_that("prep() runs for all backbone_type values", {
  for (bt in c("arrow", "solid", "dashed")) {
    el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                   backbone_type = bt)
    expect_silent(el$prep(layout_track, win))
    expect_true(nrow(el$coordCanvas) > 0)
  }
})

# ---- show_start ----

test_that("show_start defaults to FALSE", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col))
  expect_false(el$show_start)
})

test_that("show_start = TRUE accepted and stored", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 show_start = TRUE)
  expect_true(el$show_start)
})

test_that("prep() with show_start = TRUE produces tss_x column", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 show_start = TRUE)
  el$prep(layout_track, win)
  expect_true("tss_x" %in% names(el$coordCanvas))
  expect_true(all(is.finite(el$coordCanvas$tss_x)))
})

test_that("tss_position overrides auto TSS for named gene", {
  el <- seq_gene(gr_genes,
                 map(group = gene_id, strand = strand_col),
                 show_start   = TRUE,
                 tss_position = list(geneA = c(150, 200)))
  el$prep(layout_track, win)
  rows_a <- el$coordCanvas[el$coordCanvas$gene == "geneA", ]
  # TSS NPC for geneA should correspond to genomic 150, not 100
  lm <- layout_track[[1]]
  expected_npc <- (150 - lm$xscale[1]) / diff(lm$xscale) *
                  (lm$inner$x1 - lm$inner$x0) + lm$inner$x0
  expect_equal(rows_a$tss_x[1], expected_npc, tolerance = 0.001)
})

# ---- separate_strands ----

test_that("separate_strands defaults to FALSE", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col))
  expect_false(el$separate_strands)
})

test_that("separate_strands = TRUE accepted", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 separate_strands = TRUE)
  expect_true(el$separate_strands)
})

test_that("separate_strands = TRUE: plus-strand ymid > minus-strand ymid", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 separate_strands = TRUE)
  el$prep(layout_track, win)
  df <- el$coordCanvas
  plus_ymid  <- unique(df$ymid[df$strand == "+"])
  minus_ymid <- unique(df$ymid[df$strand == "-"])
  expect_true(all(plus_ymid > minus_ymid))
})

test_that("separate_strands silently ignored when only one strand present", {
  gr_plus <- gr_genes
  S4Vectors::mcols(gr_plus)$strand_col <- "+"
  lt <- make_layout(gr_plus)
  el <- seq_gene(gr_plus, map(group = gene_id, strand = strand_col),
                 separate_strands = TRUE)
  expect_silent(el$prep(lt, win))
  expect_true(nrow(el$coordCanvas) > 0)
})

# ---- draw() runs without error ----

test_that("draw() with backbone_type = 'arrow' runs without error", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 backbone_type = "arrow")
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with backbone_type = 'dashed' runs without error", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 backbone_type = "dashed")
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with show_start = TRUE runs without error", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 show_start = TRUE)
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with separate_strands = TRUE runs without error", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 separate_strands = TRUE)
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with all features enabled runs without error", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 backbone_type    = "solid",
                 show_start       = TRUE,
                 tss_position     = list(geneA = c(120, 160)),
                 separate_strands = TRUE)
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

# ---- style_type ----

test_that("style_type defaults to 'exon'", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col))
  expect_equal(el$style_type, "exon")
})

test_that("style_type = 'gene' accepted and stored", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 style_type = "gene")
  expect_equal(el$style_type, "gene")
})

test_that("style_type = 'point' accepted and stored", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
                 style_type = "point")
  expect_equal(el$style_type, "point")
})

test_that("style_type = 'invalid' errors via match.arg", {
  expect_error(
    seq_gene(gr_genes, map(group = gene_id, strand = strand_col),
             style_type = "invalid"),
    "arg"
  )
})

test_that("style_type = 'gene' emits one row per gene", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type = "gene")
  el$prep(layout_track, win)
  df <- el$coordCanvas
  expect_equal(nrow(df), length(unique(df$gene)))
  expect_equal(sort(unique(df$gene)), sort(unique(as.character(gr_genes$gene_id))))
  # chevron_npc column populated with finite, non-negative widths
  expect_true("chevron_npc" %in% names(df))
  expect_true(all(is.finite(df$chevron_npc)))
  expect_true(all(df$chevron_npc >= 0))
})

test_that("style_type = 'point' emits one row per gene with valid tss_x", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type = "point")
  el$prep(layout_track, win)
  df <- el$coordCanvas
  expect_equal(nrow(df), length(unique(df$gene)))
  expect_true(all(is.finite(df$tss_x)))
  # No exon boxes drawn in point mode
  expect_true(all(!df$draw_box))
})

test_that("draw() runs without error for style_type = 'exon'", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type = "exon")
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() runs without error for style_type = 'gene'", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type = "gene")
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() runs without error for style_type = 'point'", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type = "point")
  el$prep(layout_track, win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("style_type = 'gene' composes with separate_strands", {
  el <- seq_gene(gr_genes, map(group = gene_id, strand = strand_col,
                                label = gene_name, type = feat_type),
                 style_type       = "gene",
                 separate_strands = TRUE)
  el$prep(layout_track, win)
  df <- el$coordCanvas
  plus_ymid  <- unique(df$ymid[df$strand == "+"])
  minus_ymid <- unique(df$ymid[df$strand == "-"])
  expect_true(all(plus_ymid > minus_ymid))
})
