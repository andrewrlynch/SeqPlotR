test_that("%||% returns lhs when not NULL", {
  expect_equal("a" %||% "b", "a")
})

test_that("%||% returns rhs when lhs is NULL", {
  expect_equal(NULL %||% "fallback", "fallback")
})

test_that("flexoki_palette returns correct length", {
  expect_length(flexoki_palette(5), 5)
  expect_length(flexoki_palette(9), 9)
})

test_that("flexoki_palette interpolates beyond 9", {
  expect_length(flexoki_palette(15), 15)
})

test_that("flexoki_palette returns character hex colors", {
  pal <- flexoki_palette(3)
  expect_type(pal, "character")
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that(".stop_if_not_granges errors on non-GRanges", {
  expect_error(SeqPlotR:::.stop_if_not_granges("not_a_gr", "windows"), "GRanges")
})

test_that(".stop_if_not_granges passes on GRanges", {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 100))
  expect_silent(SeqPlotR:::.stop_if_not_granges(gr, "windows"))
})

# ── Genome window helpers ─────────────────────────────────────────────────────

test_that("default_genome_windows returns 24 chromosomes", {
  w <- default_genome_windows()
  expect_s4_class(w, "GRanges")
  expect_equal(length(w), 24L)
  expect_true(all(startsWith(as.character(GenomicRanges::seqnames(w)), "chr")))
})

test_that("default_genome_windows(add_chr=FALSE) strips prefix", {
  w <- default_genome_windows(add_chr = FALSE)
  expect_false(any(startsWith(as.character(GenomicRanges::seqnames(w)), "chr")))
})

test_that("create_genome_windows expands a bare chromosome", {
  w <- create_genome_windows("chr1")
  expect_equal(length(w), 1L)
  expect_equal(BiocGenerics::start(w), 1L)
  expect_equal(BiocGenerics::end(w), 248956422L)
})

test_that("create_genome_windows parses chr:start-end", {
  w <- create_genome_windows("chr2:1000000-2000000")
  expect_equal(BiocGenerics::start(w), 1000000L)
  expect_equal(BiocGenerics::end(w), 2000000L)
  expect_equal(as.character(GenomicRanges::seqnames(w)), "chr2")
})

test_that("create_genome_windows tolerates commas and missing prefix", {
  w <- create_genome_windows("2:1,000,000-2,000,000")
  expect_equal(BiocGenerics::start(w), 1000000L)
  expect_equal(BiocGenerics::end(w), 2000000L)
})

test_that("create_genome_windows applies padding and clips to 1", {
  w <- create_genome_windows("chr3:500-1000", padding = 1000)
  expect_equal(BiocGenerics::start(w), 1L)
  expect_equal(BiocGenerics::end(w), 2000L)
})

test_that("create_genome_windows merges overlapping ranges", {
  w <- create_genome_windows(c("chr1:100-500", "chr1:400-800"))
  expect_equal(length(w), 1L)
  expect_equal(BiocGenerics::start(w), 100L)
  expect_equal(BiocGenerics::end(w), 800L)
})

test_that("create_genome_windows errors on unknown chromosome", {
  expect_error(create_genome_windows("chrZ"), "Unknown chromosome")
})

# ── load_cytobands ────────────────────────────────────────────────────────────

test_that("load_cytobands() returns a GRanges by default", {
  cb <- load_cytobands()
  expect_s4_class(cb, "GRanges")
  expect_gt(length(cb), 1000L)
  expect_true("gieStain" %in% colnames(S4Vectors::mcols(cb)))
  expect_true("name"     %in% colnames(S4Vectors::mcols(cb)))
})

test_that("load_cytobands(as_granges=FALSE) returns the raw data frame", {
  cb <- load_cytobands(as_granges = FALSE)
  expect_s3_class(cb, "data.frame")
  expect_true(all(c("chrom", "chromStart", "chromEnd", "name", "gieStain")
                  %in% colnames(cb)))
})

# ── seq_ideogram ──────────────────────────────────────────────────────────────

test_that("seq_ideogram instantiates and inherits SeqElement", {
  expect_no_error(seq_ideogram())
  expect_true(inherits(seq_ideogram(), "SeqElement"))
})

test_that("seq_ideogram prep + draw on a real chromosome", {
  cb  <- load_cytobands()
  win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 2.5e8))
  p <- seq_plot() %|%
    seq_track(track_id = "Ideo", windows = win) %+%
    seq_ideogram(data = cb)
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_ideogram coordCanvas excludes acen bands", {
  cb  <- load_cytobands()
  win <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 2.5e8))
  ideo <- seq_ideogram(data = cb)
  p <- seq_plot() %|% seq_track(track_id = "Ideo", windows = win) %+% ideo
  p$layoutGrid()
  panels <- p$layout$panelBounds[[1]]
  ideo$prep(panels, win)
  non_acen <- ideo$coordCanvas[[1]]
  expect_s3_class(non_acen, "data.frame")
  expect_gt(nrow(non_acen), 0L)
  expect_false("#FF0000" %in% non_acen$fill)
  expect_length(ideo$centroPolys[[1]], 2L)
})

test_that(".ideogram_fill_colors maps gieStain codes correctly", {
  cols <- SeqPlotR:::.ideogram_fill_colors(
    c("gneg", "gpos100", "acen", "stalk", "gvar", "weird"))
  expect_equal(cols[[1]], "#FFFFFF")
  expect_equal(cols[[3]], "#FF0000")
  expect_equal(cols[[4]], "#7AC0CF")
  expect_equal(cols[[5]], "#CCF5FF")
  expect_equal(cols[[6]], "#CCCCCC")
  # gpos100 should be pure black
  expect_equal(cols[[2]], grDevices::grey(0))
})
