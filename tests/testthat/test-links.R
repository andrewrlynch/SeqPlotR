library(SeqPlotR)
library(GenomicRanges)

# ── Shared helpers ─────────────────────────────────────────────────────────────

make_win <- function() GRanges("chr1", IRanges(1, 1000))

# Single-locus GRanges: anchor 0 = start, anchor 1 = end (via map(x1 = end)).
make_link_gr <- function(n = 4) {
  GRanges("chr1",
          IRanges(start = c(100, 300, 500, 700)[seq_len(n)], width = 10),
          score = c(0.2, 0.8, 0.5, 0.3)[seq_len(n)])
}

# BEDPE-like data.frame for the cross-locus tests.
make_bedpe_df <- function() {
  data.frame(
    chr1   = "chr1",
    start1 = c(100, 300, 500, 700),
    chr2   = "chr1",
    start2 = c(150, 400, 600, 800),
    score  = c(0.2, 0.8, 0.5, 0.3),
    strand1 = c("+", "-", "+", "-"),
    strand2 = c("+", "-", "-", "+"),
    stringsAsFactors = FALSE
  )
}

# Minimal two-track layout for link testing
make_two_track_layout <- function() {
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win)
  p$layoutGrid()
  list(plot = p, layout = p$layout)
}

named_layout <- function(tl) {
  setNames(tl$layout$panelBounds,
           vapply(tl$plot$allTracks(),
                  function(t) t$track_id, character(1)))
}

named_windows <- function(tl) {
  setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
           vapply(tl$plot$allTracks(),
                  function(t) t$track_id, character(1)))
}

# ── drawSeqArch helper ────────────────────────────────────────────────────────

test_that("drawSeqArch runs without error", {
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(
    SeqPlotR:::drawSeqArch(x0 = 0.2, y0 = 0.1, x1 = 0.8, y1 = 0.1,
                           top0 = 0.5, top1 = 0.5,
                           orientation = "+", curve = "length",
                           stemWidth = 1, arcWidth = 1,
                           arcColor = "black", stemColor = "black")
  )
  dev.off()
})

# ── seq_arc ───────────────────────────────────────────────────────────────────

test_that("seq_arc instantiates", {
  expect_no_error(seq_arc())
  expect_no_error(seq_arc(map(x0 = start, x1 = end, y0 = score, height = score)))
})

test_that("seq_arc inherits from SeqLink and SeqElement", {
  arc <- seq_arc()
  expect_true(inherits(arc, "SeqLink"))
  expect_true(inherits(arc, "SeqElement"))
})

test_that("seq_arc t0/t1 are locked when added inside a seq_track", {
  arc <- seq_arc(map(x0 = start, x1 = end))
  trk <- seq_track(track_id = "A", windows = make_win())
  trk %+% arc
  expect_equal(arc$t0, "A")
  expect_equal(arc$t1, "A")
})

test_that("seq_arc resolve() builds anchor0_gr and anchor1_gr", {
  gr <- make_link_gr()
  arc <- seq_arc(data = gr,
                 mapping = map(x0 = start, x1 = end,
                               y0 = score, height = score))
  arc$resolve()
  expect_s4_class(arc$anchor0_gr, "GRanges")
  expect_s4_class(arc$anchor1_gr, "GRanges")
  expect_equal(length(arc$anchor0_gr), length(gr))
  expect_equal(BiocGenerics::start(arc$anchor0_gr), BiocGenerics::start(gr))
  expect_equal(BiocGenerics::start(arc$anchor1_gr), BiocGenerics::end(gr))
})

test_that("seq_arc prep() populates coordCanvas", {
  tl <- make_two_track_layout()
  gr <- make_link_gr()
  arc <- seq_arc(data = gr,
                 mapping = map(x0 = start, x1 = end,
                               y0 = score, height = score),
                 t0 = "A", t1 = "A")
  expect_no_error(arc$prep(named_layout(tl), named_windows(tl), "A"))
  expect_true(length(arc$coordCanvas) > 0)
})

test_that("seq_arc draw() runs without error", {
  tl <- make_two_track_layout()
  gr <- make_link_gr()
  arc <- seq_arc(data = gr,
                 mapping = map(x0 = start, x1 = end,
                               y0 = score, height = score),
                 t0 = "A", t1 = "A")
  arc$prep(named_layout(tl), named_windows(tl), "A")
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(arc$draw())
  dev.off()
})

test_that("seq_arc end-to-end via seq_plot", {
  gr  <- make_link_gr()
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", data = gr, mapping = map(x = start, y = score),
              windows = win) %+%
    seq_arc(map(x0 = start, x1 = end, y0 = score, height = score))
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── BEDPE data.frame input ────────────────────────────────────────────────────

test_that("seq_arch resolves a BEDPE-like data.frame", {
  bedpe <- make_bedpe_df()
  arch  <- seq_arch(data = bedpe,
                    mapping = map(x0 = start1, x1 = start2,
                                  chrom0 = chr1, chrom1 = chr2,
                                  y0 = score))
  arch$resolve()
  expect_s4_class(arch$anchor0_gr, "GRanges")
  expect_equal(length(arch$anchor0_gr), nrow(bedpe))
  expect_equal(BiocGenerics::start(arch$anchor0_gr), bedpe$start1)
  expect_equal(BiocGenerics::start(arch$anchor1_gr), bedpe$start2)
  expect_equal(as.character(GenomicRanges::seqnames(arch$anchor1_gr)),
               bedpe$chr2)
})

test_that("seq_link with data.frame errors when chrom0 is missing", {
  bedpe <- make_bedpe_df()
  arch  <- seq_arch(data = bedpe,
                    mapping = map(x0 = start1, x1 = start2))
  expect_error(arch$resolve(),
               "requires chrom0 in map\\(\\)")
})

test_that("seq_link errors when x0 or x1 is missing", {
  arc <- seq_arc(data = make_link_gr(),
                 mapping = map(x0 = start))
  expect_error(arc$resolve(), "must define both x0 and x1")
})

# ── seq_arch ──────────────────────────────────────────────────────────────────

test_that("seq_arch instantiates", {
  expect_no_error(seq_arch())
})

test_that("seq_arch inherits from SeqLink", {
  expect_true(inherits(seq_arch(), "SeqLink"))
})

test_that("seq_arch prep() populates coordCanvas", {
  tl <- make_two_track_layout()
  gr <- make_link_gr()
  arch <- seq_arch(data = gr,
                   mapping = map(x0 = start, x1 = end,
                                 y0 = score, height = score),
                   t0 = "A", t1 = "A")
  expect_no_error(arch$prep(named_layout(tl), named_windows(tl), "A"))
})

test_that("seq_arch draw() runs without error", {
  tl <- make_two_track_layout()
  gr <- make_link_gr()
  arch <- seq_arch(data = gr,
                   mapping = map(x0 = start, x1 = end,
                                 y0 = score, height = score),
                   t0 = "A", t1 = "A")
  arch$prep(named_layout(tl), named_windows(tl), "A")
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(arch$draw())
  dev.off()
})

# ── seq_recon ─────────────────────────────────────────────────────────────────

# Single GRanges with bedpe-like mcols carrying anchor-1 strand + position.
make_sv_gr <- function() {
  GRanges("chr1", IRanges(c(100, 300, 500, 700, 900), width = 10),
          strand   = c("+", "-", "+", "-", "+"),
          chr2     = c("chr1", "chr1", "chr1", "chr1", "chr2"),
          start2   = c(400, 600, 200, 800, 100),
          strand2  = c("+", "-", "-", "+", "+"))
}

recon_mapping <- function() {
  map(x0 = start,    x1 = start2,
      chrom0 = seqnames, chrom1 = chr2,
      strand0 = strand, strand1 = strand2)
}

test_that("seq_recon instantiates", {
  expect_no_error(seq_recon(data = make_sv_gr(), mapping = recon_mapping()))
})

test_that("seq_recon inherits from SeqArch and SeqLink", {
  rc <- seq_recon(data = make_sv_gr(), mapping = recon_mapping())
  expect_true(inherits(rc, "SeqArch"))
  expect_true(inherits(rc, "SeqLink"))
})

test_that("seq_recon errors when strand0 or strand1 is absent from map()", {
  rc <- seq_recon(data = make_sv_gr(),
                  mapping = map(x0 = start, x1 = start2,
                                chrom0 = seqnames, chrom1 = chr2),
                  t0 = "A", t1 = "A")
  tl <- make_two_track_layout()
  expect_error(rc$prep(named_layout(tl), named_windows(tl), "A"),
               "strand0, strand1")
})

run_recon_one <- function(strand0, strand1, chr2 = "chr1") {
  gr <- GRanges("chr1", IRanges(100, width = 10),
                strand  = strand0,
                chr2    = chr2,
                start2  = 500,
                strand2 = strand1)
  rc <- seq_recon(data = gr, mapping = recon_mapping(),
                  t0 = "A", t1 = "A")
  tl <- make_two_track_layout()
  rc$prep(named_layout(tl), named_windows(tl), "A")
  rc
}

test_that("seq_recon classifies H2H correctly (+/+)", {
  rc <- run_recon_one("+", "+")
  expect_equal(rc$aesthetics$arcColor[1], rc$col_h2h)
})

test_that("seq_recon classifies T2T correctly (-/-)", {
  rc <- run_recon_one("-", "-")
  expect_equal(rc$aesthetics$arcColor[1], rc$col_t2t)
})

test_that("seq_recon classifies Dup correctly (-/+)", {
  rc <- run_recon_one("-", "+")
  expect_equal(rc$aesthetics$arcColor[1], rc$col_dup)
})

test_that("seq_recon classifies Del correctly (+/-)", {
  rc <- run_recon_one("+", "-")
  expect_equal(rc$aesthetics$arcColor[1], rc$col_del)
})

test_that("seq_recon classifies Translocation correctly (different chr)", {
  rc <- run_recon_one("+", "+", chr2 = "chr2")
  expect_equal(rc$aesthetics$arcColor[1], rc$col_trans)
})

test_that("seq_recon draw() runs without error", {
  rc <- seq_recon(data = make_sv_gr(), mapping = recon_mapping(),
                  t0 = "A", t1 = "A")
  tl <- make_two_track_layout()
  rc$prep(named_layout(tl), named_windows(tl), "A")
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(rc$draw())
  dev.off()
})

# ── Shared helpers (Batch 5B) ─────────────────────────────────────────────────

make_string_data <- function() {
  # Single data object with both anchor positions
  GRanges("chr1",
          IRanges(c(100, 400), width = 1),
          x1pos  = c(500L, 700L),
          chrom1 = c("chr1", "chr1"),
          s0     = c("+", "-"),
          s1     = c("-", "+"),
          score  = c(0.8, 0.6))
}

make_synteny_data <- function() {
  data.frame(
    c0    = "chr1", p0    = c(100L, 400L), p0end = c(200L, 500L),
    c1    = "chr1", p1    = c(150L, 450L), p1end = c(250L, 550L),
    score = c(0.8, 0.6),
    stringsAsFactors = FALSE
  )
}

# ── seq_string ────────────────────────────────────────────────────────────────

test_that("seq_string instantiates", {
  gr <- make_string_data()
  expect_no_error(
    seq_string(data = gr,
               map(x0 = start, x1 = x1pos, chrom0 = seqnames, chrom1 = chrom1,
                   strand0 = s0, strand1 = s1, y0 = score, y1 = score),
               t0 = "A", t1 = "B")
  )
})

test_that("seq_string inherits from SeqLink", {
  expect_true(inherits(seq_string(), "SeqLink"))
})

test_that("seq_string t0/t1 locked when added inside seq_track", {
  lnk <- seq_string(data = make_string_data(),
                    map(x0 = start, x1 = x1pos))
  trk <- seq_track(track_id = "A", windows = make_win())
  trk %+% lnk
  expect_equal(lnk$t0, "A")
  expect_equal(lnk$t1, "A")
})

test_that("seq_string prep() runs without error cross-track", {
  tl  <- make_two_track_layout()
  gr  <- make_string_data()
  lnk <- seq_string(data = gr,
                    map(x0 = start, x1 = x1pos, chrom0 = seqnames, chrom1 = chrom1,
                        strand0 = s0, strand1 = s1, y0 = score, y1 = score),
                    t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  expect_no_error(lnk$prep(layout_named, windows_named))
})

test_that("seq_string draw() runs without error", {
  tl  <- make_two_track_layout()
  gr  <- make_string_data()
  lnk <- seq_string(data = gr,
                    map(x0 = start, x1 = x1pos, chrom0 = seqnames, chrom1 = chrom1,
                        strand0 = s0, strand1 = s1, y0 = score, y1 = score),
                    t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named)
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(lnk$draw())
  dev.off()
})

test_that("seq_string type auto-infers from strand (+/-) -> s curve", {
  gr <- GRanges("chr1", IRanges(100, width = 1),
                x1pos = 500L, chrom1 = "chr1", s0 = "+", s1 = "-", score = 0.8)
  lnk <- seq_string(data = gr,
                    map(x0 = start, x1 = x1pos, strand0 = s0, strand1 = s1),
                    t0 = "A", t1 = "A")
  tl <- make_two_track_layout()
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named, "A")
  expect_equal(lnk$coordCanvas[[1]]$type %||% lnk$coordGrid$type[1], "s")
})

test_that("seq_string end-to-end as plot-level link", {
  gr  <- make_string_data()
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win) %+%
    seq_string(data = gr,
               map(x0 = start, x1 = x1pos, chrom0 = seqnames, chrom1 = chrom1,
                   strand0 = s0, strand1 = s1, y0 = score, y1 = score),
               t0 = "A", t1 = "B")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_synteny ───────────────────────────────────────────────────────────────

test_that("seq_synteny instantiates", {
  df <- make_synteny_data()
  expect_no_error(
    seq_synteny(data = df,
                map(x0 = p0, x0_end = p0end, x1 = p1, x1_end = p1end,
                    chrom0 = c0, chrom1 = c1),
                t0 = "A", t1 = "B")
  )
})

test_that("seq_synteny inherits from SeqLink", {
  expect_true(inherits(seq_synteny(), "SeqLink"))
})

test_that("seq_synteny prep() runs without error", {
  tl  <- make_two_track_layout()
  df  <- make_synteny_data()
  lnk <- seq_synteny(data = df,
                     map(x0 = p0, x0_end = p0end, x1 = p1, x1_end = p1end,
                         chrom0 = c0, chrom1 = c1),
                     t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  expect_no_error(lnk$prep(layout_named, windows_named))
})

test_that("seq_synteny coordCanvas produces 4-point polygons", {
  tl  <- make_two_track_layout()
  df  <- make_synteny_data()
  lnk <- seq_synteny(data = df,
                     map(x0 = p0, x0_end = p0end, x1 = p1, x1_end = p1end,
                         chrom0 = c0, chrom1 = c1),
                     t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named)
  visible_polys <- Filter(Negate(is.null), lnk$coordCanvas)
  if (length(visible_polys) > 0)
    expect_equal(length(visible_polys[[1]]$x), 4)
})

test_that("seq_synteny draw() runs without error", {
  tl  <- make_two_track_layout()
  df  <- make_synteny_data()
  lnk <- seq_synteny(data = df,
                     map(x0 = p0, x0_end = p0end, x1 = p1, x1_end = p1end,
                         chrom0 = c0, chrom1 = c1),
                     t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named)
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(lnk$draw())
  dev.off()
})

# ── seq_zoom ──────────────────────────────────────────────────────────────────

make_zoom_data <- function() {
  GRanges("chr1", IRanges(c(200, 600), width = 100))
}

test_that("seq_zoom instantiates", {
  gr <- make_zoom_data()
  expect_no_error(
    seq_zoom(data = gr, map(x0 = start, x0_end = end), t0 = "A", t1 = "B")
  )
})

test_that("seq_zoom inherits from SeqLink", {
  expect_true(inherits(seq_zoom(), "SeqLink"))
})

test_that("seq_zoom prep() runs without error", {
  tl  <- make_two_track_layout()
  gr  <- make_zoom_data()
  lnk <- seq_zoom(data = gr, map(x0 = start, x0_end = end), t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  expect_no_error(lnk$prep(layout_named, windows_named))
})

test_that("seq_zoom produces a 4-column coordCanvas data.frame", {
  tl  <- make_two_track_layout()
  gr  <- make_zoom_data()
  lnk <- seq_zoom(data = gr, map(x0 = start, x0_end = end), t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named)
  if (!is.null(lnk$coordCanvas) && nrow(lnk$coordCanvas) > 0) {
    expect_true(all(c("x00", "x10", "x01", "x11", "y0", "y1") %in% names(lnk$coordCanvas)))
  }
})

test_that("seq_zoom draw() runs without error", {
  tl  <- make_two_track_layout()
  gr  <- make_zoom_data()
  lnk <- seq_zoom(data = gr, map(x0 = start, x0_end = end), t0 = "A", t1 = "B")
  layout_named  <- setNames(tl$layout$panelBounds,
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  windows_named <- setNames(lapply(tl$plot$allTracks(), function(t) t$windows),
                             vapply(tl$plot$allTracks(), function(t) t$track_id, character(1)))
  lnk$prep(layout_named, windows_named)
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(lnk$draw())
  dev.off()
})

test_that("seq_zoom end-to-end as plot-level link", {
  gr  <- make_zoom_data()
  win <- make_win()
  p <- seq_plot() %|%
    seq_track(track_id = "A", windows = win) %__%
    seq_track(track_id = "B", windows = win) %+%
    seq_zoom(data = gr, map(x0 = start, x0_end = end), t0 = "A", t1 = "B")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})
