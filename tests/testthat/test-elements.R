library(SeqPlotR)
library(GenomicRanges)

# ── Shared test data ──────────────────────────────────────────────────────────

make_gr <- function() {
  GRanges("chr1", IRanges(c(100, 300, 500, 700), width = 50),
          score = c(0.2, 0.8, 0.5, 0.3),
          group = c("A", "B", "A", "B"),
          label = c("a", "b", "c", "d"))
}

make_panel_meta <- function(xmin = 1, xmax = 1000, ymin = 0, ymax = 1) {
  list(list(
    xscale = c(xmin, xmax),
    yscale = c(ymin, ymax),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = make_gr(),
    track_mapping = map(x = start, y = score)
  ))
}

make_pm  <- make_panel_meta
make_win <- function() GRanges("chr1", IRanges(1, 1000))

# ── seq_point ─────────────────────────────────────────────────────────────────

test_that("seq_point instantiates", {
  expect_no_error(seq_point())
  expect_no_error(seq_point(map(x = start, y = score)))
})

test_that("seq_point resolve() populates resolved$x and resolved$y", {
  gr <- make_gr()
  pt <- seq_point()
  pt$resolve(track_data = gr, track_mapping = map(x = start, y = score))
  expect_equal(pt$resolved$x, BiocGenerics::start(gr))
  expect_equal(pt$resolved$y, gr$score)
})

test_that("seq_point prep() populates coordCanvas with npc values", {
  win <- GRanges("chr1", IRanges(1, 1000))
  pt  <- seq_point()
  pt$prep(make_panel_meta(), win)
  coords <- pt$coordCanvas[[1]]
  expect_false(is.null(coords))
  expect_true(all(coords$x >= 0 & coords$x <= 1))
  expect_true(all(coords$y >= 0 & coords$y <= 1))
})

test_that("seq_point draw() runs without error", {
  win <- GRanges("chr1", IRanges(1, 1000))
  pt  <- seq_point()
  pt$prep(make_panel_meta(), win)
  pdf(tempfile())
  grid::grid.newpage()
  expect_no_error(pt$draw())
  dev.off()
})

# ── seq_line ──────────────────────────────────────────────────────────────────

test_that("seq_line instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  ln  <- seq_line()
  expect_no_error(ln$prep(make_panel_meta(), win))
  expect_false(is.null(ln$coordCanvas))
})

test_that("seq_line draw() runs without error", {
  win <- GRanges("chr1", IRanges(1, 1000))
  ln  <- seq_line()
  ln$prep(make_panel_meta(), win)
  pdf(tempfile())
  grid::grid.newpage()
  expect_no_error(ln$draw())
  dev.off()
})

# ── seq_segment ───────────────────────────────────────────────────────────────

test_that("seq_segment instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  seg <- seq_segment(map(x = start, x_end = end, y = score, y_end = score))
  pm  <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data = gr,
    track_mapping = map(x = start, x_end = end, y = score, y_end = score)
  ))
  expect_no_error(seg$prep(pm, win))
})

test_that("seq_segment draw() runs without error", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  seg <- seq_segment(map(x = start, x_end = end, y = score, y_end = score))
  pm  <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data = gr,
    track_mapping = map(x = start, x_end = end, y = score, y_end = score)
  ))
  seg$prep(pm, win)
  pdf(tempfile())
  grid::grid.newpage()
  expect_no_error(seg$draw())
  dev.off()
})

# ── seq_curve ─────────────────────────────────────────────────────────────────

test_that("seq_curve instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  cv  <- seq_curve(map(x = start, y = score, x_end = end, y_end = score))
  pm  <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data = gr,
    track_mapping = map(x = start, y = score, x_end = end, y_end = score)
  ))
  expect_no_error(cv$prep(pm, win))
})

# ── seq_area ──────────────────────────────────────────────────────────────────

test_that("seq_area instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  ar  <- seq_area()
  expect_no_error(ar$prep(make_panel_meta(), win))
})

test_that("seq_area coordCanvas forms a closed polygon", {
  win <- GRanges("chr1", IRanges(1, 1000))
  ar  <- seq_area()
  ar$prep(make_panel_meta(), win)
  coords <- ar$coordCanvas[[1]]
  expect_false(is.null(coords))
  # Polygon should have 2x as many points as data (data line + baseline return)
  n_data <- length(make_gr())  # 4 points in window
  expect_equal(length(coords$x), n_data * 2)
})

test_that("seq_area draw() runs without error", {
  win <- GRanges("chr1", IRanges(1, 1000))
  ar  <- seq_area()
  ar$prep(make_panel_meta(), win)
  pdf(tempfile())
  grid::grid.newpage()
  expect_no_error(ar$draw())
  dev.off()
})

# ── seq_poly ──────────────────────────────────────────────────────────────────

test_that("seq_poly instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  pl  <- seq_poly()
  expect_no_error(pl$prep(make_panel_meta(), win))
})

# ── seq_path ──────────────────────────────────────────────────────────────────

test_that("seq_path instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  pa  <- seq_path()
  expect_no_error(pa$prep(make_panel_meta(), win))
})

# ── seq_text ──────────────────────────────────────────────────────────────────

test_that("seq_text instantiates and preps", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  tx  <- seq_text(map(x = start, y = score, label = label))
  pm  <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data = gr,
    track_mapping = map(x = start, y = score, label = label)
  ))
  expect_no_error(tx$prep(pm, win))
})

test_that("seq_text draw() runs without error", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  tx  <- seq_text(map(x = start, y = score, label = label))
  pm  <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data = gr,
    track_mapping = map(x = start, y = score, label = label)
  ))
  tx$prep(pm, win)
  pdf(tempfile())
  grid::grid.newpage()
  expect_no_error(tx$draw())
  dev.off()
})

# ── End-to-end ────────────────────────────────────────────────────────────────

test_that("seq_point renders end-to-end via seq_plot", {
  win <- GRanges("chr1", IRanges(1, 1000))
  gr  <- make_gr()
  p <- seq_plot() %|%
    seq_track(data = gr, mapping = map(x = start, y = score), windows = win) %+%
    seq_point()
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_bar ───────────────────────────────────────────────────────────────────

test_that("seq_bar instantiates", {
  expect_no_error(seq_bar())
  expect_no_error(seq_bar(map(x = start, y = score)))
})

test_that("seq_bar prep() populates coordCanvas", {
  bar <- seq_bar()
  expect_no_error(bar$prep(make_pm(), make_win()))
  expect_false(is.null(bar$coordCanvas[[1]]))
})

test_that("seq_bar produces correct bar count", {
  bar <- seq_bar()
  bar$prep(make_pm(), make_win())
  expect_equal(nrow(bar$coordCanvas[[1]]), length(make_gr()))
})

test_that("seq_bar bars are geometrically valid (x0<=x1, y0<=y1)", {
  bar <- seq_bar()
  bar$prep(make_pm(), make_win())
  coords <- bar$coordCanvas[[1]]
  expect_true(all(coords$x0 <= coords$x1))
  expect_true(all(coords$y0 <= coords$y1))
})

test_that("seq_bar draw() runs without error", {
  bar <- seq_bar()
  bar$prep(make_pm(), make_win())
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(bar$draw())
  dev.off()
})

test_that("seq_bar end-to-end via seq_plot", {
  p <- seq_plot() %|%
    seq_track(data = make_gr(),
              mapping = map(x = start, y = score),
              windows = make_win()) %+%
    seq_bar()
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_ribbon ────────────────────────────────────────────────────────────────

test_that("seq_ribbon instantiates", {
  expect_no_error(seq_ribbon(map(x = start, y_min = score, y_max = score)))
})

test_that("seq_ribbon polygon has 2x data points", {
  gr <- make_gr()
  pm <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(x = start, y_min = score, y_max = score)
  ))
  rb <- seq_ribbon(map(x = start, y_min = score, y_max = score))
  rb$prep(pm, make_win())
  expect_equal(length(rb$coordCanvas[[1]]$x), length(gr) * 2)
})

test_that("seq_ribbon draw() runs without error", {
  gr <- make_gr()
  pm <- list(list(
    xscale = c(1, 1000), yscale = c(0, 1),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(x = start, y_min = score, y_max = score)
  ))
  rb <- seq_ribbon(map(x = start, y_min = score, y_max = score))
  rb$prep(pm, make_win()); pdf(tempfile()); grid::grid.newpage()
  expect_no_error(rb$draw())
  dev.off()
})

# ── seq_density ───────────────────────────────────────────────────────────────

test_that("seq_density instantiates and preps", {
  expect_no_error(seq_density())
  dn <- seq_density(); dn$prep(make_pm(), make_win())
  expect_false(is.null(dn$coordCanvas[[1]]))
})

test_that("seq_density polygon size matches stats::density output", {
  dn <- seq_density()
  gr <- make_gr()
  dn$resolve(track_data = gr, track_mapping = map(x = start, y = score))
  ref <- stats::density(gr$score, bw = "nrd0")
  dn$prep(make_pm(), make_win())
  coords <- dn$coordCanvas[[1]]
  expect_equal(length(coords$x), length(ref$x) * 2)
})

test_that("seq_density draw() runs without error", {
  dn <- seq_density(); dn$prep(make_pm(), make_win())
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(dn$draw())
  dev.off()
})

# ── seq_tile ──────────────────────────────────────────────────────────────────

make_tile_gr <- function() {
  GRanges("chr1", IRanges(c(100, 300, 500), width = 100),
          color = c("red", "blue", "green"))
}

test_that("seq_tile unrotated instantiates and preps", {
  gr <- make_tile_gr()
  pm <- list(list(
    xscale = c(1, 1000), yscale = c(0.5, 1.5),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(x = start, fill = color)
  ))
  tl <- seq_tile(aesthetics = aes(rotate = FALSE))
  expect_no_error(tl$prep(pm, make_win()))
  expect_false(is.null(tl$coordCanvas[[1]]))
})

test_that("seq_tile rotated preps without error", {
  gr_x <- make_tile_gr()
  gr_y <- GRanges("chr1", IRanges(c(100, 300, 500), width = 100))
  tl   <- seq_tile(data2 = gr_y, aesthetics = aes(rotate = TRUE))
  pm   <- list(list(
    xscale = c(1, 1000), yscale = c(0, 500),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr_x,
    track_mapping = map(x = start, fill = color)
  ))
  expect_no_error(tl$prep(pm, make_win()))
})

test_that("seq_tile rotated produces different x0 from unrotated", {
  gr_x <- make_tile_gr()
  gr_y <- GRanges("chr1", IRanges(c(200, 400, 600), width = 100))
  pm_flat <- list(list(
    xscale = c(1, 1000), yscale = c(0.5, 1.5),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr_x,
    track_mapping = map(x = start, fill = color)
  ))
  pm_rot <- list(list(
    xscale = c(1, 1000), yscale = c(0, 500),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr_x,
    track_mapping = map(x = start, fill = color)
  ))
  tl_flat <- seq_tile(aesthetics = aes(rotate = FALSE))
  tl_rot  <- seq_tile(data2 = gr_y, aesthetics = aes(rotate = TRUE))
  tl_flat$prep(pm_flat, make_win())
  tl_rot$prep(pm_rot,  make_win())
  x_flat <- tl_flat$coordCanvas[[1]]$x0[1]
  x_rot  <- tl_rot$coordCanvas[[1]]$x0[1]
  expect_false(isTRUE(all.equal(x_flat, x_rot)))
})

test_that("seq_tile draw() unrotated runs without error", {
  gr <- make_tile_gr()
  pm <- list(list(
    xscale = c(1, 1000), yscale = c(0.5, 1.5),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(x = start, fill = color)
  ))
  tl <- seq_tile(aesthetics = aes(rotate = FALSE))
  tl$prep(pm, make_win()); pdf(tempfile()); grid::grid.newpage()
  expect_no_error(tl$draw())
  dev.off()
})

# ── seq_lollipop ──────────────────────────────────────────────────────────────

test_that("seq_lollipop instantiates and preps", {
  lp <- seq_lollipop()
  expect_no_error(lp$prep(make_pm(), make_win()))
  expect_true(all(c("x", "y0", "y1") %in% names(lp$coordCanvas[[1]])))
})

test_that("seq_lollipop y1 >= y0 for all stems", {
  lp <- seq_lollipop()
  lp$prep(make_pm(), make_win())
  coords <- lp$coordCanvas[[1]]
  expect_true(all(coords$y1 >= coords$y0))
})

test_that("seq_lollipop draw() runs without error", {
  lp <- seq_lollipop(); lp$prep(make_pm(), make_win())
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(lp$draw())
  dev.off()
})

# ── seq_gene ──────────────────────────────────────────────────────────────────

make_gene_gr <- function() {
  GRanges("chr1", IRanges(c(100, 200, 500, 600), width = 80),
          gene = c("G1", "G1", "G2", "G2"),
          feat = c("exon", "exon", "exon", "UTR"),
          str  = c("+", "+", "-", "-"),
          nm   = c("GeneA", "GeneA", "GeneB", "GeneB"),
          col  = c("blue", "blue", "red", "red"))
}
make_gene_pm <- function() {
  gr <- make_gene_gr()
  list(list(
    xscale = c(1, 1000), yscale = c(0, 4),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(group = gene, type = feat,
                        strand = str, label = nm, color = col)
  ))
}

test_that("seq_gene instantiates with arbitrary column names", {
  expect_no_error(
    seq_gene(map(group = gene, type = feat, strand = str, label = nm))
  )
})

test_that("seq_gene preps without error", {
  gn <- seq_gene(map(group = gene, type = feat, strand = str, label = nm))
  expect_no_error(gn$prep(make_gene_pm(), make_win()))
})

test_that("seq_gene coordCanvas is a non-empty data.frame", {
  gn <- seq_gene(map(group = gene, type = feat, strand = str, label = nm))
  gn$prep(make_gene_pm(), make_win())
  expect_true(is.data.frame(gn$coordCanvas))
  expect_gt(nrow(gn$coordCanvas), 0)
})

test_that("seq_gene draw() runs without error", {
  gn <- seq_gene(map(group = gene, type = feat, strand = str, label = nm))
  gn$prep(make_gene_pm(), make_win())
  pdf(tempfile()); grid::grid.newpage()
  expect_no_error(gn$draw())
  dev.off()
})

test_that("seq_gene works with no strand mapping", {
  gr <- GRanges("chr1", IRanges(c(100, 200), width = 80),
                gene = c("G1", "G1"), nm = c("GeneA", "GeneA"))
  pm <- list(list(
    xscale = c(1, 1000), yscale = c(0, 2),
    inner  = list(x0 = 0.1, x1 = 0.9, y0 = 0.1, y1 = 0.9),
    track_data    = gr,
    track_mapping = map(group = gene, label = nm)
  ))
  gn <- seq_gene(map(group = gene, label = nm))
  expect_no_error(gn$prep(pm, make_win()))
})

test_that("seq_gene end-to-end via seq_plot", {
  gr <- make_gene_gr()
  p <- seq_plot() %|%
    seq_track(data = gr,
              mapping = map(group = gene, type = feat,
                            strand = str, label = nm),
              windows = make_win()) %+%
    seq_gene(map(group = gene, type = feat, strand = str, label = nm))
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})
