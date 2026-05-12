library(SeqPlotR)
library(GenomicRanges)

# ── Shared synthetic data ──────────────────────────────────────────────────────

make_windows <- function() {
  GRanges(c("chr1", "chr2", "chr3"),
          IRanges(c(1, 1, 1), c(1e6, 2e6, 1.5e6)))
}

make_cn_gr <- function(n = 200) {
  set.seed(42)
  chrs <- sample(c("chr1", "chr2", "chr3"), n, replace = TRUE)
  GRanges(chrs,
          IRanges(start = sample(1:2e6, n), width = 5000),
          cn        = sample(0:6, n, replace = TRUE,
                             prob = c(.05, .1, .5, .15, .1, .05, .05)),
          log2ratio = rnorm(n, mean = 0, sd = 0.3))
}

# Hi-C simulator ported from THEfunc/inst/examples/08_rotated_hic_heatmap.R:
# upper-triangular bin pairs with exponential distance decay × lognormal
# noise, then mirrored to produce the symmetric (j, i) contacts.
.generate_hic_matrix <- function(n_bins, decay_rate = 0.2) {
  ij    <- expand.grid(i = seq_len(n_bins), j = seq_len(n_bins))
  upper <- ij$j >= ij$i
  ij    <- ij[upper, , drop = FALSE]
  d     <- ij$j - ij$i
  strength <- pmax(0.01,
                   exp(-d * decay_rate) *
                   stats::rlnorm(nrow(ij), meanlog = 0, sdlog = 0.3))
  upper_df <- data.frame(bin_i = ij$i, bin_j = ij$j, strength = strength)
  off_diag <- d > 0
  lower_df <- data.frame(bin_i = ij$j[off_diag],
                         bin_j = ij$i[off_diag],
                         strength = strength[off_diag])
  rbind(upper_df, lower_df)
}

make_hic_gr <- function(n_bins = 30, bin_size = 1e4, decay_rate = 0.2,
                        seqname = "chr1", region_start = 1L) {
  set.seed(123)
  mat <- .generate_hic_matrix(n_bins, decay_rate = decay_rate)
  bin_starts <- region_start + (seq_len(n_bins) - 1L) * bin_size
  GRanges(seqname,
          IRanges(start = bin_starts[mat$bin_i], width = bin_size),
          i_start = bin_starts[mat$bin_i],
          i_end   = bin_starts[mat$bin_i] + bin_size,
          j_start = bin_starts[mat$bin_j],
          j_end   = bin_starts[mat$bin_j] + bin_size,
          score   = mat$strength)
}

make_chip_gr <- function(n = 300, sample_name = "S1") {
  set.seed(42)
  GRanges("chr1",
          IRanges(sort(sample(1:1e6, n)), width = 500),
          score  = rexp(n, rate = 0.2),
          sample = sample_name)
}

make_peaks_gr <- function(n = 30) {
  GRanges("chr1",
          IRanges(sort(sample(1:1e6, n)), width = 2000))
}

# ── seq_copynumber ────────────────────────────────────────────────────────────

test_that("seq_copynumber returns a seq_plot", {
  p <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                      cn_col = "cn", ratio_col = "log2ratio")
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_copynumber plot() renders without error", {
  p <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                      cn_col = "cn", ratio_col = "log2ratio")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_copynumber auto-detects columns", {
  gr <- make_cn_gr()
  names(mcols(gr))[names(mcols(gr)) == "cn"]        <- "copy_number"
  names(mcols(gr))[names(mcols(gr)) == "log2ratio"] <- "logR"
  p <- seq_copynumber(data = gr, windows = make_windows())
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_copynumber with segment_data renders", {
  cn  <- make_cn_gr()
  seg <- GRanges(c("chr1", "chr2"),
                 IRanges(c(1, 1), c(1e6, 2e6)),
                 seg_mean = c(0.1, -0.4))
  p <- seq_copynumber(data = cn, windows = make_windows(),
                      cn_col = "cn", ratio_col = "log2ratio",
                      segment_data = seg, segment_col = "seg_mean")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_cn_heatmap ────────────────────────────────────────────────────────────

test_that("seq_cn_heatmap returns a seq_plot", {
  gr <- make_cn_gr(500)
  mcols(gr)$sample <- sample(paste0("S", 1:5), length(gr), replace = TRUE)
  p  <- seq_cn_heatmap(data = gr, windows = make_windows(),
                       sample_col = "sample", cn_col = "cn")
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_cn_heatmap renders without error", {
  gr <- make_cn_gr(500)
  mcols(gr)$sample <- sample(paste0("S", 1:5), length(gr), replace = TRUE)
  p  <- seq_cn_heatmap(data = gr, windows = make_windows(),
                       sample_col = "sample", cn_col = "cn")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_hic ───────────────────────────────────────────────────────────────────

test_that("seq_hic default style (triangle) returns seq_plot", {
  p <- seq_hic(data = make_hic_gr(),
               windows = GRanges("chr1", IRanges(1, 1e6)))
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_hic full style returns seq_plot", {
  p <- seq_hic(data = make_hic_gr(),
               windows = GRanges("chr1", IRanges(1, 1e6)),
               style   = "full")
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_hic triangle style renders without error", {
  p <- seq_hic(data = make_hic_gr(),
               windows = GRanges("chr1", IRanges(1, 1e6)),
               style   = "triangle")
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_hic rectangle style with max_dist renders without error", {
  p <- seq_hic(data = make_hic_gr(),
               windows  = GRanges("chr1", IRanges(1, 1e6)),
               style    = "rectangle",
               max_dist = 3e5)
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_hic rectangle errors without max_dist", {
  expect_error(
    seq_hic(data   = make_hic_gr(),
            windows = GRanges("chr1", IRanges(1, 1e6)),
            style  = "rectangle"),
    "max_dist is required"
  )
})

test_that("seq_hic errors on missing required columns", {
  bad_gr <- GRanges("chr1", IRanges(1:5 * 1e5, width = 1e4), score = 1:5)
  expect_error(
    seq_hic(data = bad_gr, windows = GRanges("chr1", IRanges(1, 1e6))),
    "missing required columns"
  )
})

test_that("seq_hic two styles combined via seq_resolve renders", {
  win  <- GRanges("chr1", IRanges(1, 1e6))
  tri  <- seq_hic(data = make_hic_gr(), windows = win, style = "triangle",
                  track_id = "tri")
  full <- seq_hic(data = make_hic_gr(), windows = win, style = "full",
                  track_id = "full")
  fig  <- seq_resolve(seq_plot(), tri, full)
  pdf(tempfile())
  expect_no_error(fig$plot())
  dev.off()
})

# ── seq_chip ──────────────────────────────────────────────────────────────────

test_that("seq_chip with named list returns seq_plot", {
  sig <- list(Rad21   = make_chip_gr(300, "Rad21"),
              NippedB = make_chip_gr(300, "NippedB"))
  p   <- seq_chip(data = sig, windows = GRanges("chr1", IRanges(1, 1e6)))
  expect_true(inherits(p, "SeqPlot"))
})

test_that("seq_chip renders without error", {
  sig <- list(S1 = make_chip_gr(), S2 = make_chip_gr())
  p   <- seq_chip(data = sig, windows = GRanges("chr1", IRanges(1, 1e6)))
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_chip with peaks renders without error", {
  sig   <- list(S1 = make_chip_gr())
  peaks <- list(S1 = make_peaks_gr())
  p     <- seq_chip(data = sig, peaks = peaks,
                    windows = GRanges("chr1", IRanges(1, 1e6)))
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

test_that("seq_chip with single GRanges and sample_col renders", {
  gr <- c(make_chip_gr(150, "S1"), make_chip_gr(150, "S2"))
  mcols(gr)$sample_id <- c(rep("S1", 150), rep("S2", 150))
  p  <- seq_chip(data = gr, sample_col = "sample_id",
                 windows = GRanges("chr1", IRanges(1, 1e6)))
  pdf(tempfile())
  expect_no_error(p$plot())
  dev.off()
})

# ── seq_resolve ───────────────────────────────────────────────────────────────

test_that("seq_resolve adds child tracks to parent", {
  cn <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                       cn_col = "cn", ratio_col = "log2ratio")
  parent <- seq_plot()
  result <- seq_resolve(parent, cn)
  expect_gt(length(result$allTracks()), 0)
})

test_that("seq_resolve with multiple children stacks correctly", {
  cn1 <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                        cn_col = "cn", ratio_col = "log2ratio",
                        track_id = "CN1")
  cn2 <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                        cn_col = "cn", ratio_col = "log2ratio",
                        track_id = "CN2")
  result <- seq_resolve(seq_plot(), cn1, cn2)
  expect_gte(length(result$allTracks()), 2)
})

test_that("seq_resolve combined plot renders without error", {
  cn <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                       cn_col = "cn", ratio_col = "log2ratio",
                       track_id = "cn")
  hic <- seq_hic(data = make_hic_gr(),
                 windows = GRanges("chr1", IRanges(1, 1e6)),
                 style   = "triangle",
                 track_id = "hic")
  fig <- seq_resolve(seq_plot(), cn, hic)
  pdf(tempfile())
  expect_no_error(fig$plot())
  dev.off()
})

test_that("seq_resolve errors on non-seq_plot parent", {
  expect_error(
    seq_resolve("not_a_plot",
                seq_copynumber(make_cn_gr(), make_windows(),
                               cn_col = "cn", ratio_col = "log2ratio")),
    "parent must be a seq_plot"
  )
})

test_that("seq_resolve errors on duplicate track_ids", {
  cn1 <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                        cn_col = "cn", ratio_col = "log2ratio",
                        track_id = "dup")
  cn2 <- seq_copynumber(data = make_cn_gr(), windows = make_windows(),
                        cn_col = "cn", ratio_col = "log2ratio",
                        track_id = "dup")
  expect_error(seq_resolve(seq_plot(), cn1, cn2), "duplicate track_id")
})
