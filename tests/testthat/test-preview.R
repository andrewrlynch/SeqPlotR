library(SeqPlotR)
library(GenomicRanges)

# ── seq_preview_circos ────────────────────────────────────────────────────────

make_circos_plot <- function() {
  seq_plot() %|%
    seq_track(track_id = "Chr1", track_width = 3,
              windows = GRanges("chr1", IRanges(1, 1000))) %|%
    seq_track(track_id = "Chr2", track_width = 2,
              windows = GRanges("chr2", IRanges(1, 1000))) %|%
    seq_track(track_id = "Chr3", track_width = 1,
              windows = GRanges("chr3", IRanges(1, 1000))) %__%
    seq_track(track_id = "Signal",
              windows = GRanges("chr1", IRanges(1, 1000))) %__%
    seq_track(track_id = "CopyNum",
              windows = GRanges("chr1", IRanges(1, 1000)))
}

test_that("seq_preview_circos runs without error", {
  p <- make_circos_plot()
  pdf(tempfile())
  expect_no_error(seq_preview_circos(plot_obj = p))
  dev.off()
})

test_that("seq_preview_circos returns polar bounds invisibly", {
  p   <- make_circos_plot()
  pdf(tempfile())
  res <- seq_preview_circos(plot_obj = p)
  dev.off()
  expect_type(res, "list")
  expect_named(res, c("Chr1", "Chr2", "Chr3", "Signal", "CopyNum"))
})

test_that("seq_preview_circos polar bounds have correct fields", {
  p   <- make_circos_plot()
  pdf(tempfile())
  res <- seq_preview_circos(plot_obj = p)
  dev.off()
  for (id in names(res)) {
    expect_true(all(c("theta0", "theta1", "r0", "r1") %in% names(res[[id]])))
  }
})

test_that("seq_preview_circos theta spans sum to ~360 degrees", {
  p   <- make_circos_plot()
  pdf(tempfile())
  res <- seq_preview_circos(plot_obj = p, gap_degrees = 0)
  dev.off()
  outer_span <- sum(vapply(c("Chr1", "Chr2", "Chr3"),
                           function(id) abs(res[[id]]$theta0 - res[[id]]$theta1),
                           numeric(1)))
  expect_equal(outer_span, 360, tolerance = 0.5)
})

test_that("seq_preview_circos track_width ratio is respected", {
  p   <- make_circos_plot()
  pdf(tempfile())
  res <- seq_preview_circos(plot_obj = p, gap_degrees = 0)
  dev.off()
  span1 <- abs(res$Chr1$theta0 - res$Chr1$theta1)
  span2 <- abs(res$Chr2$theta0 - res$Chr2$theta1)
  span3 <- abs(res$Chr3$theta0 - res$Chr3$theta1)
  expect_equal(span1 / span2, 1.5, tolerance = 0.01)
  expect_equal(span2 / span3, 2.0, tolerance = 0.01)
})

test_that("seq_preview_circos rings are ordered outer-to-inner", {
  p   <- make_circos_plot()
  pdf(tempfile())
  res <- seq_preview_circos(plot_obj = p)
  dev.off()
  expect_gt(res$Chr1$r1,   res$Signal$r1)
  expect_gt(res$Signal$r1, res$CopyNum$r1)
})

test_that(".arc_polygon produces 2*n_pts vertices", {
  poly <- SeqPlotR:::.arc_polygon(90, -270, r0 = 0.1, r1 = 0.4, n_pts = 30)
  expect_length(poly$x, 60)
  expect_length(poly$y, 60)
})

test_that(".polar_to_npc top (90 degrees) maps to top of canvas", {
  pt <- SeqPlotR:::.polar_to_npc(r = 0.4, theta_deg = 90, cx = 0.5, cy = 0.5)
  expect_equal(pt$x, 0.5, tolerance = 1e-10)
  expect_gt(pt$y, 0.5)
})
