library(testthat)

test_that(".expand_limits applies multiplicative + additive padding", {
  expect_equal(.expand_limits(c(0, 100), c(0.05, 0)), c(-5, 105))
  expect_equal(.expand_limits(c(0, 100), c(0, 2)),    c(-2, 102))
  expect_equal(.expand_limits(c(0, 100), c(0.1, 1)),  c(-11, 111))
  expect_equal(.expand_limits(c(0, 100), 0.05),       c(-5, 105))
  expect_equal(.expand_limits(c(0, 100), c(0, 0)),    c(0, 100))
})

test_that(".expand_limits is safe for degenerate input", {
  expect_null(.expand_limits(NULL, c(0.1, 0)))
  r <- c(1, 1)
  expect_equal(.expand_limits(r, c(0.1, 0)), c(1, 1))
})

test_that(".compute_scale_breaks continuous with cap='capped'", {
  sc <- seq_scale_continuous(limits = c(0, 100), n_breaks = 5,
                             expand = c(0, 0), cap = "capped")
  meta <- .compute_scale_breaks(sc, c(0, 100))
  expect_gte(length(meta$breaks), 3L)
  expect_true(0   %in% meta$breaks)
  expect_true(100 %in% meta$breaks)
  expect_equal(meta$axis_range, range(meta$breaks))
})

test_that(".compute_scale_breaks cap='full' spans expanded range", {
  sc <- seq_scale_continuous(limits = c(0, 100), n_breaks = 5,
                             expand = c(0.05, 0), cap = "full")
  meta <- .compute_scale_breaks(sc, c(0, 100))
  expect_equal(meta$axis_range, meta$plot_range)
  expect_equal(meta$plot_range, c(-5, 105))
})

test_that(".compute_scale_breaks cap='exact' uses unexpanded data range", {
  sc <- seq_scale_continuous(limits = c(0, 100), n_breaks = 5,
                             expand = c(0.05, 0), cap = "exact")
  meta <- .compute_scale_breaks(sc, c(0, 100))
  expect_equal(meta$axis_range, c(0, 100))
})

test_that(".compute_scale_breaks cap='ticks' suppresses axis line", {
  sc <- seq_scale_continuous(limits = c(0, 100), n_breaks = 5,
                             expand = c(0.05, 0), cap = "ticks")
  meta <- .compute_scale_breaks(sc, c(0, 100))
  expect_null(meta$axis_range)
})

test_that(".compute_scale_breaks honours explicit breaks", {
  sc <- seq_scale_continuous(limits = c(0, 100),
                             breaks = c(10, 50, 90),
                             expand = c(0, 0), cap = "capped")
  meta <- .compute_scale_breaks(sc, c(0, 100))
  expect_equal(meta$breaks, c(10, 50, 90))
  expect_equal(meta$axis_range, c(10, 90))
})

test_that(".compute_minor_breaks subdivides between majors", {
  mb <- .compute_minor_breaks(3, c(0, 10, 20), c(0, 20))
  # 3 interior points between each pair of majors ⇒ 6 total.
  expect_length(mb, 6L)
  expect_true(all(mb > 0 & mb < 20))
})

test_that(".compute_minor_breaks filters vector input to range", {
  mb <- .compute_minor_breaks(c(-5, 5, 15, 25), c(0, 20), c(0, 20))
  expect_equal(mb, c(5, 15))
})

test_that(".merge_scale_with_theme builds a scale from theme when none supplied", {
  flat <- .default_theme()
  flat[["axis.x1.scale.limits"]] <- c(1, 99)
  spec <- .build_axis_spec(flat, "x1")
  sc   <- .merge_scale_with_theme(NULL, spec)
  expect_false(is.null(sc))
  expect_equal(sc$limits, c(1, 99))
})

test_that(".merge_scale_with_theme preserves explicit scale fields", {
  flat <- .default_theme()
  flat[["axis.x1.scale.limits"]] <- c(1, 99)
  flat[["axis.x1.scale.breaks"]] <- c(25, 50, 75)
  spec <- .build_axis_spec(flat, "x1")
  user_sc <- seq_scale_continuous(limits = c(0, 100))
  merged  <- .merge_scale_with_theme(user_sc, spec)
  # Explicit user limits stand.
  expect_equal(merged$limits, c(0, 100))
  # NULL breaks filled from theme.
  expect_equal(merged$breaks, c(25, 50, 75))
})
