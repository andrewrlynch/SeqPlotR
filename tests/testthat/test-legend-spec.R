library(testthat)

k1 <- LegendKey(label = "A", color = "red")
k2 <- LegendKey(label = "B", color = "blue")

# --- Construction ---

test_that("seq_legend accepts a single LegendKey", {
  spec <- seq_legend(k1)
  expect_s3_class(spec, "SeqLegendSpec")
  expect_length(spec$keys, 1)
})

test_that("seq_legend accepts a list of LegendKeys", {
  spec <- seq_legend(list(k1, k2))
  expect_length(spec$keys, 2)
})

test_that("seq_legend stores title", {
  spec <- seq_legend(k1, title = "Marks")
  expect_equal(spec$title, "Marks")
})

test_that("seq_legend title defaults to NULL", {
  spec <- seq_legend(k1)
  expect_null(spec$title)
})

test_that("seq_legend position defaults to 'inside'", {
  spec <- seq_legend(k1)
  expect_equal(spec$position, "inside")
})

test_that("seq_legend stores x, y, hjust", {
  spec <- seq_legend(k1, x = 0.1, y = 0.9, hjust = 1)
  expect_equal(spec$x, 0.1)
  expect_equal(spec$y, 0.9)
  expect_equal(spec$hjust, 1)
})

# --- Orientation inference ---

test_that("orientation defaults to 'horizontal' for top/bottom side", {
  spec <- seq_legend(k1, position = "track_margin", side = "top")
  expect_equal(spec$orientation, "horizontal")
})

test_that("orientation defaults to 'horizontal' for 'inside' (no side)", {
  spec <- seq_legend(k1, position = "inside")
  expect_equal(spec$orientation, "horizontal")
})

test_that("orientation defaults to 'vertical' for left side", {
  spec <- seq_legend(k1, position = "track_margin", side = "left")
  expect_equal(spec$orientation, "vertical")
})

test_that("orientation defaults to 'vertical' for right side", {
  spec <- seq_legend(k1, position = "canvas_margin", side = "right")
  expect_equal(spec$orientation, "vertical")
})

test_that("explicit orientation overrides inference", {
  spec <- seq_legend(k1, position = "track_margin", side = "left",
                     orientation = "horizontal")
  expect_equal(spec$orientation, "horizontal")
})

# --- Side default ---

test_that("side defaults to NULL for position='inside'", {
  spec <- seq_legend(k1, position = "inside")
  expect_null(spec$side)
})

test_that("side defaults to 'top' for position='track_margin'", {
  spec <- seq_legend(k1, position = "track_margin")
  expect_equal(spec$side, "top")
})

test_that("side defaults to 'top' for position='canvas_margin'", {
  spec <- seq_legend(k1, position = "canvas_margin")
  expect_equal(spec$side, "top")
})

# --- nrow / ncol stored as-is ---

test_that("nrow and ncol are stored", {
  spec <- seq_legend(list(k1, k2), nrow = 2, ncol = 1)
  expect_equal(spec$nrow, 2)
  expect_equal(spec$ncol, 1)
})

# --- Validation ---

test_that("seq_legend errors on empty keys list", {
  expect_error(seq_legend(list()), "non-empty")
})

test_that("seq_legend errors on non-LegendKey in list", {
  expect_error(seq_legend(list(k1, "bad")), "LegendKey")
})

test_that("seq_legend errors on invalid position", {
  expect_error(seq_legend(k1, position = "outer"), "position")
})

test_that("seq_legend errors on x out of range", {
  expect_error(seq_legend(k1, x = 1.5), "\\[0, 1\\]")
})

test_that("seq_legend errors on y out of range", {
  expect_error(seq_legend(k1, y = -0.1), "\\[0, 1\\]")
})

test_that("seq_legend errors on hjust out of range", {
  expect_error(seq_legend(k1, hjust = 2), "\\[0, 1\\]")
})

test_that("seq_legend errors on invalid orientation", {
  expect_error(seq_legend(k1, orientation = "diagonal"), "orientation")
})

test_that("seq_legend errors on invalid side", {
  expect_error(seq_legend(k1, side = "center"), "side")
})

# --- is_seq_legend_spec ---

test_that("is_seq_legend_spec returns TRUE for SeqLegendSpec", {
  expect_true(is_seq_legend_spec(seq_legend(k1)))
})

test_that("is_seq_legend_spec returns FALSE for other objects", {
  expect_false(is_seq_legend_spec(k1))
  expect_false(is_seq_legend_spec(list()))
})

# --- print method ---

test_that("print.SeqLegendSpec produces output without error", {
  spec <- seq_legend(list(k1, k2), title = "Marks", position = "track_margin",
                     side = "bottom")
  expect_output(print(spec), "SeqLegendSpec")
  expect_output(print(spec), "track_margin")
  expect_output(print(spec), "Marks")
})
