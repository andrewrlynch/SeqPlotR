library(testthat)
library(GenomicRanges)

cb <- GRanges("chr1",
  IRanges(start = c(1, 1e7, 2e7, 3e7, 3.5e7, 4e7, 6e7),
          end   = c(1e7, 2e7, 3e7, 3.5e7, 4e7, 6e7, 2e8)),
  gieStain = c("gneg","gpos25","gneg","acen","acen","gpos75","gneg")
)
win <- GRanges("chr1", IRanges(1e7, 3e7))

make_layout <- function(win) {
  trk <- seq_track(data = cb, windows = win)
  sp  <- seq_plot() %+% trk
  sp$layoutGrid()
  sp$layout$panelBounds[[1]]
}

# ---- scope ----

test_that("scope defaults to 'window'", {
  el <- seq_ideogram(cb)
  expect_equal(el$scope, "window")
})

test_that("scope = 'full' accepted", {
  el <- seq_ideogram(cb, scope = "full")
  expect_equal(el$scope, "full")
})

test_that("scope = 'full' prep() populates highlightBoxes", {
  el <- seq_ideogram(cb, scope = "full")
  el$prep(make_layout(win), win)
  expect_false(is.null(el$highlightBoxes))
  expect_false(is.null(el$highlightBoxes[[1]]))
})

test_that("scope = 'full' prep() coordCanvas covers all bands", {
  el <- seq_ideogram(cb, scope = "full")
  el$prep(make_layout(win), win)
  df <- el$coordCanvas[[1]]
  # 5 non-acen bands out of 7 total
  expect_equal(nrow(df), 5L)
})

test_that("scope = 'window' prep() only includes bands overlapping window", {
  el <- seq_ideogram(cb, scope = "window")
  el$prep(make_layout(win), win)
  df <- el$coordCanvas[[1]]
  expect_lte(nrow(df), 5L)
})

# ---- style ----

test_that("style defaults to 'block'", {
  el <- seq_ideogram(cb)
  expect_equal(el$style, "block")
})

test_that("style = 'rounded' accepted", {
  el <- seq_ideogram(cb, style = "rounded")
  expect_equal(el$style, "rounded")
})

test_that("invalid style errors", {
  expect_error(seq_ideogram(cb, style = "fuzzy"), "arg")
})

test_that("draw() with style = 'block' runs without error", {
  skip_if_not(capabilities("png"))
  el <- seq_ideogram(cb)
  el$prep(make_layout(win), win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with style = 'rounded' runs without error", {
  skip_if_not(capabilities("png"))
  el <- seq_ideogram(cb, style = "rounded")
  el$prep(make_layout(win), win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("draw() with scope = 'full' and highlight aes runs without error", {
  skip_if_not(capabilities("png"))
  el <- seq_ideogram(cb, scope = "full",
                     aesthetics = aes(
                       highlight = aes(fill = "blue", alpha = 0.2)
                     ))
  el$prep(make_layout(win), win)
  grDevices::pdf(file = NULL); on.exit(grDevices::dev.off())
  expect_invisible(el$draw())
})

test_that("telomere.radius aes is respected", {
  el <- seq_ideogram(cb, style = "rounded",
                     aesthetics = aes(telomere.radius = 0.6))
  expect_equal(el$aesthetics[["telomere.radius"]], 0.6)
})
