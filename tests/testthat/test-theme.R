library(testthat)

test_that(".flatten_theme returns identical output for flat and nested forms", {
  flat <- aes(axis.x.position = "top")
  nested <- aes(axis.x = aes(position = "top"))
  expect_identical(.flatten_theme(flat), .flatten_theme(nested))
  expect_identical(.flatten_theme(flat)[["axis.x.position"]], "top")
})

test_that(".flatten_theme handles deep nesting", {
  a <- aes(axis = aes(x = aes(line = aes(col = "red"))))
  out <- .flatten_theme(a)
  expect_identical(out[["axis.x.line.col"]], "red")
})

test_that(".flatten_theme preserves vector leaves", {
  a <- aes(axis.x1.scale.limits = c(1, 1000),
           axis.x1.scale.expand = c(0.05, 0))
  out <- .flatten_theme(a)
  expect_identical(out[["axis.x1.scale.limits"]], c(1, 1000))
  expect_identical(out[["axis.x1.scale.expand"]], c(0.05, 0))
})

test_that(".flatten_theme last-write wins on collision", {
  a <- aes(axis.x.line.col = "red",
           axis.x = aes(line = aes(col = "blue")))
  out <- .flatten_theme(a)
  expect_identical(out[["axis.x.line.col"]], "blue")
})

test_that(".axis_inheritance_chain walks specific to general", {
  chain <- .axis_inheritance_chain("axis.x1.line.col")
  expect_identical(chain,
                   c("axis.x1.line.col", "axis.x.line.col", "axis.line.col"))
})

test_that(".axis_inheritance_chain handles the general-form side", {
  chain <- .axis_inheritance_chain("axis.x.line.col")
  expect_identical(chain, c("axis.x.line.col", "axis.line.col"))
})

test_that(".resolve_theme walks the inheritance chain most-specific first", {
  theme <- list("axis.line.col" = "A",
                "axis.x.line.col" = "B",
                "axis.x1.line.col" = "C")
  expect_identical(.resolve_theme(theme, "axis.x1.line.col"), "C")
  expect_identical(.resolve_theme(theme, "axis.y1.line.col"), "A")
  expect_identical(.resolve_theme(theme, "axis.x2.line.col"), "B")
  expect_null(.resolve_theme(theme, "axis.y1.nothing"))
})

test_that(".merge_themes gives track-level priority over plot-level", {
  plot_theme  <- list("axis.x.line.col" = "plotval",
                      "axis.y.line.col" = "yellow")
  track_theme <- list("axis.x.line.col" = "trackval")
  merged <- .merge_themes(plot_theme, track_theme)
  expect_identical(merged[["axis.x.line.col"]], "trackval")
  expect_identical(merged[["axis.y.line.col"]], "yellow")
})

test_that(".default_theme yields sensible axis defaults", {
  d <- .default_theme()
  expect_identical(d[["axis.line.col"]], "#1C1B1A")
  expect_identical(d[["axis.x1.position"]], "bottom")
  expect_identical(d[["axis.x2.position"]], "top")
  expect_identical(d[["axis.y1.position"]], "left")
  expect_identical(d[["axis.y2.position"]], "right")
  expect_true(d[["axis.line.visible"]])
})

test_that(".build_axis_spec builds a resolved per-side spec", {
  theme <- .merge_themes(.default_theme(),
                         list("axis.x1.line.col" = "red",
                              "axis.x.line.lwd" = 2))
  spec <- .build_axis_spec(theme, "x1")
  expect_identical(spec$side, "x1")
  expect_identical(spec$axis_dim, "x")
  expect_identical(spec$axis_index, 1L)
  expect_identical(spec$position, "bottom")
  expect_identical(spec$line$col, "red")
  expect_identical(spec$line$lwd, 2)
})

test_that(".build_resolved_theme assembles all four axes", {
  flat <- .default_theme()
  rt <- .build_resolved_theme(flat)
  expect_true(all(c("x1","x2","y1","y2") %in% names(rt$axes)))
  expect_identical(rt$axes$x1$position, "bottom")
  expect_identical(rt$axes$x2$position, "top")
  expect_false(rt$y_per_window)
  expect_true(is.list(rt$chrome))
})
