# ── LegendKey ─────────────────────────────────────────────────────────────────
#
# S3 constructor for a single legend-key entry. Stores the visual properties
# and label for one row in a legend group.  Rendering is handled elsewhere.

#' Legend Key
#'
#' Constructs a single legend key entry for a SeqPlotR element.  A `LegendKey`
#' records the visual properties (colour, shape, fill, etc.), the display
#' label, and an optional group title for one row in a rendered legend.
#'
#' @param label Character or `NULL`. Text label shown beside the key glyph.
#' @param title Character or `NULL`. Legend group title this key belongs to.
#' @param color Character. Stroke / point colour. Default `"#1C1B1A"`.
#' @param shape Character. Glyph shape code (e.g. `"-"`, `"circle"`).
#'   Default `"-"`.
#' @param size Numeric. Relative glyph size. Default `1`.
#' @param alpha Numeric in `[0, 1]`. Opacity. Default `1`.
#' @param fill Character or `NULL`. Fill colour for filled glyphs.
#'   Default `NULL`.
#' @param lty Integer or character. Line type. Default `1`.
#' @param ... Additional fields stored verbatim in the `extra` sub-list.
#'
#' @return A `LegendKey` S3 object (a named list with class `"LegendKey"`).
#' @export
LegendKey <- function(label = NULL,
                      title = NULL,
                      color = "#1C1B1A",
                      shape = "-",
                      size  = 1,
                      alpha = 1,
                      fill  = NULL,
                      lty   = 1,
                      ...) {
  key <- list(label = label, title = title, color = color, shape = shape,
              size  = size,  alpha = alpha, fill  = fill,  lty   = lty,
              extra = list(...))
  class(key) <- "LegendKey"
  key
}
