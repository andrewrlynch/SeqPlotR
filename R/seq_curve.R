# ── SeqCurveR6 ────────────────────────────────────────────────────────────────
#
# Primitive element: cubic bezier curve from (x, y) to (x_end, y_end) with a
# y-offset `curvature` controlling arch height.

#' SeqCurve R6 class
#'
#' Internal R6 generator backing [seq_curve()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqCurveR6 <- R6::R6Class("SeqCurve",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqCurveR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `x`, `y`, `x_end`, `y_end`.
    #' @param aesthetics Optional `SeqAes`. `curvature` (default 0.3) sets
    #'   the fractional y-offset of the bezier control points.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, compute bezier grobs per window.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      x0_vals <- self$resolved$x     %||% BiocGenerics::start(eff_data)
      x1_vals <- self$resolved$x_end %||% BiocGenerics::end(eff_data)
      y0_vals <- self$resolved$y     %||% rep(0.5, length(eff_data))
      y1_vals <- self$resolved$y_end %||% rep(0.5, length(eff_data))

      curvature <- self$aesthetics$curvature %||% 0.3

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]
        idx <- qh[mask]

        xpr <- panel_meta$xplot_range %||% panel_meta$xscale
        ypr <- panel_meta$yplot_range %||% panel_meta$yscale
        u0 <- pmax(pmin((x0_vals[idx] - xpr[1]) / diff(xpr), 1), 0)
        u1 <- pmax(pmin((x1_vals[idx] - xpr[1]) / diff(xpr), 1), 0)
        v0 <- pmax(pmin((y0_vals[idx] - ypr[1]) / diff(ypr), 1), 0)
        v1 <- pmax(pmin((y1_vals[idx] - ypr[1]) / diff(ypr), 1), 0)

        xw <- panel_meta$inner$x1 - panel_meta$inner$x0
        yw <- panel_meta$inner$y1 - panel_meta$inner$y0
        x0c <- panel_meta$inner$x0 + u0 * xw
        x1c <- panel_meta$inner$x0 + u1 * xw
        y0c <- panel_meta$inner$y0 + v0 * yw
        y1c <- panel_meta$inner$y0 + v1 * yw

        cp1_x <- x0c + (x1c - x0c) / 3
        cp2_x <- x0c + 2 * (x1c - x0c) / 3
        cp1_y <- pmin(pmax(y0c + curvature * yw, 0), 1)
        cp2_y <- pmin(pmax(y1c + curvature * yw, 0), 1)

        self$coordCanvas[[w]] <- list(
          x0 = x0c, x1 = x1c, y0 = y0c, y1 = y1c,
          cp1_x = cp1_x, cp2_x = cp2_x,
          cp1_y = cp1_y, cp2_y = cp2_y
        )
      }
      invisible()
    },

    #' @description Draw bezier curves via `grid::bezierGrob()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords)) next
        n <- length(coords$x0)
        for (i in seq_len(n)) {
          g <- grid::bezierGrob(
            x  = grid::unit(c(coords$x0[i], coords$cp1_x[i],
                              coords$cp2_x[i], coords$x1[i]), "npc"),
            y  = grid::unit(c(coords$y0[i], coords$cp1_y[i],
                              coords$cp2_y[i], coords$y1[i]), "npc"),
            gp = gp
          )
          grid::grid.draw(g)
        }
      }
    }
  )
)

#' Draw bezier curves
#'
#' Primitive element drawing a cubic bezier curve from `(x, y)` to
#' `(x_end, y_end)` with a y-offset controlled by `curvature`.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y`, `x_end`,
#'   `y_end`.
#' @param aesthetics Optional `SeqAes`. Supports `curvature` (default 0.3),
#'   `color`, `linewidth`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqCurveR6` instance.
#' @examples
#' seq_curve(map(x = start, y = score, x_end = end, y_end = score))
#' @export
seq_curve <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqCurveR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
