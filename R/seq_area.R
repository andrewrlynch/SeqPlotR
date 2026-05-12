# ── SeqAreaR6 ─────────────────────────────────────────────────────────────────
#
# Primitive element: filled polygon formed by the data line (x, y) followed by
# a baseline return. No stacking in this batch — grouped/stacked areas are
# deferred. `baseline` is a constant from aes() (default 0).

#' SeqArea R6 class
#'
#' Internal R6 generator backing [seq_area()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqAreaR6 <- R6::R6Class("SeqArea",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqAreaR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `x`, `y`.
    #' @param aesthetics Optional `SeqAes`. `baseline` sets the closing
    #'   y-value (default 0).
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip to windows, transform to npc, then
    #'   build a closed polygon by appending a reversed baseline path.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      x_vals <- self$resolved$x %||%
        ((BiocGenerics::start(eff_data) + BiocGenerics::end(eff_data)) / 2)
      y_vals <- self$resolved$y %||% rep(0.5, length(eff_data))

      baseline <- self$aesthetics$baseline %||% 0

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]
        idx <- qh[mask]
        if (length(idx) == 0) next

        # Sort along the genomic (positional) axis so the outline traces
        # the data curve monotonically. Falls back to x for conventional
        # tracks.
        y_is_genomic <- isTRUE(panel_meta$y_is_genomic)
        sort_axis <- if (y_is_genomic) "y" else "x"
        ord <- if (sort_axis == "y") order(y_vals[idx]) else order(x_vals[idx])
        idx <- idx[ord]

        u <- pmax(pmin((x_vals[idx] - panel_meta$xscale[1]) /
                         diff(panel_meta$xscale), 1), 0)
        v <- pmax(pmin((y_vals[idx] - panel_meta$yscale[1]) /
                         diff(panel_meta$yscale), 1), 0)
        xw <- panel_meta$inner$x1 - panel_meta$inner$x0
        yw <- panel_meta$inner$y1 - panel_meta$inner$y0
        x_canvas <- panel_meta$inner$x0 + u * xw
        y_canvas <- panel_meta$inner$y0 + v * yw

        if (sort_axis == "y") {
          # Flipped orientation: baseline is a value on the scalar x-axis,
          # closing the polygon with a vertical return path.
          base_u <- pmax(pmin((baseline - panel_meta$xscale[1]) /
                                diff(panel_meta$xscale), 1), 0)
          base_canvas <- panel_meta$inner$x0 + base_u * xw
          x_poly <- c(x_canvas, rep(base_canvas, length(x_canvas)))
          y_poly <- c(y_canvas, rev(y_canvas))
        } else {
          # Conventional orientation: baseline on the scalar y-axis,
          # horizontal return path.
          base_v <- pmax(pmin((baseline - panel_meta$yscale[1]) /
                                diff(panel_meta$yscale), 1), 0)
          base_canvas <- panel_meta$inner$y0 + base_v * yw
          x_poly <- c(x_canvas, rev(x_canvas))
          y_poly <- c(y_canvas, rep(base_canvas, length(y_canvas)))
        }

        self$coordCanvas[[w]] <- list(x = x_poly, y = y_poly)
      }

      # Auto-legend for fill or color mapping
      if (is.null(self$legend) && isTRUE(self$show_legend)) {
        if (!is.null(self$resolved$fill)) {
          fill_nm     <- .map_col_name(self$resolved$.mapping, "fill")
          fill_result <- .auto_scale_colors(self$resolved$fill, col_name = fill_nm,
                                            aes_name = "fill")
          self$auto_legend <- fill_result$legend
        } else if (!is.null(self$resolved$color)) {
          col_nm     <- .map_col_name(self$resolved$.mapping, "color")
          col_result <- .auto_scale_colors(self$resolved$color, col_name = col_nm,
                                           aes_name = "color")
          self$auto_legend <- col_result$legend
        }
      }

      invisible()
    },

    #' @description Draw filled polygons via `grid::grid.polygon()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) < 3) next
        grid::grid.polygon(
          x  = grid::unit(coords$x, "npc"),
          y  = grid::unit(coords$y, "npc"),
          gp = gp
        )
      }
    }
  )
)

#' Draw filled area under a curve
#'
#' Primitive element drawing a filled polygon formed by the data line
#' `(x, y)` and a return path along `baseline`.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y`.
#' @param aesthetics Optional `SeqAes`. Supports `baseline` (default 0),
#'   `fill`, `color`, `alpha`, `linewidth`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqAreaR6` instance.
#' @examples
#' seq_area(map(x = mid, y = score))
#' @export
seq_area <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqAreaR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
