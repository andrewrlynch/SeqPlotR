# ── SeqSegmentR6 ──────────────────────────────────────────────────────────────
#
# Primitive element: straight line segment from (x, y) to (x_end, y_end).
# Required map() fields: x, x_end, y, y_end. Defaults: x -> start, x_end ->
# end, y0 -> 0, y1 -> resolved y.

#' SeqSegment R6 class
#'
#' Internal R6 generator backing [seq_segment()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqSegmentR6 <- R6::R6Class("SeqSegment",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqSegmentR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required fields: `x`, `x_end`, `y`,
    #'   `y_end`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics)
    },

    #' @description Resolve mappings, apply auto-scale for non-concrete color,
    #'   clip to windows, transform to npc.
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
      y1_vals <- self$resolved$y     %||% rep(0.5, length(eff_data))
      y0_vals <- self$resolved$y_end %||% rep(0,   length(eff_data))

      # ---- Auto-scale color when legend is not manually set ----
      col_result <- list(colors = NULL, legend = NULL)
      if (is.null(self$legend) && isTRUE(self$show_legend) &&
          !is.null(self$resolved$color)) {
        col_nm     <- .map_col_name(self$resolved$.mapping, "color")
        col_result <- .auto_scale_colors(self$resolved$color, col_name = col_nm,
                                         aes_name = "color")
        self$auto_legend <- col_result$legend
      }

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
        mode <- panel_meta$x_oob %||% "exclude"
        # A segment is OOB iff *both* anchors fall outside the plot range.
        x_in <- (pmax(x0_vals[idx], x1_vals[idx]) >= xpr[1]) &
                (pmin(x0_vals[idx], x1_vals[idx]) <= xpr[2])
        y_in <- (pmax(y0_vals[idx], y1_vals[idx]) >= ypr[1]) &
                (pmin(y0_vals[idx], y1_vals[idx]) <= ypr[2])
        in_range <- x_in & y_in
        n_oob <- sum(!in_range, na.rm = TRUE)
        if (n_oob > 0L && !isTRUE(getOption("seqplotr.suppress_oob", FALSE))) {
          verb <- if (mode == "exclude") "excluded" else "plotted"
          message(n_oob, " out-of-bounds data points ", verb, "! (seq_segment)")
        }
        if (mode == "exclude") {
          idx <- idx[in_range]
        }
        x0i <- x0_vals[idx]; x1i <- x1_vals[idx]
        y0i <- y0_vals[idx]; y1i <- y1_vals[idx]
        u0 <- pmax(pmin((x0i - xpr[1]) / diff(xpr), 1), 0)
        u1 <- pmax(pmin((x1i - xpr[1]) / diff(xpr), 1), 0)
        v0 <- pmax(pmin((y0i - ypr[1]) / diff(ypr), 1), 0)
        v1 <- pmax(pmin((y1i - ypr[1]) / diff(ypr), 1), 0)

        xw <- panel_meta$inner$x1 - panel_meta$inner$x0
        yw <- panel_meta$inner$y1 - panel_meta$inner$y0
        self$coordCanvas[[w]] <- list(
          x0    = panel_meta$inner$x0 + u0 * xw,
          x1    = panel_meta$inner$x0 + u1 * xw,
          y0    = panel_meta$inner$y0 + v0 * yw,
          y1    = panel_meta$inner$y0 + v1 * yw,
          color = if (!is.null(col_result$colors))
                    as.character(col_result$colors[idx])
                  else NULL
        )
      }
      invisible()
    },

    #' @description Draw segments via `grid::grid.segments()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp_base <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords)) next
        # Per-segment colors override the constant aesthetic color
        col_use <- coords$color %||% gp_base$col
        gp <- grid::gpar(
          col   = col_use,
          lwd   = gp_base$lwd,
          lty   = gp_base$lty,
          alpha = gp_base$alpha
        )
        grid::grid.segments(
          x0 = grid::unit(coords$x0, "npc"),
          x1 = grid::unit(coords$x1, "npc"),
          y0 = grid::unit(coords$y0, "npc"),
          y1 = grid::unit(coords$y1, "npc"),
          gp = gp
        )
      }
    }
  )
)

#' Draw line segments
#'
#' Primitive element drawing a straight segment from `(x, y)` to
#' `(x_end, y_end)` for each observation. When `map(color = col)` is used
#' with a non-color column, colors are auto-scaled (discrete → flexoki palette;
#' numeric → viridis gradient) and an auto-legend is generated.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `x_end`, `y`,
#'   `y_end`. Defaults: `x = start`, `x_end = end`, `y0 = 0`,
#'   `y1 = resolved y`. `color` is auto-scaled when mapped to a non-color
#'   column.
#' @param aesthetics Optional `SeqAes`. Supports `color`, `linewidth`,
#'   `linetype`, `alpha`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqSegmentR6` instance.
#' @examples
#' seq_segment(map(x = start, x_end = end, y = score, y_end = score))
#' @export
seq_segment <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqSegmentR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
