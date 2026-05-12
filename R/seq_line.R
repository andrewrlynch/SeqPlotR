# ── SeqLineR6 ─────────────────────────────────────────────────────────────────
#
# Primitive element: polyline through (x, y). x defaults to interval midpoint
# `(start + end) / 2`; y defaults to 0.5. `aes(type = "step")` renders a
# step-line (x = start of each interval, duplicated with rep-each-2 trick
# ported from THEfunc).

#' SeqLine R6 class
#'
#' Internal R6 generator backing [seq_line()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqLineR6 <- R6::R6Class("SeqLine",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqLineR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`. `type = "step"` (or `"s"`) enables
    #'   step-line rendering.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip to windows, transform to npc.
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

      is_step <- isTRUE(self$aesthetics$type %in% c("s", "step"))

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]

        qh_w <- qh[mask]
        x_sub <- x_vals[qh_w]
        y_sub <- y_vals[qh_w]

        if (is_step) {
          # Port from THEfunc SeqLine: left edge of each interval, duplicated
          # so consecutive pairs form horizontal treads.
          x_sub <- BiocGenerics::start(eff_data)[qh_w]
          ord <- order(x_sub)
          x_sub <- x_sub[ord]
          y_sub <- y_sub[ord]
          x_sub <- rep(x_sub, each = 2)[-1]
          n2    <- length(y_sub) * 2
          y_sub <- rep(y_sub, each = 2)[seq_len(max(n2 - 1, 0))]
        }

        xpr <- panel_meta$xplot_range %||% panel_meta$xscale
        ypr <- panel_meta$yplot_range %||% panel_meta$yscale
        oob <- .apply_oob(x_sub, y_sub, xpr, ypr,
                          mode = panel_meta$x_oob %||% "exclude",
                          label = "seq_line")
        x_sub <- oob$x; y_sub <- oob$y
        u <- (x_sub - xpr[1]) / diff(xpr)
        v <- (y_sub - ypr[1]) / diff(ypr)
        self$coordCanvas[[w]] <- list(
          x = panel_meta$inner$x0 + u *
            (panel_meta$inner$x1 - panel_meta$inner$x0),
          y = panel_meta$inner$y0 + v *
            (panel_meta$inner$y1 - panel_meta$inner$y0)
        )
      }

      # Auto-legend for color mapping
      if (is.null(self$legend) && isTRUE(self$show_legend) &&
          !is.null(self$resolved$color)) {
        col_nm     <- .map_col_name(self$resolved$.mapping, "color")
        col_result <- .auto_scale_colors(self$resolved$color, col_name = col_nm,
                                         aes_name = "color")
        self$auto_legend <- col_result$legend
      }

      invisible()
    },

    #' @description Draw the polyline with `grid::grid.lines()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) < 2) next
        grid::grid.lines(
          x  = grid::unit(coords$x, "npc"),
          y  = grid::unit(coords$y, "npc"),
          gp = gp
        )
      }
    }
  )
)

#' Draw a line
#'
#' Primitive element drawing a polyline through `(x, y)` for each interval.
#' `x` defaults to the interval midpoint; `y` defaults to 0.5.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Recognized fields: `x`, `y`.
#' @param aesthetics Optional `SeqAes`. Recognized fields: `color`,
#'   `linewidth`, `linetype`, `alpha`, and `type` (`"step"` or `"s"` for
#'   step lines).
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqLineR6` instance.
#' @examples
#' seq_line(map(x = mid, y = score))
#' @export
seq_line <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqLineR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
