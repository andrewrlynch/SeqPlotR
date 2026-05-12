# ── SeqScalebarR6 ─────────────────────────────────────────────────────────────
#
# Decorative reference scalebar element. Independent of any axis — useful
# when the x-axis is hidden but the user still wants a visual length cue.

#' Auto-format a scalebar label from a base-pair length
#' @param bp Numeric. Length in bp.
#' @return Character label such as `"50 kb"`.
#' @keywords internal
.scalebar_auto_label <- function(bp) {
  if      (bp >= 1e6) sprintf("%g Mb", bp / 1e6)
  else if (bp >= 1e3) sprintf("%g kb", bp / 1e3)
  else                sprintf("%g bp", bp)
}

#' SeqScalebar R6 class
#'
#' Internal R6 generator backing [seq_scalebar()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqScalebarR6 <- R6::R6Class("SeqScalebar",
  inherit = SeqElementR6,
  public = list(
    #' @field length_bp Length of the bar in base pairs.
    length_bp   = NULL,
    #' @field label Character label rendered above the bar.
    label       = NULL,
    #' @field hjust Horizontal anchor of the bar's right edge in panel NPC.
    hjust       = 0.95,
    #' @field vjust Vertical anchor of the bar in panel NPC.
    vjust       = 0.05,
    #' @field ticks Whether to draw tick caps at the bar's ends.
    ticks       = TRUE,
    #' @field tick_height Tick height as a fraction of the panel height.
    tick_height = 0.15,
    #' @field bar_lwd Line width of the bar and ticks.
    bar_lwd     = 1.2,
    #' @field bar_col Bar / tick color.
    bar_col     = "#1C1B1A",
    #' @field label_cex Label `cex`.
    label_cex   = 0.65,
    #' @field label_col Label color.
    label_col   = "#1C1B1A",
    #' @field label_pad Vertical NPC gap between bar top and label bottom.
    label_pad   = 0.005,

    #' @description Construct a SeqScalebarR6.
    #' @param length_bp Numeric, positive.
    #' @param label Optional character label. When `NULL`, auto-formatted.
    #' @param hjust Numeric in 0--1. Default `0.95`.
    #' @param vjust Numeric in 0--1. Default `0.05`.
    #' @param ticks Logical. Default `TRUE`.
    #' @param tick_height Numeric. Default `0.15`.
    #' @param bar_lwd Numeric. Default `1.2`.
    #' @param bar_col Character. Default `"#1C1B1A"`.
    #' @param label_cex Numeric. Default `0.65`.
    #' @param label_col Character. Default `"#1C1B1A"`.
    #' @param label_pad Numeric. Default `0.005`.
    #' @param ... Reserved.
    initialize = function(length_bp,
                          label       = NULL,
                          hjust       = 0.95,
                          vjust       = 0.05,
                          ticks       = TRUE,
                          tick_height = 0.15,
                          bar_lwd     = 1.2,
                          bar_col     = "#1C1B1A",
                          label_cex   = 0.65,
                          label_col   = "#1C1B1A",
                          label_pad   = 0.005,
                          ...) {
      super$initialize()
      if (!is.numeric(length_bp) || length(length_bp) != 1L ||
          !is.finite(length_bp) || length_bp <= 0)
        stop("`length_bp` must be a positive number.", call. = FALSE)
      self$length_bp   <- length_bp
      self$label       <- label %||% .scalebar_auto_label(length_bp)
      self$hjust       <- hjust
      self$vjust       <- vjust
      self$ticks       <- isTRUE(ticks)
      self$tick_height <- tick_height
      self$bar_lwd     <- bar_lwd
      self$bar_col     <- bar_col
      self$label_cex   <- label_cex
      self$label_col   <- label_col
      self$label_pad   <- label_pad
    },

    #' @description Compute per-window NPC coordinates for the bar, ticks,
    #'   and label.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      panels <- vector("list", length(track_windows))

      for (w in seq_along(track_windows)) {
        pm     <- layout_track[[w]]
        xscale <- pm$xscale
        win_bp <- diff(xscale)
        pan_w  <- pm$inner$x1 - pm$inner$x0
        pan_h  <- pm$inner$y1 - pm$inner$y0

        bar_frac <- self$length_bp / win_bp
        bar_w    <- bar_frac * pan_w

        x_right <- pm$inner$x0 + self$hjust * pan_w
        x_left  <- x_right - bar_w
        y_bar   <- pm$inner$y0 + self$vjust * pan_h

        tick_h  <- self$tick_height * pan_h / 2

        panels[[w]] <- list(
          x0      = x_left,   x1 = x_right,
          y       = y_bar,
          yt0     = y_bar - tick_h, yt1 = y_bar + tick_h,
          label_x = (x_left + x_right) / 2,
          label_y = y_bar + tick_h + self$label_pad
        )
      }
      self$coordCanvas <- panels
      invisible()
    },

    #' @description Draw the bar, tick caps, and label.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())

      for (coords in self$coordCanvas) {
        if (is.null(coords)) next

        grid::grid.lines(
          x  = grid::unit(c(coords$x0, coords$x1), "npc"),
          y  = grid::unit(c(coords$y, coords$y),   "npc"),
          gp = grid::gpar(col = self$bar_col, lwd = self$bar_lwd,
                          lineend = "square")
        )

        if (isTRUE(self$ticks)) {
          for (tx in c(coords$x0, coords$x1)) {
            grid::grid.lines(
              x  = grid::unit(c(tx, tx), "npc"),
              y  = grid::unit(c(coords$yt0, coords$yt1), "npc"),
              gp = grid::gpar(col = self$bar_col, lwd = self$bar_lwd,
                              lineend = "square")
            )
          }
        }

        if (!is.null(self$label) && nzchar(self$label)) {
          grid::grid.text(
            label = self$label,
            x     = grid::unit(coords$label_x, "npc"),
            y     = grid::unit(coords$label_y, "npc"),
            just  = c("center", "bottom"),
            gp    = grid::gpar(cex = self$label_cex, col = self$label_col)
          )
        }
      }
      invisible()
    }
  )
)

#' Reference scalebar element
#'
#' Draws a horizontal scalebar of a specified genomic length with tick
#' caps and a label. Useful as a visual length reference when the x-axis
#' is hidden.
#'
#' @param length_bp Numeric. Length of the scalebar in base pairs.
#'   Required.
#' @param label Character or `NULL`. Bar label. When `NULL`, auto-
#'   formatted from `length_bp` (e.g. `"50 kb"`).
#' @param hjust Numeric in 0–1. Horizontal position of the bar's right
#'   edge within the panel. Default `0.95`.
#' @param vjust Numeric in 0–1. Vertical position of the bar within the
#'   panel. Default `0.05` (near the bottom).
#' @param ticks Logical. Draw tick caps at each end. Default `TRUE`.
#' @param tick_height Numeric. Tick height as a fraction of panel height.
#'   Default `0.15`.
#' @param bar_lwd Numeric. Line width of the bar and ticks. Default `1.2`.
#' @param bar_col Character. Bar / tick color. Default `"#1C1B1A"`.
#' @param label_cex Numeric. Label character expansion. Default `0.65`.
#' @param label_col Character. Label color. Default `"#1C1B1A"`.
#' @param label_pad Numeric. NPC gap between bar top and label bottom.
#'   Default `0.005`.
#' @return A `SeqScalebarR6` instance.
#' @export
seq_scalebar <- function(length_bp,
                         label       = NULL,
                         hjust       = 0.95,
                         vjust       = 0.05,
                         ticks       = TRUE,
                         tick_height = 0.15,
                         bar_lwd     = 1.2,
                         bar_col     = "#1C1B1A",
                         label_cex   = 0.65,
                         label_col   = "#1C1B1A",
                         label_pad   = 0.005) {
  SeqScalebarR6$new(
    length_bp   = length_bp,
    label       = label,
    hjust       = hjust,
    vjust       = vjust,
    ticks       = ticks,
    tick_height = tick_height,
    bar_lwd     = bar_lwd,
    bar_col     = bar_col,
    label_cex   = label_cex,
    label_col   = label_col,
    label_pad   = label_pad
  )
}
