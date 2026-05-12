# ── SeqLollipopR6 ─────────────────────────────────────────────────────────────
#
# Composite element: vertical stem from a baseline to `y`, terminated by
# a point at `y`.

#' SeqLollipop R6 class
#'
#' Internal R6 generator backing [seq_lollipop()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqLollipopR6 <- R6::R6Class("SeqLollipop",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqLollipopR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Recognised: `x`, `y`.
    #' @param aesthetics Optional `SeqAes`. `baseline` sets the stem's lower
    #'   y-value (default 0). Also supports `color`, `linewidth`, `size`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip to windows, transform to npc. Builds
    #'   per-window `list(x, y0, y1)` used by both the stem and the head.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      n        <- length(eff_data)
      x_vals   <- self$resolved$x %||% BiocGenerics::start(eff_data)
      y_vals   <- self$resolved$y %||% rep(0.5, n)
      baseline <- self$aesthetics$baseline %||% 0

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        idx  <- qh[mask]
        if (length(idx) == 0) next
        pm <- layout_track[[w]]

        xs <- x_vals[idx]
        ys <- y_vals[idx]

        u  <- pmax(pmin((xs       - pm$xscale[1]) / diff(pm$xscale), 1), 0)
        v1 <- pmax(pmin((ys       - pm$yscale[1]) / diff(pm$yscale), 1), 0)
        v0 <- pmax(pmin((baseline - pm$yscale[1]) / diff(pm$yscale), 1), 0)

        xw <- pm$inner$x1 - pm$inner$x0
        yw <- pm$inner$y1 - pm$inner$y0
        x_c  <- pm$inner$x0 + u  * xw
        y1_c <- pm$inner$y0 + v1 * yw
        y0_c <- pm$inner$y0 + v0 * yw

        ylow  <- pmin(y0_c, y1_c)
        yhigh <- pmax(y0_c, y1_c)

        self$coordCanvas[[w]] <- list(x = x_c, y0 = ylow, y1 = yhigh)
      }
      invisible()
    },

    #' @description Draw stems via `grid::grid.segments()` and heads via
    #'   `grid::grid.points()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp  <- .aes_to_gpar(self$aesthetics)
      pch <- self$aesthetics$shape %||% 16
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) == 0) next
        grid::grid.segments(
          x0 = grid::unit(coords$x,  "npc"),
          x1 = grid::unit(coords$x,  "npc"),
          y0 = grid::unit(coords$y0, "npc"),
          y1 = grid::unit(coords$y1, "npc"),
          gp = gp
        )
        grid::grid.points(
          x   = grid::unit(coords$x,  "npc"),
          y   = grid::unit(coords$y1, "npc"),
          pch = pch,
          gp  = gp
        )
      }
    }
  )
)

#' Draw lollipops
#'
#' Composite element drawing a vertical stem from `baseline` (default 0) to
#' `y`, with a point placed at `y`.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Recognised fields: `x`, `y`.
#' @param aesthetics Optional `SeqAes`. Supports `baseline`, `color`,
#'   `linewidth`, `size`, `shape`, `alpha`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqLollipopR6` instance.
#' @examples
#' seq_lollipop(map(x = start, y = score))
#' @export
seq_lollipop <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqLollipopR6$new(data = data, mapping = mapping,
                    aesthetics = aesthetics, ...)
}
