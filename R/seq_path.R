# ── SeqPathR6 ─────────────────────────────────────────────────────────────────
#
# Primitive element: polyline connecting observations in order, optionally
# split into separate paths via a grouping column. Within each group, points
# are drawn in input order.

#' SeqPath R6 class
#'
#' Internal R6 generator backing [seq_path()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqPathR6 <- R6::R6Class("SeqPath",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqPathR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `x`, `y`. Optional: `group`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip, transform, partition by `group`.
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
      g_vals <- self$resolved$group %||% rep(1L, length(eff_data))

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]
        idx <- qh[mask]

        x_sub <- x_vals[idx]
        y_sub <- y_vals[idx]
        g_sub <- g_vals[idx]

        u <- pmax(pmin((x_sub - panel_meta$xscale[1]) /
                         diff(panel_meta$xscale), 1), 0)
        v <- pmax(pmin((y_sub - panel_meta$yscale[1]) /
                         diff(panel_meta$yscale), 1), 0)
        x_canvas <- panel_meta$inner$x0 + u *
          (panel_meta$inner$x1 - panel_meta$inner$x0)
        y_canvas <- panel_meta$inner$y0 + v *
          (panel_meta$inner$y1 - panel_meta$inner$y0)

        grp_f <- factor(g_sub, levels = unique(g_sub))
        id    <- as.integer(grp_f)

        self$coordCanvas[[w]] <- list(
          x = x_canvas, y = y_canvas, id = id
        )
      }
      invisible()
    },

    #' @description Draw polylines via `grid::grid.polyline()` with `id`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) < 2) next
        grid::grid.polyline(
          x  = grid::unit(coords$x, "npc"),
          y  = grid::unit(coords$y, "npc"),
          id = coords$id,
          gp = gp
        )
      }
    }
  )
)

#' Draw connected paths
#'
#' Primitive element drawing a polyline connecting observations in order,
#' optionally partitioned into separate paths by a grouping column.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y`. Optional:
#'   `group` (column whose values split observations into separate paths).
#' @param aesthetics Optional `SeqAes`. Supports `color`, `linewidth`,
#'   `linetype`, `alpha`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqPathR6` instance.
#' @examples
#' seq_path(map(x = mid, y = score))
#' @export
seq_path <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqPathR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
