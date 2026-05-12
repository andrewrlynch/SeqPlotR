# ── SeqPolyR6 ─────────────────────────────────────────────────────────────────
#
# Primitive element: filled polygon(s). Vertices come from (x, y); optional
# `group` mapping splits vertices into multiple polygons passed together via
# grid.polygon()'s `id` argument.

#' SeqPoly R6 class
#'
#' Internal R6 generator backing [seq_poly()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqPolyR6 <- R6::R6Class("SeqPoly",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqPolyR6.
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

        xpr <- panel_meta$xplot_range %||% panel_meta$xscale
        ypr <- panel_meta$yplot_range %||% panel_meta$yscale
        oob <- .apply_oob(x_sub, y_sub, xpr, ypr,
                          mode = panel_meta$x_oob %||% "exclude",
                          label = "seq_poly")
        x_sub <- oob$x; y_sub <- oob$y
        g_sub <- g_sub[oob$keep]
        u <- (x_sub - xpr[1]) / diff(xpr)
        v <- (y_sub - ypr[1]) / diff(ypr)
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

    #' @description Draw polygons via `grid::grid.polygon()` with `id`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp <- .aes_to_gpar(self$aesthetics)
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) < 3) next
        grid::grid.polygon(
          x  = grid::unit(coords$x, "npc"),
          y  = grid::unit(coords$y, "npc"),
          id = coords$id,
          gp = gp
        )
      }
    }
  )
)

#' Draw filled polygons
#'
#' Primitive element drawing one or more filled polygons. Vertices come from
#' `(x, y)`; an optional `group` mapping partitions vertices across multiple
#' polygons.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y`. Optional:
#'   `group`.
#' @param aesthetics Optional `SeqAes`. Supports `fill`, `color`, `alpha`,
#'   `linewidth`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqPolyR6` instance.
#' @examples
#' seq_poly(map(x = start, y = score))
#' @export
seq_poly <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqPolyR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
