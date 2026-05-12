# ── SeqPointR6 ────────────────────────────────────────────────────────────────
#
# Primitive element: one glyph per observation at (x, y). Draws via
# grid::grid.points(). x defaults to start(data); y defaults to 0.5.

#' SeqPoint R6 class
#'
#' Internal R6 generator backing [seq_point()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqPointR6 <- R6::R6Class("SeqPoint",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqPointR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, apply auto-scales for non-concrete
    #'   color/fill/shape, clip to windows, transform genomic coordinates to
    #'   canvas npc.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      x_vals <- self$resolved$x %||% BiocGenerics::start(eff_data)
      y_vals <- self$resolved$y %||% rep(0.5, length(eff_data))

      # ---- Auto-scale color / fill / shape when legend is not manually set ----
      col_result  <- list(colors = NULL, legend = NULL)
      fill_result <- list(colors = NULL, legend = NULL)
      shp_result  <- list(shapes = NULL, legend = NULL)

      if (is.null(self$legend) && isTRUE(self$show_legend)) {
        raw_mapping <- self$resolved$.mapping

        if (!is.null(self$resolved$color)) {
          col_nm <- .map_col_name(raw_mapping, "color")
          col_result <- .auto_scale_colors(self$resolved$color, col_name = col_nm,
                                           aes_name = "color")
        }
        if (!is.null(self$resolved$fill)) {
          fill_nm <- .map_col_name(raw_mapping, "fill")
          fill_result <- .auto_scale_colors(self$resolved$fill, col_name = fill_nm,
                                            aes_name = "fill")
        }
        if (!is.null(self$resolved$shape)) {
          shp_nm <- .map_col_name(raw_mapping, "shape")
          shp_result <- .auto_scale_shapes(self$resolved$shape, col_name = shp_nm)
        }

        # Collect non-NULL legends in a list; store as single spec when only one
        auto_legs <- Filter(Negate(is.null),
                            list(col_result$legend, fill_result$legend,
                                 shp_result$legend))
        self$auto_legend <- if (length(auto_legs) == 1L) auto_legs[[1L]]
                            else if (length(auto_legs) > 1L) auto_legs
                            else NULL
      }

      col_vals  <- col_result$colors
      fill_vals <- fill_result$colors
      shp_vals  <- shp_result$shapes

      # ---- Convert shape strings to pch integers for grid.points() ----
      shape_to_pch <- c(circle = 16L, square = 15L, triangle = 17L, diamond = 18L)
      pch_vals <- if (!is.null(shp_vals)) {
        pch <- shape_to_pch[shp_vals]
        pch[is.na(pch)] <- 16L
        unname(pch)
      } else NULL

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]
        idx   <- qh[mask]
        x_sub <- x_vals[idx]
        y_sub <- y_vals[idx]
        u <- pmax(pmin((x_sub - panel_meta$xscale[1]) /
                         diff(panel_meta$xscale), 1), 0)
        v <- pmax(pmin((y_sub - panel_meta$yscale[1]) /
                         diff(panel_meta$yscale), 1), 0)
        self$coordCanvas[[w]] <- list(
          x    = panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0),
          y    = panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0),
          color = if (!is.null(col_vals))  as.character(col_vals[idx])  else NULL,
          fill  = if (!is.null(fill_vals)) as.character(fill_vals[idx]) else NULL,
          pch   = if (!is.null(pch_vals))  pch_vals[idx]                else NULL
        )
      }
      invisible()
    },

    #' @description Draw points via `grid::grid.points()` for each window.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      pch_default <- self$aesthetics$shape %||% 16
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) == 0L) next
        pch_use <- coords$pch %||% pch_default
        gp <- grid::gpar(
          col      = coords$color %||% self$aesthetics$color %||%
                     self$aesthetics$col,
          fill     = coords$fill  %||% self$aesthetics$fill,
          lwd      = self$aesthetics$linewidth %||% self$aesthetics$lwd,
          lty      = self$aesthetics$linetype,
          alpha    = self$aesthetics$alpha,
          cex      = self$aesthetics$size,
          fontsize = self$aesthetics$fontsize
        )
        grid::grid.points(
          x   = grid::unit(coords$x, "npc"),
          y   = grid::unit(coords$y, "npc"),
          pch = pch_use,
          gp  = gp
        )
      }
    }
  )
)

#' Draw points
#'
#' Primitive element drawing one glyph per observation at `(x, y)`.
#'
#' @param data Optional `GRanges`. Falls back to the parent track's data.
#' @param mapping Optional `SeqMap` from [map()]. Required fields: `x`, `y`.
#'   `x` defaults to `start(data)`; `y` defaults to 0.5 when unmapped.
#'   `color`, `fill`, and `shape` are auto-scaled when mapped to non-color
#'   columns: discrete columns get the flexoki palette; numeric columns get
#'   a viridis gradient legend.
#' @param aesthetics Optional `SeqAes` from [aes()]. Supports `color`, `fill`,
#'   `size`, `shape`, `alpha`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqPointR6` instance (S3 class `"SeqPoint" / "SeqElement"`).
#' @examples
#' seq_point(map(x = start, y = score))
#' seq_point(map(x = start, y = score, color = type))  # auto-legend
#' @export
seq_point <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqPointR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
