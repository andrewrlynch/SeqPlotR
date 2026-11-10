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
    #'   per-window `list(x, y0, y1, y_head, color, pch)`. `y0`/`y1` span the
    #'   stem (baseline to value); `y_head` is the value's own npc position so
    #'   the head renders at the value even when it is below the baseline.
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

      # ---- Auto-scale color / fill / shape (mirrors seq_point) ----
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
        auto_legs <- Filter(Negate(is.null),
                            list(col_result$legend, fill_result$legend,
                                 shp_result$legend))
        self$auto_legend <- if (length(auto_legs) == 1L) auto_legs[[1L]]
                            else if (length(auto_legs) > 1L) auto_legs
                            else NULL
      }

      col_vals  <- col_result$colors
      # When color maps directly to a concrete colour (no auto-scale), carry it
      # through per-point so map(color = <colour column>) still works.
      if (is.null(col_vals) && !is.null(self$resolved$color) &&
          is.character(self$resolved$color))
        col_vals <- self$resolved$color
      shp_vals  <- shp_result$shapes

      shape_to_pch <- c(circle = 16L, square = 15L, triangle = 17L, diamond = 18L)
      pch_vals <- if (!is.null(shp_vals)) {
        pch <- shape_to_pch[shp_vals]; pch[is.na(pch)] <- 16L; unname(pch)
      } else NULL

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

        xpr <- pm$xplot_range %||% pm$xscale
        ypr <- pm$yplot_range %||% pm$yscale
        u  <- pmax(pmin((xs       - xpr[1]) / diff(xpr), 1), 0)
        v1 <- pmax(pmin((ys       - ypr[1]) / diff(ypr), 1), 0)
        v0 <- pmax(pmin((baseline - ypr[1]) / diff(ypr), 1), 0)

        xw <- pm$inner$x1 - pm$inner$x0
        yw <- pm$inner$y1 - pm$inner$y0
        x_c  <- pm$inner$x0 + u  * xw
        y1_c <- pm$inner$y0 + v1 * yw   # value's npc position (head sits here)
        y0_c <- pm$inner$y0 + v0 * yw   # baseline's npc position

        # Stem spans baseline<->value; the head is anchored at the value.
        ylow  <- pmin(y0_c, y1_c)
        yhigh <- pmax(y0_c, y1_c)

        self$coordCanvas[[w]] <- list(
          x      = x_c,
          y0     = ylow,
          y1     = yhigh,
          y_head = y1_c,
          color  = if (!is.null(col_vals)) as.character(col_vals[idx]) else NULL,
          pch    = if (!is.null(pch_vals)) pch_vals[idx]               else NULL
        )
      }
      invisible()
    },

    #' @description Draw stems via `grid::grid.segments()` and heads via
    #'   `grid::grid.points()`, honouring any per-point `color`/`shape`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      pch_def   <- self$aesthetics$shape %||% 16
      col_def   <- self$aesthetics$color %||% self$aesthetics$col
      lwd_stem  <- self$aesthetics$linewidth %||% self$aesthetics$lwd
      cex_head  <- self$aesthetics$size
      alpha     <- self$aesthetics$alpha
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) == 0) next
        pt_col  <- coords$color %||% col_def
        pt_pch  <- coords$pch   %||% pch_def
        grid::grid.segments(
          x0 = grid::unit(coords$x,  "npc"),
          x1 = grid::unit(coords$x,  "npc"),
          y0 = grid::unit(coords$y0, "npc"),
          y1 = grid::unit(coords$y1, "npc"),
          gp = grid::gpar(col = pt_col, lwd = lwd_stem, alpha = alpha)
        )
        grid::grid.points(
          x   = grid::unit(coords$x,      "npc"),
          y   = grid::unit(coords$y_head, "npc"),
          pch = pt_pch,
          gp  = grid::gpar(col = pt_col, fill = pt_col,
                           cex = cex_head, alpha = alpha)
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
#' @param mapping Optional `SeqMap`. Recognised fields: `x`, `y`, and
#'   `color`/`fill`/`shape` (auto-scaled to a palette with a legend when
#'   mapped to a non-colour column, exactly as in [seq_point()]). Values
#'   below `baseline` place the head at the value, not at the baseline.
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
