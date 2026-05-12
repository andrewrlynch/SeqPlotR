# ── SeqTextR6 ─────────────────────────────────────────────────────────────────
#
# Primitive element: text labels at (x, y). label may come from map(label =
# col) (per-observation) or aes(label = "fixed") (recycled to every point).

#' SeqText R6 class
#'
#' Internal R6 generator backing [seq_text()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqTextR6 <- R6::R6Class("SeqText",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqTextR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `x`, `y`, `label`.
    #' @param aesthetics Optional `SeqAes`. Supports `size`, `color`, `angle`,
    #'   `hjust`, `vjust`, and a constant `label`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip, transform to npc.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      n <- length(eff_data)
      x_vals <- self$resolved$x %||% BiocGenerics::start(eff_data)
      y_vals <- self$resolved$y %||% rep(0.5, n)

      lbl_vals <- self$resolved$label
      if (is.null(lbl_vals) && !is.null(self$aesthetics$label))
        lbl_vals <- rep(as.character(self$aesthetics$label), length.out = n)
      if (is.null(lbl_vals))
        lbl_vals <- rep("", n)
      lbl_vals <- as.character(lbl_vals)

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        panel_meta <- layout_track[[w]]
        idx <- qh[mask]

        u <- pmax(pmin((x_vals[idx] - panel_meta$xscale[1]) /
                         diff(panel_meta$xscale), 1), 0)
        v <- pmax(pmin((y_vals[idx] - panel_meta$yscale[1]) /
                         diff(panel_meta$yscale), 1), 0)
        self$coordCanvas[[w]] <- list(
          x     = panel_meta$inner$x0 + u *
                    (panel_meta$inner$x1 - panel_meta$inner$x0),
          y     = panel_meta$inner$y0 + v *
                    (panel_meta$inner$y1 - panel_meta$inner$y0),
          label = lbl_vals[idx]
        )
      }
      invisible()
    },

    #' @description Draw text via `grid::grid.text()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      gp    <- .aes_to_gpar(self$aesthetics)
      angle <- self$aesthetics$angle %||% 0
      hjust <- self$aesthetics$hjust %||% 0.5
      vjust <- self$aesthetics$vjust %||% 0.5
      for (coords in self$coordCanvas) {
        if (is.null(coords) || length(coords$x) == 0) next
        grid::grid.text(
          label = coords$label,
          x     = grid::unit(coords$x, "npc"),
          y     = grid::unit(coords$y, "npc"),
          just  = c(hjust, vjust),
          rot   = angle,
          gp    = gp
        )
      }
    }
  )
)

#' Draw text labels
#'
#' Primitive element drawing a text label at each observation's `(x, y)`.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y`, `label`.
#' @param aesthetics Optional `SeqAes`. Supports `size`, `color`, `angle`,
#'   `hjust`, `vjust`, and a constant `label` when not supplied via `map()`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqTextR6` instance.
#' @examples
#' seq_text(map(x = start, y = score, label = name))
#' @export
seq_text <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqTextR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
