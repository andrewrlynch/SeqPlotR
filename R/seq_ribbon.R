# ── SeqRibbonR6 ───────────────────────────────────────────────────────────────
#
# Composite element: fills the band between y_min and y_max along x.

#' SeqRibbon R6 class
#'
#' Internal R6 generator backing [seq_ribbon()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqRibbonR6 <- R6::R6Class("SeqRibbon",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqRibbonR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `x`, `y_min`, `y_max`.
    #' @param aesthetics Optional `SeqAes`. Defaults to grey fill, alpha 0.8.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      if (is.null(aesthetics) || length(aesthetics) == 0) {
        aesthetics <- aes(fill = "grey60", alpha = 0.8)
      } else {
        if (is.null(aesthetics$fill))  aesthetics$fill  <- "grey60"
        if (is.null(aesthetics$alpha)) aesthetics$alpha <- 0.8
      }
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip to windows, transform to npc, and
    #'   build a closed polygon bounded above by `y_max` and below by `y_min`.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      n      <- length(eff_data)
      x_vals <- self$resolved$x     %||% BiocGenerics::start(eff_data)
      y_min  <- self$resolved$y_min %||% rep(0, n)
      y_max  <- self$resolved$y_max %||% rep(1, n)

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
        ylo <- y_min[idx]
        yhi <- y_max[idx]
        ord <- order(xs)
        xs  <- xs[ord]
        ylo <- ylo[ord]
        yhi <- yhi[ord]

        u    <- pmax(pmin((xs  - pm$xscale[1]) / diff(pm$xscale), 1), 0)
        v_lo <- pmax(pmin((ylo - pm$yscale[1]) / diff(pm$yscale), 1), 0)
        v_hi <- pmax(pmin((yhi - pm$yscale[1]) / diff(pm$yscale), 1), 0)

        xw <- pm$inner$x1 - pm$inner$x0
        yw <- pm$inner$y1 - pm$inner$y0
        x_c   <- pm$inner$x0 + u    * xw
        ylo_c <- pm$inner$y0 + v_lo * yw
        yhi_c <- pm$inner$y0 + v_hi * yw

        self$coordCanvas[[w]] <- list(
          x = c(x_c, rev(x_c)),
          y = c(yhi_c, rev(ylo_c))
        )
      }
      invisible()
    },

    #' @description Draw the ribbon polygon via `grid::grid.polygon()`.
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

#' Draw a filled ribbon between two y series
#'
#' Composite element filling the band between `y_min` and `y_max` at each
#' `x` value. Useful for confidence intervals and uncertainty bands.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required fields: `x`, `y_min`, `y_max`.
#' @param aesthetics Optional `SeqAes`. Supports `fill`, `color`, `alpha`,
#'   `linewidth`. Defaults to `aes(fill = "grey60", alpha = 0.8)`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqRibbonR6` instance.
#' @examples
#' seq_ribbon(map(x = start, y_min = lo, y_max = hi))
#' @export
seq_ribbon <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqRibbonR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
