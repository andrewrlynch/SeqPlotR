# ── SeqDensityR6 ──────────────────────────────────────────────────────────────
#
# Composite element: computes a 1D kernel density from the resolved y-values
# and renders it as a filled area. The density x-axis (the distribution of
# y-values) is drawn horizontally across the panel and mapped via `yscale`.

#' SeqDensity R6 class
#'
#' Internal R6 generator backing [seq_density()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqDensityR6 <- R6::R6Class("SeqDensity",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqDensityR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Required: `y`.
    #' @param aesthetics Optional `SeqAes`. `bw` controls the kernel
    #'   bandwidth (passed to [stats::density()]). Default fill `"grey60"`,
    #'   alpha 0.8.
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

    #' @description Compute `stats::density()` of the resolved y, then build a
    #'   closed polygon. Density evaluation points are mapped to canvas x via
    #'   the panel's `yscale`; heights are normalised to `[0, 1]` of panel
    #'   height.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      y_raw <- self$resolved$y
      if (is.null(y_raw)) return(invisible())
      y_raw <- y_raw[!is.na(y_raw)]
      if (length(y_raw) < 2) return(invisible())

      bw <- self$aesthetics$bw %||% "nrd0"
      dens <- stats::density(y_raw, bw = bw)

      self$coordCanvas <- vector("list", length(track_windows))
      for (w in seq_along(track_windows)) {
        pm <- layout_track[[w]]
        if (is.null(pm)) next

        u_d <- pmax(pmin((dens$x - pm$yscale[1]) / diff(pm$yscale), 1), 0)
        xw  <- pm$inner$x1 - pm$inner$x0
        yw  <- pm$inner$y1 - pm$inner$y0
        x_c <- pm$inner$x0 + u_d * xw

        dmax <- max(dens$y)
        v_d  <- if (dmax > 0) dens$y / dmax else rep(0, length(dens$y))
        y_c  <- pm$inner$y0 + v_d * yw

        self$coordCanvas[[w]] <- list(
          x = c(x_c, rev(x_c)),
          y = c(y_c, rep(pm$inner$y0, length(y_c)))
        )
      }
      invisible()
    },

    #' @description Draw the density polygon via `grid::grid.polygon()`.
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

#' Draw a kernel density estimate as a filled area
#'
#' Composite element that computes [stats::density()] of the resolved `y`
#' values and renders the distribution as a filled area. The density
#' evaluation axis (the distribution of y values) is drawn horizontally
#' across the panel using the track's `yscale`; densities are normalised
#' to the panel height.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Required field: `y`.
#' @param aesthetics Optional `SeqAes`. Supports `fill`, `color`, `alpha`,
#'   `linewidth`, and `bw` (kernel bandwidth, default `"nrd0"`).
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqDensityR6` instance.
#' @examples
#' seq_density(map(y = score))
#' @export
seq_density <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqDensityR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
