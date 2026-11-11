# TWO-CHANNEL TILE ELEMENT

# Mix two channel intensities in [0, 1] into one colour per tile.
# Hue is interpolated by b / (a + b), depth from white by max(a, b).
.mix_colors <- function(a, b, col_a, col_b) {
  ca <- grDevices::col2rgb(col_a)[, 1] / 255
  cb <- grDevices::col2rgb(col_b)[, 1] / 255
  a  <- pmax(pmin(a, 1), 0)
  b  <- pmax(pmin(b, 1), 0)
  w  <- ifelse(a + b > 0, b / (a + b), 0)
  depth <- pmax(a, b)
  hue <- function(i) (1 - w) * ca[i] + w * cb[i]
  grDevices::rgb(
    1 - depth * (1 - hue(1)),
    1 - depth * (1 - hue(2)),
    1 - depth * (1 - hue(3))
  )
}

# Element drawing rectangles [x0, x1] x [y_min, y_max] coloured by two intensities.
SeqTileMixR6 <- R6::R6Class("SeqTileMix",
  inherit = SeqElementR6,
  public = list(
    y_fields   = c("y_min", "y_max"),
    legend_arg = NULL,

    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          legend = NULL, show_legend = TRUE, ...) {
      super$initialize(data, mapping, aesthetics, legend = NULL,
                       show_legend = show_legend)
      self$legend_arg <- legend
    },

    .infer_scale_y = function() {
      r <- self$resolved
      if (is.null(r$y_min) || is.null(r$y_max)) return(NULL)
      seq_scale_continuous(limits = range(c(r$y_min, r$y_max), na.rm = TRUE))
    },

    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      col_a <- self$aesthetics$color_a %||% "#D62728"
      col_b <- self$aesthetics$color_b %||% "#2B59C3"

      # legend: NULL gives the default 2D key, FALSE omits it
      la <- self$legend_arg
      self$legend <- if (isFALSE(la) || !isTRUE(self$show_legend)) NULL
                     else la %||% seq_legend_2d(col_a, col_b)

      self$coordCanvas <- vector("list", length(track_windows))
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())
      if (!inherits(eff_data, "GRanges"))
        stop("seq_tile_mix: data must be a GRanges.", call. = FALSE)

      r  <- self$resolved
      n  <- length(eff_data)
      a  <- as.numeric(r$a %||% rep(NA_real_, n))
      b  <- as.numeric(r$b %||% rep(NA_real_, n))
      x0 <- r$x0 %||% BiocGenerics::start(eff_data)
      x1 <- r$x1 %||% BiocGenerics::end(eff_data)
      ylo <- r$y_min
      yhi <- r$y_max
      if (is.null(ylo) || is.null(yhi))
        stop("seq_tile_mix: map(y_min = , y_max = ) are required.", call. = FALSE)

      # zero or missing tiles are not drawn
      keep <- !is.na(a) & !is.na(b) & pmax(a, b) > 0 &
              !is.na(ylo) & !is.na(yhi)
      fill <- rep(NA_character_, n)
      fill[keep] <- .mix_colors(a[keep], b[keep], col_a, col_b)

      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)
      for (w in unique(sh)) {
        idx <- qh[sh == w]
        idx <- idx[keep[idx]]
        if (length(idx) == 0L) next
        pm  <- layout_track[[w]]
        xpr <- pm$xplot_range %||% pm$xscale
        ypr <- pm$yplot_range %||% pm$yscale
        xw  <- pm$inner$x1 - pm$inner$x0
        yw  <- pm$inner$y1 - pm$inner$y0
        u0 <- pmax(pmin((pmin(x0[idx], x1[idx]) - xpr[1]) / diff(xpr), 1), 0)
        u1 <- pmax(pmin((pmax(x0[idx], x1[idx]) - xpr[1]) / diff(xpr), 1), 0)
        v0 <- pmax(pmin((pmin(ylo[idx], yhi[idx]) - ypr[1]) / diff(ypr), 1), 0)
        v1 <- pmax(pmin((pmax(ylo[idx], yhi[idx]) - ypr[1]) / diff(ypr), 1), 0)
        self$coordCanvas[[w]] <- data.frame(
          row  = idx,
          x0   = pm$inner$x0 + u0 * xw,
          x1   = pm$inner$x0 + u1 * xw,
          y0   = pm$inner$y0 + v0 * yw,
          y1   = pm$inner$y0 + v1 * yw,
          fill = fill[idx],
          stringsAsFactors = FALSE
        )
      }
      invisible()
    },

    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      seam <- self$aesthetics$seam %||% 0.25
      for (w in seq_along(self$coordCanvas)) {
        d <- self$coordCanvas[[w]]
        if (is.null(d) || nrow(d) == 0L) next
        grid::grid.rect(
          x      = grid::unit(d$x0, "npc"),
          y      = grid::unit(d$y0, "npc"),
          width  = grid::unit(d$x1 - d$x0, "npc"),
          height = grid::unit(d$y1 - d$y0, "npc"),
          just   = c("left", "bottom"),
          gp     = grid::gpar(fill = d$fill, col = d$fill, lwd = seam),
          name   = sprintf("seqplotr.tile_mix.%d", w)
        )
      }
      invisible()
    }
  )
)

#' Draw two-channel colour tiles
#'
#' Draws one rectangle per row over `[x0, x1] x [y_min, y_max]`. Two
#' intensities in `[0, 1]`, `a` and `b`, set the colour: the hue is
#' interpolated between `color_a` and `color_b` by `b / (a + b)` and the depth
#' from white is `max(a, b)`. Rows with both intensities zero or missing are
#' not drawn. Tiles carry a hairline border in their own colour so adjacent
#' tiles show no seams. Draw bound lines with [seq_segment()].
#'
#' @param data A `GRanges` with one row per tile.
#' @param mapping A `SeqMap`. Fields: `x0` (default `start`), `x1` (default
#'   `end`), `y_min`, `y_max`, `a`, `b`.
#' @param aesthetics A `SeqAes`. Recognised: `color_a` (default `"#D62728"`),
#'   `color_b` (default `"#2B59C3"`), `seam` (border line width, default
#'   `0.25`).
#' @param legend `NULL` for the default 2D key from [seq_legend_2d()], `FALSE`
#'   for no key, or a `SeqLegend2DSpec` used as given.
#' @param show_legend Logical. `FALSE` suppresses the key.
#' @param ... Unused.
#' @return A `SeqTileMixR6` instance.
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 101), c(100, 200)),
#'   lo = c(0, 0.2), hi = c(0.2, 0.4), g = c(1, 0.5), l = c(0, 0.5))
#' seq_tile_mix(gr, map(y_min = lo, y_max = hi, a = g, b = l))
#' @export
seq_tile_mix <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                         legend = NULL, show_legend = TRUE, ...) {
  SeqTileMixR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                   legend = legend, show_legend = show_legend, ...)
}
