# TWO-DIMENSIONAL COLOUR KEY

#' Two-channel colour key specification
#'
#' Describes a smooth colour square over `[0, 1] x [0, 1]` for
#' [seq_tile_mix()], with channel `a` along x and channel `b` along y. Titles,
#' tick positions and tick labels are supplied by the caller.
#'
#' @param color_a,color_b Colours of the two channels.
#' @param x_title,y_title Axis titles, or `NULL`.
#' @param x_breaks,y_breaks Tick positions in `[0, 1]`.
#' @param x_labels,y_labels Tick labels. Default to the breaks.
#' @param note Optional text drawn below the key.
#' @param size Side length of the key in inches.
#' @param resolution Number of raster cells per side.
#' @param position One of `"canvas_margin"`, `"track_margin"`, `"inside"`.
#' @param side Margin side for margin positions.
#' @param x Horizontal anchor as a fraction of the target width.
#' @param y Vertical anchor of the key top as a fraction of the target height.
#'   In a canvas margin the key is placed below the other legend keys of the
#'   same side when they exist, and `y` is used otherwise.
#' @return A `SeqLegend2DSpec` object.
#' @examples
#' seq_legend_2d("#D62728", "#2B59C3", x_title = "A", y_title = "B")
#' @export
seq_legend_2d <- function(color_a, color_b, x_title = NULL, y_title = NULL,
                          x_breaks = c(0, 1), x_labels = NULL,
                          y_breaks = c(0, 1), y_labels = NULL,
                          note = NULL, size = 0.85, resolution = 128,
                          position = c("canvas_margin", "track_margin", "inside"),
                          side = "right", x = 0.08, y = 0.5) {
  position <- match.arg(position)
  side <- match.arg(side, c("top", "bottom", "left", "right"))
  for (cl in list(color_a, color_b))
    tryCatch(grDevices::col2rgb(cl),
             error = function(e) stop("seq_legend_2d: invalid colour '", cl, "'.",
                                      call. = FALSE))
  if (!is.numeric(size) || length(size) != 1L || size <= 0)
    stop("seq_legend_2d: size must be a positive number of inches.", call. = FALSE)
  x_labels <- x_labels %||% as.character(x_breaks)
  y_labels <- y_labels %||% as.character(y_breaks)
  if (length(x_labels) != length(x_breaks) || length(y_labels) != length(y_breaks))
    stop("seq_legend_2d: labels must match breaks in length.", call. = FALSE)
  structure(
    list(color_a = color_a, color_b = color_b, x_title = x_title,
         y_title = y_title, x_breaks = x_breaks, x_labels = x_labels,
         y_breaks = y_breaks, y_labels = y_labels, note = note, size = size,
         resolution = as.integer(resolution), position = position,
         side = side, x = x, y = y),
    class = "SeqLegend2DSpec"
  )
}

# Draw a 2D key into a target rect (named list x0, x1, y0, y1 in canvas npc).
# top is the key top in npc and overrides the anchor y when not NULL.
.draw_legend_2d <- function(spec, rect, top = NULL, cex = 0.65) {
  if (!inherits(spec, "SeqLegend2DSpec"))
    stop("spec must be a SeqLegend2DSpec.", call. = FALSE)
  if ((rect$x1 - rect$x0) <= 0 || (rect$y1 - rect$y0) <= 0)
    return(invisible())
  n <- spec$resolution
  u <- seq(0, 1, length.out = n)
  v <- seq(1, 0, length.out = n)
  img <- matrix(
    .mix_colors(rep(u, each = n), rep(v, times = n), spec$color_a, spec$color_b),
    nrow = n, ncol = n
  )
  gp_txt <- grid::gpar(cex = cex, col = "#1C1B1A")
  key_in <- grid::unit(spec$size, "inches")
  top_y  <- top %||% (rect$y0 + spec$y * (rect$y1 - rect$y0))

  # block left edge; the y title and tick labels sit left of the key
  blk <- grid::unit(rect$x0 + spec$x * (rect$x1 - rect$x0), "npc")
  ylab_w <- grid::grobWidth(grid::textGrob(spec$y_labels, gp = gp_txt))
  left <- blk + ylab_w + grid::unit(1.5, "mm") +
          (if (is.null(spec$y_title)) grid::unit(0, "mm") else grid::unit(4, "mm"))
  bottom <- grid::unit(top_y, "npc") - key_in

  grid::grid.raster(
    img, x = left, y = bottom, width = key_in, height = key_in,
    just = c("left", "bottom"), interpolate = TRUE,
    name = "seqplotr.legend_2d.raster"
  )
  grid::grid.rect(
    x = left, y = bottom, width = key_in, height = key_in,
    just = c("left", "bottom"),
    gp = grid::gpar(col = "grey40", fill = NA),
    name = "seqplotr.legend_2d.border"
  )
  grid::grid.text(
    spec$x_labels, x = left + grid::unit(spec$x_breaks * spec$size, "inches"),
    y = bottom - grid::unit(1, "mm"), just = c("center", "top"), gp = gp_txt,
    name = "seqplotr.legend_2d.xlabels"
  )
  grid::grid.text(
    spec$y_labels, x = left - grid::unit(1, "mm"),
    y = bottom + grid::unit(spec$y_breaks * spec$size, "inches"),
    just = c("right", "center"), gp = gp_txt,
    name = "seqplotr.legend_2d.ylabels"
  )
  if (!is.null(spec$x_title))
    grid::grid.text(
      spec$x_title, x = blk, y = bottom - grid::unit(4.5, "mm"),
      just = c("left", "top"), gp = gp_txt,
      name = "seqplotr.legend_2d.xtitle"
    )
  if (!is.null(spec$y_title))
    grid::grid.text(
      spec$y_title, x = blk, y = bottom + key_in * 0.5, rot = 90,
      just = c("center", "top"), gp = gp_txt,
      name = "seqplotr.legend_2d.ytitle"
    )
  if (!is.null(spec$note))
    grid::grid.text(
      spec$note, x = blk, y = bottom - grid::unit(8.5, "mm"),
      just = c("left", "top"), gp = grid::gpar(cex = cex, col = "grey35"),
      name = "seqplotr.legend_2d.note"
    )
  invisible()
}

# Collect canvas-margin 2D specs from the elements of a list of tracks.
.collect_canvas_2d_specs <- function(tracks) {
  out <- list()
  for (ti in seq_along(tracks)) {
    trk <- tracks[[ti]]
    if (!isTRUE(trk$show_legend)) next
    for (el in trk$elements) {
      if (!isTRUE(el$show_legend)) next
      spec <- el$legend
      if (!inherits(spec, "SeqLegend2DSpec")) next
      if (spec$position != "canvas_margin") next
      out <- c(out, list(list(spec = spec, track_idx = ti)))
    }
  }
  out
}
