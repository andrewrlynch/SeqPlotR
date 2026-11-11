# CANVAS MAPPING

#' Map data coordinates into a panel's canvas npc coordinates
#'
#' Maps genomic `x` and data `y` values into one panel of a track, the same
#' way the built-in elements do: `x` through the panel's x plot range and `y`
#' through its y plot range (falling back to the scale limits), into the
#' panel's `inner` rectangle. Use it in `prep()` of custom elements and links,
#' including links that place values on another track via
#' `layout_all_tracks[[track_id]]`.
#'
#' @param x Numeric vector of genomic (or x-axis) positions.
#' @param y Numeric vector of y values in the track's axis units.
#' @param panel One panel-metadata list, e.g. `layout_track[[w]]` in an
#'   element's `prep()` or `layout_all_tracks[[track_id]][[w]]` in a link's.
#' @param apply_transform Logical. Apply the panel's `y_transform` (set by
#'   transformed scales such as [seq_scale_count()]) to `y` first. Use
#'   `FALSE` for values that were already transformed by `resolve()`.
#' @param clamp Logical. Clamp the result to the panel's inner rectangle.
#' @return A list with numeric vectors `x` and `y` in canvas npc units.
#'   `NA` inputs stay `NA`.
#' @examples
#' panel <- list(xscale = c(0, 100), yscale = c(0, 1),
#'               inner = list(x0 = 0.1, x1 = 0.9, y0 = 0.2, y1 = 0.8))
#' seq_to_canvas(c(25, 50), c(0, 0.5), panel)
#' @export
seq_to_canvas <- function(x, y, panel, apply_transform = TRUE, clamp = TRUE) {
  if (isTRUE(apply_transform) && is.function(panel$y_transform))
    y <- panel$y_transform(y)
  xpr <- panel$xplot_range %||% panel$xscale
  ypr <- panel$yplot_range %||% panel$yscale
  u <- (as.numeric(x) - xpr[1]) / diff(xpr)
  v <- (as.numeric(y) - ypr[1]) / diff(ypr)
  if (isTRUE(clamp)) {
    u <- pmin(pmax(u, 0), 1)
    v <- pmin(pmax(v, 0), 1)
  }
  inner <- panel$inner
  list(x = inner$x0 + u * (inner$x1 - inner$x0),
       y = inner$y0 + v * (inner$y1 - inner$y0))
}

# SHARED DRAWING HELPERS

# per-element alpha on concrete colours (adjustcolor takes one alpha only)
.col_alpha <- function(col, alpha) {
  col <- as.character(col)
  col[is.na(col)] <- "grey70"
  m <- grDevices::col2rgb(col, alpha = TRUE)
  a <- rep_len(as.numeric(alpha), ncol(m))
  a[is.na(a)] <- 1
  grDevices::rgb(m[1, ], m[2, ], m[3, ], pmin(pmax(m[4, ] * a, 0), 255),
                 maxColorValue = 255)
}

# grob name seqplotr.<element>.<part>.<window>
.grob_name <- function(prefix, part, suffix) {
  if (is.null(prefix)) return(NULL)
  paste(c(prefix, part, if (length(suffix) && nzchar(suffix)) suffix),
        collapse = ".")
}

# Interval glyphs at canvas npc x/y/lo/hi; glyph in point, latest, earliest, interval, untimed.
.draw_interval_glyphs <- function(x, y, lo, hi, glyph, color, alpha = 1,
                                  size = 2, bar_lwd = 4, bar_alpha = 1,
                                  cap_mm = 2.5, points = TRUE,
                                  name = "seqplotr.interval_glyph",
                                  suffix = "") {
  n <- length(x)
  if (n == 0L) return(invisible())
  glyph <- rep_len(as.character(glyph), n)
  glyph[is.na(glyph)] <- "untimed"
  y <- rep_len(as.numeric(y), n)
  lo <- rep_len(as.numeric(lo), n)
  hi <- rep_len(as.numeric(hi), n)
  alpha <- rep_len(as.numeric(alpha), n)
  alpha[is.na(alpha)] <- 1
  col <- .col_alpha(rep_len(color, n), alpha)
  bar_col <- .col_alpha(rep_len(color, n), alpha * bar_alpha)
  ok_y <- is.finite(y)
  ok_b <- is.finite(lo) & is.finite(hi)
  half <- grid::unit(cap_mm / 2, "mm")

  # solid bars for point glyphs
  s <- which(glyph == "point" & ok_b)
  if (length(s))
    grid::grid.segments(grid::unit(x[s], "npc"), grid::unit(lo[s], "npc"),
                        grid::unit(x[s], "npc"), grid::unit(hi[s], "npc"),
                        name = .grob_name(name, "bar", suffix),
                        gp = grid::gpar(col = bar_col[s], lwd = bar_lwd,
                                        lineend = "butt"))

  # open bars: latest runs from y down to bar_lo, earliest from y up to bar_hi
  top <- ifelse(glyph == "latest", ifelse(ok_y, y, hi), hi)
  bot <- ifelse(glyph == "earliest", ifelse(ok_y, y, lo), lo)
  s <- which(glyph %in% c("latest", "earliest") & is.finite(top) &
               is.finite(bot))
  if (length(s))
    grid::grid.rect(x = grid::unit(x[s], "npc"),
                    y = grid::unit((top[s] + bot[s]) / 2, "npc"),
                    width = grid::unit(bar_lwd / 96, "inches"),
                    height = grid::unit(abs(top[s] - bot[s]), "npc"),
                    name = .grob_name(name, "open_bar", suffix),
                    gp = grid::gpar(col = bar_col[s], fill = NA, lwd = 1))

  # caps at y for one-sided bounds
  s <- which(glyph %in% c("latest", "earliest") & ok_y)
  if (length(s))
    grid::grid.segments(grid::unit(x[s], "npc") - half, grid::unit(y[s], "npc"),
                        grid::unit(x[s], "npc") + half, grid::unit(y[s], "npc"),
                        name = .grob_name(name, "cap", suffix),
                        gp = grid::gpar(col = col[s], lwd = 1.5,
                                        lineend = "butt"))

  # brackets: vertical line from bar_lo to bar_hi with end ticks
  s <- which(glyph == "interval" & ok_b)
  if (length(s)) {
    grid::grid.segments(grid::unit(x[s], "npc"), grid::unit(lo[s], "npc"),
                        grid::unit(x[s], "npc"), grid::unit(hi[s], "npc"),
                        name = .grob_name(name, "bracket", suffix),
                        gp = grid::gpar(col = col[s], lwd = 1.2))
    xs <- grid::unit(c(x[s], x[s]), "npc")
    ys <- grid::unit(c(lo[s], hi[s]), "npc")
    grid::grid.segments(xs - half, ys, xs + half, ys,
                        name = .grob_name(name, "bracket_tick", suffix),
                        gp = grid::gpar(col = c(col[s], col[s]), lwd = 1.2))
  }

  # points at y for point glyphs and for intervals with a y
  s <- which(glyph %in% c("point", "interval") & ok_y)
  if (isTRUE(points) && length(s))
    grid::grid.points(grid::unit(x[s], "npc"), grid::unit(y[s], "npc"),
                      pch = 21, size = grid::unit(size, "mm"),
                      name = .grob_name(name, "point", suffix),
                      gp = grid::gpar(col = col[s], fill = col[s]))
  invisible()
}

# Open caution triangle with "!"; base-left corner at canvas (x, y) shifted by dx_mm/dy_mm.
.draw_caution_marks <- function(x, y, mm = 1.6, dx_mm = 0.6, dy_mm = 0,
                                col = "#1C1B1A", name = NULL) {
  n <- length(x)
  if (n == 0L) return(invisible())
  px <- grid::unit(rep(x, each = 3), "npc") +
    grid::unit(dx_mm + rep(c(0, mm, mm / 2), n), "mm")
  py <- grid::unit(rep(y, each = 3), "npc") +
    grid::unit(dy_mm + rep(c(0, 0, mm * 0.9), n), "mm")
  grid::grid.polygon(px, py, id = rep(seq_len(n), each = 3),
                     name = if (!is.null(name)) paste0(name, ".triangle"),
                     gp = grid::gpar(col = col, fill = "white", lwd = 0.8))
  grid::grid.text("!", grid::unit(x, "npc") + grid::unit(dx_mm + mm / 2, "mm"),
                  grid::unit(y, "npc") + grid::unit(dy_mm + mm * 0.34, "mm"),
                  name = if (!is.null(name)) paste0(name, ".text"),
                  gp = grid::gpar(col = col, fontsize = mm * 3.6,
                                  fontface = "bold"))
  invisible()
}

# legend key drawing the caution mark centred in its cell
.caution_legend_key <- function(label = "caution", mm = 1.6,
                                col = "#1C1B1A") {
  LegendKey(label = label, color = col, shape = "custom",
            draw = function(x0, x1, y, height)
              .draw_caution_marks((x0 + x1) / 2, y, mm = mm, dx_mm = -mm / 2,
                                  dy_mm = -mm * 0.45, col = col))
}

# Legend for elements taking legend = NULL (auto keys) / FALSE / spec / key list.
.glyph_element_legend <- function(auto_keys, legend_arg, legend_labels = NULL,
                                  legend_extra = NULL,
                                  position = "canvas_margin", side = "right") {
  if (isFALSE(legend_arg)) return(NULL)
  if (inherits(legend_extra, "LegendKey")) legend_extra <- list(legend_extra)
  if (inherits(legend_arg, "SeqLegendSpec")) {
    if (length(legend_extra)) legend_arg$keys <- c(legend_arg$keys,
                                                   legend_extra)
    return(legend_arg)
  }
  keys <- if (inherits(legend_arg, "LegendKey")) list(legend_arg)
          else if (is.list(legend_arg) && length(legend_arg)) legend_arg
          else auto_keys
  if (length(legend_labels) && length(keys) && !is.null(names(keys))) {
    for (nm in intersect(names(legend_labels), names(keys)))
      keys[[nm]]$label <- unname(legend_labels[[nm]])
  }
  keys <- unname(c(keys, legend_extra))
  if (!length(keys)) return(NULL)
  is_margin <- position %in% c("track_margin", "canvas_margin")
  seq_legend(keys, position = position, side = if (is_margin) side,
             orientation = "vertical", x = 0.05, y = 0.95)
}
