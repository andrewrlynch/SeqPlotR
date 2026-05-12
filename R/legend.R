# ── SeqLegendSpec ─────────────────────────────────────────────────────────────
#
# Lightweight S3 placement-spec for a legend group. Records all layout intent
# (position, side, orientation, anchor) without performing any rendering.
# Consumed by rendering agents in later batches.

#' Create a legend placement specification
#'
#' @description
#' Constructs a `SeqLegendSpec` object that describes where and how a legend
#' group should be rendered. The spec captures placement intent only — no
#' rendering logic is executed here. It is consumed by the legend-rendering
#' layer in a later batch.
#'
#' @param keys A single `LegendKey`, or a named or unnamed list of `LegendKey`
#'   objects. Required.
#' @param title Optional character string. The legend group title drawn above
#'   or beside the key block.
#' @param position Character. One of `"inside"`, `"track_margin"`, or
#'   `"canvas_margin"`. Controls which area the legend occupies.
#'   Default `"inside"`.
#' @param x Numeric in \[0, 1\]. Horizontal position of the legend anchor
#'   within the target area. Default `0.5`.
#' @param y Numeric in \[0, 1\]. Vertical position of the legend anchor
#'   within the target area. Default `0.5`.
#' @param hjust Numeric in \[0, 1\]. Horizontal justification of the legend
#'   content block relative to the anchor. `0` = left-aligned, `1` =
#'   right-aligned. Default `0`.
#' @param orientation Character. One of `"horizontal"` or `"vertical"`.
#'   Controls whether keys are laid out as a row or a column. When `NULL`
#'   (default), orientation is inferred from `side`: `"vertical"` for
#'   `side %in% c("left", "right")`, `"horizontal"` otherwise.
#' @param nrow Integer or `NULL`. Number of rows in the key grid. When both
#'   `nrow` and `ncol` are `NULL`, defaults to 1 row (all keys in one row for
#'   horizontal orientation) or 1 column (all keys in one column for vertical
#'   orientation).
#' @param ncol Integer or `NULL`. Number of columns in the key grid.
#' @param side Character or `NULL`. For `position %in% c("track_margin",
#'   "canvas_margin")`, which margin to target. One of `"top"`, `"bottom"`,
#'   `"left"`, `"right"`. When `NULL`, defaults to `"top"` for margin
#'   positions and is left `NULL` for `"inside"`.
#'
#' @return A `SeqLegendSpec` S3 object (a named list with class
#'   `"SeqLegendSpec"`).
#'
#' @examples
#' k1 <- LegendKey(label = "H3K27ac", color = "firebrick")
#' k2 <- LegendKey(label = "H3K4me3", color = "steelblue")
#'
#' # Inside the track, top-left
#' spec <- seq_legend(list(k1, k2), title = "Marks", x = 0.02, y = 0.95)
#'
#' # In the outer track margin, bottom
#' spec <- seq_legend(list(k1, k2), position = "track_margin", side = "bottom")
#'
#' @export
seq_legend <- function(keys,
                       title       = NULL,
                       position    = "inside",
                       x           = 0.5,
                       y           = 0.5,
                       hjust       = 0,
                       orientation = NULL,
                       nrow        = NULL,
                       ncol        = NULL,
                       side        = NULL) {

  # --- Validate keys ---
  if (inherits(keys, "LegendKey")) {
    keys <- list(keys)
  }
  if (!is.list(keys) || length(keys) == 0L) {
    stop("`keys` must be a LegendKey or a non-empty list of LegendKey objects.")
  }
  for (i in seq_along(keys)) {
    if (!inherits(keys[[i]], "LegendKey")) {
      stop(sprintf("`keys[[%d]]` must be a LegendKey object.", i))
    }
  }

  # --- Validate position ---
  valid_positions <- c("inside", "track_margin", "canvas_margin")
  if (!position %in% valid_positions) {
    stop(sprintf(
      "`position` must be one of %s.",
      paste(sprintf('"%s"', valid_positions), collapse = ", ")
    ))
  }

  # --- Validate x / y ---
  if (!is.numeric(x) || length(x) != 1L || x < 0 || x > 1) {
    stop("`x` must be a single numeric value in [0, 1].")
  }
  if (!is.numeric(y) || length(y) != 1L || y < 0 || y > 1) {
    stop("`y` must be a single numeric value in [0, 1].")
  }

  # --- Validate hjust ---
  if (!is.numeric(hjust) || length(hjust) != 1L || hjust < 0 || hjust > 1) {
    stop("`hjust` must be a single numeric value in [0, 1].")
  }

  # --- Validate orientation ---
  if (!is.null(orientation) && !orientation %in% c("horizontal", "vertical")) {
    stop('`orientation` must be "horizontal", "vertical", or NULL.')
  }

  # --- Validate side ---
  valid_sides <- c("top", "bottom", "left", "right")
  if (!is.null(side) && !side %in% valid_sides) {
    stop(sprintf(
      "`side` must be one of %s, or NULL.",
      paste(sprintf('"%s"', valid_sides), collapse = ", ")
    ))
  }

  # --- Infer orientation when NULL ---
  resolved_orientation <- orientation
  if (is.null(resolved_orientation)) {
    if (!is.null(side) && side %in% c("left", "right")) {
      resolved_orientation <- "vertical"
    } else {
      resolved_orientation <- "horizontal"
    }
  }

  # --- Infer side default for margin positions ---
  resolved_side <- side
  if (is.null(resolved_side) && position %in% c("track_margin", "canvas_margin")) {
    resolved_side <- "top"
  }

  spec <- list(
    keys        = keys,
    title       = title,
    position    = position,
    x           = x,
    y           = y,
    hjust       = hjust,
    orientation = resolved_orientation,
    nrow        = nrow,
    ncol        = ncol,
    side        = resolved_side
  )
  class(spec) <- "SeqLegendSpec"
  spec
}


#' Test if an object is a SeqLegendSpec
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_seq_legend_spec <- function(x) inherits(x, "SeqLegendSpec")


#' @export
print.SeqLegendSpec <- function(x, ...) {
  cat(sprintf(
    "<SeqLegendSpec> position=%s  side=%s  orientation=%s  keys=%d%s\n",
    x$position,
    if (is.null(x$side)) "NULL" else x$side,
    x$orientation,
    length(x$keys),
    if (!is.null(x$title)) sprintf('  title="%s"', x$title) else ""
  ))
  invisible(x)
}


# ── Internal rendering helpers ─────────────────────────────────────────────────

# Draw a single legend glyph for one LegendKey at a given position.
#
# Arguments:
#   key    A LegendKey object.
#   x0     Left edge of the key symbol area (NPC).
#   x1     Right edge of the key symbol area (NPC).
#   y      Vertical centre of the key row (NPC).
#   height Height to allocate for the symbol (NPC).
.SeqLegend_drawKey <- function(key, x0, x1, y, height) {
  col  <- key$color %||% "#1C1B1A"
  fill <- key$fill  %||% col
  alp  <- key$alpha %||% 1
  lty  <- key$lty   %||% 1
  shp  <- key$shape %||% "-"

  if (identical(shp, "-") || identical(shp, "line")) {
    grid::grid.lines(
      x  = grid::unit(c(x0, x1), "npc"),
      y  = grid::unit(c(y,  y),  "npc"),
      gp = grid::gpar(col = col, alpha = alp, lty = lty, lwd = 1.5)
    )
  } else if (identical(shp, "rect")) {
    grid::grid.rect(
      x      = grid::unit((x0 + x1) / 2, "npc"),
      y      = grid::unit(y, "npc"),
      width  = grid::unit(x1 - x0, "npc"),
      height = grid::unit(max(height * 0.7, 0.01), "npc"),
      gp     = grid::gpar(col = col, fill = fill, alpha = alp)
    )
  } else {
    pch <- switch(as.character(shp),
      "circle"   = 21L,
      "square"   = 22L,
      "triangle" = 24L,
      "diamond"  = 23L,
      "point"    = 21L,
      21L
    )
    grid::grid.points(
      x    = grid::unit((x0 + x1) / 2, "npc"),
      y    = grid::unit(y, "npc"),
      pch  = pch,
      size = grid::unit(max(height * 0.8, 0.01), "npc"),
      gp   = grid::gpar(col = col, fill = fill, alpha = alp)
    )
  }
}


# Compute per-key cell geometry from a SeqLegendSpec and a bounding rect.
#
# Arguments:
#   spec    SeqLegendSpec produced by seq_legend().
#   bbox    Named list: x0, x1, y0, y1  (NPC coords of the target area).
#   pad     Fractional padding inside bbox (applied uniformly). Default 0.04.
#   key_w   Width (NPC) allocated for the key symbol. Default 0.04.
#   key_gap Gap between key symbol right edge and label left edge. Default 0.006.
#   row_h   Row height as a fraction of bbox height. Default 0.16.
#
# Returns a list with:
#   $title_x, $title_y  NPC coords for the group title (NULL when spec$title is NULL).
#   $cells              List of per-key lists:
#                         key, label, key_x0, key_x1, y, text_x, text_y.
#   $bbox               The bbox passed in.
.legend_layout_cells <- function(spec,
                                 bbox,
                                 pad     = 0.04,
                                 key_w   = 0.04,
                                 key_gap = 0.006,
                                 row_h   = 0.16) {
  keys    <- spec$keys
  n_items <- length(keys)
  if (n_items == 0L) return(NULL)

  bw <- bbox$x1 - bbox$x0
  bh <- bbox$y1 - bbox$y0

  # Resolve nrow / ncol from spec or orientation default
  ncol <- spec$ncol
  nrow <- spec$nrow

  if (is.null(ncol) && is.null(nrow)) {
    if (spec$orientation == "vertical") {
      ncol <- 1L; nrow <- n_items
    } else {
      nrow <- 1L; ncol <- n_items
    }
  } else if (is.null(ncol)) {
    ncol <- max(1L, ceiling(n_items / nrow))
  } else if (is.null(nrow)) {
    nrow <- max(1L, ceiling(n_items / ncol))
  }
  ncol <- max(1L, as.integer(ncol))
  nrow <- max(1L, as.integer(nrow))

  has_title <- !is.null(spec$title)

  # Fixed row height: row_h fraction of bbox height, capped so keys stay visible.
  # This prevents cells from consuming the entire panel when bh is large.
  cell_h <- max(min(row_h * bh, 0.07), 0.008)

  # Usable area inside padding
  inner_x0 <- bbox$x0 + pad * bw
  inner_x1 <- bbox$x1 - pad * bw
  inner_w   <- inner_x1 - inner_x0
  cell_w    <- inner_w / ncol

  # Horizontal: x_block_left derived from spec$x (fraction of inner width)
  x_block_left <- inner_x0 + spec$x * inner_w - (ncol * cell_w) * spec$hjust

  # Vertical: spec$y is the TOP of the entire legend block (title included).
  # The legend extends downward from that anchor point.
  y_legend_top <- bbox$y0 + spec$y * bh

  # Title sits in the first row above the key rows
  title_x <- NULL; title_y <- NULL
  if (has_title) {
    title_x     <- x_block_left
    title_y     <- y_legend_top - cell_h / 2
    y_block_top <- y_legend_top - cell_h
  } else {
    y_block_top <- y_legend_top
  }

  cells <- vector("list", n_items)
  for (i in seq_len(n_items)) {
    r  <- ((i - 1L) %/% ncol) + 1L
    cl <- ((i - 1L) %%  ncol) + 1L

    x0 <- x_block_left + (cl - 1L) * cell_w
    x1 <- x0 + cell_w
    y1 <- y_block_top  - (r  - 1L) * cell_h
    y0 <- y1 - cell_h

    key <- keys[[i]]
    xk1 <- min(x1, x0 + key_w)

    cells[[i]] <- list(
      key    = key,
      label  = key$label,
      key_x0 = x0,
      key_x1 = xk1,
      y      = (y0 + y1) / 2,
      text_x = min(inner_x1, xk1 + key_gap),
      text_y = (y0 + y1) / 2
    )
  }

  list(title_x = title_x, title_y = title_y, cells = cells, bbox = bbox)
}


# Draw a SeqLegendSpec with position="inside" into the data panel of one window.
#
# Arguments:
#   spec        SeqLegendSpec (position must equal "inside").
#   panel_meta  One entry from layout$panelBounds[[track_key]][[window_idx]];
#               must have $full (list: x0,x1,y0,y1) and optionally $inner.
#   key_cex     Character size for key labels. Default 0.85.
#   title_cex   Character size for the group title. Default 0.90.
#   title_col   Colour for title text. Default "#1C1B1A".
#   label_col   Colour for key label text. Default "#1C1B1A".
.draw_legend_inside <- function(spec,
                                panel_meta,
                                key_cex   = 0.85,
                                title_cex = 0.90,
                                title_col = "#1C1B1A",
                                label_col = "#1C1B1A") {
  if (!inherits(spec, "SeqLegendSpec")) stop("spec must be a SeqLegendSpec.")
  if (spec$position != "inside") return(invisible())

  # Use inner panel area (plot area after margins), falling back to full.
  p    <- if (!is.null(panel_meta$inner)) panel_meta$inner else panel_meta$full
  bbox <- list(x0 = p$x0, x1 = p$x1, y0 = p$y0, y1 = p$y1)

  lay <- .legend_layout_cells(spec, bbox)
  if (is.null(lay)) return(invisible())

  # Title
  if (!is.null(lay$title_x)) {
    grid::grid.text(
      label = spec$title,
      x     = grid::unit(lay$title_x, "npc"),
      y     = grid::unit(lay$title_y, "npc"),
      just  = c("left", "center"),
      gp    = grid::gpar(col = title_col, cex = title_cex, fontface = "bold")
    )
  }

  # Keys and labels
  key_height <- min(0.05, (p$y1 - p$y0) * 0.08)
  for (cell in lay$cells) {
    .SeqLegend_drawKey(
      key    = cell$key,
      x0     = cell$key_x0,
      x1     = cell$key_x1,
      y      = cell$y,
      height = key_height
    )
    if (!is.null(cell$label)) {
      grid::grid.text(
        label = cell$label,
        x     = grid::unit(cell$text_x, "npc"),
        y     = grid::unit(cell$text_y, "npc"),
        just  = c("left", "center"),
        gp    = grid::gpar(col = label_col, cex = key_cex)
      )
    }
  }

  invisible()
}


# ---- .draw_legend_track_margin -----------------------------------------------
#
# Draw a SeqLegendSpec with position="track_margin" into the outer margin of
# the specified track.
#
# Arguments:
#   spec              SeqLegendSpec (position must equal "track_margin")
#   track_margin_rect Named list x0, x1, y0, y1 — the margin rect for the
#                     chosen side, as produced by trackMarginBounds[[t]][[side]].
#   key_cex           Label character size. Default 0.85.
#   title_cex         Title character size. Default 0.90.
#   title_col         Title colour. Default "#1C1B1A".
#   label_col         Label colour. Default "#1C1B1A".
.draw_legend_track_margin <- function(spec,
                                      track_margin_rect,
                                      key_cex   = 0.85,
                                      title_cex = 0.90,
                                      title_col = "#1C1B1A",
                                      label_col = "#1C1B1A") {

  if (!inherits(spec, "SeqLegendSpec"))  stop("spec must be a SeqLegendSpec.")
  if (spec$position != "track_margin")   return(invisible())

  r <- track_margin_rect
  if ((r$x1 - r$x0) <= 0 || (r$y1 - r$y0) <= 0) return(invisible())

  layout <- .legend_layout_cells(spec, r)
  if (is.null(layout)) return(invisible())

  key_height <- min(0.04, (r$y1 - r$y0) * 0.35)

  if (!is.null(layout$title_x)) {
    grid::grid.text(
      label = spec$title,
      x     = grid::unit(layout$title_x, "npc"),
      y     = grid::unit(layout$title_y, "npc"),
      just  = c("left", "center"),
      gp    = grid::gpar(col = title_col, cex = title_cex, fontface = "bold")
    )
  }

  for (cell in layout$cells) {
    .SeqLegend_drawKey(
      key    = cell$key,
      x0     = cell$key_x0,
      x1     = cell$key_x1,
      y      = cell$y,
      height = key_height
    )
    grid::grid.text(
      label = cell$label,
      x     = grid::unit(cell$text_x, "npc"),
      y     = grid::unit(cell$text_y, "npc"),
      just  = c("left", "center"),
      gp    = grid::gpar(col = label_col, cex = key_cex)
    )
  }

  invisible()
}


# ---- .collect_canvas_legend_specs --------------------------------------------
#
# Walk all tracks and collect SeqLegendSpec objects with
# position = "canvas_margin".
#
# Arguments:
#   tracks   Flat list of SeqTrack objects (from sp$allTracks()).
#
# Returns a list of lists, each with:
#   $spec        The SeqLegendSpec
#   $track_idx   Integer, 1-based track index
.collect_canvas_legend_specs <- function(tracks) {
  out <- list()
  for (ti in seq_along(tracks)) {
    trk <- tracks[[ti]]
    if (!isTRUE(trk$show_legend)) next
    for (el in trk$elements) {
      if (!isTRUE(el$show_legend)) next
      if (is.null(el$legend))      next
      spec <- el$legend
      if (inherits(spec, "LegendKey") ||
          (is.list(spec) && !inherits(spec, "SeqLegendSpec"))) next
      if (!inherits(spec, "SeqLegendSpec"))  next
      if (spec$position != "canvas_margin")  next
      out <- c(out, list(list(spec = spec, track_idx = ti)))
    }
  }
  out
}


# ---- .merge_canvas_specs -----------------------------------------------------
#
# Merge a list of SeqLegendSpec entries (same side) into a single flat key
# list. Title changes between consecutive entries are represented by inserting
# a synthetic title key (shape = "none") as a visual group header.
#
# Arguments:
#   entries   List produced by .collect_canvas_legend_specs(), filtered to
#             a single side.
#   side      Character; the target side string.
#
# Returns a SeqLegendSpec or NULL when entries is empty.
.merge_canvas_specs <- function(entries, side) {
  if (length(entries) == 0L) return(NULL)

  merged_keys <- list()
  last_title  <- NA_character_

  for (entry in entries) {
    sp  <- entry$spec
    ttl <- if (is.null(sp$title)) NA_character_ else sp$title

    if (!identical(ttl, last_title) && !is.na(ttl)) {
      merged_keys <- c(merged_keys, list(LegendKey(label = ttl, shape = "none")))
      last_title  <- ttl
    }

    merged_keys <- c(merged_keys, sp$keys)
  }

  if (length(merged_keys) == 0L) return(NULL)

  first <- entries[[1L]]$spec
  seq_legend(
    keys        = merged_keys,
    title       = NULL,
    position    = "canvas_margin",
    x           = first$x,
    y           = first$y,
    hjust       = first$hjust,
    orientation = first$orientation,
    nrow        = first$nrow,
    ncol        = first$ncol,
    side        = side
  )
}


# ---- .draw_legend_canvas_margin ----------------------------------------------
#
# Draw a merged SeqLegendSpec into a canvas margin rect.
#
# Arguments:
#   spec              SeqLegendSpec (position = "canvas_margin")
#   canvas_margin_rect Named list x0, x1, y0, y1
#   key_cex, title_cex, title_col, label_col  — same as other draw helpers
.draw_legend_canvas_margin <- function(spec,
                                       canvas_margin_rect,
                                       key_cex   = 0.85,
                                       title_cex = 0.90,
                                       title_col = "#1C1B1A",
                                       label_col = "#1C1B1A") {

  if (!inherits(spec, "SeqLegendSpec"))  stop("spec must be a SeqLegendSpec.")
  if (spec$position != "canvas_margin")  return(invisible())

  r <- canvas_margin_rect
  if ((r$x1 - r$x0) <= 0 || (r$y1 - r$y0) <= 0) return(invisible())

  layout <- .legend_layout_cells(spec, r)
  if (is.null(layout)) return(invisible())

  key_height <- min(0.04, (r$y1 - r$y0) * 0.4)

  for (cell in layout$cells) {
    if (identical(cell$key$shape, "none")) {
      grid::grid.text(
        label = cell$label,
        x     = grid::unit(cell$key_x0, "npc"),
        y     = grid::unit(cell$y, "npc"),
        just  = c("left", "center"),
        gp    = grid::gpar(col = title_col, cex = title_cex, fontface = "bold")
      )
      next
    }

    .SeqLegend_drawKey(
      key    = cell$key,
      x0     = cell$key_x0,
      x1     = cell$key_x1,
      y      = cell$y,
      height = key_height
    )
    grid::grid.text(
      label = cell$label,
      x     = grid::unit(cell$text_x, "npc"),
      y     = grid::unit(cell$text_y, "npc"),
      just  = c("left", "center"),
      gp    = grid::gpar(col = label_col, cex = key_cex)
    )
  }

  invisible()
}


# ── GradientLegendSpec ─────────────────────────────────────────────────────────

#' Create a color gradient (color-bar) legend specification
#'
#' @description
#' Constructs a `GradientLegendSpec` describing a continuous color-scale legend
#' (color bar). When `breaks` is supplied the bar is rendered as discrete
#' `LegendKey` entries instead — useful for heatmaps that need only a handful
#' of labelled stops.
#'
#' @param palette One of `"viridis"`, `"plasma"`, `"magma"`, `"blues"`,
#'   `"reds"`. Default `"viridis"`.
#' @param limits Numeric vector of length 2. The data range the gradient spans.
#'   Default `c(0, 1)`.
#' @param title Optional character string. Legend group title.
#' @param position One of `"inside"`, `"track_margin"`, `"canvas_margin"`.
#' @param x,y Numeric in \[0, 1\]. Anchor of the bar within the target area.
#'   Default `0.5`.
#' @param hjust Horizontal justification of the bar relative to the anchor.
#'   Default `0`.
#' @param orientation `"horizontal"` or `"vertical"`. Bar direction.
#'   Inferred from `side` when `NULL`.
#' @param side For margin positions: one of `"top"`, `"bottom"`, `"left"`,
#'   `"right"`. Defaults to `"top"` for margin positions.
#' @param breaks `NULL` (continuous color bar), a single positive integer
#'   (that many evenly-spaced discrete keys), or a numeric vector of explicit
#'   break values.
#'
#' @return A `GradientLegendSpec` S3 object.
#' @examples
#' # Continuous color bar
#' seq_gradient_legend(palette = "viridis", limits = c(0, 100), title = "Score")
#'
#' # Five discrete keys
#' seq_gradient_legend(palette = "plasma", limits = c(-2, 2), breaks = 5)
#' @export
seq_gradient_legend <- function(palette     = "viridis",
                                 limits      = c(0, 1),
                                 title       = NULL,
                                 position    = "inside",
                                 x           = 0.5,
                                 y           = 0.5,
                                 hjust       = 0,
                                 orientation = NULL,
                                 side        = NULL,
                                 breaks      = NULL) {

  valid_positions <- c("inside", "track_margin", "canvas_margin")
  if (!position %in% valid_positions)
    stop(sprintf('`position` must be one of %s.',
                 paste(sprintf('"%s"', valid_positions), collapse = ", ")))

  valid_palettes <- c("viridis", "plasma", "magma", "blues", "reds")
  if (!palette %in% valid_palettes)
    stop(sprintf('`palette` must be one of %s.',
                 paste(sprintf('"%s"', valid_palettes), collapse = ", ")))

  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)))
    stop("`limits` must be a finite numeric vector of length 2.")

  if (!is.null(orientation) && !orientation %in% c("horizontal", "vertical"))
    stop('`orientation` must be "horizontal", "vertical", or NULL.')

  if (is.null(orientation)) {
    if (!is.null(side) && side %in% c("left", "right")) {
      orientation <- "vertical"
    } else {
      orientation <- "horizontal"
    }
  }
  if (is.null(side) && position %in% c("track_margin", "canvas_margin"))
    side <- "top"

  spec <- list(
    palette     = palette,
    limits      = limits,
    title       = title,
    position    = position,
    x           = x,
    y           = y,
    hjust       = hjust,
    orientation = orientation,
    side        = side,
    breaks      = breaks
  )
  class(spec) <- "GradientLegendSpec"
  spec
}

#' @export
print.GradientLegendSpec <- function(x, ...) {
  brk_str <- if (is.null(x$breaks)) "NULL (continuous)"
             else if (length(x$breaks) == 1L) sprintf("n=%d discrete", as.integer(x$breaks))
             else sprintf("c(%s)", paste(formatC(x$breaks, digits = 2L, format = "g"),
                                         collapse = ", "))
  cat(sprintf(
    "<GradientLegendSpec> palette=%s  limits=[%.3g, %.3g]  position=%s  side=%s  breaks=%s%s\n",
    x$palette, x$limits[1], x$limits[2], x$position,
    if (is.null(x$side)) "NULL" else x$side,
    brk_str,
    if (!is.null(x$title)) sprintf('  title="%s"', x$title) else ""
  ))
  invisible(x)
}


# ── Element → legend key shape lookup ─────────────────────────────────────────

# Named vector: R6 class name → legend key shape string
.element_key_shape <- c(
  SeqPoint   = "point",
  SeqBar     = "rect",
  SeqLine    = "line",
  SeqArea    = "rect",
  SeqSegment = "line",
  SeqTile    = "rect",
  SeqGene    = "rect",
  SeqArch    = "line",
  SeqRecon   = "line"
)

# Return the default key shape for an element, falling back to "rect".
.element_shape_for <- function(element) {
  cls <- class(element)[1L]
  val <- .element_key_shape[cls]
  if (is.na(val)) "rect" else unname(val)
}


# Format a number for axis tick labels (compact g-format).
.axis_format_num <- function(x) {
  formatC(x, digits = 3L, format = "g")
}


# Return the stop colors for a named palette.
.palette_stops <- function(palette_nm) {
  switch(palette_nm %||% "viridis",
    viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
    plasma  = c("#0d0887", "#cc4778", "#f0f921"),
    magma   = c("#000004", "#b63679", "#fcfdbf"),
    blues   = c("#f7fbff", "#2171b5", "#08306b"),
    reds    = c("#fff5f0", "#ef3b2c", "#67000d"),
    c("#440154", "#31688e", "#35b779", "#fde725")  # viridis fallback
  )
}


# Map t_vals in [0,1] to hex colors from the named palette.
.palette_color_at <- function(t_vals, palette_nm) {
  ramp <- grDevices::colorRamp(.palette_stops(palette_nm))
  t_c  <- pmax(0, pmin(1, t_vals))
  cols <- ramp(t_c)
  grDevices::rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)
}


# Convert a GradientLegendSpec with `breaks` to a list of LegendKey objects.
# Returns NULL when spec$breaks is NULL.
.gradient_spec_to_keys <- function(spec) {
  brks <- spec$breaks
  if (is.null(brks)) return(NULL)

  if (length(brks) == 1L) {
    n    <- max(1L, as.integer(brks))
    brks <- seq(spec$limits[1], spec$limits[2], length.out = n + 2L)[-c(1L, n + 2L)]
  } else {
    brks <- as.numeric(brks)
  }
  if (length(brks) == 0L) return(NULL)

  t_vals <- pmax(0, pmin(1, (brks - spec$limits[1]) / diff(spec$limits)))
  colors <- .palette_color_at(t_vals, spec$palette)
  labels <- formatC(brks, digits = 3L, format = "g")
  Map(function(lbl, col) LegendKey(label = lbl, color = col, fill = col),
      labels, as.list(colors))
}


# Draw a GradientLegendSpec rendered as discrete LegendKey rows.
.draw_gradient_as_keys <- function(spec, panel_meta, ...) {
  keys <- .gradient_spec_to_keys(spec)
  if (is.null(keys) || length(keys) == 0L) return(invisible())
  fake_leg <- seq_legend(
    keys        = keys,
    title       = spec$title,
    position    = "inside",
    x           = spec$x,
    y           = spec$y,
    hjust       = spec$hjust,
    orientation = spec$orientation %||% "vertical",
    side        = NULL
  )
  .draw_legend_inside(fake_leg, panel_meta, ...)
}


# Draw a GradientLegendSpec as a continuous color bar.
.draw_gradient_colorbar <- function(spec, panel_meta,
                                     key_cex   = 0.75,
                                     title_cex = 0.85,
                                     title_col = "#1C1B1A",
                                     label_col = "#1C1B1A") {

  p  <- if (!is.null(panel_meta$inner)) panel_meta$inner else panel_meta$full
  bw <- p$x1 - p$x0
  bh <- p$y1 - p$y0
  if (bw <= 0 || bh <= 0) return(invisible())

  is_horiz <- identical(spec$orientation, "horizontal")

  if (is_horiz) {
    bar_w_npc <- min(0.22, bw * 0.30)
    bar_h_npc <- min(0.018, bh * 0.05)
  } else {
    bar_w_npc <- min(0.018, bw * 0.05)
    bar_h_npc <- min(0.22, bh * 0.30)
  }
  gap <- min(0.012, bh * 0.025)

  x_anchor <- p$x0 + spec$x * bw - bar_w_npc * (spec$hjust %||% 0)
  y_anchor <- p$y0 + spec$y * bh   # top of legend block

  title_h <- 0
  if (!is.null(spec$title)) {
    grid::grid.text(
      label = spec$title,
      x     = grid::unit(x_anchor, "npc"),
      y     = grid::unit(y_anchor, "npc"),
      just  = c("left", "top"),
      gp    = grid::gpar(col = title_col, cex = title_cex, fontface = "bold")
    )
    # Reserve enough vertical space for the title text so the bar starts
    # below it. Use convertHeight() for the actual rendered line height.
    line_h  <- as.numeric(grid::convertHeight(
      grid::unit(title_cex * 1.2, "lines"), "npc", valueOnly = TRUE
    ))
    title_h <- line_h + gap
  }

  bar_y_top <- y_anchor - title_h
  bar_cx    <- x_anchor + bar_w_npc / 2
  bar_cy    <- bar_y_top - bar_h_npc / 2

  n_steps   <- 128L
  grad_cols <- grDevices::colorRampPalette(.palette_stops(spec$palette))(n_steps)
  if (is_horiz) {
    img <- matrix(grad_cols, nrow = 1L, ncol = n_steps)
  } else {
    img <- matrix(rev(grad_cols), nrow = n_steps, ncol = 1L)
  }

  grid::grid.raster(
    image       = img,
    x           = grid::unit(bar_cx, "npc"),
    y           = grid::unit(bar_cy, "npc"),
    width       = grid::unit(bar_w_npc, "npc"),
    height      = grid::unit(bar_h_npc, "npc"),
    interpolate = TRUE
  )

  # Explicit corner coords used for tick placement
  bar_x0 <- bar_cx - bar_w_npc / 2
  bar_x1 <- bar_cx + bar_w_npc / 2
  bar_y0 <- bar_cy - bar_h_npc / 2  # bottom (low value for vertical)
  bar_y1 <- bar_cy + bar_h_npc / 2  # top    (high value for vertical)

  # --- Axis ticks and labels ---
  tick_len  <- 0.008
  label_gap <- 0.004

  lims <- spec$limits %||% c(0, 1)
  brks <- spec$breaks
  if (is.null(brks)) {
    brks <- pretty(lims, n = 5L)
  } else if (length(brks) == 1L && is.numeric(brks)) {
    brks <- pretty(lims, n = as.integer(brks))
  } else {
    brks <- as.numeric(brks)
  }
  brks <- brks[is.finite(brks) & brks >= lims[1] & brks <= lims[2]]

  span <- diff(lims)
  if (span == 0) span <- 1

  for (b in brks) {
    t <- max(0, min(1, (b - lims[1]) / span))

    if (is_horiz) {
      tick_x <- bar_x0 + t * (bar_x1 - bar_x0)
      grid::grid.lines(
        x  = grid::unit(c(tick_x, tick_x), "npc"),
        y  = grid::unit(c(bar_y0 - tick_len, bar_y0), "npc"),
        gp = grid::gpar(col = label_col, lwd = 0.8)
      )
      grid::grid.text(
        label = .axis_format_num(b),
        x     = grid::unit(tick_x, "npc"),
        y     = grid::unit(bar_y0 - tick_len - label_gap, "npc"),
        just  = c("center", "top"),
        gp    = grid::gpar(col = label_col, cex = key_cex)
      )
    } else {
      tick_y <- bar_y0 + t * (bar_y1 - bar_y0)
      grid::grid.lines(
        x  = grid::unit(c(bar_x1, bar_x1 + tick_len), "npc"),
        y  = grid::unit(c(tick_y, tick_y), "npc"),
        gp = grid::gpar(col = label_col, lwd = 0.8)
      )
      grid::grid.text(
        label = .axis_format_num(b),
        x     = grid::unit(bar_x1 + tick_len + label_gap, "npc"),
        y     = grid::unit(tick_y, "npc"),
        just  = c("left", "center"),
        gp    = grid::gpar(col = label_col, cex = key_cex)
      )
    }
  }

  invisible()
}


# Dispatch: draw a GradientLegendSpec into the panel given by panel_meta.
# Always renders as a continuous color bar with a tick+label axis.
# spec$breaks controls tick positions (NULL = pretty(), integer = n ticks,
# numeric vector = explicit positions).
.draw_gradient_legend <- function(spec, panel_meta, ...) {
  if (!inherits(spec, "GradientLegendSpec")) return(invisible())
  .draw_gradient_colorbar(spec, panel_meta, ...)
}
