# seq_preview_layout ----

#' Preview the grid layout of a SeqPlotR plot
#'
#' @description
#' Renders a schematic toy representation of the intended plot layout — no data
#' required. Each track region is drawn as a filled, labeled rectangle.
#' Useful for verifying a layout string or operator chain before adding data.
#'
#' @param plot_obj A `seq_plot` object built with the operator chain, or `NULL`.
#'   If provided, the layout is extracted from the object.
#' @param layout A multiline layout string (e.g. `"##AA\n##AA\nBBBC\nBBBD"`),
#'   or `NULL`. Used when `plot_obj` is not provided.
#' @param labels Logical. If `TRUE` (default), draw the track ID label centered
#'   in each region.
#' @param colors Optional named character vector mapping track IDs to fill
#'   colors (e.g. `c(A = "#FF0000", B = "#00FF00")`). Unspecified IDs receive
#'   automatic colors.
#' @param margins Logical. If `TRUE` (default), draw a thin border rectangle
#'   indicating the outer canvas margin.
#' @param margin_size Numeric in the range 0 to 0.5. Fractional canvas margin
#'   on each side. Default `0.05`.
#'
#' @return A named list of npc bounding boxes (invisibly):
#'   each element is `list(x0, x1, y0, y1)` keyed by track ID.
#'
#' @examples
#' \dontrun{
#' # Preview a patchwork layout string directly
#' layout <- "
#' ##AA
#' ##AA
#' BBBC
#' BBBD
#' "
#' seq_preview_layout(layout = layout)
#'
#' # Preview a layout built with operators
#' p <- seq_plot() |>
#'   (\(x) x %|% seq_track(track_id = "Signal",  track_width = 3))() |>
#'   (\(x) x %|% seq_track(track_id = "Ideogram", track_width = 1))() |>
#'   (\(x) x %__% seq_track(track_id = "Genes"))()
#' seq_preview_layout(plot_obj = p)
#' }
#'
#' @export
seq_preview_layout <- function(plot_obj = NULL,
                                layout   = NULL,
                                labels   = TRUE,
                                colors   = NULL,
                                margins  = TRUE,
                                margin_size = 0.05) {

  # ── 1. Obtain region definitions ──────────────────────────────────────────
  if (!is.null(plot_obj)) {
    if (!inherits(plot_obj, "SeqPlot"))
      stop("plot_obj must be a seq_plot object (class 'SeqPlot').")
    regions <- .extract_regions_from_plot(plot_obj, margin_size)
  } else if (!is.null(layout)) {
    if (!is.character(layout) || length(layout) != 1)
      stop("layout must be a single character string.")
    parsed  <- .parse_layout_string(layout)
    regions <- .regions_from_parse(parsed, margin_size)
  } else {
    stop("Either plot_obj or layout must be provided.")
  }

  if (length(regions) == 0)
    stop("No track regions found. Check that track_ids or layout letters are specified.")

  # ── 2. Assign colors ──────────────────────────────────────────────────────
  track_ids  <- names(regions)
  fill_colors <- .preview_colors(track_ids, colors)

  # ── 3. Render ─────────────────────────────────────────────────────────────
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(name = "preview_root"))

  # Outer margin border
  if (isTRUE(margins)) {
    grid::grid.rect(
      x      = grid::unit(0.5, "npc"),
      y      = grid::unit(0.5, "npc"),
      width  = grid::unit(1 - 2 * margin_size, "npc"),
      height = grid::unit(1 - 2 * margin_size, "npc"),
      gp     = grid::gpar(fill = NA, col = "grey70", lwd = 0.8, lty = "dashed")
    )
  }

  # Draw blank (#) cells as light grey background first
  if (!is.null(attr(regions, "blank_cells"))) {
    for (cell in attr(regions, "blank_cells")) {
      grid::grid.rect(
        x      = grid::unit((cell$x0 + cell$x1) / 2, "npc"),
        y      = grid::unit((cell$y0 + cell$y1) / 2, "npc"),
        width  = grid::unit(cell$x1 - cell$x0 - 0.005, "npc"),
        height = grid::unit(cell$y1 - cell$y0 - 0.005, "npc"),
        just   = "centre",
        gp     = grid::gpar(fill = "#F0F0F0", col = NA)
      )
    }
  }

  # Draw track regions
  for (id in track_ids) {
    r   <- regions[[id]]
    col <- fill_colors[[id]]

    # Filled rectangle with small inset gap
    gap <- 0.004
    grid::grid.rect(
      x      = grid::unit((r$x0 + r$x1) / 2, "npc"),
      y      = grid::unit((r$y0 + r$y1) / 2, "npc"),
      width  = grid::unit(r$x1 - r$x0 - gap, "npc"),
      height = grid::unit(r$y1 - r$y0 - gap, "npc"),
      just   = "centre",
      gp     = grid::gpar(fill = col, col = "white", lwd = 1.5)
    )

    # Label
    if (isTRUE(labels)) {
      # Choose text color based on background luminance
      text_col <- .contrast_color(col)
      # Scale font size to region size (clamped between 0.5 and 1.2 relative cex)
      region_size <- min(r$x1 - r$x0, r$y1 - r$y0)
      cex <- max(0.5, min(1.2, region_size * 5))

      grid::grid.text(
        label = id,
        x     = grid::unit((r$x0 + r$x1) / 2, "npc"),
        y     = grid::unit((r$y0 + r$y1) / 2, "npc"),
        just  = "centre",
        gp    = grid::gpar(col = text_col, cex = cex, fontface = "bold")
      )
    }
  }

  grid::popViewport()

  invisible(regions)
}


# ── Internal helpers ──────────────────────────────────────────────────────────
#
# `.parse_layout_string()` is defined in R/layout.R as the canonical parser.
# Helpers below convert the parsed structure into preview-specific npc bounds.

#' Convert parsed layout to npc bounding boxes
#'
#' @param parsed Output of `.parse_layout_string()`.
#' @param margin Numeric fractional margin on each side.
#' @return Named list of `list(x0, x1, y0, y1)` per track ID.
#' @keywords internal
.regions_from_parse <- function(parsed, margin = 0.05) {
  nrow   <- parsed$nrow
  ncol   <- parsed$ncol
  canvas <- list(
    x0 = margin, x1 = 1 - margin,
    y0 = margin, y1 = 1 - margin
  )
  cw <- canvas$x1 - canvas$x0
  ch <- canvas$y1 - canvas$y0

  regions <- lapply(parsed$regions, function(reg) {
    x0 <- canvas$x0 + (reg$c0 - 1) / ncol * cw
    x1 <- canvas$x0 +  reg$c1       / ncol * cw
    # Layout string is top-to-bottom; npc y is bottom-to-top
    y0 <- canvas$y0 + (1 - reg$r1 / nrow) * ch
    y1 <- canvas$y0 + (1 - (reg$r0 - 1) / nrow) * ch
    list(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
  })

  # Attach blank cell npc boxes as attribute
  blank_npc <- lapply(parsed$blank_cells, function(bc) {
    x0 <- canvas$x0 + (bc$col - 1) / ncol * cw
    x1 <- canvas$x0 +  bc$col       / ncol * cw
    y0 <- canvas$y0 + (1 - bc$row / nrow) * ch
    y1 <- canvas$y0 + (1 - (bc$row - 1) / nrow) * ch
    list(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
  })
  attr(regions, "blank_cells") <- blank_npc

  regions
}


#' Extract layout regions from a seq_plot object
#'
#' Handles both patchwork (layout_str) and positional (rows) modes.
#'
#' @param plot_obj A `SeqPlot` R6 object.
#' @param margin Numeric fractional margin.
#' @return Named list of npc bounding boxes per track ID.
#' @keywords internal
.extract_regions_from_plot <- function(plot_obj, margin = 0.05) {

  # Patchwork mode
  if (!is.null(plot_obj$layout_str)) {
    parsed <- .parse_layout_string(plot_obj$layout_str)
    return(.regions_from_parse(parsed, margin))
  }

  # Positional mode — reconstruct from rows
  rows <- plot_obj$rows
  if (is.null(rows) || length(rows) == 0)
    stop("seq_plot has no tracks and no layout string.")

  canvas <- list(x0 = margin, x1 = 1 - margin,
                 y0 = margin, y1 = 1 - margin)
  cw <- canvas$x1 - canvas$x0
  ch <- canvas$y1 - canvas$y0

  nrows <- length(rows)

  # Row heights: max track_height per row, normalized
  row_heights <- vapply(rows, function(row) {
    max(vapply(row, function(trk) trk$track_height %||% 1, numeric(1)))
  }, numeric(1))
  row_heights_norm <- row_heights / sum(row_heights)

  regions <- list()
  y_top <- canvas$y1

  for (ri in seq_len(nrows)) {
    row  <- rows[[ri]]
    rh   <- row_heights_norm[ri] * ch
    y_bottom <- y_top - rh

    # Track widths within this row
    tw_raw  <- vapply(row, function(trk) trk$track_width %||% 1, numeric(1))
    tw_norm <- tw_raw / sum(tw_raw)

    x_left <- canvas$x0
    for (ti in seq_along(row)) {
      trk  <- row[[ti]]
      tw   <- tw_norm[ti] * cw
      id   <- trk$track_id %||% paste0("R", ri, "T", ti)

      regions[[id]] <- list(
        x0 = x_left,
        x1 = x_left + tw,
        y0 = y_bottom,
        y1 = y_top
      )
      x_left <- x_left + tw
    }
    y_top <- y_bottom
  }

  regions
}


#' Assign preview fill colors to track IDs
#'
#' @param ids Character vector of track IDs.
#' @param user_colors Optional named character vector overriding specific IDs.
#' @return Named character vector of hex colors.
#' @keywords internal
.preview_colors <- function(ids, user_colors = NULL) {
  palette <- c(
    "#A8C5DA", "#F4C6A0", "#B5D5C5", "#F2E0B6",
    "#C9B8D8", "#F7B7B7", "#B8D4E8", "#D4E8B8",
    "#E8D4B8", "#D8B8C9", "#C5DAA8", "#DAC5A8"
  )
  n <- length(ids)
  colors <- setNames(
    palette[(seq_len(n) - 1) %% length(palette) + 1],
    ids
  )
  if (!is.null(user_colors)) {
    valid <- names(user_colors) %in% ids
    if (!all(valid))
      warning("colors contains unknown track IDs: ",
              paste(names(user_colors)[!valid], collapse = ", "))
    colors[names(user_colors)[valid]] <- user_colors[valid]
  }
  colors
}


#' Choose black or white text for contrast against a background color
#'
#' Uses the WCAG relative luminance formula.
#'
#' @param hex_color A single hex color string (e.g. `"#A8C5DA"`).
#' @return `"#1C1B1A"` (dark) or `"#FFFFFF"` (light).
#' @keywords internal
.contrast_color <- function(hex_color) {
  tryCatch({
    rgb_vals <- grDevices::col2rgb(hex_color) / 255
    # Linearize
    lin <- ifelse(rgb_vals <= 0.04045,
                  rgb_vals / 12.92,
                  ((rgb_vals + 0.055) / 1.055) ^ 2.4)
    L <- 0.2126 * lin[1] + 0.7152 * lin[2] + 0.0722 * lin[3]
    if (L > 0.179) "#1C1B1A" else "#FFFFFF"
  }, error = function(e) "#1C1B1A")
}


# seq_preview_circos ----

#' Preview a SeqPlotR plot as a circos layout
#'
#' @description
#' Renders a schematic circular layout — no data required. Each track is shown
#' as a colored arc sector. Rows (defined via `direction = "under"` / `%__%`)
#' map to concentric rings from outer to inner; tracks within a row (defined
#' via `direction = "right"` / `%|%`) map to arc sectors with angular width
#' proportional to `track_width`. Ring radial thickness is proportional to the
#' maximum `track_height` in that row.
#'
#' Only positional layouts are supported — a plot built with a patchwork
#' layout string has no natural circos interpretation and is rejected.
#'
#' @param plot_obj A `seq_plot` object built with the operator chain. Required.
#' @param labels Logical. Draw `track_id` labels. Default `TRUE`.
#' @param colors Optional named character vector mapping track IDs to fill
#'   colors. Unspecified IDs receive automatic colors.
#' @param start_angle Numeric clock-face degrees where the first sector starts.
#'   Default `90` (top). Sectors sweep clockwise (decreasing angle).
#' @param end_angle Numeric clock-face degrees where the last sector ends.
#'   Default `-270` (full circle = `start_angle - 360`).
#' @param gap_degrees Numeric blank gap in degrees between sectors of the same
#'   ring. Default `2`.
#' @param ring_gap Numeric npc blank gap between concentric rings. Default
#'   `0.02`.
#' @param outer_radius Numeric npc radius of the outermost ring's outer edge.
#'   Default `0.45`.
#' @param inner_radius Numeric npc radius of the innermost ring's inner edge.
#'   Default `0.08`.
#' @param cx,cy Numeric npc coordinates of the circle centre. Default
#'   `(0.5, 0.5)`.
#'
#' @return Named list of polar bounding boxes (invisibly). Each entry is
#'   `list(theta0, theta1, r0, r1)` in clock-face degrees and npc radii,
#'   keyed by track ID.
#'
#' @examples
#' \dontrun{
#' p <- seq_plot() %|%
#'   seq_track(track_id = "Chr1", track_width = 3) %|%
#'   seq_track(track_id = "Chr2", track_width = 2) %|%
#'   seq_track(track_id = "Chr3", track_width = 1) %__%
#'   seq_track(track_id = "Signal") %__%
#'   seq_track(track_id = "CopyNum")
#' seq_preview_circos(plot_obj = p)
#' }
#'
#' @export
seq_preview_circos <- function(plot_obj     = NULL,
                               labels       = TRUE,
                               colors       = NULL,
                               start_angle  = 90,
                               end_angle    = -270,
                               gap_degrees  = 2,
                               ring_gap     = 0.02,
                               outer_radius = 0.45,
                               inner_radius = 0.08,
                               cx           = 0.5,
                               cy           = 0.5) {

  # ── 1. Validate + extract rows ──────────────────────────────────────────────
  if (is.null(plot_obj))
    stop("plot_obj must be provided for seq_preview_circos().",
         call. = FALSE)
  if (!inherits(plot_obj, "SeqPlot"))
    stop("plot_obj must be a seq_plot object (class 'SeqPlot').",
         call. = FALSE)
  if (!is.null(plot_obj$layout_str))
    stop("seq_preview_circos() does not support patchwork layout strings. ",
         "Use a positional layout (%|% / %__%) to build the plot.",
         call. = FALSE)

  rows <- plot_obj$rows
  if (is.null(rows) || length(rows) == 0)
    stop("seq_plot has no tracks. Add tracks with %|% or %__%.",
         call. = FALSE)

  # ── 2. Compute radial ring spans ────────────────────────────────────────────
  n_rings     <- length(rows)
  row_heights <- vapply(rows, function(row)
    max(vapply(row, function(trk) trk$track_height %||% 1, numeric(1))),
    numeric(1))
  usable_r    <- (outer_radius - inner_radius) -
                 ring_gap * max(0, n_rings - 1)
  rh_norm     <- row_heights / sum(row_heights)
  ring_heights_npc <- rh_norm * usable_r

  # ── 3. Compute angular sector spans ─────────────────────────────────────────
  total_degrees <- start_angle - end_angle  # 360 for a full circle

  polar_bounds <- list()
  all_ids      <- character(0)
  r_cursor     <- outer_radius

  for (ri in seq_len(n_rings)) {
    row <- rows[[ri]]
    r1  <- r_cursor
    r0  <- r_cursor - ring_heights_npc[ri]
    r_cursor <- r0 - ring_gap

    n_sectors  <- length(row)
    tw_raw     <- vapply(row, function(trk) trk$track_width %||% 1,
                         numeric(1))
    tw_norm    <- tw_raw / sum(tw_raw)
    usable_deg <- total_degrees - gap_degrees * max(0, n_sectors - 1)
    sector_deg <- tw_norm * usable_deg

    theta_cursor <- start_angle
    for (ti in seq_along(row)) {
      trk    <- row[[ti]]
      id     <- trk$track_id %||% paste0("R", ri, "T", ti)
      theta0 <- theta_cursor
      theta1 <- theta_cursor - sector_deg[ti]
      theta_cursor <- theta1 - gap_degrees

      polar_bounds[[id]] <- list(theta0 = theta0, theta1 = theta1,
                                 r0     = r0,    r1     = r1)
      all_ids <- c(all_ids, id)
    }
  }

  # ── 4. Assign fill colors ───────────────────────────────────────────────────
  fill_colors <- .preview_colors(all_ids, colors)

  # ── 5. Render ───────────────────────────────────────────────────────────────
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(name = "circos_preview_root"))

  # Outer boundary circle (dashed)
  theta_seq <- seq(0, 360, length.out = 361)
  outer_xs  <- cx + outer_radius * cos(theta_seq * pi / 180)
  outer_ys  <- cy + outer_radius * sin(theta_seq * pi / 180)
  grid::grid.lines(
    x  = grid::unit(outer_xs, "npc"),
    y  = grid::unit(outer_ys, "npc"),
    gp = grid::gpar(col = "grey70", lwd = 0.8, lty = "dashed")
  )

  # Arc polygons
  for (id in all_ids) {
    pb   <- polar_bounds[[id]]
    poly <- .arc_polygon(pb$theta0, pb$theta1, pb$r0, pb$r1, cx, cy)
    grid::grid.polygon(
      x  = grid::unit(poly$x, "npc"),
      y  = grid::unit(poly$y, "npc"),
      gp = grid::gpar(fill = fill_colors[[id]], col = "white", lwd = 1.2)
    )
  }

  # Labels
  if (isTRUE(labels)) {
    for (id in all_ids) {
      pb        <- polar_bounds[[id]]
      theta_mid <- (pb$theta0 + pb$theta1) / 2
      r_mid     <- (pb$r0 + pb$r1) / 2
      arc_span  <- abs(pb$theta0 - pb$theta1)
      if (arc_span < 3) next  # too narrow to label

      lpos     <- .polar_to_npc(r_mid, theta_mid, cx, cy)
      text_col <- .contrast_color(fill_colors[[id]])

      arc_len_approx <- r_mid * arc_span * pi / 180
      cex <- max(0.35, min(0.9, arc_len_approx * 4))

      # Tangent to the arc (perpendicular to radius); flip to keep upright.
      rot <- theta_mid + 90
      if (rot %% 360 > 180) rot <- rot - 180

      grid::grid.text(
        label = id,
        x     = grid::unit(lpos$x, "npc"),
        y     = grid::unit(lpos$y, "npc"),
        rot   = rot,
        gp    = grid::gpar(col = text_col, cex = cex, fontface = "bold")
      )
    }
  }

  grid::popViewport()
  invisible(polar_bounds)
}


# ── Circos-specific helpers ──────────────────────────────────────────────────

#' Convert polar coordinates to npc Cartesian
#'
#' Mathematical convention: 0 degrees = 3 o'clock, counter-clockwise positive.
#'
#' @param r Numeric npc radius.
#' @param theta_deg Numeric angle in degrees.
#' @param cx,cy Numeric npc coordinates of the circle centre.
#' @return `list(x, y)` in npc.
#' @keywords internal
.polar_to_npc <- function(r, theta_deg, cx = 0.5, cy = 0.5) {
  theta_rad <- theta_deg * pi / 180
  list(x = cx + r * cos(theta_rad),
       y = cy + r * sin(theta_rad))
}

#' Build a filled annular-sector polygon as xy vertex vectors
#'
#' Samples `n_pts` points along the outer arc (radius `r1`) from `theta0_deg`
#' to `theta1_deg`, then `n_pts` points along the inner arc (radius `r0`)
#' reversed to close the polygon.
#'
#' @param theta0_deg,theta1_deg Sector angular bounds (degrees).
#' @param r0,r1 Inner and outer radii (npc).
#' @param cx,cy Circle centre (npc).
#' @param n_pts Samples per arc edge. Default `60`.
#' @return `list(x, y)` of length `2 * n_pts` each.
#' @keywords internal
.arc_polygon <- function(theta0_deg, theta1_deg, r0, r1,
                         cx = 0.5, cy = 0.5, n_pts = 60) {
  thetas  <- seq(theta0_deg, theta1_deg, length.out = n_pts)
  outer   <- lapply(thetas,      function(th) .polar_to_npc(r1, th, cx, cy))
  inner   <- lapply(rev(thetas), function(th) .polar_to_npc(r0, th, cx, cy))
  all_pts <- c(outer, inner)
  list(
    x = vapply(all_pts, `[[`, numeric(1), "x"),
    y = vapply(all_pts, `[[`, numeric(1), "y")
  )
}
