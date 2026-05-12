# ── Layout engine ────────────────────────────────────────────────────────────
#
# Three internal entry points used by SeqPlotR6$layoutGrid():
#
#   .parse_layout_string()      — parse a multiline patchwork string
#   .build_positional_layout()  — build panel/track bounds from row structure
#   .build_patchwork_layout()   — build panel/track bounds from a layout string
#
# Both layout builders return a list with the same shape as THEfunc's
# layoutGrid() output:
#
#   list(
#     panelBounds = <list of per-track lists of per-window panel metadata>,
#     trackBounds = <list of per-track bounding boxes>
#   )
#
# panelBounds entries (per window) include `full`, `inner`, `xscale`, `yscale`,
# `window`, `track`, `xScaleFactor`, `y_scale_type`, and `y_sub_panels`. The
# field set is intentionally a strict superset of what the layout tests inspect
# so that downstream element prep() methods (Batch 4) can rely on the same
# panel metadata structure.

# ── .parse_layout_string ─────────────────────────────────────────────────────

#' Parse a patchwork layout string
#'
#' Converts a multiline character string into a structured layout
#' description. Each non-`#` letter denotes a region; `#` cells are
#' rendered as blank.
#'
#' @param s A single multiline character string. Rows are split on `"\n"`,
#'   trimmed, and dropped if empty. All rows must have the same number of
#'   characters.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`nrow`}{number of rows}
#'     \item{`ncol`}{number of columns}
#'     \item{`regions`}{named list, one entry per unique non-`#` letter, each
#'       a list `list(r0, r1, c0, c1)` (1-indexed, inclusive)}
#'     \item{`blank_cells`}{list of `list(row, col)` for every `#` cell}
#'   }
#'
#' @keywords internal
.parse_layout_string <- function(s) {
  lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L)
    stop("Layout string is empty.", call. = FALSE)

  widths <- nchar(lines)
  if (length(unique(widths)) != 1L)
    stop("Layout string rows have unequal length: ",
         paste(widths, collapse = ", "), ".", call. = FALSE)

  nrow_ <- length(lines)
  ncol_ <- widths[1]
  m <- matrix(unlist(strsplit(lines, "", fixed = TRUE)),
              nrow = nrow_, ncol = ncol_, byrow = TRUE)

  blank_cells <- list()
  for (i in seq_len(nrow_)) for (j in seq_len(ncol_)) {
    if (m[i, j] == "#")
      blank_cells <- c(blank_cells, list(list(row = i, col = j)))
  }

  letters_ <- setdiff(unique(as.vector(m)), "#")
  regions  <- list()
  for (L in letters_) {
    pos <- which(m == L, arr.ind = TRUE)
    r0 <- min(pos[, 1]); r1 <- max(pos[, 1])
    c0 <- min(pos[, 2]); c1 <- max(pos[, 2])
    sub <- m[r0:r1, c0:c1, drop = FALSE]
    if (!all(sub == L))
      stop("Layout region '", L,
           "' is not rectangular - every cell in its bounding box must be '",
           L, "'.", call. = FALSE)
    regions[[L]] <- list(r0 = r0, r1 = r1, c0 = c0, c1 = c1)
  }

  list(nrow = nrow_, ncol = ncol_,
       regions = regions, blank_cells = blank_cells)
}

# ── shared helpers ───────────────────────────────────────────────────────────

#' Compute the y-axis scale for a track
#'
#' Used by the layout builders. Dispatches on (1) genomic y-windows, (2)
#' explicit `scale_y` type, and (3) a `c(0, 1)` placeholder fallback.
#'
#' @param track A `SeqTrackR6` instance.
#' @return Numeric length-2 vector `c(y_min, y_max)`.
#' @keywords internal
.compute_track_yscale <- function(track) {
  if (isTRUE(track$uses_genomic_y) && !is.null(track$y_windows)) {
    return(c(min(BiocGenerics::start(track$y_windows)),
             max(BiocGenerics::end(track$y_windows))))
  }
  sy <- track$scale_y
  if (is.null(sy)) return(c(0, 1))
  if (!is.null(sy$type)) {
    if (sy$type == "genomic" && !is.null(sy$windows)) {
      return(c(min(BiocGenerics::start(sy$windows)),
               max(BiocGenerics::end(sy$windows))))
    }
    if (sy$type == "discrete" && !is.null(sy$levels)) {
      return(c(0.5, length(sy$levels) + 0.5))
    }
    if (sy$type == "continuous" && !is.null(sy$limits)) {
      return(as.numeric(sy$limits))
    }
  }
  c(0, 1)
}

#' Infer a shared x-axis scale factor for a track
#'
#' Uses the narrowest window to determine the appropriate unit (Mb/kb/bp)
#' so that all windows in the track share the same scale.
#'
#' @param widths Integer vector of genomic window widths in bp.
#' @return A single numeric scale factor: `1e-6` (Mb), `1e-3` (kb), or `1`
#'   (bp).
#' @keywords internal
.infer_track_scale_factor <- function(widths) {
  narrowest <- min(widths, na.rm = TRUE)
  if      (narrowest >= 1e6) 1e-6   # >= 1 Mb   -> Mb
  else if (narrowest >= 1e2) 1e-3   # >= 100 bp -> kb
  else                       1      # < 100 bp  -> bp
}

#' Build per-window relative widths and resolve scale factors
#'
#' @param win A `GRanges` object of windows.
#' @param track Optional `SeqTrackR6`. When provided, `window_scale` is
#'   read from the track and the track-level inferred unit is used as the
#'   default. When `NULL`, falls back to `1e-6` (original behaviour).
#' @return A list with `rel` and `scale` numeric vectors of length
#'   `length(win)`.
#' @keywords internal
.window_relative_widths <- function(win, track = NULL) {
  raw <- BiocGenerics::width(win)
  n   <- length(win)

  # Priority 1: explicit mcols$scale on the GRanges (user-set, wins always)
  if ("scale" %in% names(S4Vectors::mcols(win))) {
    scale_factor <- S4Vectors::mcols(win)$scale
    eff <- raw * scale_factor
    return(list(rel = eff / sum(eff), scale = scale_factor))
  }

  # Priority 2: window_scale on the SeqTrack
  if (!is.null(track) && !is.null(track$window_scale)) {
    ws <- track$window_scale
    if (length(ws) == 1L) {
      scale_factor <- rep(ws, n)
    } else if (length(ws) == n) {
      scale_factor <- ws
    } else {
      warning(
        "`window_scale` has length ", length(ws), " but track has ", n,
        " window(s). Recycling with rep_len().",
        call. = FALSE
      )
      scale_factor <- rep_len(ws, n)
    }
    eff <- raw * scale_factor
    return(list(rel = eff / sum(eff), scale = scale_factor))
  }

  # Priority 3: infer from the narrowest window (new default)
  if (!is.null(track)) {
    sf <- .infer_track_scale_factor(raw)
    scale_factor <- rep(sf, n)
    eff <- raw * scale_factor
    return(list(rel = eff / sum(eff), scale = scale_factor))
  }

  # Priority 4: legacy fallback (no track provided - keep original behaviour)
  scale_factor <- rep(1e-6, n)
  eff <- raw * scale_factor
  list(rel = eff / sum(eff), scale = scale_factor)
}

#' Build the panel metadata list for one track in one (x, y) bounding box
#'
#' Splits the track's x-band into per-window panels and packages each one
#' with `full`, `inner`, `xscale`, `yscale`, and other downstream-required
#' fields. Used by both layout builders.
#'
#' @param track A `SeqTrackR6` instance.
#' @param x0,x1 Numeric npc x bounds for the track.
#' @param y0,y1 Numeric npc y bounds for the track.
#' @param window_gap Numeric npc gap between adjacent window panels.
#' @param track_key The value to store in each panel's `track` field — integer
#'   index in positional mode, `track_id` string in patchwork mode.
#' @return A list of per-window panel metadata.
#' @keywords internal
.build_track_panels <- function(track, x0, x1, y0, y1,
                                window_gap, track_key) {
  # combine_windows: collapse multiple genomic windows into a single
  # virtual panel. Each window keeps its own genomic scale internally
  # (via the virtual map) so axes and per-window labels still render,
  # but cross-window data (e.g. inter-chromosomal Hi-C contacts) lays
  # out continuously within one track inner area.
  vmap_x <- NULL
  if (isTRUE(track$combine_windows) && length(track$windows) > 1L) {
    vmap_x        <- .build_virtual_map(track$windows, virtual_gap = 0)
    win           <- vmap_x$combined_window
  } else {
    win           <- track$windows
  }
  ww   <- .window_relative_widths(win, track = track)
  nWin <- length(win)

  # Margin shrink helper
  shrink <- function(box, m) {
    r <- list(
      x0 = box$x0 + m$left,   x1 = box$x1 - m$right,
      y0 = box$y0 + m$bottom, y1 = box$y1 - m$top
    )
    if (r$x1 < r$x0) r$x1 <- r$x0
    if (r$y1 < r$y0) r$y1 <- r$y0
    r
  }

  # Track-level zones: full → outer → inner (holds all windows)
  tom <- .normalize_margin(track$track_outer_margin)
  tim <- .normalize_margin(track$track_inner_margin)
  track_full  <- list(x0 = x0, x1 = x1, y0 = y0, y1 = y1)
  track_outer <- shrink(track_full,  tom)   # axis titles live between full & outer
  track_inner <- shrink(track_outer, tim)   # contains all windows

  # Window gap: prefer the track's resolved (plot+track) flat theme so
  # aes("window.gap.width") set at the track level overrides the plot
  # level. Falls back to the plot-level `window_gap` arg.
  trk_flat <- track$resolved_theme$flat %||% list()
  wgap <- trk_flat[["window.gap.width"]] %||%
          trk_flat$window_gaps           %||% window_gap
  # `track$window_margin` is deprecated and intentionally not consulted —
  # see the warning in SeqTrackR6$initialize().
  win_gap_total <- if (nWin > 1) wgap * (nWin - 1) else 0
  plot_w        <- (track_inner$x1 - track_inner$x0) - win_gap_total
  if (plot_w < 0) plot_w <- 0

  wom <- .normalize_margin(track$window_outer_margin)
  wim <- .normalize_margin(track$window_inner_margin)

  yscale  <- .compute_track_yscale(track)
  yscale2 <- .scale_to_yrange(track$scale_y2, track$y_windows2)
  y_scale_type <- if (isTRUE(track$uses_genomic_y)) "genomic"
                  else if (!is.null(track$scale_y))  track$scale_y$type
                  else                                "continuous"

  # Build y sub-panels when the track has multiple y-windows on the
  # genomic y axis (and combine_y_windows is FALSE). Each sub-panel
  # gets its own npc band and per-window yscale; tiles route into
  # whichever sub-panel matches their j-position seqname.
  y_sub_panels_default <- NULL
  vmap_y <- NULL
  if (isTRUE(track$uses_genomic_y) &&
      !is.null(track$y_windows) &&
      length(track$y_windows) > 1L &&
      !isTRUE(track$combine_y_windows)) {
    y_win    <- track$y_windows
    n_yw     <- length(y_win)
    y_w      <- BiocGenerics::width(y_win)
    y_rel    <- y_w / sum(y_w)
    y_gap    <- wgap   # same window.gap.width applies to y sub-panels
    h_total  <- (track_inner$y1 - track_inner$y0) - y_gap * (n_yw - 1L)
    if (h_total < 0) h_total <- 0
    y_sub_panels_default <- vector("list", n_yw)
    y_cur <- track_inner$y0
    for (k in seq_len(n_yw)) {
      if (k > 1L) y_cur <- y_cur + y_gap
      y0_sub <- y_cur
      y1_sub <- y_cur + y_rel[k] * h_total
      y_cur  <- y1_sub
      ys_k <- c(BiocGenerics::start(y_win[k]),
                BiocGenerics::end(y_win[k]))
      y_sub_panels_default[[k]] <- list(
        y0      = y0_sub,
        y1      = y1_sub,
        yscale  = ys_k,
        yplot_range = if (!is.null(track$scale_y))
          .compute_scale_breaks(track$scale_y, ys_k)$plot_range else ys_k,
        seqname = as.character(GenomicRanges::seqnames(y_win[k]))
      )
    }
    # Track-wide yscale spans the union for any single-axis fallbacks.
    yscale <- c(min(BiocGenerics::start(y_win)),
                max(BiocGenerics::end(y_win)))
  }
  # combine_y_windows: collapse multi-y-windows into a single virtual
  # y axis spanning the concatenation. Tiles' j positions need to be
  # virtualized (the seq_hic wrapper does this when combine_windows
  # is TRUE on a square matrix). The yscale becomes the virtual range
  # and `virtual_map_y` is stashed on each panel for axis drawing.
  if (isTRUE(track$uses_genomic_y) &&
      !is.null(track$y_windows) &&
      length(track$y_windows) > 1L &&
      isTRUE(track$combine_y_windows)) {
    vmap_y <- .build_virtual_map(track$y_windows, virtual_gap = 0)
    yscale <- c(1, vmap_y$virtual_total)
  }

  panels <- vector("list", nWin)
  x_cur  <- track_inner$x0
  # Pre-compute oob modes (default exclude) from the scales.
  x_oob  <- if (!is.null(track$scale_x)  && !is.null(track$scale_x$oob))
              track$scale_x$oob  else "exclude"
  y_oob  <- if (!is.null(track$scale_y)  && !is.null(track$scale_y$oob))
              track$scale_y$oob  else "exclude"
  x_oob2 <- if (!is.null(track$scale_x2) && !is.null(track$scale_x2$oob))
              track$scale_x2$oob else "exclude"
  y_oob2 <- if (!is.null(track$scale_y2) && !is.null(track$scale_y2$oob))
              track$scale_y2$oob else "exclude"
  for (w in seq_len(nWin)) {
    w_width <- ww$rel[w] * plot_w
    px0 <- x_cur
    px1 <- x_cur + w_width
    window_full  <- list(x0 = px0, x1 = px1,
                         y0 = track_inner$y0, y1 = track_inner$y1)
    window_outer <- shrink(window_full,  wom)  # after window_outer_margin
    plot_area    <- shrink(window_outer, wim)  # after window_inner_margin
    xscale_w  <- .compute_track_xscale(track, win, w)
    xscale2_w <- .scale_to_xrange(track$scale_x2, win, w)

    # Compute the expanded plot range that elements should map against.
    xplot_range_w <- if (!is.null(track$scale_x))
      .compute_scale_breaks(track$scale_x, xscale_w)$plot_range
      else xscale_w
    yplot_range_w <- if (!is.null(track$scale_y))
      .compute_scale_breaks(track$scale_y, yscale)$plot_range
      else yscale
    xplot_range2_w <- if (!is.null(track$scale_x2) && !is.null(xscale2_w))
      .compute_scale_breaks(track$scale_x2, xscale2_w)$plot_range
      else xscale2_w
    yplot_range2_w <- if (!is.null(track$scale_y2) && !is.null(yscale2))
      .compute_scale_breaks(track$scale_y2, yscale2)$plot_range
      else yscale2

    panels[[w]] <- list(
      track               = track_key,
      window              = w,
      full                = window_full,
      inner               = plot_area,
      window_outer        = window_outer,
      track_full          = track_full,
      track_outer         = track_outer,
      track_inner         = track_inner,
      track_mapping       = track$mapping,
      track_outer_margin  = tom,
      track_inner_margin  = tim,
      window_outer_margin = wom,
      window_inner_margin = wim,
      xscale              = xscale_w,
      yscale              = yscale,
      xscale2             = xscale2_w,
      yscale2             = yscale2,
      xplot_range         = xplot_range_w,
      yplot_range         = yplot_range_w,
      xplot_range2        = xplot_range2_w,
      yplot_range2        = yplot_range2_w,
      x_oob               = x_oob,
      y_oob               = y_oob,
      x_oob2              = x_oob2,
      y_oob2              = y_oob2,
      data_x              = xscale_w,
      data_y              = yscale,
      xScaleFactor        = ww$scale[[w]],
      yScaleFactor        = NULL,
      x_scale_type        = if (!is.null(track$scale_x)) track$scale_x$type else "genomic",
      y_scale_type        = y_scale_type,
      y_levels            = if (!is.null(track$scale_y) && identical(track$scale_y$type, "discrete"))
                              track$scale_y$levels else NULL,
      y_labels            = if (!is.null(track$scale_y) && identical(track$scale_y$type, "discrete"))
                              (track$scale_y$labels %||% track$scale_y$levels) else NULL,
      y_is_genomic        = isTRUE(track$uses_genomic_y),
      y_sub_panels        = y_sub_panels_default,
      # combine_windows metadata: virtualization maps for the x and y
      # axes, carried so axis drawing can render per-original-window
      # labels and titles, and elements/wrappers can virtualize
      # companion coordinates (e.g. seq_tile data2 j-positions).
      virtual_map_x       = vmap_x,
      virtual_map_y       = vmap_y,
      # Axis flip flags propagated from the track. Elements that
      # respect them mirror their rendered npc coordinates around
      # the panel inner centre; axis drawing uses them to flip
      # tick / label positions so labels follow the data.
      flip_x              = isTRUE(track$flip_x),
      flip_y              = isTRUE(track$flip_y)
    )
    x_cur <- px1 + wgap
  }

  # For multi-window tracks with genomic x scales, compute a shared break
  # step from the widest window so all windows show the same tick interval.
  if (nWin > 1L) {
    xscales <- lapply(panels, `[[`, "xscale")
    valid   <- Filter(function(s) !is.null(s) && length(s) == 2L && diff(s) > 0, xscales)
    if (length(valid) >= 2L) {
      widths_geo   <- sapply(valid, diff)
      widest_scale <- valid[[which.max(widths_geo)]]
      ref_breaks   <- pretty(widest_scale, n = 5L)
      if (length(ref_breaks) >= 2L) {
        step   <- diff(ref_breaks)[1L]
        panels <- lapply(panels, function(p) { p$x_break_step <- step; p })
      }
    }
  }

  # When the scale has an explicit unit, override the auto-inferred scale
  # factors so axis label formatting uses the requested unit.
  if (!is.null(track$scale_x) && inherits(track$scale_x, "SeqScaleGenomic") &&
      !is.null(track$scale_x$unit)) {
    unit_sf <- switch(track$scale_x$unit, Mb = 1e-6, Kb = 1e-3, bp = 1)
    for (k in seq_along(panels)) panels[[k]]$xScaleFactor <- unit_sf
  }
  if (!is.null(track$scale_y) && inherits(track$scale_y, "SeqScaleGenomic") &&
      !is.null(track$scale_y$unit)) {
    unit_sf <- switch(track$scale_y$unit, Mb = 1e-6, Kb = 1e-3, bp = 1)
    for (k in seq_along(panels)) panels[[k]]$yScaleFactor <- unit_sf
  }

  panels
}

#' Compute the bounding box of a track from its panel list
#'
#' @param panels Output of `.build_track_panels()`.
#' @return A list with `x0`, `x1`, `y0`, `y1`.
#' @keywords internal
.panels_to_track_bounds <- function(panels) {
  # Prefer the track-wide `track_full` field (the raw allocated cell,
  # pre-margin). Fall back to aggregating per-window `full` boxes.
  if (length(panels) && !is.null(panels[[1]]$track_full)) {
    return(panels[[1]]$track_full)
  }
  x0s <- vapply(panels, function(p) p$full$x0, numeric(1))
  x1s <- vapply(panels, function(p) p$full$x1, numeric(1))
  y0s <- vapply(panels, function(p) p$full$y0, numeric(1))
  y1s <- vapply(panels, function(p) p$full$y1, numeric(1))
  list(x0 = min(x0s), x1 = max(x1s),
       y0 = min(y0s), y1 = max(y1s))
}

# ── .build_positional_layout ────────────────────────────────────────────────

#' Build layout bounds for the row/column positional mode
#'
#' Iterates `rows` (a list of lists of `SeqTrackR6`), assigning each row a
#' y-band and each track within a row its proportional x-band. Within each
#' track, the genomic windows are split into per-window panels.
#'
#' @param rows A list of lists of `SeqTrackR6` objects.
#' @param aesthetics Merged plot-level aesthetics (already containing
#'   `margins`, `trackGaps`, `windowGaps`).
#' @return A list with `panelBounds` (integer-indexed) and `trackBounds`.
#' @keywords internal
.build_positional_layout <- function(rows, aesthetics) {
  margins    <- aesthetics$margins %||% list(top = 0, right = 0,
                                             bottom = 0, left = 0)
  trackGap   <- aesthetics$track_gaps  %||% 0.01
  windowGap  <- aesthetics[["window.gap.width"]] %||%
                aesthetics$window_gaps           %||% 0.01

  rows <- Filter(function(r) length(r) > 0, rows)
  nRows <- length(rows)
  if (nRows == 0L) {
    return(list(panelBounds = list(), trackBounds = list()))
  }

  availWidth  <- 1 - (margins$left + margins$right)
  availHeight <- 1 - (margins$top  + margins$bottom)

  row_heights <- vapply(rows, function(r) {
    max(vapply(r, function(t) t$track_height %||% 1, numeric(1)))
  }, numeric(1))
  totalRowGap     <- if (nRows > 1) trackGap * (nRows - 1) else 0
  availRowHeight  <- availHeight - totalRowGap
  rel_row_heights <- row_heights / sum(row_heights)
  rowHeights_npc  <- rel_row_heights * availRowHeight

  panelBounds <- list()
  trackBounds <- list()
  track_count <- 0L
  cursor_y <- 1 - margins$top

  for (r in seq_len(nRows)) {
    if (r > 1L) cursor_y <- cursor_y - trackGap
    y_top <- cursor_y
    cursor_y <- cursor_y - rowHeights_npc[r]
    y_bot <- cursor_y

    tracks_in_row <- rows[[r]]
    nT <- length(tracks_in_row)

    track_widths <- vapply(tracks_in_row,
                           function(t) t$track_width %||% 1, numeric(1))
    rel_track_widths <- track_widths / sum(track_widths)
    nGaps <- nT - 1
    totalGapW <- if (nGaps > 0) trackGap * nGaps else 0
    availTrackW <- availWidth - totalGapW
    trackWidths_npc <- rel_track_widths * availTrackW

    cursor_x <- margins$left
    for (ti in seq_len(nT)) {
      if (ti > 1L) cursor_x <- cursor_x + trackGap
      x_left  <- cursor_x
      x_right <- cursor_x + trackWidths_npc[ti]
      cursor_x <- x_right

      track <- tracks_in_row[[ti]]
      track_count <- track_count + 1L
      panels <- .build_track_panels(
        track, x_left, x_right, y_bot, y_top,
        window_gap = windowGap, track_key = track_count
      )
      panelBounds[[track_count]] <- panels
      trackBounds[[track_count]] <- .panels_to_track_bounds(panels)
    }
  }

  list(panelBounds = panelBounds, trackBounds = trackBounds)
}

# ── .build_patchwork_layout ─────────────────────────────────────────────────

#' Build layout bounds for the patchwork-string mode
#'
#' Looks up each track's `track_id` in the parsed layout regions; tracks not
#' present in the layout string are silently skipped.
#'
#' @param tracks Flat list of `SeqTrackR6` objects.
#' @param layout_str The raw layout string.
#' @param aesthetics Merged plot-level aesthetics.
#' @return A list with `panelBounds` and `trackBounds`, each a *named* list
#'   keyed by `track_id`.
#' @keywords internal
.build_patchwork_layout <- function(tracks, layout_str, aesthetics) {
  parsed     <- .parse_layout_string(layout_str)
  margins    <- aesthetics$margins %||% list(top = 0, right = 0,
                                             bottom = 0, left = 0)
  windowGap  <- aesthetics[["window.gap.width"]] %||%
                aesthetics$window_gaps           %||% 0.01

  availWidth  <- 1 - (margins$left + margins$right)
  availHeight <- 1 - (margins$top  + margins$bottom)
  ncol_       <- parsed$ncol
  nrow_       <- parsed$nrow

  panelBounds <- list()
  trackBounds <- list()

  for (track in tracks) {
    tid <- track$track_id
    if (is.null(tid) || is.na(tid) || !tid %in% names(parsed$regions)) next
    region <- parsed$regions[[tid]]

    x0_track <- margins$left + (region$c0 - 1) / ncol_ * availWidth
    x1_track <- margins$left +  region$c1     / ncol_ * availWidth
    y_top_canvas <- 1 - margins$top
    y0_track <- y_top_canvas -  region$r1     / nrow_ * availHeight
    y1_track <- y_top_canvas - (region$r0 - 1) / nrow_ * availHeight

    panels <- .build_track_panels(
      track, x0_track, x1_track, y0_track, y1_track,
      window_gap = windowGap, track_key = tid
    )
    panelBounds[[tid]] <- panels
    trackBounds[[tid]] <- .panels_to_track_bounds(panels)
  }

  list(panelBounds = panelBounds, trackBounds = trackBounds)
}
