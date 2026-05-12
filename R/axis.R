# ── Panel per-element view ───────────────────────────────────────────────────

#' Remap a per-window panel list so xscale/yscale reflect the element's
#' primary/secondary axis selection
#'
#' When an element sets `map(axis.x = 2)` or `map(axis.y = 2)`, its
#' coordinate transforms must use `xscale2` / `yscale2` instead of the
#' primary scales. Every primitive's `prep()` reads `panel$xscale` /
#' `panel$yscale` directly, so this helper builds a per-element shallow
#' copy of the panel list with those fields swapped.
#'
#' @param layout_track List of per-window panel metadata.
#' @param resolved Element `resolved` list (contains `axis_x`/`axis_y`).
#' @return A layout_track with `xscale` / `yscale` rewritten if needed.
#' @keywords internal
.panels_for_element <- function(layout_track, resolved) {
  if (is.null(layout_track)) return(layout_track)
  ax <- resolved$axis_x %||% 1L
  ay <- resolved$axis_y %||% 1L
  if (ax == 1L && ay == 1L) return(layout_track)
  lapply(layout_track, function(panel) {
    if (ax == 2L && !is.null(panel$xscale2)) panel$xscale <- panel$xscale2
    if (ay == 2L && !is.null(panel$yscale2)) panel$yscale <- panel$yscale2
    panel
  })
}

# ── Axis & chrome drawing ────────────────────────────────────────────────────
#
# These helpers render track chrome (backgrounds / borders) and up to four
# axes per track (x1, x2, y1, y2). Each axis is laid out in one of two
# "slots" on its side; when x1 and x2 share a position (e.g. both "top"),
# x1 sits in the slot adjacent to the panel and x2 sits outside. All
# styling is driven by the track's `resolved_theme` (built in
# SeqPlotR6$layoutGrid()).

#' Map a data value to a canvas npc coordinate inside a panel range
#'
#' @param val Numeric values to map.
#' @param scale_lim Length-2 data-range (from which `val` comes).
#' @param npc_lo,npc_hi Canvas npc endpoints of the panel span.
#' @param flip Logical. When `TRUE`, mirror the npc output around the
#'   midpoint of `[npc_lo, npc_hi]` — i.e. low data values render at
#'   `npc_hi` and high values at `npc_lo`.
#' @return Numeric vector of canvas npc coordinates.
#' @keywords internal
.axis_map_npc <- function(val, scale_lim, npc_lo, npc_hi, flip = FALSE) {
  span <- diff(scale_lim)
  if (!is.finite(span) || span == 0) return(rep(npc_lo, length(val)))
  npc <- npc_lo + (val - scale_lim[1]) / span * (npc_hi - npc_lo)
  if (isTRUE(flip)) npc <- (npc_lo + npc_hi) - npc
  npc
}

#' Format a numeric tick label
#'
#' Uses decimal-comma-grouped formatting, no scientific notation, and
#' trims whitespace. A separate helper so the axis draw code reads clean.
#'
#' @param x Numeric vector of tick values.
#' @return Character vector of formatted labels.
#' @keywords internal
.axis_fmt <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

# ── Axis band geometry ───────────────────────────────────────────────────────

#' Layout the two axis slots available on each side of a track
#'
#' For each side (top / bottom / left / right) there are up to two axes
#' (x1 & x2 or y1 & y2) that could live there. The slot nearer the panel
#' is slot 1; the outer slot is slot 2. When both axes with the same
#' dimension share a side, x1 takes slot 1 and x2 takes slot 2; otherwise
#' the single axis takes slot 1 (the full band).
#'
#' @param resolved_theme The track's resolved theme.
#' @param first_panel The track's first panel (source of band bounds).
#' @return A named list with one entry per side (`bottom`, `top`,
#'   `left`, `right`). Each entry is a list of up to two slots; each
#'   slot has `anchor` (the npc coordinate of the panel-adjacent edge),
#'   `outer` (the outer edge of the slot), and `band_outer` (the outer
#'   edge of the whole track_inner margin band).
#' @keywords internal
.compute_axis_bands <- function(resolved_theme, first_panel) {
  inner   <- first_panel$inner
  tr_in   <- first_panel$track_inner  # container for axes + panel
  win_im  <- first_panel$window_inner_margin %||%
    list(top = 0, right = 0, bottom = 0, left = 0)

  # The band for each side is the window_inner_margin between the panel
  # edge and the track_inner edge.
  sides <- list(
    bottom = list(anchor = inner$y0, band_outer = inner$y0 - win_im$bottom,
                  axis = "x"),
    top    = list(anchor = inner$y1, band_outer = inner$y1 + win_im$top,
                  axis = "x"),
    left   = list(anchor = inner$x0, band_outer = inner$x0 - win_im$left,
                  axis = "y"),
    right  = list(anchor = inner$x1, band_outer = inner$x1 + win_im$right,
                  axis = "y")
  )

  # Assign x1/x2/y1/y2 to a side based on their resolved position.
  per_side <- list(bottom = character(0), top = character(0),
                   left   = character(0), right = character(0))
  for (side in c("x1","x2","y1","y2")) {
    pos <- resolved_theme$axes[[side]]$position
    if (is.null(pos)) next
    if (!(pos %in% names(per_side))) next
    # Only x axes can sit top/bottom; only y axes can sit left/right.
    is_x <- substr(side, 1L, 1L) == "x"
    is_h <- pos %in% c("top", "bottom")
    if (is_x != is_h) next
    per_side[[pos]] <- c(per_side[[pos]], side)
  }
  # Stable order: x1 before x2, y1 before y2 — so x1 always sits in
  # slot 1 (adjacent to panel) when they share a position.
  for (side in names(per_side)) {
    per_side[[side]] <- sort(per_side[[side]])
  }

  # Turn per-side axis lists into concrete slots.
  slots <- list()
  for (side in names(sides)) {
    axes <- per_side[[side]]
    if (length(axes) == 0L) { slots[[side]] <- list(); next }
    anchor     <- sides[[side]]$anchor
    band_outer <- sides[[side]]$band_outer
    total_span <- band_outer - anchor
    n <- length(axes)
    per <- total_span / n
    out <- vector("list", n)
    for (i in seq_along(axes)) {
      out[[i]] <- list(
        side       = side,
        axis       = axes[i],
        slot_index = i,
        anchor     = anchor + (i - 1L) * per,
        outer      = anchor + i * per,
        band_outer = band_outer
      )
    }
    names(out) <- axes
    slots[[side]] <- out
  }
  slots
}

# ── Active-axis rule ─────────────────────────────────────────────────────────

#' Decide whether a given axis should render
#'
#' @param track The `SeqTrackR6`.
#' @param resolved_theme Track resolved theme.
#' @param side One of `"x1"`, `"x2"`, `"y1"`, `"y2"`.
#' @return Logical.
#' @keywords internal
.axis_active <- function(track, resolved_theme, side) {
  spec <- resolved_theme$axes[[side]]
  if (!is.null(spec$visible) && !is.na(spec$visible))
    return(isTRUE(spec$visible))
  # Default: primary axes visible; secondary only when active.
  if (side == "x1") {
    return(!is.null(track$scale_x))
  } else if (side == "y1") {
    return(!is.null(track$scale_y))
  } else if (side == "x2") {
    return(isTRUE(track$has_axis_x2) && !is.null(track$scale_x2))
  } else if (side == "y2") {
    return(isTRUE(track$has_axis_y2) && !is.null(track$scale_y2))
  }
  FALSE
}

# ── combine_windows axis support ────────────────────────────────────────────

#' Expand a combined virtual panel into per-original-window sub-panels
#'
#' When a track was laid out with `combine_windows = TRUE`, the layout
#' returns one panel whose xscale is virtual `c(1, virtual_total)`.
#' For axis drawing we want each original window to render its own
#' ticks, labels, and title at the corresponding virtual sub-range —
#' so produce a list of pseudo-panels covering each `vmap_x` entry,
#' with `xscale` set to the original genomic range and `inner.x0/x1`
#' narrowed to the virtual sub-range's npc extent.
#'
#' @param panels The track's panel list as returned by
#'   `.build_track_panels()`.
#' @return Either the input `panels` unchanged (when no virtual map is
#'   present) or a list of virtual sub-panels suitable for x-axis
#'   drawing.
#' @keywords internal
.expand_panels_for_combined_axis <- function(panels) {
  if (length(panels) == 0L) return(panels)
  vmap <- panels[[1]]$virtual_map_x
  if (is.null(vmap)) return(panels)
  if (length(vmap$seqnames) <= 1L) return(panels)

  base <- panels[[1]]
  # Map a virtual-x value to the panel's inner-x npc.
  to_npc_x <- function(v) {
    span <- diff(base$xscale)
    if (!is.finite(span) || span == 0) return(base$inner$x0)
    base$inner$x0 + (v - base$xscale[1]) / span *
      (base$inner$x1 - base$inner$x0)
  }
  # Default scale factor: Mb when scale type is genomic, else 1.
  default_sf <- if (identical(base$x_scale_type, "genomic")) 1e-6 else 1

  lapply(seq_along(vmap$seqnames), function(k) {
    win <- base
    win$xscale  <- c(vmap$genomic_start[k], vmap$genomic_end[k])
    win$data_x  <- win$xscale
    win$inner$x0 <- to_npc_x(vmap$virtual_start[k])
    win$inner$x1 <- to_npc_x(vmap$virtual_end[k])
    win$xScaleFactor <- default_sf
    win$.combined_window_label <- paste0(
      vmap$seqnames[k], ": ",
      format(vmap$genomic_start[k], big.mark = ",", scientific = FALSE), "-",
      format(vmap$genomic_end[k],   big.mark = ",", scientific = FALSE))
    win
  })
}

#' Expand a multi-y-window panel into one pseudo-panel per y sub-panel
#'
#' Mirror of [.expand_panels_for_combined_axis] but for the y-axis.
#' Two modes:
#' \itemize{
#'   \item `y_sub_panels` (vertically stacked sub-panels): each entry
#'     gets its own npc band with the corresponding window's yscale.
#'   \item `virtual_map_y` (combined virtual y axis): the single y
#'     band is split into per-original-window sub-ranges along the
#'     virtual axis, each carrying the original genomic yscale for
#'     label rendering.
#' }
#'
#' @param panels The track's panel list.
#' @return A list of pseudo-panels (one per y sub-panel); or the input
#'   panels unchanged when no multi-y data is present.
#' @keywords internal
.expand_panels_for_multi_y_axis <- function(panels) {
  if (length(panels) == 0L) return(panels)
  base <- panels[[1]]
  default_sf <- if (identical(base$y_scale_type, "genomic")) 1e-6 else 1

  y_subs <- base$y_sub_panels
  if (!is.null(y_subs) && length(y_subs) > 0L) {
    return(lapply(seq_along(y_subs), function(k) {
      sub <- y_subs[[k]]
      p   <- base
      p$yscale <- sub$yscale
      p$data_y <- sub$yscale
      p$inner$y0 <- sub$y0
      p$inner$y1 <- sub$y1
      p$yScaleFactor <- default_sf
      p$.combined_window_label <- paste0(
        sub$seqname, ": ",
        format(sub$yscale[1], big.mark = ",", scientific = FALSE), "-",
        format(sub$yscale[2], big.mark = ",", scientific = FALSE))
      p
    }))
  }

  vmap <- base$virtual_map_y
  if (!is.null(vmap) && length(vmap$seqnames) > 1L) {
    span <- diff(base$yscale)
    if (!is.finite(span) || span == 0) return(panels)
    to_npc_y <- function(v)
      base$inner$y0 + (v - base$yscale[1]) / span *
        (base$inner$y1 - base$inner$y0)
    return(lapply(seq_along(vmap$seqnames), function(k) {
      p <- base
      p$yscale <- c(vmap$genomic_start[k], vmap$genomic_end[k])
      p$data_y <- p$yscale
      p$inner$y0 <- to_npc_y(vmap$virtual_start[k])
      p$inner$y1 <- to_npc_y(vmap$virtual_end[k])
      p$yScaleFactor <- default_sf
      p$.combined_window_label <- paste0(
        vmap$seqnames[k], ": ",
        format(vmap$genomic_start[k], big.mark = ",", scientific = FALSE), "-",
        format(vmap$genomic_end[k],   big.mark = ",", scientific = FALSE))
      p
    }))
  }

  panels
}

# ── One-axis draw ────────────────────────────────────────────────────────────

#' Draw a single axis (line, ticks, minor ticks, labels, title)
#'
#' @param side One of `"x1"`, `"x2"`, `"y1"`, `"y2"`.
#' @param track The `SeqTrackR6`.
#' @param panels The track's panel list.
#' @param slot The slot from `.compute_axis_bands()` for this axis.
#' @return Invisibly `NULL`.
#' @keywords internal
.draw_one_axis <- function(side, track, panels, slot) {
  spec <- track$resolved_theme$axes[[side]]
  is_x <- substr(side, 1L, 1L) == "x"
  sec  <- substr(side, 2L, 2L) == "2"

  scale  <- if (is_x) (if (sec) track$scale_x2 else track$scale_x)
            else     (if (sec) track$scale_y2 else track$scale_y)
  if (is.null(scale)) return(invisible())

  # combine_windows: when the track was rendered as a single virtual
  # panel, expand it into per-original-window sub-panels for axis
  # drawing so each window gets its own ticks, labels, and title.
  if (is_x && !sec) {
    panels <- .expand_panels_for_combined_axis(panels)
  }
  # Multi-y-window: when a panel carries y_sub_panels, expand into
  # per-y-window pseudo panels so y-axis labels and titles draw per
  # sub-panel.
  if (!is_x && !sec) {
    panels <- .expand_panels_for_multi_y_axis(panels)
  }

  # For x-axes, iterate per-window; for y-axes, draw once per track
  # unless the theme asks for per-window or the panels have been
  # expanded for multi-y (each pseudo-panel = one y sub-window).
  per_window <- is_x ||
                isTRUE(track$resolved_theme$y_per_window) ||
                (length(panels) > 1L &&
                 !is.null(panels[[1]]$.combined_window_label))

  tick_dir <- if (side == "x1")        -1
              else if (side == "x2")    1
              else if (side == "y1")   -1
              else                      1

  for (w in seq_along(panels)) {
    if (!per_window && w > 1L) next
    win    <- panels[[w]]
    p      <- win$inner
    scale_lim <- if (is_x) (if (sec) win$xscale2 else win$xscale)
                 else      (if (sec) win$yscale2 else win$yscale)
    if (is.null(scale_lim)) next
    # When the panel carries a track-level break step (set by layout for
    # multi-window tracks), force all windows to use the same tick interval.
    # Only applies to genomic x axes where no explicit breaks are set.
    if (is_x && !sec &&
        !is.null(win$x_break_step) &&
        is.null(scale$breaks) &&
        identical(scale$type, "genomic")) {
      step     <- win$x_break_step
      br_start <- ceiling(scale_lim[1] / step) * step
      br_end   <- floor(scale_lim[2]   / step) * step
      forced   <- if (br_start <= br_end) seq(br_start, br_end, by = step) else numeric(0)
      # Prepend 0 when the window starts below the first step — keeps 0 as
      # the natural axis origin; scale.R will retain it via its own snap.
      if (scale_lim[1] > 0 && scale_lim[1] < step && length(forced) > 0)
        forced <- c(0, forced)
      sc_tmp        <- scale
      sc_tmp$breaks <- forced
      meta <- .compute_scale_breaks(sc_tmp, scale_lim, plot_range = scale_lim)
    } else {
      meta <- .compute_scale_breaks(scale, scale_lim, plot_range = scale_lim)
    }
    if (length(meta$breaks) == 0L) next

    # Panel npc endpoints along the axis direction.
    if (is_x) { lo <- p$x0; hi <- p$x1 } else { lo <- p$y0; hi <- p$y1 }

    # Flip flag for this axis. Honour the panel's flip_x / flip_y so
    # tick label positions follow whatever the elements rendered.
    flip <- if (is_x) isTRUE(win$flip_x) else isTRUE(win$flip_y)

    # Anchor is the npc coord on the opposite dimension.
    anchor <- slot$anchor
    tick_len <- spec$ticks$length %||% 0.005

    # --- Axis line ---
    if (isTRUE(spec$line$visible) && !is.null(meta$axis_range)) {
      ar <- meta$axis_range
      ar_lo <- .axis_map_npc(ar[1], scale_lim, lo, hi, flip = flip)
      ar_hi <- .axis_map_npc(ar[2], scale_lim, lo, hi, flip = flip)
      gp_line <- grid::gpar(col = spec$line$col %||% "#1C1B1A",
                            lwd = spec$line$lwd %||% 1,
                            alpha = spec$line$alpha %||% 1)
      if (is_x) {
        grid::grid.lines(x = grid::unit(c(ar_lo, ar_hi), "npc"),
                         y = grid::unit(c(anchor, anchor), "npc"),
                         gp = gp_line)
      } else {
        grid::grid.lines(x = grid::unit(c(anchor, anchor), "npc"),
                         y = grid::unit(c(ar_lo, ar_hi), "npc"),
                         gp = gp_line)
      }
    }

    # --- Range-style label (early exit for ticks/labels block) ---
    label_style <- spec$text$style    %||% "tick"
    label_pos   <- spec$text$position %||% "axis"
    if (isTRUE(spec$text$visible) && identical(label_style, "range")) {
      sf <- if (is_x) {
        if (identical(scale$type, "genomic")) (win$xScaleFactor %||% 1) else 1
      } else {
        if (identical(scale$type, "genomic")) (win$yScaleFactor %||% 1e-6) else 1
      }
      lo_val <- scale_lim[1] * sf
      hi_val <- scale_lim[2] * sf
      range_label <- sprintf("[%s\u2013%s]",
                             .axis_fmt(lo_val), .axis_fmt(hi_val))

      offset <- spec$text$offset %||% (if (is_x) 0.015 else 0.010)

      if (is.numeric(label_pos) && length(label_pos) == 2L) {
        # Panel-relative NPC -> canvas NPC.
        lx  <- p$x0 + label_pos[1] * (p$x1 - p$x0)
        ly  <- p$y0 + label_pos[2] * (p$y1 - p$y0)
        rot <- 0
      } else {
        if (is_x) {
          lx <- (p$x0 + p$x1) / 2
          ly <- anchor + tick_dir * offset
          rot <- 0
        } else {
          lx <- anchor + tick_dir * offset
          ly <- (p$y0 + p$y1) / 2
          rot <- 90
        }
      }

      grid::grid.text(
        label = range_label,
        x     = grid::unit(lx, "npc"),
        y     = grid::unit(ly, "npc"),
        rot   = rot,
        just  = c(spec$text$hjust %||% "center",
                  spec$text$vjust %||% "center"),
        gp    = grid::gpar(col = spec$text$col %||% "#1C1B1A",
                            cex = spec$text$size %||% 0.6)
      )

      # Still draw ticks (just not the per-tick labels) when ticks visible.
      if (isTRUE(spec$ticks$visible)) {
        gp_tick <- grid::gpar(col = spec$ticks$col %||% spec$line$col,
                              lwd = spec$ticks$lwd %||% spec$line$lwd,
                              alpha = spec$ticks$alpha %||% 1)
        for (i in seq_along(meta$breaks)) {
          b   <- meta$breaks[i]
          pos <- .axis_map_npc(b, scale_lim, lo, hi, flip = flip)
          if (is_x) {
            grid::grid.lines(x = grid::unit(c(pos, pos), "npc"),
                             y = grid::unit(c(anchor,
                                              anchor + tick_dir * tick_len),
                                            "npc"),
                             gp = gp_tick)
          } else {
            grid::grid.lines(x = grid::unit(c(anchor,
                                              anchor + tick_dir * tick_len),
                                            "npc"),
                             y = grid::unit(c(pos, pos), "npc"),
                             gp = gp_tick)
          }
        }
      }
      next
    }

    # --- Major ticks & labels ---
    if (isTRUE(spec$ticks$visible) || isTRUE(spec$text$visible)) {
      gp_tick <- grid::gpar(col = spec$ticks$col %||% spec$line$col,
                            lwd = spec$ticks$lwd %||% spec$line$lwd,
                            alpha = spec$ticks$alpha %||% 1)
      for (i in seq_along(meta$breaks)) {
        b   <- meta$breaks[i]
        pos <- .axis_map_npc(b, scale_lim, lo, hi, flip = flip)
        if (isTRUE(spec$ticks$visible)) {
          if (is_x) {
            grid::grid.lines(x = grid::unit(c(pos, pos), "npc"),
                             y = grid::unit(c(anchor,
                                              anchor + tick_dir * tick_len),
                                            "npc"),
                             gp = gp_tick)
          } else {
            grid::grid.lines(x = grid::unit(c(anchor,
                                              anchor + tick_dir * tick_len),
                                            "npc"),
                             y = grid::unit(c(pos, pos), "npc"),
                             gp = gp_tick)
          }
        }

        if (isTRUE(spec$text$visible)) {
          lbl <- meta$labels[i]
          if (is.numeric(lbl) || (is.character(lbl) && !is.na(suppressWarnings(as.numeric(lbl))))) {
            scaleF <- if (is_x) {
              if (identical(scale$type, "genomic")) (win$xScaleFactor %||% 1) else 1
            } else {
              if (identical(scale$type, "genomic")) (win$yScaleFactor %||% 1e-6) else 1
            }
            lbl <- .axis_fmt(as.numeric(lbl) * scaleF)
          } else {
            lbl <- as.character(lbl)
          }
          offset <- spec$text$offset %||% (if (is_x) 0.015 else 0.010)
          gp_text <- grid::gpar(col = spec$text$col %||% "#1C1B1A",
                                cex = spec$text$size %||% 0.6)
          if (is_x) {
            just_v <- if (side == "x1") "top" else "bottom"
            grid::grid.text(
              label = lbl,
              x = grid::unit(pos, "npc"),
              y = grid::unit(anchor + tick_dir * offset, "npc"),
              just = c("center", just_v),
              rot  = spec$text$angle %||% 0,
              gp   = gp_text
            )
          } else {
            just_h <- if (side == "y1") "right" else "left"
            grid::grid.text(
              label = lbl,
              x = grid::unit(anchor + tick_dir * offset, "npc"),
              y = grid::unit(pos, "npc"),
              just = c(just_h, "center"),
              rot  = spec$text$angle %||% 0,
              gp   = gp_text
            )
          }
        }
      }

      # Minor ticks (no labels)
      if (isTRUE(spec$ticks$visible) && length(meta$minor_breaks) > 0L) {
        minor_len <- tick_len * 0.6
        for (mb in meta$minor_breaks) {
          pos <- .axis_map_npc(mb, scale_lim, lo, hi, flip = flip)
          if (is_x) {
            grid::grid.lines(x = grid::unit(c(pos, pos), "npc"),
                             y = grid::unit(c(anchor,
                                              anchor + tick_dir * minor_len),
                                            "npc"),
                             gp = gp_tick)
          } else {
            grid::grid.lines(x = grid::unit(c(anchor,
                                              anchor + tick_dir * minor_len),
                                            "npc"),
                             y = grid::unit(c(pos, pos), "npc"),
                             gp = gp_tick)
          }
        }
      }
    }
  }

  # --- Axis title.
  # Default: one title per track, centered in the track_outer margin.
  # With combine_windows on the x-axis: one title per *original* window,
  # centered over that window's virtual npc sub-range, naming the chrom
  # and genomic span (e.g. "chr14: 98,310,000-99,310,000").
  if (isTRUE(spec$title$visible)) {
    first <- panels[[1]]
    tm    <- first$track_mapping
    default_text <- if (!is.null(tm) && !is.null(tm[[if (is_x) "x" else "y"]]))
      paste(deparse(tm[[if (is_x) "x" else "y"]]), collapse = " ")
    else NULL
    user_title <- spec$title$text
    gp_title <- grid::gpar(col = spec$title$col %||% "#1C1B1A",
                           cex = spec$title$size %||% 0.8)

    # In-panel placement: if title$position is c(x_npc, y_npc), draw the
    # title once at those panel-relative NPC coords (mapped to canvas NPC
    # via the first panel's inner box), bypassing the margin-band logic.
    title_pos <- spec$title$position %||% "axis"
    if (!identical(title_pos, "axis") &&
        is.numeric(title_pos) && length(title_pos) == 2L) {
      title_text <- user_title %||% default_text
      if (!is.null(title_text) && nzchar(title_text)) {
        in_box <- first$inner
        tx <- in_box$x0 + title_pos[1] * (in_box$x1 - in_box$x0)
        ty <- in_box$y0 + title_pos[2] * (in_box$y1 - in_box$y0)
        rot <- if (!is_x) 90 else 0
        grid::grid.text(
          label = title_text,
          x     = grid::unit(tx, "npc"),
          y     = grid::unit(ty, "npc"),
          rot   = rot,
          just  = c(spec$title$hjust %||% "center",
                    spec$title$vjust %||% "center"),
          gp    = gp_title
        )
      }
      return(invisible())
    }

    is_per_window <- length(panels) > 1L &&
                     !is.null(panels[[1]]$.combined_window_label)

    if (is_per_window && is_x) {
      ty <- if (side == "x1")
              (first$track_full$y0 + first$track_outer$y0) / 2
            else
              (first$track_full$y1 + first$track_outer$y1) / 2
      for (w in seq_along(panels)) {
        win    <- panels[[w]]
        title_text <- user_title %||% win$.combined_window_label %||%
                      default_text
        if (is.null(title_text) || !nzchar(title_text)) next
        tx <- (win$inner$x0 + win$inner$x1) / 2
        grid::grid.text(label = title_text,
                        x = grid::unit(tx, "npc"),
                        y = grid::unit(ty, "npc"),
                        gp = gp_title)
      }
    } else if (is_per_window && !is_x) {
      tx <- if (side == "y1")
              (first$track_full$x0 + first$track_outer$x0) / 2
            else
              (first$track_full$x1 + first$track_outer$x1) / 2
      rot <- if (side == "y1") 90 else -90
      for (w in seq_along(panels)) {
        win <- panels[[w]]
        title_text <- user_title %||% win$.combined_window_label %||%
                      default_text
        if (is.null(title_text) || !nzchar(title_text)) next
        ty <- (win$inner$y0 + win$inner$y1) / 2
        grid::grid.text(label = title_text,
                        x = grid::unit(tx, "npc"),
                        y = grid::unit(ty, "npc"),
                        rot = rot,
                        gp = gp_title)
      }
    } else {
      title_text <- user_title %||% default_text
      if (!is.null(title_text) && nzchar(title_text)) {
        if (is_x) {
          tx <- (first$track_outer$x0 + first$track_outer$x1) / 2
          ty <- if (side == "x1")
                  (first$track_full$y0  + first$track_outer$y0) / 2
                else
                  (first$track_full$y1  + first$track_outer$y1) / 2
          grid::grid.text(label = title_text,
                          x = grid::unit(tx, "npc"),
                          y = grid::unit(ty, "npc"),
                          gp = gp_title)
        } else {
          ty <- (first$track_outer$y0 + first$track_outer$y1) / 2
          tx <- if (side == "y1")
                  (first$track_full$x0  + first$track_outer$x0) / 2
                else
                  (first$track_full$x1  + first$track_outer$x1) / 2
          rot <- if (side == "y1") 90 else -90
          grid::grid.text(label = title_text,
                          x = grid::unit(tx, "npc"),
                          y = grid::unit(ty, "npc"),
                          rot = rot,
                          gp = gp_title)
        }
      }
    }
  }

  invisible()
}

# ── Top-level axis dispatcher ────────────────────────────────────────────────

#' Draw all active axes for one track
#'
#' @param track A `SeqTrackR6`.
#' @param panels Panel list (from `layoutGrid()`).
#' @return Invisibly `NULL`.
#' @keywords internal
.draw_track_axes <- function(track, panels) {
  rt <- track$resolved_theme
  if (is.null(rt)) return(invisible())
  first <- panels[[1]]
  bands <- .compute_axis_bands(rt, first)
  for (side in c("x1","x2","y1","y2")) {
    if (!.axis_active(track, rt, side)) next
    pos <- rt$axes[[side]]$position
    slot <- bands[[pos]][[side]]
    if (is.null(slot)) next
    .draw_one_axis(side, track, panels, slot)
  }
  invisible()
}

# ── Track chrome ─────────────────────────────────────────────────────────────

#' Draw the track background / border and per-window panel chrome
#'
#' Reads from `track$resolved_theme$chrome`.
#'
#' @param track The `SeqTrackR6`.
#' @param panels Panel list for this track.
#' @param tb Track bounds list (`x0`, `x1`, `y0`, `y1`).
#' @return Invisibly `NULL`.
#' @keywords internal
.draw_track_chrome <- function(track, panels, tb) {
  rt <- track$resolved_theme
  if (is.null(rt) || is.null(tb)) return(invisible())
  chrome <- rt$chrome

  .rect <- function(box, fill, col, lwd, alpha) {
    if (is.null(box)) return(invisible())
    grid::grid.rect(
      x = grid::unit(box$x0, "npc"),
      y = grid::unit(box$y0, "npc"),
      width  = grid::unit(box$x1 - box$x0, "npc"),
      height = grid::unit(box$y1 - box$y0, "npc"),
      just = c("left", "bottom"),
      gp   = grid::gpar(fill = fill, col = col, lwd = lwd %||% 0.5,
                        alpha = alpha %||% 1)
    )
  }

  # Track background + border — one rect covering the full track cell.
  .rect(tb,
        fill  = chrome$background$fill,
        col   = chrome$border$col,
        lwd   = chrome$border$lwd,
        alpha = chrome$background$alpha)
  if (!is.na(chrome$border$col)) {
    .rect(tb,
          fill  = NA,
          col   = chrome$border$col,
          lwd   = chrome$border$lwd,
          alpha = chrome$border$alpha)
  }

  # Per-window panel background + border — the `inner` box (plot area).
  for (win in panels) {
    p <- win$inner
    .rect(p,
          fill  = chrome$window$background$fill,
          col   = NA,
          lwd   = 0.5,
          alpha = chrome$window$background$alpha)
    if (!is.na(chrome$window$border$col %||% NA)) {
      .rect(p,
            fill  = NA,
            col   = chrome$window$border$col,
            lwd   = chrome$window$border$lwd,
            alpha = chrome$window$border$alpha)
    }
  }

  # combine_windows: draw a thin vertical separator at each interior
  # window boundary (between virtual sub-windows). Helps the eye pick
  # out which side of the panel belongs to which original window.
  vmap <- if (length(panels)) panels[[1]]$virtual_map_x else NULL
  if (!is.null(vmap) && length(vmap$seqnames) > 1L) {
    base <- panels[[1]]
    span <- diff(base$xscale)
    if (is.finite(span) && span > 0) {
      to_npc_x <- function(v)
        base$inner$x0 + (v - base$xscale[1]) / span *
          (base$inner$x1 - base$inner$x0)
      sep_col <- chrome$window$border$col %||% "grey60"
      if (is.na(sep_col)) sep_col <- "grey60"
      sep_lwd <- chrome$window$border$lwd %||% 0.6
      for (k in seq_len(length(vmap$seqnames) - 1L)) {
        x_npc <- to_npc_x(vmap$virtual_end[k])
        grid::grid.lines(
          x = grid::unit(c(x_npc, x_npc), "npc"),
          y = grid::unit(c(base$inner$y0, base$inner$y1), "npc"),
          gp = grid::gpar(col = sep_col, lwd = sep_lwd)
        )
      }
    }
  }

  # Multi-y-window: draw a thin horizontal separator at each interior
  # y-window boundary.
  y_subs <- if (length(panels)) panels[[1]]$y_sub_panels else NULL
  if (!is.null(y_subs) && length(y_subs) > 1L) {
    sep_col <- chrome$window$border$col %||% "grey60"
    if (is.na(sep_col)) sep_col <- "grey60"
    sep_lwd <- chrome$window$border$lwd %||% 0.6
    base <- panels[[1]]
    for (k in seq_len(length(y_subs) - 1L)) {
      y_npc <- (y_subs[[k]]$y1 + y_subs[[k + 1L]]$y0) / 2
      grid::grid.lines(
        x = grid::unit(c(base$inner$x0, base$inner$x1), "npc"),
        y = grid::unit(c(y_npc, y_npc), "npc"),
        gp = grid::gpar(col = sep_col, lwd = sep_lwd)
      )
    }
  }
  invisible()
}
