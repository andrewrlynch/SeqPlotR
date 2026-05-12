# ── Coordinate operations ────────────────────────────────────────────────────
#
# Internal helpers that move data between the three coordinate systems used by
# SeqPlotR:
#
#   data  -- raw genomic / numeric data values
#   npc   -- [0, 1] panel-relative coordinates
#   canvas-- [0, 1] page-relative npc coordinates within a panel's `inner` box
#
# Ports of THEfunc's `clipToXscale()`, `convertDataToGrid()`, `gridToCanvas()`.

#' Clip genomic ranges to an x-axis scale
#'
#' Returns the subset of paired ranges that intersect the closed interval
#' `[xscale[1], xscale[2]]`, with each range trimmed to fit inside that
#' interval. Out-of-bounds ranges are silently dropped — no warning.
#'
#' @param x0 Numeric vector of range starts.
#' @param x1 Numeric vector of range ends (same length as `x0`).
#' @param xscale Numeric vector of length 2: `c(min, max)`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`x0`}{clipped starts}
#'     \item{`x1`}{clipped ends}
#'     \item{`mask`}{logical vector aligned with the input that marks which
#'       rows survived clipping (used to subset companion vectors).}
#'   }
#'
#' @keywords internal
.clip_to_windows <- function(x0, x1, xscale) {
  keep <- x1 >= xscale[1] & x0 <= xscale[2]
  list(
    x0   = pmax(x0[keep], xscale[1]),
    x1   = pmin(x1[keep], xscale[2]),
    mask = keep
  )
}

#' Normalise a data x value to the unit interval within a scale range
#'
#' @param x Numeric vector of data x values.
#' @param xscale Numeric vector of length 2: `c(min, max)`.
#' @return Numeric vector of the same length as `x`, clamped to `[0, 1]`.
#' @keywords internal
.data_to_npc <- function(x, xscale) {
  span <- diff(xscale)
  if (!is.finite(span) || span == 0) return(rep(0, length(x)))
  pmax(0, pmin(1, (x - xscale[1]) / span))
}

#' Normalise a data y value to the unit interval within a scale range
#'
#' @param y Numeric vector of data y values.
#' @param yscale Numeric vector of length 2: `c(min, max)`.
#' @return Numeric vector of the same length as `y`, clamped to `[0, 1]`.
#' @keywords internal
.data_to_npc_y <- function(y, yscale) {
  span <- diff(yscale)
  if (!is.finite(span) || span == 0) return(rep(0, length(y)))
  pmax(0, pmin(1, (y - yscale[1]) / span))
}

#' Map unit-interval panel coordinates to canvas npc coordinates
#'
#' Uses `panel_meta$inner` as the destination rectangle.
#'
#' @param u Numeric vector of x in `[0, 1]` panel coordinates.
#' @param v Numeric vector of y in `[0, 1]` panel coordinates.
#' @param panel_meta A panel metadata list with element `inner` containing
#'   `x0`, `x1`, `y0`, `y1`.
#' @return A list with `x` and `y` numeric vectors of canvas npc coordinates.
#' @keywords internal
.npc_to_canvas <- function(u, v, panel_meta) {
  x <- panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0)
  y <- panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)
  list(x = x, y = y)
}

# ── combine_windows helpers ──────────────────────────────────────────────────

#' Build a virtualization map for a multi-region GRanges
#'
#' Concatenates the windows into a single virtual coordinate system in
#' insertion order. Each window contributes `width(window)` to the
#' virtual axis, with optional gaps (`virtual_gap`) inserted between
#' windows so a separator can be drawn at the boundaries.
#'
#' @param windows A `GRanges` of one or more windows.
#' @param virtual_gap Numeric (bp) gap inserted between consecutive
#'   windows in virtual coordinates. Default `0`.
#' @return A list with components:
#' \describe{
#'   \item{`seqnames`}{character vector of original seqnames, length
#'     `length(windows)`.}
#'   \item{`genomic_start`, `genomic_end`}{numeric vectors of original
#'     genomic ranges.}
#'   \item{`virtual_start`, `virtual_end`}{numeric vectors of the
#'     corresponding virtual ranges.}
#'   \item{`virtual_total`}{the total virtual extent (right edge of the
#'     last window).}
#'   \item{`virtual_gap`}{the gap size used (echoed back).}
#'   \item{`combined_window`}{a single-range `GRanges` covering the
#'     full virtual extent, used in place of `windows` for layout.}
#' }
#' @keywords internal
.build_virtual_map <- function(windows, virtual_gap = 0) {
  n <- length(windows)
  seqs   <- as.character(GenomicRanges::seqnames(windows))
  gstart <- BiocGenerics::start(windows)
  gend   <- BiocGenerics::end(windows)
  widths <- BiocGenerics::width(windows)

  vstart <- numeric(n)
  vend   <- numeric(n)
  cur    <- 1
  for (k in seq_len(n)) {
    vstart[k] <- cur
    vend[k]   <- cur + widths[k] - 1
    cur       <- vend[k] + 1 + virtual_gap
  }
  total <- if (n == 0L) 0 else vend[n]

  combined <- GenomicRanges::GRanges(
    "__combined__",
    IRanges::IRanges(start = 1, end = max(total, 1))
  )

  list(
    seqnames        = seqs,
    genomic_start   = as.numeric(gstart),
    genomic_end     = as.numeric(gend),
    virtual_start   = vstart,
    virtual_end     = vend,
    virtual_total   = total,
    virtual_gap     = virtual_gap,
    combined_window = combined
  )
}

#' Convert genomic (seqname, position) pairs to virtual positions
#'
#' Uses the map returned by [.build_virtual_map()]. Each position is
#' resolved against the window whose seqname matches and whose genomic
#' range contains the position. Positions that match no window get
#' `NA_real_`.
#'
#' @param seqnames Character vector of seqnames.
#' @param positions Numeric vector of genomic positions, same length.
#' @param vmap A virtualization map.
#' @return Numeric vector of virtual positions, `NA` for unmatched.
#' @keywords internal
.virtualize_positions <- function(seqnames, positions, vmap) {
  out <- rep(NA_real_, length(positions))
  for (k in seq_along(vmap$seqnames)) {
    mask <- seqnames == vmap$seqnames[k] &
            positions >= vmap$genomic_start[k] &
            positions <= vmap$genomic_end[k]
    out[mask] <- vmap$virtual_start[k] +
                 (positions[mask] - vmap$genomic_start[k])
  }
  out
}

#' Virtualize a `GRanges` according to a virtualization map
#'
#' Returns a new `GRanges` with seqnames replaced by the combined
#' sentinel, and ranges replaced by virtual coordinates. Rows whose
#' original (seqname, position) pair does not fall in any window are
#' silently dropped. Optional companion mcols giving a *secondary*
#' position pair (e.g. `j_start`/`j_end` with an optional `j_chrom`
#' column for Hi-C) are also virtualized in place.
#'
#' @param gr A `GRanges`.
#' @param vmap A virtualization map.
#' @param j_start_col,j_end_col,j_chrom_col Optional mcols column names
#'   for a secondary position pair (with optional companion seqnames).
#'   When `j_chrom_col` is `NULL`, the j seqname is taken to be the
#'   same as the primary GRanges seqname.
#' @return A virtualized `GRanges`, possibly shorter than `gr`.
#' @keywords internal
.virtualize_granges <- function(gr, vmap,
                                j_start_col = NULL, j_end_col = NULL,
                                j_chrom_col = NULL) {
  if (length(gr) == 0L) return(gr)
  i_seq   <- as.character(GenomicRanges::seqnames(gr))
  i_start <- BiocGenerics::start(gr)
  i_end   <- BiocGenerics::end(gr)

  i_vstart <- .virtualize_positions(i_seq, i_start, vmap)
  i_vend   <- .virtualize_positions(i_seq, i_end,   vmap)

  keep <- !is.na(i_vstart) & !is.na(i_vend)

  if (!is.null(j_start_col) && !is.null(j_end_col) &&
      j_start_col %in% names(S4Vectors::mcols(gr)) &&
      j_end_col   %in% names(S4Vectors::mcols(gr))) {
    j_seq <- if (!is.null(j_chrom_col) &&
                 j_chrom_col %in% names(S4Vectors::mcols(gr)))
               as.character(S4Vectors::mcols(gr)[[j_chrom_col]])
             else
               i_seq
    j_start <- as.numeric(S4Vectors::mcols(gr)[[j_start_col]])
    j_end   <- as.numeric(S4Vectors::mcols(gr)[[j_end_col]])
    j_vstart <- .virtualize_positions(j_seq, j_start, vmap)
    j_vend   <- .virtualize_positions(j_seq, j_end,   vmap)
    keep     <- keep & !is.na(j_vstart) & !is.na(j_vend)
  } else {
    j_vstart <- NULL
    j_vend   <- NULL
  }

  gr2 <- gr[keep]
  if (length(gr2) == 0L) {
    return(GenomicRanges::GRanges(
      "__combined__",
      IRanges::IRanges(integer(0), integer(0))
    ))
  }

  out <- GenomicRanges::GRanges(
    rep("__combined__", sum(keep)),
    IRanges::IRanges(start = i_vstart[keep], end = i_vend[keep])
  )
  S4Vectors::mcols(out) <- S4Vectors::mcols(gr2)
  if (!is.null(j_vstart)) {
    S4Vectors::mcols(out)[[j_start_col]] <- j_vstart[keep]
    S4Vectors::mcols(out)[[j_end_col]]   <- j_vend[keep]
  }
  out
}

#' Genomic xscale for a single window
#'
#' @param window_gr A `GRanges` object containing one or more windows.
#' @param w Integer index of the window to use.
#' @return Numeric length-2 vector: `c(start, end)` of the selected window.
#' @keywords internal
.genomic_xscale <- function(window_gr, w) {
  c(BiocGenerics::start(window_gr)[w], BiocGenerics::end(window_gr)[w])
}

#' Infer a continuous x-range for a track from its elements' resolved data
#'
#' Used when the track has a non-genomic x mapping but no explicit
#' `scale_x`. Tries each element in turn: resolves against the track's
#' data + mapping and returns the range of the first numeric `x` vector.
#' Returns `NULL` when nothing usable is found.
#'
#' @param track A `SeqTrackR6` instance.
#' @return Numeric length-2 vector or `NULL`.
#' @keywords internal
.infer_x_range <- function(track) {
  for (elem in track$elements) {
    tryCatch(
      elem$resolve(track_data = track$data, track_mapping = track$mapping),
      error = function(e) NULL
    )
    x <- elem$resolved$x
    if (is.numeric(x) && length(x) > 0) {
      r <- range(x, na.rm = TRUE)
      if (all(is.finite(r))) return(as.numeric(r))
    }
  }
  NULL
}

#' Is a `map()` expression a genomic special (bare symbol)?
#'
#' Returns `TRUE` when `expr` is a bare symbol matching one of `start`,
#' `end`, `mid`, or `width` — the coordinate specials injected into the
#' mapping eval environment. Used to auto-detect genomic axes.
#'
#' @param expr An R language object from a `SeqMap`.
#' @return Logical scalar.
#' @keywords internal
.is_genomic_special <- function(expr) {
  if (is.null(expr)) return(FALSE)
  if (!is.symbol(expr)) return(FALSE)
  as.character(expr) %in% c("start", "end", "mid", "width")
}

#' Compute the xscale for a single track window
#'
#' Dispatches on `track$scale_x`: continuous uses `limits` (or a
#' `c(0, 1)` placeholder when `limits` is NULL), discrete uses
#' `c(0.5, n_levels + 0.5)`, and everything else (including an explicit
#' `seq_scale_genomic()`) falls back to the genomic window range.
#'
#' @param track A `SeqTrackR6` instance.
#' @param window_gr The track's `windows` `GRanges`.
#' @param w Integer index of the window.
#' @return Numeric length-2 vector: `c(x_min, x_max)`.
#' @keywords internal
.compute_track_xscale <- function(track, window_gr, w) {
  .scale_to_xrange(track$scale_x, window_gr, w)
}

#' Compute an x-range from a scale object (helper shared by x1 and x2)
#'
#' @param sx A `SeqPositionScale` or `NULL`.
#' @param window_gr A `GRanges` — used as the genomic fallback.
#' @param w Integer index of the window.
#' @return Numeric length-2 vector, or `NULL` when `sx` is NULL and no
#'   genomic fallback is appropriate (caller decides).
#' @keywords internal
.scale_to_xrange <- function(sx, window_gr, w) {
  if (!is.null(sx) && !is.null(sx$type)) {
    if (sx$type == "continuous")
      return(as.numeric(sx$limits %||% c(0, 1)))
    if (sx$type == "discrete" && !is.null(sx$levels))
      return(c(0.5, length(sx$levels) + 0.5))
    # genomic falls through to the window range
  }
  if (is.null(sx)) return(NULL)
  .genomic_xscale(window_gr, w)
}

#' Compute a y-range from a scale object (helper shared by y1 and y2)
#'
#' @param sy A `SeqPositionScale` or `NULL`.
#' @param y_windows Optional `GRanges` — used when the scale is genomic.
#' @return Numeric length-2 vector, or `NULL` when nothing can be derived.
#' @keywords internal
.scale_to_yrange <- function(sy, y_windows = NULL) {
  if (is.null(sy)) return(NULL)
  if (!is.null(sy$type)) {
    if (sy$type == "genomic") {
      win <- sy$windows %||% y_windows
      if (!is.null(win))
        return(c(min(BiocGenerics::start(win)), max(BiocGenerics::end(win))))
    }
    if (sy$type == "discrete" && !is.null(sy$levels))
      return(c(0.5, length(sy$levels) + 0.5))
    if (sy$type == "continuous" && !is.null(sy$limits))
      return(as.numeric(sy$limits))
  }
  NULL
}

#' Genomic yscale for a track
#'
#' Returns the yscale that should be used for a track's y-axis. When the track
#' uses a genomic y-axis (`uses_genomic_y == TRUE`), the range spans the full
#' extent of `track$y_windows`. Otherwise returns `c(0, 1)` as a placeholder
#' that callers may override using element data.
#'
#' @param track A `SeqTrackR6` instance.
#' @return Numeric length-2 vector: `c(y_min, y_max)`.
#' @keywords internal
.genomic_yscale <- function(track) {
  if (isTRUE(track$uses_genomic_y) && !is.null(track$y_windows)) {
    return(c(min(BiocGenerics::start(track$y_windows)),
             max(BiocGenerics::end(track$y_windows))))
  }
  c(0, 1)
}
