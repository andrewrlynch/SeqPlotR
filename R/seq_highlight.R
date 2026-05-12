# ── seq_highlight — multi-track region highlight band ────────────────────────

#' SeqHighlight R6 class
#'
#' Internal R6 generator backing [seq_highlight()]. Inherits from
#' [`SeqLinkR6`]. Draws a continuous filled highlight band that passes
#' through every track from `t0` down to `t1` (inclusive). Per-track
#' widths in NPC follow each track's own genomic scale, so the band
#' fans / compresses across scale changes. Adjacent tracks are bridged
#' by trapezoids whose top / bottom edges match the corresponding
#' track's projected `xL` / `xR`.
#'
#' @keywords internal
SeqHighlightR6 <- R6::R6Class("SeqHighlight",
  inherit = SeqLinkR6,
  public = list(
    #' @description Construct a SeqHighlightR6.
    #' @param data Optional `GRanges` or `data.frame` carrying highlight
    #'   region(s).
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x0_end`. Optional:
    #'   `chrom0` (auto-derived from GRanges seqnames). The same genomic
    #'   region is projected onto every track in `[t0..t1]`.
    #' @param t0,t1 Track identifiers (string or integer index) bracketing
    #'   the inclusive run of tracks the band passes through. When
    #'   `t1` is `NULL`, the highlight is restricted to `t0` only.
    #' @param aesthetics Optional `SeqAes`. Recognised: `fill`
    #'   (default `"grey50"`), `alpha` (default `0.25`), `color`
    #'   (border, default `NA`), `linewidth` (default `0.5`).
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve the mapping. Like [`SeqZoomR6`], this only
    #'   requires `x0` / `x0_end` — there is no `x1` / `x1_end`, since
    #'   the same genomic region is reused on every track in the band.
    #' @param track_data Optional `GRanges`/`data.frame` from the parent
    #'   track (used as a fallback when the link has no own `data`).
    #' @param track_mapping Optional `SeqMap` from the parent track.
    #' @return The link, invisibly.
    resolve = function(track_data = NULL, track_mapping = NULL) {
      eff_data <- self$data %||% track_data
      eff_mapping <- if (is.null(self$mapping) && is.null(track_mapping)) {
        NULL
      } else if (is.null(self$mapping)) {
        track_mapping
      } else if (is.null(track_mapping)) {
        self$mapping
      } else {
        merged <- as.list(track_mapping)
        for (nm in names(self$mapping)) merged[[nm]] <- self$mapping[[nm]]
        structure(merged, class = "SeqMap")
      }

      self$resolved          <- .resolve_mapping(eff_data, eff_mapping)
      self$resolved$.data    <- eff_data
      self$resolved$.mapping <- eff_mapping
      invisible(self)
    },

    #' @description Project each highlight region onto every track in
    #'   `[t0..t1]` and stash a per-rectangle `data.frame` in
    #'   `self$coordCanvas` with columns
    #'   `region_id`, `track_pos`, `window_idx`, `xL`, `xR`,
    #'   `y0_npc`, `y1_npc`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Unused — `seq_highlight` requires explicit
    #'   `t0` (and optionally `t1`).
    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      tid0 <- self$t0 %||% plot_track_index
      tid1 <- self$t1 %||% tid0
      if (is.null(tid0))
        stop("seq_highlight requires t0 to be set.", call. = FALSE)

      pos0 <- .highlight_track_pos(tid0, layout_all_tracks)
      pos1 <- .highlight_track_pos(tid1, layout_all_tracks)
      pos_lo <- min(pos0, pos1)
      pos_hi <- max(pos0, pos1)
      track_positions <- seq.int(pos_lo, pos_hi)

      # Resolve mapping using the first referenced track's data fallback.
      first_layout <- layout_all_tracks[[pos_lo]]
      self$resolve(
        track_data    = first_layout[[1]]$track_data,
        track_mapping = first_layout[[1]]$track_mapping
      )

      self$coordCanvas <- NULL

      eff_data <- self$resolved$.data
      if (is.null(eff_data)) return(invisible())
      n <- if (is.data.frame(eff_data)) nrow(eff_data) else length(eff_data)
      if (n == 0L) return(invisible())

      x0_v  <- self$resolved$x0
      x0e_v <- self$resolved$x0_end
      if (is.null(x0_v) || is.null(x0e_v))
        stop("seq_highlight requires x0 and x0_end in map().", call. = FALSE)

      is_gr <- inherits(eff_data, "GRanges")
      chrom0_v <- self$resolved$chrom0 %||%
        (if (is_gr) as.character(GenomicRanges::seqnames(eff_data))
         else stop("seq_highlight with data.frame data requires chrom0 in map().",
                   call. = FALSE))

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      x0_v     <- .vec(x0_v,     n); x0e_v <- .vec(x0e_v, n)
      chrom0_v <- .vec(chrom0_v, n)

      region_gr <- GenomicRanges::GRanges(
        seqnames = chrom0_v,
        ranges   = IRanges::IRanges(start = pmin(x0_v, x0e_v),
                                    end   = pmax(x0_v, x0e_v))
      )

      rows <- list()
      for (pos in track_positions) {
        layout_t <- layout_all_tracks[[pos]]
        windows_t <- track_windows_list[[pos]]
        if (is.null(layout_t) || is.null(windows_t)) next

        ov <- suppressWarnings(
          GenomicRanges::findOverlaps(region_gr, windows_t)
        )
        if (length(ov) == 0L) next
        qh <- S4Vectors::queryHits(ov)
        sh <- S4Vectors::subjectHits(ov)
        for (k in seq_along(qh)) {
          i <- qh[k]; w <- sh[k]
          pm <- layout_t[[w]]
          if (is.null(pm)) next
          proj <- .highlight_project(x0_v[i], x0e_v[i], pm)
          rows[[length(rows) + 1L]] <- data.frame(
            region_id  = i,
            track_pos  = pos,
            window_idx = w,
            xL         = proj$xL,
            xR         = proj$xR,
            y0_npc     = proj$y0,
            y1_npc     = proj$y1,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(rows) == 0L) return(invisible())
      self$coordCanvas <- do.call(rbind, rows)
      invisible()
    },

    #' @description Render the highlight band: one filled rectangle per
    #'   per-track / per-window slice, and a connecting trapezoid in the
    #'   gap between adjacent tracks (matched by `window_idx`).
    draw = function() {
      df <- self$coordCanvas
      if (is.null(df) || nrow(df) == 0L) return(invisible())
      fill  <- self$aesthetics$fill      %||% "grey50"
      alpha <- self$aesthetics$alpha     %||% 0.25
      col   <- self$aesthetics$color     %||% NA
      lwd   <- self$aesthetics$linewidth %||% 0.5
      fill2 <- grDevices::adjustcolor(fill, alpha.f = alpha)
      gp <- grid::gpar(fill = fill2, col = col, lwd = lwd)

      regions <- split(df, df$region_id)
      for (rdf in regions) {
        # Per-track rectangles inside each panel
        for (i in seq_len(nrow(rdf))) {
          grid::grid.polygon(
            x  = grid::unit(c(rdf$xL[i], rdf$xR[i],
                              rdf$xR[i], rdf$xL[i]), "npc"),
            y  = grid::unit(c(rdf$y0_npc[i], rdf$y0_npc[i],
                              rdf$y1_npc[i], rdf$y1_npc[i]), "npc"),
            gp = gp
          )
        }

        # Bridging trapezoids between adjacent tracks (matched by window_idx).
        # Order tracks top-to-bottom by y1_npc descending so the upper track
        # comes first in each pair.
        ord <- order(rdf$y1_npc, decreasing = TRUE)
        sorted <- rdf[ord, , drop = FALSE]
        positions <- unique(sorted$track_pos)
        if (length(positions) < 2L) next
        for (p_i in seq_len(length(positions) - 1L)) {
          upper <- sorted[sorted$track_pos == positions[p_i],   , drop = FALSE]
          lower <- sorted[sorted$track_pos == positions[p_i + 1L], , drop = FALSE]
          if (nrow(upper) == 0L || nrow(lower) == 0L) next
          # Match by window_idx; if exactly one in each, just connect.
          if (nrow(upper) == 1L && nrow(lower) == 1L) {
            pairs <- list(list(u = upper[1, ], l = lower[1, ]))
          } else {
            common <- intersect(upper$window_idx, lower$window_idx)
            pairs <- lapply(common, function(w) list(
              u = upper[upper$window_idx == w, ][1, ],
              l = lower[lower$window_idx == w, ][1, ]
            ))
          }
          for (pr in pairs) {
            u <- pr$u; l <- pr$l
            # Upper rectangle bottom edge -> lower rectangle top edge
            y_top <- min(u$y0_npc, u$y1_npc)
            y_bot <- max(l$y0_npc, l$y1_npc)
            grid::grid.polygon(
              x  = grid::unit(c(u$xL, u$xR, l$xR, l$xL), "npc"),
              y  = grid::unit(c(y_top, y_top, y_bot, y_bot), "npc"),
              gp = gp
            )
          }
        }
      }
      invisible()
    }
  )
)

# ── helpers ──────────────────────────────────────────────────────────────────

# Convert a track reference (string or integer) to its integer position in
# `layout_all_tracks`. Mirrors `.resolve_track_ref()` but returns the position
# rather than the panel-bounds list.
.highlight_track_pos <- function(ref, layout_all_tracks) {
  if (is.null(ref))
    stop("seq_highlight reference is NULL - t0/t1 must be set.", call. = FALSE)
  if (is.character(ref)) {
    pos <- match(ref, names(layout_all_tracks))
    if (is.na(pos))
      stop("seq_highlight references unknown track_id '", ref, "'.",
           call. = FALSE)
    return(pos)
  }
  if (is.numeric(ref)) {
    idx <- as.integer(ref)
    if (idx < 1L || idx > length(layout_all_tracks))
      stop("seq_highlight references track index ", idx,
           " which is out of bounds (1..", length(layout_all_tracks), ").",
           call. = FALSE)
    return(idx)
  }
  stop("seq_highlight reference must be a track_id string or integer index.",
       call. = FALSE)
}

# Project a genomic [x_left, x_right] range onto a panel's npc inner box.
# Returns the four NPC coordinates (xL, xR, y0, y1) describing the rect.
.highlight_project <- function(x_left, x_right, pm) {
  xs <- pm$xscale
  x_lc <- pmax(pmin(x_left,  xs[2]), xs[1])
  x_rc <- pmax(pmin(x_right, xs[2]), xs[1])
  u0 <- (x_lc - xs[1]) / diff(xs)
  u1 <- (x_rc - xs[1]) / diff(xs)
  box <- pm$inner
  list(xL = box$x0 + u0 * (box$x1 - box$x0),
       xR = box$x0 + u1 * (box$x1 - box$x0),
       y0 = box$y0, y1 = box$y1)
}

#' Cross-track region highlight band
#'
#' Draws a continuous filled highlight band that passes through every
#' track from `t0` down to `t1` (inclusive), with per-track widths
#' determined by each track's own genomic scale. Adjacent tracks are
#' bridged by trapezoids in the inter-track gap, so the band fans /
#' compresses smoothly across tracks with different windows. Useful for
#' ChIP-seq style stacks where a single locus or region of interest
#' should be highlighted across many panels at once, or for an overview
#' track stacked above a zoomed detail.
#'
#' Always used as a plot-level link — `t0` must be specified explicitly.
#' Set `t1 = NULL` (or omit it) to highlight a single track only.
#'
#' Map vocabulary:
#' \describe{
#'   \item{`x0`, `x0_end`}{required — genomic edges of the highlight
#'     region. The same genomic region is projected onto every track in
#'     `[t0..t1]` and naturally compresses or expands per-track based on
#'     each track's window range.}
#'   \item{`chrom0`}{optional for GRanges (default `seqnames(data)`);
#'     required for `data.frame`.}
#' }
#'
#' @param data Optional `GRanges` or `data.frame` of one or more
#'   highlight regions.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers (string `track_id` or integer index)
#'   bracketing the inclusive run of tracks. `t0` is required;
#'   `t1` defaults to `t0` (single-track highlight).
#' @param aesthetics Optional [aes()]: `fill` (default `"grey50"`),
#'   `alpha` (default `0.25`), `color` (border; default `NA`),
#'   `linewidth` (default `0.5`).
#' @param ... Reserved.
#' @return A `SeqHighlightR6` instance.
#' @export
seq_highlight <- function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
  SeqHighlightR6$new(data = data, mapping = mapping,
                     t0 = t0, t1 = t1,
                     aesthetics = aesthetics, ...)
}
