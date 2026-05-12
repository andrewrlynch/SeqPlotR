# ── seq_zoom — cross-track zoom/highlight polygon ────────────────────────────

#' SeqZoom R6 class
#'
#' Internal R6 generator backing [seq_zoom()]. Inherits from [`SeqLinkR6`].
#' Draws a filled quadrilateral projecting a genomic region in track `t0`
#' onto the same (or a different) region in track `t1`. Typically used to
#' connect an overview / ideogram track to a zoomed detail track.
#'
#' @keywords internal
SeqZoomR6 <- R6::R6Class("SeqZoom",
  inherit = SeqLinkR6,
  public = list(
    #' @description Construct a SeqZoomR6.
    #' @param data Optional `GRanges` or `data.frame` carrying the region(s).
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x0_end`. Optional:
    #'   `x1`, `x1_end` (default to `x0`, `x0_end` when absent), `chrom0`,
    #'   `chrom1`.
    #' @param t0,t1 Track identifiers. Must be specified explicitly for
    #'   plot-level placement.
    #' @param aesthetics Optional `SeqAes`. Recognised: `fill`, `color`,
    #'   `alpha`, `linewidth`, `stemOffset` (npc, default `0.01`).
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve the mapping and populate region fields. Unlike
    #'   other `SeqLink` subclasses, `seq_zoom` only requires `x0` / `x0_end`
    #'   — the `x1` / `x1_end` defaults mirror the `x0` region when absent.
    #' @param track_data Optional `GRanges`/`data.frame` from the parent track.
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

    #' @description Project the region(s) onto both tracks and build a
    #'   `data.frame` of 4-corner polygon coordinates in `self$coordCanvas`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Unused — `seq_zoom` requires explicit `t0`/`t1`.
    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      tid0 <- self$t0 %||% plot_track_index
      tid1 <- self$t1 %||% plot_track_index

      layout_t0  <- self$.resolve_track_ref(tid0, layout_all_tracks)
      layout_t1  <- self$.resolve_track_ref(tid1, layout_all_tracks)
      windows_t0 <- track_windows_list[[tid0]]
      windows_t1 <- track_windows_list[[tid1]]

      self$resolve(
        track_data    = layout_t0[[1]]$track_data,
        track_mapping = layout_t0[[1]]$track_mapping
      )

      self$coordCanvas <- NULL

      eff_data <- self$resolved$.data
      if (is.null(eff_data)) return(invisible())
      n <- if (is.data.frame(eff_data)) nrow(eff_data) else length(eff_data)
      if (n == 0L) return(invisible())

      x0_v  <- self$resolved$x0
      x0e_v <- self$resolved$x0_end
      if (is.null(x0_v) || is.null(x0e_v))
        stop("seq_zoom requires x0 and x0_end in map().", call. = FALSE)

      x1_v  <- self$resolved$x1     %||% x0_v
      x1e_v <- self$resolved$x1_end %||% x0e_v

      is_gr <- inherits(eff_data, "GRanges")
      chrom0_v <- self$resolved$chrom0 %||%
        (if (is_gr) as.character(GenomicRanges::seqnames(eff_data))
         else stop("seq_zoom with data.frame data requires chrom0 in map().",
                   call. = FALSE))
      chrom1_v <- self$resolved$chrom1 %||% chrom0_v

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      x0_v     <- .vec(x0_v,     n); x0e_v <- .vec(x0e_v, n)
      x1_v     <- .vec(x1_v,     n); x1e_v <- .vec(x1e_v, n)
      chrom0_v <- .vec(chrom0_v, n); chrom1_v <- .vec(chrom1_v, n)

      region0_gr <- GenomicRanges::GRanges(
        seqnames = chrom0_v,
        ranges   = IRanges::IRanges(start = pmin(x0_v, x0e_v),
                                    end   = pmax(x0_v, x0e_v))
      )
      region1_gr <- GenomicRanges::GRanges(
        seqnames = chrom1_v,
        ranges   = IRanges::IRanges(start = pmin(x1_v, x1e_v),
                                    end   = pmax(x1_v, x1e_v))
      )

      ov0 <- suppressWarnings(
        GenomicRanges::findOverlaps(region0_gr, windows_t0)
      )
      ov1 <- suppressWarnings(
        GenomicRanges::findOverlaps(region1_gr, windows_t1)
      )
      if (length(ov0) == 0L || length(ov1) == 0L) return(invisible())

      # Auto-detect attachment edges by comparing track y positions.
      p0_first <- layout_t0[[1]]
      p1_first <- layout_t1[[1]]
      t0_above <- p0_first$inner$y0 > p1_first$inner$y0
      stemOffset <- self$aesthetics$stemOffset %||% 0.01

      .project <- function(gr_row, x_left, x_right, pm, attach_edge,
                           offset_sign) {
        xs <- pm$xscale
        x_lc <- pmax(pmin(x_left,  xs[2]), xs[1])
        x_rc <- pmax(pmin(x_right, xs[2]), xs[1])
        u0 <- (x_lc - xs[1]) / diff(xs)
        u1 <- (x_rc - xs[1]) / diff(xs)
        box <- pm$inner
        xL <- box$x0 + u0 * (box$x1 - box$x0)
        xR <- box$x0 + u1 * (box$x1 - box$x0)
        yEdge <- if (attach_edge == "top") box$y1 else box$y0
        list(xL = xL, xR = xR, y = yEdge + offset_sign * stemOffset)
      }

      edge0 <- if (t0_above) "bottom" else "top"
      edge1 <- if (t0_above) "top"    else "bottom"
      sign0 <- if (t0_above) -1       else  1
      sign1 <- if (t0_above)  1       else -1

      qh0 <- S4Vectors::queryHits(ov0);   sh0 <- S4Vectors::subjectHits(ov0)
      qh1 <- S4Vectors::queryHits(ov1);   sh1 <- S4Vectors::subjectHits(ov1)

      rows_a <- lapply(seq_along(qh0), function(k) {
        i <- qh0[k]; w <- sh0[k]; pm <- layout_t0[[w]]
        proj <- .project(region0_gr[i], x0_v[i], x0e_v[i], pm, edge0, sign0)
        data.frame(idx = i, x00 = proj$xL, x10 = proj$xR, y0 = proj$y,
                   stringsAsFactors = FALSE)
      })
      rows_b <- lapply(seq_along(qh1), function(k) {
        i <- qh1[k]; w <- sh1[k]; pm <- layout_t1[[w]]
        proj <- .project(region1_gr[i], x1_v[i], x1e_v[i], pm, edge1, sign1)
        data.frame(idx = i, x01 = proj$xL, x11 = proj$xR, y1 = proj$y,
                   stringsAsFactors = FALSE)
      })

      a <- do.call(rbind, rows_a)
      b <- do.call(rbind, rows_b)
      merged <- merge(a, b, by = "idx")
      if (nrow(merged) == 0) return(invisible())

      self$coordCanvas <- merged[, c("x00", "x10", "x01", "x11", "y0", "y1")]
      invisible()
    },

    #' @description Draw each prepared polygon with [grid::grid.polygon()].
    draw = function() {
      df <- self$coordCanvas
      if (is.null(df) || nrow(df) == 0) return(invisible())
      fill  <- self$aesthetics$fill      %||% "grey50"
      alpha <- self$aesthetics$alpha     %||% 0.15
      col   <- self$aesthetics$color     %||% NA
      lwd   <- self$aesthetics$linewidth %||% 0.5
      fill2 <- grDevices::adjustcolor(fill, alpha.f = alpha)
      for (i in seq_len(nrow(df))) {
        grid::grid.polygon(
          x  = grid::unit(c(df$x00[i], df$x10[i], df$x11[i], df$x01[i]),
                          "npc"),
          y  = grid::unit(c(df$y0[i],  df$y0[i],  df$y1[i],  df$y1[i]),
                          "npc"),
          gp = grid::gpar(fill = fill2, col = col, lwd = lwd)
        )
      }
      invisible()
    }
  )
)

#' Zoom / highlight polygon between two tracks
#'
#' Projects a genomic region onto two tracks and connects the two projections
#' with a filled quadrilateral. Typically used to link an overview track to
#' a zoomed detail track (e.g. ideogram → locus). Always used as a plot-level
#' link — `t0` and `t1` must be specified explicitly.
#'
#' Map vocabulary:
#' \describe{
#'   \item{`x0`, `x0_end`}{required — genomic edges of the region in `t0`.}
#'   \item{`x1`, `x1_end`}{optional — genomic edges of the region in `t1`;
#'     default to `x0`, `x0_end` (same region projected onto both tracks).}
#'   \item{`chrom0`, `chrom1`}{optional for GRanges (default `seqnames(data)`);
#'     required for `data.frame`.}
#' }
#'
#' Attachment edges auto-select so that the upper track attaches from its
#' bottom and the lower track from its top.
#'
#' @param data Optional `GRanges` or `data.frame`.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers (required).
#' @param aesthetics Optional [aes()]: `fill` (default `"grey50"`),
#'   `alpha` (default `0.15`), `color` (border; default `NA`),
#'   `linewidth` (default `0.5`), `stemOffset` (default `0.01`).
#' @param ... Reserved.
#' @return A `SeqZoomR6` instance.
#' @export
seq_zoom <- function(data = NULL, mapping = NULL,
                     t0 = NULL, t1 = NULL,
                     aesthetics = aes(), ...) {
  SeqZoomR6$new(data = data, mapping = mapping,
                t0 = t0, t1 = t1,
                aesthetics = aesthetics, ...)
}
