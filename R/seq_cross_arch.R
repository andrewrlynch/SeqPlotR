# CROSS-TRACK ARCH

# R6 generator backing seq_cross_arch(); stems on the source track, tops on the target track.
SeqCrossArchR6 <- R6::R6Class("SeqCrossArch",
  inherit = SeqArchR6,
  public = list(
    target          = NULL,
    legend_arg      = NULL,
    legend_extra    = NULL,
    legend_position = "canvas_margin",
    legend_side     = "right",
    y_fields        = c("y0", "marker"),

    initialize = function(data = NULL, mapping = NULL, target = NULL,
                          t0 = NULL, t1 = NULL, aesthetics = aes(),
                          legend = NULL, legend_extra = NULL,
                          show_legend = TRUE,
                          legend_position = "canvas_margin",
                          legend_side = "right", ...) {
      if (inherits(data, "SeqMap") && is.null(mapping)) {
        mapping <- data
        data    <- NULL
      }
      if (is.null(target) || length(target) != 1L || is.na(target))
        stop("seq_cross_arch needs a single `target` track_id.", call. = FALSE)
      # x0 defaults to start; an unmapped x1 means single stems; chrom1 defaults to chrom0
      m <- unclass(mapping %||% list())
      if (is.null(m$x0)) m$x0 <- quote(start)
      if (is.null(m$x1)) m$x1 <- NA
      if (is.null(m$chrom1) && !is.null(m$chrom0)) m$chrom1 <- m$chrom0
      super$initialize(data = data, mapping = structure(m, class = "SeqMap"),
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
      self$target          <- target
      self$show_legend     <- show_legend
      self$legend_arg      <- legend
      self$legend_extra    <- legend_extra
      self$legend_position <- legend_position
      self$legend_side     <- legend_side
    },

    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      src <- self$t0 %||% plot_track_index
      tgt <- self$target
      src_panels <- self$.resolve_track_ref(src, layout_all_tracks)
      tgt_panels <- self$.resolve_track_ref(tgt, layout_all_tracks)
      .id <- function(ref) if (is.numeric(ref))
        names(layout_all_tracks)[as.integer(ref)] else ref
      if (identical(.id(src), .id(tgt)))
        stop("seq_cross_arch: target '", tgt, "' is the source track; ",
             "use seq_arch() for arches within one track.", call. = FALSE)
      win_s <- track_windows_list[[src]]
      win_t <- track_windows_list[[tgt]]
      same <- length(win_s) == length(win_t) &&
        all(as.character(GenomicRanges::seqnames(win_s)) ==
              as.character(GenomicRanges::seqnames(win_t))) &&
        all(BiocGenerics::start(win_s) == BiocGenerics::start(win_t)) &&
        all(BiocGenerics::end(win_s) == BiocGenerics::end(win_t))
      if (!same)
        stop("seq_cross_arch: target track '", tgt, "' must use the same ",
             "windows as the source track '", .id(src), "'.", call. = FALSE)

      self$resolve(track_data    = src_panels[[1]]$track_data,
                   track_mapping = src_panels[[1]]$track_mapping)
      self$coordCanvas <- list()
      self$stubs       <- list()
      self$legend <- .glyph_element_legend(
        NULL, self$legend_arg, legend_extra = self$legend_extra,
        position = self$legend_position, side = self$legend_side)
      if (is.null(self$anchor0_gr) || length(self$anchor0_gr) == 0L)
        return(invisible())

      n <- length(self$anchor0_gr)
      r <- self$resolved
      a <- self$aesthetics
      .col <- function(v, default) if (is.null(v)) rep(default, n)
                                   else rep_len(v, n)
      single <- self$anchor1_missing %||% rep(FALSE, n)
      ids    <- if (is.null(r$id)) NULL else as.character(.col(r$id, NA))

      ov0 <- suppressWarnings(GenomicRanges::findOverlaps(
        self$anchor0_gr, win_s, select = "first", ignore.strand = TRUE))
      ov1 <- suppressWarnings(GenomicRanges::findOverlaps(
        self$anchor1_gr, win_s, select = "first", ignore.strand = TRUE))
      bad <- which(is.na(ov0) | is.na(ov1))
      if (length(bad)) {
        i  <- bad[1]
        gr <- if (is.na(ov0[i])) self$anchor0_gr[i] else self$anchor1_gr[i]
        stop(sprintf("seq_cross_arch: row %d%s anchor %s:%s is outside every window",
                     i, if (is.null(ids)) "" else sprintf(" (%s)", ids[i]),
                     as.character(GenomicRanges::seqnames(gr)),
                     format(BiocGenerics::start(gr), scientific = FALSE)),
             call. = FALSE)
      }

      x0g    <- BiocGenerics::start(self$anchor0_gr)
      x1g    <- BiocGenerics::start(self$anchor1_gr)
      y0     <- as.numeric(.col(r$y0, 0))
      top    <- as.numeric(.col(r$top, NA))
      lo     <- as.numeric(.col(r$bar_lo, NA))
      hi     <- as.numeric(.col(r$bar_hi, NA))
      marker <- as.numeric(.col(r$marker, NA))
      glyph  <- as.character(.col(r$glyph, NA))
      caut   <- as.logical(.col(r$caution, FALSE))
      caut[is.na(caut)] <- FALSE
      col    <- .col(r$color, a$color %||% "#1C1B1A")
      if (!.looks_like_color(col)) col <- .auto_scale_colors(col)$colors
      col    <- .col_alpha(col, .col(r$alpha, a$alpha %||% 1))
      lwd    <- a$linewidth %||% 1.2

      # tops and bars map through the target panel and its transform; stem values are already transformed
      for (i in seq_len(n)) {
        ps <- src_panels[[ov0[i]]]
        pt <- tgt_panels[[ov0[i]]]
        b0 <- seq_to_canvas(x0g[i], y0[i], ps, apply_transform = FALSE)
        t0 <- if (is.na(top[i]))
          seq_to_canvas(x0g[i], pt$yscale[2], pt, apply_transform = FALSE)$y
        else seq_to_canvas(x0g[i], top[i], pt)$y
        m0 <- seq_to_canvas(x0g[i], marker[i], ps, apply_transform = FALSE)$y
        bl0 <- seq_to_canvas(x0g[i], c(lo[i], hi[i]), pt)$y
        rec <- list(row = i, x0 = b0$x, y0 = b0$y, x1 = NA_real_, y1 = NA_real_,
                    top0 = t0, top1 = NA_real_, single = single[i],
                    color = col[i], lwd = lwd, glyph = glyph[i],
                    lo0 = bl0[1], hi0 = bl0[2], lo1 = NA_real_, hi1 = NA_real_,
                    ytop = t0, marker0 = m0, marker1 = NA_real_,
                    caution = caut[i])
        if (!single[i]) {
          ps1 <- src_panels[[ov1[i]]]
          pt1 <- tgt_panels[[ov1[i]]]
          b1  <- seq_to_canvas(x1g[i], y0[i], ps1, apply_transform = FALSE)
          t1  <- if (is.na(top[i]))
            seq_to_canvas(x1g[i], pt1$yscale[2], pt1, apply_transform = FALSE)$y
          else seq_to_canvas(x1g[i], top[i], pt1)$y
          bl1 <- seq_to_canvas(x1g[i], c(lo[i], hi[i]), pt1)$y
          rec$x1 <- b1$x; rec$y1 <- b1$y; rec$top1 <- t1
          rec$lo1 <- bl1[1]; rec$hi1 <- bl1[2]; rec$ytop <- max(t0, t1)
          rec$marker1 <- seq_to_canvas(x1g[i], marker[i], ps1,
                                       apply_transform = FALSE)$y
        }
        self$coordCanvas[[i]] <- rec
      }
      invisible()
    },

    draw = function() {
      rows <- self$coordCanvas
      if (!length(rows)) return(invisible())
      a <- self$aesthetics
      curve <- a$curve %||% 0.04
      smm <- a$single_mm %||% 1.6
      for (r in rows) {
        if (!r$single) {
          drawSeqArch(x0 = r$x0, y0 = r$y0, x1 = r$x1, y1 = r$y1,
                      top0 = r$top0, top1 = r$top1, orientation = "*",
                      curve = curve, stemWidth = r$lwd, arcWidth = r$lwd,
                      arcColor = r$color, stemColor = r$color)
        } else {
          # single anchor: one stem with a pennant at its top
          grid::grid.segments(r$x0, r$y0, r$x0, r$top0,
                              name = paste0("seqplotr.cross_arch.stem.", r$row),
                              gp = grid::gpar(col = r$color, lwd = r$lwd))
          grid::grid.polygon(
            grid::unit(r$x0, "npc") + grid::unit(c(0, 0, smm), "mm"),
            grid::unit(r$top0, "npc") + grid::unit(c(0, -0.8, -0.4) * smm, "mm"),
            name = paste0("seqplotr.cross_arch.pennant.", r$row),
            gp = grid::gpar(fill = r$color, col = NA))
        }
      }

      .v <- function(f) vapply(rows, function(r) as.numeric(r[[f]]), numeric(1))
      dbl <- !vapply(rows, `[[`, logical(1), "single")
      cols <- vapply(rows, `[[`, character(1), "color")
      glyph <- vapply(rows, function(r) as.character(r$glyph), character(1))
      gx <- c(.v("x0"), .v("x1")[dbl])
      gy <- c(.v("top0"), .v("top1")[dbl])
      .draw_interval_glyphs(gx, gy, c(.v("lo0"), .v("lo1")[dbl]),
                            c(.v("hi0"), .v("hi1")[dbl]),
                            c(glyph, glyph[dbl]), c(cols, cols[dbl]),
                            alpha = 1, size = a$size %||% 2,
                            bar_lwd = a$bar_linewidth %||% 4,
                            bar_alpha = a$bar_alpha %||% 1,
                            cap_mm = a$cap_mm %||% 2.5,
                            points = a$glyph_points %||% TRUE,
                            name = "seqplotr.cross_arch.glyph")

      my <- c(.v("marker0"), .v("marker1")[dbl])
      s  <- which(is.finite(my))
      if (length(s))
        grid::grid.points(grid::unit(gx[s], "npc"), grid::unit(my[s], "npc"),
                          pch = a$marker_shape %||% 16,
                          name = "seqplotr.cross_arch.marker",
                          gp = grid::gpar(col = c(cols, cols[dbl])[s],
                                          cex = a$marker_size %||% 0.8))

      # caution marks above the arch peak, or beside the pennant of a single stem
      cmm <- a$caution_mm %||% 1.6
      peak <- if (is.numeric(curve)) 0.75 * curve else 0
      for (r in rows) {
        if (!isTRUE(r$caution)) next
        if (r$single)
          .draw_caution_marks(r$x0, r$top0, mm = cmm, dx_mm = smm + 0.5,
                              dy_mm = -0.8 * smm,
                              name = paste0("seqplotr.cross_arch.caution.", r$row))
        else
          .draw_caution_marks((r$x0 + r$x1) / 2, r$ytop + peak, mm = cmm,
                              dx_mm = -cmm / 2, dy_mm = 0.5,
                              name = paste0("seqplotr.cross_arch.caution.", r$row))
      }
      invisible()
    }
  )
)

#' Arches whose stems rise from one track to tops on another
#'
#' Draws, per row, vertical stems at one or two genomic anchors starting at
#' `y0` on the track the element is added to, joined by a Bezier arch whose
#' tops sit at `top` on the `target` track. Both tracks must use the same
#' windows: x is shared and only y maps across tracks, through the target
#' track's panels (and its y transform, e.g. [seq_scale_count()]). Rows
#' without a second anchor draw a single stem with a pennant. Optional
#' interval glyphs (see [seq_interval_glyph()]) are drawn on both stems at
#' the top, and an optional marker point on each stem at `marker` in source
#' track units. An anchor outside every window is an error naming the row.
#'
#' `coordCanvas` holds one list per row with canvas npc values `row, x0, y0,
#' x1, y1, top0, top1, single, color, lwd, glyph, lo0, hi0, lo1, hi1, ytop,
#' marker0, marker1, caution`.
#'
#' @param data `GRanges` or `data.frame` with one row per arch.
#' @param mapping `SeqMap` from [map()]. Fields: `x0` (default `start`),
#'   `chrom0` (default `seqnames` for GRanges), `x1` and `chrom1` (second
#'   anchor; `NA` or unmapped gives a single stem; `chrom1` defaults to
#'   `chrom0`), `top` (arch top in target track units; `NA` peaks at the
#'   target axis top), `y0` (stem base, default 0), `color`, `alpha`,
#'   `glyph`, `bar_lo`, `bar_hi` (glyph type and bar in target track units),
#'   `marker` (point on each stem), `caution` (logical; draws a caution mark)
#'   and `id` (row label used in error messages).
#' @param target `track_id` of the track holding the arch tops.
#' @param aesthetics `SeqAes` from [aes()]: `linewidth` (1.2), `curve`
#'   (0.04 npc), `color`, `alpha`, `bar_linewidth` (4), `bar_alpha` (1),
#'   `glyph_points` (`TRUE`), `size` (glyph point size in mm, 2), `cap_mm`
#'   (2.5), `marker_size` (0.8), `marker_shape` (16), `single_mm` (pennant
#'   width, 1.6), `caution_mm` (1.6).
#' @param legend `NULL` (no keys of its own), `FALSE`, a `SeqLegendSpec` or a
#'   list of `LegendKey` objects used as given.
#' @param legend_extra Optional `LegendKey` or list of them appended.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the built legend.
#' @param ... Passed to the generator (e.g. `t0` for a plot-level link).
#' @return A `SeqCrossArch` link element; add it to the source track.
#' @examples
#' library(GenomicRanges)
#' win <- GRanges("chr1", IRanges(1, 1e6))
#' sv <- data.frame(chr = "chr1", pos1 = 2e5, pos2 = 6e5, t = 0.4)
#' p <- seq_plot() %+%
#'   seq_track(track_id = "top", windows = win,
#'             scale_y = seq_scale_continuous(limits = c(0, 1))) %+%
#'   seq_track(track_id = "base", windows = win,
#'             scale_y = seq_scale_continuous(limits = c(0, 4))) %+%
#'   seq_cross_arch(sv, map(chrom0 = chr, x0 = pos1, x1 = pos2, top = t),
#'                  target = "top")
#' @export
seq_cross_arch <- function(data = NULL, mapping = NULL, target,
                           aesthetics = aes(), legend = NULL,
                           legend_extra = NULL, show_legend = TRUE,
                           legend_position = "canvas_margin",
                           legend_side = "right", ...) {
  if (missing(target))
    stop("seq_cross_arch needs a `target` track_id.", call. = FALSE)
  SeqCrossArchR6$new(data = data, mapping = mapping, target = target,
                     aesthetics = aesthetics, legend = legend,
                     legend_extra = legend_extra, show_legend = show_legend,
                     legend_position = legend_position,
                     legend_side = legend_side, ...)
}
