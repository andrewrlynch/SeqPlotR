# ── seq_arc — within-track Bezier arch (no stems) ────────────────────────────

#' SeqArc R6 class
#'
#' Internal R6 generator backing [seq_arc()]. Inherits from [`SeqLinkR6`].
#' Draws a single Bezier arch between two genomic loci within one track.
#' Both anchors live on the same track (`t0 == t1`), and `%+%` locks them
#' to the parent track when added inside a [seq_track()].
#'
#' @keywords internal
SeqArcR6 <- R6::R6Class("SeqArc",
  inherit = SeqLinkR6,
  public = list(
    #' @description Construct a SeqArcR6.
    #' @param data Optional `GRanges` or `data.frame` carrying both anchors.
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x1`. Optional:
    #'   `chrom0`, `chrom1`, `y0`, `y1`, `height`.
    #' @param t0,t1 Track identifiers. Both default to the parent track when
    #'   the arc is added via `%+%` inside a [seq_track()].
    #' @param aesthetics Optional `SeqAes`. Recognised: `color`, `linewidth`,
    #'   `orientation`, `curve`, `height`.
    #' @param ... Unused — accepted for forward compatibility.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve the arc's mapping against the referenced track's
    #'   panels, find overlaps with the track's windows, and store per-link
    #'   canvas coordinates for `draw()`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Fallback track reference when `t0` is `NULL`.
    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      tid       <- self$t0 %||% plot_track_index
      layout_t  <- self$.resolve_track_ref(tid, layout_all_tracks)
      windows_t <- track_windows_list[[tid]]

      self$resolve(
        track_data    = layout_t[[1]]$track_data,
        track_mapping = layout_t[[1]]$track_mapping
      )

      self$coordCanvas <- list()
      if (is.null(self$anchor0_gr) || length(self$anchor0_gr) == 0L)
        return(invisible())

      n     <- length(self$anchor0_gr)
      x0_g  <- BiocGenerics::start(self$anchor0_gr)
      x1_g  <- BiocGenerics::start(self$anchor1_gr)
      y0_v  <- self$resolved$y0     %||% rep(0, n)
      y1_v  <- self$resolved$y1     %||% y0_v
      h_v   <- self$resolved$height %||% (self$aesthetics$height %||% 1)

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      y0_v <- .vec(y0_v, n); y1_v <- .vec(y1_v, n); h_v <- .vec(h_v, n)

      ov <- suppressWarnings(
        GenomicRanges::findOverlaps(self$anchor0_gr, windows_t)
      )
      if (length(ov) == 0L) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      orient_a <- self$aesthetics$orientation %||% "*"
      curve_a  <- self$aesthetics$curve       %||% "length"
      color_a  <- self$aesthetics$color       %||% "#1C1B1A"
      width_a  <- self$aesthetics$linewidth   %||% 1

      for (w in unique(sh)) {
        mask <- sh == w
        pm   <- layout_t[[w]]
        idx  <- qh[mask]

        p0   <- .to_canvas(x0_g[idx], y0_v[idx], pm)
        p1   <- .to_canvas(x1_g[idx], y1_v[idx], pm)
        ptop <- .to_canvas(x0_g[idx], h_v[idx],  pm)$y

        new_rows <- lapply(seq_along(idx), function(k) {
          list(x0 = p0$x[k], y0 = p0$y[k],
               x1 = p1$x[k], y1 = p1$y[k],
               top0 = ptop[k], top1 = ptop[k],
               orientation = orient_a,
               curve       = curve_a,
               arcColor    = color_a,
               arcWidth    = width_a)
        })
        self$coordCanvas <- append(self$coordCanvas, new_rows)
      }
      invisible()
    },

    #' @description Draw each prepared arch with `drawSeqArch()`. No stems
    #'   are rendered.
    draw = function() {
      if (is.null(self$coordCanvas) || length(self$coordCanvas) == 0L)
        return(invisible())
      for (r in self$coordCanvas) {
        drawSeqArch(
          x0 = r$x0,    y0 = r$y0,
          x1 = r$x1,    y1 = r$y1,
          top0 = r$top0, top1 = r$top1,
          orientation = r$orientation,
          curve       = r$curve,
          stemWidth = 0, stemColor = r$arcColor,
          arcWidth  = r$arcWidth, arcColor = r$arcColor
        )
      }
      invisible()
    }
  )
)

#' Bezier arch link between two genomic loci on the same track
#'
#' Single Bezier arch with no stems. Anchors live on the same track:
#' when added inside a [seq_track()] via `%+%`, both `t0` and `t1` are
#' locked to the parent track.
#'
#' `data` may be a `GRanges` (with anchor-1 columns in `mcols`) or a
#' `data.frame` (BEDPE-like). `map()` field names are explicit:
#' \describe{
#'   \item{`x0`, `x1`}{required — genomic position of each anchor.}
#'   \item{`chrom0`, `chrom1`}{optional for GRanges (default `seqnames(data)`);
#'     required for `data.frame`.}
#'   \item{`y0`, `y1`}{optional baselines (default `0`).}
#'   \item{`height`}{optional arch peak height in data-scale units (default `1`).}
#' }
#'
#' @param data Optional `GRanges` or `data.frame`. Falls back to the parent
#'   track's data.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers. Locked to the parent track when added
#'   inside a [seq_track()].
#' @param aesthetics Optional [aes()]: `color`, `linewidth`, `orientation`
#'   (`"+"` / `"-"` / `"*"`), `curve` (`"length"`, `"equal"`, or numeric).
#' @param ... Reserved.
#' @return A `SeqArcR6` instance.
#' @export
seq_arc <- function(data = NULL, mapping = NULL,
                    t0 = NULL, t1 = NULL,
                    aesthetics = aes(), ...) {
  SeqArcR6$new(data = data, mapping = mapping,
               t0 = t0, t1 = t1,
               aesthetics = aesthetics, ...)
}
