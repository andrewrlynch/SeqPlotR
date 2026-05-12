# ── seq_string — Bezier string link (within- or cross-track) ─────────────────

#' SeqString R6 class
#'
#' Internal R6 generator backing [seq_string()]. Inherits from [`SeqLinkR6`].
#' Draws a smooth cubic-Bezier "string" between two anchors, typically used
#' to connect breakpoints between two tracks. C vs. S curve shape is inferred
#' from the resolved `strand0`/`strand1` fields when `aes(type = "auto")`
#' (the default); `aes(type = "c")` or `aes(type = "s")` forces the shape.
#'
#' @keywords internal
SeqStringR6 <- R6::R6Class("SeqString",
  inherit = SeqLinkR6,
  public = list(
    #' @description Construct a SeqStringR6.
    #' @param data Optional `GRanges` or `data.frame` carrying both anchors.
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x1`. Optional:
    #'   `chrom0`, `chrom1`, `strand0`, `strand1`, `y0`, `y1`.
    #' @param t0,t1 Track identifiers. Locked to the parent track when added
    #'   inside a [seq_track()] via `%+%`.
    #' @param aesthetics Optional `SeqAes`. Recognised: `color`, `linewidth`,
    #'   `alpha`, `type` (`"auto"`, `"c"`, `"s"`), `bulge`, `orientation`.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve both anchors, find overlaps with each referenced
    #'   track's windows, and store per-link canvas coordinates + resolved
    #'   curve type in `self$coordCanvas`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Fallback track reference when `t0`/`t1` are `NULL`.
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

      self$coordCanvas <- list()
      if (is.null(self$anchor0_gr) || length(self$anchor0_gr) == 0L)
        return(invisible())

      n    <- length(self$anchor0_gr)
      x0_g <- BiocGenerics::start(self$anchor0_gr)
      x1_g <- BiocGenerics::start(self$anchor1_gr)
      y0_v <- self$resolved$y0 %||% rep(0, n)
      y1_v <- self$resolved$y1 %||% y0_v
      s0_v <- as.character(BiocGenerics::strand(self$anchor0_gr))
      s1_v <- as.character(BiocGenerics::strand(self$anchor1_gr))

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      y0_v <- .vec(y0_v, n); y1_v <- .vec(y1_v, n)

      ov0_all <- rep(NA_integer_, n)
      ov1_all <- rep(NA_integer_, n)
      ov0 <- suppressWarnings(
        GenomicRanges::findOverlaps(self$anchor0_gr, windows_t0)
      )
      if (length(ov0) > 0L)
        ov0_all[S4Vectors::queryHits(ov0)] <- S4Vectors::subjectHits(ov0)
      ov1 <- suppressWarnings(
        GenomicRanges::findOverlaps(self$anchor1_gr, windows_t1)
      )
      if (length(ov1) > 0L)
        ov1_all[S4Vectors::queryHits(ov1)] <- S4Vectors::subjectHits(ov1)

      type_aes <- tolower(as.character(
        self$aesthetics$type %||% "auto"
      ))
      col_a   <- self$aesthetics$color     %||% "purple"
      lwd_a   <- self$aesthetics$linewidth %||% 1.5
      alpha_a <- self$aesthetics$alpha     %||% 0.8
      bulge_a <- self$aesthetics$bulge     %||% 0.04
      orient_a <- self$aesthetics$orientation %||% "*"

      valid <- !is.na(ov0_all) & !is.na(ov1_all)
      if (!any(valid)) return(invisible())

      for (i in which(valid)) {
        pm0 <- layout_t0[[ov0_all[i]]]
        pm1 <- layout_t1[[ov1_all[i]]]

        p0 <- .to_canvas(x0_g[i], y0_v[i], pm0)
        p1 <- .to_canvas(x1_g[i], y1_v[i], pm1)

        type_i <- if (type_aes %in% c("c", "s")) {
          type_aes
        } else {
          .string_type_from_strand(s0_v[i], s1_v[i], default = "c")
        }

        self$coordCanvas[[length(self$coordCanvas) + 1L]] <- list(
          x0          = p0$x,
          y0          = p0$y,
          x1          = p1$x,
          y1          = p1$y,
          strand0     = s0_v[i],
          strand1     = s1_v[i],
          type        = type_i,
          orientation = orient_a,
          col         = col_a,
          lwd         = lwd_a,
          alpha       = alpha_a,
          bulge       = bulge_a
        )
      }
      invisible()
    },

    #' @description Draw each prepared string with [drawSeqString()].
    draw = function() {
      if (is.null(self$coordCanvas) || length(self$coordCanvas) == 0L)
        return(invisible())
      for (r in self$coordCanvas) {
        drawSeqString(
          x0          = r$x0, y0 = r$y0,
          x1          = r$x1, y1 = r$y1,
          strand1     = r$strand0,
          strand2     = r$strand1,
          orientation = r$orientation,
          type        = r$type,
          bulge       = r$bulge,
          lwd         = r$lwd,
          col         = r$col,
          alpha       = r$alpha
        )
      }
      invisible()
    }
  )
)

#' Cubic-Bezier string link between two genomic loci
#'
#' Draws a smooth Bezier curve connecting two anchors. When `aes(type = "auto")`
#' (the default), the C vs. S shape is inferred from the resolved
#' `strand0`/`strand1` fields — opposing strands (`+/-`, `-/+`) produce an
#' S-curve, matching strands produce a C-curve.
#'
#' `data` may be a `GRanges` (anchor-1 columns in `mcols`) or a `data.frame`
#' (BEDPE-like). Required `map()` fields are `x0` and `x1`; `chrom0` and
#' `chrom1` are required for `data.frame` input and optional for `GRanges`
#' (default `seqnames(data)`). `y0` and `y1` default to `0`.
#'
#' @param data Optional `GRanges` or `data.frame`.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers. Locked to the parent track when added
#'   inside a [seq_track()].
#' @param aesthetics Optional [aes()]: `color` (default `"purple"`),
#'   `linewidth` (default `1.5`), `alpha` (default `0.8`), `type`
#'   (`"auto"`, `"c"`, or `"s"`; default `"auto"`), `bulge`, `orientation`.
#' @param ... Reserved.
#' @return A `SeqStringR6` instance.
#' @export
seq_string <- function(data = NULL, mapping = NULL,
                       t0 = NULL, t1 = NULL,
                       aesthetics = aes(), ...) {
  SeqStringR6$new(data = data, mapping = mapping,
                  t0 = t0, t1 = t1,
                  aesthetics = aesthetics, ...)
}
