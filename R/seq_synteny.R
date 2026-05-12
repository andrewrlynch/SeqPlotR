# ── seq_synteny — filled trapezoid between two tracks ────────────────────────

#' SeqSynteny R6 class
#'
#' Internal R6 generator backing [seq_synteny()]. Inherits from [`SeqLinkR6`].
#' Connects a genomic region in track `t0` to a homologous region in track
#' `t1` with a filled quadrilateral. Each region is defined by a pair of
#' edges (`x0`/`x0_end` and `x1`/`x1_end`) plus an optional data-scale y for
#' each side.
#'
#' @keywords internal
SeqSyntenyR6 <- R6::R6Class("SeqSynteny",
  inherit = SeqLinkR6,
  public = list(
    #' @description Construct a SeqSyntenyR6.
    #' @param data Optional `GRanges` or `data.frame` carrying both region
    #'   pairs.
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x1` (and, for
    #'   `data.frame`, `chrom0`, `chrom1`). Optional: `x0_end`, `x1_end`,
    #'   `y0`, `y1`, `color`, `fill`.
    #' @param t0,t1 Track identifiers. Locked to the parent track when added
    #'   inside a [seq_track()] via `%+%`.
    #' @param aesthetics Optional `SeqAes`. Recognised: `fill`, `color`,
    #'   `alpha`, `linewidth`.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve both region endpoints, find which windows each
    #'   region's left anchor falls into, and store a 4-corner polygon per
    #'   region in `self$coordCanvas`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Fallback track reference.
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

      n      <- length(self$anchor0_gr)
      x0_v   <- BiocGenerics::start(self$anchor0_gr)
      x1_v   <- BiocGenerics::start(self$anchor1_gr)
      x0e_v  <- self$resolved$x0_end %||% (x0_v + 1L)
      x1e_v  <- self$resolved$x1_end %||% (x1_v + 1L)
      y0_map <- self$resolved$y0
      y1_map <- self$resolved$y1
      fill_v <- self$resolved$fill %||%
                self$resolved$color %||%
                (self$aesthetics$fill  %||%
                 self$aesthetics$color %||% "grey50")

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      x0e_v  <- .vec(x0e_v,  n)
      x1e_v  <- .vec(x1e_v,  n)
      fill_v <- .vec(fill_v, n)
      if (!is.null(y0_map)) y0_map <- .vec(y0_map, n)
      if (!is.null(y1_map)) y1_map <- .vec(y1_map, n)

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

      valid <- !is.na(ov0_all) & !is.na(ov1_all)

      .x_npc <- function(x_gen, pm) {
        xpr <- pm$xplot_range %||% pm$xscale
        u <- pmax(pmin((x_gen - xpr[1]) / diff(xpr), 1), 0)
        pm$inner$x0 + u * (pm$inner$x1 - pm$inner$x0)
      }
      .y_npc <- function(y_val, pm) {
        ypr <- pm$yplot_range %||% pm$yscale
        v <- pmax(pmin((y_val - ypr[1]) / diff(ypr), 1), 0)
        pm$inner$y0 + v * (pm$inner$y1 - pm$inner$y0)
      }

      self$coordCanvas <- vector("list", n)
      for (i in seq_len(n)) {
        if (!valid[i]) next

        pm0 <- layout_t0[[ov0_all[i]]]
        pm1 <- layout_t1[[ov1_all[i]]]

        xBL <- .x_npc(x0_v[i],  pm0)
        xBR <- .x_npc(x0e_v[i], pm0)
        xTL <- .x_npc(x1_v[i],  pm1)
        xTR <- .x_npc(x1e_v[i], pm1)

        yB <- if (!is.null(y0_map)) .y_npc(y0_map[i], pm0) else pm0$inner$y0
        yT <- if (!is.null(y1_map)) .y_npc(y1_map[i], pm1) else pm1$inner$y1

        self$coordCanvas[[i]] <- list(
          x    = c(xBL, xBR, xTR, xTL),
          y    = c(yB,  yB,  yT,  yT),
          fill = fill_v[i]
        )
      }
      invisible()
    },

    #' @description Draw each prepared polygon with [grid::grid.polygon()].
    draw = function() {
      if (is.null(self$coordCanvas) || length(self$coordCanvas) == 0L)
        return(invisible())
      alpha <- self$aesthetics$alpha     %||% 0.4
      col   <- self$aesthetics$color     %||% NA
      lwd   <- self$aesthetics$linewidth %||% 0.5
      for (poly in self$coordCanvas) {
        if (is.null(poly)) next
        grid::grid.polygon(
          x  = grid::unit(poly$x, "npc"),
          y  = grid::unit(poly$y, "npc"),
          gp = grid::gpar(
            fill = grDevices::adjustcolor(poly$fill, alpha.f = alpha),
            col  = col,
            lwd  = lwd
          )
        )
      }
      invisible()
    }
  )
)

#' Filled trapezoid connecting two tracks
#'
#' Draws a filled quadrilateral whose bottom edge spans a region in track
#' `t0` and whose top edge spans a homologous region in track `t1`. Use to
#' highlight syntenic blocks, orthologous segments, or otherwise paired
#' intervals across two stacked tracks.
#'
#' `data` may be a `GRanges` or a `data.frame`. Required `map()` fields are
#' `x0` and `x1`; `x0_end` and `x1_end` default to one base past `x0`/`x1`
#' when absent. `y0` and `y1` default to the bottom edge of t0's inner panel
#' and the top edge of t1's inner panel respectively — no scale conversion is
#' applied in that case.
#'
#' @param data Optional `GRanges` or `data.frame`.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers. Locked to the parent track when added
#'   inside a [seq_track()].
#' @param aesthetics Optional [aes()]: `fill` / `color` (fill color),
#'   `alpha` (default `0.4`), `linewidth` (border width).
#' @param ... Reserved.
#' @return A `SeqSyntenyR6` instance.
#' @export
seq_synteny <- function(data = NULL, mapping = NULL,
                        t0 = NULL, t1 = NULL,
                        aesthetics = aes(), ...) {
  SeqSyntenyR6$new(data = data, mapping = mapping,
                   t0 = t0, t1 = t1,
                   aesthetics = aesthetics, ...)
}
