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
    #' @field height Optional block height (fraction of the panel inner
    #'   height). When set, ribbons anchor at the facing block edges.
    height = NULL,
    #' @field smooth Logical. Draw smoothstep sides instead of straight ones.
    smooth = FALSE,

    #' @description Construct a SeqSyntenyR6.
    #' @param data Optional `GRanges` or `data.frame` carrying both region
    #'   pairs.
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x1` (and, for
    #'   `data.frame`, `chrom0`, `chrom1`). Optional: `x0_end`, `x1_end`,
    #'   `y0`, `y1`, `color`, `fill`, `twisted`.
    #' @param t0,t1 Track identifiers. Locked to the parent track when added
    #'   inside a [seq_track()] via `%+%`.
    #' @param aesthetics Optional `SeqAes`. Recognised: `fill`, `color`,
    #'   `alpha`, `linewidth`.
    #' @param height Optional block height (fraction of the panel inner
    #'   height) used to anchor at the facing block edges.
    #' @param smooth Logical. Smoothstep sides.
    #' @param layer `"above"` or `"below"`.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), height = NULL,
                          smooth = FALSE, layer = c("above", "below"), ...) {
      if (!is.null(height) &&
          (!is.numeric(height) || length(height) != 1L || !is.finite(height)))
        stop("seq_synteny height must be NULL or a single number.",
             call. = FALSE)
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics,
                       layer = match.arg(layer))
      self$height <- height
      self$smooth <- isTRUE(smooth)
    },

    #' @description Resolve both region endpoints, find which windows each
    #'   region's left anchor falls into, and store one polygon per region
    #'   in `self$coordCanvas`.
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

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      n <- length(self$anchor0_gr)
      # exact numeric edges where available, anchor starts otherwise
      x0_v <- if (is.numeric(self$resolved$x0)) .vec(self$resolved$x0, n)
              else BiocGenerics::start(self$anchor0_gr)
      x1_v <- if (is.numeric(self$resolved$x1)) .vec(self$resolved$x1, n)
              else BiocGenerics::start(self$anchor1_gr)
      x0e_v  <- self$resolved$x0_end %||% (x0_v + 1L)
      x1e_v  <- self$resolved$x1_end %||% (x1_v + 1L)
      y0_map <- self$resolved$y0
      y1_map <- self$resolved$y1
      fill_v <- self$resolved$fill %||%
                self$resolved$color %||%
                (self$aesthetics$fill  %||%
                 self$aesthetics$color %||% "grey50")
      tw_v <- as.logical(self$resolved$twisted %||% FALSE)

      x0e_v  <- .vec(x0e_v,  n)
      x1e_v  <- .vec(x1e_v,  n)
      fill_v <- .vec(fill_v, n)
      tw_v   <- .vec(tw_v,   n)
      tw_v[is.na(tw_v)] <- FALSE
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
      .yc <- function(pm) (pm$inner$y0 + pm$inner$y1) / 2
      .hh <- function(pm) self$height * (pm$inner$y1 - pm$inner$y0) / 2
      tt <- seq(0, 1, length.out = 16L)
      sm <- 3 * tt^2 - 2 * tt^3

      self$coordCanvas <- vector("list", n)
      for (i in seq_len(n)) {
        if (!valid[i]) next

        pm0 <- layout_t0[[ov0_all[i]]]
        pm1 <- layout_t1[[ov1_all[i]]]

        xBL <- .x_npc(x0_v[i],  pm0)
        xBR <- .x_npc(x0e_v[i], pm0)
        xTL <- .x_npc(x1_v[i],  pm1)
        xTR <- .x_npc(x1e_v[i], pm1)
        # twisted ribbons pair each t0 edge with the opposite t1 edge
        if (tw_v[i]) { tmp <- xTL; xTL <- xTR; xTR <- tmp }

        if (!is.null(self$height)) {
          up <- .yc(pm0) > .yc(pm1)
          yB_e <- if (up) .yc(pm0) - .hh(pm0) else .yc(pm0) + .hh(pm0)
          yT_e <- if (up) .yc(pm1) + .hh(pm1) else .yc(pm1) - .hh(pm1)
        } else {
          yB_e <- pm0$inner$y0
          yT_e <- pm1$inner$y1
        }
        yB <- if (!is.null(y0_map)) .y_npc(y0_map[i], pm0) else yB_e
        yT <- if (!is.null(y1_map)) .y_npc(y1_map[i], pm1) else yT_e

        if (self$smooth) {
          ys <- yB + (yT - yB) * tt
          x  <- c(xBL + (xTL - xBL) * sm, rev(xBR + (xTR - xBR) * sm))
          y  <- c(ys, rev(ys))
        } else {
          x <- c(xBL, xBR, xTR, xTL)
          y <- c(yB,  yB,  yT,  yT)
        }
        self$coordCanvas[[i]] <- list(
          x       = x,
          y       = y,
          fill    = fill_v[i],
          twisted = tw_v[i],
          x0_npc  = c(xBL, xBR),
          x1_npc  = c(xTL, xTR)
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
          ),
          name = "seqsynteny"
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
#' when absent. Numeric edges are used exactly (no rounding to whole bases).
#' `y0` and `y1` default to the bottom edge of t0's inner panel and the top
#' edge of t1's inner panel respectively — no scale conversion is applied in
#' that case. With `height` set, the defaults become the facing edges of
#' blocks of that height centred in each panel (the same geometry as
#' [seq_block()] and [seq_connector()]).
#'
#' The optional logical `map()` field `twisted` pairs `x0` with `x1_end` and
#' `x0_end` with `x1`, so the two sides cross (an inverted copy). With
#' `smooth = TRUE` each side follows a smoothstep curve (16 points per side),
#' which turns twisted ribbons into S shapes.
#'
#' @param data Optional `GRanges` or `data.frame`.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers. Locked to the parent track when added
#'   inside a [seq_track()].
#' @param aesthetics Optional [aes()]: `fill` / `color` (fill color),
#'   `alpha` (default `0.4`), `linewidth` (border width).
#' @param height Optional block height as a fraction of the panel inner
#'   height. `NULL` (default) anchors at the panel edges.
#' @param smooth Logical. Smoothstep sides instead of straight ones.
#' @param layer `"above"` (default) draws after track elements; `"below"`
#'   draws before them, so ribbons sit beneath blocks.
#' @param ... Reserved.
#' @return A `SeqSyntenyR6` instance.
#' @export
seq_synteny <- function(data = NULL, mapping = NULL,
                        t0 = NULL, t1 = NULL,
                        aesthetics = aes(), height = NULL,
                        smooth = FALSE, layer = c("above", "below"), ...) {
  SeqSyntenyR6$new(data = data, mapping = mapping,
                   t0 = t0, t1 = t1,
                   aesthetics = aesthetics, height = height,
                   smooth = smooth, layer = layer, ...)
}
