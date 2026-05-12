# ── seq_arch — Bezier arch with vertical stems and partial-window stubs ──────

#' SeqArch R6 class
#'
#' Internal R6 generator backing [seq_arch()]. Inherits from [`SeqLinkR6`].
#' Same anchor pattern as [`SeqArcR6`] but draws full vertical stems from the
#' baseline to the arch endpoints and renders partial-window stubs when only
#' one anchor is visible.
#'
#' @keywords internal
SeqArchR6 <- R6::R6Class("SeqArch",
  inherit = SeqLinkR6,
  public = list(
    #' @field stubs List of partial-window stub descriptions populated by
    #'   `prep()`. Each entry: `list(x, y0, y1, dir, partner, color, width)`.
    stubs = list(),

    #' @description Construct a SeqArchR6.
    #' @param data Optional `GRanges` or `data.frame` with both anchors.
    #' @param mapping Optional `SeqMap`. Required: `x0`, `x1`. Optional:
    #'   `chrom0`, `chrom1`, `strand0`, `strand1`, `y0`, `y1`, `height`.
    #' @param t0,t1 Track identifiers. Locked to the parent track when added
    #'   inside a [seq_track()] via `%+%`.
    #' @param aesthetics Optional `SeqAes`. Recognised: `color`, `linewidth`,
    #'   `arcColor`, `stemColor`, `arcWidth`, `stemWidth`, `orientation`,
    #'   `curve`, `height`, `plotStubs`, `stubAngle`, `stubLength`.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(), ...) {
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
      self$stubs <- list()
    },

    #' @description Resolve both anchors against their referenced tracks,
    #'   assign each anchor to a window via `findOverlaps()`, and build
    #'   `coordCanvas` for fully-visible links plus `stubs` for half-visible
    #'   links.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Fallback identifier for `t0` / `t1` when both
    #'   are `NULL` (within-track placement).
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
      self$stubs       <- list()
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

      ov0_all <- rep(NA_integer_, n)
      ov1_all <- rep(NA_integer_, n)
      # findOverlaps warns when seqlevels do not intersect (cross-chromosome
      # links to a windowed track) — that's the very case `seq_recon`
      # classifies as a translocation, so the warning is expected noise.
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

      orient_a   <- .vec(self$aesthetics$orientation %||% "*",      n)
      curve_a    <- .vec(self$aesthetics$curve       %||% "length", n)
      arc_col_a  <- .vec(self$aesthetics$arcColor    %||%
                          (self$aesthetics$color %||% "#1C1B1A"),   n)
      stem_col_a <- .vec(self$aesthetics$stemColor   %||%
                          (self$aesthetics$color %||% "#1C1B1A"),   n)
      arc_w_a    <- .vec(self$aesthetics$arcWidth    %||%
                          (self$aesthetics$linewidth %||% 1),       n)
      stem_w_a   <- .vec(self$aesthetics$stemWidth   %||%
                          (self$aesthetics$linewidth %||% 1),       n)

      # Stubs ---------------------------------------------------------------
      left_only  <- which(!is.na(ov0_all) &  is.na(ov1_all))
      right_only <- which( is.na(ov0_all) & !is.na(ov1_all))

      if (isTRUE(self$aesthetics$plotStubs %||% TRUE)) {
        for (i in c(left_only, right_only)) {
          use1      <- i %in% left_only
          tid_local <- if (use1) tid0 else tid1
          win_local <- if (use1) ov0_all[i] else ov1_all[i]
          pm_local  <- layout_all_tracks[[tid_local]][[win_local]]

          x_local <- if (use1) x0_g[i] else x1_g[i]
          y_local <- if (use1) y0_v[i] else y1_v[i]

          p_base <- .to_canvas(x_local, y_local, pm_local)
          p_top  <- .to_canvas(x_local, h_v[i], pm_local)

          partner_chr <- if (use1)
            as.character(GenomicRanges::seqnames(self$anchor1_gr[i]))
          else
            as.character(GenomicRanges::seqnames(self$anchor0_gr[i]))

          self$stubs[[length(self$stubs) + 1L]] <- list(
            x       = p_base$x,
            y0      = p_base$y,
            y1      = p_top$y,
            dir     = if (use1) 1L else -1L,
            partner = partner_chr,
            color   = stem_col_a[i],
            width   = stem_w_a[i]
          )
        }
      }

      # Full arches ---------------------------------------------------------
      valid <- !is.na(ov0_all) & !is.na(ov1_all)
      if (!any(valid)) return(invisible())

      for (i in which(valid)) {
        pm0  <- layout_t0[[ov0_all[i]]]
        pm1  <- layout_t1[[ov1_all[i]]]

        p0   <- .to_canvas(x0_g[i], y0_v[i], pm0)
        p1   <- .to_canvas(x1_g[i], y1_v[i], pm1)
        top0 <- .to_canvas(x0_g[i], h_v[i],  pm0)$y
        top1 <- .to_canvas(x1_g[i], h_v[i],  pm1)$y

        self$coordCanvas[[length(self$coordCanvas) + 1L]] <- list(
          x0 = p0$x,    y0 = p0$y,
          x1 = p1$x,    y1 = p1$y,
          top0 = top0,  top1 = top1,
          orientation = orient_a[i],
          curve       = curve_a[i],
          arcColor    = arc_col_a[i],
          stemColor   = stem_col_a[i],
          arcWidth    = arc_w_a[i],
          stemWidth   = stem_w_a[i]
        )
      }
      invisible()
    },

    #' @description Draw arches with stems, then partial-window stubs (when
    #'   `aes(plotStubs = TRUE)`).
    draw = function() {
      if (!is.null(self$coordCanvas) && length(self$coordCanvas) > 0L) {
        for (r in self$coordCanvas) {
          drawSeqArch(
            x0 = r$x0, y0 = r$y0, x1 = r$x1, y1 = r$y1,
            top0 = r$top0, top1 = r$top1,
            orientation = r$orientation, curve = r$curve,
            stemWidth = r$stemWidth, stemColor = r$stemColor,
            arcWidth  = r$arcWidth,  arcColor  = r$arcColor
          )
        }
      }

      if (isTRUE(self$aesthetics$plotStubs %||% TRUE) &&
          length(self$stubs) > 0L) {
        angle_deg <- self$aesthetics$stubAngle  %||% 45
        L_npc     <- self$aesthetics$stubLength %||% 0.02
        theta     <- angle_deg * pi / 180
        dx        <- L_npc * cos(theta)
        dy        <- L_npc * sin(theta)

        xs   <- vapply(self$stubs, `[[`, numeric(1),   "x")
        y0s  <- vapply(self$stubs, `[[`, numeric(1),   "y0")
        y1s  <- vapply(self$stubs, `[[`, numeric(1),   "y1")
        dirs <- vapply(self$stubs, `[[`, integer(1),   "dir")
        cols <- vapply(self$stubs, `[[`, character(1), "color")
        wds  <- vapply(self$stubs, `[[`, numeric(1),   "width")

        grid::grid.segments(
          x0 = grid::unit(xs,  "npc"),
          y0 = grid::unit(y0s, "npc"),
          x1 = grid::unit(xs,  "npc"),
          y1 = grid::unit(y1s, "npc"),
          gp = grid::gpar(
            col = grDevices::adjustcolor(cols, alpha.f = 0.5),
            lwd = wds
          )
        )

        grid::grid.segments(
          x0 = grid::unit(xs,             "npc"),
          y0 = grid::unit(y1s,            "npc"),
          x1 = grid::unit(xs + dirs * dx, "npc"),
          y1 = grid::unit(y1s +       dy, "npc"),
          gp = grid::gpar(
            col = grDevices::adjustcolor(cols, alpha.f = 0.5),
            lwd = wds
          ),
          arrow = grid::arrow(type = "open", length = grid::unit(1, "mm"))
        )

        partners <- gsub("chr", "",
                         vapply(self$stubs, `[[`, character(1), "partner"))
        grid::grid.text(
          label = partners,
          x = grid::unit(xs,         "npc"),
          y = grid::unit(y1s + 0.01, "npc"),
          just = c("center", "bottom"),
          gp   = grid::gpar(col = cols, fontsize = 7)
        )
      }
      invisible()
    }
  )
)

#' Bezier arch with vertical stems
#'
#' Same anchor / `t0`-`t1` rules as [seq_arc()], but draws stems from the
#' baseline up to the arch endpoints and renders stubs (with a partner
#' chromosome label) for half-visible links. Stub rendering is controlled
#' by `aes(plotStubs = TRUE)` (default).
#'
#' @inheritParams seq_arc
#' @return A `SeqArchR6` instance.
#' @export
seq_arch <- function(data = NULL, mapping = NULL,
                     t0 = NULL, t1 = NULL,
                     aesthetics = aes(), ...) {
  SeqArchR6$new(data = data, mapping = mapping,
                t0 = t0, t1 = t1,
                aesthetics = aesthetics, ...)
}
