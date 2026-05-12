# ── SeqTileR6 ─────────────────────────────────────────────────────────────────
#
# Composite element: one filled shape per genomic range.
# - Unrotated (default): rectangle at (x0, y_idx-0.5) -- (x1, y_idx+0.5).
# - Rotated: diamond built via a linear coordinate transform in genomic space.
#   The diamond's four corners are then mapped to canvas npc through the
#   usual xscale / yscale pipeline.

#' SeqTile R6 class
#'
#' Internal R6 generator backing [seq_tile()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqTileR6 <- R6::R6Class("SeqTile",
  inherit = SeqElementR6,
  public = list(
    #' @field data2 Optional `GRanges` giving the y-axis genomic coordinates
    #'   for rotated mode (one range per row of the primary `data`).
    data2 = NULL,

    #' @description Construct a SeqTileR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Recognised: `x`, `y`, `fill`.
    #' @param aesthetics Optional `SeqAes`. `rotate` toggles diamond mode;
    #'   `fill` sets a constant fill color when no `fill` mapping is given.
    #' @param data2 Optional `GRanges` for rotated mode (y genomic ranges).
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          data2 = NULL, ...) {
      super$initialize(data, mapping, aesthetics, ...)
      self$data2 <- data2
    },

    #' @description Resolve mappings, clip to windows, transform to npc.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      rotate <- isTRUE(self$aesthetics$rotate)
      n      <- length(eff_data)

      fill_map   <- self$resolved$fill
      const_fill <- self$aesthetics$fill %||% "grey60"

      # Auto-scale fill when it is mapped but not already concrete colors
      fill_scaled <- NULL
      if (!is.null(fill_map) && is.null(self$legend) && isTRUE(self$show_legend)) {
        fill_nm <- .map_col_name(self$resolved$.mapping, "fill")
        fill_res <- .auto_scale_colors(fill_map, col_name = fill_nm,
                                        aes_name = "fill")
        fill_scaled   <- fill_res$colors
        self$auto_legend <- fill_res$legend
      } else if (!is.null(fill_map) && .looks_like_color(fill_map)) {
        fill_scaled <- as.character(fill_map)
      }

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        mask <- sh == w
        idx  <- qh[mask]
        if (length(idx) == 0) next
        pm <- layout_track[[w]]

        x0_g <- BiocGenerics::start(eff_data)[idx]
        x1_g <- BiocGenerics::end(eff_data)[idx]

        fill_vals <- if (!is.null(fill_scaled))
                       fill_scaled[idx]
                     else
                       rep(const_fill, length(idx))

        if (!rotate) {
          # When `data2` is supplied, each tile spans a genomic y-range.
          # Otherwise y is a unit-height row indexed by `resolved$y`.
          if (!is.null(self$data2) &&
              length(self$data2) == length(eff_data)) {
            y0_gu <- BiocGenerics::start(self$data2)[idx]
            y1_gu <- BiocGenerics::end(self$data2)[idx]
            j_seq <- as.character(GenomicRanges::seqnames(self$data2))[idx]
          } else {
            y_idx <- if (!is.null(self$resolved$y))
                       self$resolved$y[idx]
                     else
                       rep(1, length(idx))
            y0_gu <- y_idx - 0.5
            y1_gu <- y_idx + 0.5
            j_seq <- NULL
          }

          xpr <- pm$xplot_range %||% pm$xscale
          u0 <- pmax(pmin((x0_g  - xpr[1]) / diff(xpr), 1), 0)
          u1 <- pmax(pmin((x1_g  - xpr[1]) / diff(xpr), 1), 0)
          xw <- pm$inner$x1 - pm$inner$x0
          x0c <- pm$inner$x0 + u0 * xw
          x1c <- pm$inner$x0 + u1 * xw

          # Multi-y-window: route each tile into the y sub-panel matching
          # its j-bin seqname. Falls back to a single-band layout when
          # y_sub_panels is absent.
          y_subs <- pm$y_sub_panels
          if (!is.null(y_subs) && length(y_subs) > 0L && !is.null(j_seq)) {
            frames <- vector("list", length(y_subs))
            for (yk in seq_along(y_subs)) {
              sub  <- y_subs[[yk]]
              ymask <- j_seq == sub$seqname
              if (!any(ymask)) next
              ys <- sub$yplot_range %||% sub$yscale
              v0 <- pmax(pmin((y0_gu[ymask] - ys[1]) / diff(ys), 1), 0)
              v1 <- pmax(pmin((y1_gu[ymask] - ys[1]) / diff(ys), 1), 0)
              y0c <- sub$y0 + v0 * (sub$y1 - sub$y0)
              y1c <- sub$y0 + v1 * (sub$y1 - sub$y0)
              keep <- (x1_g[ymask] >= pm$xscale[1]) &
                      (x0_g[ymask] <= pm$xscale[2])
              if (!any(keep)) next
              frames[[yk]] <- data.frame(
                x0 = x0c[ymask][keep], x1 = x1c[ymask][keep],
                y0 = y0c[keep],         y1 = y1c[keep],
                fill = fill_vals[ymask][keep],
                stringsAsFactors = FALSE
              )
            }
            frames <- Filter(Negate(is.null), frames)
            df_w <- if (length(frames)) do.call(rbind, frames) else NULL
          } else {
            ypr <- pm$yplot_range %||% pm$yscale
            v0 <- pmax(pmin((y0_gu - ypr[1]) / diff(ypr), 1), 0)
            v1 <- pmax(pmin((y1_gu - ypr[1]) / diff(ypr), 1), 0)
            yw <- pm$inner$y1 - pm$inner$y0
            y0c <- pm$inner$y0 + v0 * yw
            y1c <- pm$inner$y0 + v1 * yw
            keep <- (x1_g >= pm$xscale[1]) & (x0_g <= pm$xscale[2])
            df_w <- data.frame(
              x0 = x0c[keep], x1 = x1c[keep],
              y0 = y0c[keep], y1 = y1c[keep],
              fill = fill_vals[keep],
              stringsAsFactors = FALSE
            )
          }
          # Apply axis flips by mirroring rectangle bounds around the
          # panel inner centre. Swap x0 <-> x1 / y0 <-> y1 to keep the
          # df_w$x0 < df_w$x1 invariant; this matters for grid::grid.rect.
          if (!is.null(df_w) && nrow(df_w) > 0L) {
            if (isTRUE(pm$flip_x)) {
              c_x <- pm$inner$x0 + pm$inner$x1
              tmp <- df_w$x0; df_w$x0 <- c_x - df_w$x1; df_w$x1 <- c_x - tmp
            }
            if (isTRUE(pm$flip_y)) {
              c_y <- pm$inner$y0 + pm$inner$y1
              tmp <- df_w$y0; df_w$y0 <- c_y - df_w$y1; df_w$y1 <- c_y - tmp
            }
          }
          self$coordCanvas[[w]] <- df_w
        } else {
          if (is.null(self$data2) || length(self$data2) != length(eff_data)) {
            stop("seq_tile(rotate = TRUE) requires a `data2` GRanges of the ",
                 "same length as `data`.", call. = FALSE)
          }
          y0_g <- BiocGenerics::start(self$data2)[idx]
          y1_g <- BiocGenerics::end(self$data2)[idx]

          # Genomic-space rotation. x_rot = (x + y) / 2 gives the
          # position midpoint; y_rot = (y - x) gives the full genomic
          # distance (bp). Using the full distance (not halved) means
          # the panel's y-axis reads as real interaction distance and
          # `max_dist`/`scale_y` limits map directly to bp.
          x0r <- (x0_g + y0_g) / 2
          x1r <- (x1_g + y1_g) / 2
          y0r <-  y0_g - x1_g                  # bottom corner y (bp)
          y1r <-  y1_g - x0_g                  # top    corner y (bp)
          xcr <- (x0r + x1r) / 2
          ycr <- (y0r + y1r) / 2

          xw <- pm$inner$x1 - pm$inner$x0
          yw <- pm$inner$y1 - pm$inner$y0
          # Unclamped NPC transforms — rotated tiles straddling the
          # panel edge must keep their true corner positions so the
          # diamond shape survives. Clipping happens at draw() time via
          # a clipped viewport.
          xpr <- pm$xplot_range %||% pm$xscale
          ypr <- pm$yplot_range %||% pm$yscale
          to_x_npc <- function(v)
            pm$inner$x0 + (v - xpr[1]) / diff(xpr) * xw
          to_y_npc <- function(v)
            pm$inner$y0 + (v - ypr[1]) / diff(ypr) * yw

          # Skip tiles whose bounding box is entirely outside the panel
          # in either axis — they contribute nothing visible.
          vis <- (x1r >= pm$xscale[1]) & (x0r <= pm$xscale[2]) &
                 (y1r >= pm$yscale[1]) & (y0r <= pm$yscale[2])
          x0r <- x0r[vis]; x1r <- x1r[vis]; xcr <- xcr[vis]
          y0r <- y0r[vis]; y1r <- y1r[vis]; ycr <- ycr[vis]
          fill_vals <- fill_vals[vis]

          if (length(x0r) == 0) next

          # Four diamond corners per tile: left, bottom, right, top.
          poly_x <- c(rbind(to_x_npc(x0r), to_x_npc(xcr),
                            to_x_npc(x1r), to_x_npc(xcr)))
          poly_y <- c(rbind(to_y_npc(ycr), to_y_npc(y0r),
                            to_y_npc(ycr), to_y_npc(y1r)))

          # Apply axis flips by mirroring corner positions around the
          # panel inner centre.
          if (isTRUE(pm$flip_x))
            poly_x <- (pm$inner$x0 + pm$inner$x1) - poly_x
          if (isTRUE(pm$flip_y))
            poly_y <- (pm$inner$y0 + pm$inner$y1) - poly_y

          self$coordCanvas[[w]] <- list(
            rotated    = TRUE,
            x          = poly_x,
            y          = poly_y,
            id.lengths = rep(4L, length(x0r)),
            fill       = fill_vals,
            # Inner-panel NPC bounds for clip viewport at draw time.
            clip       = list(x0 = pm$inner$x0, x1 = pm$inner$x1,
                              y0 = pm$inner$y0, y1 = pm$inner$y1)
          )
        }
      }
      invisible()
    },

    #' @description Draw tiles. Rectangles via `grid::grid.rect()` for the
    #'   unrotated mode; diamond polygons via
    #'   `grid::grid.polygon(id.lengths = ...)` for the rotated mode.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      rotate <- isTRUE(self$aesthetics$rotate)
      border <- self$aesthetics$color %||% self$aesthetics$col %||% NA
      lwd    <- self$aesthetics$linewidth %||% self$aesthetics$lwd %||% 0.5
      alpha  <- self$aesthetics$alpha

      for (coords in self$coordCanvas) {
        if (is.null(coords)) next

        if (!rotate) {
          if (!is.data.frame(coords) || nrow(coords) == 0) next
          grid::grid.rect(
            x      = grid::unit((coords$x0 + coords$x1) / 2, "npc"),
            y      = grid::unit((coords$y0 + coords$y1) / 2, "npc"),
            width  = grid::unit(abs(coords$x1 - coords$x0), "npc"),
            height = grid::unit(abs(coords$y1 - coords$y0), "npc"),
            gp = grid::gpar(
              fill  = coords$fill,
              col   = border,
              lwd   = lwd,
              alpha = alpha
            )
          )
        } else {
          if (length(coords$x) == 0) next
          # Push a clipped viewport exactly covering the inner panel
          # so diamond polygons that straddle the edge are clipped to
          # the visible area instead of being drawn outside it.
          cl  <- coords$clip
          vp  <- grid::viewport(
            x      = grid::unit(cl$x0, "npc"),
            y      = grid::unit(cl$y0, "npc"),
            width  = grid::unit(cl$x1 - cl$x0, "npc"),
            height = grid::unit(cl$y1 - cl$y0, "npc"),
            just   = c("left", "bottom"),
            xscale = c(cl$x0, cl$x1),
            yscale = c(cl$y0, cl$y1),
            clip   = "on"
          )
          grid::pushViewport(vp)
          grid::grid.polygon(
            x          = grid::unit(coords$x, "native"),
            y          = grid::unit(coords$y, "native"),
            id.lengths = coords$id.lengths,
            gp = grid::gpar(
              fill  = coords$fill,
              col   = border,
              lwd   = lwd,
              alpha = alpha
            )
          )
          grid::popViewport()
        }
      }
    }
  )
)

#' Draw tiles (rectangles or rotated diamonds)
#'
#' Composite element drawing one filled shape per genomic range. In the
#' default unrotated mode, each tile is a rectangle spanning `start..end`
#' on the x-axis at the row indicated by `y` (default 1). With
#' `aes(rotate = TRUE)` and a `data2` GRanges giving per-tile y-axis
#' coordinates, tiles are rendered as diamonds via a linear
#' coordinate transform in genomic space.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Recognised fields: `x`, `y`, `fill`.
#' @param aesthetics Optional `SeqAes`. `rotate` toggles diamond mode;
#'   `fill` sets a constant fill color.
#' @param data2 Optional `GRanges` providing y-axis genomic ranges for
#'   rotated mode. Must match `data` in length.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqTileR6` instance.
#' @examples
#' seq_tile(map(x = start, fill = color))
#' @export
seq_tile <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                     data2 = NULL, ...) {
  SeqTileR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                data2 = data2, ...)
}
