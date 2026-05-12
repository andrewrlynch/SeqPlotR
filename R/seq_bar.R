# ── SeqBarR6 ──────────────────────────────────────────────────────────────────
#
# Composite element: one filled rectangle per genomic range, spanning
# start..end on the x-axis with height y. Optional stacking by group.

#' SeqBar R6 class
#'
#' Internal R6 generator backing [seq_bar()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqBarR6 <- R6::R6Class("SeqBar",
  inherit = SeqElementR6,
  public = list(
    #' @description Construct a SeqBarR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Recognised: `x`, `y`, `group`, `fill`.
    #' @param aesthetics Optional `SeqAes`. Supports `fill`, `color`,
    #'   `linewidth`, `alpha`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
      super$initialize(data, mapping, aesthetics, ...)
    },

    #' @description Resolve mappings, clip to windows, transform to npc. When
    #'   a `group` mapping is present, bars stack at identical x positions.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())

      n        <- length(eff_data)
      y_vals   <- self$resolved$y     %||% rep(1, n)
      g_vals   <- self$resolved$group %||% rep("default", n)
      g_vals   <- as.character(g_vals)
      g_levels <- unique(g_vals)
      has_group <- !is.null(self$resolved$group)

      fill_map <- self$resolved$fill
      const_fill <- self$aesthetics$fill

      # Scale fill to hex colors; set auto_legend only when no manual legend
      if (!is.null(fill_map)) {
        fill_nm     <- .map_col_name(self$resolved$.mapping, "fill")
        fill_result <- .auto_scale_colors(fill_map, col_name = fill_nm,
                                          aes_name = "fill")
        fill_map <- fill_result$colors %||% as.character(fill_map)
        if (is.null(self$legend) && isTRUE(self$show_legend))
          self$auto_legend <- fill_result$legend
      }

      start_pos <- BiocGenerics::start(eff_data)
      end_pos   <- BiocGenerics::end(eff_data)

      self$coordCanvas <- vector("list", length(track_windows))
      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      if (length(ov) == 0) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      palette <- flexoki_palette(max(length(g_levels), 1))
      names(palette) <- g_levels

      for (w in unique(sh)) {
        mask <- sh == w
        idx  <- qh[mask]
        pm   <- layout_track[[w]]

        x0 <- start_pos[idx]
        x1 <- end_pos[idx]
        clip <- .clip_to_windows(x0, x1, pm$xscale)
        if (length(clip$x0) == 0) next

        keep <- clip$mask
        idx_k <- idx[keep]

        df <- data.frame(
          x0    = clip$x0,
          x1    = clip$x1,
          x     = (clip$x0 + clip$x1) / 2,
          y     = y_vals[idx_k],
          group = g_vals[idx_k],
          stringsAsFactors = FALSE
        )
        if (!is.null(fill_map)) {
          df$fill_per <- as.character(fill_map[idx_k])
        } else {
          df$fill_per <- NA_character_
        }

        df$group <- factor(df$group, levels = g_levels)
        df <- df[order(df$x, df$group), , drop = FALSE]

        if (has_group) {
          df$y0 <- 0
          df$y1 <- 0
          for (xv in unique(df$x)) {
            ii <- which(df$x == xv)
            heights <- df$y[ii]
            df$y0[ii] <- cumsum(c(0, utils::head(heights, -1)))
            df$y1[ii] <- cumsum(heights)
          }
        } else {
          df$y0 <- 0
          df$y1 <- df$y
        }

        u0 <- pmax(pmin((df$x0 - pm$xscale[1]) / diff(pm$xscale), 1), 0)
        u1 <- pmax(pmin((df$x1 - pm$xscale[1]) / diff(pm$xscale), 1), 0)
        v0 <- pmax(pmin((df$y0 - pm$yscale[1]) / diff(pm$yscale), 1), 0)
        v1 <- pmax(pmin((df$y1 - pm$yscale[1]) / diff(pm$yscale), 1), 0)

        xw <- pm$inner$x1 - pm$inner$x0
        yw <- pm$inner$y1 - pm$inner$y0
        x0c <- pm$inner$x0 + u0 * xw
        x1c <- pm$inner$x0 + u1 * xw
        y0c <- pm$inner$y0 + v0 * yw
        y1c <- pm$inner$y0 + v1 * yw

        if (!is.null(fill_map)) {
          fill_vec <- df$fill_per
        } else if (!is.null(const_fill)) {
          fill_vec <- rep(const_fill, nrow(df))
        } else if (has_group) {
          fill_vec <- unname(palette[as.character(df$group)])
        } else {
          fill_vec <- rep("grey60", nrow(df))
        }

        self$coordCanvas[[w]] <- data.frame(
          x0   = x0c,
          x1   = x1c,
          y0   = y0c,
          y1   = y1c,
          fill = fill_vec,
          stringsAsFactors = FALSE
        )
      }
      invisible()
    },

    #' @description Draw bars via `grid::grid.rect()`.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      color <- self$aesthetics$color %||% self$aesthetics$col %||% NA
      lwd   <- self$aesthetics$linewidth %||% self$aesthetics$lwd %||% 1
      alpha <- self$aesthetics$alpha
      for (coords in self$coordCanvas) {
        if (!is.data.frame(coords) || nrow(coords) == 0) next
        grid::grid.rect(
          x      = grid::unit((coords$x0 + coords$x1) / 2, "npc"),
          y      = grid::unit((coords$y0 + coords$y1) / 2, "npc"),
          width  = grid::unit(abs(coords$x1 - coords$x0), "npc"),
          height = grid::unit(abs(coords$y1 - coords$y0), "npc"),
          gp = grid::gpar(
            fill  = coords$fill,
            col   = color,
            lwd   = lwd,
            alpha = alpha
          )
        )
      }
    }
  )
)

#' Draw bars
#'
#' Composite element drawing one filled rectangle per genomic range, spanning
#' `start..end` on the x-axis with height `y`. Supply a `group` mapping to
#' stack bars vertically at identical x positions.
#'
#' @param data Optional `GRanges`. Falls back to the parent track's data.
#' @param mapping Optional `SeqMap` from [map()]. Recognised fields: `x`,
#'   `y` (bar height, default 1), `group` (stacking), `fill` (per-bar color).
#' @param aesthetics Optional `SeqAes` from [aes()]. Supports `fill`,
#'   `color`, `linewidth`, `alpha`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqBarR6` instance.
#' @examples
#' seq_bar(map(x = start, y = score))
#' @export
seq_bar <- function(data = NULL, mapping = NULL, aesthetics = aes(), ...) {
  SeqBarR6$new(data = data, mapping = mapping, aesthetics = aesthetics, ...)
}
