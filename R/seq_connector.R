# SEQ CONNECTOR

# default line styles by key
.connector_styles <- function() {
  list(junction = list(col = "#3B82D6", lty = 1, lwd = 1.3),
       inferred = list(col = "#D33A2C", lty = "22", lwd = 1.3),
       unused   = list(col = "grey60", lty = 1, lwd = 1.1))
}

#' SeqConnector R6 class
#'
#' Internal R6 generator backing [seq_connector()]. Inherits from
#' [`SeqLinkR6`]. Draws a cubic Bezier from one edge of an interval to one
#' edge of another, on the same or different tracks.
#'
#' @keywords internal
SeqConnectorR6 <- R6::R6Class("SeqConnector",
  inherit = SeqLinkR6,
  public = list(
    #' @field styles Named list of `list(col, lty, lwd)` keyed by style.
    styles       = NULL,
    #' @field height Block height as a fraction of the panel inner height.
    height       = 0.5,
    #' @field handle Handle length as a fraction of the track inner width.
    handle       = 0.04,
    #' @field bulge_height Same-track bulge as a fraction of the panel inner
    #'   height.
    bulge_height = 0.45,

    #' @description Construct a SeqConnectorR6.
    #' @param data Optional `GRanges` or `data.frame`.
    #' @param mapping Optional `SeqMap`.
    #' @param t0,t1 Track identifiers.
    #' @param aesthetics Optional `SeqAes`.
    #' @param styles Optional named list of `list(col, lty, lwd)`.
    #' @param height,handle,bulge_height Geometry fractions.
    #' @param layer `"above"` or `"below"`.
    #' @param legend Optional legend.
    #' @param show_legend Logical.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL, t0 = NULL, t1 = NULL,
                          aesthetics = aes(), styles = NULL, height = 0.5,
                          handle = 0.04, bulge_height = 0.45,
                          layer = c("above", "below"), legend = NULL,
                          show_legend = TRUE, ...) {
      if (!is.null(data) && !inherits(data, c("GRanges", "SeqMap")) &&
          !is.data.frame(data))
        stop("seq_connector requires GRanges or data.frame data.", call. = FALSE)
      for (nm in c("height", "handle", "bulge_height")) {
        v <- get(nm)
        if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
          stop("seq_connector ", nm, " must be a single number.", call. = FALSE)
      }
      if (!is.null(styles) && (!is.list(styles) || is.null(names(styles)) ||
                               !all(vapply(styles, is.list, logical(1)))))
        stop("seq_connector styles must be a named list of lists.", call. = FALSE)
      super$initialize(data = data, mapping = mapping, t0 = t0, t1 = t1,
                       aesthetics = aesthetics, legend = legend,
                       show_legend = show_legend, layer = match.arg(layer))
      self$styles       <- styles
      self$height       <- height
      self$handle       <- handle
      self$bulge_height <- bulge_height
    },

    #' @description Resolve both anchors, map the chosen edges to canvas npc
    #'   and store four Bezier control points per connector in
    #'   `self$coordCanvas`.
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

      self$resolve(track_data    = layout_t0[[1]]$track_data,
                   track_mapping = layout_t0[[1]]$track_mapping)
      self$coordCanvas <- list()
      self$auto_legend <- NULL
      if (is.null(self$anchor0_gr) || length(self$anchor0_gr) == 0L)
        return(invisible())

      r <- self$resolved
      n <- length(self$anchor0_gr)
      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      .side <- function(v, def) {
        v <- tolower(as.character(.vec(v %||% def, n)))
        v[is.na(v)] <- def
        if (!all(v %in% c("left", "right")))
          stop("seq_connector sides must be 'left' or 'right'.", call. = FALSE)
        v
      }
      x0  <- .vec(as.numeric(r$x0), n)
      x1  <- .vec(as.numeric(r$x1), n)
      x0e <- .vec(as.numeric(r$x0_end %||% r$x0), n)
      x1e <- .vec(as.numeric(r$x1_end %||% r$x1), n)
      s0  <- .side(r$side0, "right")
      s1  <- .side(r$side1, "left")
      bul <- .vec(as.logical(r$bulge %||% TRUE), n)
      sty <- if (is.null(r$style)) rep(NA_character_, n)
             else .vec(as.character(r$style), n)

      # first window hit per anchor
      .first_hit <- function(gr, win) {
        out <- rep(NA_integer_, n)
        ov <- suppressWarnings(GenomicRanges::findOverlaps(gr, win, select = "first"))
        out[!is.na(ov)] <- ov[!is.na(ov)]
        out
      }
      w0 <- .first_hit(self$anchor0_gr, windows_t0)
      w1 <- .first_hit(self$anchor1_gr, windows_t1)

      .x_npc <- function(v, pm) {
        xpr <- pm$xplot_range %||% pm$xscale
        u <- pmax(pmin((v - xpr[1]) / diff(xpr), 1), 0)
        pm$inner$x0 + u * (pm$inner$x1 - pm$inner$x0)
      }
      fallback <- list(col = self$aesthetics$color     %||% "#3B82D6",
                       lty = self$aesthetics$linetype  %||% 1,
                       lwd = self$aesthetics$linewidth %||% 1.3)
      all_sty <- utils::modifyList(.connector_styles(), self$styles %||% list())
      .style <- function(k) {
        if (is.na(k) || is.null(all_sty[[k]])) return(fallback)
        utils::modifyList(fallback, all_sty[[k]])
      }
      same <- identical(tid0, tid1)

      for (i in which(!is.na(w0) & !is.na(w1))) {
        pm0 <- layout_t0[[w0[i]]]
        pm1 <- layout_t1[[w1[i]]]
        xa <- .x_npc(if (s0[i] == "left") x0[i] else x0e[i], pm0)
        xb <- .x_npc(if (s1[i] == "left") x1[i] else x1e[i], pm1)
        d0 <- if (s0[i] == "left") -1 else 1
        d1 <- if (s1[i] == "left") -1 else 1
        ti <- pm0$track_inner %||% pm0$inner
        hl <- self$handle * (ti$x1 - ti$x0)
        yc0 <- (pm0$inner$y0 + pm0$inner$y1) / 2
        yc1 <- (pm1$inner$y0 + pm1$inner$y1) / 2
        if (same) {
          s   <- if (isFALSE(bul[i])) -1 else 1
          ya  <- yc0 + s * self$height * (pm0$inner$y1 - pm0$inner$y0) / 2
          yb  <- yc1 + s * self$height * (pm1$inner$y1 - pm1$inner$y0) / 2
          bh  <- self$bulge_height * (pm0$inner$y1 - pm0$inner$y0)
          x   <- c(xa, xa + d0 * 0.3 * hl, xb + d1 * 0.3 * hl, xb)
          y   <- c(ya, ya + s * bh, yb + s * bh, yb)
        } else {
          x <- c(xa, xa + d0 * hl, xb + d1 * hl, xb)
          y <- c(yc0, yc0, yc1, yc1)
        }
        st <- .style(sty[i])
        self$coordCanvas[[length(self$coordCanvas) + 1L]] <- list(
          i = i, x = x, y = y, col = st$col, lty = st$lty, lwd = st$lwd,
          style = sty[i]
        )
      }

      if (!is.null(r$style) && isTRUE(self$show_legend) &&
          length(self$coordCanvas) > 0L) {
        keys <- unique(vapply(self$coordCanvas, function(cc) cc$style,
                              character(1)))
        keys <- keys[!is.na(keys)]
        if (length(keys)) {
          self$auto_legend <- seq_legend(keys = lapply(keys, function(k) {
            st <- .style(k)
            LegendKey(label = k, color = st$col, lty = st$lty, shape = "-")
          }))
        }
      }
      invisible()
    },

    #' @description Draw each prepared connector with [grid::grid.bezier()].
    draw = function() {
      for (cc in self$coordCanvas) {
        grid::grid.bezier(grid::unit(cc$x, "npc"), grid::unit(cc$y, "npc"),
                          gp = grid::gpar(col = cc$col, lty = cc$lty,
                                          lwd = cc$lwd),
                          name = "seqconnector")
      }
      invisible()
    }
  )
)

#' Bezier connector between interval edges
#'
#' Draws a cubic Bezier from the left or right edge of one interval to the
#' left or right edge of another, with handles pointing away from each
#' edge. Intervals may sit on the same track or on different tracks. Use it
#' to join consecutive pieces of an assembly, a rearranged chromosome or a
#' read path.
#'
#' Same-track connectors leave from the top (or bottom) edge of blocks of
#' `height` centred in the panel and bulge up (or down); connectors between
#' tracks leave from the block centres with horizontal handles.
#'
#' `map()` fields: `x0`, `x1` (required: the left edges of interval 0 and
#' interval 1, in track coordinates), `x0_end`, `x1_end` (right edges;
#' default the left edges), `chrom0`, `chrom1` (window seqnames; required
#' for `data.frame` data), `side0`, `side1` (`"left"` or `"right"`; default
#' `"right"` and `"left"`), `style` (a key into `styles`) and `bulge`
#' (same-track only; `TRUE` or `NA` up, `FALSE` down). Numeric edges are used
#' exactly.
#'
#' Default styles: `junction` solid blue, `inferred` red dashed, `unused`
#' grey. Entries of `styles` override these by name. A missing or unknown
#' key is drawn with `aes(color, linetype, linewidth)` (defaults
#' `"#3B82D6"`, `1`, `1.3`).
#'
#' @param data Optional `GRanges` or `data.frame`.
#' @param mapping Optional [map()].
#' @param t0,t1 Track identifiers. Locked to the parent track when added
#'   inside a [seq_track()].
#' @param aesthetics Optional [aes()]: `color`, `linetype`, `linewidth` for
#'   unstyled connectors.
#' @param styles Optional named list of `list(col, lty, lwd)`.
#' @param height Block height as a fraction of the panel inner height
#'   (default `0.5`).
#' @param handle Handle length as a fraction of the track inner width
#'   (default `0.04`; same-track handles use 30 percent of it).
#' @param bulge_height Same-track bulge as a fraction of the panel inner
#'   height (default `0.45`).
#' @param layer `"above"` (default) or `"below"` the track elements.
#' @param legend Optional `LegendKey`, list of keys or [seq_legend()].
#' @param show_legend Logical. When `style` is mapped, an automatic legend
#'   lists the keys present.
#' @param ... Reserved.
#' @return A `SeqConnectorR6` instance.
#' @examples
#' w <- GenomicRanges::GRanges("chrA", IRanges::IRanges(0, 100))
#' df <- data.frame(chrom0 = "chrA", x0 = 10, x0_end = 30,
#'                  chrom1 = "chrA", x1 = 60, x1_end = 90)
#' p <- seq_plot() %+%
#'   seq_track(track_id = "a", windows = w) %+%
#'   seq_connector(df, map(x0 = x0, x0_end = x0_end, chrom0 = chrom0,
#'                         x1 = x1, x1_end = x1_end, chrom1 = chrom1),
#'                 t0 = "a", t1 = "a")
#' @export
seq_connector <- function(data = NULL, mapping = NULL, t0 = NULL, t1 = NULL,
                          aesthetics = aes(), styles = NULL, height = 0.5,
                          handle = 0.04, bulge_height = 0.45,
                          layer = c("above", "below"), legend = NULL,
                          show_legend = TRUE, ...) {
  SeqConnectorR6$new(data = data, mapping = mapping, t0 = t0, t1 = t1,
                     aesthetics = aesthetics, styles = styles, height = height,
                     handle = handle, bulge_height = bulge_height,
                     layer = layer, legend = legend,
                     show_legend = show_legend, ...)
}
