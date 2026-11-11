# LABELLED SPAN STRIP ELEMENT

# Coloured bars over genomic spans with a text label per span.
SeqSpanStripR6 <- R6::R6Class("SeqSpanStrip",
  inherit = SeqElementR6,
  public = list(
    legend_arg = NULL,
    legend_extra = NULL,
    legend_position = "canvas_margin",
    legend_side = "right",

    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          legend = NULL, legend_extra = NULL, show_legend = TRUE,
                          legend_position = "canvas_margin",
                          legend_side = "right", ...) {
      super$initialize(data, mapping, aesthetics, legend = NULL,
                       show_legend = show_legend)
      self$legend_arg <- legend
      self$legend_extra <- legend_extra
      self$legend_position <- legend_position
      self$legend_side <- legend_side
    },

    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      # no auto keys; legend and legend_extra come from the caller
      self$legend <- if (isTRUE(self$show_legend))
        .element_legend_spec(self$legend_arg, list(), self$legend_extra,
                             self$legend_position, self$legend_side)
      else NULL

      eff_data <- self$resolved$.data
      self$coordCanvas <- vector("list", length(track_windows))
      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())
      if (!inherits(eff_data, "GRanges"))
        stop("seq_span_strip: data must be a GRanges.", call. = FALSE)

      r   <- self$resolved
      n   <- length(eff_data)
      x0  <- r$x0 %||% BiocGenerics::start(eff_data)
      x1  <- r$x1 %||% BiocGenerics::end(eff_data)
      lab <- as.character(r$label %||% rep(NA_character_, n))
      fill <- as.character(r$fill %||% rep(NA_character_, n))
      fill[is.na(fill)] <- "grey70"
      a <- self$aesthetics
      alpha <- a$alpha %||% 1
      bar_h <- a$bar_height %||% 0.35
      bar_y <- a$bar_y %||% 0.5

      ov <- GenomicRanges::findOverlaps(eff_data, track_windows)
      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)
      for (w in unique(sh)) {
        idx <- qh[sh == w]
        pm  <- layout_track[[w]]
        xpr <- pm$xplot_range %||% pm$xscale
        xw  <- pm$inner$x1 - pm$inner$x0
        yw  <- pm$inner$y1 - pm$inner$y0
        u0 <- pmax(pmin((pmin(x0[idx], x1[idx]) - xpr[1]) / diff(xpr), 1), 0)
        u1 <- pmax(pmin((pmax(x0[idx], x1[idx]) - xpr[1]) / diff(xpr), 1), 0)
        cx0 <- pm$inner$x0 + u0 * xw
        cx1 <- pm$inner$x0 + u1 * xw
        cy0 <- pm$inner$y0 + (bar_y - bar_h / 2) * yw
        cy1 <- pm$inner$y0 + (bar_y + bar_h / 2) * yw
        # label anchor: above the bar top, or on the bar line when bar_height is 0
        yl <- if (identical(a$label_position, "inside") || bar_h == 0)
                pm$inner$y0 + bar_y * yw
              else cy1
        self$coordCanvas[[w]] <- data.frame(
          row = idx, x0 = cx0, x1 = cx1, y0 = cy0, y1 = cy1,
          xl = (cx0 + cx1) / 2, yl = yl, fill = fill[idx], label = lab[idx],
          stringsAsFactors = FALSE
        )
      }
      invisible()
    },

    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      a <- self$aesthetics
      alpha <- a$alpha %||% 1
      bar_h <- a$bar_height %||% 0.35
      cex   <- a$label_cex %||% 0.6
      above <- !identical(a$label_position, "inside") && bar_h > 0
      for (w in seq_along(self$coordCanvas)) {
        d <- self$coordCanvas[[w]]
        if (is.null(d) || nrow(d) == 0L) next
        if (bar_h > 0)
          grid::grid.rect(
            x = grid::unit(d$x0, "npc"), y = grid::unit(d$y0, "npc"),
            width = grid::unit(d$x1 - d$x0, "npc"),
            height = grid::unit(d$y1 - d$y0, "npc"),
            just = c("left", "bottom"),
            gp = grid::gpar(fill = .with_alpha(d$fill, alpha), col = NA),
            name = sprintf("seqplotr.span_strip.rect.%d", w)
          )
        ok <- !is.na(d$label)
        if (any(ok))
          grid::grid.text(
            d$label[ok], x = grid::unit(d$xl[ok], "npc"),
            y = grid::unit(d$yl[ok], "npc") +
                (if (above) grid::unit(0.6, "mm") else grid::unit(0, "mm")),
            just = c("center", if (above) "bottom" else "center"),
            gp = grid::gpar(cex = cex, col = "#1C1B1A"),
            name = sprintf("seqplotr.span_strip.text.%d", w)
          )
      }
      invisible()
    }
  )
)

#' Draw labelled genomic spans
#'
#' Draws a coloured bar over each genomic span with a text label at the span
#' centre. Spans are clipped to each window. With `bar_height = 0` only the
#' labels are drawn.
#'
#' @param data A `GRanges` with one row per span.
#' @param mapping A `SeqMap`. Fields: `x0` (default `start`), `x1` (default
#'   `end`), `label`, `fill` (concrete colours; `NA` gives `grey70`).
#' @param aesthetics A `SeqAes`. Recognised: `bar_height` (fraction of the
#'   panel height, default 0.35), `bar_y` (bar centre as a panel fraction,
#'   default 0.5), `label_cex` (0.6), `label_position` (`"above"` or
#'   `"inside"`), `alpha` (1).
#' @param legend `NULL` for no key, `FALSE` for none, or a `SeqLegendSpec` or
#'   list of `LegendKey` used as given.
#' @param legend_extra List of `LegendKey` appended to the legend.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the generated legend.
#' @param ... Unused.
#' @return A `SeqSpanStripR6` instance.
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 500), c(400, 900)),
#'   name = c("a", "b"), col = c("tomato", "steelblue"))
#' seq_span_strip(gr, map(label = name, fill = col))
#' @export
seq_span_strip <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                           legend = NULL, legend_extra = NULL,
                           show_legend = TRUE,
                           legend_position = "canvas_margin",
                           legend_side = "right", ...) {
  SeqSpanStripR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                     legend = legend, legend_extra = legend_extra,
                     show_legend = show_legend, legend_position = legend_position,
                     legend_side = legend_side, ...)
}
