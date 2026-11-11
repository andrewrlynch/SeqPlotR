# Paired ribbon element ----

#' SeqPairedRibbon R6 class
#'
#' R6 generator backing [seq_paired_ribbon()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqPairedRibbonR6 <- R6::R6Class("SeqPairedRibbon",
  inherit = SeqElementR6,
  public = list(
    #' @field y_fields Resolved fields in y axis units.
    y_fields        = c("y_up", "y_down"),
    #' @field legend_arg The `legend` argument as given to the constructor.
    legend_arg      = NULL,
    #' @field legend_labels Character relabelling the two auto keys.
    legend_labels   = NULL,
    #' @field legend_extra List of `LegendKey` appended to the auto keys.
    legend_extra    = NULL,
    #' @field legend_position Position of the element legend.
    legend_position = "canvas_margin",
    #' @field legend_side Side of the element legend.
    legend_side     = "right",

    #' @description Construct a SeqPairedRibbonR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param legend `NULL` (auto keys), `FALSE`, a `SeqLegendSpec`, or a list
    #'   of `LegendKey`.
    #' @param legend_labels Character relabelling the two auto keys.
    #' @param legend_extra List of `LegendKey` appended to the legend keys.
    #' @param legend_position,legend_side Placement of the element legend.
    #' @param show_legend Logical.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          legend = NULL, legend_labels = NULL,
                          legend_extra = NULL, legend_position = "canvas_margin",
                          legend_side = "right", show_legend = TRUE, ...) {
      super$initialize(data, mapping, aesthetics, legend = NULL,
                       show_legend = show_legend)
      self$legend_arg      <- legend
      self$legend_labels   <- legend_labels
      self$legend_extra    <- legend_extra
      self$legend_position <- legend_position
      self$legend_side     <- legend_side
    },

    #' @description Resolve mappings, clip spans to windows and fill
    #'   `coordCanvas` with one data.frame per window.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    #' @return `coordCanvas`, invisibly.
    prep = function(layout_track, track_windows) {
      self$resolve(track_data    = layout_track[[1]]$track_data,
                   track_mapping = layout_track[[1]]$track_mapping)
      self$coordCanvas <- vector("list", length(track_windows))
      a <- self$aesthetics
      fill_up   <- a$fill_up   %||% "#D62728"
      fill_down <- a$fill_down %||% "#1F77B4"
      ll <- self$legend_labels
      if (length(ll) && is.null(names(ll))) names(ll) <- c("up", "down")[seq_along(ll)]
      self$legend <- .glyph_element_legend(
        list(up = LegendKey(label = "first value", shape = "rect",
                            color = fill_up, fill = fill_up),
             down = LegendKey(label = "second value", shape = "rect",
                              color = fill_down, fill = fill_down)),
        self$legend_arg, ll, self$legend_extra,
        self$legend_position, self$legend_side)

      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0L) return(invisible(self$coordCanvas))
      r  <- self$resolved
      n  <- length(eff_data)
      x0 <- .span_field(r, "x0", BiocGenerics::start(eff_data), n)
      x1 <- .span_field(r, "x1", BiocGenerics::end(eff_data), n)
      yu <- .span_field(r, "y_up", NA_real_, n)
      yd <- .span_field(r, "y_down", NA_real_, n)

      hits <- .span_window_hits(eff_data, x0, x1, track_windows)
      for (w in unique(hits$window)) {
        m  <- hits$window == w
        i  <- hits$row[m]
        pm <- layout_track[[w]]
        t  <- (a$thickness %||% 0.012) * pm$yscale[2]
        k  <- length(i)
        cx <- seq_to_canvas(c(hits$x0[m], hits$x1[m]), rep(0, 2L * k), pm,
                            apply_transform = FALSE)$x
        cy <- function(v) seq_to_canvas(hits$x0[m], v, pm, apply_transform = FALSE)$y
        self$coordCanvas[[w]] <- data.frame(
          row = i, x0 = cx[seq_len(k)], x1 = cx[k + seq_len(k)],
          up0 = cy(yu[i]), up1 = cy(yu[i] + t),
          down0 = cy(yd[i] - t), down1 = cy(yd[i]))
      }
      invisible(self$coordCanvas)
    },

    #' @description Draw all second-value rects, then all first-value rects.
    draw = function() {
      a <- self$aesthetics
      parts <- list(down = c("down0", "down1", a$fill_down %||% "#1F77B4"),
                    up   = c("up0",   "up1",   a$fill_up   %||% "#D62728"))
      for (p in names(parts)) {
        for (w in seq_along(self$coordCanvas)) {
          d <- self$coordCanvas[[w]]
          if (is.null(d) || nrow(d) == 0L) next
          y0 <- d[[parts[[p]][1]]]
          y1 <- d[[parts[[p]][2]]]
          ok <- is.finite(y0) & is.finite(y1)
          if (!any(ok)) next
          grid::grid.rect(
            x = grid::unit(d$x0[ok], "npc"), y = grid::unit(y0[ok], "npc"),
            width = grid::unit(d$x1[ok] - d$x0[ok], "npc"),
            height = grid::unit(y1[ok] - y0[ok], "npc"),
            just = c("left", "bottom"),
            gp = grid::gpar(fill = parts[[p]][3], col = NA, alpha = a$alpha),
            name = paste0("seqplotr.paired_ribbon.", p, ".", w))
        }
      }
      invisible()
    }
  )
)

#' Draw paired ribbons for two values per span
#'
#' Draws two thin bands per genomic span: one from the first value `y_up`
#' upward and one from the second value `y_down` downward, so equal values
#' touch rather than overlap. All second-value bands are drawn before all
#' first-value bands. Band thickness is a fraction of the top of the track's
#' y range, so it is the same in axis units for every row.
#'
#' @param data Optional `GRanges`. Falls back to the parent track's data.
#' @param mapping Optional `SeqMap`. Fields: `x0` (default `start`), `x1`
#'   (default `end`), `y_up`, `y_down` (y axis units).
#' @param aesthetics Optional `SeqAes`. `fill_up` (default `"#D62728"`),
#'   `fill_down` (`"#1F77B4"`), `thickness` (fraction of the y range top,
#'   default 0.012), `alpha`.
#' @param legend `NULL` for two automatic rect keys (`"first value"`,
#'   `"second value"`), `FALSE` for none, or a `SeqLegendSpec` or list of
#'   `LegendKey` used as given.
#' @param legend_labels Character relabelling the automatic keys, named `up`
#'   and `down` or given in that order.
#' @param legend_extra List of `LegendKey` appended to the legend keys.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the legend. Default
#'   `"canvas_margin"`, `"right"`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqPairedRibbonR6` instance.
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 5e5), width = 5e5),
#'                              a = c(2, 3), b = c(1, 1))
#' seq_paired_ribbon(gr, map(y_up = a, y_down = b),
#'                   legend_labels = c(up = "A", down = "B"))
#' @export
seq_paired_ribbon <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                              legend = NULL, legend_labels = NULL,
                              legend_extra = NULL, show_legend = TRUE,
                              legend_position = "canvas_margin",
                              legend_side = "right", ...) {
  SeqPairedRibbonR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                        legend = legend, legend_labels = legend_labels,
                        legend_extra = legend_extra, show_legend = show_legend,
                        legend_position = legend_position,
                        legend_side = legend_side, ...)
}
