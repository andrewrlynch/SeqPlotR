# Step ribbon element ----

# resolved per-row field with a default, recycled to n
.span_field <- function(resolved, name, default, n) {
  rep_len(resolved[[name]] %||% default, n)
}

# clip spans to windows; returns overlap hits with clipped genomic x0/x1 (NA spans skipped)
.span_window_hits <- function(data, x0, x1, windows) {
  ok <- which(!is.na(x0) & !is.na(x1))
  q <- GenomicRanges::GRanges(GenomicRanges::seqnames(data)[ok],
                              IRanges::IRanges(pmin(x0, x1)[ok], pmax(x0, x1)[ok]))
  ov <- GenomicRanges::findOverlaps(q, windows)
  qh <- ok[S4Vectors::queryHits(ov)]
  sh <- S4Vectors::subjectHits(ov)
  list(row = qh, window = sh,
       x0 = pmax(pmin(x0, x1)[qh], BiocGenerics::start(windows)[sh]),
       x1 = pmin(pmax(x0, x1)[qh], BiocGenerics::end(windows)[sh]))
}

#' SeqStepRibbon R6 class
#'
#' R6 generator backing [seq_step_ribbon()]. Inherits from [`SeqElementR6`].
#'
#' @keywords internal
SeqStepRibbonR6 <- R6::R6Class("SeqStepRibbon",
  inherit = SeqElementR6,
  public = list(
    #' @field y_fields Resolved fields in y axis units.
    y_fields        = c("y", "bar_lo", "bar_hi"),
    #' @field legend_arg The `legend` argument as given to the constructor.
    legend_arg      = NULL,
    #' @field legend_labels Named character relabelling auto keys by thickness.
    legend_labels   = NULL,
    #' @field legend_extra List of `LegendKey` appended to the auto keys.
    legend_extra    = NULL,
    #' @field legend_position Position of the element legend.
    legend_position = "canvas_margin",
    #' @field legend_side Side of the element legend.
    legend_side     = "right",

    #' @description Construct a SeqStepRibbonR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param legend `NULL` (auto keys), `FALSE`, a `SeqLegendSpec`, or a list
    #'   of `LegendKey`.
    #' @param legend_labels Named character relabelling auto keys.
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
      self$legend <- NULL
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0L) return(invisible(self$coordCanvas))

      r   <- self$resolved
      a   <- self$aesthetics
      n   <- length(eff_data)
      x0  <- .span_field(r, "x0", BiocGenerics::start(eff_data), n)
      x1  <- .span_field(r, "x1", BiocGenerics::end(eff_data), n)
      y   <- .span_field(r, "y", NA_real_, n)
      th  <- .span_field(r, "thickness", 1, n)
      col <- as.character(.span_field(r, "color", a$color %||% "grey30", n))
      alp <- .span_field(r, "alpha", a$alpha %||% 1, n)
      lab <- as.character(.span_field(r, "label", NA_character_, n))
      blo <- .span_field(r, "bar_lo", NA_real_, n)
      bhi <- .span_field(r, "bar_hi", NA_real_, n)
      pnt <- as.character(.span_field(r, "pointing", NA_character_, n))
      hh  <- a$half_height %||% 0.006

      hits <- .span_window_hits(eff_data, x0, x1, track_windows)
      keep <- !is.na(y[hits$row])
      for (w in unique(hits$window[keep])) {
        m  <- keep & hits$window == w
        i  <- hits$row[m]
        pm <- layout_track[[w]]
        cx <- seq_to_canvas(c(hits$x0[m], hits$x1[m], (hits$x0[m] + hits$x1[m]) / 2),
                            rep(y[i], 3L), pm, apply_transform = FALSE)$x
        k  <- length(i)
        cy <- function(v) seq_to_canvas(hits$x0[m], v, pm, apply_transform = FALSE)$y
        self$coordCanvas[[w]] <- data.frame(
          row = i, x0 = cx[seq_len(k)], x1 = cx[k + seq_len(k)],
          y0 = cy(y[i] - hh * th[i]), y1 = cy(y[i] + hh * th[i]),
          xm = cx[2L * k + seq_len(k)], ym = cy(y[i]),
          lo = cy(blo[i]), hi = cy(bhi[i]), pointing = pnt[i],
          color = col[i], alpha = alp[i], label = lab[i],
          stringsAsFactors = FALSE)
      }

      present <- sort(unique(th[!is.na(y)]))
      auto <- lapply(present, function(k)
        LegendKey(label = paste("thickness", k), shape = "line",
                  color = "grey40", lwd = 2 * k))
      names(auto) <- as.character(present)
      self$legend <- .glyph_element_legend(auto, self$legend_arg,
                                           self$legend_labels, self$legend_extra,
                                           self$legend_position, self$legend_side)
      invisible(self$coordCanvas)
    },

    #' @description Draw rects, interval lines, triangles and labels.
    draw = function() {
      a     <- self$aesthetics
      th_mm <- a$tri_half_mm   %||% 1.1
      tv_mm <- a$tri_height_mm %||% 1.9
      for (w in seq_along(self$coordCanvas)) {
        d <- self$coordCanvas[[w]]
        if (is.null(d) || nrow(d) == 0L) next
        fill <- .col_alpha(d$color, d$alpha)
        grid::grid.rect(
          x = grid::unit(d$x0, "npc"), y = grid::unit(d$y0, "npc"),
          width = grid::unit(d$x1 - d$x0, "npc"),
          height = grid::unit(d$y1 - d$y0, "npc"), just = c("left", "bottom"),
          gp = grid::gpar(fill = fill, col = NA),
          name = paste0("seqplotr.step_ribbon.rect.", w))
        b <- is.finite(d$lo) & is.finite(d$hi)
        if (any(b))
          grid::grid.segments(
            x0 = grid::unit(d$xm[b], "npc"), y0 = grid::unit(d$lo[b], "npc"),
            x1 = grid::unit(d$xm[b], "npc"), y1 = grid::unit(d$hi[b], "npc"),
            gp = grid::gpar(col = fill[b], lwd = a$line_width %||% 1),
            name = paste0("seqplotr.step_ribbon.bar.", w))
        t <- which(d$pointing %in% c("up", "down"))
        if (length(t)) {
          dir <- ifelse(d$pointing[t] == "up", 1, -1)
          grid::grid.polygon(
            x = grid::unit(rep(d$xm[t], each = 3L), "npc") +
              grid::unit(rep(c(-th_mm, th_mm, 0), length(t)), "mm"),
            y = grid::unit(rep(d$ym[t], each = 3L), "npc") +
              grid::unit(as.vector(rbind(0, 0, dir * tv_mm)), "mm"),
            id = rep(seq_along(t), each = 3L),
            gp = grid::gpar(fill = fill[t], col = NA),
            name = paste0("seqplotr.step_ribbon.triangle.", w))
        }
        l <- !is.na(d$label)
        if (any(l))
          grid::grid.text(
            d$label[l],
            x = grid::unit(d$x1[l], "npc") + grid::unit(a$label_mm %||% 1, "mm"),
            y = grid::unit(d$ym[l], "npc"), just = c("left", "centre"),
            gp = grid::gpar(cex = a$label_cex %||% 0.55, col = d$color[l]),
            name = paste0("seqplotr.step_ribbon.label.", w))
      }
      invisible()
    }
  )
)

#' Draw step ribbons with interval bars and direction triangles
#'
#' Draws each row as a thin horizontal band over its genomic span at height
#' `y`, with an optional vertical interval line from `bar_lo` to `bar_hi` at
#' the span midpoint, a fixed-size triangle at the midpoint pointing up or
#' down, and an optional label right of the span. Band thickness scales with
#' the mapped `thickness`; triangle size is in mm, independent of span width
#' and axis range. A row overlapping several windows is drawn in each,
#' clipped to the window.
#'
#' @param data Optional `GRanges`. Falls back to the parent track's data.
#' @param mapping Optional `SeqMap`. Fields: `x0` (default `start`), `x1`
#'   (default `end`), `y` (required), `thickness` (default 1), `color`
#'   (concrete colours, default `"grey30"`), `alpha` (default 1), `label`,
#'   `bar_lo`, `bar_hi`, `pointing` (`"up"`, `"down"` or `NA`). `y`, `bar_lo`
#'   and `bar_hi` are in y axis units.
#' @param aesthetics Optional `SeqAes`. `half_height` (axis units per
#'   thickness unit, default 0.006), `tri_half_mm` (1.1), `tri_height_mm`
#'   (1.9), `line_width` (1), `label_cex` (0.55), `label_mm` (1).
#' @param legend `NULL` for automatic keys (one line key per thickness value
#'   present, labelled `"thickness <k>"`), `FALSE` for none, or a
#'   `SeqLegendSpec` or list of `LegendKey` used as given.
#' @param legend_labels Named character vector relabelling the automatic keys,
#'   named by thickness value (e.g. `c("1" = "one copy")`).
#' @param legend_extra List of `LegendKey` appended to the legend keys.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the legend. Default
#'   `"canvas_margin"`, `"right"`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqStepRibbonR6` instance.
#' @examples
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1e5, 5e5), width = 2e5),
#'                              t = c(0.3, 0.6), k = c(1L, 3L))
#' seq_step_ribbon(gr, map(y = t, thickness = k, pointing = "up"))
#' @export
seq_step_ribbon <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                            legend = NULL, legend_labels = NULL,
                            legend_extra = NULL, show_legend = TRUE,
                            legend_position = "canvas_margin",
                            legend_side = "right", ...) {
  SeqStepRibbonR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                      legend = legend, legend_labels = legend_labels,
                      legend_extra = legend_extra, show_legend = show_legend,
                      legend_position = legend_position,
                      legend_side = legend_side, ...)
}
