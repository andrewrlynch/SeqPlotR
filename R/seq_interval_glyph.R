# INTERVAL GLYPH

.glyph_types <- c("point", "latest", "earliest", "interval", "untimed")
.glyph_default_labels <- c(point = "point estimate", latest = "upper bound",
                           earliest = "lower bound", interval = "interval",
                           untimed = "no estimate")

# R6 generator backing seq_interval_glyph(); one glyph per row at (x, y) with bars from bar_lo to bar_hi.
SeqIntervalGlyphR6 <- R6::R6Class("SeqIntervalGlyph",
  inherit = SeqElementR6,
  public = list(
    legend_arg      = NULL,
    legend_labels   = NULL,
    legend_all      = FALSE,
    legend_extra    = NULL,
    legend_position = "canvas_margin",
    legend_side     = "right",
    y_fields        = c("y", "bar_lo", "bar_hi"),

    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          legend = NULL, legend_labels = NULL,
                          legend_all = FALSE, legend_extra = NULL,
                          show_legend = TRUE,
                          legend_position = "canvas_margin",
                          legend_side = "right", ...) {
      super$initialize(data, mapping, aesthetics, show_legend = show_legend)
      self$legend_arg      <- legend
      self$legend_labels   <- legend_labels
      self$legend_all      <- isTRUE(legend_all)
      self$legend_extra    <- legend_extra
      self$legend_position <- legend_position
      self$legend_side     <- legend_side
    },

    prep = function(layout_track, track_windows) {
      self$resolve(track_data    = layout_track[[1]]$track_data,
                   track_mapping = layout_track[[1]]$track_mapping)
      self$coordCanvas <- vector("list", length(track_windows))
      self$legend <- NULL
      eff <- self$resolved$.data
      n <- if (is.null(eff)) 0L else length(eff)
      if (n > 0L && !inherits(eff, "GRanges"))
        stop("seq_interval_glyph needs GRanges data.", call. = FALSE)

      r <- self$resolved
      a <- self$aesthetics
      .col <- function(v, default) {
        if (is.null(v)) return(rep(default, n))
        rep_len(v, n)
      }
      x     <- as.numeric(.col(r$x %||% if (n) BiocGenerics::start(eff), NA))
      y     <- as.numeric(.col(r$y, NA))
      lo    <- as.numeric(.col(r$bar_lo, NA))
      hi    <- as.numeric(.col(r$bar_hi, NA))
      glyph <- as.character(.col(r$glyph, "point"))
      glyph[is.na(glyph)] <- "point"
      col   <- .col(r$color, a$color %||% "#1C1B1A")
      if (!.looks_like_color(col))
        col <- .auto_scale_colors(col)$colors
      alpha <- as.numeric(.col(r$alpha, a$alpha %||% 1))
      label <- as.character(.col(r$label, NA))
      caut  <- as.logical(.col(r$caution, FALSE))
      caut[is.na(caut)] <- FALSE

      bad <- setdiff(unique(glyph), c(.glyph_types, "hidden"))
      if (length(bad))
        stop("seq_interval_glyph: unknown glyph value(s) ",
             paste(bad, collapse = ", "), ".", call. = FALSE)

      keep <- which(glyph != "hidden" & !is.na(x))
      present <- unique(glyph[keep])
      self$legend <- .glyph_element_legend(
        auto_keys = seq_glyph_legend_keys(present, labels = self$legend_labels,
                                          all = self$legend_all),
        legend_arg = self$legend_arg, legend_extra = self$legend_extra,
        position = self$legend_position, side = self$legend_side)
      if (!length(keep)) return(invisible())

      gr <- GenomicRanges::GRanges(
        as.character(GenomicRanges::seqnames(eff))[keep],
        IRanges::IRanges(start = x[keep], width = 1L))
      ov <- GenomicRanges::findOverlaps(gr, track_windows, ignore.strand = TRUE)
      if (!length(ov)) return(invisible())
      qh <- keep[S4Vectors::queryHits(ov)]
      sh <- S4Vectors::subjectHits(ov)
      for (w in unique(sh)) {
        i <- qh[sh == w]
        pm <- layout_track[[w]]
        p  <- seq_to_canvas(x[i], y[i], pm, apply_transform = FALSE)
        self$coordCanvas[[w]] <- data.frame(
          row = i, glyph = glyph[i], x = p$x, y = p$y,
          lo = seq_to_canvas(x[i], lo[i], pm, apply_transform = FALSE)$y,
          hi = seq_to_canvas(x[i], hi[i], pm, apply_transform = FALSE)$y,
          color = as.character(col[i]), alpha = alpha[i], label = label[i],
          caution = caut[i], stringsAsFactors = FALSE)
      }
      invisible()
    },

    draw = function() {
      a <- self$aesthetics
      for (w in seq_along(self$coordCanvas)) {
        d <- self$coordCanvas[[w]]
        if (is.null(d) || !nrow(d)) next
        .draw_interval_glyphs(d$x, d$y, d$lo, d$hi, d$glyph, d$color, d$alpha,
                              size = a$size %||% 2,
                              bar_lwd = a$bar_linewidth %||% 4,
                              bar_alpha = a$bar_alpha %||% 1,
                              cap_mm = a$cap_mm %||% 2.5,
                              points = a$glyph_points %||% TRUE,
                              name = "seqplotr.interval_glyph",
                              suffix = as.character(w))
        # labels sit right of the glyph at y, or at the bar centre when y is NA
        ly <- ifelse(is.finite(d$y), d$y, (d$lo + d$hi) / 2)
        s <- which(!is.na(d$label) & is.finite(ly))
        if (length(s))
          grid::grid.text(d$label[s],
                          grid::unit(d$x[s], "npc") +
                            grid::unit(a$label_mm %||% 1.5, "mm"),
                          grid::unit(ly[s], "npc"), just = c("left", "center"),
                          name = paste0("seqplotr.interval_glyph.label.", w),
                          gp = grid::gpar(col = .col_alpha(d$color[s], 1),
                                          cex = a$label_cex %||% 0.6))
        # caution marks sit left of the glyph so they never cover the label
        s <- which(d$caution & is.finite(ly))
        cmm <- a$caution_mm %||% 1.6
        if (length(s))
          .draw_caution_marks(d$x[s], ly[s], mm = cmm,
                              dx_mm = -((a$cap_mm %||% 2.5) / 2 + cmm + 0.3),
                              dy_mm = -0.45 * cmm,
                              name = paste0("seqplotr.interval_glyph.caution.",
                                            w))
      }
      invisible()
    },

    .infer_scale_y = function() {
      r <- self$resolved
      v <- suppressWarnings(as.numeric(c(r$y, r$bar_lo, r$bar_hi)))
      v <- v[is.finite(v)]
      if (!length(v)) return(NULL)
      seq_scale_continuous(limits = range(v))
    }
  )
)

#' Draw interval glyphs: point estimates, one-sided bounds and intervals
#'
#' Draws one glyph per row at genomic position `x` and value `y`, with a
#' vertical bar from `bar_lo` to `bar_hi`. The glyph type is read per row
#' from the `glyph` mapping:
#'
#' * `"point"`: filled point at `y` and a solid bar from `bar_lo` to `bar_hi`.
#' * `"latest"`: cap at `y` and an open bar from `y` down to `bar_lo`.
#' * `"earliest"`: cap at `y` and an open bar from `y` up to `bar_hi`.
#' * `"interval"`: bracket from `bar_lo` to `bar_hi`, plus a point at `y`
#'   when `y` is not `NA`.
#' * `"untimed"`: nothing is drawn.
#' * `"hidden"`: the row is dropped.
#'
#' A row with `NA` `y` gets no point and no cap. All values are drawn as
#' given. `coordCanvas[[w]]` holds one data.frame per window with columns
#' `row, glyph, x, y, lo, hi, color, alpha, label, caution` (canvas npc).
#'
#' @param data `GRanges`. Falls back to the parent track's data.
#' @param mapping `SeqMap` from [map()]. Fields: `x` (default `start`), `y`,
#'   `bar_lo`, `bar_hi`, `glyph` (default `"point"`), `color`, `alpha`,
#'   `label` (printed right of the glyph when not `NA`) and `caution`
#'   (logical; draws a small caution triangle).
#' @param aesthetics `SeqAes` from [aes()]: `size` (point size in mm,
#'   default 2), `bar_linewidth` (4), `bar_alpha` (1; multiplies the row
#'   alpha for bars), `cap_mm` (2.5), `glyph_points` (`TRUE`), `label_cex`
#'   (0.6), `label_mm` (1.5), `caution_mm` (1.6), `color`, `alpha`.
#' @param legend `NULL` (keys for the glyph types present, from
#'   [seq_glyph_legend_keys()]), `FALSE` (no legend), a `SeqLegendSpec`, or a
#'   list of `LegendKey` objects used as given.
#' @param legend_labels Named character relabelling the auto keys by glyph
#'   type, e.g. `c(point = "estimate")`.
#' @param legend_all Logical. Show keys for all five drawn glyph types.
#' @param legend_extra Optional `LegendKey` or list of them appended to the
#'   legend.
#' @param show_legend Logical. `FALSE` suppresses the legend.
#' @param legend_position,legend_side Placement of the built legend.
#'   Default `"canvas_margin"`, `"right"`.
#' @param ... Unused.
#' @return A `SeqIntervalGlyph` element.
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(c(2e5, 5e5, 8e5), width = 1),
#'               y = c(0.4, 0.6, NA), lo = c(0.3, 0.2, 0.5),
#'               hi = c(0.5, 0.6, 0.9),
#'               g = c("point", "latest", "interval"))
#' el <- seq_interval_glyph(gr, map(y = y, bar_lo = lo, bar_hi = hi,
#'                                  glyph = g))
#' @export
seq_interval_glyph <- function(data = NULL, mapping = NULL,
                               aesthetics = aes(), legend = NULL,
                               legend_labels = NULL, legend_all = FALSE,
                               legend_extra = NULL, show_legend = TRUE,
                               legend_position = "canvas_margin",
                               legend_side = "right", ...) {
  SeqIntervalGlyphR6$new(data = data, mapping = mapping,
                         aesthetics = aesthetics, legend = legend,
                         legend_labels = legend_labels,
                         legend_all = legend_all, legend_extra = legend_extra,
                         show_legend = show_legend,
                         legend_position = legend_position,
                         legend_side = legend_side, ...)
}

#' Legend keys for interval glyph types
#'
#' Builds one `LegendKey` per glyph type, each drawing a miniature of the
#' glyph through the key's `draw` function. Use it to label
#' [seq_interval_glyph()] or [seq_cross_arch()] glyphs in a legend.
#'
#' @param glyphs Character vector of glyph types (`"point"`, `"latest"`,
#'   `"earliest"`, `"interval"`, `"untimed"`); others, including `"hidden"`,
#'   are ignored. Keys follow that fixed order.
#' @param labels Optional named character overriding the default labels
#'   (`point = "point estimate"`, `latest = "upper bound"`,
#'   `earliest = "lower bound"`, `interval = "interval"`,
#'   `untimed = "no estimate"`).
#' @param all Logical. Return keys for all five types regardless of
#'   `glyphs`.
#' @param color Glyph colour in the keys.
#' @return A named list of `LegendKey` objects (names are glyph types).
#' @examples
#' keys <- seq_glyph_legend_keys(c("point", "latest"))
#' names(keys)
#' @export
seq_glyph_legend_keys <- function(glyphs = character(), labels = NULL,
                                  all = FALSE, color = "grey30") {
  types <- if (isTRUE(all)) .glyph_types
           else .glyph_types[.glyph_types %in% as.character(glyphs)]
  lab <- .glyph_default_labels
  if (length(labels)) {
    nm <- intersect(names(labels), .glyph_types)
    lab[nm] <- as.character(labels[nm])
  }
  keys <- lapply(types, function(g) {
    LegendKey(label = unname(lab[[g]]), color = color, shape = "custom",
              glyph = g, draw = function(x0, x1, y, height) {
                xm <- (x0 + x1) / 2
                h  <- max(height, 0.004) * 0.5
                if (g == "untimed") {
                  grid::grid.segments(xm, y - h, xm, y + h,
                                      gp = grid::gpar(col = .col_alpha(color, 0.5),
                                                      lwd = 1))
                  return(invisible())
                }
                gy <- switch(g, latest = y + h * 0.6, earliest = y - h * 0.6,
                             interval = NA, y)
                .draw_interval_glyphs(xm, gy, y - h, y + h, g, color,
                                      size = 1.6, bar_lwd = 3, bar_alpha = 0.5,
                                      cap_mm = 1.8, name = NULL)
              })
  })
  stats::setNames(keys, types)
}
