# ── seq_ideogram — chromosome ideogram from cytogenetic bands ────────────────
#
# Ported from THEfunc `SeqIdeogram.R`. Non-centromeric bands render as filled
# rectangles; each pair of `acen` bands renders as two inward-pointing red
# triangles (the centromere). Two new arguments extend the original behaviour:
#
#   scope  — "window" (default) or "full" (whole-chromosome view with the
#            window highlighted)
#   style  — "block" (default) or "rounded" (telomere-end rounded caps)

#' Giemsa stain code → fill color
#'
#' Maps UCSC `gieStain` codes (`gneg`, `gpos25`..`gpos100`, `acen`, `stalk`,
#' `gvar`) to fill colors. Unknown codes fall through to `"#CCCCCC"`.
#'
#' @param stain Character vector of `gieStain` codes.
#' @return Character vector of fill colors, same length as `stain`.
#' @keywords internal
.ideogram_fill_colors <- function(stain) {
  vapply(stain, function(s) {
    if (is.na(s))                 "#CCCCCC"
    else if (s == "gneg")         "#FFFFFF"
    else if (startsWith(s, "gpos")) {
      pct <- 1 - suppressWarnings(as.numeric(sub("gpos", "", s))) / 100
      if (!is.finite(pct)) "#CCCCCC" else grDevices::grey(pct)
    }
    else if (s == "acen")         "#FF0000"
    else if (s == "stalk")        "#7AC0CF"
    else if (s == "gvar")         "#CCF5FF"
    else                          "#CCCCCC"
  }, character(1))
}

#' SeqIdeogram R6 class
#'
#' Internal R6 generator backing [seq_ideogram()]. Inherits from
#' [`SeqElementR6`]. Draws chromosome ideograms from a `GRanges` of
#' cytogenetic bands. Each band becomes a filled rectangle; consecutive
#' `acen` bands within a window render as two red triangles meeting at the
#' centromere.
#'
#' @keywords internal
SeqIdeogramR6 <- R6::R6Class("SeqIdeogram",
  inherit = SeqElementR6,
  public = list(
    #' @field centroPolys Per-window list of centromere triangle polygons
    #'   populated by `prep()`. Each entry is `NULL` or
    #'   `list(list(x, y), list(x, y))`.
    centroPolys = NULL,
    #' @field highlightBoxes Per-window list of highlight rectangles
    #'   `list(x0, x1, y0, y1)`. Populated only when `scope = "full"`.
    highlightBoxes = NULL,
    #' @field scope Character. `"window"` (default) or `"full"`.
    scope = "window",
    #' @field style Character. `"block"` (default) or `"rounded"`.
    style = "block",
    #' @field highlight_range Optional `GRanges`. When set (and
    #'   `scope = "full"`), the highlight rectangle covers this range
    #'   instead of the parent track's windows.
    highlight_range = NULL,

    #' @description Construct a SeqIdeogramR6.
    #' @param data Optional `GRanges` of cytogenetic bands with a `gieStain`
    #'   mcol (or a `stain` mapping in `mapping`). Falls back to the parent
    #'   track's data.
    #' @param mapping Optional `SeqMap`. Recognised: `stain` (defaults to
    #'   the `gieStain` mcol when absent).
    #' @param aesthetics Optional `SeqAes`. Recognised: `color` (band
    #'   border, default `"black"`), `linewidth` (band border width,
    #'   default `0.1`), `outline` (nested aes for the chromosome's
    #'   perimeter outline; sub-keys: `col`, `lwd`, `visible`),
    #'   `highlight` (nested aes for the `scope = "full"` highlight box;
    #'   sub-keys: `fill`, `col`, `lwd`, `alpha`),
    #'   `telomere.radius` (corner radius as a fraction of band height for
    #'   `style = "rounded"`; default `0.4`).
    #' @param scope Character. `"window"` (default — only bands overlapping
    #'   the track windows are drawn) or `"full"` (the whole chromosome is
    #'   drawn rescaled to the panel; the current window is overlaid as a
    #'   highlight rectangle).
    #' @param style Character. `"block"` (default — rectangular bands) or
    #'   `"rounded"` (rounded caps on the leftmost and rightmost bands).
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          aesthetics = aes(),
                          scope = "window",
                          style = "block",
                          highlight_range = NULL,
                          ...) {
      super$initialize(data, mapping, aesthetics, ...)
      self$scope          <- match.arg(scope, c("window", "full"))
      self$style          <- match.arg(style, c("block", "rounded"))
      if (!is.null(highlight_range) &&
          !inherits(highlight_range, "GRanges"))
        stop("`highlight_range` must be a GRanges or NULL.", call. = FALSE)
      self$highlight_range <- highlight_range
      self$centroPolys     <- NULL
      self$highlightBoxes  <- NULL
    },

    #' @description Resolve mapping, find band overlaps with each window,
    #'   and populate `coordCanvas` (non-centromeric bands) and
    #'   `centroPolys` (per-window centromere triangles). When
    #'   `scope = "full"`, the whole chromosome is mapped into the panel
    #'   and `highlightBoxes` is populated with one entry per window.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      n_panels <- length(track_windows)
      self$coordCanvas    <- vector("list", n_panels)
      self$centroPolys    <- vector("list", n_panels)
      self$highlightBoxes <- vector("list", n_panels)

      if (is.null(eff_data) || length(eff_data) == 0) return(invisible())
      if (!inherits(eff_data, "GRanges"))
        stop("seq_ideogram requires GRanges `data`.", call. = FALSE)

      stain_all <- self$resolved$stain %||%
        S4Vectors::mcols(eff_data)$gieStain
      if (is.null(stain_all))
        stop("seq_ideogram needs a `gieStain` mcol or a `stain` mapping.",
             call. = FALSE)
      stain_all <- as.character(stain_all)

      # ── scope = "full" path: whole-chromosome view ───────────────────────
      if (self$scope == "full") {
        chr <- as.character(GenomicRanges::seqnames(track_windows[1]))
        keep <- as.character(GenomicRanges::seqnames(eff_data)) == chr
        all_bands <- eff_data[keep]
        stain_chr <- stain_all[keep]
        if (length(all_bands) == 0L) return(invisible())

        chr_start <- min(BiocGenerics::start(all_bands))
        chr_end   <- max(BiocGenerics::end(all_bands))
        chr_span  <- chr_end - chr_start
        if (!is.finite(chr_span) || chr_span <= 0) return(invisible())

        for (w in seq_len(n_panels)) {
          pm    <- layout_track[[w]]
          win_w <- track_windows[w]

          u0 <- (BiocGenerics::start(all_bands) - chr_start) / chr_span
          u1 <- (BiocGenerics::end(all_bands)   - chr_start) / chr_span
          x0c <- pm$inner$x0 + u0 * (pm$inner$x1 - pm$inner$x0)
          x1c <- pm$inner$x0 + u1 * (pm$inner$x1 - pm$inner$x0)
          y0c <- pm$inner$y0; y1c <- pm$inner$y1

          fill_cols  <- .ideogram_fill_colors(stain_chr)
          non_acen   <- stain_chr != "acen"

          self$coordCanvas[[w]] <- data.frame(
            x0 = x0c[non_acen], x1 = x1c[non_acen],
            y0 = y0c, y1 = y1c,
            fill = fill_cols[non_acen],
            is_tel_left  = seq_along(x0c)[non_acen] == 1L,
            is_tel_right = seq_along(x0c)[non_acen] == length(non_acen),
            stringsAsFactors = FALSE
          )

          # Centromere triangles
          acen_idx <- which(stain_chr == "acen")
          if (length(acen_idx) == 2L) {
            x0a <- x0c[acen_idx[1]]; x1a <- x1c[acen_idx[1]]
            x0b <- x0c[acen_idx[2]]; x1b <- x1c[acen_idx[2]]
            ym  <- (y0c + y1c) / 2
            self$centroPolys[[w]] <- list(
              list(x = c(x0a, x1a, x0a), y = c(y0c, ym, y1c)),
              list(x = c(x1b, x0b, x1b), y = c(y0c, ym, y1c))
            )
          }

          # Highlight box: prefer explicit highlight_range when given;
          # otherwise fall back to the current window region.
          hl_w <- if (!is.null(self$highlight_range))
                    self$highlight_range else win_w
          hl_chr <- as.character(GenomicRanges::seqnames(hl_w))
          if (any(hl_chr == chr) &&
              BiocGenerics::start(hl_w[hl_chr == chr][1]) <= chr_end &&
              BiocGenerics::end  (hl_w[hl_chr == chr][1]) >= chr_start) {
            hl_w  <- hl_w[hl_chr == chr][1]
            win_u0 <- (BiocGenerics::start(hl_w) - chr_start) / chr_span
            win_u1 <- (BiocGenerics::end(hl_w)   - chr_start) / chr_span
            win_u0 <- pmax(pmin(win_u0, 1), 0)
            win_u1 <- pmax(pmin(win_u1, 1), 0)
            win_x0 <- pm$inner$x0 + win_u0 * (pm$inner$x1 - pm$inner$x0)
            win_x1 <- pm$inner$x0 + win_u1 * (pm$inner$x1 - pm$inner$x0)
            # Suppress when highlight covers the entire chromosome.
            if (win_x1 - win_x0 < (pm$inner$x1 - pm$inner$x0) - 1e-6) {
              self$highlightBoxes[[w]] <- list(x0 = win_x0, x1 = win_x1,
                                               y0 = y0c,    y1 = y1c)
            }
          }
        }
        return(invisible())
      }

      # ── scope = "window" path: existing overlap-only logic ───────────────
      ov <- suppressWarnings(
        GenomicRanges::findOverlaps(eff_data, track_windows)
      )
      if (length(ov) == 0L) return(invisible())

      qh <- S4Vectors::queryHits(ov)
      sh <- S4Vectors::subjectHits(ov)

      for (w in unique(sh)) {
        pm    <- layout_track[[w]]
        idx   <- qh[sh == w]
        bands <- eff_data[idx]
        stain <- stain_all[idx]

        u0 <- (BiocGenerics::start(bands) - pm$xscale[1]) / diff(pm$xscale)
        u1 <- (BiocGenerics::end(bands)   - pm$xscale[1]) / diff(pm$xscale)
        u0 <- pmax(pmin(u0, 1), 0); u1 <- pmax(pmin(u1, 1), 0)

        x0c <- pm$inner$x0 + u0 * (pm$inner$x1 - pm$inner$x0)
        x1c <- pm$inner$x0 + u1 * (pm$inner$x1 - pm$inner$x0)
        y0c <- pm$inner$y0
        y1c <- pm$inner$y1

        fill_cols <- .ideogram_fill_colors(stain)

        non_acen <- stain != "acen"
        self$coordCanvas[[w]] <- data.frame(
          x0 = x0c[non_acen], x1 = x1c[non_acen],
          y0 = y0c,           y1 = y1c,
          fill = fill_cols[non_acen],
          is_tel_left  = seq_along(x0c)[non_acen] == 1L,
          is_tel_right = seq_along(x0c)[non_acen] == length(non_acen),
          stringsAsFactors = FALSE
        )

        acen_idx <- which(stain == "acen")
        if (length(acen_idx) == 2L) {
          x0a <- x0c[acen_idx[1]]; x1a <- x1c[acen_idx[1]]
          x0b <- x0c[acen_idx[2]]; x1b <- x1c[acen_idx[2]]
          ym  <- (y0c + y1c) / 2
          self$centroPolys[[w]] <- list(
            list(x = c(x0a, x1a, x0a), y = c(y0c, ym, y1c)),
            list(x = c(x1b, x0b, x1b), y = c(y0c, ym, y1c))
          )
        }
      }
      invisible()
    },

    #' @description Draw non-centromeric bands with `grid::grid.rect()` and
    #'   centromere triangles with `grid::grid.polygon()`. With
    #'   `style = "rounded"`, the leftmost / rightmost bands draw with
    #'   rounded telomere caps. With `scope = "full"`, the current window
    #'   region is overlaid as a highlight rectangle.
    draw = function() {
      border_col <- self$aesthetics$color     %||% "black"
      border_lwd <- self$aesthetics$linewidth %||% 0.1
      radius_frac <- self$aesthetics[["telomere.radius"]] %||% 0.4

      n_arc <- 32L

      if (!is.null(self$coordCanvas)) {
        for (coords in self$coordCanvas) {
          if (is.null(coords) || nrow(coords) == 0L) next

          band_h <- coords$y1[1] - coords$y0[1]
          # Compute cap horizontal radius in physical units so the cap
          # is aspect-ratio-correct. `telomere.radius = 1.0` gives a
          # full semicircle whose depth equals band_height / 2 inches,
          # independent of canvas pixel aspect ratio.
          band_h_in <- tryCatch(
            grid::convertHeight(grid::unit(band_h, "npc"),
                                 "inches", valueOnly = TRUE),
            error = function(e) band_h
          )
          rx_inches <- max(0, radius_frac) * band_h_in / 2
          rx_max <- tryCatch(
            grid::convertWidth(grid::unit(rx_inches, "inches"),
                                "npc", valueOnly = TRUE),
            error = function(e) rx_inches
          )

          do_round <- identical(self$style, "rounded") && rx_max > 0

          for (j in seq_len(nrow(coords))) {
            bx0 <- coords$x0[j]; bx1 <- coords$x1[j]
            by0 <- coords$y0[j]; by1 <- coords$y1[j]
            bfill <- coords$fill[j]
            cy <- (by0 + by1) / 2
            ry <- (by1 - by0) / 2

            is_left  <- do_round && isTRUE(coords$is_tel_left[j])
            is_right <- do_round && isTRUE(coords$is_tel_right[j])

            # Cap depth, clamped so the cap never extends past the
            # band's own right (or left) edge.
            cap_w <- min(rx_max, bx1 - bx0)

            if (!is_left && !is_right) {
              grid::grid.rect(
                x      = grid::unit((bx0 + bx1) / 2, "npc"),
                y      = grid::unit(cy, "npc"),
                width  = grid::unit(bx1 - bx0, "npc"),
                height = grid::unit(by1 - by0, "npc"),
                gp = grid::gpar(fill = bfill,
                                col  = border_col,
                                lwd  = border_lwd)
              )
              next
            }

            # Trace the band's visible region — the intersection of
            # the band rectangle with the local stadium-shaped cap.
            cx_l <- bx0 + cap_w
            cx_r <- bx1 - cap_w
            ys <- seq(by0, by1, length.out = n_arc)
            arc_l <- if (is_left) {
              cx_l - cap_w * sqrt(pmax(0, 1 - ((ys - cy) / ry)^2))
            } else rep(-Inf, n_arc)
            arc_r <- if (is_right) {
              cx_r + cap_w * sqrt(pmax(0, 1 - ((ys - cy) / ry)^2))
            } else rep( Inf, n_arc)

            left_edge  <- pmax(bx0, arc_l)
            right_edge <- pmin(bx1, arc_r)
            vis <- left_edge <= right_edge
            if (!any(vis)) next

            xs <- c(left_edge[vis], rev(right_edge[vis]))
            ys <- c(ys[vis],        rev(ys[vis]))

            grid::grid.polygon(
              x  = grid::unit(xs, "npc"),
              y  = grid::unit(ys, "npc"),
              gp = grid::gpar(fill = bfill,
                              col  = border_col,
                              lwd  = border_lwd)
            )
          }
        }
      }

      if (!is.null(self$centroPolys)) {
        for (polys in self$centroPolys) {
          if (is.null(polys)) next
          for (tri in polys) {
            grid::grid.polygon(
              x  = grid::unit(tri$x, "npc"),
              y  = grid::unit(tri$y, "npc"),
              gp = grid::gpar(fill = "#FF0000",
                              col  = border_col,
                              lwd  = border_lwd)
            )
          }
        }
      }

      # ── Outer perimeter outline ─────────────────────────────────────
      # Trace a single closed polygon around the outside of the
      # chromosome (rectangular for `block`, stadium-shaped for
      # `rounded`) and stroke it on top of the bands.
      ol_aes     <- self$aesthetics[["outline"]] %||% list()
      ol_visible <- if (is.null(ol_aes$visible)) TRUE
                    else isTRUE(ol_aes$visible)

      if (ol_visible && !is.null(self$coordCanvas)) {
        ol_col <- ol_aes$col %||% border_col
        ol_lwd <- ol_aes$lwd %||% max(border_lwd, 1)

        for (coords in self$coordCanvas) {
          if (is.null(coords) || nrow(coords) == 0L) next

          left_idx  <- which(coords$is_tel_left  %in% TRUE)
          right_idx <- which(coords$is_tel_right %in% TRUE)
          if (length(left_idx) == 0L || length(right_idx) == 0L) next

          by0 <- coords$y0[1]; by1 <- coords$y1[1]
          cy  <- (by0 + by1) / 2
          ry  <- (by1 - by0) / 2

          chr_left  <- coords$x0[left_idx[1]]
          chr_right <- coords$x1[right_idx[1]]

          if (identical(self$style, "rounded")) {
            band_h_in <- tryCatch(
              grid::convertHeight(grid::unit(by1 - by0, "npc"),
                                   "inches", valueOnly = TRUE),
              error = function(e) by1 - by0
            )
            rx_inches <- max(0, radius_frac) * band_h_in / 2
            rx_max <- tryCatch(
              grid::convertWidth(grid::unit(rx_inches, "inches"),
                                  "npc", valueOnly = TRUE),
              error = function(e) rx_inches
            )
            cap_w_l <- min(rx_max,
                           coords$x1[left_idx[1]]  - coords$x0[left_idx[1]])
            cap_w_r <- min(rx_max,
                           coords$x1[right_idx[1]] - coords$x0[right_idx[1]])
          } else {
            cap_w_l <- 0
            cap_w_r <- 0
          }

          if (cap_w_l > 0 || cap_w_r > 0) {
            n_arc <- 32L
            # Left half-ellipse: top → leftmost → bottom (CCW)
            t_l <- seq(pi / 2, 3 * pi / 2, length.out = n_arc)
            ax_l <- (chr_left + cap_w_l) + cap_w_l * cos(t_l)
            ay_l <- cy + ry * sin(t_l)
            # Right half-ellipse: bottom → rightmost → top (CCW)
            t_r <- seq(-pi / 2, pi / 2, length.out = n_arc)
            ax_r <- (chr_right - cap_w_r) + cap_w_r * cos(t_r)
            ay_r <- cy + ry * sin(t_r)
            xs <- c(ax_l, ax_r)
            ys <- c(ay_l, ay_r)
          } else {
            xs <- c(chr_left, chr_right, chr_right, chr_left)
            ys <- c(by1,      by1,       by0,       by0)
          }

          grid::grid.polygon(
            x  = grid::unit(xs, "npc"),
            y  = grid::unit(ys, "npc"),
            gp = grid::gpar(fill = NA, col = ol_col, lwd = ol_lwd)
          )
        }
      }

      # Highlight rectangles (only populated when scope = "full").
      if (!is.null(self$highlightBoxes)) {
        hl_aes   <- self$aesthetics[["highlight"]] %||% list()
        hl_fill  <- hl_aes$fill  %||% "red"
        hl_col   <- hl_aes$col   %||% NA
        hl_lwd   <- hl_aes$lwd   %||% 0.8
        hl_alpha <- hl_aes$alpha %||% 0.15

        for (hb in self$highlightBoxes) {
          if (is.null(hb)) next
          if (!is.finite(hb$x1 - hb$x0) || hb$x1 <= hb$x0) next
          grid::grid.rect(
            x      = grid::unit((hb$x0 + hb$x1) / 2, "npc"),
            y      = grid::unit((hb$y0 + hb$y1) / 2, "npc"),
            width  = grid::unit(hb$x1 - hb$x0, "npc"),
            height = grid::unit(hb$y1 - hb$y0, "npc"),
            gp = grid::gpar(fill  = hl_fill,
                            col   = hl_col,
                            lwd   = hl_lwd,
                            alpha = hl_alpha)
          )
        }
      }
      invisible()
    }
  )
)

#' Chromosome ideogram from cytogenetic bands
#'
#' Draws a chromosome ideogram: each cytoband becomes a filled rectangle
#' shaded by its Giemsa stain (`gpos25` darkens through `gpos100`; `gneg`
#' is white; `stalk` and `gvar` carry their conventional colors). Paired
#' `acen` bands within a window render as two inward-pointing red
#' triangles, marking the centromere.
#'
#' The simplest call supplies a `GRanges` of cytobands — use
#' [load_cytobands()] to load the bundled hg38 table:
#' ```r
#' cb <- load_cytobands()
#' seq_plot() %|%
#'   seq_track(track_id = "Ideo",
#'             windows = default_genome_windows()) %+%
#'   seq_ideogram(data = cb)
#' ```
#'
#' Map a non-standard stain column with `map(stain = my_col)`.
#'
#' @param data Optional `GRanges` of cytobands. Falls back to the parent
#'   track's data. Must carry a `gieStain` mcol unless a `stain` mapping
#'   is supplied.
#' @param mapping Optional [map()]. Recognised: `stain`.
#' @param aesthetics Optional [aes()]: `color` (band border color,
#'   default `"black"`), `linewidth` (band border width, default
#'   `0.1`), `outline` (nested `aes()` controlling the chromosome's
#'   outer perimeter outline: sub-keys `col`, `lwd`, `visible`),
#'   `highlight` (nested `aes()` controlling the `scope = "full"`
#'   highlight rectangle: sub-keys `fill`, `col`, `lwd`, `alpha`),
#'   `telomere.radius` (numeric; corner radius as a fraction of band
#'   height for `style = "rounded"`; `1.0` = full half-circle cap).
#' @param scope Character. One of `"window"` (default — only bands
#'   overlapping the track windows are drawn) or `"full"` (the whole
#'   chromosome is drawn rescaled to fill the panel; the current window
#'   region is overlaid as a translucent highlight rectangle).
#' @param style Character. One of `"block"` (default — rectangular
#'   bands) or `"rounded"` (rounded telomere caps on the leftmost and
#'   rightmost bands).
#' @param highlight_range Optional `GRanges`. Only honoured when
#'   `scope = "full"`. When set, the highlight rectangle marks this
#'   range instead of the parent track's `windows`. Lets you set the
#'   track's `windows` to span the full chromosome (so the x-axis
#'   reads chromosome coordinates) while still highlighting a
#'   sub-range.
#' @param ... Reserved.
#' @return A `SeqIdeogramR6` instance.
#' @export
seq_ideogram <- function(data = NULL, mapping = NULL,
                         aesthetics = aes(),
                         scope = "window",
                         style = "block",
                         highlight_range = NULL,
                         ...) {
  SeqIdeogramR6$new(data = data, mapping = mapping,
                    aesthetics = aesthetics,
                    scope = scope, style = style,
                    highlight_range = highlight_range, ...)
}
