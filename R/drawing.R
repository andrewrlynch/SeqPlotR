# ── Internal drawing helpers for link elements ───────────────────────────────
#
# Ported from THEfunc `drawing.R`. `scales::alpha()` (a non-dependency) is
# replaced with `grDevices::adjustcolor(..., alpha.f=)` everywhere it appears.

#' Draw a single cubic-Bezier arch and its vertical stems
#'
#' Internal helper used by `seq_arc`, `seq_arch`, and `seq_recon`. Draws the
#' Bezier arch connecting two stem tops at canvas npc coordinates, then the
#' two vertical stems. Pass `stemWidth = 0` to suppress stem rendering.
#'
#' @param x0,y0 Canvas npc coordinates of the first stem base.
#' @param x1,y1 Canvas npc coordinates of the second stem base.
#' @param top0,top1 Canvas npc y-coordinates where the arch meets the stem
#'   tops at `x0` and `x1`.
#' @param orientation Character. Controls which side the arch bulges to.
#'   `"+"` / `"*"` curves up; `"-"` curves down. Mixed strings such as
#'   `"+/-"` control each end independently.
#' @param curve `"length"` (offset scales with span), `"equal"` (fixed
#'   `0.2` offset), or a numeric value used directly as the offset.
#' @param stemWidth,stemColor Stem aesthetic.
#' @param arcWidth,arcColor Arch aesthetic.
#' @return Invisible `NULL`.
#' @keywords internal
drawSeqArch <- function(x0, y0, x1, y1, top0, top1,
                        orientation = "*", curve = "length",
                        stemWidth = 1, arcWidth = 1,
                        arcColor = "black", stemColor = "black") {

  span <- abs(x1 - x0)

  if (is.numeric(curve)) {
    curve_offset <- curve
  } else if (identical(curve, "equal")) {
    curve_offset <- 0.2
  } else if (identical(curve, "length")) {
    curve_offset <- span * 0.2
  } else {
    warning("Unknown curve value in drawSeqArch(); using default 0.2")
    curve_offset <- span * 0.2
  }

  curve_offset1 <- ifelse(grepl("^\\*|^\\+", orientation),
                          curve_offset, curve_offset * -1)
  curve_offset2 <- ifelse(grepl("\\*$|\\+$", orientation),
                          curve_offset, curve_offset * -1)

  ctrl_spread <- 0
  dx          <- abs(x1 - x0)
  ctrl_dx     <- dx * ctrl_spread

  mid_y   <- max(top0, top1)
  ctrl_x1 <- x0 + ctrl_dx
  ctrl_x2 <- x1 - ctrl_dx
  ctrl_y1 <- mid_y + curve_offset1
  ctrl_y2 <- mid_y + curve_offset2

  P0 <- c(x0,      top0)
  P1 <- c(ctrl_x1, ctrl_y1)
  P2 <- c(ctrl_x2, ctrl_y2)
  P3 <- c(x1,      top1)

  t <- seq(0, 1, length.out = 100)
  bez_x <- (1 - t)^3 * P0[1] + 3 * (1 - t)^2 * t * P1[1] +
    3 * (1 - t) * t^2 * P2[1] + t^3 * P3[1]
  bez_y <- (1 - t)^3 * P0[2] + 3 * (1 - t)^2 * t * P1[2] +
    3 * (1 - t) * t^2 * P2[2] + t^3 * P3[2]

  grid::grid.lines(bez_x, bez_y,
                   gp = grid::gpar(col = arcColor, lwd = arcWidth))

  if (isTRUE(stemWidth > 0)) {
    grid::grid.segments(x0, y0, x0, top0,
                        gp = grid::gpar(col = stemColor, lwd = stemWidth))
    grid::grid.segments(x1, y1, x1, top1,
                        gp = grid::gpar(col = stemColor, lwd = stemWidth))
  }

  invisible(NULL)
}

#' Draw a single cubic-Bezier string between two anchors
#'
#' Internal helper used by `seq_string`. Draws a smooth horizontal Bezier
#' curve connecting `(x0, y0)` and `(x1, y1)` in canvas npc coordinates. The
#' shape (C vs. S) is determined by the strand pair — matching strands bend
#' both control points the same direction, opposing strands bend them in
#' opposite directions.
#'
#' @param x0,y0 Canvas npc coordinates of the first anchor.
#' @param x1,y1 Canvas npc coordinates of the second anchor.
#' @param strand1,strand2 Strand character for each anchor (`"+"`, `"-"`, `"*"`).
#' @param orientation Optional explicit direction (`"+"` / `"-"`); when
#'   unset, direction is inferred from `strand1`.
#' @param type `"c"` or `"s"`. Not consulted directly — retained for API
#'   compatibility; the C vs. S shape is driven by the strand pair.
#' @param bulge Horizontal control-point offset in npc units. Clamped to
#'   `[0, 0.35]`.
#' @param lwd,col,alpha Line aesthetics.
#' @return Invisible `NULL`.
#' @keywords internal
drawSeqString <- function(x0, y0, x1, y1,
                          strand1 = "*",
                          strand2 = "*",
                          orientation = "*",
                          type = c("c", "s"),
                          bulge = 0.04,
                          lwd = 1.5,
                          col = "red",
                          alpha = 1) {

  if (length(x0) == 0 || length(y0) == 0 ||
      length(x1) == 0 || length(y1) == 0)
    return(invisible(NULL))

  x0 <- as.numeric(x0)[1]; y0 <- as.numeric(y0)[1]
  x1 <- as.numeric(x1)[1]; y1 <- as.numeric(y1)[1]
  if (!is.finite(x0) || !is.finite(y0) ||
      !is.finite(x1) || !is.finite(y1))
    return(invisible(NULL))

  strand1     <- as.character(strand1[1])
  strand2     <- as.character(strand2[1])
  orientation <- tolower(as.character(orientation[1]))

  type <- tolower(as.character(type[1]))
  if (!type %in% c("c", "s")) type <- "c"

  # keep stable left->right; if we swap endpoints, swap strand roles too
  if (x1 < x0) {
    tx <- x0; x0 <- x1; x1 <- tx
    ty <- y0; y0 <- y1; y1 <- ty
    ts <- strand1; strand1 <- strand2; strand2 <- ts
  }

  dx <- x1 - x0
  if (!is.finite(dx) || dx <= 0) return(invisible(NULL))

  bx <- as.numeric(bulge)[1]
  if (!is.finite(bx)) bx <- 0.07
  bx <- max(0, min(0.35, bx))

  dir <- 1
  if (orientation %in% c("-", "down", "left", "neg")) {
    dir <- -1
  } else if (orientation %in% c("+", "up", "right", "pos")) {
    dir <- 1
  } else {
    if (strand1 == "-") dir <- -1
    if (strand1 == "+") dir <- 1
  }

  yP1 <- y0
  yP2 <- y1

  same_strand <- (strand1 %in% c("+", "-")) &&
                 (strand2 %in% c("+", "-")) &&
                 (strand1 == strand2)

  if (same_strand) {
    xP1 <- x0 + dir * bx
    xP2 <- x1 + dir * bx
  } else {
    xP1 <- x0 + dir * bx
    xP2 <- x1 - dir * bx
  }

  xP1 <- max(0, min(1, xP1))
  xP2 <- max(0, min(1, xP2))

  gp <- grid::gpar(
    col = grDevices::adjustcolor(col, alpha.f = alpha),
    lwd = lwd
  )

  grid::grid.draw(grid::bezierGrob(
    x = c(x0, xP1, xP2, x1),
    y = c(y0, yP1, yP2, y1),
    default.units = "npc",
    gp = gp
  ))

  invisible(NULL)
}

#' Infer C vs. S curve type from a pair of strands
#'
#' Opposing strands (`+/-` or `-/+`) produce an S curve; matching or
#' unknown strands fall back to the supplied default (typically `"c"`).
#'
#' @param strand0,strand1 Character vectors of strands (`"+"`, `"-"`, `"*"`).
#' @param default Fallback curve type when strand info is missing or matching.
#' @return Character vector of `"c"` or `"s"` values.
#' @keywords internal
.string_type_from_strand <- function(strand0, strand1, default = "c") {
  code <- paste0(as.character(strand0), "/", as.character(strand1))
  ifelse(code %in% c("+/-", "-/+"), "s", default)
}

#' Convert genomic x and data y to canvas npc coordinates for a single panel
#'
#' Internal helper used by link prep methods. Clamps to the panel's `inner`
#' rectangle so any value outside the data ranges lands on the boundary.
#'
#' @param x_gen Numeric vector of genomic x positions.
#' @param y_data Numeric vector of data-scale y values (same length as `x_gen`).
#' @param panel_meta Panel metadata list (with `xscale`, `yscale`, and
#'   `inner` containing `x0`, `x1`, `y0`, `y1`).
#' @return List with `x` and `y` numeric vectors of canvas npc coordinates.
#' @keywords internal
.to_canvas <- function(x_gen, y_data, panel_meta) {
  u <- pmax(pmin((x_gen  - panel_meta$xscale[1]) /
                   diff(panel_meta$xscale), 1), 0)
  v <- pmax(pmin((y_data - panel_meta$yscale[1]) /
                   diff(panel_meta$yscale), 1), 0)
  list(
    x = panel_meta$inner$x0 + u * (panel_meta$inner$x1 - panel_meta$inner$x0),
    y = panel_meta$inner$y0 + v * (panel_meta$inner$y1 - panel_meta$inner$y0)
  )
}
