# PLOT-LEVEL ANNOTATIONS

#' SeqAnnotation R6 base class
#'
#' Internal R6 generator for plot-level annotations added to a `SeqPlot`
#' with `%+%`. Subclasses implement `prep(layout)` and `draw()`; `allowance`
#' is the bottom space (in lines) reserved for them by `layoutGrid()`.
#'
#' @keywords internal
SeqAnnotationR6 <- R6::R6Class("SeqAnnotation",
  public = list(
    #' @field allowance Bottom space reserved, in lines.
    allowance   = 0,
    #' @field coordCanvas Canvas positions set by `prep()`.
    coordCanvas = NULL,

    #' @description Override in subclasses. Default implementation errors.
    #' @param layout The plot's `layout` list.
    prep = function(layout) {
      stop("prep() must be implemented by ", class(self)[1], call. = FALSE)
    },

    #' @description Override in subclasses. Default implementation errors.
    draw = function() {
      stop("draw() must be implemented by ", class(self)[1], call. = FALSE)
    }
  )
)

#' SeqCaption R6 class
#'
#' Internal R6 generator backing [seq_caption()].
#'
#' @keywords internal
SeqCaptionR6 <- R6::R6Class("SeqCaption",
  inherit = SeqAnnotationR6,
  public = list(
    #' @field text Caption text.
    text       = NULL,
    #' @field banner Optional banner text.
    banner     = NULL,
    #' @field aesthetics `SeqAes` with caption and banner styling.
    aesthetics = NULL,

    #' @description Construct a SeqCaptionR6.
    #' @param text Caption text.
    #' @param banner Optional banner text.
    #' @param allowance Bottom space reserved, in lines.
    #' @param aesthetics Optional `SeqAes`.
    initialize = function(text, banner = NULL, allowance = 2.8,
                          aesthetics = aes()) {
      .one <- function(v, nm) {
        if (is.null(v)) return(NULL)
        if (length(v) != 1L)
          stop("seq_caption ", nm, " must be a single string or NULL.",
               call. = FALSE)
        as.character(v)
      }
      if (!is.numeric(allowance) || length(allowance) != 1L ||
          !is.finite(allowance) || allowance < 0)
        stop("seq_caption allowance must be a single non-negative number.",
             call. = FALSE)
      self$text       <- .one(text, "text")
      self$banner     <- .one(banner, "banner")
      self$allowance  <- allowance
      self$aesthetics <- if (is.null(aesthetics)) aes() else aesthetics
    },

    #' @description Place the caption at the centre top of the bottom canvas
    #'   margin.
    #' @param layout The plot's `layout` list (uses `canvasMarginBounds`).
    prep = function(layout) {
      b <- layout$canvasMarginBounds$bottom
      self$coordCanvas <- list(x = (b$x0 + b$x1) / 2, y = b$y1,
                               text = self$text, banner = self$banner)
      invisible()
    },

    #' @description Draw the caption line and the banner under it.
    draw = function() {
      cc <- self$coordCanvas
      if (is.null(cc)) return(invisible())
      a <- self$aesthetics
      .ok <- function(v) !is.null(v) && !is.na(v) && nzchar(v)
      if (.ok(cc$text))
        grid::grid.text(cc$text, x = grid::unit(cc$x, "npc"),
                        y = grid::unit(cc$y, "npc") - grid::unit(0.9, "lines"),
                        gp = grid::gpar(cex = a$caption.size %||% 0.6,
                                        col = a$caption.col %||% "grey30"),
                        name = "seqcaption.text")
      if (.ok(cc$banner))
        grid::grid.text(cc$banner, x = grid::unit(cc$x, "npc"),
                        y = grid::unit(cc$y, "npc") - grid::unit(1.9, "lines"),
                        gp = grid::gpar(cex = a$banner.size %||% 0.7,
                                        col = a$banner.col %||% "#D33A2C",
                                        fontface = a$banner.face %||% "bold"),
                        name = "seqcaption.banner")
      invisible()
    }
  )
)

#' Caption and banner under the figure
#'
#' Adds a caption line centred under the whole figure and an optional
#' banner line beneath it (for example a warning that the figure is
#' approximate). Both strings are drawn exactly as given. The bottom canvas
#' margin grows by a fixed `allowance` (in lines) whether or not a banner is
#' given, so track heights never depend on the caption content.
#'
#' Add to a plot with `%+%`; a caption cannot be added to a track.
#'
#' @param text Caption text (a single string). `NULL`, `NA` or `""` draws
#'   nothing.
#' @param banner Optional banner text, drawn one line under the caption.
#' @param allowance Bottom space reserved, in lines (default `2.8`).
#' @param aesthetics Optional [aes()]: `caption.size` (cex, `0.6`),
#'   `caption.col` (`"grey30"`), `banner.size` (`0.7`), `banner.col`
#'   (`"#D33A2C"`), `banner.face` (`"bold"`).
#' @return A `SeqCaptionR6` instance.
#' @examples
#' w <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1e6))
#' p <- seq_plot() %+% seq_track(windows = w) %+%
#'   seq_caption("Model output, not truth", banner = "APPROXIMATE")
#' @export
seq_caption <- function(text, banner = NULL, allowance = 2.8,
                        aesthetics = aes()) {
  SeqCaptionR6$new(text = text, banner = banner, allowance = allowance,
                   aesthetics = aesthetics)
}
