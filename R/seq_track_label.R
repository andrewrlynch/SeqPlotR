# SeqTrackLabelR6 class ----

#' SeqTrackLabel R6 class
#'
#' Internal R6 generator backing [seq_track_label()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqTrackLabelR6 <- R6::R6Class("SeqTrackLabel",
  inherit = SeqElementR6,
  public = list(
    #' @field label Row label text.
    label = NULL,
    #' @field sublabel Optional second-line text.
    sublabel = NULL,
    #' @field side `"left"` or `"right"`.
    side = "left",

    #' @description Construct a SeqTrackLabelR6.
    #' @param label Label text.
    #' @param sublabel Optional second line.
    #' @param side `"left"` or `"right"`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param ... Unused.
    initialize = function(label, sublabel = NULL, side = c("left", "right"),
                          aesthetics = aes(), ...) {
      super$initialize(data = NULL, mapping = NULL, aesthetics = aesthetics,
                       show_legend = FALSE)
      self$label    <- label
      self$sublabel <- sublabel
      self$side     <- match.arg(side)
    },

    #' @description Place the label beside the track's window area.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges (unused).
    #' @return The `coordCanvas` list, invisibly.
    prep = function(layout_track, track_windows) {
      self$coordCanvas <- NULL
      if (length(layout_track) == 0L) return(invisible(NULL))
      lab <- as.character(self$label)[1]
      sub <- if (is.null(self$sublabel)) NA_character_ else as.character(self$sublabel)[1]
      if (is.na(lab) || !nzchar(lab)) return(invisible(NULL))
      p1 <- layout_track[[1]]
      x  <- if (self$side == "left") p1$track_inner$x0
            else layout_track[[length(layout_track)]]$track_inner$x1
      self$coordCanvas <- list(
        x = x, y = (p1$inner$y0 + p1$inner$y1) / 2, label = lab,
        sublabel = if (is.na(sub) || !nzchar(sub)) NULL else sub, side = self$side
      )
      invisible(self$coordCanvas)
    },

    #' @description Draw the label and optional sublabel.
    draw = function() {
      cc <- self$coordCanvas
      if (is.null(cc)) return(invisible())
      a   <- self$aesthetics
      off <- a$offset %||% 0.6
      sgn <- if (cc$side == "left") -1 else 1
      x   <- grid::unit(cc$x, "npc") + sgn * grid::unit(off, "lines")
      just <- if (cc$side == "left") "right" else "left"
      has_sub <- !is.null(cc$sublabel)
      grid::grid.text(
        cc$label, x = x,
        y = grid::unit(cc$y, "npc") + grid::unit(if (has_sub) 0.45 else 0, "lines"),
        just = just, name = "seqtracklabel.label",
        gp = grid::gpar(cex = a$size %||% 0.75, col = a$col %||% "#1C1B1A")
      )
      if (has_sub)
        grid::grid.text(
          cc$sublabel, x = x, y = grid::unit(cc$y, "npc") - grid::unit(0.5, "lines"),
          just = just, name = "seqtracklabel.sublabel",
          gp = grid::gpar(cex = a$sublabel.size %||% 0.65,
                          col = a$sublabel.col %||% "grey35",
                          fontface = a$sublabel.face %||% "italic")
        )
      invisible()
    }
  )
)

#' Track label with an optional second line
#'
#' Draws a text label beside a track's window area, vertically centred on the
#' track, with an optional smaller second line beneath it. Useful for naming
#' rows in multi-track figures (samples, haplotypes, assembled sequences)
#' without relying on axis titles.
#'
#' @param label Label text (a single string). `NA` or `""` draws nothing.
#' @param sublabel Optional second line (a single string). `NULL` or `NA`
#'   draws no second line.
#' @param side `"left"` (default) places the label right-aligned just left of
#'   the window area; `"right"` places it left-aligned just right of it.
#' @param aesthetics Optional `SeqAes`. Keys: `size` (cex, default 0.75),
#'   `col` (default `"#1C1B1A"`), `offset` (gap from the window area in lines,
#'   default 0.6), `sublabel.size` (default 0.65), `sublabel.col` (default
#'   `"grey35"`), `sublabel.face` (default `"italic"`).
#' @return A `SeqTrackLabelR6` instance.
#' @examples
#' seq_track_label("sample 1", sublabel = "haplotype A")
#' @export
seq_track_label <- function(label, sublabel = NULL, side = c("left", "right"),
                            aesthetics = aes()) {
  SeqTrackLabelR6$new(label = label, sublabel = sublabel, side = side,
                      aesthetics = aesthetics)
}
