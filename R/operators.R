# %+% — primary operator ------------------------------------------------------

#' Add a track, element, or plot-level feature to a SeqPlot or SeqTrack
#'
#' The single composition operator for SeqPlotR. Dispatches on the class of
#' the right-hand side:
#'
#' - `SeqTrack`      -> added to the plot layout
#' - `SeqElement`    -> added to the current (last) track
#' - `SeqLink`       -> stored in `plot$plot_links` (plot-level, deferred)
#' - `SeqAnnotation` -> stored in `plot$plot_annotations` (plot-level, deferred)
#'
#' When the LHS is a `SeqTrack`, only `SeqElement` and `SeqLink` are accepted.
#' Links added via `%+%` on a `SeqTrack` have their `t0` and `t1` locked to
#' the parent track's `track_id` — no override is permitted.
#'
#' @param e1 A `SeqPlot` or `SeqTrack` object.
#' @param e2 A `SeqTrack`, `SeqElement`, `SeqLink`, or `SeqAnnotation`.
#' @return `e1`, invisibly modified in place (R6 reference semantics).
#' @name op-plus
#' @usage e1 \%+\% e2
#' @export
`%+%` <- function(e1, e2) {

  # ── SeqTrack as LHS ─────────────────────────────────────────────────────────
  if (inherits(e1, "SeqTrack")) {
    if (!inherits(e2, c("SeqElement", "SeqLink")))
      stop("%+% cannot add object of class '", class(e2)[1], "' to a SeqTrack. ",
           "Only SeqElement or SeqLink objects can be added to a track.",
           call. = FALSE)
    if (inherits(e2, "SeqLink")) {
      # Lock t0/t1 to the parent track — no override permitted
      e2$t0 <- e1$track_id
      e2$t1 <- e1$track_id
    }
    e1$addElement(e2)
    return(invisible(e1))
  }

  # ── SeqPlot as LHS ──────────────────────────────────────────────────────────
  if (inherits(e1, "SeqPlot")) {

    # SeqTrack: add to layout
    if (inherits(e2, "SeqTrack")) {
      e1$addTrack(e2)
      return(invisible(e1))
    }

    # SeqLink: when neither t0 nor t1 is set, treat as a within-track link
    # added to the most recently added track (mirrors the SeqElement path).
    # Otherwise validate the references and store as a deferred plot-level link.
    if (inherits(e2, "SeqLink")) {
      if (is.null(e2$t0) && is.null(e2$t1)) {
        all_tracks <- e1$allTracks()
        n <- length(all_tracks)
        if (n == 0)
          stop("No tracks in plot yet. Add a seq_track() before adding links.",
               call. = FALSE)
        last_track <- all_tracks[[n]]
        e2$t0 <- last_track$track_id
        e2$t1 <- last_track$track_id
        last_track$addElement(e2)
        return(invisible(e1))
      }

      existing_ids <- e1$trackIds()
      for (ref in c(e2$t0, e2$t1)) {
        if (!is.null(ref) && !is.na(ref) && !ref %in% existing_ids) {
          stop("seq_link references track_id '", ref, "' which has not been added ",
               "to the plot yet. Define all referenced tracks before the link.",
               call. = FALSE)
        }
      }
      e1$plot_links <- append(e1$plot_links, list(e2))
      return(invisible(e1))
    }

    # SeqAnnotation: store deferred
    if (inherits(e2, "SeqAnnotation")) {
      e1$plot_annotations <- append(e1$plot_annotations, list(e2))
      return(invisible(e1))
    }

    # SeqElement: add to last track
    if (inherits(e2, "SeqElement")) {
      all_tracks <- e1$allTracks()
      n <- length(all_tracks)
      if (n == 0)
        stop("No tracks in plot yet. Add a seq_track() before adding elements.",
             call. = FALSE)
      all_tracks[[n]]$addElement(e2)
      return(invisible(e1))
    }

    # SeqBlank: no-op
    if (inherits(e2, "SeqBlank")) return(invisible(e1))

    stop("%+% cannot add object of class '", class(e2)[1], "' to a SeqPlot.",
         call. = FALSE)
  }

  stop("No %+% method for left-hand side of class '", class(e1)[1], "'.",
       call. = FALSE)
}

# %|% — convenience alias for direction = "right" -----------------------------

#' Append a track to the current row (horizontal stacking)
#'
#' Convenience alias. Equivalent to adding a [seq_track()] with
#' `direction = "right"`. All logic is in `%+%`.
#'
#' @param e1 A `SeqPlot` object.
#' @param e2 A `SeqTrack` (direction forced to `"right"`) or any valid `%+%` RHS.
#' @return `e1`, invisibly modified in place.
#' @name op-horizontal
#' @usage e1 \%|\% e2
#' @export
`%|%` <- function(e1, e2) {
  if (inherits(e2, "SeqTrack")) e2$direction <- "right"
  e1 %+% e2
}

# %__% — convenience alias for direction = "under" ----------------------------

#' Start a new row (vertical stacking)
#'
#' Convenience alias. Equivalent to adding a [seq_track()] with
#' `direction = "under"`. All logic is in `%+%`.
#'
#' @param e1 A `SeqPlot` object.
#' @param e2 A `SeqTrack` (direction forced to `"under"`) or any valid `%+%` RHS.
#' @return `e1`, invisibly modified in place.
#' @name op-vertical
#' @usage e1 \%__\% e2
#' @export
`%__%` <- function(e1, e2) {
  if (inherits(e2, "SeqTrack")) e2$direction <- "under"
  e1 %+% e2
}

#' Return a blank plot-level placeholder
#'
#' Used to occupy a named cell in a patchwork layout without rendering anything.
#'
#' @return A list with class `"SeqBlank"`.
#' @export
seq_blank <- function() {
  structure(list(), class = "SeqBlank")
}
