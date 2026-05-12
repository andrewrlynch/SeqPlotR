# ── seq_resolve — compose wrapper-built seq_plot objects into a parent ──────
#
# The key composition helper for the wrapper family. A call to
# `seq_copynumber()` / `seq_hic()` / `seq_chip()` / etc. returns a
# `seq_plot`. `seq_resolve()` extracts the tracks (and any plot-level
# links/annotations) from those child plots and appends them to a
# parent `seq_plot`, preserving operator-chain semantics.

#' Compose wrapper-built plots into a parent plot
#'
#' Each wrapper (e.g. [seq_copynumber()], [seq_hic()], [seq_chip()])
#' returns a `seq_plot`. `seq_resolve()` unpacks those child plots into
#' a parent `seq_plot`, letting you combine heterogeneous track stacks
#' in a single figure.
#'
#' Tracks from each child are appended to `parent` respecting their
#' original row grouping. The `direction` argument controls where the
#' **first** row of each child lands relative to the current bottom of
#' `parent`:
#' \itemize{
#'   \item `"under"` — start a new row below whatever is already in
#'     the parent (the common case).
#'   \item `"right"` — append to the current row (put the child beside
#'     what is already there).
#' }
#'
#' Subsequent rows within a child are always placed `"under"` (i.e. the
#' child's internal row structure is preserved). Multiple children
#' supplied in one call are stacked in order, each with `direction`
#' applied between children.
#'
#' Plot-level links and annotations carried by the children are
#' transferred to `parent` in document order.
#'
#' Duplicate `track_id`s across children cause an error — pass unique
#' `track_id` values to each wrapper call to resolve the conflict.
#'
#' @param parent A `seq_plot` object to add tracks into.
#' @param ... One or more `seq_plot` objects produced by wrapper
#'   functions.
#' @param direction `"under"` (default) or `"right"`. See Details.
#' @return The modified `parent`, invisibly (R6 reference semantics).
#' @examples
#' \dontrun{
#' win <- GRanges("chr1", IRanges(1, 1e6))
#' cn  <- seq_copynumber(cn_gr, windows = win, track_id = "CN")
#' hic <- seq_hic(hic_gr, windows = win, style = "triangle",
#'                track_id = "HiC")
#' seq_resolve(seq_plot(), cn, hic)$plot()
#' }
#' @export
seq_resolve <- function(parent, ..., direction = "under") {
  if (!inherits(parent, "SeqPlot"))
    stop("parent must be a seq_plot object.", call. = FALSE)
  direction <- match.arg(direction, c("under", "right"))

  children <- list(...)
  if (length(children) == 0L) return(invisible(parent))

  for (ci in seq_along(children)) {
    child <- children[[ci]]
    if (!inherits(child, "SeqPlot"))
      stop("All ... arguments to seq_resolve() must be seq_plot objects ",
           "(child #", ci, " is '", class(child)[1], "').", call. = FALSE)

    # Extract the child's row structure. Positional children have
    # $rows; patchwork children keep tracks flat in $tracks and the
    # layout string decides placement. For patchwork children we
    # flatten into a single row list so each track lands one per row.
    child_rows <- child$rows
    if (is.null(child_rows))
      child_rows <- lapply(child$tracks %||% list(), list)

    existing_ids <- parent$trackIds()
    existing_ids <- existing_ids[!is.na(existing_ids)]

    first_row <- TRUE
    for (row in child_rows) {
      if (length(row) == 0L) next
      for (ti in seq_along(row)) {
        trk <- row[[ti]]
        tid <- trk$track_id
        if (!is.null(tid) && !is.na(tid) && tid %in% existing_ids)
          stop("seq_resolve(): duplicate track_id '", tid, "' from child #", ci,
               ". Pass a unique 'track_id' to each wrapper call to resolve.",
               call. = FALSE)
        if (!is.null(tid) && !is.na(tid))
          existing_ids <- c(existing_ids, tid)

        # First track of the very first row inherits `direction`;
        # subsequent track rows within this child always go "under";
        # tracks within a row go "right".
        if (first_row && ti == 1L) {
          trk$direction <- direction
        } else if (ti == 1L) {
          trk$direction <- "under"
        } else {
          trk$direction <- "right"
        }

        parent <- parent %+% trk
      }
      first_row <- FALSE
    }

    if (length(child$plot_links) > 0L)
      parent$plot_links <- append(parent$plot_links, child$plot_links)
    if (length(child$plot_annotations) > 0L)
      parent$plot_annotations <- append(parent$plot_annotations,
                                        child$plot_annotations)
  }

  invisible(parent)
}
