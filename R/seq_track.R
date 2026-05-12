# ── SeqTrackR6 ────────────────────────────────────────────────────────────────

#' SeqTrack R6 class
#'
#' Internal R6 class backing [seq_track()]. Public users should go through the
#' snake_case constructor; this class is exported only so that other internal
#' files can reference it.
#'
#' @keywords internal
SeqTrackR6 <- R6::R6Class("SeqTrack",
  public = list(
    #' @field data GRanges default data for elements in this track.
    data                 = NULL,
    #' @field mapping SeqMap default mapping for elements in this track.
    mapping              = NULL,
    #' @field aesthetics SeqAes default constant aesthetics.
    aesthetics           = NULL,
    #' @field windows GRanges windows for this track (only place windows are set).
    windows              = NULL,
    #' @field track_id Character unique identifier for this track.
    track_id             = NULL,
    #' @field track_width Relative width unit within its row.
    track_width          = 1,
    #' @field track_height Relative height unit within its column.
    track_height         = 1,
    #' @field direction One of "right" or "under".
    direction            = "right",
    #' @field elements List of SeqElement / SeqLink added via `%+%`.
    elements             = list(),
    #' @field scale_x Optional position scale for the primary x axis.
    scale_x              = NULL,
    #' @field scale_y Optional position scale for the primary y axis.
    scale_y              = NULL,
    #' @field scale_x2 Optional position scale for the secondary x axis.
    #'   Populated automatically at `layoutGrid()` time from elements
    #'   whose `map(axis.x = 2)` targets the secondary x axis.
    scale_x2             = NULL,
    #' @field scale_y2 Optional position scale for the secondary y axis.
    scale_y2             = NULL,
    #' @field uses_genomic_y TRUE when mapping$y resolves to a genomic special.
    uses_genomic_y       = FALSE,
    #' @field y_windows Optional GRanges windows for the primary y axis.
    y_windows            = NULL,
    #' @field y_windows2 Optional GRanges windows for the secondary y axis.
    y_windows2           = NULL,
    #' @field has_axis_x2 TRUE when any element targets axis.x = 2, when
    #'   scale_x2 is explicitly set, or when the theme forces visibility.
    #'   Set at layoutGrid() time.
    has_axis_x2          = FALSE,
    #' @field has_axis_y2 Symmetric to `has_axis_x2` for the y direction.
    has_axis_y2          = FALSE,
    #' @field resolved_theme Populated at layoutGrid(): nested list with
    #'   per-axis specs (x1/x2/y1/y2), track chrome, and per-window flag.
    resolved_theme       = NULL,
    #' @field track_outer_margin Named list of track-level outer margin (npc)
    #'   — reserves the band where axis titles sit.
    track_outer_margin   = NULL,
    #' @field track_inner_margin Named list of track-level inner margin (npc)
    #'   — separates axis titles from the window row.
    track_inner_margin   = NULL,
    #' @field window_outer_margin Named list of per-window outer margin (npc).
    window_outer_margin  = NULL,
    #' @field window_inner_margin Named list of per-window inner margin (npc)
    #'   — holds per-window axis ticks and labels.
    window_inner_margin  = NULL,
    #' @field window_margin Deprecated. Set
    #'   `aes("window.gap.width" = <value>)` on the plot or track instead.
    #'   Stored for backward compatibility but has no effect on layout.
    window_margin        = NULL,
    #' @field combine_windows Logical. When `TRUE`, multi-region windows
    #'   are concatenated into a single virtual panel, with per-original-
    #'   window axis labels and a thin separator at each window boundary.
    #'   Used to draw cross-window data (e.g. inter-chromosomal Hi-C
    #'   contacts) within one continuous track.
    combine_windows      = FALSE,
    #' @field combine_y_windows Logical. Symmetric to `combine_windows`
    #'   for the y-axis (genomic y tracks only).
    combine_y_windows    = FALSE,
    #' @field flip_x Logical. When `TRUE`, mirror the x-axis: low data
    #'   values render at the right edge of the panel and high values
    #'   at the left. Tick labels follow the same orientation.
    flip_x               = FALSE,
    #' @field flip_y Logical. Symmetric to `flip_x` for the y-axis.
    #'   For Hi-C `triangle` style this produces a downward-pointing
    #'   triangle; for `diagonal` it shows the lower diagonal.
    flip_y               = FALSE,
    #' @field show_legend Logical. When `FALSE`, this track contributes no
    #'   legend keys regardless of the `legend` fields on its elements.
    #'   Default `TRUE`.
    show_legend          = TRUE,
    #' @field window_scale Numeric vector or `NULL`. Per-window x-axis scale
    #'   factors (e.g. `1e-6` for Mb, `1e-3` for kb, `1` for bp). When `NULL`
    #'   (default), the scale is inferred from the narrowest window. When length
    #'   1, the value is applied to all windows. When length equals the number of
    #'   windows, values are applied positionally. Any other length triggers a
    #'   warning and recycles with `rep_len()`.
    window_scale         = NULL,

    #' @description Construct a SeqTrackR6.
    #' @param data GRanges default data.
    #' @param mapping SeqMap default mapping.
    #' @param aesthetics SeqAes default constant aesthetics.
    #' @param windows GRanges windows for this track.
    #' @param track_id Character unique identifier.
    #' @param direction One of "right" or "under".
    #' @param track_width Relative width unit.
    #' @param track_height Relative height unit.
    #' @param scale_x Optional primary x position scale.
    #' @param scale_y Optional primary y position scale.
    #' @param scale_x2 Optional secondary x position scale.
    #' @param scale_y2 Optional secondary y position scale.
    #' @param y_windows Optional GRanges y-axis windows for the primary axis.
    #' @param y_windows2 Optional GRanges y-axis windows for the secondary axis.
    #' @param track_outer_margin Scalar, length-4 `c(bottom, left, top,
    #'   right)`, or named list. Default `0.02`. Band where axis titles
    #'   (derived from the track `mapping`) are drawn.
    #' @param track_inner_margin Same form. Default `0.02`. Separates the
    #'   title zone from the window row.
    #' @param window_outer_margin Same form. Default `0`. Spacer around
    #'   each window inside the track plot region.
    #' @param window_inner_margin Same form. Default `0.02`. Holds each
    #'   window's axis ticks and labels.
    #' @param window_margin Deprecated. Use
    #'   `aes("window.gap.width" = <value>)` on the plot or track aesthetics
    #'   to control inter-window gap. Passing a non-`NULL` value emits a
    #'   deprecation warning and is otherwise ignored.
    #' @param elements Optional list of `SeqElement` / `SeqLink` objects to
    #'   pre-populate this track. Elements can also be added later via
    #'   `addElement()` or the `%+%` operator.
    #' @param show_legend Logical. When `FALSE`, this track contributes no
    #'   legend keys regardless of element `legend` fields. Default `TRUE`.
    #' @param window_scale Numeric vector or `NULL`. Per-window x-axis scale
    #'   factors (e.g. `1e-6` for Mb, `1e-3` for kb, `1` for bp). When `NULL`
    #'   (default), the scale is inferred from the narrowest window. Length 1
    #'   applies the value to all windows; length equal to the number of windows
    #'   uses values positionally; any other length triggers a warning and
    #'   recycles with `rep_len()`.
    #' @param ... Additional arguments (currently ignored).
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          windows = NULL, track_id = NULL,
                          direction = "right",
                          track_width = 1, track_height = 1,
                          scale_x = NULL, scale_y = NULL,
                          scale_x2 = NULL, scale_y2 = NULL,
                          y_windows = NULL, y_windows2 = NULL,
                          track_outer_margin = 0.02,
                          track_inner_margin = 0.02,
                          window_outer_margin = 0,
                          window_inner_margin = 0.02,
                          window_margin = NULL,
                          combine_windows = FALSE,
                          combine_y_windows = FALSE,
                          flip_x = FALSE, flip_y = FALSE,
                          elements = list(),
                          show_legend = TRUE,
                          window_scale = NULL, ...) {
      self$data                 <- data
      self$mapping              <- mapping
      self$aesthetics           <- aesthetics
      self$windows              <- windows
      self$track_id             <- track_id
      self$direction            <- direction
      self$track_width          <- track_width
      self$track_height         <- track_height
      self$scale_x              <- scale_x
      self$scale_y              <- scale_y
      self$scale_x2             <- scale_x2
      self$scale_y2             <- scale_y2
      self$y_windows            <- y_windows
      self$y_windows2           <- y_windows2
      # Auto-detect genomic y when user supplies y_windows or a genomic y scale.
      if (!is.null(y_windows)) {
        self$uses_genomic_y <- TRUE
        # Auto-build a genomic scale_y so axis drawing has something to
        # render against. Users can override by passing scale_y.
        if (is.null(self$scale_y))
          self$scale_y <- seq_scale_genomic(y_windows)
      }
      if (inherits(scale_y, "SeqScaleGenomic")) {
        self$uses_genomic_y <- TRUE
        if (is.null(self$y_windows)) self$y_windows <- scale_y$windows
      }
      self$track_outer_margin   <- .normalize_margin(track_outer_margin)
      self$track_inner_margin   <- .normalize_margin(track_inner_margin)
      self$window_outer_margin  <- .normalize_margin(window_outer_margin)
      self$window_inner_margin  <- .normalize_margin(window_inner_margin)
      self$window_margin        <- window_margin
      if (!is.null(window_margin)) {
        warning(
          "`window_margin` is deprecated and has no effect. ",
          "Use aes(\"window.gap.width\" = <value>) instead.",
          call. = FALSE
        )
      }
      self$combine_windows      <- isTRUE(combine_windows)
      self$combine_y_windows    <- isTRUE(combine_y_windows)
      self$flip_x               <- isTRUE(flip_x)
      self$flip_y               <- isTRUE(flip_y)
      if (length(elements) > 0L) self$elements <- elements
      self$show_legend          <- isTRUE(show_legend)
      self$window_scale         <- window_scale
    },

    #' @description Append an element (SeqElement or SeqLink) to this track.
    #' @param elem The element to add.
    addElement = function(elem) {
      self$elements <- append(self$elements, list(elem))
      invisible(self)
    },

    #' @description
    #' Collect legend keys from all elements in this track.
    #'
    #' Iterates every element (including `SeqLink` subclasses) and calls
    #' `$collect_legend_keys()` on each. Returns a flat list of entries, where
    #' each entry is a named list with fields `title`, `key`, and
    #' `element_class` (as produced by `SeqElement$collect_legend_keys()`).
    #' Returns `NULL` when `show_legend` is `FALSE` or no element contributes
    #' any keys.
    #'
    #' @return A list of legend-key entries, or `NULL`.
    collect_legend_keys = function() {
      if (!isTRUE(self$show_legend)) return(NULL)

      out <- list()
      for (el in self$elements) {
        if (is.function(el$collect_legend_keys)) {
          keys <- el$collect_legend_keys()
          if (!is.null(keys)) out <- c(out, keys)
        }
      }

      if (length(out) == 0L) return(NULL)
      out
    }
  )
)

# ── seq_track() constructor ──────────────────────────────────────────────────

#' Create a new seq_track
#'
#' A track cell is partitioned into five nested zones:
#' \enumerate{
#'   \item **Track outer margin** — outermost band on the track cell;
#'     reserved for axis titles derived from the track `mapping`.
#'   \item **Track inner margin** — separates titles from the window row.
#'   \item **Window outer margin** — per-window spacer inside the track
#'     plot region.
#'   \item **Window inner margin** — per-window band holding axis lines,
#'     ticks, and tick labels.
#'   \item **Plot area** — where elements render.
#' }
#'
#' Axis lines sit at the boundary between the plot area and the window
#' inner margin. Tick labels and axis titles pick up their text from the
#' track's `mapping`: the x-axis title is the expression assigned to `x`,
#' and the y-axis title is the expression assigned to `y`.
#'
#' @param data A GRanges object providing the default data for elements in this track.
#' @param mapping A SeqMap object from [map()] providing default aesthetic mappings.
#' @param aesthetics A SeqAes object from [aes()] providing constant aesthetics.
#' @param windows A GRanges object defining the genomic windows for this track.
#'   Elements cannot define their own windows — this is the only place windows are set.
#' @param track_id Character. Unique identifier for this track, used for patchwork
#'   layout matching and cross-track link references.
#' @param direction One of `"right"` (append to current row) or `"under"` (start new row).
#'   Ignored for the first track (always top-left) and when [seq_plot()] is given a
#'   layout string.
#' @param track_width Relative width unit within its row. Default 1.
#' @param track_height Relative height unit within its column. Default 1.
#' @param scale_x Optional position scale for the primary x axis — a
#'   `seq_scale_genomic()` (default behaviour when `NULL`),
#'   `seq_scale_continuous(limits = ...)` for scalar data x, or
#'   `seq_scale_discrete(levels = ...)` for categorical x.
#' @param scale_y Optional position scale for the primary y axis. Pass
#'   `seq_scale_genomic(y_windows)` to flip the track so genomic
#'   position runs along y; continuous / discrete scales are supported too.
#' @param scale_x2,scale_y2 Optional position scales for the secondary x
#'   and y axes. Elements routed to the secondary axis via
#'   `map(axis.x = 2)` or `map(axis.y = 2)` are positioned against
#'   these scales. When `NULL`, secondary scales are auto-inferred
#'   from the contributing elements at `layoutGrid()` time.
#' @param y_windows Optional `GRanges` used as the genomic y-axis range
#'   when flipping the track. Setting this (or a `SeqScaleGenomic`
#'   `scale_y`) auto-enables `uses_genomic_y`.
#' @param y_windows2 Optional `GRanges` for a secondary genomic y axis.
#' @param track_outer_margin Per-track outer margin in npc units. Scalar,
#'   length-4 `c(bottom, left, top, right)` (base-R `par(mar = ...)`
#'   order), or a named list with any of `top`, `right`, `bottom`,
#'   `left`. Default `0.02`. Axis titles draw here.
#' @param track_inner_margin Per-track inner margin. Same form as
#'   `track_outer_margin`. Default `0.02`.
#' @param window_outer_margin Per-window outer margin. Same form. Default
#'   `0`. Optional spacer around each window.
#' @param window_inner_margin Per-window inner margin. Same form. Default
#'   `0.02`. Axis ticks and tick labels draw here.
#' @param window_margin Deprecated. Use
#'   `aes("window.gap.width" = <value>)` on the plot or track aesthetics
#'   to control the inter-window gap. Passing a non-`NULL` value emits a
#'   deprecation warning and is otherwise ignored.
#' @param combine_windows Logical. When `TRUE`, multi-region `windows`
#'   are concatenated into a single virtual panel — useful for drawing
#'   data that spans regions (e.g. inter-chromosomal Hi-C contacts) in
#'   one continuous track. Each original window's axis labels and
#'   title are still drawn, separated by a thin boundary marker.
#'   Default `FALSE` (multi-region windows render as separate panels).
#' @param combine_y_windows Logical. Symmetric to `combine_windows`
#'   for the genomic y-axis (only relevant for tracks with multiple
#'   `y_windows`). Default `FALSE`.
#' @param flip_x,flip_y Logical. When `TRUE`, mirror the x or y axis
#'   so low data values render at the high end of the panel and vice
#'   versa. Tick labels follow the same orientation.
#' @param elements Optional list of `SeqElement` / `SeqLink` objects to
#'   pre-populate this track. Elements can also be added via `addElement()`
#'   or the `%+%` operator.
#' @param show_legend Logical. When `FALSE`, this track contributes no legend
#'   keys regardless of element `legend` fields. Default `TRUE`.
#' @param window_scale Numeric vector or `NULL`. Per-window x-axis scale
#'   factors (e.g. `1e-6` for Mb, `1e-3` for kb, `1` for bp). When `NULL`
#'   (default), the scale is inferred from the narrowest window. Length 1
#'   applies the value to all windows; length equal to the number of windows
#'   uses values positionally; any other length triggers a warning and
#'   recycles with `rep_len()`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqTrackR6` instance (S3 class `"SeqTrack"`).
#' @examples
#' seq_track(track_id = "A")
#' @export
seq_track <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                      windows = NULL, track_id = NULL,
                      direction = "right",
                      track_width = 1, track_height = 1,
                      scale_x = NULL, scale_y = NULL,
                      scale_x2 = NULL, scale_y2 = NULL,
                      y_windows = NULL, y_windows2 = NULL,
                      track_outer_margin = 0.02,
                      track_inner_margin = 0.02,
                      window_outer_margin = 0,
                      window_inner_margin = 0.02,
                      window_margin = NULL,
                      combine_windows = FALSE,
                      combine_y_windows = FALSE,
                      flip_x = FALSE, flip_y = FALSE,
                      elements = list(),
                      show_legend = TRUE,
                      window_scale = NULL, ...) {
  SeqTrackR6$new(
    data = data, mapping = mapping, aesthetics = aesthetics,
    windows = windows, track_id = track_id, direction = direction,
    track_width = track_width, track_height = track_height,
    scale_x = scale_x, scale_y = scale_y,
    scale_x2 = scale_x2, scale_y2 = scale_y2,
    y_windows = y_windows, y_windows2 = y_windows2,
    track_outer_margin  = track_outer_margin,
    track_inner_margin  = track_inner_margin,
    window_outer_margin = window_outer_margin,
    window_inner_margin = window_inner_margin,
    window_margin       = window_margin,
    combine_windows     = combine_windows,
    combine_y_windows   = combine_y_windows,
    flip_x              = flip_x,
    flip_y              = flip_y,
    elements            = elements,
    show_legend         = show_legend,
    window_scale        = window_scale, ...
  )
}
