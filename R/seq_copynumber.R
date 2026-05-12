# ── seq_copynumber — single-sample copy-number scatter ───────────────────────
#
# Wrapper that builds a `seq_plot` with one track rendering per-bin
# log-ratio / ratio values as coloured points, with optional dotted
# reference lines at integer copy-number states and an optional segment
# overlay.

#' Default CN state colour palette
#'
#' Diverging palette for integer copy-number states: blue at CN=0, grey
#' at the diploid value (2), and increasingly warm colours above. Values
#' beyond the range fall back to `"grey80"`.
#'
#' @return A named character vector keyed by CN state (as a string).
#' @keywords internal
.default_cn_state_colors <- function() {
  c("0" = "#205EA6",
    "1" = "#66A0C8",
    "2" = "#B7B5AC",
    "3" = "#F9AE77",
    "4" = "#DA702C",
    "5" = "#D14D41",
    "6" = "#AF3029",
    "7" = "#8B1C3B",
    "8" = "#A02C6D")
}

#' Build a `SeqMap` from a named list of column names
#'
#' Helper used by the wrapper functions to construct mappings
#' programmatically — `map()` captures unevaluated expressions, which
#' is awkward when the column name is determined at runtime. This
#' builds the same object shape directly from string column names (or
#' already-unquoted symbols).
#'
#' @param ... Named arguments. Strings are converted to bare R symbols;
#'   language / call objects are passed through unchanged.
#' @return A `SeqMap` object.
#' @keywords internal
.map_from_names <- function(...) {
  xs <- list(...)
  out <- lapply(xs, function(v) {
    if (is.character(v)) as.name(v)
    else if (is.language(v) || is.symbol(v)) v
    else v
  })
  structure(out, class = "SeqMap")
}

#' Detect a column name from a list of candidates
#'
#' @param gr A `GRanges` object.
#' @param hint Optional explicit column name.
#' @param candidates Character vector of candidate names, in preference
#'   order.
#' @param predicate A function `function(x)` accepting an mcols column
#'   and returning `TRUE` when the column is usable.
#' @param what Label for the error/warning message (e.g. `"CN state"`).
#' @return The detected column name.
#' @keywords internal
.detect_col <- function(gr, hint, candidates, predicate, what) {
  mc  <- S4Vectors::mcols(gr)
  nms <- names(mc)
  if (!is.null(hint)) {
    if (!hint %in% nms)
      stop("Column '", hint, "' not found in data mcols.", call. = FALSE)
    return(hint)
  }
  found <- intersect(candidates, nms)
  if (length(found) > 0L) return(found[1])
  usable <- nms[vapply(nms, function(n) predicate(mc[[n]]), logical(1))]
  if (length(usable) == 0L)
    stop("No column suitable for ", what, " found in data.", call. = FALSE)
  warning("No standard ", what, " column found; using '", usable[1], "'.",
          call. = FALSE)
  usable[1]
}

#' Copy-number scatter plot
#'
#' Builds a [seq_plot()] with one track showing per-bin copy-number data:
#' a scatter of continuous ratio / log-ratio values, coloured by integer
#' CN state. Optional dotted reference lines mark integer CN values
#' present in the data, and an optional segmentation overlay draws a
#' horizontal line segment for each called segment.
#'
#' The CN state and ratio column names are auto-detected when not
#' specified. For CN state the search order is
#' `c("cn", "copy_number", "CN", "state", "integer_cn")`; for ratio it
#' is `c("log2ratio", "logR", "log2R", "ratio", "log2_ratio")`. When no
#' match is found the first integer / numeric column is used with a
#' warning.
#'
#' @param data A `GRanges` with one row per bin and mcols columns for
#'   the CN state and the continuous ratio.
#' @param windows A `GRanges` of genomic windows defining the view.
#' @param cn_col Name of the mcols column giving integer CN state. Auto-
#'   detected when `NULL`.
#' @param ratio_col Name of the mcols column giving the continuous
#'   ratio / log-ratio. Auto-detected when `NULL`.
#' @param state_colors Named character vector keyed by CN state (as a
#'   string). Defaults to a diverging blue -> grey -> orange palette.
#' @param segment_data Optional `GRanges` of segmentation calls. Drawn as
#'   horizontal line segments over the scatter.
#' @param segment_col Column in `segment_data` giving the per-segment y
#'   value (e.g. segment mean). Auto-detected from the same candidates as
#'   `ratio_col` when `NULL`.
#' @param show_reference_lines Logical; draw dotted horizontal lines at
#'   each integer CN value present in the data range. Default `TRUE`.
#' @param reference_line_col Colour for the reference lines.
#' @param track_height Relative track height.
#' @param track_id Character `track_id` for the generated track. Defaults
#'   to `"copynumber"`.
#' @param legend A `LegendKey` or `SeqLegendSpec` attached to the scatter
#'   element. `NULL` (default) produces no legend entry.
#' @param show_legend Logical. When `FALSE`, the scatter element contributes no
#'   legend. Default `TRUE`.
#' @param ... Additional arguments forwarded to [seq_track()].
#' @return A `SeqPlot` object composable via `%+%`, `%|%`, `%__%`, and
#'   [seq_resolve()].
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(seq(1, 1e6, by = 1e4), width = 5000),
#'               cn        = sample(0:4, 100, replace = TRUE),
#'               log2ratio = rnorm(100, 0, 0.3))
#' win <- GRanges("chr1", IRanges(1, 1e6))
#' seq_copynumber(gr, windows = win)
#' @export
seq_copynumber <- function(data,
                           windows,
                           cn_col               = NULL,
                           ratio_col            = NULL,
                           state_colors         = NULL,
                           segment_data         = NULL,
                           segment_col          = NULL,
                           show_reference_lines = TRUE,
                           reference_line_col   = "grey50",
                           track_height         = 1,
                           track_id             = NULL,
                           legend               = NULL,
                           show_legend          = TRUE,
                           ...) {
  .stop_if_not_granges(data,    "data")
  .stop_if_not_granges(windows, "windows")

  cn_col <- .detect_col(
    data, cn_col,
    candidates = c("cn", "copy_number", "CN", "state", "integer_cn"),
    predicate  = function(x) is.numeric(x) || is.integer(x),
    what       = "CN state"
  )
  ratio_col <- .detect_col(
    data, ratio_col,
    candidates = c("log2ratio", "logR", "log2R", "ratio", "log2_ratio"),
    predicate  = function(x) is.numeric(x),
    what       = "ratio"
  )

  palette   <- state_colors %||% .default_cn_state_colors()
  cn_vals   <- as.character(S4Vectors::mcols(data)[[cn_col]])
  color_vec <- palette[cn_vals]
  color_vec[is.na(color_vec)] <- "grey80"
  S4Vectors::mcols(data)$.point_color <- unname(color_vec)

  p <- seq_plot() %+%
    seq_track(
      data         = data,
      mapping      = .map_from_names(x = "start", y = ratio_col),
      windows      = windows,
      track_height = track_height,
      track_id     = track_id %||% "copynumber",
      ...
    ) %+%
    seq_point(
      mapping     = .map_from_names(color = ".point_color"),
      aesthetics  = aes(size = 0.4),
      legend      = legend,
      show_legend = show_legend
    )

  if (isTRUE(show_reference_lines)) {
    cn_int  <- suppressWarnings(as.integer(cn_vals))
    present <- sort(unique(cn_int[is.finite(cn_int)]))
    if (length(present) > 0L) {
      nwin   <- length(windows)
      starts <- BiocGenerics::start(windows)
      ends   <- BiocGenerics::end(windows)
      seqs   <- as.character(GenomicRanges::seqnames(windows))
      ref_gr <- GenomicRanges::GRanges(
        seqnames = rep(seqs,   each = length(present)),
        ranges   = IRanges::IRanges(
          start  = rep(starts, each = length(present)),
          end    = rep(ends,   each = length(present))
        ),
        ref_y    = rep(as.numeric(present), nwin)
      )
      p <- p %+% seq_segment(
        data       = ref_gr,
        mapping    = .map_from_names(x = "start", x_end = "end",
                                     y = "ref_y", y_end = "ref_y"),
        aesthetics = aes(color     = reference_line_col,
                         linetype  = "dotted",
                         linewidth = 0.4)
      )
    }
  }

  if (!is.null(segment_data)) {
    .stop_if_not_granges(segment_data, "segment_data")
    segment_col <- .detect_col(
      segment_data, segment_col,
      candidates = c("seg_mean", "seg.mean", "mean", "log2ratio",
                     "logR", "log2R", "ratio"),
      predicate  = function(x) is.numeric(x),
      what       = "segment mean"
    )
    p <- p %+% seq_segment(
      data       = segment_data,
      mapping    = .map_from_names(x = "start", x_end = "end",
                                   y = segment_col, y_end = segment_col),
      aesthetics = aes(color = "black", linewidth = 1.5)
    )
  }

  p
}
