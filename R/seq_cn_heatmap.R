# ── seq_cn_heatmap — multi-sample CN tile heatmap ────────────────────────────
#
# Wrapper that arranges per-sample CN calls into a heatmap: samples on
# the y-axis, genome on the x-axis, tile fill encoding CN state. Uses
# `seq_tile()` with `rotate = FALSE` and a discrete y-axis.

#' Multi-sample copy-number heatmap
#'
#' Builds a [seq_plot()] with a single track arranging per-sample CN
#' calls as a heatmap: samples are placed on the y-axis (top to bottom
#' according to `sample_order`), genomic windows along the x-axis, and
#' each bin is a tile coloured by CN state.
#'
#' Data may be supplied either as a long-format `GRanges` with one row
#' per (sample, bin) carrying `sample` and `cn` mcols columns, or as a
#' numeric matrix (rows = samples, cols = bins) plus a `bins` `GRanges`
#' giving the genomic position of each column (passed via `...`).
#' Matrix input requires the matrix argument to have sample names as
#' `rownames`.
#'
#' @param data A `GRanges` (long format) or numeric matrix.
#' @param windows A `GRanges` of genomic view windows.
#' @param sample_col Name of the mcols column giving sample identity.
#'   Auto-detected from `c("sample", "sample_id", "Sample", "id")` when
#'   `NULL`.
#' @param cn_col Name of the mcols column giving integer CN state. Auto-
#'   detected from the same candidates as [seq_copynumber()].
#' @param state_colors Named character vector keyed by CN state string.
#'   Defaults to the [seq_copynumber()] palette.
#' @param sample_order Character vector of sample names in display order
#'   (top to bottom). When `NULL`, samples are sorted alphabetically.
#' @param bins Optional `GRanges` giving the genomic position of each
#'   matrix column; required when `data` is a matrix.
#' @param track_height Relative track height.
#' @param track_id Character `track_id` for the generated track.
#' @param ... Additional arguments forwarded to [seq_track()].
#' @return A `SeqPlot` object.
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges("chr1", IRanges(seq(1, 1e6, by = 2e4), width = 2e4),
#'               sample = rep(c("S1","S2","S3"), length.out = 50),
#'               cn     = sample(0:4, 50, replace = TRUE))
#' win <- GRanges("chr1", IRanges(1, 1e6))
#' seq_cn_heatmap(gr, windows = win)
#' @export
seq_cn_heatmap <- function(data,
                           windows,
                           sample_col   = NULL,
                           cn_col       = NULL,
                           state_colors = NULL,
                           sample_order = NULL,
                           bins         = NULL,
                           track_height = 3,
                           track_id     = NULL,
                           ...) {
  .stop_if_not_granges(windows, "windows")

  if (is.matrix(data)) {
    if (is.null(bins))
      stop("When 'data' is a matrix, 'bins' (a GRanges matching matrix columns) ",
           "must be supplied.", call. = FALSE)
    .stop_if_not_granges(bins, "bins")
    if (ncol(data) != length(bins))
      stop("ncol(data) must equal length(bins).", call. = FALSE)
    sample_names <- rownames(data) %||% paste0("S", seq_len(nrow(data)))
    long_gr <- rep(bins, nrow(data))
    S4Vectors::mcols(long_gr)$sample <- rep(sample_names, each = ncol(data))
    S4Vectors::mcols(long_gr)$cn     <- as.integer(t(data))
    data        <- long_gr
    sample_col  <- sample_col %||% "sample"
    cn_col      <- cn_col     %||% "cn"
  }

  .stop_if_not_granges(data, "data")

  sample_col <- .detect_col(
    data, sample_col,
    candidates = c("sample", "sample_id", "Sample", "id"),
    predicate  = function(x) is.character(x) || is.factor(x) || is.integer(x),
    what       = "sample identity"
  )
  cn_col <- .detect_col(
    data, cn_col,
    candidates = c("cn", "copy_number", "CN", "state", "integer_cn"),
    predicate  = function(x) is.numeric(x) || is.integer(x),
    what       = "CN state"
  )

  palette   <- state_colors %||% .default_cn_state_colors()
  sample_v  <- as.character(S4Vectors::mcols(data)[[sample_col]])
  sample_lv <- sample_order %||% sort(unique(sample_v))

  y_idx <- match(sample_v, sample_lv)
  cn_v  <- as.character(S4Vectors::mcols(data)[[cn_col]])
  fill  <- palette[cn_v]
  fill[is.na(fill)] <- "grey80"

  S4Vectors::mcols(data)$.y_idx      <- y_idx
  S4Vectors::mcols(data)$.fill_color <- unname(fill)

  keep <- !is.na(y_idx)
  data <- data[keep]

  p <- seq_plot() %+%
    seq_track(
      data         = data,
      mapping      = .map_from_names(x = "start",
                                     y = ".y_idx",
                                     fill = ".fill_color"),
      windows      = windows,
      track_height = track_height,
      track_id     = track_id %||% "cn_heatmap",
      scale_y      = seq_scale_discrete(levels = sample_lv),
      ...
    ) %+%
    seq_tile(aesthetics = aes(rotate = FALSE))

  p
}
