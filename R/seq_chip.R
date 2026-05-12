# ── seq_chip — ChIP-style multi-track signal + peaks wrapper ────────────────
#
# Stacks one `seq_area` (coverage signal) track plus an optional
# `seq_bar` (peak calls) track per sample. Samples are placed
# vertically; each sample inherits a colour from `colors` or the
# `flexoki_palette()`.

#' Detect the signal column in a GRanges
#'
#' @param gr A `GRanges` with a numeric signal column.
#' @param hint Optional explicit column name.
#' @return The detected column name.
#' @keywords internal
.detect_signal_col <- function(gr, hint = NULL) {
  mc  <- S4Vectors::mcols(gr)
  nms <- names(mc)
  if (!is.null(hint)) {
    if (!hint %in% nms)
      stop("Column '", hint, "' not found in signal GRanges.", call. = FALSE)
    return(hint)
  }
  candidates <- c("score", "signal", "coverage", "count",
                  "rpm", "rpkm", "fpkm", "cpm")
  found <- intersect(candidates, nms)
  if (length(found) > 0L) return(found[1])
  num_cols <- nms[vapply(nms, function(n) is.numeric(mc[[n]]), logical(1))]
  if (length(num_cols) == 0L)
    stop("No numeric signal column found in data for sample.", call. = FALSE)
  warning("No standard signal column found; using '", num_cols[1], "'",
          call. = FALSE)
  num_cols[1]
}

#' Split a GRanges by sample column into a named list
#'
#' @keywords internal
.split_gr_by_sample <- function(gr, sample_col) {
  vals <- as.character(S4Vectors::mcols(gr)[[sample_col]])
  lvls <- unique(vals)
  out  <- lapply(lvls, function(v) gr[vals == v])
  stats::setNames(out, lvls)
}

#' ChIP-style multi-track signal + peaks plot
#'
#' For each sample, `seq_chip()` stacks a [seq_area()] coverage track
#' on top of an optional [seq_bar()] peak track, colouring both from
#' the sample's colour. Tracks flow top-to-bottom (each one uses
#' `direction = "under"`), so the returned `SeqPlot` is a single
#' column of stacked tracks.
#'
#' @param data Either a named `list` of `GRanges` (one per sample), or a
#'   single `GRanges` with a sample column (pass its name via
#'   `sample_col`). Each signal `GRanges` must carry a numeric signal
#'   column (auto-detected).
#' @param windows `GRanges` defining the view region.
#' @param sample_col Column in `data` giving sample identity when
#'   `data` is a single `GRanges`. Ignored for list input.
#' @param signal_col Explicit signal column name; auto-detected per
#'   sample when `NULL`.
#' @param peaks Optional peak calls, mirroring `data`: either a named
#'   list of `GRanges` with the same names, or a single `GRanges` with
#'   `sample_col`.
#' @param peak_col Optional column in peaks used for bar height; default
#'   renders uniform-height bars.
#' @param colors Named character vector mapping sample name to colour.
#'   Defaults to cycling the [flexoki_palette()].
#' @param scale_max Numeric scalar or named vector capping the signal y-
#'   axis per sample. `NULL` (default) autoscales.
#' @param signal_height Relative height of each signal track.
#' @param peak_height Relative height of each peak track.
#' @param show_genes Optional `GRanges` for gene annotation — adds a
#'   final [seq_gene()] track beneath the sample tracks.
#' @param track_id_prefix Prefix prepended to all auto-generated
#'   `track_id`s. Useful when composing multiple `seq_chip()` calls via
#'   [seq_resolve()].
#' @param legend A `LegendKey` or `SeqLegendSpec` forwarded to each signal
#'   area element. `NULL` (default) produces no legend entry.
#' @param show_legend Logical. When `FALSE`, signal area elements contribute no
#'   legend. Default `TRUE`.
#' @param ... Additional arguments forwarded to [seq_track()].
#' @return A `SeqPlot` with per-sample signal (and optionally peak) tracks.
#' @examples
#' library(GenomicRanges)
#' set.seed(1)
#' make_sig <- function() GRanges("chr1",
#'   IRanges(sort(sample(1:1e6, 200)), width = 500),
#'   score = rexp(200, rate = 0.2))
#' sigs <- list(S1 = make_sig(), S2 = make_sig())
#' seq_chip(sigs, windows = GRanges("chr1", IRanges(1, 1e6)))
#' @export
seq_chip <- function(data,
                     windows,
                     sample_col      = NULL,
                     signal_col      = NULL,
                     peaks           = NULL,
                     peak_col        = NULL,
                     colors          = NULL,
                     scale_max       = NULL,
                     signal_height   = 1,
                     peak_height     = 0.25,
                     show_genes      = NULL,
                     track_id_prefix = "",
                     legend          = NULL,
                     show_legend     = TRUE,
                     ...) {
  .stop_if_not_granges(windows, "windows")

  if (inherits(data, "GRanges")) {
    if (is.null(sample_col))
      stop("When 'data' is a single GRanges, 'sample_col' must name the ",
           "column giving sample identity.", call. = FALSE)
    if (!sample_col %in% names(S4Vectors::mcols(data)))
      stop("sample_col '", sample_col, "' not found in data mcols.",
           call. = FALSE)
    signal_list <- .split_gr_by_sample(data, sample_col)
  } else if (is.list(data)) {
    if (is.null(names(data)) || any(!nzchar(names(data))))
      stop("'data' list must have non-empty names (one per sample).",
           call. = FALSE)
    signal_list <- data
  } else {
    stop("'data' must be a GRanges or a named list of GRanges.",
         call. = FALSE)
  }

  peaks_list <- NULL
  if (!is.null(peaks)) {
    if (inherits(peaks, "GRanges")) {
      if (is.null(sample_col))
        stop("When 'peaks' is a single GRanges, 'sample_col' must name the ",
             "sample column.", call. = FALSE)
      peaks_list <- .split_gr_by_sample(peaks, sample_col)
    } else if (is.list(peaks)) {
      peaks_list <- peaks
    } else {
      stop("'peaks' must be a GRanges, a named list of GRanges, or NULL.",
           call. = FALSE)
    }
  }

  sample_names <- names(signal_list)
  if (is.null(colors) || length(colors) == 0L) {
    pal <- flexoki_palette(max(length(sample_names), 1L))
    track_colors <- stats::setNames(pal[seq_along(sample_names)], sample_names)
  } else {
    track_colors <- colors
    missing_cols <- setdiff(sample_names, names(track_colors))
    if (length(missing_cols) > 0L) {
      pal <- flexoki_palette(length(missing_cols))
      track_colors[missing_cols] <- pal
    }
  }

  p <- seq_plot()

  for (nm in sample_names) {
    sig_gr <- signal_list[[nm]]
    .stop_if_not_granges(sig_gr, paste0("signal for sample '", nm, "'"))
    this_sig_col <- .detect_signal_col(sig_gr, signal_col)
    col          <- unname(track_colors[[nm]])

    y_max <- if (is.null(scale_max)) NULL
             else if (length(scale_max) == 1L) scale_max
             else if (nm %in% names(scale_max)) scale_max[[nm]]
             else NULL

    sig_scale <- if (!is.null(y_max))
                   seq_scale_continuous(limits = c(0, y_max))
                 else NULL

    sig_track <- seq_track(
      data         = sig_gr,
      mapping      = .map_from_names(x = "start", y = this_sig_col),
      windows      = windows,
      track_id     = paste0(track_id_prefix, nm, "_signal"),
      track_height = signal_height,
      scale_y      = sig_scale,
      direction    = "under",
      ...
    )
    p <- p %+% sig_track %+%
      seq_area(aesthetics  = aes(fill = col, color = col, baseline = 0),
               legend      = legend,
               show_legend = show_legend)

    if (!is.null(peaks_list) && nm %in% names(peaks_list)) {
      pk_gr <- peaks_list[[nm]]
      .stop_if_not_granges(pk_gr, paste0("peaks for sample '", nm, "'"))
      pk_mapping <- if (!is.null(peak_col) &&
                        peak_col %in% names(S4Vectors::mcols(pk_gr)))
                      .map_from_names(x = "start", y = peak_col)
                    else
                      .map_from_names(x = "start")
      pk_track <- seq_track(
        data         = pk_gr,
        mapping      = pk_mapping,
        windows      = windows,
        track_id     = paste0(track_id_prefix, nm, "_peaks"),
        track_height = peak_height,
        direction    = "under",
        ...
      )
      p <- p %+% pk_track %+%
        seq_bar(aesthetics = aes(fill = col, color = col))
    }
  }

  if (!is.null(show_genes)) {
    .stop_if_not_granges(show_genes, "show_genes")
    gene_track <- seq_track(
      data         = show_genes,
      mapping      = .map_from_names(x = "start"),
      windows      = windows,
      track_id     = paste0(track_id_prefix, "genes"),
      track_height = 0.5,
      direction    = "under"
    )
    p <- p %+% gene_track %+% seq_gene(mapping = .map_from_names(x = "start"))
  }

  p
}
