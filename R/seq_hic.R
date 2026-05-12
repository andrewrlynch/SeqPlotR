# ── seq_hic — Hi-C contact matrix wrapper ────────────────────────────────────
#
# Single-style Hi-C visualisation. One `seq_hic()` call produces one
# `seq_plot` rendered in the chosen `style`. To combine styles (for
# example a full square beside a triangle), call `seq_hic()` multiple
# times and compose via `seq_resolve()` or `%|%`.

#' Coerce Hi-C input into a sparse contact `GRanges`
#'
#' Accepts either:
#' * a `GRanges` whose mcols already carry `i_start`, `i_end`,
#'   `j_start`, `j_end`, `score`;
#' * a numeric matrix / data.frame whose `rownames` and `colnames` name
#'   bin positions as `"chr:start-end"` or integer bin indices.
#'
#' @param data Input as described above.
#' @param windows A `GRanges` giving the genomic view (used to resolve
#'   integer-indexed matrix rownames / colnames).
#' @return A `GRanges` with `i_start`, `i_end`, `j_start`, `j_end`,
#'   `score` mcols.
#' @keywords internal
.hic_to_granges <- function(data, windows) {
  if (is.matrix(data) || is.data.frame(data)) {
    m   <- as.matrix(data)
    rn  <- rownames(m) %||% as.character(seq_len(nrow(m)))
    cn  <- colnames(m) %||% as.character(seq_len(ncol(m)))
    parse_bin <- function(labels) {
      if (all(grepl("^[A-Za-z0-9_]+:[0-9]+-[0-9]+$", labels))) {
        parsed <- utils::strcapture(
          pattern = "^([^:]+):([0-9]+)-([0-9]+)$",
          x       = labels,
          proto   = list(chrom = character(),
                         start = integer(),
                         end   = integer())
        )
        list(chrom = parsed$chrom,
             start = parsed$start,
             end   = parsed$end)
      } else {
        # Integer bin indices — treat the first window in `windows` as
        # the reference and place bins at start + (index-1)*step.
        idx   <- as.integer(labels)
        nbin  <- length(idx)
        w     <- windows[1]
        step  <- max(1L, as.integer(BiocGenerics::width(w) / nbin))
        chrom <- rep(as.character(GenomicRanges::seqnames(w)), nbin)
        list(
          chrom = chrom,
          start = BiocGenerics::start(w) + (idx - 1L) * step,
          end   = BiocGenerics::start(w) + idx * step - 1L
        )
      }
    }
    rb <- parse_bin(rn)
    cb <- parse_bin(cn)
    idx_r <- rep(seq_len(nrow(m)), times = ncol(m))
    idx_c <- rep(seq_len(ncol(m)), each  = nrow(m))
    scr   <- as.numeric(m)
    keep  <- is.finite(scr) & scr != 0
    if (!any(keep))
      return(GenomicRanges::GRanges(
        seqnames = character(0),
        ranges   = IRanges::IRanges(integer(0), integer(0))))
    idx_r <- idx_r[keep]; idx_c <- idx_c[keep]; scr <- scr[keep]
    return(GenomicRanges::GRanges(
      seqnames = rb$chrom[idx_r],
      ranges   = IRanges::IRanges(start = rb$start[idx_r],
                                  end   = rb$end[idx_r]),
      i_start  = rb$start[idx_r], i_end = rb$end[idx_r],
      j_start  = cb$start[idx_c], j_end = cb$end[idx_c],
      score    = scr
    ))
  }

  .stop_if_not_granges(data, "data")
  required <- c("i_start", "i_end", "j_start", "j_end", "score")
  missing  <- setdiff(required, names(S4Vectors::mcols(data)))
  if (length(missing) > 0L)
    stop("Hi-C data missing required columns: ",
         paste(missing, collapse = ", "),
         "\nExpected: i_start, i_end, j_start, j_end, score",
         call. = FALSE)
  data
}

#' Hi-C contact matrix
#'
#' Builds a [seq_plot()] rendering a Hi-C contact matrix in one of four
#' styles:
#' \describe{
#'   \item{`"full"`}{Symmetric square heatmap with genomic position on
#'     both x- and y-axes.}
#'   \item{`"diagonal"`}{Same coordinate system as `"full"` (kept as a
#'     separate keyword so call sites can switch styles without changing
#'     shape).}
#'   \item{`"triangle"`}{Rotated 45 degrees; upper triangle only, y-axis
#'     is interaction distance in base pairs.}
#'   \item{`"rectangle"`}{Same rotation as `"triangle"`, but y-axis is
#'     capped at `max_dist`, yielding a rectangle.}
#' }
#'
#' To show multiple styles side-by-side, call `seq_hic()` multiple times
#' and combine via [seq_resolve()] or `%|%`.
#'
#' @param data Sparse `GRanges` (mcols: `i_start`, `i_end`, `j_start`,
#'   `j_end`, `score`) or a numeric matrix / data.frame whose row/column
#'   names encode bin positions. Each row must describe a well-formed
#'   contact where `i_end - i_start` and `j_end - j_start` both equal
#'   the Hi-C bin width (they define the tile's footprint on the
#'   position and distance axes). For cross-chromosomal contacts,
#'   include an optional `j_chrom` mcols column giving the j-bin's
#'   chromosome (the GRanges's `seqnames` is taken to be the i-bin's
#'   chromosome). When absent, both bins are assumed to live on the
#'   same chromosome.
#' @param windows `GRanges` defining the genomic region(s) to display
#'   on the x-axis. Multiple ranges produce side-by-side panels (one
#'   per range), useful for comparing several regions.
#' @param style One of `"full"`, `"diagonal"`, `"triangle"`,
#'   `"rectangle"`. Default `"triangle"`.
#' @param max_dist For `style = "rectangle"` only: cap the distance axis
#'   at this value (bp). Required for `"rectangle"`.
#' @param palette Colour scale palette for the tile fill. Passed to
#'   [seq_scale_color_continuous()].
#' @param na_color Colour for zero/NA contacts.
#' @param y_windows Optional `GRanges` for the genomic y-axis range in
#'   `"full"` / `"diagonal"` styles. Defaults to `windows` (square
#'   matrix). Pass a different `GRanges` to display an asymmetric
#'   region pair. Ignored for rotated styles.
#' @param combine_windows Logical; when `TRUE`, multi-region `windows`
#'   are concatenated into a single virtual track so cross-window data
#'   (e.g. inter-chromosomal contacts) renders continuously in one
#'   panel. Default `FALSE`.
#' @param combine_y_windows Symmetric to `combine_windows` for
#'   multi-region `y_windows` in the `full` / `diagonal` styles.
#' @param flip_x,flip_y Logical. Mirror the x or y axis. For the
#'   `triangle` style `flip_y = TRUE` produces a downward-pointing
#'   triangle; for `diagonal` it switches to the lower diagonal; for
#'   `full` it flips the matrix vertically (y) or horizontally (x).
#'   Tick labels follow the same orientation. Default `FALSE`.
#' @param track_height Relative track height.
#' @param track_id `track_id` for the generated track. Defaults to
#'   `paste0("hic_", style)`.
#' @param legend A `LegendKey` or `SeqLegendSpec` forwarded to the tile
#'   element. `NULL` (default) produces no legend entry.
#' @param show_legend Logical. When `FALSE`, the tile element contributes no
#'   legend. Default `TRUE`.
#' @param ... Additional arguments forwarded to [seq_track()].
#' @return A `SeqPlot` with a single Hi-C track.
#' @examples
#' library(GenomicRanges)
#' set.seed(1)
#' n  <- 80
#' st <- sort(sample(seq(1, 1e6, by = 1e4), n))
#' gr <- GRanges("chr1", IRanges(st, width = 1e4),
#'               i_start = st, i_end = st + 1e4,
#'               j_start = st + sample(0:5e5, n, replace = TRUE),
#'               j_end   = st + sample(0:5e5, n, replace = TRUE) + 1e4,
#'               score   = rexp(n, rate = 0.5))
#' win <- GRanges("chr1", IRanges(1, 1e6))
#' seq_hic(gr, windows = win, style = "triangle")
#' @export
seq_hic <- function(data,
                    windows,
                    style             = "triangle",
                    max_dist          = NULL,
                    palette           = "blues",
                    na_color          = "#FFFFD9",
                    y_windows         = NULL,
                    combine_windows   = FALSE,
                    combine_y_windows = FALSE,
                    flip_x            = FALSE,
                    flip_y            = FALSE,
                    track_height      = 1,
                    track_id          = NULL,
                    legend            = NULL,
                    show_legend       = TRUE,
                    ...) {
  style <- match.arg(style, c("full", "diagonal", "triangle", "rectangle"))
  if (style == "rectangle" && is.null(max_dist))
    stop("max_dist is required for style = 'rectangle'.", call. = FALSE)
  .stop_if_not_granges(windows, "windows")
  if (!is.null(y_windows)) .stop_if_not_granges(y_windows, "y_windows")

  sparse_gr <- .hic_to_granges(data, windows)
  if (length(sparse_gr) == 0L) {
    # Build a non-empty placeholder so the plot can still be laid out.
    sparse_gr <- GenomicRanges::GRanges(
      seqnames = as.character(GenomicRanges::seqnames(windows))[1],
      ranges   = IRanges::IRanges(BiocGenerics::start(windows)[1], width = 1L),
      i_start  = BiocGenerics::start(windows)[1],
      i_end    = BiocGenerics::start(windows)[1],
      j_start  = BiocGenerics::start(windows)[1],
      j_end    = BiocGenerics::start(windows)[1],
      score    = 0
    )
  }

  S4Vectors::mcols(sparse_gr)$.log_score <- log1p(S4Vectors::mcols(sparse_gr)$score)

  rotate <- style %in% c("triangle", "rectangle")

  # Canonicalise i/j-bin ranges so start <= end (Hi-C inputs sometimes
  # encode upper-triangle-only data with j < i).
  mc_raw  <- S4Vectors::mcols(sparse_gr)
  is_lo   <- pmin(mc_raw$i_start, mc_raw$i_end)
  is_hi   <- pmax(mc_raw$i_start, mc_raw$i_end)
  js_lo   <- pmin(mc_raw$j_start, mc_raw$j_end)
  js_hi   <- pmax(mc_raw$j_start, mc_raw$j_end)

  # Per-bin chromosome assignment. The GRanges's seqnames is the i-bin
  # chromosome; an optional `j_chrom` mcols column gives the j-bin's
  # chromosome (defaults to the i-bin chromosome when absent).
  i_chrom <- as.character(GenomicRanges::seqnames(sparse_gr))
  j_chrom <- if ("j_chrom" %in% names(mc_raw))
               as.character(mc_raw$j_chrom)
             else
               i_chrom
  same_chrom <- i_chrom == j_chrom

  # For triangle / rectangle / diagonal, keep upper-triangle tiles
  # only (j >= i) for *intra-chromosomal* contacts. Inter-chromosomal
  # contacts have no diagonal symmetry to deduplicate, so keep them all.
  if (style %in% c("diagonal", "triangle", "rectangle")) {
    keep <- (!same_chrom) | (js_lo >= is_lo)
    if (!any(keep)) keep[1] <- TRUE
    sparse_gr <- sparse_gr[keep]
    mc_raw    <- mc_raw[keep, , drop = FALSE]
    is_lo <- is_lo[keep]; is_hi <- is_hi[keep]
    js_lo <- js_lo[keep]; js_hi <- js_hi[keep]
    i_chrom <- i_chrom[keep]; j_chrom <- j_chrom[keep]
  }

  # `seq_tile` expects per-tile x/y genomic ranges. Build data2 (y-axis
  # genomic coordinates) from the j-bin columns, carrying the j-bin's
  # chromosome so the layout can route inter-chrom tiles correctly.
  tile_data2 <- GenomicRanges::GRanges(
    seqnames = j_chrom,
    ranges   = IRanges::IRanges(start = js_lo, end = js_hi)
  )
  # Swap the primary GRanges ranges over the i-bin columns so that
  # `start(gr)`/`end(gr)` drive the x extent consistently regardless of
  # the caller's original ranges.
  old_mc    <- mc_raw
  sparse_gr <- GenomicRanges::GRanges(
    seqnames = i_chrom,
    ranges   = IRanges::IRanges(start = is_lo, end = is_hi)
  )
  S4Vectors::mcols(sparse_gr) <- old_mc

  # Genomic y axis only applies to full/diagonal. Default to mirroring
  # the x windows (square matrix); a user-supplied `y_windows` lets the
  # x and y ranges differ (e.g., one chromosome on x, another on y).
  y_windows <- if (!rotate) (y_windows %||% windows) else NULL

  # Virtualize j positions in tile_data2 when the y axis is combined:
  #   - rotated styles (triangle/rectangle): combine_windows virtualizes
  #     against the x-window map (the y-axis is distance, derived from
  #     virtualized i and j on the same combined axis).
  #   - full / diagonal styles: combine_y_windows virtualizes against
  #     the y_windows map (which may differ from the x map). For
  #     symmetric square matrices y_windows == windows so the two maps
  #     are identical.
  vmap_for_j <- NULL
  if (rotate && isTRUE(combine_windows) && length(windows) > 1L) {
    vmap_for_j <- .build_virtual_map(windows, virtual_gap = 0)
  } else if (!rotate && isTRUE(combine_y_windows) &&
             !is.null(y_windows) && length(y_windows) > 1L) {
    vmap_for_j <- .build_virtual_map(y_windows, virtual_gap = 0)
  }
  if (!is.null(vmap_for_j)) {
    j_seq    <- as.character(GenomicRanges::seqnames(tile_data2))
    j_vstart <- .virtualize_positions(j_seq,
                                       BiocGenerics::start(tile_data2),
                                       vmap_for_j)
    j_vend   <- .virtualize_positions(j_seq,
                                       BiocGenerics::end(tile_data2),
                                       vmap_for_j)
    keep_j   <- !is.na(j_vstart) & !is.na(j_vend)
    if (any(!keep_j)) {
      sparse_gr <- sparse_gr[keep_j]
      tile_data2 <- tile_data2[keep_j]
      j_vstart  <- j_vstart[keep_j]
      j_vend    <- j_vend[keep_j]
      js_lo     <- js_lo[keep_j]; js_hi <- js_hi[keep_j]
      is_lo     <- is_lo[keep_j]; is_hi <- is_hi[keep_j]
      i_chrom   <- i_chrom[keep_j]; j_chrom <- j_chrom[keep_j]
    }
    tile_data2 <- GenomicRanges::GRanges(
      "__combined__",
      IRanges::IRanges(start = j_vstart, end = j_vend)
    )
  }

  # For triangle/rectangle, the rotated y-axis represents interaction
  # distance (full bp). Derive the data-driven cap from the largest
  # virtualized j-i span when no explicit `max_dist` is given.
  if (rotate) {
    if (isTRUE(combine_windows) && length(windows) > 1L) {
      vmap_x_local <- .build_virtual_map(windows, virtual_gap = 0)
      i_vstart <- .virtualize_positions(i_chrom, is_lo, vmap_x_local)
      i_vend   <- .virtualize_positions(i_chrom, is_hi, vmap_x_local)
      j_vstart_full <- .virtualize_positions(j_chrom, js_lo, vmap_x_local)
      j_vend_full   <- .virtualize_positions(j_chrom, js_hi, vmap_x_local)
      max_full_dist <- max(c(j_vend_full - i_vstart,
                             i_vend - j_vstart_full), 0, na.rm = TRUE)
    } else {
      max_full_dist <- max(js_hi - is_lo, 0, na.rm = TRUE)
    }
  }
  scale_y <- if (style == "rectangle")
               seq_scale_continuous(limits = c(0, max_dist))
             else if (rotate)
               seq_scale_continuous(limits = c(0, max_full_dist))
             else
               NULL

  fill_scale <- seq_scale_fill_continuous(palette = palette,
                                          na_value = na_color)
  # The tile fill is a per-tile colour derived from log_score via a
  # continuous colour ramp. We resolve the palette here (so the wrapper
  # is self-contained) and attach the resulting colours to each tile
  # via `.fill_color` BEFORE the GRanges goes into the track — R's
  # copy-on-modify semantics mean any later mcols<- creates a new
  # GRanges that the track wouldn't see.
  logs <- S4Vectors::mcols(sparse_gr)$.log_score
  lims <- range(logs, na.rm = TRUE)
  if (!all(is.finite(lims)) || diff(lims) == 0) lims <- c(0, 1)
  t_vals <- pmax(0, pmin(1, (logs - lims[1]) / diff(lims)))
  stops  <- switch(palette,
    viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
    plasma  = c("#0d0887", "#cc4778", "#f0f921"),
    magma   = c("#000004", "#b63679", "#fcfdbf"),
    blues   = c("#f7fbff", "#2171b5", "#08306b"),
    reds    = c("#fff5f0", "#ef3b2c", "#67000d"),
    c("#f7fbff", "#2171b5", "#08306b"))
  ramp  <- grDevices::colorRamp(stops)
  cols  <- ramp(t_vals)
  hex   <- grDevices::rgb(cols[, 1], cols[, 2], cols[, 3],
                          maxColorValue = 255)
  hex[!is.finite(logs) | logs == 0] <- na_color
  S4Vectors::mcols(sparse_gr)$.fill_color <- hex

  trk <- seq_track(
    data              = sparse_gr,
    mapping           = .map_from_names(x = "start"),
    windows           = windows,
    y_windows         = y_windows,
    track_height      = track_height,
    track_id          = track_id %||% paste0("hic_", style),
    scale_y           = scale_y,
    combine_windows   = combine_windows,
    combine_y_windows = combine_y_windows,
    flip_x            = flip_x,
    flip_y            = flip_y,
    ...
  )

  # Don't set the tile element's `data` directly — let it inherit from
  # the track via `track_data`. This matters when `combine_windows` is
  # TRUE: the track's data is virtualized at draw time, so inheriting
  # gives the tile the virtualized GRanges. Setting element-level data
  # would shadow the virtualized version with the original genomic
  # GRanges and produce an empty plot.
  tile_elem <- seq_tile(
    mapping     = .map_from_names(x = "start", fill = ".fill_color"),
    data2       = tile_data2,
    aesthetics  = aes(rotate = rotate, na_color = na_color),
    legend      = legend,
    show_legend = show_legend
  )

  p <- seq_plot() %+% trk %+% tile_elem
  p
}
