# map() column names evaluated against element data ----
utils::globalVariables(c("start", "end", "cn", "x0", "x1", "x0_end", "x1_end", "ref_chrom", "ref_start", "ref_end",
                         "flipped", "label", "fill_key", "tint", "dashed", "centromere_start", "centromere_end",
                         "chrom0", "chrom1", "side0", "side1", "kind", "bulge", "twisted"))

#' Plot a derivative-chromosome reconstruction
#'
#' Draws a `plot_derivatives` table set (schema version 1), as written by the
#' ROSE package, with [seq_block()], [seq_connector()], [seq_synteny()],
#' [seq_track_label()] and [seq_caption()]. The adapter only selects, renames
#' and maps table columns to element arguments. One track is drawn per row:
#' a reference row, derivative rows, deleted copies and candidate strips. All
#' rows share one base-pair scale. Clusters with `cn_panels` also get copy
#' number tracks above the rows.
#'
#' @param x A directory holding `plot_derivatives.*.tsv` files, or a list of
#'   data frames with a `schema_version`. See [seq_read_table_set()].
#' @param fill `"cytoband"` (default) or `"solid"`.
#' @param palette `NULL` (default hues), `"soft"`, or a character vector of
#'   colours, optionally named by `fill_key`.
#' @param block_colors Optional named colours overriding block fills by label,
#'   then by key.
#' @param synteny Logical. Draw ribbons between the reference and derivative
#'   rows.
#' @param label_inside Logical. Place block labels inside blocks when they fit.
#' @param flip_style How flipped blocks are marked: `"rotate"`, `"overline"`
#'   or `"minus"`.
#' @param show_loose Logical. Draw the loose derivative rows and the rows that
#'   continue them.
#' @param cytobands Cytoband table passed to [seq_block()]. `NULL` uses the
#'   bundled hg38 bands.
#' @param annotations Optional data frame for [seq_block()] annotations, with
#'   an optional `track` column naming the row label to draw on. Rows without
#'   `track` go on the reference and derivative rows.
#' @param ... Arguments merged over the `seq_plot()` defaults.
#' @return A `SeqPlot`, extendable with `%+%`.
#' @export
seq_rose_derivatives <- function(x, fill = c("cytoband", "solid"), palette = NULL,
                                 block_colors = NULL, synteny = FALSE, label_inside = FALSE,
                                 flip_style = c("rotate", "overline", "minus"), show_loose = TRUE,
                                 cytobands = NULL, annotations = NULL, ...) {
  fill <- match.arg(fill)
  flip_style <- match.arg(flip_style)
  req <- list(
    rows = c("row_id", "kind", "y_order", "label", "sublabel", "dashed", "candidate_id", "derivative_id"),
    blocks = c("row_id", "block_id", "fragment_id", "label", "chrom", "start", "end", "x0", "x1",
               "flipped", "fill_key", "tint", "centromere_start", "centromere_end", "label_inside"),
    connectors = c("row_from", "block_from", "side_from", "row_to", "block_to", "side_to", "kind", "same_row_up"),
    ribbons = c("fragment_id", "row_id", "block_id", "twisted"),
    windows = c("row_id", "chrom", "start", "end"),
    caption = c("text", "approximate", "approximate_reason")
  )
  s <- seq_read_table_set(x, "plot_derivatives", supported = 1L, required = req)
  rows <- s$rows[order(s$rows$y_order), , drop = FALSE]
  blocks <- s$blocks; conns <- s$connectors; ribs <- s$ribbons; wins <- s$windows

  # drop loose rows and their continuation rows ----
  if (!isTRUE(show_loose)) {
    key <- paste(rows$kind, rows$candidate_id, rows$derivative_id)
    loose <- !is.na(rows$sublabel) & rows$sublabel == "loose"
    drop <- rows$row_id[key %in% key[loose]]
    rows <- rows[!rows$row_id %in% drop, , drop = FALSE]
    blocks <- blocks[!blocks$row_id %in% drop, , drop = FALSE]
    wins <- wins[!wins$row_id %in% drop, , drop = FALSE]
    conns <- conns[!conns$row_from %in% drop & !conns$row_to %in% drop, , drop = FALSE]
    ribs <- ribs[!ribs$row_id %in% drop, , drop = FALSE]
  }

  # palette over fill keys in data order ----
  keys <- unique(blocks$fill_key)
  pal <- .block_palette(palette)
  pal_named <- if (!is.null(names(pal))) pal else stats::setNames(pal[(seq_along(keys) - 1L) %% length(pal) + 1L], keys)
  tints <- c(0.15, 0.5, 0.3, 0.65, 0.4)

  # window name of each block ----
  blocks$win <- vapply(seq_len(nrow(blocks)), function(i) {
    w <- wins[wins$row_id == blocks$row_id[i] & wins$start <= blocks$x0[i] & blocks$x1[i] <= wins$end + 1, , drop = FALSE]
    if (nrow(w) == 0L) w <- wins[wins$row_id == blocks$row_id[i], , drop = FALSE]
    w$chrom[1]
  }, character(1))

  ref_row <- rows$row_id[rows$kind == "reference"][1]
  p <- do.call(seq_plot, utils::modifyList(list(plot_margin = c(0.02, 0.02, 0.02, 0.02)), list(...)))
  left_right <- list(left = 0.13, right = 0.02)

  # copy number tracks of cluster plots ----
  if (!is.null(s$cn_panels) && nrow(s$cn_panels) > 0L && !is.null(s$cn_segments)) {
    rb <- blocks[blocks$row_id == ref_row, , drop = FALSE]
    rb$col <- .block_fill_colors(rb$fill_key, rb$tint, rb$label, pal_named, tints, block_colors)
    pn <- s$cn_panels[order(s$cn_panels$y_order), , drop = FALSE]
    for (k in seq_len(nrow(pn))) {
      sg <- s$cn_segments[s$cn_segments$panel_id == pn$panel_id[k], , drop = FALSE]
      sg$col <- ifelse(!is.na(sg$flank) & sg$flank, "grey75", rb$col[match(sg$colour_key, rb$label)])
      sg$col[is.na(sg$col)] <- "grey75"
      p <- p %+% seq_track(
        track_id = paste0("cn", pn$panel_id[k]), direction = "under",
        windows = GenomicRanges::GRanges(pn$chrom[k], IRanges::IRanges(pn$win_start[k], pn$win_end[k])),
        track_height = 2, track_outer_margin = left_right
      )
      if (nrow(sg) > 0L) {
        p <- p %+% seq_segment(
          data = GenomicRanges::GRanges(sg$chrom, IRanges::IRanges(sg$start, sg$end), cn = sg$cn, col = sg$col),
          mapping = map(x = start, x_end = end, y = cn, y_end = cn, color = col),
          aesthetics = aes(linewidth = 3.2)
        )
      }
      p <- p %+% seq_track_label(label = pn$label[k], sublabel = "Copy number")
    }
  }

  # one track per row ----
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    w <- wins[wins$row_id == r$row_id, , drop = FALSE]
    b <- blocks[blocks$row_id == r$row_id, , drop = FALSE]
    multi_ref <- r$kind == "reference" && nrow(w) > 1L
    outer <- list(top = 0, bottom = if (multi_ref) 0.03 else 0, left = 0.13, right = 0.02)
    A <- if (multi_ref) aes(axis.x.line = NA, axis.x.ticks = NA, axis.x.text = NA) else aes(axis.x.visible = FALSE)
    p <- p %+% seq_track(
      track_id = paste0("row", r$row_id), direction = "under",
      windows = GenomicRanges::GRanges(w$chrom, IRanges::IRanges(w$start, w$end)),
      x_span = "shared", track_height = 1, track_outer_margin = outer,
      track_inner_margin = 0, window_outer_margin = 0, window_inner_margin = 0,
      aesthetics = A
    )
    if (nrow(b) > 0L) {
      gr <- GenomicRanges::GRanges(
        b$win, IRanges::IRanges(b$x0, b$x1), block_id = b$block_id, fragment_id = b$fragment_id,
        label = b$label, ref_chrom = b$chrom, ref_start = b$start, ref_end = b$end,
        x0 = b$x0, x1 = b$x1, flipped = b$flipped, fill_key = b$fill_key, tint = b$tint,
        label_inside = isTRUE(label_inside) & b$label_inside, dashed = rep(isTRUE(r$dashed), nrow(b)),
        centromere_start = b$centromere_start, centromere_end = b$centromere_end
      )
      strip <- r$kind == "candidate_strip"
      ann <- annotations
      if (!is.null(ann)) {
        if ("track" %in% names(ann)) {
          keep <- if (r$kind %in% c("reference", "derivative")) is.na(ann$track) | ann$track == r$label else !is.na(ann$track) & ann$track == r$label
          ann <- ann[keep, , drop = FALSE]
        } else if (!r$kind %in% c("reference", "derivative")) {
          ann <- NULL
        }
        if (!is.null(ann) && nrow(ann) == 0L) ann <- NULL
      }
      p <- p %+% seq_block(
        data = gr,
        mapping = map(x = x0, x_end = x1, chrom = ref_chrom, ref_start = ref_start, ref_end = ref_end,
                      flipped = flipped, label = label, key = fill_key, tint = tint,
                      label_inside = label_inside, dashed = dashed,
                      centromere_start = centromere_start, centromere_end = centromere_end),
        fill = if (strip) "solid" else fill, shape = if (strip) "arrow" else "rect",
        palette = pal_named, block_colors = block_colors, flip_style = flip_style,
        cytobands = cytobands, label_side = if (r$kind == "reference") "above" else "below",
        aesthetics = if (r$kind == "deleted") aes(fade = 0.55) else aes(),
        annotations = ann, show_legend = FALSE
      )
    }
    p <- p %+% seq_track_label(label = r$label, sublabel = r$sublabel)
  }

  # junction connectors by row pair ----
  if (nrow(conns) > 0L) {
    f <- blocks[match(conns$block_from, blocks$block_id), ]
    g <- blocks[match(conns$block_to, blocks$block_id), ]
    df <- data.frame(chrom0 = f$win, x0 = f$x0, x0_end = f$x1, side0 = conns$side_from,
                     chrom1 = g$win, x1 = g$x0, x1_end = g$x1, side1 = conns$side_to,
                     kind = conns$kind, bulge = conns$same_row_up, stringsAsFactors = FALSE)
    grp <- paste(conns$row_from, conns$row_to)
    for (gk in unique(grp)) {
      ix <- which(grp == gk)
      p <- p %+% seq_connector(
        data = df[ix, , drop = FALSE],
        mapping = map(x0 = x0, x0_end = x0_end, chrom0 = chrom0, side0 = side0, x1 = x1,
                      x1_end = x1_end, chrom1 = chrom1, side1 = side1, style = kind, bulge = bulge),
        t0 = paste0("row", conns$row_from[ix[1]]), t1 = paste0("row", conns$row_to[ix[1]]),
        show_legend = FALSE
      )
    }
  }

  # ribbons from the reference row to each derivative row ----
  if (isTRUE(synteny) && nrow(ribs) > 0L) {
    rb <- blocks[blocks$row_id == ref_row, , drop = FALSE]
    r <- rb[match(ribs$fragment_id, rb$fragment_id), ]
    d <- blocks[match(ribs$block_id, blocks$block_id), ]
    df <- data.frame(chrom0 = r$win, x0 = r$x0, x0_end = r$x1, chrom1 = d$win, x1 = d$x0, x1_end = d$x1,
                     fill = .block_fill_colors(d$fill_key, d$tint, d$label, pal_named, tints, block_colors),
                     twisted = ribs$twisted, stringsAsFactors = FALSE)
    alpha <- if (length(unique(ribs$ref_chrom)) > 1L) 0.45 else 0.6
    for (rid in unique(ribs$row_id)) {
      ix <- which(ribs$row_id == rid & !is.na(df$chrom0))
      if (length(ix) == 0L) next
      p <- p %+% seq_synteny(
        data = df[ix, , drop = FALSE],
        mapping = map(x0 = x0, x0_end = x0_end, x1 = x1, x1_end = x1_end, chrom0 = chrom0,
                      chrom1 = chrom1, fill = fill, twisted = twisted),
        t0 = paste0("row", ref_row), t1 = paste0("row", rid), height = 0.5, smooth = TRUE,
        layer = "below", aesthetics = aes(alpha = alpha, color = NA)
      )
    }
  }

  # caption ----
  cap <- s$caption
  if (!is.null(cap) && nrow(cap) > 0L) {
    banner <- if (isTRUE(cap$approximate[1])) {
      if (is.na(cap$approximate_reason[1])) "APPROXIMATE" else cap$approximate_reason[1]
    }
    p <- p %+% seq_caption(text = cap$text[1], banner = banner)
  }
  p
}
