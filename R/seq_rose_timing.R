# ROSE timing adapter ----

# default palettes, merged with the user palette element-wise
.timing_palette <- function() {
  list(
    sv = c(DEL = "#E69F00", DUP = "#009E73", h2hINV = "#7B3294", t2tINV = "#C2A5CF",
           TRA = "#8C564B", SGL = "#7F7F7F", untimed = "grey70"),
    step = c(gain = "#A50F15", loss = "#08519C", wgd = "#000000"),
    cn = c(major = "#D62728", minor = "#1F77B4"),
    sbs = c("C>A" = "#1EBFF0", "C>G" = "#050708", "C>T" = "#E62725",
            "T>A" = "#CBCACB", "T>C" = "#A1CE63", "T>G" = "#EDB6C2"),
    heatmap = c(gain = "#D62728", loss = "#2B59C3"),
    glyph = c(snv = "#1C1B1A", proxy = "grey55"),
    annotation = NULL,
    node = c(SV = "#7F7F7F", GAIN = "#A50F15", LOSS = "#08519C", WGD = "#000000"),
    edge = list(snv = list(color = "grey10", linetype = 1, linewidth = 1.4),
                mixed = list(color = "grey45", linetype = 1, linewidth = 1.2),
                structure = list(color = "grey65", linetype = 2, linewidth = 1.0),
                other = list(color = "grey60", linetype = 2, linewidth = 0.9))
  )
}

# required columns of each table in the plot_timing set
.timing_required <- list(
  windows = c("chr", "start", "end", "order"),
  cn_bands = c("chr", "start", "end", "cn_major", "cn_minor"),
  snv_points = c("chr", "pos", "multiplicity", "sbs"),
  sv_arches = c("arch_id", "chr1", "pos1", "chr2", "pos2", "sv_type", "multiplicity",
                "glyph", "y", "bar_lo", "bar_hi", "layer", "faded", "alpha",
                "ecdna_caution", "colour_key"),
  cn_steps = c("chr", "start", "end", "direction", "k", "label", "credit", "glyph",
               "y", "bar_lo", "bar_hi", "faded", "alpha", "colour_key"),
  glyphs = c("glyph_id", "kind", "chr", "start", "end", "x", "glyph", "y", "bar_lo",
             "bar_hi", "faded", "alpha", "label", "ecdna_caution", "colour_key"),
  heatmap = c("chr", "start", "end", "tau_lo", "tau_hi", "dens_gain", "dens_loss",
              "intensity_gain", "intensity_loss"),
  heatmap_bounds = c("chr", "start", "end", "direction", "y"),
  poset_nodes = c("event_id", "shape_key", "x", "y", "bar_lo", "bar_hi", "alpha", "label"),
  poset_edges = c("from_event", "to_event", "kind", "arrow", "style_key", "label", "conflict"),
  poset_groups = c("group_id", "x", "y", "label"),
  annotations = c("chr", "start", "end", "label", "colour_key"),
  settings = c("smooth_sd", "d_max", "cn_log_axis")
)

#' Draw a ROSE timing plot table set
#'
#' Thin adapter that turns a `plot_timing` table set written by ROSE into a
#' [seq_plot()] built from the generic timing elements ([seq_cross_arch()],
#' [seq_step_ribbon()], [seq_paired_ribbon()], [seq_interval_glyph()],
#' [seq_tile_mix()], [seq_node_link()], [seq_span_strip()]). The tables carry
#' every position, interval, opacity, label and colour category, so the adapter
#' only selects rows, renames columns, looks up palette colours and sets axes.
#' Proxy rows are order bounds, not times. Glyph rows without a position, such as
#' a whole-genome doubling, are drawn across every window at their given `y`.
#'
#' @param x A directory holding the `plot_timing.*` files, or a list of data
#'   frames with a `schema_version` element or attribute. Read with
#'   [seq_read_table_set()].
#' @param layout One of `"arches"` (SV arches over copy number with step
#'   ribbons), `"order"` (arches at their poset layer), `"heatmap"` (arches
#'   over a two-channel timing heatmap) or `"poset"` (node-link view).
#' @param poset_view For `layout = "poset"`: `"poset"` (all edges), `"hasse"`
#'   (reduction edges only) or `"groups"` (group boxes only).
#' @param palette Named list merged over the default palettes. Entries are
#'   `sv`, `step`, `cn`, `sbs`, `heatmap`, `glyph` and `node` (named colour
#'   vectors), `annotation` (named colours by `colour_key`) and `edge` (named
#'   lists of `color`, `linetype` and `linewidth`).
#' @param show_hidden Logical. Add a strip for glyph rows marked `hidden`.
#' @param plain Logical. Draw a genomic x axis on the last timing track instead
#'   of an ideogram.
#' @param cytobands Cytoband `GRanges` for the ideogram. Defaults to
#'   [load_cytobands()].
#' @param legend_width Right margin of the page, as a fraction, left for the
#'   legends.
#' @param ... Passed to [seq_plot()].
#' @return A `SeqPlot` that can be extended with `%+%`.
#' @examples
#' d <- system.file("extdata", "rose", "timing", "sim_bfb_amplicon", package = "SeqPlotR")
#' if (nzchar(d)) {
#'   p <- seq_rose_timing(d, layout = "arches", plain = TRUE)
#' }
#' @export
seq_rose_timing <- function(x, layout = c("arches", "order", "heatmap", "poset"),
                            poset_view = c("poset", "hasse", "groups"),
                            palette = NULL, show_hidden = FALSE, plain = FALSE,
                            cytobands = NULL, legend_width = 0.2, ...) {
  layout <- match.arg(layout)
  poset_view <- match.arg(poset_view)
  tb <- seq_read_table_set(x, "plot_timing", supported = 1L, required = .timing_required)
  tb <- lapply(tb[names(.timing_required)], as.data.frame, stringsAsFactors = FALSE)
  pal <- .timing_palette()
  # named vectors merge element-wise over the defaults, other entries replace them
  for (k in names(palette)) {
    u <- palette[[k]]
    pal[[k]] <- if (is.atomic(u) && !is.null(names(u)) && is.atomic(pal[[k]]) && !is.null(names(pal[[k]]))) {
      replace(pal[[k]], names(u), u)
    } else if (is.list(u) && is.list(pal[[k]])) utils::modifyList(pal[[k]], u) else u
  }
  st <- tb$settings[1, , drop = FALSE]
  na0 <- function(v) v[!is.na(v)]

  # shared helpers ----
  pick <- function(vec, key, fallback) {
    out <- unname(vec[as.character(key)])
    out[is.na(out)] <- fallback
    out
  }
  to_gr <- function(df, chr = "chr", s = "start", e = "end") {
    m <- df[setdiff(names(df), c("chr", "start", "end", "width", "strand", "seqnames"))]
    gr <- GenomicRanges::GRanges(df[[chr]], IRanges::IRanges(df[[s]], df[[e]]))
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(m)
    gr
  }
  leg <- function(keys) {
    if (!length(keys)) return(NULL)
    seq_legend(unname(keys), position = "canvas_margin", side = "right",
               orientation = "vertical", row_height = 1.1, x = 0.05, y = 0.95)
  }
  k_rect <- function(label, col, alpha = 1) LegendKey(label, color = col, fill = col, shape = "rect", alpha = alpha)
  k_line <- function(label, col, lwd = 2, lty = 1) LegendKey(label, color = col, shape = "-", lwd = lwd, lty = lty)
  k_pt <- function(label, col) LegendKey(label, color = col, fill = col, shape = "circle")
  # column-name mapping built without non-standard evaluation
  mp <- function(...) structure(lapply(list(...), as.name), class = "SeqMap")

  # windows ----
  w <- tb$windows[order(tb$windows$order), , drop = FALSE]
  wins <- GenomicRanges::GRanges(w$chr, IRanges::IRanges(w$start, w$end))
  tau_scale <- seq_scale_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

  # poset layout ----
  if (layout == "poset") {
    nd <- tb$poset_nodes
    ed <- tb$poset_edges
    gp <- tb$poset_groups
    if (poset_view == "hasse") ed <- ed[ed$kind == "hasse", , drop = FALSE]
    if (poset_view == "groups") { nd <- nd[0, , drop = FALSE]; ed <- ed[0, , drop = FALSE] }
    else gp <- gp[0, , drop = FALSE]
    shp <- c(SV = "circle", GAIN = "triangle_up", LOSS = "triangle_down", WGD = "diamond")
    nodes <- if (nrow(nd)) data.frame(
      id = nd$event_id, x = nd$x, y = if ("anchor_y" %in% names(nd)) ifelse(is.na(nd$anchor_y), nd$y, nd$anchor_y) else nd$y, bar_lo = nd$bar_lo, bar_hi = nd$bar_hi,
      shape = pick(shp, nd$shape_key, "circle"), fill = pick(pal$node, nd$shape_key, "grey60"),
      alpha = nd$alpha, label = nd$label, stringsAsFactors = FALSE)
    edge_style <- function(kind, style) {
      s <- ifelse(kind == "hasse" & style %in% names(pal$edge), style, "other")
      lapply(c("color", "linetype", "linewidth"), function(f) vapply(s, function(z) pal$edge[[z]][[f]], pal$edge$other[[f]]))
    }
    edges <- if (nrow(ed)) {
      es <- edge_style(ed$kind, ed$style_key)
      data.frame(from = ed$from_event, to = ed$to_event, arrow = ed$arrow,
                 color = unname(es[[1]]), linetype = unname(es[[2]]), linewidth = unname(es[[3]]),
                 label = ed$label, flag = ed$conflict, stringsAsFactors = FALSE)
    }
    groups <- if (nrow(gp)) data.frame(x = gp$x, y = gp$y, label = gp$label, stringsAsFactors = FALSE)
    xs <- c(nodes$x, groups$x)
    keys <- list()
    if (!is.null(nodes)) {
      nlab <- c(circle = "SV event", triangle_up = "Gain", triangle_down = "Loss", diamond = "WGD")
      for (s in intersect(names(nlab), nodes$shape)) {
        fc <- nodes$fill[nodes$shape == s][1]
        keys <- c(keys, list(LegendKey(nlab[[s]], shape = c(circle = "circle", triangle_up = "triangle",
                                                              triangle_down = "triangle", diamond = "diamond")[[s]],
                                       color = "black", fill = fc,
                                       pch = c(circle = 21L, triangle_up = 24L, triangle_down = 25L, diamond = 23L)[[s]])))
      }
    }
    if (!is.null(edges)) {
      sl <- c(snv = "SNV order", mixed = "Mixed order", structure = "Structural order")
      for (i in match(unique(paste(ed$kind, ed$style_key)), paste(ed$kind, ed$style_key))) {
        lab <- if (ed$kind[i] == "hasse" && ed$style_key[i] %in% names(sl)) sl[[ed$style_key[i]]] else ed$kind[i]
        keys <- c(keys, list(k_line(lab, edges$color[i], lwd = 1.5 * edges$linewidth[i], lty = edges$linetype[i])))
      }
      if (any(edges$flag)) keys <- c(keys, list(LegendKey("Order conflict", shape = "diamond",
                                                         color = "#D62728", fill = "#D62728")))
    }
    el <- if (!is.null(nodes) || !is.null(groups))
      seq_node_link(nodes = nodes, edges = edges, groups = groups, legend = leg(keys))
    trk <- seq_track(track_id = "poset", windows = GenomicRanges::GRanges("events", IRanges::IRanges(1, 1)),
                     direction = "under",
                     scale_x = seq_scale_continuous(limits = if (length(xs)) range(xs) + c(-0.6, 0.6) else c(0, 1)),
                     scale_y = seq_scale_continuous(limits = c(-0.03, 1.03), breaks = seq(0, 1, 0.2)),
                     aesthetics = aes("axis.x.visible" = FALSE, "axis.y.title.text" = "Event time (tau)"),
                     elements = if (!is.null(el)) list(el) else list())
    p <- seq_plot(aesthetics = aes(margins = list(top = 0, right = legend_width, bottom = 0.02, left = 0)), ...)
    return(p %+% trk)
  }

  # selections ----
  g <- tb$glyphs
  g_ok <- g[g$kind %in% c("amplicon", "driver") & g$glyph != "hidden" & !is.na(g$x), , drop = FALSE]
  g_span <- g_ok[!is.na(g_ok$start), , drop = FALSE]
  g_hid <- g[g$glyph == "hidden" & !is.na(g$x), , drop = FALSE]
  an <- tb$annotations
  av <- tb$sv_arches
  steps <- tb$cn_steps
  # genome-wide event rows (NA position) are drawn across every window
  gw <- g[!is.na(g$glyph) & g$glyph != "hidden" & is.na(g$chr) & is.na(g$x), , drop = FALSE]
  if (layout == "arches") gw <- gw[!gw$source_id %in% steps$event_id, , drop = FALSE]
  gw_key <- steps$colour_key[match(gw$source_id, steps$event_id)]
  gw_key <- ifelse(is.na(gw_key), tolower(tb$poset_nodes$shape_key[match(gw$source_id, tb$poset_nodes$event_id)]), gw_key)
  gw_layer <- tb$poset_nodes$layer[match(gw$source_id, tb$poset_nodes$event_id)]
  if (layout == "order") { gw <- gw[!is.na(gw_layer), , drop = FALSE]; gw_key <- gw_key[!is.na(gw_layer)]; gw_layer <- gw_layer[!is.na(gw_layer)] }
  step_lab <- c(gain = "Gain", loss = "Loss", wgd = "WGD")
  gw_els <- function() {
    if (!nrow(gw)) return(list())
    ord <- layout == "order"
    ix <- rep(seq_len(nrow(gw)), each = nrow(w))
    iw <- rep(seq_len(nrow(w)), times = nrow(gw))
    col <- ifelse(is.na(gw_key), pick(pal$glyph, gw$colour_key, "grey55"), pick(pal$step, gw_key, "grey40"))
    lab <- ifelse(gw_key %in% names(step_lab), step_lab[gw_key], ifelse(is.na(gw_key), gw$colour_key, gw_key))
    ukeys <- !duplicated(lab)
    list(seq_step_ribbon(
      to_gr(data.frame(chr = w$chr[iw], start = w$start[iw], end = w$end[iw],
                       y = if (ord) gw_layer[ix] else gw$y[ix],
                       lo = if (ord) NA_real_ else gw$bar_lo[ix], hi = if (ord) NA_real_ else gw$bar_hi[ix],
                       col = col[ix], alpha = gw$alpha[ix], label = NA_character_, stringsAsFactors = FALSE)),
      mp(x0 = "start", x1 = "end", y = "y", label = "label", bar_lo = "lo", bar_hi = "hi", color = "col", alpha = "alpha"),
      legend = leg(Map(k_rect, lab[ukeys], col[ukeys]))))
  }
  glyph_lab <- c(point = "SNV time", latest = "Proxy bound: before", earliest = "Proxy bound: after",
                 interval = "Proxy bound: between", untimed = "SV, no timing")

  # annotation palette ----
  ann_keys <- unique(an$colour_key[an$colour_key != "grey"])
  ann_pal <- c(stats::setNames(flexoki_palette(max(length(ann_keys), 1))[seq_along(ann_keys)], ann_keys),
               grey = "grey70")
  ann_pal[names(pal$annotation)] <- unlist(pal$annotation)

  # sv track ----
  sv_els <- list()
  lay_max <- max(c(0, av$layer, gw_layer), na.rm = TRUE)
  sv_ylim <- if (layout == "order") c(0, lay_max + 0.3) else c(0, 1)
  sv_scale <- if (layout == "order") {
    seq_scale_continuous(limits = sv_ylim, breaks = 0:lay_max)
  } else tau_scale
  gtypes <- c(av$glyph[layout != "order"], g_ok$glyph)
  gkeys <- seq_glyph_legend_keys(unique(gtypes), labels = glyph_lab)
  has_arch <- nrow(av) > 0
  if (layout != "order" && nrow(g_ok)) {
    sv_els <- c(sv_els, list(seq_interval_glyph(
      to_gr(data.frame(chr = g_ok$chr, start = g_ok$x, end = g_ok$x, y = g_ok$y, lo = g_ok$bar_lo,
                       hi = g_ok$bar_hi, glyph = g_ok$glyph, col = pick(pal$glyph, g_ok$colour_key, "grey55"),
                       alpha = g_ok$alpha, label = g_ok$label, caution = g_ok$ecdna_caution,
                       stringsAsFactors = FALSE)),
      mp(y = "y", bar_lo = "lo", bar_hi = "hi", glyph = "glyph", color = "col", alpha = "alpha",
          label = "label", caution = "caution"),
      legend = if (has_arch) FALSE else leg(unname(gkeys)))))
  }

  if (layout %in% c("arches", "order")) sv_els <- c(sv_els, gw_els())

  # cn track ----
  cn_els <- list(seq_paired_ribbon(
    to_gr(tb$cn_bands), mp(x0 = "start", x1 = "end", y_up = "cn_major", y_down = "cn_minor"),
    aesthetics = aes(fill_up = pal$cn[["major"]], fill_down = pal$cn[["minor"]]),
    legend = leg(list(k_rect("Major CN", pal$cn[["major"]]), k_rect("Minor CN", pal$cn[["minor"]])))))
  sn <- tb$snv_points
  if (nrow(sn)) {
    sbs_present <- unique(sn$sbs)
    cn_els <- c(cn_els, list(seq_point(
      to_gr(data.frame(chr = sn$chr, start = sn$pos, end = sn$pos, mult = sn$multiplicity,
                       col = pick(pal$sbs, sn$sbs, "grey50"), stringsAsFactors = FALSE)),
      mp(y = "mult", color = "col"), aes(size = 0.25, alpha = 0.6),
      legend = leg(lapply(sbs_present, function(s) k_pt(s, pal$sbs[[s]]))))))
  }
  if (has_arch) {
    ord <- layout == "order"
    ad <- data.frame(chr1 = av$chr1, pos1 = av$pos1, chr2 = av$chr2, pos2 = av$pos2,
                     top = if (ord) av$layer else av$y, glyph = av$glyph, lo = av$bar_lo, hi = av$bar_hi,
                     mult = av$multiplicity, col = pick(pal$sv, av$colour_key, pal$sv[["untimed"]]),
                     alpha = av$alpha, caution = av$ecdna_caution, id = av$arch_id, stringsAsFactors = FALSE)
    present <- intersect(names(pal$sv), unique(av$colour_key))
    present <- present[present != "untimed"]
    akeys <- lapply(present, function(s) k_line(s, pal$sv[[s]]))
    if (any(av$colour_key == "untimed")) akeys <- c(akeys, list(k_line("SV, no timing", pal$sv[["untimed"]])))
    if (!ord) akeys <- c(akeys, unname(gkeys[setdiff(names(gkeys), "untimed")]))
    if (any(av$ecdna_caution, na.rm = TRUE)) {
      akeys <- c(akeys, list(.caution_legend_key("ecDNA: timing not directly comparable with the chromosomal clock")))
    }
    cn_els <- c(cn_els, list(
      if (ord) seq_cross_arch(ad, mp(chrom0 = "chr1", x0 = "pos1", chrom1 = "chr2", x1 = "pos2", top = "top",
                                      marker = "mult", color = "col", alpha = "alpha", caution = "caution", id = "id"),
                              target = "sv", legend = leg(akeys))
      else seq_cross_arch(ad, mp(chrom0 = "chr1", x0 = "pos1", chrom1 = "chr2", x1 = "pos2", top = "top",
                                  glyph = "glyph", bar_lo = "lo", bar_hi = "hi", marker = "mult", color = "col",
                                  alpha = "alpha", caution = "caution", id = "id"),
                          target = "sv", aesthetics = aes(bar_alpha = 0.35, glyph_points = FALSE),
                          legend = leg(akeys))))
  }
  cn_max <- max(c(0, tb$cn_bands$cn_major, av$multiplicity), na.rm = TRUE)
  cn_log <- isTRUE(as.logical(st$cn_log_axis))

  # steps track ----
  steps_els <- list()
  if (nrow(steps)) {
    pt <- c(gain = "up", wgd = "up", loss = "down")
    ks <- sort(unique(steps$k))
    skeys <- list()
    if (any(steps$k > 1, na.rm = TRUE)) {
      tl <- c("1" = "1 copy", "2" = "2 copies", "3" = "3+ copies")
      skeys <- lapply(ks, function(k) k_line(if (as.character(k) %in% names(tl)) tl[[as.character(k)]] else paste(k, "copies"),
                                             "grey40", lwd = 2 * k))
    }
    dl <- c(gain = "Gain", loss = "Loss", wgd = "WGD")
    for (d in intersect(names(dl), unique(steps$colour_key))) skeys <- c(skeys, list(k_rect(dl[[d]], pal$step[[d]])))
    if (any(steps$credit < 1, na.rm = TRUE)) skeys <- c(skeys, list(k_rect("Partial credit", "grey40", 0.5)))
    if (any(steps$faded, na.rm = TRUE)) skeys <- c(skeys, list(k_rect("Proxy bound", "grey40", 0.3)))
    steps_els <- list(seq_step_ribbon(
      to_gr(data.frame(chr = steps$chr, start = steps$start, end = steps$end, y = steps$y,
                       k = steps$k, label = steps$label, lo = steps$bar_lo, hi = steps$bar_hi,
                       col = pick(pal$step, steps$colour_key, "grey40"), alpha = steps$alpha,
                       pointing = unname(pt[steps$direction]), stringsAsFactors = FALSE)),
      mp(x0 = "start", x1 = "end", y = "y", thickness = "k", label = "label", bar_lo = "lo", bar_hi = "hi",
          color = "col", alpha = "alpha", pointing = "pointing"),
      legend = leg(skeys)))
  }

  # heat track ----
  heat_els <- list()
  hm <- tb$heatmap
  if (layout == "heatmap") {
    if (nrow(hm)) {
      has_dens <- any(hm$dens_gain > 0 | hm$dens_loss > 0, na.rm = TRUE)
      dmax <- sprintf("%.1f", st$d_max)
      key2d <- if (has_dens) seq_legend_2d(
        pal$heatmap[["gain"]], pal$heatmap[["loss"]],
        x_title = "Gain density (per tau)", y_title = "Loss density (per tau)",
        x_labels = c("0", dmax), y_labels = c("0", dmax),
        note = sprintf("smoothing sd = %.3g tau", st$smooth_sd))
      heat_els <- list(seq_tile_mix(
        to_gr(hm), mp(x0 = "start", x1 = "end", y_min = "tau_lo", y_max = "tau_hi", a = "intensity_gain", b = "intensity_loss"),
        aesthetics = aes(color_a = pal$heatmap[["gain"]], color_b = pal$heatmap[["loss"]]),
        legend = if (has_dens) key2d else FALSE))
    }
    hb <- tb$heatmap_bounds
    if (nrow(hb)) {
      bl <- c(gain = "Gain, proxy bound", loss = "Loss, proxy bound")
      hkeys <- lapply(intersect(names(bl), unique(hb$direction)), function(d) k_line(bl[[d]], pal$heatmap[[d]], lwd = 2.5))
      heat_els <- c(heat_els, list(seq_segment(
        to_gr(data.frame(chr = hb$chr, start = hb$start, end = hb$end, y = hb$y,
                         col = pick(pal$heatmap, hb$direction, "grey40"), stringsAsFactors = FALSE)),
        mp(x = "start", x_end = "end", y = "y", y_end = "y", color = "col"), aes(linewidth = 1.2),
        legend = leg(hkeys))))
    }
  }

  if (layout == "heatmap") heat_els <- c(heat_els, gw_els())

  # labels and hidden tracks ----
  lab_els <- list()
  if (nrow(g_span)) {
    lab_els <- c(lab_els, list(seq_span_strip(
      to_gr(data.frame(chr = g_span$chr, start = g_span$start, end = g_span$end, label = g_span$label,
                       stringsAsFactors = FALSE)),
      mp(x0 = "start", x1 = "end", label = "label"), aes(bar_y = 0.2))))
  }
  if (nrow(an)) {
    lab_els <- c(lab_els, list(seq_span_strip(
      to_gr(data.frame(chr = an$chr, start = an$start, end = an$end, label = an$label,
                       fill = pick(ann_pal, an$colour_key, "grey70"), stringsAsFactors = FALSE)),
      mp(x0 = "start", x1 = "end", label = "label", fill = "fill"), aes(bar_y = 0.68),
      legend = leg(lapply(ann_keys, function(k) k_rect(k, pick(ann_pal, k, "grey70")))))))
  }
  hid_els <- list()
  if (show_hidden && nrow(g_hid)) {
    hid_els <- list(seq_span_strip(
      to_gr(data.frame(chr = g_hid$chr, start = g_hid$x, end = g_hid$x, label = g_hid$label, stringsAsFactors = FALSE)),
      mp(x0 = "start", x1 = "end", label = "label"), aes(bar_height = 0)))
  }

  # assemble tracks ----
  cn_title <- paste0("Copy number / multiplicity", if (cn_log) " (log scale)" else "")
  sv_title <- if (layout == "order") "SV event order (poset layer; 0 = unordered)" else "SV event time (tau)"
  unit_y <- seq_scale_continuous(limits = c(0, 1))
  spec <- list()
  if (length(lab_els)) spec$labels <- list(els = lab_els, scale = unit_y, title = NULL, height = 0.3)
  spec$sv <- list(els = sv_els, scale = sv_scale, title = sv_title, height = 1)
  spec$cn <- list(els = cn_els, scale = seq_scale_count(cn_max, log = cn_log), title = cn_title, height = 1)
  if (layout == "arches") spec$steps <- list(els = steps_els, scale = tau_scale, title = "CN event time (tau)", height = 1)
  if (layout == "heatmap") spec$heat <- list(els = heat_els, scale = tau_scale, title = "CN event time (tau)", height = 1.6)
  if (length(hid_els)) spec$hidden <- list(els = hid_els, scale = unit_y, title = NULL, height = 0.12)
  last <- names(spec)[length(spec)]
  # hidden y axis parts keep the panel x extent aligned across tracks
  blank_y <- list("axis.y.line.visible" = FALSE, "axis.y.ticks.visible" = FALSE,
                  "axis.y.text.visible" = FALSE, "axis.y.title.text" = " ")
  flat_y <- list(top = 0, bottom = 0, left = 0.02, right = 0.02)
  trk <- lapply(names(spec), function(id) {
    s <- spec[[id]]
    a <- list("axis.x1.visible" = plain && id == last)
    a <- c(a, if (is.null(s$title)) blank_y else list("axis.y.title.text" = s$title))
    seq_track(track_id = id, windows = wins, direction = "under", scale_y = s$scale,
              track_height = s$height, aesthetics = do.call(aes, a), elements = s$els,
              track_outer_margin = if (is.null(s$title)) flat_y else 0.02,
              track_inner_margin = if (is.null(s$title)) flat_y else 0.02,
              window_inner_margin = if (is.null(s$title)) flat_y else 0.02)
  })
  if (!plain) {
    cb <- cytobands %||% load_cytobands()
    ide_m <- list(top = 0, bottom = 0.02, left = 0.02, right = 0.02)
    trk <- c(trk, list(seq_track(track_id = "ideogram", windows = wins, direction = "under", track_height = 0.45,
                                 track_outer_margin = ide_m, track_inner_margin = ide_m, window_inner_margin = ide_m,
                                 aesthetics = do.call(aes, blank_y), elements = list(seq_ideogram(data = cb)))))
  }
  p <- seq_plot(aesthetics = aes(margins = list(top = if (length(lab_els)) 0.04 else 0, right = legend_width, bottom = 0.02, left = 0)), ...)
  for (t in trk) p <- p %+% t
  p
}
