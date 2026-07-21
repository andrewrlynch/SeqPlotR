# ── SeqGeneR6 ─────────────────────────────────────────────────────────────────
#
# Composite element: gene models with backbone lines, exon/UTR boxes,
# directional arrows, and labels. Format-agnostic — all column references
# come from the user's `map()` call; nothing is hard-coded.

# Parse a `gene.label$position` value into x_pos / y_pos slots.
#
# `position` may be NULL, a single string from c("start","end","top","bottom"),
# or a length-2 vector pairing one x value with one y value. Missing axes
# fall back to "strand" (x) and "center" (y).
.parse_gene_label_pos <- function(position) {
  if (is.null(position)) return(list(x_pos = "strand", y_pos = "center"))
  x_vals <- c("start", "end")
  y_vals <- c("top", "bottom")
  x_pos  <- intersect(position, x_vals)
  y_pos  <- intersect(position, y_vals)
  list(
    x_pos = if (length(x_pos)) x_pos[1] else "strand",
    y_pos = if (length(y_pos)) y_pos[1] else "center"
  )
}

# Measure the rendered width of each label, in npc, at expansion `cex`.
#
# `grid::stringWidth()` returns a "strwidth" unit whose stored numeric part is
# a constant multiplier (always 1 per string); the physical width only
# materialises when the unit is *converted*. So `as.numeric()` on it yields 1
# for every string regardless of length — `grid::convertWidth()` is required.
# The result is npc relative to whatever viewport is active at prep() time (the
# default `pushViewport(viewport())` here does not rescale from its parent),
# which is device-level npc — consistent with how `pm$inner$x0/x1` are used.
# The tryCatch fallback only fires on a genuine device error (prep() with no
# open device auto-opens a null device rather than throwing).
.measure_label_npc <- function(labels, cex = 1) {
  if (length(labels) == 0L) return(numeric(0))
  tryCatch({
    grid::pushViewport(grid::viewport())
    on.exit(grid::popViewport())
    grid::convertWidth(grid::stringWidth(labels) * cex, "npc",
                       valueOnly = TRUE)
  }, error = function(e) nchar(labels) * 0.01 * cex)
}

#' SeqGene R6 class
#'
#' Internal R6 generator backing [seq_gene()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqGeneR6 <- R6::R6Class("SeqGene",
  inherit = SeqElementR6,
  public = list(
    #' @field exon_height Proportion of per-tier height allocated to the
    #'   exon box.
    exon_height   = 0.8,
    #' @field label_pad Extra padding (in bp) added, on top of the measured
    #'   rendered label width, to each gene's packing interval on the label
    #'   side. Only affects `label_position = "flush"`. Default `50000`.
    label_pad     = 50000,

    #' @field label_position Character. One of `"flush"`, `"above"`, `"below"`.
    #'   Controls where the gene label is drawn. Default `"flush"`.
    label_position = "flush",

    #' @field label_row_frac Numeric in `[0, 1)`. Fraction of each tier row
    #'   reserved for the label when `label_position` is `"above"`/`"below"`;
    #'   the exon/gene glyph height is shrunk accordingly. Default `0.35`.
    label_row_frac = 0.35,

    #' @field pack_direction Character. One of `"bottom_to_top"`,
    #'   `"top_to_bottom"`. Physical edge tier 1 is drawn from. Default
    #'   `"bottom_to_top"`.
    pack_direction = "bottom_to_top",

    #' @field expanded Logical. When `TRUE`, every gene gets its own tier
    #'   (one gene per row) regardless of genomic overlap, consistently across
    #'   all windows in the track. Default `FALSE`.
    expanded = FALSE,
    #' @field label_cex Expansion factor for gene label text.
    label_cex     = 0.6,
    #' @field label_offset Horizontal offset (in npc units) applied when
    #'   placing labels flush against the gene backbone.
    label_offset  = 0.01,

    #' @field backbone_type Character. One of `"arrow"`, `"solid"`, `"dashed"`.
    #'   Controls how the gene backbone line is rendered. Default `"arrow"`.
    backbone_type = "arrow",

    #' @field show_start Logical. When `TRUE`, a TSS flag arrow is drawn above
    #'   the first exon (by row order) of each gene. Default `FALSE`.
    show_start = FALSE,

    #' @field tss_position Named list of length-2 numeric vectors. Keys are gene
    #'   IDs (matching `group` values); values are `c(start, end)` genomic
    #'   coordinates overriding the auto-detected first-exon TSS position.
    #'   Default `NULL` (auto from first exon by row order).
    tss_position = NULL,

    #' @field separate_strands Logical. When `TRUE`, genes are placed in two
    #'   horizontal sub-bands by strand (`"+"` top, `"-"` bottom), each labelled.
    #'   Silently ignored when all strands are `"*"` or only one strand is
    #'   present. Default `FALSE`.
    separate_strands = FALSE,

    #' @field style_type Character. One of `"exon"`, `"gene"`, `"point"`.
    #'   Selects the per-gene rendering style: `"exon"` draws backbone +
    #'   exon/UTR boxes (default); `"gene"` draws a single chevron-shaped
    #'   polygon spanning the full gene extent; `"point"` draws a single
    #'   filled circle at the TSS. Default `"exon"`.
    style_type = "exon",

    #' @field backbone_lwd Numeric. Backbone line width; also the default
    #'   stroke width of the directional chevron arrows. Default `1`.
    backbone_lwd = 1,

    #' @field exon_lwd Numeric or `NA`. Exon/UTR box border width. When `NA`
    #'   (default) boxes draw with no border (`col = NA`, current behaviour);
    #'   when set, the border is stroked.
    exon_lwd = NA,

    #' @field exon_border_col Border colour for exon/UTR boxes, used only when
    #'   `exon_lwd` is not `NA`. Default `NA` (inherits the gene colour).
    exon_border_col = NA,

    #' @field arrow_size Numeric. Chevron arrowhead length in millimetres.
    #'   Default `1.5`.
    arrow_size = 1.5,

    #' @field tss_arrow_size Numeric. TSS marker arrowhead (filled triangle)
    #'   size in millimetres. Default `3`.
    tss_arrow_size = 3,

    #' @field tss_stem_lwd Numeric. Line width of the TSS marker stem/arm.
    #'   Default `1.2`.
    tss_stem_lwd = 1.2,

    #' @field tss_rise Numeric or `NULL`. Vertical stem rise of the TSS marker
    #'   in npc. When `NULL` (default) it is derived from the gene box height
    #'   (`0.9 *` box height, `0.02` fallback), matching prior behaviour.
    tss_rise = NULL,

    #' @field tss_run Numeric. Horizontal arm run of the TSS marker in npc.
    #'   Default `0.018`.
    tss_run = 0.018,

    #' @field oob Character. Out-of-bounds handling for genes straddling a
    #'   window edge: `"keep"`, `"truncate"` (default), or `"remove"`.
    oob = "truncate",

    #' @description Construct a SeqGeneR6.
    #' @param data Optional `GRanges`.
    #' @param mapping Optional `SeqMap`. Recognised: `group`, `strand`,
    #'   `label`, `type`, `color`.
    #' @param aesthetics Optional `SeqAes`. Supports `color` (default
    #'   `"gray30"` when no `color` mapping is given), `linewidth`, `alpha`.
    #' @param backbone_type Character. One of `"arrow"`, `"solid"`, `"dashed"`.
    #'   Default `"arrow"`. Ignored when `style_type` is `"gene"` or `"point"`.
    #' @param show_start Logical. Draw TSS flag arrow above first exon. Default
    #'   `FALSE`. Ignored when `style_type = "point"` (the point itself marks
    #'   the TSS).
    #' @param tss_position Named list overriding per-gene TSS genomic positions.
    #'   Default `NULL`.
    #' @param separate_strands Logical. Split track into `"+"` and `"-"`
    #'   sub-bands. Default `FALSE`.
    #' @param style_type Character. One of `"exon"` (default), `"gene"`, or
    #'   `"point"`. Controls per-gene rendering: full backbone + exon boxes,
    #'   single chevron polygon, or a single TSS point.
    #' @param label_position Character. One of `"flush"` (default), `"above"`,
    #'   `"below"`. Where the gene label is drawn.
    #' @param label_row_frac Numeric in `[0, 1)`. Fraction of each tier row
    #'   reserved for `"above"`/`"below"` labels. Default `0.35`.
    #' @param pack_direction Character. One of `"bottom_to_top"` (default) or
    #'   `"top_to_bottom"`. Physical edge tier 1 is drawn from.
    #' @param expanded Logical. When `TRUE`, one gene per tier regardless of
    #'   overlap, consistent across all windows. Default `FALSE`.
    #' @param label_pad Numeric. Extra bp padding beyond the measured label
    #'   width, applied on the label side in `"flush"` mode. Default `50000`.
    #' @param backbone_lwd Numeric. Backbone (and chevron) line width. Default `1`.
    #' @param exon_lwd Numeric or `NA`. Exon/UTR box border width. Default `NA`.
    #' @param exon_border_col Exon/UTR border colour when `exon_lwd` is set.
    #'   Default `NA` (inherits gene colour).
    #' @param arrow_size Numeric. Chevron arrowhead length (mm). Default `1.5`.
    #' @param tss_arrow_size Numeric. TSS triangle size (mm). Default `3`.
    #' @param tss_stem_lwd Numeric. TSS stem/arm line width. Default `1.2`.
    #' @param tss_rise Numeric or `NULL`. TSS stem rise (npc); `NULL` derives it
    #'   from box height. Default `NULL`.
    #' @param tss_run Numeric. TSS arm run (npc). Default `0.018`.
    #' @param oob Character. Out-of-bounds gene handling: `"keep"`,
    #'   `"truncate"` (default), or `"remove"`.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          backbone_type    = "arrow",
                          show_start       = FALSE,
                          tss_position     = NULL,
                          separate_strands = FALSE,
                          style_type       = c("exon", "gene", "point"),
                          label_position   = c("flush", "above", "below"),
                          label_row_frac   = 0.35,
                          pack_direction   = c("bottom_to_top", "top_to_bottom"),
                          expanded         = FALSE,
                          label_pad        = 50000,
                          backbone_lwd     = 1,
                          exon_lwd         = NA,
                          exon_border_col  = NA,
                          arrow_size       = 1.5,
                          tss_arrow_size   = 3,
                          tss_stem_lwd     = 1.2,
                          tss_rise         = NULL,
                          tss_run          = 0.018,
                          oob              = c("truncate", "keep", "remove"),
                          ...) {
      super$initialize(data, mapping, aesthetics)
      backbone_type  <- match.arg(backbone_type, c("arrow", "solid", "dashed"))
      style_type     <- match.arg(style_type)
      label_position <- match.arg(label_position)
      pack_direction <- match.arg(pack_direction)
      oob            <- match.arg(oob)
      self$backbone_type    <- backbone_type
      self$show_start       <- isTRUE(show_start)
      self$tss_position     <- tss_position
      self$separate_strands <- isTRUE(separate_strands)
      self$style_type       <- style_type
      self$label_position   <- label_position
      self$label_row_frac   <- label_row_frac
      self$pack_direction   <- pack_direction
      self$expanded         <- isTRUE(expanded)
      self$label_pad        <- label_pad
      self$backbone_lwd     <- backbone_lwd
      self$exon_lwd         <- exon_lwd
      self$exon_border_col  <- exon_border_col
      self$arrow_size       <- arrow_size
      self$tss_arrow_size   <- tss_arrow_size
      self$tss_stem_lwd     <- tss_stem_lwd
      self$tss_rise         <- tss_rise
      self$tss_run          <- tss_run
      self$oob              <- oob
    },

    #' @description Resolve mappings, stack genes into non-overlapping tiers,
    #'   and compute canvas coordinates for backbones, exons/UTRs, arrows,
    #'   and labels. Populates `self$coordCanvas` with a single
    #'   `data.frame` whose rows are per-feature draw primitives.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$resolve(
        track_data    = layout_track[[1]]$track_data,
        track_mapping = layout_track[[1]]$track_mapping
      )
      eff_data <- self$resolved$.data
      if (is.null(eff_data) || length(eff_data) == 0) {
        self$coordCanvas <- data.frame()
        return(invisible())
      }

      n <- length(eff_data)
      # Treat function-valued resolutions as "not mapped" — this happens when
      # a map() expression like `strand = str` unintentionally resolves to a
      # base R function because no matching mcols column is present.
      .usable <- function(v) !is.null(v) && !is.function(v)
      pick <- function(v, fallback) if (.usable(v)) v else fallback

      group_vals  <- as.character(pick(self$resolved$group,  seq_along(eff_data)))
      strand_vals <- as.character(pick(self$resolved$strand, rep("+", n)))
      label_vals  <- as.character(pick(self$resolved$label,  group_vals))
      type_vals   <- as.character(pick(self$resolved$type,   rep("exon", n)))
      color_vals  <- as.character(pick(
        self$resolved$color,
        rep(self$aesthetics$color %||% "gray30", n)
      ))

      # Normalise strand values outside of "+/-" to "+".
      strand_vals[is.na(strand_vals) | !strand_vals %in% c("+", "-")] <- "+"

      all_rows <- list()
      nWin    <- length(layout_track)

      # ── Track-level packing setup ───────────────────────────────────────
      # Strand separation is decided once for the whole track (not per window)
      # so that band geometry and row heights are consistent across windows.
      use_sep_strands <- isTRUE(self$separate_strands) &&
        length(unique(strand_vals[strand_vals %in% c("+", "-")])) > 1L

      # Global gene table — one entry per unique gene id in the whole track.
      # Provides each gene's full genomic extent (used for oob tests) and
      # label, plus a once-measured rendered label width in npc.
      all_gids   <- unique(group_vals)
      g_all_start <- BiocGenerics::start(eff_data)
      g_all_end   <- BiocGenerics::end(eff_data)
      g_start  <- integer(length(all_gids))
      g_end    <- integer(length(all_gids))
      g_strand <- character(length(all_gids))
      g_label  <- character(length(all_gids))
      for (i in seq_along(all_gids)) {
        m <- group_vals == all_gids[i]
        g_start[i]  <- min(g_all_start[m])
        g_end[i]    <- max(g_all_end[m])
        g_strand[i] <- strand_vals[m][1]
        g_label[i]  <- label_vals[m][1]
      }

      # A1: rendered label width per gene (npc), measured once. Only used to
      # pad packing intervals in `label_position = "flush"`; the width is
      # device-relative npc (window-independent) — converted to bp per window.
      label_npc_w <- setNames(
        .measure_label_npc(g_label, self$label_cex), all_gids)

      # Greedy interval-scheduling tier assignment for a group of genes.
      # `idx` indexes into p0/p1; returns tiers aligned with `idx`.
      .pack_greedy <- function(idx, p0, p1) {
        out       <- integer(length(idx))
        ends_last <- numeric(0)
        for (k in order(p0[idx])) {
          s <- p0[idx[k]]; e <- p1[idx[k]]
          w2 <- which(ends_last < s)
          if (length(w2) == 0) {
            out[k]    <- length(ends_last) + 1L
            ends_last <- c(ends_last, e)
          } else {
            out[k]           <- w2[1]
            ends_last[w2[1]] <- e
          }
        }
        out
      }

      # ── Pass 1: per-window tier assignment (no canvas geometry yet) ──────
      win_ctx <- vector("list", nWin)
      for (w in seq_len(nWin)) {
        pm  <- layout_track[[w]]
        win <- track_windows[w]

        hits <- GenomicRanges::findOverlaps(eff_data, win)
        if (length(hits) == 0) next
        keep <- unique(S4Vectors::queryHits(hits))

        exons_start <- g_all_start[keep]
        exons_end   <- g_all_end[keep]
        gid_all     <- group_vals[keep]
        label_all   <- label_vals[keep]
        type_all    <- type_vals[keep]

        gids <- unique(gid_all)
        if (length(gids) == 0) next

        gi_global  <- match(gids, all_gids)
        gene_start <- g_start[gi_global]
        gene_end   <- g_end[gi_global]
        strand_g   <- g_strand[gi_global]
        label_g    <- g_label[gi_global]
        color_g    <- vapply(gids, function(gg) color_vals[group_vals == gg][1],
                             character(1))

        # A6: out-of-bounds handling relative to this window, applied *before*
        # packing so tiers reflect the clipped/removed set.
        wstart <- BiocGenerics::start(win)
        wend   <- BiocGenerics::end(win)
        if (identical(self$oob, "remove")) {
          inside <- gene_start >= wstart & gene_end <= wend
          if (!any(inside)) next
          gids       <- gids[inside];       gene_start <- gene_start[inside]
          gene_end   <- gene_end[inside];   strand_g   <- strand_g[inside]
          label_g    <- label_g[inside];    color_g    <- color_g[inside]
        } else if (identical(self$oob, "truncate")) {
          gene_start  <- pmax(gene_start, wstart)
          gene_end    <- pmin(gene_end,   wend)
          exons_start <- pmax(exons_start, wstart)
          exons_end   <- pmin(exons_end,   wend)
        }
        # "keep": leave extents untouched (genes may run past the window edge).

        xpr <- pm$xplot_range %||% pm$xscale

        # Packing intervals start from the (possibly clipped) body extent.
        p0 <- gene_start
        p1 <- gene_end
        pad_bp <- numeric(length(gids))
        # A1: in flush mode, widen each gene's packing interval on the label
        # side by its rendered label width (converted to bp) plus label_pad, so
        # genes whose *labels* collide are pushed to separate tiers. The body
        # extents used for drawing (gene_start/gene_end) stay unpadded.
        if (identical(self$label_position, "flush")) {
          panel_w    <- pm$inner$x1 - pm$inner$x0
          bp_per_npc <- if (isTRUE(panel_w > 0)) diff(xpr) / panel_w else 0
          pad_bp     <- as.numeric(label_npc_w[gids]) * bp_per_npc + self$label_pad
          is_minus   <- strand_g == "-"
          # +/* genes draw their flush label on the left (5'); - on the right.
          p0[!is_minus] <- p0[!is_minus] - pad_bp[!is_minus]
          p1[ is_minus] <- p1[ is_minus] + pad_bp[ is_minus]
        }

        tiers <- integer(length(gids))
        if (isTRUE(self$expanded)) {
          # A2: per-window ranking — one gene per tier, ordered by start (ties
          # by label), over the genes present in THIS window only. A gene
          # absent from a window does not reserve a row there.
          if (use_sep_strands) {
            for (st in c("+", "-")) {
              sel <- which(strand_g == st)
              if (!length(sel)) next
              ordv <- sel[order(gene_start[sel], label_g[sel])]
              tiers[ordv] <- seq_along(ordv)
            }
            ntiers_plus  <- max(1L, sum(strand_g == "+"))
            ntiers_minus <- max(1L, sum(strand_g == "-"))
            ntiers_flat  <- max(tiers)
          } else {
            ordv <- order(gene_start, label_g)
            tiers[ordv]  <- seq_along(ordv)
            ntiers_plus  <- length(gids)
            ntiers_minus <- length(gids)
            ntiers_flat  <- length(gids)
          }
        } else if (use_sep_strands) {
          for (band_gids in list(gids[strand_g == "+"], gids[strand_g == "-"])) {
            if (length(band_gids) == 0) next
            idx <- match(band_gids, gids)
            tiers[idx] <- .pack_greedy(idx, p0, p1)
          }
          ntiers_plus  <- if (any(strand_g == "+"))
            max(tiers[strand_g == "+"]) else 1L
          ntiers_minus <- if (any(strand_g == "-"))
            max(tiers[strand_g == "-"]) else 1L
          ntiers_flat  <- max(tiers)
        } else {
          tiers        <- .pack_greedy(seq_along(gids), p0, p1)
          ntiers_plus  <- max(tiers)
          ntiers_minus <- max(tiers)
          ntiers_flat  <- max(tiers)
        }

        u0 <- (gene_start - xpr[1]) / diff(xpr)
        u1 <- (gene_end   - xpr[1]) / diff(xpr)
        x0c <- pm$inner$x0 + u0 * (pm$inner$x1 - pm$inner$x0)
        x1c <- pm$inner$x0 + u1 * (pm$inner$x1 - pm$inner$x0)

        # TSS detection for show_start.
        tss_genomic <- setNames(numeric(length(gids)), gids)
        for (i in seq_along(gids)) {
          gid <- gids[i]
          m   <- gid_all == gid
          if (!is.null(self$tss_position) &&
              gid %in% names(self$tss_position)) {
            tss_genomic[gid] <- self$tss_position[[gid]][1]
          } else if (strand_g[i] == "-") {
            tss_genomic[gid] <- max(exons_end[m])
          } else {
            tss_genomic[gid] <- min(exons_start[m])
          }
        }
        tss_npc <- (tss_genomic - xpr[1]) / diff(xpr) *
                   (pm$inner$x1 - pm$inner$x0) + pm$inner$x0

        win_ctx[[w]] <- list(
          pm = pm, xpr = xpr, gids = gids, strand_g = strand_g,
          label_g = label_g, color_g = color_g, tiers = tiers, pad_bp = pad_bp,
          ntiers_plus = ntiers_plus, ntiers_minus = ntiers_minus,
          ntiers_flat = ntiers_flat, x0c = x0c, x1c = x1c, tss_npc = tss_npc,
          gid_all = gid_all, type_all = type_all,
          exons_start = exons_start, exons_end = exons_end
        )
      }

      # ── Global tier counts (constant row height across the track) ────────
      # Always the max per-window tier count — including the `expanded` path
      # (A2), which now ranks per window. Glyph height stays constant across
      # windows while packing/rows are window-local.
      present <- Filter(Negate(is.null), win_ctx)
      if (length(present) == 0L) {
        self$coordCanvas <- data.frame()
        return(invisible())
      }
      ntiers_global       <- max(vapply(present, function(c) c$ntiers_flat,  integer(1)))
      ntiers_plus_global  <- max(vapply(present, function(c) c$ntiers_plus,  integer(1)))
      ntiers_minus_global <- max(vapply(present, function(c) c$ntiers_minus, integer(1)))

      # Glyph vertical geometry for one gene, using the *global* tier count so
      # thickness is identical across windows. Honours pack_direction and,
      # for non-flush labels, reserves label_row_frac of the row.
      .glyph_geom <- function(band_y0, band_y1, nt, tier) {
        nt    <- max(1L, nt)
        row_h <- (band_y1 - band_y0) / nt
        flush <- identical(self$label_position, "flush")
        exon_scale <- self$exon_height *
          (if (flush) 1 else (1 - self$label_row_frac))
        exon_h <- row_h * exon_scale
        if (identical(self$pack_direction, "top_to_bottom")) {
          row_hi <- band_y1 - (tier - 1) * row_h
          row_lo <- row_hi - row_h
        } else {
          row_lo <- band_y0 + (tier - 1) * row_h
          row_hi <- row_lo + row_h
        }
        if (flush) {
          ymid <- if (identical(self$pack_direction, "top_to_bottom"))
            row_hi - exon_h / 2 else row_lo + exon_h / 2
          laby <- NA_real_
        } else if (identical(self$label_position, "above")) {
          frac <- self$label_row_frac
          ymid <- row_lo + (1 - frac) * row_h / 2
          laby <- row_hi - frac * row_h / 2
        } else {  # "below"
          frac <- self$label_row_frac
          ymid <- row_hi - (1 - frac) * row_h / 2
          laby <- row_lo + frac * row_h / 2
        }
        list(ymid = ymid, exon_h = exon_h, laby = laby)
      }

      # ── Pass 2: build draw rows against the shared global row height ──────
      lpos <- self$label_position
      for (w in seq_len(nWin)) {
        ctx <- win_ctx[[w]]
        if (is.null(ctx)) next
        pm       <- ctx$pm
        xpr      <- ctx$xpr
        gids     <- ctx$gids
        strand_g <- ctx$strand_g
        label_g  <- ctx$label_g
        color_g  <- ctx$color_g
        tiers    <- ctx$tiers
        pad_vec  <- ctx$pad_bp
        x0c      <- ctx$x0c
        x1c      <- ctx$x1c
        tss_npc  <- ctx$tss_npc
        gid_all  <- ctx$gid_all
        type_all <- ctx$type_all
        exons_start <- ctx$exons_start
        exons_end   <- ctx$exons_end

        track_h <- pm$inner$y1 - pm$inner$y0
        if (use_sep_strands) {
          plus_y0  <- pm$inner$y0 + track_h / 2; plus_y1  <- pm$inner$y1
          minus_y0 <- pm$inner$y0;               minus_y1 <- pm$inner$y0 + track_h / 2
        } else {
          plus_y0  <- pm$inner$y0; plus_y1  <- pm$inner$y1
          minus_y0 <- pm$inner$y0; minus_y1 <- pm$inner$y1
        }

        ymid       <- numeric(length(gids))
        exon_h_vec <- numeric(length(gids))
        laby_vec   <- numeric(length(gids))
        for (i in seq_along(gids)) {
          g <- if (use_sep_strands && strand_g[i] == "+")
                 .glyph_geom(plus_y0,  plus_y1,  ntiers_plus_global,  tiers[i])
               else if (use_sep_strands && strand_g[i] == "-")
                 .glyph_geom(minus_y0, minus_y1, ntiers_minus_global, tiers[i])
               else
                 .glyph_geom(pm$inner$y0, pm$inner$y1, ntiers_global, tiers[i])
          ymid[i]       <- g$ymid
          exon_h_vec[i] <- g$exon_h
          laby_vec[i]   <- g$laby
        }

        for (i in seq_along(gids)) {
          gene_id <- gids[i]
          m       <- gid_all == gene_id
          type_sub  <- type_all[m]
          start_sub <- exons_start[m]
          end_sub   <- exons_end[m]

          # Height computed in Pass 2 against the shared global row count.
          exon_h <- exon_h_vec[i]
          laby_i <- laby_vec[i]
          pad_i  <- pad_vec[i]

          # ── style_type = "gene": one chevron polygon row per gene ───────
          if (self$style_type == "gene") {
            x0b_i <- x0c[i]
            x1b_i <- x1c[i]
            chev  <- min(0.012, max(0, (x1b_i - x0b_i) * 0.3))
            ey0   <- ymid[i] - exon_h / 2
            ey1   <- ymid[i] + exon_h / 2
            all_rows[[length(all_rows) + 1]] <- data.frame(
              gene            = gene_id,
              x0b             = x0b_i, x1b = x1b_i,
              ymid            = ymid[i],
              dir             = ifelse(strand_g[i] == "-", -1, 1),
              label           = label_g[i],
              color           = color_g[i],
              exon_x0         = NA_real_,
              exon_x1         = NA_real_,
              exon_y0         = ey0,
              exon_y1         = ey1,
              draw_box        = FALSE,
              strand          = strand_g[i],
              tss_x           = tss_npc[gene_id],
              chevron_npc     = chev,
              use_sep_strands = use_sep_strands,
              label_pos       = lpos,
              label_y         = laby_i,
              pad_bp          = pad_i,
              panel_y0        = pm$inner$y0,
              panel_y1        = pm$inner$y1,
              stringsAsFactors = FALSE,
              row.names       = NULL
            )
            next
          }

          # ── style_type = "point": one TSS point per gene ────────────────
          if (self$style_type == "point") {
            all_rows[[length(all_rows) + 1]] <- data.frame(
              gene            = gene_id,
              x0b             = x0c[i], x1b = x1c[i],
              ymid            = ymid[i],
              dir             = ifelse(strand_g[i] == "-", -1, 1),
              label           = label_g[i],
              color           = color_g[i],
              exon_x0         = NA_real_,
              exon_x1         = NA_real_,
              exon_y0         = NA_real_,
              exon_y1         = NA_real_,
              draw_box        = FALSE,
              strand          = strand_g[i],
              tss_x           = tss_npc[gene_id],
              chevron_npc     = NA_real_,
              use_sep_strands = use_sep_strands,
              label_pos       = lpos,
              label_y         = laby_i,
              pad_bp          = pad_i,
              panel_y0        = pm$inner$y0,
              panel_y1        = pm$inner$y1,
              stringsAsFactors = FALSE,
              row.names       = NULL
            )
            next
          }

          # ── style_type = "exon" (default): per-feature rows ─────────────
          # Per-feature height scaling: exon full, UTR 80%, other skipped.
          height_scale <- rep(NA_real_, length(type_sub))
          height_scale[type_sub == "exon"] <- 1
          height_scale[type_sub == "UTR"]  <- 0.8
          draw_mask <- !is.na(height_scale)
          if (!any(draw_mask)) {
            # Still emit the gene's backbone row with a placeholder exon
            # entry (zero-height) so draw() sees it.
            all_rows[[length(all_rows) + 1]] <- data.frame(
              gene            = gene_id,
              x0b             = x0c[i], x1b = x1c[i],
              ymid            = ymid[i],
              dir             = ifelse(strand_g[i] == "-", -1, 1),
              label           = label_g[i],
              color           = color_g[i],
              exon_x0         = NA_real_,
              exon_x1         = NA_real_,
              exon_y0         = NA_real_,
              exon_y1         = NA_real_,
              draw_box        = FALSE,
              strand          = strand_g[i],
              tss_x           = tss_npc[gene_id],
              chevron_npc     = NA_real_,
              use_sep_strands = use_sep_strands,
              label_pos       = lpos,
              label_y         = laby_i,
              pad_bp          = pad_i,
              panel_y0        = pm$inner$y0,
              panel_y1        = pm$inner$y1,
              stringsAsFactors = FALSE,
              row.names  = NULL
            )
            next
          }

          start_sub <- start_sub[draw_mask]
          end_sub   <- end_sub[draw_mask]
          height_scale <- height_scale[draw_mask]

          ux0 <- (start_sub - xpr[1]) / diff(xpr)
          ux1 <- (end_sub   - xpr[1]) / diff(xpr)
          ex0 <- pm$inner$x0 + ux0 * (pm$inner$x1 - pm$inner$x0)
          ex1 <- pm$inner$x0 + ux1 * (pm$inner$x1 - pm$inner$x0)
          half_h <- exon_h / 2 * height_scale
          ey0 <- ymid[i] - half_h
          ey1 <- ymid[i] + half_h

          all_rows[[length(all_rows) + 1]] <- data.frame(
            gene            = gene_id,
            x0b             = x0c[i], x1b = x1c[i],
            ymid            = ymid[i],
            dir             = ifelse(strand_g[i] == "-", -1, 1),
            label           = label_g[i],
            color           = color_g[i],
            exon_x0         = ex0,
            exon_x1         = ex1,
            exon_y0         = ey0,
            exon_y1         = ey1,
            draw_box        = TRUE,
            strand          = strand_g[i],
            tss_x           = tss_npc[gene_id],
            chevron_npc     = NA_real_,
            use_sep_strands = use_sep_strands,
            label_pos       = lpos,
            label_y         = laby_i,
            pad_bp          = pad_i,
            panel_y0        = pm$inner$y0,
            panel_y1        = pm$inner$y1,
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
      }

      self$coordCanvas <- if (length(all_rows) > 0)
        do.call(rbind, all_rows)
      else
        data.frame()
      invisible()
    },

    #' @description Draw backbones, exons/UTRs, directional arrows, and
    #'   labels.
    draw = function() {
      df <- self$coordCanvas
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
        return(invisible())

      by_gene <- split(df, df$gene)
      for (sub in by_gene) {
        col <- sub$color[1]
        dir <- sub$dir[1]
        lbl <- sub$label[1]
        ym  <- sub$ymid[1]
        x0b <- sub$x0b[1]
        x1b <- sub$x1b[1]

        # Glyph border (exon boxes and the "gene" chevron). NA exon_lwd keeps
        # the borderless default; a set exon_lwd strokes exon_border_col
        # (falling back to the gene colour).
        if (is.na(self$exon_lwd)) {
          box_border <- NA; box_lwd <- 0.5
        } else {
          box_lwd    <- self$exon_lwd
          box_border <- if (is.na(self$exon_border_col)) col
                        else self$exon_border_col
        }

        # ── style_type = "exon" (default): backbone + boxes + arrows ──────
        if (self$style_type == "exon") {
          # Backbone line — lty depends on backbone_type.
          lty_val <- switch(self$backbone_type,
            arrow  = 1,
            solid  = 1,
            dashed = 2
          )
          grid::grid.lines(
            x  = grid::unit(c(x0b, x1b), "npc"),
            y  = grid::unit(c(ym,  ym),  "npc"),
            gp = grid::gpar(col = col, lwd = self$backbone_lwd,
                            lty = lty_val, lineend = "butt")
          )

          # Chevron arrows only for backbone_type == "arrow".
          if (self$backbone_type == "arrow") {
            spacing_npc   <- 0.02
            arrow_len_npc <- 0.008   # nonzero so arrow has a defined direction
            # Use a panel-global NPC grid (0.02, 0.04, ...) so all genes share
            # the same arrow phase regardless of where they start.
            all_xs <- seq(spacing_npc, 1 - spacing_npc / 2, by = spacing_npc)
            xs <- all_xs[all_xs > x0b & all_xs < x1b]
            for (x in xs) {
              grid::grid.segments(
                x0 = grid::unit(x - dir * arrow_len_npc, "npc"),
                x1 = grid::unit(x, "npc"),
                y0 = grid::unit(ym, "npc"),
                y1 = grid::unit(ym, "npc"),
                gp = grid::gpar(col = col, lwd = self$backbone_lwd),
                arrow = grid::arrow(type = "open", angle = 45,
                                    length = grid::unit(self$arrow_size, "mm"))
              )
            }
          }

          # Exon/UTR boxes use the shared glyph border computed above.
          for (j in seq_len(nrow(sub))) {
            if (!isTRUE(sub$draw_box[j])) next
            grid::grid.rect(
              x      = grid::unit((sub$exon_x0[j] + sub$exon_x1[j]) / 2, "npc"),
              y      = grid::unit((sub$exon_y0[j] + sub$exon_y1[j]) / 2, "npc"),
              width  = grid::unit(sub$exon_x1[j] - sub$exon_x0[j], "npc"),
              height = grid::unit(sub$exon_y1[j] - sub$exon_y0[j], "npc"),
              gp     = grid::gpar(fill = col, col = box_border, lwd = box_lwd)
            )
          }
        }

        # ── style_type = "gene": single chevron polygon per gene ──────────
        if (self$style_type == "gene") {
          ey0  <- sub$exon_y0[1]
          ey1  <- sub$exon_y1[1]
          chev <- sub$chevron_npc[1]
          if (!is.finite(chev)) chev <- 0
          # Clamp chevron to body width so a tiny gene doesn't get a
          # pointed tip wider than itself.
          chev <- min(chev, max(0, x1b - x0b))

          if (dir > 0) {
            xs <- c(x0b, x1b - chev, x1b,        x1b - chev, x0b)
            ys <- c(ey0, ey0,        ym,         ey1,        ey1)
          } else {
            xs <- c(x0b + chev, x1b, x1b, x0b + chev, x0b)
            ys <- c(ey0,        ey0, ey1, ey1,        ym)
          }
          grid::grid.polygon(
            x  = grid::unit(xs, "npc"),
            y  = grid::unit(ys, "npc"),
            gp = grid::gpar(fill = col, col = box_border, lwd = box_lwd)
          )
        }

        # ── style_type = "point": single TSS dot per gene ─────────────────
        if (self$style_type == "point") {
          tss_x <- sub$tss_x[1]
          if (is.finite(tss_x)) {
            grid::grid.points(
              x    = grid::unit(tss_x, "npc"),
              y    = grid::unit(ym,    "npc"),
              pch  = 16,
              size = grid::unit(2, "mm"),
              gp   = grid::gpar(col = col)
            )
          }
        }

        # TSS start-site marker: up-then-over stem terminated by a solid
        # filled triangle (a grid.polygon, not a grid arrow — grid's
        # `type = "closed"` arrows render as an unfilled outline on many
        # devices). Skipped in "point" mode. Strand-aware: `+` up-right,
        # `-` up-left, with the triangle leading slightly ahead of the arm.
        if (isTRUE(self$show_start) && self$style_type != "point") {
          tss_x <- sub$tss_x[1]
          rise  <- self$tss_rise
          if (is.null(rise)) {
            rise <- (sub$exon_y1[1] - sub$exon_y0[1]) * 0.9
            if (!(is.finite(rise) && rise > 0)) rise <- 0.02
          }
          run     <- self$tss_run
          stem_y1 <- ym + rise
          arm_end <- tss_x + dir * run

          # Vertical stem
          grid::grid.lines(
            x  = grid::unit(c(tss_x, tss_x), "npc"),
            y  = grid::unit(c(ym, stem_y1), "npc"),
            gp = grid::gpar(col = col, lwd = self$tss_stem_lwd)
          )
          # Horizontal arm
          grid::grid.segments(
            x0 = grid::unit(tss_x,   "npc"), x1 = grid::unit(arm_end, "npc"),
            y0 = grid::unit(stem_y1, "npc"), y1 = grid::unit(stem_y1, "npc"),
            gp = grid::gpar(col = col, lwd = self$tss_stem_lwd)
          )
          # Solid filled triangle arrowhead, led slightly ahead of the arm.
          tri_w <- grid::convertWidth(
            grid::unit(self$tss_arrow_size, "mm"), "npc", valueOnly = TRUE)
          tri_h <- grid::convertHeight(
            grid::unit(self$tss_arrow_size, "mm"), "npc", valueOnly = TRUE)
          base_x <- arm_end + dir * tri_w * 0.25   # small forward offset
          apex_x <- base_x  + dir * tri_w
          grid::grid.polygon(
            x  = grid::unit(c(base_x, base_x, apex_x), "npc"),
            y  = grid::unit(c(stem_y1 - tri_h / 2,
                              stem_y1 + tri_h / 2, stem_y1), "npc"),
            gp = grid::gpar(fill = col, col = col)
          )
        }

        gl_aes  <- self$aesthetics[["gene.label"]]
        lbl_col <- (if (!is.null(gl_aes)) gl_aes$color else NULL) %||% col
        lbl_cex <- (if (!is.null(gl_aes)) gl_aes$size  else NULL) %||%
                   self$label_cex
        lbl_hjust <- if (!is.null(gl_aes)) gl_aes$hjust else NULL
        lbl_vjust <- if (!is.null(gl_aes)) gl_aes$vjust else NULL

        lp <- if (!is.null(gl_aes) && !is.null(gl_aes$position))
          .parse_gene_label_pos(gl_aes$position)
        else
          list(x_pos = "strand", y_pos = "center")

        # label_position = "above"/"below" centres the label horizontally over
        # the gene body at the reserved-row y computed in prep(). An explicit
        # gene.label$position still wins (handled by the else branch).
        stack_lbl <- self$label_position %in% c("above", "below") &&
          (is.null(gl_aes) || is.null(gl_aes$position))

        if (stack_lbl) {
          labx <- (x0b + x1b) / 2
          laby <- sub$label_y[1]
          if (!is.finite(laby))
            laby <- if (self$label_position == "above") ym + 0.02 else ym - 0.02
          default_hjust_str <- "center"
          default_vjust_str <- "center"
        } else {
          labx <- switch(lp$x_pos,
            start  = x0b,
            end    = x1b,
            strand = if (dir > 0) x0b - self$label_offset
                     else         x1b + self$label_offset
          )

          default_hjust_str <- switch(lp$x_pos,
            start  = "left",
            end    = "right",
            strand = if (dir > 0) "right" else "left"
          )

          # Resolve y from exon bounds when available
          exon_rows <- sub[!is.na(sub$exon_y0) &
                           (sub$draw_box %in% TRUE), ]
          if (nrow(exon_rows) > 0) {
            exon_top    <- max(exon_rows$exon_y1, na.rm = TRUE)
            exon_bottom <- min(exon_rows$exon_y0, na.rm = TRUE)
          } else {
            exon_top    <- ym + 0.01
            exon_bottom <- ym - 0.01
          }

          laby <- switch(lp$y_pos,
            top    = exon_top    + 0.005,
            bottom = exon_bottom - 0.005,
            center = ym
          )

          default_vjust_str <- switch(lp$y_pos,
            top    = "bottom",
            bottom = "top",
            center = "center"
          )
        }

        just_h <- lbl_hjust %||% default_hjust_str
        just_v <- lbl_vjust %||% default_vjust_str

        grid::grid.text(
          label = lbl,
          x     = grid::unit(labx, "npc"),
          y     = grid::unit(laby, "npc"),
          just  = c(just_h, just_v),
          gp    = grid::gpar(cex = lbl_cex, col = lbl_col)
        )
      }

      # Strand band divider + labels when separate_strands is active.
      if (nrow(df) > 0 && isTRUE(df$use_sep_strands[1])) {
        py0   <- df$panel_y0[1]
        py1   <- df$panel_y1[1]
        mid_y <- (py0 + py1) / 2

        # Dashed divider between + and - bands
        grid::grid.lines(
          x  = grid::unit(c(0, 1), "npc"),
          y  = grid::unit(c(mid_y, mid_y), "npc"),
          gp = grid::gpar(col = "grey60", lwd = 0.5, lty = 2)
        )

        # + label at 0.75 of track height (centre of plus band)
        grid::grid.text(
          label = "+",
          x     = grid::unit(0.01, "npc"),
          y     = grid::unit(py0 + 0.75 * (py1 - py0), "npc"),
          just  = c("left", "center"),
          gp    = grid::gpar(cex = 0.7, col = "grey40", fontface = "bold")
        )

        # − label at 0.25 of track height (centre of minus band)
        grid::grid.text(
          label = "−",
          x     = grid::unit(0.01, "npc"),
          y     = grid::unit(py0 + 0.25 * (py1 - py0), "npc"),
          just  = c("left", "center"),
          gp    = grid::gpar(cex = 0.7, col = "grey40", fontface = "bold")
        )
      }
    }
  )
)

#' Draw gene models
#'
#' Composite element rendering gene models with backbone lines, exon/UTR
#' boxes, directional arrows, and labels. Entirely format-agnostic —
#' supply the relevant column names via `map()`.
#'
#' @param data Optional `GRanges`.
#' @param mapping Optional `SeqMap`. Recognised fields: `group`
#'   (feature-to-gene grouping; features sharing a value are one gene),
#'   `strand` (default `"+"`), `label` (defaults to the group value),
#'   `type` (`"exon"` full-height box, `"UTR"` 80-percent-height box,
#'   anything else no box), and `color` (per-feature color).
#' @param aesthetics Optional `SeqAes`. Supports `color` (default
#'   `"gray30"`), `linewidth`, `alpha`.
#' @param backbone_type Character. One of `"arrow"` (default), `"solid"`, or
#'   `"dashed"`. Controls backbone line style and whether chevron arrows are
#'   drawn. Ignored when `style_type` is `"gene"` or `"point"` (no backbone
#'   exists in those modes).
#' @param show_start Logical. When `TRUE`, a TSS flag arrow is drawn above the
#'   first exon of each gene. Default `FALSE`. Ignored when
#'   `style_type = "point"` (the point itself marks the TSS).
#' @param tss_position Named list overriding per-gene TSS genomic positions.
#'   Keys are gene IDs; values are `c(start, end)` genomic coordinates.
#'   Default `NULL` (auto-detected from first exon by row order).
#' @param separate_strands Logical. When `TRUE`, genes are split into `"+"`
#'   (top) and `"-"` (bottom) sub-bands, each labelled. Silently ignored when
#'   only one strand is present. Default `FALSE`.
#' @param style_type Character. One of `"exon"` (default), `"gene"`, or
#'   `"point"`. Selects the per-gene rendering style: `"exon"` draws the full
#'   backbone with exon/UTR boxes; `"gene"` draws a single chevron-shaped
#'   polygon spanning the gene extent (no exon detail); `"point"` draws a
#'   single filled circle at the TSS. Labels are drawn in all three modes.
#' @param label_position Character. One of `"flush"` (default), `"above"`, or
#'   `"below"`. `"flush"` places the label inline against the gene's 5'/3' end
#'   (current behaviour); `"above"`/`"below"` centre the label horizontally
#'   over the gene body and offset it vertically within the gene's tier row.
#' @param label_row_frac Numeric in `[0, 1)`. Fraction of each tier row
#'   reserved for the label text when `label_position` is `"above"`/`"below"`;
#'   the glyph height shrinks accordingly so labels do not collide with the
#'   next tier. Default `0.35`. Ignored for `"flush"`.
#' @param pack_direction Character. One of `"bottom_to_top"` (default) or
#'   `"top_to_bottom"`. `"bottom_to_top"` draws tier 1 at the track bottom and
#'   stacks upward; `"top_to_bottom"` draws tier 1 at the track top and stacks
#'   downward. Tier *numbers* keep meaning packing order in both cases.
#' @param expanded Logical. When `TRUE`, the greedy overlap packing is skipped
#'   and every gene is given its own tier (one gene per row), ordered by
#'   genomic start (ties broken by label). Ranking is per window (a gene absent
#'   from a window does not reserve a row there); glyph height stays constant
#'   across windows. Default `FALSE`.
#' @param label_pad Numeric. In `label_position = "flush"`, extra bp padding
#'   added beyond the measured rendered label width when packing tiers, so
#'   genes whose labels would collide land on separate tiers. Default `50000`.
#' @param backbone_lwd Numeric. Backbone (and chevron) line width. Default `1`.
#' @param exon_lwd Numeric or `NA`. Exon/UTR box border width. `NA` (default)
#'   draws unbordered boxes as before; a number strokes the border.
#' @param exon_border_col Border colour used only when `exon_lwd` is set.
#'   Default `NA` (falls back to the gene colour).
#' @param arrow_size Numeric. Chevron arrowhead length in mm. Default `1.5`.
#' @param tss_arrow_size Numeric. TSS marker triangle size in mm. Default `3`.
#' @param tss_stem_lwd Numeric. TSS marker stem/arm line width. Default `1.2`.
#' @param tss_rise Numeric or `NULL`. TSS stem rise in npc; `NULL` (default)
#'   derives it from the gene box height.
#' @param tss_run Numeric. TSS arm run in npc. Default `0.018`.
#' @param oob Character. Out-of-bounds gene handling: `"truncate"` (default,
#'   clip straddling genes to the window), `"keep"` (draw out-of-window
#'   portions — the prior behaviour), or `"remove"` (drop any gene not entirely
#'   within the window).
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqGeneR6` instance.
#' @examples
#' seq_gene(map(group = gene_id, strand = strand, label = gene_name))
#' @export
seq_gene <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                     backbone_type    = "arrow",
                     show_start       = FALSE,
                     tss_position     = NULL,
                     separate_strands = FALSE,
                     style_type       = c("exon", "gene", "point"),
                     label_position   = c("flush", "above", "below"),
                     label_row_frac   = 0.35,
                     pack_direction   = c("bottom_to_top", "top_to_bottom"),
                     expanded         = FALSE,
                     label_pad        = 50000,
                     backbone_lwd     = 1,
                     exon_lwd         = NA,
                     exon_border_col  = NA,
                     arrow_size       = 1.5,
                     tss_arrow_size   = 3,
                     tss_stem_lwd     = 1.2,
                     tss_rise         = NULL,
                     tss_run          = 0.018,
                     oob              = c("truncate", "keep", "remove"),
                     ...) {
  SeqGeneR6$new(
    data             = data,
    mapping          = mapping,
    aesthetics       = aesthetics,
    backbone_type    = backbone_type,
    show_start       = show_start,
    tss_position     = tss_position,
    separate_strands = separate_strands,
    style_type       = style_type,
    label_position   = label_position,
    label_row_frac   = label_row_frac,
    pack_direction   = pack_direction,
    expanded         = expanded,
    label_pad        = label_pad,
    backbone_lwd     = backbone_lwd,
    exon_lwd         = exon_lwd,
    exon_border_col  = exon_border_col,
    arrow_size       = arrow_size,
    tss_arrow_size   = tss_arrow_size,
    tss_stem_lwd     = tss_stem_lwd,
    tss_rise         = tss_rise,
    tss_run          = tss_run,
    oob              = oob,
    ...
  )
}
