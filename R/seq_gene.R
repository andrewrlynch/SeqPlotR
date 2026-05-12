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
    #' @field label_pad Minimum padding (in bp) around each gene's extent to
    #'   keep labels from overlapping neighbouring genes.
    label_pad     = 50000,
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
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          backbone_type    = "arrow",
                          show_start       = FALSE,
                          tss_position     = NULL,
                          separate_strands = FALSE,
                          style_type       = c("exon", "gene", "point"),
                          ...) {
      super$initialize(data, mapping, aesthetics)
      backbone_type <- match.arg(backbone_type, c("arrow", "solid", "dashed"))
      style_type    <- match.arg(style_type)
      self$backbone_type   <- backbone_type
      self$show_start      <- isTRUE(show_start)
      self$tss_position    <- tss_position
      self$separate_strands <- isTRUE(separate_strands)
      self$style_type      <- style_type
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

      for (w in seq_len(nWin)) {
        pm  <- layout_track[[w]]
        win <- track_windows[w]

        hits <- GenomicRanges::findOverlaps(eff_data, win)
        if (length(hits) == 0) next
        keep <- unique(S4Vectors::queryHits(hits))

        exons_start <- BiocGenerics::start(eff_data)[keep]
        exons_end   <- BiocGenerics::end(eff_data)[keep]
        gid_all     <- group_vals[keep]
        strand_all  <- strand_vals[keep]
        label_all   <- label_vals[keep]
        type_all    <- type_vals[keep]
        color_all   <- color_vals[keep]

        gids <- unique(gid_all)
        if (length(gids) == 0) next

        gene_start <- integer(length(gids))
        gene_end   <- integer(length(gids))
        strand_g   <- character(length(gids))
        label_g    <- character(length(gids))
        color_g    <- character(length(gids))
        for (i in seq_along(gids)) {
          m <- gid_all == gids[i]
          gene_start[i] <- min(exons_start[m])
          gene_end[i]   <- max(exons_end[m])
          strand_g[i]   <- strand_all[m][1]
          label_g[i]    <- label_all[m][1]
          color_g[i]    <- color_all[m][1]
        }

        # Pack tiers by raw genomic extent — label space is not part of
        # the overlap test (labels can extend outside the gene body without
        # forcing a new tier).
        p0 <- gene_start
        p1 <- gene_end

        # Determine whether strand separation is active for this window.
        use_sep_strands <- isTRUE(self$separate_strands) &&
          length(unique(strand_all[strand_all %in% c("+", "-")])) > 1L

        track_h <- pm$inner$y1 - pm$inner$y0

        if (use_sep_strands) {
          # Plus-strand genes occupy the top half, minus-strand the bottom half.
          plus_y0  <- pm$inner$y0 + track_h / 2
          plus_y1  <- pm$inner$y1
          minus_y0 <- pm$inner$y0
          minus_y1 <- pm$inner$y0 + track_h / 2
        } else {
          plus_y0  <- pm$inner$y0
          plus_y1  <- pm$inner$y1
          minus_y0 <- pm$inner$y0
          minus_y1 <- pm$inner$y1
        }

        # Tier stacking — separate per band when use_sep_strands.
        tiers <- integer(length(gids))
        if (use_sep_strands) {
          plus_gids  <- gids[strand_g == "+"]
          minus_gids <- gids[strand_g == "-"]

          for (band_gids in list(plus_gids, minus_gids)) {
            if (length(band_gids) == 0) next
            idx       <- match(band_gids, gids)
            ends_last <- numeric(0)
            for (i in order(p0[idx])) {
              gi <- idx[i]
              s <- p0[gi]; e <- p1[gi]
              w2 <- which(ends_last < s)
              if (length(w2) == 0) {
                tiers[gi]  <- length(ends_last) + 1
                ends_last  <- c(ends_last, e)
              } else {
                tiers[gi]        <- w2[1]
                ends_last[w2[1]] <- e
              }
            }
          }
          ntiers_plus  <- if (length(plus_gids))
            max(tiers[match(plus_gids,  gids)]) else 1L
          ntiers_minus <- if (length(minus_gids))
            max(tiers[match(minus_gids, gids)]) else 1L
        } else {
          # Original single-stack tier assignment.
          ord       <- order(p0)
          ends_last <- numeric(0)
          for (i in ord) {
            s <- p0[i]; e <- p1[i]
            w2 <- which(ends_last < s)
            if (length(w2) == 0) {
              tiers[i]   <- length(ends_last) + 1
              ends_last  <- c(ends_last, e)
            } else {
              tiers[i]            <- w2[1]
              ends_last[w2[1]]    <- e
            }
          }
          ntiers_plus  <- max(tiers)
          ntiers_minus <- max(tiers)
        }

        xpr <- pm$xplot_range %||% pm$xscale
        u0 <- (gene_start - xpr[1]) / diff(xpr)
        u1 <- (gene_end   - xpr[1]) / diff(xpr)
        x0c <- pm$inner$x0 + u0 * (pm$inner$x1 - pm$inner$x0)
        x1c <- pm$inner$x0 + u1 * (pm$inner$x1 - pm$inner$x0)

        # Per-gene ymid using correct band.
        ymid <- numeric(length(gids))
        for (i in seq_along(gids)) {
          if (use_sep_strands && strand_g[i] == "+") {
            band_y0 <- plus_y0;  band_y1 <- plus_y1;  nt <- ntiers_plus
          } else if (use_sep_strands && strand_g[i] == "-") {
            band_y0 <- minus_y0; band_y1 <- minus_y1; nt <- ntiers_minus
          } else {
            band_y0 <- pm$inner$y0; band_y1 <- pm$inner$y1
            nt <- max(tiers)
          }
          band_h_i <- band_y1 - band_y0
          row_h_i  <- band_h_i / nt
          exon_h_i <- row_h_i * self$exon_height
          ymid[i]  <- band_y0 + (tiers[i] - 1) * row_h_i + exon_h_i / 2
        }

        # TSS detection for show_start.
        tss_genomic <- setNames(numeric(length(gids)), gids)
        for (i in seq_along(gids)) {
          gid <- gids[i]
          m   <- gid_all == gid
          if (!is.null(self$tss_position) &&
              gid %in% names(self$tss_position)) {
            tp <- self$tss_position[[gid]]
            tss_genomic[gid] <- tp[1]
          } else {
            # For + strand: TSS is at the leftmost start (min genomic).
            # For - strand: TSS is at the rightmost end (max genomic).
            if (strand_g[i] == "-") {
              tss_genomic[gid] <- max(exons_end[m])
            } else {
              tss_genomic[gid] <- min(exons_start[m])
            }
          }
        }

        # Convert TSS genomic positions to NPC.
        tss_npc <- (tss_genomic - xpr[1]) / diff(xpr) *
                   (pm$inner$x1 - pm$inner$x0) + pm$inner$x0

        for (i in seq_along(gids)) {
          gene_id <- gids[i]
          m       <- gid_all == gene_id
          type_sub  <- type_all[m]
          start_sub <- exons_start[m]
          end_sub   <- exons_end[m]

          # Compute exon_h for this gene from its band.
          if (use_sep_strands && strand_g[i] == "+") {
            nt_i <- ntiers_plus
            bh_i <- plus_y1 - plus_y0
          } else if (use_sep_strands && strand_g[i] == "-") {
            nt_i <- ntiers_minus
            bh_i <- minus_y1 - minus_y0
          } else {
            nt_i <- max(tiers)
            bh_i <- pm$inner$y1 - pm$inner$y0
          }
          exon_h <- (bh_i / nt_i) * self$exon_height

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
            gp = grid::gpar(col = col, lwd = 1, lty = lty_val, lineend = "butt")
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
                gp = grid::gpar(col = col, lwd = 1),
                arrow = grid::arrow(type = "open", angle = 45,
                                    length = grid::unit(1.5, "mm"))
              )
            }
          }

          for (j in seq_len(nrow(sub))) {
            if (!isTRUE(sub$draw_box[j])) next
            grid::grid.rect(
              x      = grid::unit((sub$exon_x0[j] + sub$exon_x1[j]) / 2, "npc"),
              y      = grid::unit((sub$exon_y0[j] + sub$exon_y1[j]) / 2, "npc"),
              width  = grid::unit(sub$exon_x1[j] - sub$exon_x0[j], "npc"),
              height = grid::unit(sub$exon_y1[j] - sub$exon_y0[j], "npc"),
              gp     = grid::gpar(fill = col, col = NA)
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
            gp = grid::gpar(fill = col, col = NA)
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

        # TSS start-site flag: upward L-shaped open arrow.
        # Skipped in "point" mode — the point itself marks the TSS.
        if (isTRUE(self$show_start) && self$style_type != "point") {
          tss_x   <- sub$tss_x[1]
          # Use exon/gene box height when present; fall back when NA.
          arrow_h <- (sub$exon_y1[1] - sub$exon_y0[1]) * 0.9
          arrow_h <- if (is.finite(arrow_h) && arrow_h > 0) arrow_h else 0.02

          stem_y0 <- ym
          stem_y1 <- ym + arrow_h

          # Vertical stem
          grid::grid.lines(
            x  = grid::unit(c(tss_x, tss_x), "npc"),
            y  = grid::unit(c(stem_y0, stem_y1), "npc"),
            gp = grid::gpar(col = col, lwd = 1.2)
          )
          # Horizontal arm with open arrowhead at the tip
          flag_len <- 0.018
          flag_x1  <- tss_x + dir * flag_len
          grid::grid.segments(
            x0 = grid::unit(tss_x,   "npc"),
            x1 = grid::unit(flag_x1, "npc"),
            y0 = grid::unit(stem_y1, "npc"),
            y1 = grid::unit(stem_y1, "npc"),
            gp = grid::gpar(col = col, lwd = 1.2),
            arrow = grid::arrow(type = "closed", angle = 25,
                                length = grid::unit(3, "mm"),
                                ends   = "last")
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
    ...
  )
}
