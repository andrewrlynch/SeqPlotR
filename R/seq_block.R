# Block palettes and colours ----

# reference hues by default, a softer set for "soft", any other vector as given
.block_palette <- function(palette = NULL) {
  if (is.null(palette))
    return(c("#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948", "#9C755F", "#E15759"))
  if (identical(palette, "soft"))
    return(c("#C49A6C", "#8DB0D6", "#9CC79A", "#D6A0C4", "#E2C46F", "#A9A3D6"))
  palette
}

# block fill colours: key colour or palette entry, tinted toward white, then block_colors overrides
.block_fill_colors <- function(key, tint = 1L, label = NULL, palette = NULL,
                               tints = c(0.15, 0.5, 0.3, 0.65, 0.4), block_colors = NULL) {
  n <- length(key)
  if (n == 0L) return(character(0))
  key <- as.character(key)
  pal <- .block_palette(palette)
  base <- if (!anyNA(key) && .looks_like_color(key)) {
    key
  } else if (!is.null(names(pal))) {
    unname(pal[key])
  } else {
    u <- unique(key)
    pal[(match(key, u) - 1L) %% length(pal) + 1L]
  }
  base[is.na(base)] <- "grey80"
  tint <- rep_len(suppressWarnings(as.integer(tint)), n)
  tint[is.na(tint)] <- 1L
  t <- if (length(tints)) tints[(tint - 1L) %% length(tints) + 1L] else rep(0, n)
  m <- grDevices::col2rgb(base) / 255
  out <- grDevices::rgb(m[1, ] + (1 - m[1, ]) * t, m[2, ] + (1 - m[2, ]) * t,
                        m[3, ] + (1 - m[3, ]) * t)
  if (length(block_colors) && !is.null(names(block_colors))) {
    bc <- unlist(block_colors)
    by_label <- if (is.null(label)) rep(NA_character_, n) else unname(bc[as.character(label)])
    by_key <- unname(bc[key])
    out <- ifelse(!is.na(by_label), by_label, ifelse(!is.na(by_key), by_key, out))
  }
  out
}

# cytoband stain colours
.block_stain_colors <- function(stain) {
  pal <- c(gneg = "white", gpos25 = "grey80", gpos33 = "grey72", gpos50 = "grey60",
           gpos66 = "grey45", gpos75 = "grey40", gpos100 = "grey10", gvar = "#C7E9F1",
           stalk = "grey70", acen = "#D33A2C")
  out <- unname(pal[as.character(stain)])
  out[is.na(out)] <- "grey80"
  out
}

# SeqBlockR6 class ----

#' SeqBlock R6 class
#'
#' Internal R6 generator backing [seq_block()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqBlockR6 <- R6::R6Class("SeqBlock",
  inherit = SeqElementR6,
  public = list(
    #' @field fill `"solid"` or `"cytoband"`.
    fill = "solid",
    #' @field cytobands Cytoband table (`GRanges` or data frame) or `NULL`.
    cytobands = NULL,
    #' @field palette Palette spec passed to the fill colour helper.
    palette = NULL,
    #' @field tints Tint fractions mixed toward white.
    tints = NULL,
    #' @field block_colors Named colour overrides by label or key.
    block_colors = NULL,
    #' @field shape `"rect"` or `"arrow"`.
    shape = "rect",
    #' @field height Block height as a fraction of the panel height.
    height = 0.5,
    #' @field label_side `"below"` or `"above"`.
    label_side = "below",
    #' @field flip_style `"rotate"`, `"overline"` or `"minus"`.
    flip_style = "rotate",
    #' @field annotations Optional annotation data frame.
    annotations = NULL,

    #' @description Construct a SeqBlockR6.
    #' @param data `GRanges` of block positions in track coordinates.
    #' @param mapping Optional `SeqMap`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param fill,cytobands,palette,tints,block_colors,shape,height See [seq_block()].
    #' @param label_side,flip_style,annotations See [seq_block()].
    #' @param legend Optional legend.
    #' @param show_legend Logical.
    #' @param ... Unused.
    initialize = function(data = NULL, mapping = NULL, aesthetics = aes(),
                          fill = c("solid", "cytoband"), cytobands = NULL,
                          palette = NULL, tints = c(0.15, 0.5, 0.3, 0.65, 0.4),
                          block_colors = NULL, shape = c("rect", "arrow"),
                          height = 0.5, label_side = c("below", "above"),
                          flip_style = c("rotate", "overline", "minus"),
                          annotations = NULL, legend = NULL, show_legend = TRUE, ...) {
      super$initialize(data, mapping, aesthetics, legend = legend,
                       show_legend = show_legend)
      if (!is.null(self$data) && !inherits(self$data, "GRanges"))
        stop("seq_block requires GRanges data.", call. = FALSE)
      if (!is.null(annotations) &&
          !(is.data.frame(annotations) && all(c("chrom", "start", "label") %in% names(annotations))))
        stop("annotations must be a data.frame with columns chrom, start and label.",
             call. = FALSE)
      self$fill         <- match.arg(fill)
      self$cytobands    <- cytobands
      self$palette      <- palette
      self$tints        <- tints
      self$block_colors <- block_colors
      self$shape        <- match.arg(shape)
      self$height       <- height
      self$label_side   <- match.arg(label_side)
      self$flip_style   <- match.arg(flip_style)
      self$annotations  <- annotations
    },

    #' @description Resolve mappings, assign blocks to windows and compute npc
    #'   geometry for fills, bands, centromeres, labels and annotations.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    #' @return The `coordCanvas` list (one entry per window), invisibly.
    prep = function(layout_track, track_windows) {
      self$resolve(track_data    = layout_track[[1]]$track_data,
                   track_mapping = layout_track[[1]]$track_mapping)
      self$coordCanvas <- vector("list", length(track_windows))
      d <- self$resolved$.data
      if (is.null(d) || length(d) == 0L) return(invisible(self$coordCanvas))
      if (!inherits(d, "GRanges")) stop("seq_block requires GRanges data.", call. = FALSE)

      # resolved fields with defaults ----
      r <- self$resolved
      n <- length(d)
      vec <- function(v, default) {
        if (is.null(v)) v <- default
        if (length(v) == 1L) rep(v, n) else v
      }
      lgl <- function(v) { v <- as.logical(v); !is.na(v) & v }
      xs    <- as.numeric(vec(r$x, BiocGenerics::start(d)))
      xe    <- as.numeric(vec(r$x_end, BiocGenerics::end(d)))
      chrom <- as.character(vec(r$chrom, as.character(GenomicRanges::seqnames(d))))
      rs    <- as.numeric(vec(r$ref_start, xs))
      re    <- as.numeric(vec(r$ref_end, xe))
      flip  <- lgl(vec(r$flipped, FALSE))
      lab   <- as.character(vec(r$label, NA_character_))
      key   <- vec(r$key, chrom)
      tint  <- vec(r$tint, 1L)
      lin   <- lgl(vec(r$label_inside, FALSE))
      dash  <- lgl(vec(r$dashed, FALSE))
      cs    <- as.numeric(vec(r$centromere_start, NA_real_))
      ce    <- as.numeric(vec(r$centromere_end, NA_real_))
      fills <- .block_fill_colors(key, tint, lab, self$palette, self$tints, self$block_colors)

      a     <- self$aesthetics
      lsize <- a$label.size %||% 0.7
      ti    <- layout_track[[1]]$track_inner
      tw    <- ti$x1 - ti$x0
      hfrac <- a$arrow.head %||% 0.35
      hmax  <- a$arrow.head_max %||% 0.012
      toggle <- function(ch) ifelse(startsWith(ch, "chr"), sub("^chr", "", ch), paste0("chr", ch))

      # cytoband table as half-open chrom/start/end/stain ----
      bt <- NULL
      if (self$fill == "cytoband" && self$shape == "rect") {
        cb <- private$bands %||% (self$cytobands %||% load_cytobands())
        if (inherits(cb, "GRanges")) {
          mc <- S4Vectors::mcols(cb)
          st <- mc$stain %||% mc$gieStain
          bt <- data.frame(chrom = as.character(GenomicRanges::seqnames(cb)),
                           start = BiocGenerics::start(cb) - 1, end = BiocGenerics::end(cb),
                           stain = if (is.null(st)) NA_character_ else as.character(st),
                           stringsAsFactors = FALSE)
        } else if (is.data.frame(cb)) {
          nm <- names(cb)
          pick <- function(a1, a2) if (a1 %in% nm) a1 else if (a2 %in% nm) a2 else NA_character_
          cols <- c(pick("chrom", "chrom"), pick("start", "chromStart"),
                    pick("end", "chromEnd"), pick("stain", "gieStain"))
          if (anyNA(cols))
            stop("cytobands must have columns chrom, start, end and stain.", call. = FALSE)
          bt <- data.frame(chrom = as.character(cb[[cols[1]]]), start = as.numeric(cb[[cols[2]]]),
                           end = as.numeric(cb[[cols[3]]]), stain = as.character(cb[[cols[4]]]),
                           stringsAsFactors = FALSE)
        } else {
          stop("cytobands must be a GRanges or a data.frame.", call. = FALSE)
        }
        private$bands <- bt
      }

      # annotation columns with defaults ----
      an <- self$annotations
      if (!is.null(an) && nrow(an) > 0L) {
        na_ <- nrow(an)
        opt <- function(nm, default) if (nm %in% names(an)) an[[nm]] else rep(default, na_)
        an <- data.frame(chrom = as.character(an$chrom), start = as.numeric(an$start),
                         end = as.numeric(opt("end", NA_real_)), label = as.character(an$label),
                         side = as.character(opt("side", "above")),
                         color = as.character(opt("color", "#2336B8")),
                         italic = as.logical(opt("italic", TRUE)), dot = as.logical(opt("dot", TRUE)),
                         strand = as.character(opt("strand", NA_character_)),
                         stringsAsFactors = FALSE)
        an$end[is.na(an$end)] <- an$start[is.na(an$end)]
        an$side[is.na(an$side) | !an$side %in% c("above", "below")] <- "above"
        an$color[is.na(an$color)] <- "#2336B8"
        an$italic[is.na(an$italic)] <- TRUE
        an$dot[is.na(an$dot)] <- TRUE
      } else {
        an <- NULL
      }

      ov <- GenomicRanges::findOverlaps(d, track_windows, select = "first")
      missing_bands <- character(0)
      bind <- function(l) if (length(l)) do.call(rbind, l) else data.frame()

      for (w in sort(unique(ov[!is.na(ov)]))) {
        pm  <- layout_track[[w]]
        inn <- pm$inner
        xpr <- pm$xplot_range %||% pm$xscale
        npc <- function(v) pmin(pmax(inn$x0 + (v - xpr[1]) / diff(xpr) * (inn$x1 - inn$x0), inn$x0), inn$x1)
        yc <- (inn$y0 + inn$y1) / 2
        hh <- self$height * (inn$y1 - inn$y0) / 2
        y0 <- yc - hh
        y1 <- yc + hh
        L_blocks <- L_bands <- L_arms <- L_bow <- L_arrows <- L_labels <- L_ann <- list()

        for (i in which(ov == w)) {
          x0n <- min(npc(xs[i]), npc(xe[i]))
          x1n <- max(npc(xs[i]), npc(xe[i]))
          fi  <- flip[i]
          f <- function(g) {
            if (re[i] == rs[i]) return(rep((x0n + x1n) / 2, length(g)))
            u <- (g - rs[i]) / (re[i] - rs[i])
            if (fi) x1n - u * (x1n - x0n) else x0n + u * (x1n - x0n)
          }
          solid <- TRUE
          hd <- NA_real_
          if (self$shape == "arrow") {
            hd <- min(hfrac * (x1n - x0n), hmax * tw)
            vx <- if (fi) c(x1n, x0n + hd, x0n, x0n + hd, x1n) else c(x0n, x1n - hd, x1n, x1n - hd, x0n)
            L_arrows[[length(L_arrows) + 1L]] <- data.frame(i = i, k = 1:5, x = vx, y = c(y0, y0, yc, y1, y1))
          } else {
            # arms around a centromere, with its bowtie
            cen <- !is.na(cs[i]) && !is.na(ce[i]) && cs[i] < ce[i] &&
              cs[i] >= min(rs[i], re[i]) && ce[i] <= max(rs[i], re[i])
            arm_ref <- if (cen) Filter(function(v) v[1] < v[2], list(c(rs[i], cs[i]), c(ce[i], re[i])))
                       else list(c(rs[i], re[i]))
            for (ar in arm_ref)
              L_arms[[length(L_arms) + 1L]] <- data.frame(i = i, x0 = min(f(ar)), x1 = max(f(ar)))
            if (cen) {
              xm <- f((cs[i] + ce[i]) / 2)
              L_bow[[length(L_bow) + 1L]] <- data.frame(
                i = i, part = 1:2, xa = c(f(cs[i]), f(ce[i])), xb = c(f(ce[i]), f(cs[i])),
                xm = xm, y0 = y0, y1 = y1, yc = yc)
            }
            if (!is.null(bt)) {
              bb <- bt[bt$chrom == chrom[i], , drop = FALSE]
              if (nrow(bb) == 0L) bb <- bt[bt$chrom == toggle(chrom[i]), , drop = FALSE]
              got <- 0L
              for (ar in arm_ref) {
                hit <- bb[bb$end > ar[1] & bb$start < ar[2], , drop = FALSE]
                if (nrow(hit) == 0L) next
                ga <- pmax(hit$start, ar[1])
                gb <- pmin(hit$end, ar[2])
                L_bands[[length(L_bands) + 1L]] <- data.frame(
                  i = i, x0 = pmin(f(ga), f(gb)), x1 = pmax(f(ga), f(gb)),
                  fill = .block_stain_colors(hit$stain), stringsAsFactors = FALSE)
                got <- got + nrow(hit)
              }
              solid <- got == 0L
              if (solid) missing_bands <- c(missing_bands, chrom[i])
            }
          }
          L_blocks[[length(L_blocks) + 1L]] <- data.frame(
            i = i, x0 = x0n, x1 = x1n, y0 = y0, y1 = y1,
            fill = if (solid) fills[i] else "white", flipped = fi, dashed = dash[i],
            shape = self$shape, head = hd, solid = solid, stringsAsFactors = FALSE)

          # label text, flip marking and fit ----
          if (!is.na(lab[i]) && nzchar(lab[i])) {
            txt <- if (self$flip_style == "minus" && fi) paste0("-", lab[i]) else lab[i]
            wl <- grid::convertWidth(grid::grobWidth(grid::textGrob(txt, gp = grid::gpar(cex = lsize))),
                                     "npc", valueOnly = TRUE)
            L_labels[[length(L_labels) + 1L]] <- data.frame(
              i = i, label = txt, x = (x0n + x1n) / 2, wl = wl,
              inside = lin[i] && 1.15 * wl <= (x1n - x0n),
              rot = if (self$flip_style == "rotate" && fi) 180 else 0,
              overline = self$flip_style == "overline" && fi, stringsAsFactors = FALSE)
          }

          # annotations mapped through the block ----
          if (!is.null(an)) {
            mid <- (an$start + an$end) / 2
            hit <- which((an$chrom == chrom[i] | an$chrom == toggle(chrom[i])) &
                           mid >= min(rs[i], re[i]) & mid <= max(rs[i], re[i]))
            for (j in hit) {
              sg <- if (an$side[j] == "above") 1 else -1
              st <- an$strand[j]
              ax <- if (!is.na(st) && st == "+") f(c(an$start[j], an$end[j]))
                    else if (!is.na(st) && st == "-") f(c(an$end[j], an$start[j]))
                    else c(NA_real_, NA_real_)
              L_ann[[length(L_ann) + 1L]] <- data.frame(
                x = f(mid[j]), y = if (sg > 0) y1 else y0, dy = sg * 0.5, label = an$label[j],
                color = an$color[j], face = if (isTRUE(an$italic[j])) "italic" else "plain",
                dot = isTRUE(an$dot[j]), ax0 = ax[1], ax1 = ax[2], side = an$side[j],
                stringsAsFactors = FALSE)
            }
          }
        }

        # stagger outside labels over three levels ----
        labels <- bind(L_labels)
        if (nrow(labels)) {
          labels$level <- 0L
          out <- which(!labels$inside)
          last <- rep(-Inf, 3L)
          for (k in out[order(labels$x[out])]) {
            hw <- labels$wl[k] / 2
            lv <- which(last < labels$x[k] - hw - 0.004)[1]
            if (is.na(lv)) lv <- which.min(last)
            labels$level[k] <- lv - 1L
            last[lv] <- labels$x[k] + hw
          }
          above <- self$label_side == "above"
          sg <- if (above) 1 else -1
          labels$y <- ifelse(labels$inside, yc, if (above) y1 else y0)
          labels$dy <- ifelse(labels$inside, 0, sg * (0.3 + 0.95 * labels$level))
          vj <- if (above) 0 else 1
          labels$vjust <- ifelse(labels$inside, 0.5, ifelse(labels$rot == 180, 1 - vj, vj))
          labels <- labels[, c("i", "label", "x", "y", "dy", "level", "inside", "rot", "vjust", "overline", "wl")]
        }

        self$coordCanvas[[w]] <- list(
          blocks = bind(L_blocks), bands = bind(L_bands), arms = bind(L_arms),
          bowties = bind(L_bow), arrows = bind(L_arrows), labels = labels,
          annotations = bind(L_ann), track_inner_width = tw
        )
      }
      if (length(missing_bands))
        message("no cytobands for ", paste(unique(missing_bands), collapse = ", "),
                "; drawing solid blocks")
      invisible(self$coordCanvas)
    },

    #' @description Draw fills, bands, centromeres, fade, outlines, labels and
    #'   annotations with grid.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())
      a <- self$aesthetics
      u <- function(v) grid::unit(v, "npc")
      lsize <- a$label.size %||% 0.7
      lcol  <- a$label.col %||% "black"
      fade  <- a$fade %||% 0
      ocol  <- a$outline.col %||% "grey25"
      olwd  <- a$outline.lwd %||% 0.7
      for (cc in self$coordCanvas) {
        if (is.null(cc)) next
        b <- cc$blocks
        bi <- function(i) match(i, b$i)

        # fills ----
        am <- cc$arms
        if (nrow(am)) {
          k <- bi(am$i)
          grid::grid.rect(u(am$x0), u(b$y0[k]), u(am$x1 - am$x0), u(b$y1[k] - b$y0[k]),
                          just = c("left", "bottom"), name = "seqblock.fill",
                          gp = grid::gpar(fill = b$fill[k], col = NA))
        }
        ar <- cc$arrows
        if (nrow(ar)) {
          ids <- sort(unique(ar$i))
          grid::grid.polygon(u(ar$x), u(ar$y), id = ar$i, name = "seqblock.fill",
                             gp = grid::gpar(fill = b$fill[bi(ids)], col = NA))
        }

        # cytobands and centromere bowties ----
        bd <- cc$bands
        if (nrow(bd)) {
          k <- bi(bd$i)
          grid::grid.rect(u(bd$x0), u(b$y0[k]), u(bd$x1 - bd$x0), u(b$y1[k] - b$y0[k]),
                          just = c("left", "bottom"), name = "seqblock.band",
                          gp = grid::gpar(fill = bd$fill, col = NA))
        }
        bw <- cc$bowties
        if (nrow(bw)) {
          m <- nrow(bw)
          grid::grid.polygon(u(as.vector(rbind(bw$xa, bw$xa, bw$xm))),
                             u(as.vector(rbind(bw$y1, bw$y0, bw$yc))),
                             id = rep(seq_len(m), each = 3L), name = "seqblock.centromere",
                             gp = grid::gpar(fill = a$centromere.fill %||% "#D33A2C", col = NA))
        }

        # fade overlay and outlines ----
        if (fade > 0 && nrow(b))
          grid::grid.rect(u(b$x0), u(b$y0), u(b$x1 - b$x0), u(b$y1 - b$y0),
                          just = c("left", "bottom"), name = "seqblock.fade",
                          gp = grid::gpar(fill = "white", col = NA, alpha = fade))
        if (nrow(am)) {
          k <- bi(am$i)
          grid::grid.rect(u(am$x0), u(b$y0[k]), u(am$x1 - am$x0), u(b$y1[k] - b$y0[k]),
                          just = c("left", "bottom"), name = "seqblock.outline",
                          gp = grid::gpar(fill = NA, col = ocol, lwd = olwd,
                                          lty = ifelse(b$dashed[k], "22", "solid")))
        }
        if (nrow(ar)) {
          dids <- intersect(sort(unique(ar$i)), b$i[b$dashed])
          if (length(dids)) {
            s <- ar[ar$i %in% dids, ]
            grid::grid.polygon(u(s$x), u(s$y), id = s$i, name = "seqblock.outline",
                               gp = grid::gpar(fill = NA, col = ocol, lwd = olwd, lty = "22"))
          }
        }

        # labels ----
        lb <- cc$labels
        for (j in seq_len(nrow(lb))) {
          y <- u(lb$y[j]) + grid::unit(lb$dy[j], "lines")
          gp <- grid::gpar(cex = lsize, col = lcol)
          grid::grid.text(lb$label[j], x = u(lb$x[j]), y = y, rot = lb$rot[j],
                          hjust = 0.5, vjust = lb$vjust[j], name = "seqblock.label", gp = gp)
          if (isTRUE(lb$overline[j])) {
            tg <- grid::textGrob(lb$label[j], gp = gp)
            top <- y + (1 - lb$vjust[j]) * grid::grobHeight(tg) + grid::unit(0.12, "lines")
            hw <- 0.5 * grid::grobWidth(tg)
            grid::grid.segments(u(lb$x[j]) - hw, top, u(lb$x[j]) + hw, top,
                                name = "seqblock.overline",
                                gp = grid::gpar(col = lcol, lwd = 0.8))
          }
        }

        # annotations ----
        an <- cc$annotations
        for (j in seq_len(nrow(an))) {
          sg <- if (an$side[j] == "above") 1 else -1
          gp <- grid::gpar(col = an$color[j])
          grid::grid.text(an$label[j], x = u(an$x[j]),
                          y = u(an$y[j]) + grid::unit(an$dy[j], "lines"),
                          vjust = if (sg > 0) 0 else 1, name = "seqblock.annotation",
                          gp = grid::gpar(col = an$color[j], fontface = an$face[j], cex = lsize))
          if (isTRUE(an$dot[j]))
            grid::grid.points(u(an$x[j]), u(an$y[j]) + grid::unit(sg * 0.2, "lines"),
                              pch = 16, size = grid::unit(0.35, "char"),
                              name = "seqblock.annotation", gp = gp)
          if (!is.na(an$ax0[j]))
            grid::grid.segments(u(an$ax0[j]), u(an$y[j]) + grid::unit(sg * 1.3, "lines"),
                                u(an$ax1[j]), u(an$y[j]) + grid::unit(sg * 1.3, "lines"),
                                arrow = grid::arrow(length = grid::unit(0.07, "in"), type = "closed"),
                                name = "seqblock.annotation",
                                gp = grid::gpar(col = an$color[j], fill = an$color[j]))
        }
      }
      invisible()
    }
  ),
  private = list(bands = NULL)
)

# seq_block() constructor ----

#' Ideogram-slice blocks
#'
#' Draws one block per range: a rectangle (or an arrow strip) spanning the
#' exact track positions `[x, x_end]`, filled either with a solid colour or
#' with the cytobands of a reference interval `chrom:ref_start-ref_end` mapped
#' into the block, mirrored when the block is flipped. A centromere interval
#' is drawn as a red bowtie with the outline on the arms only. Labels sit
#' inside the block when requested and they fit, otherwise outside (above or
#' below) with a three-level stagger. Useful for synteny, rearrangement and
#' assembly figures where pieces of reference chromosomes are laid out along
#' new coordinates.
#'
#' @param data `GRanges` whose seqnames are the track windows the blocks sit
#'   in and whose ranges are the block positions in track coordinates.
#' @param mapping Optional `SeqMap`. Fields (all optional): `x`, `x_end`
#'   (exact numeric block edges; default `start`, `end`); `chrom`,
#'   `ref_start`, `ref_end` (reference interval drawn inside the block,
#'   half-open; default `seqnames`, `x`, `x_end`); `flipped` (mirror the
#'   reference interval); `label`; `key` (fill key or a colour; default
#'   `chrom`); `tint` (integer tint index; default 1); `label_inside`;
#'   `dashed` (dashed outline); `centromere_start`, `centromere_end`
#'   (reference coordinates; `NA` means none).
#' @param aesthetics Optional `SeqAes`. Keys: `outline.col` (default
#'   `"grey25"`), `outline.lwd` (0.7), `label.size` (cex, 0.7), `label.col`
#'   (`"black"`), `centromere.fill` (`"#D33A2C"`), `fade` (alpha of a white
#'   overlay, 0), `arrow.head` (maximum head length as a fraction of block
#'   width, 0.35), `arrow.head_max` (maximum head length as a fraction of the
#'   track width, 0.012).
#' @param fill `"solid"` (default) or `"cytoband"`. Blocks with no matching
#'   cytobands fall back to solid with a message.
#' @param cytobands Cytoband table: a `GRanges` as returned by
#'   [load_cytobands()] or a data frame with `chrom`, `start` (0-based),
#'   `end`, `stain` (UCSC `chromStart`, `chromEnd`, `gieStain` also work).
#'   `NULL` loads the bundled hg38 bands. Chromosome names are matched with or
#'   without a `chr` prefix.
#' @param palette `NULL` (default hues), `"soft"` (a softer set), or a
#'   character vector of colours. A named vector is looked up by key; an
#'   unnamed one is recycled over the distinct keys in data order.
#' @param tints Fractions mixed toward white, indexed by `tint`.
#' @param block_colors Optional named colours overriding the fill by label,
#'   then by key (not tinted).
#' @param shape `"rect"` (default) or `"arrow"` (solid strips pointing toward
#'   increasing reference coordinate, reversed when flipped).
#' @param height Block height as a fraction of the panel height. Default 0.5.
#' @param label_side `"below"` (default) or `"above"` for outside labels.
#' @param flip_style How flipped blocks are marked: `"rotate"` (label rotated
#'   180 degrees), `"overline"` (line over the label) or `"minus"` (a leading
#'   minus sign).
#' @param annotations Optional data frame of point or interval annotations
#'   with columns `chrom`, `start`, optional `end`, `label`, and optional
#'   `side` (`"above"`/`"below"`), `color`, `italic`, `dot`, `strand`
#'   (`"+"`/`"-"` draws a direction arrow). Each row is drawn on every block
#'   whose reference interval contains its midpoint, mapped through flips.
#' @param legend Optional `LegendKey` list or `SeqLegendSpec`. No legend is
#'   generated automatically.
#' @param show_legend Logical. Default `TRUE`.
#' @param ... Additional arguments reserved for future use.
#' @return A `SeqBlockR6` instance.
#' @examples
#' gr <- GenomicRanges::GRanges("der1", IRanges::IRanges(c(1, 3e7), c(3e7, 6e7)),
#'                              label = c("A", "B"), flipped = c(FALSE, TRUE))
#' seq_block(gr, map(label = label, flipped = flipped, chrom = "chr1"))
#' @export
seq_block <- function(data = NULL, mapping = NULL, aesthetics = aes(),
                      fill = c("solid", "cytoband"), cytobands = NULL, palette = NULL,
                      tints = c(0.15, 0.5, 0.3, 0.65, 0.4), block_colors = NULL,
                      shape = c("rect", "arrow"), height = 0.5,
                      label_side = c("below", "above"),
                      flip_style = c("rotate", "overline", "minus"),
                      annotations = NULL, legend = NULL, show_legend = TRUE, ...) {
  SeqBlockR6$new(data = data, mapping = mapping, aesthetics = aesthetics,
                 fill = fill, cytobands = cytobands, palette = palette, tints = tints,
                 block_colors = block_colors, shape = shape, height = height,
                 label_side = label_side, flip_style = flip_style,
                 annotations = annotations, legend = legend,
                 show_legend = show_legend, ...)
}
