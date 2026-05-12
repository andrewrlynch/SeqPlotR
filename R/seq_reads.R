# ── SeqReadsR6 ────────────────────────────────────────────────────────────────
#
# IGV-style read alignment track. Loads reads from an indexed BAM via
# Rsamtools / GenomicAlignments at construction time, packs them into rows,
# and renders chevron polygons (right-facing for + strand, left for − strand)
# during draw(). Mate pairs are optionally connected by a horizontal line.
#
# Inherits from SeqElementR6 so it slots into seq_track / seq_plot like any
# other element.

#' SeqReads R6 class
#'
#' Internal R6 generator backing [seq_reads()]. Inherits from
#' [`SeqElementR6`].
#'
#' @keywords internal
SeqReadsR6 <- R6::R6Class("SeqReads",
  inherit = SeqElementR6,
  public = list(
    #' @field bam Path to the BAM file.
    bam        = NULL,
    #' @field region Sorted GRanges of windows to load.
    region     = NULL,
    #' @field gr_reads GRanges of loaded reads with `qname` and `y` mcols.
    gr_reads   = NULL,
    #' @field link_spans GRanges of mate-pair spans with `y` mcols.
    link_spans = NULL,
    #' @field nrows Maximum row index used for read packing.
    nrows      = 0L,
    #' @field min_mapq Minimum mapping quality threshold.
    min_mapq   = 0L,
    #' @field max_reads Maximum number of reads to load.
    max_reads  = 20000L,
    #' @field max_width Guardrail on per-window width in bp.
    max_width  = 100000L,
    #' @field sort_by Row-packing sort order (`"insert_length"` or `"start"`).
    sort_by    = "insert_length",
    #' @field link_mates Whether to draw mate-pair link lines.
    link_mates = TRUE,
    #' @field show_strand Whether to render strand chevrons.
    show_strand = TRUE,

    #' @description Construct a SeqReadsR6 from a BAM file and region.
    #' @param bam Character. Path to an indexed BAM file.
    #' @param region `GRanges` of windows to load.
    #' @param min_mapq Integer. Minimum mapping quality. Default `0`.
    #' @param max_reads Integer. Maximum total reads to load. Default `20000`.
    #' @param max_width Integer. Maximum allowed window width in bp. Default
    #'   `100000`.
    #' @param sort_by Character. `"insert_length"` (default) or `"start"`.
    #' @param link_mates Logical. Default `TRUE`.
    #' @param show_strand Logical. Default `TRUE`.
    #' @param aesthetics Named list of optional aesthetics. See `seq_reads()`.
    #' @param ... Reserved.
    initialize = function(bam,
                          region,
                          min_mapq    = 0L,
                          max_reads   = 20000L,
                          max_width   = 100000L,
                          sort_by     = c("insert_length", "start"),
                          link_mates  = TRUE,
                          show_strand = TRUE,
                          aesthetics  = list(),
                          ...) {

      if (!inherits(region, "GRanges") || length(region) < 1L)
        stop("`region` must be a GRanges with at least one range.",
             call. = FALSE)
      wide <- which(BiocGenerics::width(region) > max_width)
      if (length(wide)) {
        stop(sprintf(
          "%d window(s) exceed max_width (%d bp). Narrow your windows or ",
          length(wide), max_width),
          "increase max_width with care (large windows may load millions of reads).",
          call. = FALSE)
      }
      if (!is.character(bam) || length(bam) != 1L || !nzchar(bam))
        stop("`bam` must be a single non-empty path.", call. = FALSE)
      if (!file.exists(bam))
        stop("BAM file not found: ", bam, call. = FALSE)
      for (pkg in c("Rsamtools", "GenomicAlignments")) {
        if (!requireNamespace(pkg, quietly = TRUE))
          stop("seq_reads() requires '", pkg, "'.", call. = FALSE)
      }

      sort_by <- match.arg(sort_by)
      default_aes <- list(
        row_gap          = 0.12,
        tip_mm           = 2.0,
        tip_min_body_mm  = 1.0,
        col              = NA,
        lwd              = 0.3,
        fill_plus        = "steelblue3",
        fill_minus       = "tomato3",
        fill_unstranded  = "grey70",
        link_col         = "grey60",
        link_lwd         = 0.6,
        link_lty         = 1
      )

      self$bam         <- bam
      self$region      <- GenomicRanges::sort(region)
      self$min_mapq    <- as.integer(min_mapq)
      self$max_reads   <- as.integer(max_reads)
      self$max_width   <- as.integer(max_width)
      self$sort_by     <- sort_by
      self$link_mates  <- isTRUE(link_mates)
      self$show_strand <- isTRUE(show_strand)
      effective_aes    <- modifyList(default_aes, aesthetics)

      self$gr_reads <- private$read_bam_window(
        bam       = self$bam,
        region    = self$region,
        min_mapq  = self$min_mapq,
        max_reads = self$max_reads,
        sort_by   = self$sort_by
      )

      if (length(self$gr_reads) == 0L) {
        self$nrows      <- 0L
        self$link_spans <- GenomicRanges::GRanges()
        super$initialize(data = self$region, aesthetics = effective_aes)
        return(invisible())
      }

      y <- S4Vectors::mcols(self$gr_reads)$y
      self$nrows <- if (length(y)) max(as.integer(y)) else 0L

      if (isTRUE(self$link_mates))
        self$link_spans <- private$make_link_spans(self$gr_reads)
      else
        self$link_spans <- GenomicRanges::GRanges()

      super$initialize(data = self$gr_reads, aesthetics = effective_aes)
    },

    #' @description Return the y-axis limits for the row-packed reads.
    getYLimits = function() {
      if (is.null(self$nrows) || self$nrows <= 0L) return(c(0, 1))
      c(0.5, self$nrows + 0.5)
    },

    #' @description Infer a continuous y scale covering the row-packed
    #'   reads. Used by `seq_plot$layoutGrid()` to pick a primary y scale
    #'   when the parent track does not define one.
    .infer_scale_y = function() {
      lims <- self$getYLimits()
      seq_scale_continuous(limits = lims)
    },

    #' @description Convert reads and link spans to canvas NPC coordinates,
    #'   one entry per window.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows Track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$coordCanvas <- vector("list", length(track_windows))
      if (is.null(self$gr_reads) || length(self$gr_reads) == 0L)
        return(invisible())

      # --- Reads (chevron polygons) ---
      ov <- GenomicRanges::findOverlaps(self$gr_reads, track_windows,
                                        ignore.strand = TRUE)
      if (length(ov) > 0L) {
        qh <- S4Vectors::queryHits(ov)
        sh <- S4Vectors::subjectHits(ov)

        x0_all <- BiocGenerics::start(self$gr_reads)[qh]
        x1_all <- BiocGenerics::end(self$gr_reads)[qh]
        y_all  <- S4Vectors::mcols(self$gr_reads)$y[qh]
        st_all <- as.character(BiocGenerics::strand(self$gr_reads)[qh])
        qn_all <- S4Vectors::mcols(self$gr_reads)$qname[qh]

        for (w in unique(sh)) {
          mask <- sh == w
          p    <- layout_track[[w]]

          x0 <- x0_all[mask]; x1 <- x1_all[mask]
          y  <- y_all[mask];  st <- st_all[mask]
          qn <- qn_all[mask]

          clip <- .clip_to_windows(x0, x1, p$xscale)
          if (length(clip$x0) == 0L) { self$coordCanvas[[w]] <- NULL; next }
          x0 <- clip$x0; x1 <- clip$x1
          y  <- y[clip$mask]; st <- st[clip$mask]; qn <- qn[clip$mask]

          u0 <- (x0 - p$xscale[1]) / diff(p$xscale)
          u1 <- (x1 - p$xscale[1]) / diff(p$xscale)
          xleft  <- p$inner$x0 + u0 * (p$inner$x1 - p$inner$x0)
          xright <- p$inner$x0 + u1 * (p$inner$x1 - p$inner$x0)

          v_center <- (y - p$yscale[1]) / diff(p$yscale)
          v_center <- pmax(pmin(v_center, 1), 0)
          row_h    <- (p$inner$y1 - p$inner$y0) /
                      max(1e-9, diff(p$yscale))
          gap   <- self$aesthetics$row_gap * row_h
          half  <- (row_h - gap) / 2
          ycenter <- p$inner$y0 + v_center * (p$inner$y1 - p$inner$y0)

          self$coordCanvas[[w]] <- data.frame(
            x0     = xleft,
            x1     = xright,
            y0     = ycenter - half,
            y1     = ycenter + half,
            yc     = ycenter,
            strand = st,
            qname  = qn,
            stringsAsFactors = FALSE
          )
        }
      }

      # --- Link spans ---
      if (isTRUE(self$link_mates) &&
          !is.null(self$link_spans) &&
          length(self$link_spans) > 0L) {
        ov2 <- GenomicRanges::findOverlaps(self$link_spans, track_windows,
                                           ignore.strand = TRUE)
        if (length(ov2) > 0L) {
          qh2 <- S4Vectors::queryHits(ov2)
          sh2 <- S4Vectors::subjectHits(ov2)

          sp_start <- BiocGenerics::start(self$link_spans)[qh2]
          sp_end   <- BiocGenerics::end(self$link_spans)[qh2]
          sp_row   <- S4Vectors::mcols(self$link_spans)$y[qh2]

          for (w in unique(sh2)) {
            mask2 <- sh2 == w
            p     <- layout_track[[w]]

            s0 <- sp_start[mask2]; s1 <- sp_end[mask2]
            ry <- sp_row[mask2]

            clip2 <- .clip_to_windows(s0, s1, p$xscale)
            if (length(clip2$x0) == 0L) next
            s0 <- clip2$x0; s1 <- clip2$x1; ry <- ry[clip2$mask]

            u0 <- (s0 - p$xscale[1]) / diff(p$xscale)
            u1 <- (s1 - p$xscale[1]) / diff(p$xscale)
            x0c <- p$inner$x0 + u0 * (p$inner$x1 - p$inner$x0)
            x1c <- p$inner$x0 + u1 * (p$inner$x1 - p$inner$x0)

            v   <- (ry - p$yscale[1]) / diff(p$yscale)
            v   <- pmax(pmin(v, 1), 0)
            yc  <- p$inner$y0 + v * (p$inner$y1 - p$inner$y0)

            link_df <- data.frame(
              x0 = x0c, x1 = x1c, y0 = yc, y1 = yc,
              stringsAsFactors = FALSE
            )

            if (is.null(self$coordCanvas[[w]])) {
              self$coordCanvas[[w]] <- list(reads = NULL, links = link_df)
            } else if (is.data.frame(self$coordCanvas[[w]])) {
              self$coordCanvas[[w]] <- list(reads = self$coordCanvas[[w]],
                                            links = link_df)
            } else {
              self$coordCanvas[[w]]$links <- link_df
            }
          }
        }
      }

      # Normalise: bare data.frame → list(reads, links=NULL) for uniformity.
      for (w in seq_along(self$coordCanvas)) {
        if (is.data.frame(self$coordCanvas[[w]])) {
          self$coordCanvas[[w]] <- list(reads = self$coordCanvas[[w]],
                                        links = NULL)
        }
      }
      invisible()
    },

    #' @description Render reads (as chevron polygons) and mate-pair links.
    draw = function() {
      if (is.null(self$coordCanvas)) return(invisible())

      tip_npc <- tryCatch(
        as.numeric(grid::convertWidth(
          grid::unit(self$aesthetics$tip_mm, "mm"), "npc",
          valueOnly = TRUE)),
        error = function(e) 0.01
      )
      min_body_npc <- tryCatch(
        as.numeric(grid::convertWidth(
          grid::unit(self$aesthetics$tip_min_body_mm, "mm"), "npc",
          valueOnly = TRUE)),
        error = function(e) 0.005
      )

      for (win_coords in self$coordCanvas) {
        if (is.null(win_coords)) next

        reads <- win_coords$reads
        links <- win_coords$links

        if (!is.null(links) && nrow(links) > 0L) {
          grid::grid.segments(
            x0 = grid::unit(links$x0, "npc"),
            x1 = grid::unit(links$x1, "npc"),
            y0 = grid::unit(links$y0, "npc"),
            y1 = grid::unit(links$y1, "npc"),
            gp = grid::gpar(col = self$aesthetics$link_col,
                            lwd = self$aesthetics$link_lwd,
                            lty = self$aesthetics$link_lty)
          )
        }

        if (is.null(reads) || nrow(reads) == 0L) next

        fill <- rep(self$aesthetics$fill_unstranded, nrow(reads))
        fill[reads$strand == "+"] <- self$aesthetics$fill_plus
        fill[reads$strand == "-"] <- self$aesthetics$fill_minus

        for (i in seq_len(nrow(reads))) {
          x0 <- reads$x0[i]; x1 <- reads$x1[i]
          y0 <- reads$y0[i]; y1 <- reads$y1[i]
          st <- reads$strand[i]
          w_npc <- abs(x1 - x0)
          tip   <- min(tip_npc, max(0, w_npc - min_body_npc))
          ym    <- (y0 + y1) / 2

          if (!isTRUE(self$show_strand) ||
              !(st %in% c("+", "-")) ||
              w_npc <= 0) {
            xs <- c(x0, x1, x1, x0)
            ys <- c(y0, y0, y1, y1)
          } else if (st == "+") {
            xb <- max(x0, x1 - tip)
            xs <- c(x0, xb, x1, xb, x0)
            ys <- c(y0, y0, ym, y1, y1)
          } else {
            xb <- min(x1, x0 + tip)
            xs <- c(x1, xb, x0, xb, x1)
            ys <- c(y0, y0, ym, y1, y1)
          }

          grid::grid.polygon(
            x  = grid::unit(xs, "npc"),
            y  = grid::unit(ys, "npc"),
            gp = grid::gpar(
              fill = fill[i],
              col  = self$aesthetics$col,
              lwd  = self$aesthetics$lwd
            )
          )
        }
      }
      invisible()
    }
  ),

  private = list(

    read_bam_window = function(bam, region, min_mapq, max_reads, sort_by) {
      bf <- Rsamtools::BamFile(bam)
      param <- Rsamtools::ScanBamParam(
        which = region,
        what  = c("qname", "flag", "pos", "qwidth",
                  "mapq", "rname", "cigar"),
        flag  = Rsamtools::scanBamFlag(
          isUnmappedQuery          = FALSE,
          isSecondaryAlignment     = FALSE,
          isSupplementaryAlignment = FALSE
        )
      )
      gal <- GenomicAlignments::readGAlignments(bf, param = param,
                                                use.names = TRUE)
      if (length(gal) == 0L) return(GenomicRanges::GRanges())

      mq   <- S4Vectors::mcols(gal)$mapq
      keep <- is.na(mq) | mq >= min_mapq
      gal  <- gal[keep]
      if (length(gal) == 0L) return(GenomicRanges::GRanges())

      # qname-aware downsampling — keep mate pairs intact.
      if (length(gal) > max_reads) {
        qn_all <- S4Vectors::mcols(gal)$qname
        tab    <- sort(table(qn_all), decreasing = TRUE)
        keep_q <- character(0); n <- 0L
        for (q in names(tab)) {
          keep_q <- c(keep_q, q)
          n <- n + as.integer(tab[[q]])
          if (n >= max_reads) break
        }
        gal <- gal[qn_all %in% keep_q]
      }

      gr <- GenomicAlignments::granges(gal)
      S4Vectors::mcols(gr)$qname <- S4Vectors::mcols(gal)$qname

      qn <- S4Vectors::mcols(gr)$qname
      df <- data.frame(
        qname = qn,
        seq   = as.character(GenomicRanges::seqnames(gr)),
        start = BiocGenerics::start(gr),
        end   = BiocGenerics::end(gr),
        stringsAsFactors = FALSE
      )
      sp <- split(df, df$qname)
      span_start <- vapply(sp, function(d) min(d$start), numeric(1))
      span_end   <- vapply(sp, function(d) max(d$end),   numeric(1))
      insert_len <- span_end - span_start
      names(insert_len) <- names(sp)

      if (sort_by == "insert_length") {
        sorted_q <- names(sort(insert_len, decreasing = TRUE))
      } else {
        sorted_q <- names(sort(span_start))
      }

      span_spans <- GenomicRanges::GRanges(
        seqnames = vapply(sp, function(d) d$seq[1], character(1)),
        ranges   = IRanges::IRanges(start = span_start, end = span_end)
      )
      names(span_spans) <- names(sp)
      ord_spans <- span_spans[sorted_q]
      span_bins <- IRanges::disjointBins(IRanges::ranges(ord_spans))
      span_row  <- as.integer(span_bins)
      names(span_row) <- sorted_q

      S4Vectors::mcols(gr)$y <- span_row[qn]
      gr
    },

    make_link_spans = function(gr) {
      qn  <- S4Vectors::mcols(gr)$qname
      tab <- table(qn)
      keep_q <- names(tab)[tab >= 2L]
      if (length(keep_q) == 0L) return(GenomicRanges::GRanges())

      sub <- gr[qn %in% keep_q]
      qn2 <- S4Vectors::mcols(sub)$qname
      df  <- data.frame(
        qname = qn2,
        seq   = as.character(GenomicRanges::seqnames(sub)),
        start = BiocGenerics::start(sub),
        end   = BiocGenerics::end(sub),
        y     = S4Vectors::mcols(sub)$y,
        stringsAsFactors = FALSE
      )
      sp <- split(df, df$qname)
      out <- GenomicRanges::GRanges(
        seqnames = vapply(sp, function(d) d$seq[1], character(1)),
        ranges   = IRanges::IRanges(
          start = vapply(sp, function(d) min(d$start), numeric(1)),
          end   = vapply(sp, function(d) max(d$end),   numeric(1))
        )
      )
      S4Vectors::mcols(out)$y     <- vapply(sp, function(d) d$y[1], numeric(1))
      S4Vectors::mcols(out)$qname <- names(sp)
      out
    }
  )
)

#' IGV-style read alignment track element
#'
#' Renders read alignments from an indexed BAM file as chevron polygons
#' (pointing right for + strand reads, left for −). Mate pairs are
#' optionally connected by a thin horizontal line. Reads are loaded once
#' at construction time and cached on the element; only the requested
#' windows are pulled.
#'
#' @param bam Character. Path to an indexed BAM file.
#' @param region `GRanges` of windows to load. Every window must be
#'   ≤ `max_width` bp.
#' @param min_mapq Integer. Minimum mapping quality. Default `0`.
#' @param max_reads Integer. Maximum total reads to load (qname-aware
#'   downsampling preserves mate pairs). Default `20000`.
#' @param max_width Integer. Maximum allowed window width in bp. Default
#'   `100000`. Prevents accidental loading of entire chromosomes.
#' @param sort_by Character. Row-packing sort order: `"insert_length"`
#'   (default — longest inserts at the top) or `"start"` (leftmost reads
#'   at the top).
#' @param link_mates Logical. Draw a horizontal line connecting each
#'   mated pair's outer extent. Default `TRUE`.
#' @param show_strand Logical. Draw chevron arrowheads indicating strand.
#'   Default `TRUE`.
#' @param aesthetics Named list. Supported keys: `fill_plus`,
#'   `fill_minus`, `fill_unstranded`, `col`, `lwd`, `row_gap`, `tip_mm`,
#'   `tip_min_body_mm`, `link_col`, `link_lwd`, `link_lty`.
#' @return A `SeqReadsR6` instance.
#' @export
seq_reads <- function(bam,
                      region,
                      min_mapq    = 0L,
                      max_reads   = 20000L,
                      max_width   = 100000L,
                      sort_by     = c("insert_length", "start"),
                      link_mates  = TRUE,
                      show_strand = TRUE,
                      aesthetics  = list()) {
  SeqReadsR6$new(
    bam         = bam,
    region      = region,
    min_mapq    = min_mapq,
    max_reads   = max_reads,
    max_width   = max_width,
    sort_by     = sort_by,
    link_mates  = link_mates,
    show_strand = show_strand,
    aesthetics  = aesthetics
  )
}
