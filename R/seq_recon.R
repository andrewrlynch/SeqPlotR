# ── seq_recon — strand-classified SV arches ──────────────────────────────────

#' SeqRecon R6 class
#'
#' Internal R6 generator backing [seq_recon()]. Inherits from [`SeqArchR6`].
#' Classifies each link by the strand pair on its two anchors:
#' \itemize{
#'   \item `+/+` head-to-head inversion
#'   \item `-/-` tail-to-tail inversion
#'   \item `-/+` tandem duplication
#'   \item `+/-` deletion
#'   \item different chromosomes — translocation
#' }
#' Each class draws on its own vertical tier with a class-specific color and
#' a guide line + label per tier.
#'
#' @keywords internal
SeqReconR6 <- R6::R6Class("SeqRecon",
  inherit = SeqArchR6,
  public = list(
    #' @field col_h2h Color for `+/+` head-to-head inversions.
    col_h2h     = NULL,
    #' @field col_t2t Color for `-/-` tail-to-tail inversions.
    col_t2t     = NULL,
    #' @field col_dup Color for `-/+` tandem duplications.
    col_dup     = NULL,
    #' @field col_del Color for `+/-` deletions.
    col_del     = NULL,
    #' @field col_trans Color for cross-chromosome translocations.
    col_trans   = NULL,
    #' @field drawClasses Class tiers to draw.
    drawClasses = c("Inversion", "Dup/Del", "Translocation"),
    #' @field tierMultipliers Per-link tier value populated by `prep()`.
    tierMultipliers = NULL,
    #' @field last_arc_track Cached panels list of the arc track, used by
    #'   `draw()` to render tier guides and class labels.
    last_arc_track  = NULL,

    #' @description Construct a SeqReconR6.
    #' @param data BEDPE-like `GRanges` or `data.frame`.
    #' @param mapping `SeqMap`. Required: `x0`, `x1`, `strand0`, `strand1`.
    #'   For data.frame `data`, also required: `chrom0`, `chrom1`.
    #' @param t0,t1 Track identifiers.
    #' @param aesthetics Optional `SeqAes`. Override defaults via
    #'   `h2hColor`, `t2tColor`, `dupColor`, `delColor`, `transColor`.
    #' @param drawClasses Character vector of tiers to render.
    #' @param ... Reserved.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(),
                          drawClasses = c("Inversion", "Dup/Del",
                                          "Translocation"), ...) {
      self$col_h2h     <- aesthetics$h2hColor   %||% flexoki_palette(9)[3]
      self$col_t2t     <- aesthetics$t2tColor   %||% flexoki_palette(9)[4]
      self$col_dup     <- aesthetics$dupColor   %||% flexoki_palette(9)[1]
      self$col_del     <- aesthetics$delColor   %||% flexoki_palette(9)[2]
      self$col_trans   <- aesthetics$transColor %||% flexoki_palette(9)[9]
      self$drawClasses <- drawClasses
      super$initialize(data = data, mapping = mapping,
                       t0 = t0, t1 = t1, aesthetics = aesthetics)
    },

    #' @description Resolve both anchors, classify each link by strand and
    #'   chromosome, set per-link tier height + color, then delegate the
    #'   coordinate work to `SeqArchR6$prep()`.
    #' @param layout_all_tracks Named list of per-track panel-bounds lists.
    #' @param track_windows_list Named list of per-track `GRanges` windows.
    #' @param plot_track_index Fallback identifier for `t0`/`t1`.
    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      tid0 <- self$t0 %||% plot_track_index
      tid1 <- self$t1 %||% plot_track_index

      layout_t0 <- self$.resolve_track_ref(tid0, layout_all_tracks)
      self$.resolve_track_ref(tid1, layout_all_tracks)  # validate t1
      self$last_arc_track <- layout_t0

      self$resolve(
        track_data    = layout_t0[[1]]$track_data,
        track_mapping = layout_t0[[1]]$track_mapping
      )

      # seq_recon requires explicit strand0 / strand1 — we cannot meaningfully
      # classify SV type without them, so reject silent defaults.
      eff_mapping <- self$resolved$.mapping
      mapping_keys <- names(eff_mapping %||% list())
      missing_keys <- setdiff(c("strand0", "strand1"), mapping_keys)
      if (length(missing_keys) > 0L)
        stop("seq_recon requires the following map() field(s): ",
             paste(missing_keys, collapse = ", "), ".", call. = FALSE)

      if (is.null(self$anchor0_gr) || length(self$anchor0_gr) == 0L) {
        self$coordCanvas <- list()
        self$stubs       <- list()
        return(invisible())
      }
      n <- length(self$anchor0_gr)

      seq1 <- as.character(GenomicRanges::seqnames(self$anchor0_gr))
      seq2 <- as.character(GenomicRanges::seqnames(self$anchor1_gr))
      s1   <- as.character(BiocGenerics::strand(self$anchor0_gr))
      s2   <- as.character(BiocGenerics::strand(self$anchor1_gr))
      s1[!s1 %in% c("+", "-")] <- "+"
      s2[!s2 %in% c("+", "-")] <- "+"
      code <- paste0(s1, "/", s2)

      tier <- numeric(n)
      cols <- character(n)
      ori  <- character(n)

      for (i in seq_len(n)) {
        if (seq1[i] != seq2[i]) {
          tier[i] <- 1;   cols[i] <- self$col_trans; ori[i] <- "+"
        } else if (code[i] == "+/+") {
          tier[i] <- 0;   cols[i] <- self$col_h2h;   ori[i] <- "+"
        } else if (code[i] == "-/-") {
          tier[i] <- 0;   cols[i] <- self$col_t2t;   ori[i] <- "-"
        } else if (code[i] == "-/+") {
          tier[i] <- 0.5; cols[i] <- self$col_dup;   ori[i] <- "+"
        } else if (code[i] == "+/-") {
          tier[i] <- 0.5; cols[i] <- self$col_del;   ori[i] <- "-"
        } else {
          tier[i] <- 1;   cols[i] <- self$col_trans; ori[i] <- "+"
        }
      }

      self$tierMultipliers       <- tier
      self$aesthetics$height     <- tier
      self$aesthetics$arcColor   <- cols
      self$aesthetics$stemColor  <- cols
      self$aesthetics$orientation <- ori

      super$prep(layout_all_tracks  = layout_all_tracks,
                 track_windows_list = track_windows_list,
                 plot_track_index   = plot_track_index)
    },

    #' @description Draw class tier guides + labels, then the arches via
    #'   `SeqArchR6$draw()` (only for the `Inversion` tier so the labels are
    #'   layered correctly above the arches, matching THEfunc `SeqRecon`).
    draw = function() {
      if (is.null(self$last_arc_track)) return(invisible())

      panels <- self$last_arc_track
      x0s <- vapply(panels, function(pm) pm$inner$x0, numeric(1))
      x1s <- vapply(panels, function(pm) pm$inner$x1, numeric(1))
      y0s <- vapply(panels, function(pm) pm$inner$y0, numeric(1))
      y1s <- vapply(panels, function(pm) pm$inner$y1, numeric(1))
      tb_x0 <- min(x0s); tb_x1 <- max(x1s)
      tb_y0 <- min(y0s); tb_y1 <- max(y1s)
      ysc   <- panels[[1]]$yscale

      drawClasses <- self$drawClasses
      drawClasses <- stats::setNames(
        seq(0, 1, length.out = max(2L, length(drawClasses))),
        drawClasses
      )

      classes <- list(
        Inversion     = list(mult = drawClasses["Inversion"],
                             text = "HH/TT", color = self$col_h2h),
        `Dup/Del`     = list(mult = drawClasses["Dup/Del"],
                             text = c("DEL", "DUP"), color = self$col_del),
        Translocation = list(mult = drawClasses["Translocation"],
                             text = "TRA", color = self$col_trans)
      )

      for (cls in rev(names(drawClasses))) {
        info <- classes[[cls]]
        if (is.null(info)) next
        v_npc <- tb_y0 + (info$mult - ysc[1]) / diff(ysc) * (tb_y1 - tb_y0)

        grid::grid.lines(
          x  = grid::unit(c(tb_x0, tb_x1), "npc"),
          y  = grid::unit(rep(v_npc, 2),   "npc"),
          gp = grid::gpar(col = "grey40", lty = 3, lwd = 0.5)
        )

        if (cls == "Translocation") {
          grid::grid.text(
            label = info$text,
            x = grid::unit(tb_x0, "npc") - grid::unit(tb_x0 * 0.015, "npc"),
            y = grid::unit(v_npc, "npc"),
            just = "right",
            gp   = grid::gpar(col = "grey40", cex = 0.5)
          )
        } else if (cls == "Dup/Del") {
          grid::grid.text(
            "DEL",
            x = grid::unit(tb_x0, "npc") - grid::unit(tb_x0 * 0.015, "npc"),
            y = grid::unit(v_npc, "npc") - grid::unit(v_npc * 0.015, "npc"),
            just = "right",
            gp   = grid::gpar(col = self$col_del, cex = 0.5)
          )
          grid::grid.text(
            "DUP",
            x = grid::unit(tb_x0, "npc") - grid::unit(tb_x0 * 0.015, "npc"),
            y = grid::unit(v_npc, "npc") + grid::unit(v_npc * 0.015, "npc"),
            just = "right",
            gp   = grid::gpar(col = self$col_dup, cex = 0.5)
          )
        } else if (cls == "Inversion") {
          grid::grid.text(
            "TT",
            x = grid::unit(tb_x0, "npc") - grid::unit(tb_x0 * 0.015, "npc"),
            y = grid::unit(v_npc, "npc") - grid::unit(v_npc * 0.015, "npc"),
            just = "right",
            gp   = grid::gpar(col = self$col_t2t, cex = 0.5)
          )
          grid::grid.text(
            "HH",
            x = grid::unit(tb_x0, "npc") - grid::unit(tb_x0 * 0.015, "npc"),
            y = grid::unit(v_npc, "npc") + grid::unit(v_npc * 0.015, "npc"),
            just = "right",
            gp   = grid::gpar(col = self$col_h2h, cex = 0.5)
          )
          super$draw()
        }
      }
      invisible()
    }
  )
)

#' Strand-classified SV reconstruction arches
#'
#' Inherits everything from [seq_arch()] and additionally classifies each
#' link by strand orientation and chromosome pair, drawing each class on a
#' fixed tier with a class-specific color. Use to summarise structural
#' variant calls.
#'
#' Default colors come from [flexoki_palette()] and can be overridden via
#' `aes(h2hColor=, t2tColor=, dupColor=, delColor=, transColor=)`.
#'
#' Both `strand0` and `strand1` must be specified in `map()` — `seq_recon`
#' errors at `prep()` time if either is absent.
#'
#' @inheritParams seq_arc
#' @param drawClasses Character vector of tier labels to render
#'   (default: `c("Inversion", "Dup/Del", "Translocation")`).
#' @return A `SeqReconR6` instance.
#' @export
seq_recon <- function(data = NULL, mapping = NULL,
                      t0 = NULL, t1 = NULL,
                      aesthetics = aes(),
                      drawClasses = c("Inversion", "Dup/Del",
                                      "Translocation"), ...) {
  SeqReconR6$new(data = data, mapping = mapping,
                 t0 = t0, t1 = t1,
                 aesthetics = aesthetics,
                 drawClasses = drawClasses, ...)
}
