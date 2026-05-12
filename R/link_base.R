# ── SeqLinkR6 ────────────────────────────────────────────────────────────────
#
# Internal R6 base class for cross-track and within-track link elements
# (arcs, arches, strings, syntenies, zooms). Inherits from SeqElementR6 and
# adds:
#   * `t0`, `t1`              — track identifiers for the two anchors
#   * `anchor0_gr`, `anchor1_gr` — synthetic point GRanges built from the
#     resolved `map()` fields, populated by `resolve()`
#
# The link API is single-data: paired loci are encoded in one `data`
# argument (a BEDPE-like data.frame, or a GRanges with bedpe-like mcols).
# The `map()` vocabulary uses suffix-`0`/`1` field names: `x0`, `x1`,
# `chrom0`, `chrom1`, `strand0`, `strand1`, `y0`, `y1`, `height`.

#' SeqLink R6 base class
#'
#' Internal R6 generator for cross-track and within-track link elements.
#' Extends [`SeqElementR6`] with track identifiers and a pair of anchor
#' GRanges synthesised from the resolved `map()` fields.
#'
#' @keywords internal
SeqLinkR6 <- R6::R6Class("SeqLink",
  inherit = SeqElementR6,
  public = list(
    #' @field t0 `track_id` (or integer index) of the first anchor's track.
    t0         = NULL,
    #' @field t1 `track_id` (or integer index) of the second anchor's track.
    t1         = NULL,
    #' @field anchor0_gr Synthetic point GRanges for anchor 0, populated by
    #'   `resolve()` from the resolved `chrom0`, `x0`, and `strand0` fields.
    anchor0_gr = NULL,
    #' @field anchor1_gr Synthetic point GRanges for anchor 1.
    anchor1_gr = NULL,

    #' @description Construct a new SeqLinkR6.
    #' @param data Optional `GRanges` or `data.frame`. A single argument
    #'   carries both anchors; there is no `data2`.
    #' @param mapping Optional `SeqMap`. Must define `x0` and `x1` (and,
    #'   for data.frame `data`, `chrom0` and `chrom1`).
    #' @param t0,t1 Track identifiers (string or integer index). Locked to
    #'   the parent track when added inside a [seq_track()] via `%+%`.
    #' @param aesthetics Optional `SeqAes`.
    #' @param legend Optional `LegendKey` or list of `LegendKey` objects.
    #' @param show_legend Logical. Set to `FALSE` to suppress legend output.
    #' @param ... Unused — accepted so subclasses can pass extra arguments.
    initialize = function(data = NULL, mapping = NULL,
                          t0 = NULL, t1 = NULL,
                          aesthetics = aes(),
                          legend = NULL, show_legend = TRUE, ...) {
      super$initialize(data, mapping, aesthetics)
      self$t0          <- t0
      self$t1          <- t1
      self$legend      <- legend
      self$show_legend <- show_legend
    },

    #' @description Resolve the mapping against the effective `data`, then
    #'   synthesise `anchor0_gr` and `anchor1_gr` from the resolved fields.
    #'   Errors if `x0` or `x1` is missing, or if `chrom0`/`chrom1` cannot
    #'   be derived (only auto-fillable for GRanges `data`).
    #' @param track_data Optional `GRanges`/`data.frame` from the parent
    #'   track (used as a fallback when the link has no own `data`).
    #' @param track_mapping Optional `SeqMap` from the parent track.
    #' @return The link, invisibly.
    resolve = function(track_data = NULL, track_mapping = NULL) {
      super$resolve(track_data, track_mapping)

      eff_data <- self$resolved$.data
      self$anchor0_gr <- NULL
      self$anchor1_gr <- NULL
      if (is.null(eff_data)) return(invisible(self))

      n <- if (is.data.frame(eff_data)) nrow(eff_data) else length(eff_data)
      if (n == 0L) return(invisible(self))

      x0_v <- self$resolved$x0
      x1_v <- self$resolved$x1
      if (is.null(x0_v) || is.null(x1_v))
        stop("seq_link mapping must define both x0 and x1.", call. = FALSE)

      is_gr <- inherits(eff_data, "GRanges")

      chrom0_v <- self$resolved$chrom0 %||%
        (if (is_gr) as.character(GenomicRanges::seqnames(eff_data))
         else stop("seq_link with data.frame data requires chrom0 in map().",
                   call. = FALSE))
      chrom1_v <- self$resolved$chrom1 %||%
        (if (is_gr) as.character(GenomicRanges::seqnames(eff_data))
         else stop("seq_link with data.frame data requires chrom1 in map().",
                   call. = FALSE))

      strand0_v <- self$resolved$strand0 %||%
        (if (is_gr) as.character(BiocGenerics::strand(eff_data))
         else rep("*", n))
      strand1_v <- self$resolved$strand1 %||%
        (if (is_gr) as.character(BiocGenerics::strand(eff_data))
         else rep("*", n))

      .vec <- function(v, n) if (length(v) == 1L) rep(v, n) else v
      x0_v      <- .vec(x0_v,      n); x1_v      <- .vec(x1_v,      n)
      chrom0_v  <- .vec(chrom0_v,  n); chrom1_v  <- .vec(chrom1_v,  n)
      strand0_v <- .vec(strand0_v, n); strand1_v <- .vec(strand1_v, n)

      self$anchor0_gr <- GenomicRanges::GRanges(
        seqnames = chrom0_v,
        ranges   = IRanges::IRanges(start = x0_v, width = 1L),
        strand   = strand0_v
      )
      self$anchor1_gr <- GenomicRanges::GRanges(
        seqnames = chrom1_v,
        ranges   = IRanges::IRanges(start = x1_v, width = 1L),
        strand   = strand1_v
      )
      invisible(self)
    },

    #' @description Override in subclasses. Default implementation errors.
    #' @param layout_all_tracks Named list of panel bounds keyed by `track_id`.
    #' @param track_windows_list Named list of `GRanges` windows keyed by
    #'   `track_id`.
    #' @param plot_track_index Optional integer index of the parent track,
    #'   used only to set defaults for within-track links.
    prep = function(layout_all_tracks, track_windows_list,
                    plot_track_index = NULL) {
      stop("prep() must be implemented by ", class(self)[1], call. = FALSE)
    },

    #' @description Look up a track reference (`track_id` string or integer
    #'   index) in `layout_all_tracks` and return its panel-bounds list.
    #'   Errors with a clear message when the reference cannot be resolved.
    #' @param ref A `track_id` (character) or integer index.
    #' @param layout_all_tracks Named list of per-track panel bounds (the
    #'   `panelBounds` produced by [`SeqPlotR6`]'s `layoutGrid()`).
    #' @return The selected track's panel-bounds list.
    .resolve_track_ref = function(ref, layout_all_tracks) {
      if (is.null(ref))
        stop("seq_link reference is NULL - t0/t1 must be set.", call. = FALSE)

      if (is.character(ref)) {
        if (ref %in% names(layout_all_tracks))
          return(layout_all_tracks[[ref]])
        stop("seq_link references unknown track_id '", ref, "'.", call. = FALSE)
      }

      if (is.numeric(ref)) {
        idx <- as.integer(ref)
        if (idx >= 1L && idx <= length(layout_all_tracks))
          return(layout_all_tracks[[idx]])
        stop("seq_link references track index ", idx,
             " which is out of bounds (1..", length(layout_all_tracks), ").",
             call. = FALSE)
      }

      stop("seq_link reference must be a track_id string or integer index.",
           call. = FALSE)
    }
  )
)
