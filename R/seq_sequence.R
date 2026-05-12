# ── SeqSequenceR6 ─────────────────────────────────────────────────────────────
#
# IGV-style nucleotide sequence element: renders per-base coloured rectangles
# (and optionally letters) for genomic windows up to 200 bp wide.

# Default UCSC nucleotide colours
.SEQ_COLORS_UCSC <- c(
  A = "#00AA00",   # green
  T = "#FF0000",   # red
  C = "#0000FF",   # blue
  G = "#FFB300",   # orange/gold
  N = "#AAAAAA"    # grey for unknown
)

#' SeqSequence R6 class
#'
#' Internal R6 generator backing [seq_sequence()]. Inherits from
#' [`SeqElementR6`]. Renders per-base coloured rectangles (and optionally
#' letters) for genomic windows up to 200 bp wide. Wider windows emit a
#' message and render nothing.
#'
#' @keywords internal
SeqSequenceR6 <- R6::R6Class("SeqSequence",
  inherit = SeqElementR6,
  public = list(

    #' @field genome Character. BSgenome package name (e.g.
    #'   `"BSgenome.Hsapiens.UCSC.hg38"`). Required when `sequence` is `NULL`.
    genome = NULL,

    #' @field sequence Character string of nucleotides. When provided,
    #'   `genome` is ignored. The string is assumed to span the first window.
    sequence = NULL,

    #' @field show_letters Logical. When `TRUE` and window width <= 80 bp,
    #'   draw the nucleotide letter centred on each rectangle. Default `FALSE`.
    show_letters = FALSE,

    #' @field rect_height Numeric in (0, 1] or `NULL`. Fraction of track
    #'   height allocated to each nucleotide rectangle. `NULL` (default)
    #'   uses the full track height.
    rect_height = NULL,

    #' @field colors Named character vector mapping nucleotide characters to
    #'   hex color strings. Defaults to UCSC standard.
    colors = NULL,

    #' @description Construct a SeqSequenceR6.
    #' @param data Ignored; sequence is fetched from `genome` or `sequence`.
    #' @param mapping Ignored.
    #' @param aesthetics Optional `SeqAes`. Supports `color` (letter color,
    #'   defaults to matching the rectangle fill), `background` (letter
    #'   background rectangle color, default `NA`).
    #' @param genome Character. BSgenome package name. Required when `sequence`
    #'   is `NULL`.
    #' @param sequence Character string of nucleotides spanning the first
    #'   window. When provided, `genome` is ignored.
    #' @param show_letters Logical. Show nucleotide letters when window width
    #'   is <= 80 bp. Default `FALSE`.
    #' @param rect_height Numeric in (0, 1] or `NULL`. Fraction of track
    #'   height used for the rectangles. Default `NULL` (full track height).
    #' @param colors Named character vector mapping nucleotide codes to hex
    #'   colors. Defaults to UCSC standard.
    #' @param ... Reserved.
    initialize = function(data        = NULL,
                          mapping     = NULL,
                          aesthetics  = aes(),
                          genome      = NULL,
                          sequence    = NULL,
                          show_letters = FALSE,
                          rect_height  = NULL,
                          colors       = NULL,
                          ...) {
      super$initialize(data, mapping, aesthetics)
      if (is.null(genome) && is.null(sequence))
        stop("`seq_sequence()` requires either `genome` or `sequence`.",
             call. = FALSE)
      self$genome       <- genome
      self$sequence     <- sequence
      self$show_letters <- isTRUE(show_letters)
      self$rect_height  <- rect_height
      self$colors       <- colors %||% .SEQ_COLORS_UCSC
    },

    #' @description Fetch or validate the nucleotide sequence, check window
    #'   length, and compute per-base canvas coordinates.
    #' @param layout_track Per-window panel metadata list.
    #' @param track_windows The current track's `windows` GRanges.
    prep = function(layout_track, track_windows) {
      self$coordCanvas <- data.frame()

      win       <- track_windows[1]
      win_start <- BiocGenerics::start(win)
      win_end   <- BiocGenerics::end(win)
      win_width <- win_end - win_start + 1L

      if (win_width > 200L) {
        message("seq_sequence: window is ", win_width,
                " bp (> 200 bp) — nucleotide display suppressed.")
        return(invisible())
      }

      # --- Obtain sequence ---
      if (!is.null(self$sequence)) {
        seq_str <- self$sequence
        if (nchar(seq_str) < win_width) {
          warning("Provided sequence is shorter than the window (",
                  nchar(seq_str), " < ", win_width, " bp); padding with N.",
                  call. = FALSE)
          seq_str <- paste0(seq_str,
                            strrep("N", win_width - nchar(seq_str)))
        }
        seq_str <- substr(seq_str, 1L, win_width)
      } else {
        # BSgenome pull
        if (!requireNamespace("BSgenome", quietly = TRUE))
          stop("Package 'BSgenome' is required for automatic sequence fetching.",
               call. = FALSE)
        bsg_pkg <- self$genome
        if (!requireNamespace(bsg_pkg, quietly = TRUE))
          stop("BSgenome package '", bsg_pkg, "' is not installed.",
               call. = FALSE)
        bsg  <- get(bsg_pkg, envir = asNamespace(bsg_pkg))
        chr  <- as.character(GenomicRanges::seqnames(win))
        seq_str <- as.character(
          BSgenome::getSeq(bsg,
                           GenomicRanges::GRanges(chr,
                             IRanges::IRanges(win_start, win_end)))
        )
      }

      seq_str  <- toupper(seq_str)
      bases    <- strsplit(seq_str, "")[[1]]
      n_bases  <- length(bases)

      pm  <- layout_track[[1]]
      x0p <- pm$inner$x0
      x1p <- pm$inner$x1
      y0p <- pm$inner$y0
      y1p <- pm$inner$y1
      pan_w <- x1p - x0p

      rect_h_frac <- self$rect_height %||% 1
      rect_h_frac <- max(0.01, min(1, rect_h_frac))
      mid_y   <- (y0p + y1p) / 2
      half_h  <- (y1p - y0p) / 2 * rect_h_frac
      rect_y0 <- mid_y - half_h
      rect_y1 <- mid_y + half_h

      base_w <- pan_w / n_bases

      colors_out <- character(n_bases)
      for (i in seq_len(n_bases)) {
        b <- bases[i]
        colors_out[i] <- self$colors[b] %||% self$colors["N"] %||% "#AAAAAA"
      }

      self$coordCanvas <- data.frame(
        base      = bases,
        x0        = x0p + (seq_len(n_bases) - 1L) * base_w,
        x1        = x0p + seq_len(n_bases) * base_w,
        y0        = rect_y0,
        y1        = rect_y1,
        mid_y     = mid_y,
        color     = colors_out,
        win_width = win_width,
        stringsAsFactors = FALSE
      )
      invisible()
    },

    #' @description Draw per-base coloured rectangles and optional letters.
    draw = function() {
      df <- self$coordCanvas
      if (is.null(df) || nrow(df) == 0) return(invisible())

      n <- nrow(df)
      # Rectangles — 0.99 width leaves a thin gap between positions
      grid::grid.rect(
        x      = grid::unit((df$x0 + df$x1) / 2, "npc"),
        y      = grid::unit((df$y0 + df$y1) / 2, "npc"),
        width  = grid::unit((df$x1 - df$x0) * 0.99, "npc"),
        height = grid::unit(df$y1 - df$y0,           "npc"),
        gp     = grid::gpar(fill = df$color, col = NA)
      )

      # Letters (only when show_letters and window <= 80 bp)
      if (isTRUE(self$show_letters) && df$win_width[1] <= 80L) {
        aes_txt  <- self$aesthetics
        text_col <- aes_txt$color %||% NULL   # NULL = per-base contrast

        # Background rect (opt-in via aes(background = "white") etc.)
        bg_col <- aes_txt$background %||% NA
        if (!is.na(bg_col)) {
          grid::grid.rect(
            x      = grid::unit((df$x0 + df$x1) / 2, "npc"),
            y      = grid::unit(df$mid_y, "npc"),
            width  = grid::unit(df$x1 - df$x0, "npc"),
            height = grid::unit(df$y1 - df$y0, "npc"),
            gp     = grid::gpar(fill = bg_col, col = NA)
          )
        }

        # Per-base letter color: use aes color if set, else match rect fill
        letter_cols <- if (!is.null(text_col))
          rep(text_col, n)
        else
          df$color

        grid::grid.text(
          label = df$base,
          x     = grid::unit((df$x0 + df$x1) / 2, "npc"),
          y     = grid::unit(df$mid_y, "npc"),
          gp    = grid::gpar(cex = 0.7, col = letter_cols, fontface = "bold")
        )
      }

      invisible()
    }
  )
)

#' IGV-style nucleotide sequence track element
#'
#' Renders per-base coloured rectangles (and optionally letters) for genomic
#' windows up to 200 bp wide. Wider windows emit a message and show nothing.
#'
#' @param data Ignored; sequence is fetched from `genome` or `sequence`.
#' @param mapping Ignored.
#' @param aesthetics Optional `SeqAes`. Supports `color` (letter color,
#'   defaults to matching the rectangle fill), `background` (letter
#'   background rectangle color, default `NA`).
#' @param genome Character. BSgenome package name
#'   (e.g. `"BSgenome.Hsapiens.UCSC.hg38"`). Required when `sequence`
#'   is `NULL`.
#' @param sequence Character string of nucleotides spanning the first window.
#'   When provided, `genome` is ignored.
#' @param show_letters Logical. Show nucleotide letters when window width is
#'   <= 80 bp. Default `FALSE`.
#' @param rect_height Numeric in (0, 1] or `NULL`. Fraction of track height
#'   used for the rectangles. Default `NULL` (full track height).
#' @param colors Named character vector mapping nucleotide codes to hex
#'   colors. Defaults to UCSC standard:
#'   `A="#00AA00"`, `T="#FF0000"`, `C="#0000FF"`, `G="#FFB300"`.
#' @param ... Reserved.
#' @return A `SeqSequenceR6` instance.
#' @examples
#' # From BSgenome (requires BSgenome package)
#' \dontrun{
#'   seq_sequence(genome = "BSgenome.Hsapiens.UCSC.hg38")
#' }
#'
#' # From a string
#' seq_sequence(sequence = "ATCGATCGATCG", show_letters = TRUE)
#' @export
seq_sequence <- function(data        = NULL,
                         mapping     = NULL,
                         aesthetics  = aes(),
                         genome      = NULL,
                         sequence    = NULL,
                         show_letters = FALSE,
                         rect_height  = NULL,
                         colors       = NULL,
                         ...) {
  SeqSequenceR6$new(
    data         = data,
    mapping      = mapping,
    aesthetics   = aesthetics,
    genome       = genome,
    sequence     = sequence,
    show_letters = show_letters,
    rect_height  = rect_height,
    colors       = colors,
    ...
  )
}
