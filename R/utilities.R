#' Null-coalescing operator
#'
#' Returns `e1` if it is not `NULL`, otherwise returns `e2`.
#'
#' @param e1 Left-hand side value.
#' @param e2 Right-hand side fallback value.
#'
#' @return `e1` if not `NULL`, else `e2`.
#'
#' @examples
#' NULL %||% "fallback"
#' "a" %||% "b"
#'
#' @name op-null-coalesce
#' @usage e1 \%||\% e2
#' @importFrom stats setNames
#' @export
`%||%` <- function(e1, e2) {
  if (!is.null(e1)) e1 else e2
}


#' Assert that an object is a GRanges
#'
#' Internal helper. Errors if `x` does not inherit from `"GRanges"`.
#'
#' @param x Object to check.
#' @param arg_name Name of the argument, used in the error message.
#'
#' @return Invisibly `TRUE` if the check passes; otherwise stops.
#' @keywords internal
.stop_if_not_granges <- function(x, arg_name) {
  if (!inherits(x, "GRanges")) {
    stop("'", arg_name, "' must be a GRanges object.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Normalise a margin specification
#'
#' Accepts a scalar (applied to all four sides), a length-4 numeric in
#' base-R order `c(bottom, left, top, right)`, or an already-normalised
#' list with fields `top`, `right`, `bottom`, `left`. Returns a list in
#' the `list(top, right, bottom, left)` canonical form.
#'
#' @param m A scalar, length-4 numeric, or list.
#' @return A list with fields `top`, `right`, `bottom`, `left`.
#' @keywords internal
.normalize_margin <- function(m) {
  if (is.null(m)) {
    return(list(top = 0, right = 0, bottom = 0, left = 0))
  }
  if (is.list(m)) {
    return(list(
      top    = m$top    %||% 0,
      right  = m$right  %||% 0,
      bottom = m$bottom %||% 0,
      left   = m$left   %||% 0
    ))
  }
  m <- as.numeric(m)
  if (length(m) == 1L) {
    return(list(top = m, right = m, bottom = m, left = m))
  }
  if (length(m) == 4L) {
    # base R par(mar = c(bottom, left, top, right)) order
    return(list(bottom = m[1], left = m[2], top = m[3], right = m[4]))
  }
  stop("margin must be a scalar, a length-4 numeric (bottom, left, top, right), ",
       "or a named list.", call. = FALSE)
}


#' Check that required mcols columns exist on a GRanges
#'
#' Internal helper. Given a GRanges `gr` and a character vector of column names
#' `cols`, errors if any are missing from `mcols(gr)`.
#'
#' @param gr A `GRanges` object.
#' @param cols Character vector of required column names in `mcols(gr)`.
#'
#' @return Invisibly `TRUE` if the check passes; otherwise stops.
#' @keywords internal
.check_mcols <- function(gr, cols) {
  .stop_if_not_granges(gr, "gr")
  have <- colnames(S4Vectors::mcols(gr))
  missing <- setdiff(cols, have)
  if (length(missing) > 0) {
    stop("Missing required mcols column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}


# ── Genome window helpers ────────────────────────────────────────────────────

#' Default whole-genome windows
#'
#' Returns a `GRanges` spanning all hg38 autosomes plus X and Y, one range
#' per chromosome, covering the full chromosome length. Useful as a default
#' `windows` argument for whole-genome ideogram or summary tracks.
#'
#' @param add_chr Logical. If `TRUE` (default), seqnames carry the `"chr"`
#'   prefix (`"chr1"`, `"chr2"`, ...). If `FALSE`, the prefix is stripped
#'   (`"1"`, `"2"`, ...).
#' @return A `GRanges` of length 24 (chr1–22, X, Y).
#' @examples
#' default_genome_windows()
#' default_genome_windows(add_chr = FALSE)
#' @export
default_genome_windows <- function(add_chr = TRUE) {
  chr_lengths <- c(
    "1"  = 248956422, "2"  = 242193529, "3"  = 198295559,
    "4"  = 190214555, "5"  = 181538259, "6"  = 170805979,
    "7"  = 159345973, "8"  = 145138636, "9"  = 138394717,
    "10" = 133797422, "11" = 135086622, "12" = 133275309,
    "13" = 114364328, "14" = 107043718, "15" = 101991189,
    "16" = 90338345,  "17" = 83257441,  "18" = 80373285,
    "19" = 58617616,  "20" = 64444167,  "21" = 46709983,
    "22" = 50818468,  "X"  = 156040895, "Y"  = 57227415
  )
  seqs <- if (isTRUE(add_chr)) paste0("chr", names(chr_lengths))
          else                 names(chr_lengths)
  GenomicRanges::GRanges(
    seqnames = seqs,
    ranges   = IRanges::IRanges(start = 1, end = chr_lengths)
  )
}


#' Build genome windows from region strings
#'
#' Parses a vector of region strings — either `"chr:start-end"` ranges or
#' bare chromosome names — into a `GRanges` of whole-chromosome or
#' sub-chromosomal windows. Bare chromosome names expand to the full length
#' of that chromosome in the supplied genome (via [default_genome_windows()]
#' when `genome = "hg38"`). Overlapping / adjacent regions are merged via
#' [GenomicRanges::reduce()] before return.
#'
#' @param regions Character vector of region strings. Each element is either
#'   `"chr1:1000-2000"` / `"1:1,000-2,000"` (commas tolerated) or a bare
#'   chromosome name (`"chr1"`, `"X"`, etc.). The `"chr"` prefix is
#'   optional on input and controlled on output by `add_chr`.
#' @param padding Integer. Base pairs to extend on each side of every
#'   region. Clipped so `start >= 1`. Default `0`.
#' @param genome Reference genome for bare-chromosome expansion. Only
#'   `"hg38"` is supported.
#' @param add_chr Logical. If `TRUE` (default), output seqnames carry the
#'   `"chr"` prefix.
#' @return A sorted, non-overlapping `GRanges` of windows.
#' @examples
#' create_genome_windows(c("chr1", "chr2:1000000-2000000"))
#' create_genome_windows("chr3", padding = 5e5)
#' @export
create_genome_windows <- function(regions, padding = 0,
                                  genome = "hg38", add_chr = TRUE) {
  if (!identical(genome, "hg38"))
    stop("Only genome = 'hg38' is supported.", call. = FALSE)

  regions <- gsub("chr|Chr|chromosome", "", as.character(regions))
  is_chr_only <- !grepl(":", regions)

  parsed <- utils::strcapture(
    pattern = "^([^:]+):([0-9,]+)-([0-9,]+)$",
    x       = regions[!is_chr_only],
    proto   = list(chrom = character(),
                   start = character(),
                   end   = character())
  )

  if (any(is_chr_only)) {
    chr_df <- data.frame(
      chrom = regions[is_chr_only],
      stringsAsFactors = FALSE
    )
    defaults <- default_genome_windows(add_chr = FALSE)
    def_chrom <- as.character(GenomicRanges::seqnames(defaults))
    m <- match(chr_df$chrom, def_chrom)
    if (any(is.na(m)))
      stop("Unknown chromosome(s): ",
           paste(chr_df$chrom[is.na(m)], collapse = ", "), call. = FALSE)
    chr_df$start <- as.character(BiocGenerics::start(defaults)[m])
    chr_df$end   <- as.character(BiocGenerics::end(defaults)[m])
    parsed <- rbind(parsed, chr_df)
  }

  parsed$start <- as.integer(gsub(",", "", parsed$start))
  parsed$end   <- as.integer(gsub(",", "", parsed$end))

  parsed$start <- parsed$start - padding
  parsed$end   <- parsed$end   + padding
  parsed$start[parsed$start < 1] <- 1

  chrom <- if (isTRUE(add_chr)) paste0("chr", parsed$chrom) else parsed$chrom

  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(start = parsed$start, end = parsed$end)
  )
  BiocGenerics::sort(GenomicRanges::reduce(gr))
}


#' Load a cytoband annotation table
#'
#' Loads UCSC-style cytoband annotations as a `GRanges`. With no arguments,
#' returns the built-in [`cytoband_hg38`] dataset. Given a `path`, reads a
#' TSV cytoband file with columns `chrom chromStart chromEnd name gieStain`
#' (the standard UCSC cytoBand layout).
#'
#' @param path Optional path to a UCSC-style cytoband TSV. If `NULL`
#'   (default), the bundled `cytoband_hg38` table is used.
#' @param as_granges Logical. If `TRUE` (default), return a `GRanges` with
#'   `name` and `gieStain` as mcols — suitable for direct use with
#'   [seq_ideogram()]. If `FALSE`, return the raw data frame.
#' @return A `GRanges` (default) or `data.frame` of cytobands.
#' @examples
#' cb <- load_cytobands()
#' head(cb)
#' @export
load_cytobands <- function(path = NULL, as_granges = TRUE) {
  if (is.null(path)) {
    env <- new.env()
    utils::data("cytoband_hg38", package = "SeqPlotR", envir = env)
    cb <- env$cytoband_hg38
  } else {
    cb <- utils::read.table(
      path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
      col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
    )
  }
  if (!isTRUE(as_granges)) return(cb)

  GenomicRanges::GRanges(
    seqnames = cb$chrom,
    ranges   = IRanges::IRanges(start = cb$chromStart + 1L,
                                end   = cb$chromEnd),
    name     = cb$name,
    gieStain = cb$gieStain
  )
}
