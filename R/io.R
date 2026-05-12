# ── SeqPlotR file connection helpers ─────────────────────────────────────────
#
# open_bigwig() / open_hic() / open_bam() / open_h5() return lightweight S3
# connection objects. Each carries enough metadata to render legends and
# validate inputs but defers the actual data load to a per-call
# $fetch(region) method, with a per-format max_fetch_bp guardrail to prevent
# runaway whole-file reads. $fetch_binned() variants aggregate fetched data
# into fixed-width bins using a configurable function.

# ── Internal binning helpers ───────────────────────────────────────────────────

#' Resolve the `fun` argument to a scalar aggregation function
#' @keywords internal
.resolve_bin_fun <- function(fun) {
  if (is.function(fun)) return(fun)
  switch(fun,
    mean   = function(x) mean(x,   na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE),
    sum    = function(x) sum(x,    na.rm = TRUE),
    stop(sprintf(
      '`fun` must be "mean", "median", "sum", or a function. Got: "%s".',
      fun), call. = FALSE)
  )
}

#' Bin a GRanges signal track into fixed-width genomic bins
#'
#' Assigns each range in `gr` (weighted by overlap fraction) to bins of width
#' `bin_size` covering `region`, then aggregates using `agg_fun`.
#'
#' @param gr      GRanges with a numeric `score` mcols column.
#' @param region  Single-range GRanges defining the window.
#' @param bin_size Integer bin width in bp.
#' @param agg_fun Resolved aggregation function (from `.resolve_bin_fun()`).
#' @return data.frame with columns: seqnames, start, end, score.
#' @keywords internal
.bin_signal_gr <- function(gr, region, bin_size, agg_fun) {
  bin_size  <- as.integer(bin_size)
  win_start <- BiocGenerics::start(region)
  win_end   <- BiocGenerics::end(region)
  chr       <- as.character(GenomicRanges::seqnames(region))

  # Bin breakpoints
  bin_starts <- seq(win_start, win_end, by = bin_size)
  bin_ends   <- pmin(bin_starts + bin_size - 1L, win_end)
  n_bins     <- length(bin_starts)

  scores <- S4Vectors::mcols(gr)$score
  gr_s   <- BiocGenerics::start(gr)
  gr_e   <- BiocGenerics::end(gr)

  bin_vals <- numeric(n_bins)
  for (i in seq_len(n_bins)) {
    bs <- bin_starts[i]; be <- bin_ends[i]
    # Overlap fraction of each signal range with this bin
    ov_start  <- pmax(gr_s, bs)
    ov_end    <- pmin(gr_e, be)
    ov_width  <- pmax(0L, ov_end - ov_start + 1L)
    sig_width <- gr_e - gr_s + 1L
    frac      <- ifelse(sig_width > 0L, ov_width / sig_width, 0)
    contrib   <- scores * frac
    keep      <- ov_width > 0L & !is.na(contrib)
    bin_vals[i] <- if (any(keep)) agg_fun(contrib[keep]) else NA_real_
  }

  data.frame(
    seqnames = chr,
    start    = bin_starts,
    end      = bin_ends,
    score    = bin_vals,
    stringsAsFactors = FALSE
  )
}

#' Rebin a contact data.frame to a coarser bin size
#'
#' Aggregates an existing BEDPE-style contact data.frame (as returned by
#' `$fetch()`) into larger bins using `agg_fun`.
#'
#' @param contacts data.frame with columns seqnames1, start1, end1,
#'   seqnames2, start2, end2, score.
#' @param bin_size Integer. New (coarser) bin size in bp.
#' @param agg_fun  Resolved aggregation function.
#' @return data.frame with the same columns, at the new resolution.
#' @keywords internal
.rebin_contacts <- function(contacts, bin_size, agg_fun) {
  if (nrow(contacts) == 0L) return(contacts)
  bin_size <- as.integer(bin_size)

  # Snap starts to coarser grid
  contacts$bin1 <- (contacts$start1 %/% bin_size) * bin_size
  contacts$bin2 <- (contacts$start2 %/% bin_size) * bin_size

  # Aggregate score per (seqnames1, bin1, seqnames2, bin2)
  key <- paste(contacts$seqnames1, contacts$bin1,
               contacts$seqnames2, contacts$bin2, sep = "\t")
  tab <- split(contacts$score, key)

  agg_scores <- vapply(tab, agg_fun, numeric(1L))
  keys_split <- strsplit(names(agg_scores), "\t")

  data.frame(
    seqnames1 = vapply(keys_split, `[[`, character(1L), 1L),
    start1    = as.integer(vapply(keys_split, `[[`, character(1L), 2L)),
    end1      = as.integer(vapply(keys_split, `[[`, character(1L), 2L)) + bin_size - 1L,
    seqnames2 = vapply(keys_split, `[[`, character(1L), 3L),
    start2    = as.integer(vapply(keys_split, `[[`, character(1L), 4L)),
    end2      = as.integer(vapply(keys_split, `[[`, character(1L), 4L)) + bin_size - 1L,
    score     = agg_scores,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ── open_bigwig ───────────────────────────────────────────────────────────────

#' Open a bigWig file connection
#'
#' Creates a lightweight connection object that probes the bigWig header
#' (sequence info only — no data) and exposes a `$fetch(region)` method
#' returning a `GRanges` of signal values restricted to the requested
#' genomic span.
#'
#' @param path Character. Path to a `.bw` or `.bigwig` file.
#' @param max_fetch_bp Integer. Maximum genomic span allowed per
#'   `$fetch()` call. Default `5e7` (50 Mb). Prevents accidental
#'   whole-genome loads.
#' @return A `SeqBigWig` S3 object with fields `path`, `seqnames`,
#'   `seqinfo`, `max_fetch_bp`, and a `fetch(region)` closure.
#' @export
open_bigwig <- function(path, max_fetch_bp = 5e7L) {
  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)
  if (!requireNamespace("rtracklayer", quietly = TRUE))
    stop("open_bigwig() requires the 'rtracklayer' package.", call. = FALSE)

  bw_info <- tryCatch(
    rtracklayer::seqinfo(rtracklayer::BigWigFile(path)),
    error = function(e) stop("Cannot read bigWig header: ",
                              conditionMessage(e), call. = FALSE)
  )

  obj <- list(
    path         = path,
    max_fetch_bp = as.integer(max_fetch_bp),
    seqinfo      = bw_info,
    seqnames     = GenomeInfoDb::seqnames(bw_info)
  )

  obj$fetch <- function(region) {
    if (!inherits(region, "GRanges") || length(region) < 1L)
      stop("`region` must be a GRanges.", call. = FALSE)
    total_bp <- sum(BiocGenerics::width(region))
    if (total_bp > obj$max_fetch_bp)
      stop(sprintf(
        "Requested region spans %s bp but max_fetch_bp = %s. ",
        format(total_bp, big.mark = ","),
        format(obj$max_fetch_bp, big.mark = ",")),
        "Narrow your windows or increase max_fetch_bp with care.",
        call. = FALSE)
    rtracklayer::import(rtracklayer::BigWigFile(obj$path),
                        which = region, as = "GRanges")
  }

  obj$fetch_binned <- function(region,
                               bin_size = 100L,
                               fun      = "mean") {
    if (!is.numeric(bin_size) || bin_size < 1L)
      stop("`bin_size` must be a positive integer.", call. = FALSE)
    agg_fun <- .resolve_bin_fun(fun)

    gr <- obj$fetch(region)     # uses existing guardrails (max_fetch_bp etc.)
    if (length(gr) == 0L) {
      return(data.frame(seqnames = character(), start = integer(),
                        end = integer(), score = numeric(),
                        stringsAsFactors = FALSE))
    }

    # `region` may span multiple ranges; bin each separately then combine
    results <- lapply(seq_len(length(region)), function(i) {
      .bin_signal_gr(gr, region[i], bin_size = bin_size, agg_fun = agg_fun)
    })
    do.call(rbind, results)
  }

  class(obj) <- c("SeqBigWig", "SeqFileConn")
  obj
}

#' @export
print.SeqBigWig <- function(x, ...) {
  cat(sprintf(
    "<SeqBigWig>  %s\n  sequences: %s\n  max_fetch_bp: %s\n  methods: $fetch(region)  $fetch_binned(region, bin_size, fun)\n",
    x$path,
    paste(utils::head(x$seqnames, 5), collapse = ", "),
    format(x$max_fetch_bp, big.mark = ",")))
  invisible(x)
}

# ── open_bam ──────────────────────────────────────────────────────────────────

#' Open a BAM file connection
#'
#' Creates a lightweight connection object that validates the BAM index,
#' probes the header for sequence names and lengths, and exposes a
#' `$fetch(region, ...)` method returning a `GRanges` of alignments.
#'
#' @param path Character. Path to an indexed BAM file (`.bam`).
#' @param max_fetch_bp Integer. Maximum genomic span per `$fetch()` call.
#'   Default `1e5` (100 kb). BAM fetches over wide regions can be very
#'   slow — keep this tight.
#' @return A `SeqBam` S3 object with `seqnames`, `seq_lengths`,
#'   `max_fetch_bp`, and a `fetch(region, min_mapq, max_reads)` closure.
#' @export
open_bam <- function(path, max_fetch_bp = 1e5L) {
  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)
  idx_paths <- c(paste0(path, ".bai"),
                 sub("\\.bam$", ".bai", path))
  if (!any(file.exists(idx_paths)))
    stop("BAM index not found. Run `samtools index` on ", path,
         call. = FALSE)
  if (!requireNamespace("Rsamtools", quietly = TRUE))
    stop("open_bam() requires 'Rsamtools'.", call. = FALSE)
  if (!requireNamespace("GenomicAlignments", quietly = TRUE))
    stop("open_bam() requires 'GenomicAlignments'.", call. = FALSE)

  hdr <- tryCatch(
    Rsamtools::scanBamHeader(Rsamtools::BamFile(path)),
    error = function(e) stop("Cannot read BAM header: ",
                              conditionMessage(e), call. = FALSE)
  )
  # scanBamHeader() returns the per-file payload directly when called with
  # a BamFile, but wraps it in a path-keyed outer list when called with a
  # path string. Handle both shapes by descending into the first list-typed
  # element until we hit the inner `targets` field.
  sq <- if (is.list(hdr) && !is.null(hdr$targets)) hdr$targets
        else if (length(hdr) > 0L) hdr[[1]]$targets
        else NULL
  if (is.null(sq))
    stop("BAM header missing `targets`: ", path, call. = FALSE)

  obj <- list(
    path         = path,
    max_fetch_bp = as.integer(max_fetch_bp),
    seqnames     = names(sq),
    seq_lengths  = as.integer(sq)
  )

  obj$fetch <- function(region,
                        min_mapq  = 0L,
                        max_reads = 20000L) {
    if (!inherits(region, "GRanges") || length(region) < 1L)
      stop("`region` must be a GRanges.", call. = FALSE)
    total_bp <- sum(BiocGenerics::width(region))
    if (total_bp > obj$max_fetch_bp)
      stop(sprintf(
        "Requested region spans %s bp but max_fetch_bp = %s. ",
        format(total_bp, big.mark = ","),
        format(obj$max_fetch_bp, big.mark = ",")),
        "Narrow your windows or increase max_fetch_bp with care.",
        call. = FALSE)
    bf <- Rsamtools::BamFile(obj$path)
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
    keep <- is.na(mq) | mq >= as.integer(min_mapq)
    gal  <- gal[keep]
    if (length(gal) > max_reads) gal <- gal[seq_len(max_reads)]
    GenomicAlignments::granges(gal)
  }
  class(obj) <- c("SeqBam", "SeqFileConn")
  obj
}

#' @export
print.SeqBam <- function(x, ...) {
  cat(sprintf("<SeqBam>  %s\n  sequences: %s\n  max_fetch_bp: %s\n",
              x$path,
              paste(utils::head(x$seqnames, 5), collapse = ", "),
              format(x$max_fetch_bp, big.mark = ",")))
  invisible(x)
}

# ── open_hic ──────────────────────────────────────────────────────────────────

#' Open a `.hic` file connection (Juicer format)
#'
#' Creates a lightweight connection object that probes the Juicer
#' header (chromosomes and resolutions) and exposes a
#' `$fetch(region, ...)` method returning a `data.frame` of contact
#' pairs. The `region` argument may contain multiple genomic ranges;
#' each range yields its own intra-region contact submatrix and the
#' results are concatenated.
#'
#' @param path Character. Path to a Juicer `.hic` file.
#' @param resolution Integer. Default bin resolution in bp. May be
#'   omitted at construction and supplied to `$fetch()`.
#' @param max_fetch_bp Integer. Maximum genomic span per range in a
#'   single `$fetch()` call. Default `2.8e8` (280 Mb). Each range in
#'   the input GRanges is checked independently.
#' @return A `SeqHic` S3 object with `chromosomes`, `resolutions`,
#'   `max_fetch_bp`, and a `fetch()` closure.
#' @export
open_hic <- function(path, resolution = NULL, max_fetch_bp = 280000000L) {
  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)
  if (!requireNamespace("strawr", quietly = TRUE))
    stop("open_hic() requires the 'strawr' package (CRAN).", call. = FALSE)

  chroms <- tryCatch(
    strawr::readHicChroms(path),
    error = function(e) stop("Cannot read .hic header: ",
                              conditionMessage(e), call. = FALSE)
  )
  resolutions <- tryCatch(
    strawr::readHicBpResolutions(path),
    error = function(e) NULL
  )

  obj <- list(
    path         = path,
    max_fetch_bp = as.integer(max_fetch_bp),
    resolution   = resolution,
    chromosomes  = chroms,
    resolutions  = resolutions
  )

  obj$fetch <- function(region,
                        resolution = obj$resolution,
                        norm       = "NONE",
                        unit       = "BP",
                        matrix     = "observed") {
    if (is.null(resolution))
      stop("A `resolution` must be specified for $fetch().", call. = FALSE)
    if (!inherits(region, "GRanges") || length(region) < 1L)
      stop("`region` must be a GRanges with at least one range.",
           call. = FALSE)

    out <- vector("list", length(region))
    for (i in seq_along(region)) {
      r_i <- region[i]
      w_i <- BiocGenerics::width(r_i)
      if (w_i > obj$max_fetch_bp)
        stop(sprintf(
          "Requested range %d spans %s bp but max_fetch_bp = %s. ",
          i, format(w_i, big.mark = ","),
          format(obj$max_fetch_bp, big.mark = ",")),
          "Narrow your windows or increase max_fetch_bp with care.",
          call. = FALSE)
      chr_i <- as.character(GenomicRanges::seqnames(r_i))
      loc_i <- sprintf("%s:%d:%d", chr_i,
                       BiocGenerics::start(r_i),
                       BiocGenerics::end(r_i))
      raw <- strawr::straw(
        norm    = norm,
        fname   = obj$path,
        chr1loc = loc_i,
        chr2loc = loc_i,
        unit    = unit,
        binsize = as.integer(resolution),
        matrix  = matrix
      )
      if (is.null(raw) || nrow(raw) == 0L) next
      out[[i]] <- data.frame(
        seqnames1 = chr_i,
        start1    = raw$x,
        end1      = raw$x + as.integer(resolution) - 1L,
        seqnames2 = chr_i,
        start2    = raw$y,
        end2      = raw$y + as.integer(resolution) - 1L,
        score     = raw$counts,
        stringsAsFactors = FALSE
      )
    }
    out <- Filter(Negate(is.null), out)
    if (length(out) == 0L) return(data.frame())
    do.call(rbind, out)
  }

  obj$fetch_binned <- function(region,
                               bin_size   = NULL,
                               fun        = "mean",
                               resolution = obj$resolution,
                               norm       = "NONE",
                               unit       = "BP",
                               matrix     = "observed") {

    if (is.null(resolution))
      stop("A `resolution` must be specified.", call. = FALSE)

    contacts <- obj$fetch(
      region     = region,
      resolution = resolution,
      norm       = norm,
      unit       = unit,
      matrix     = matrix
    )

    if (is.null(bin_size) || nrow(contacts) == 0L) return(contacts)

    bin_size <- as.integer(bin_size)
    if (bin_size <= as.integer(resolution))
      stop(sprintf(
        "`bin_size` (%d) must be larger than the native resolution (%d).",
        bin_size, as.integer(resolution)), call. = FALSE)

    agg_fun <- .resolve_bin_fun(fun)
    .rebin_contacts(contacts, bin_size, agg_fun)
  }

  class(obj) <- c("SeqHic", "SeqFileConn")
  obj
}

#' @export
print.SeqHic <- function(x, ...) {
  res_str <- if (!is.null(x$resolutions))
    paste(utils::head(x$resolutions, 5), collapse = ", ")
  else "unknown"
  n_chr <- if (is.data.frame(x$chromosomes)) nrow(x$chromosomes) else
             length(x$chromosomes)
  cat(sprintf(
    "<SeqHic>  %s\n  chromosomes: %d\n  resolutions (bp): %s\n  max_fetch_bp: %s\n  methods: $fetch(...)  $fetch_binned(..., bin_size, fun)\n",
    x$path,
    n_chr,
    res_str,
    format(x$max_fetch_bp, big.mark = ",")))
  invisible(x)
}

# ── open_h5 ───────────────────────────────────────────────────────────────────

#' Open a .cool or .mcool HDF5 contact matrix file
#'
#' Supports both `.cool` (single-resolution) and `.mcool`
#' (multi-resolution) files produced by cooler / HiCExplorer.
#'
#' @param path Character. Path to the `.cool` or `.mcool` file.
#' @param resolution Integer or `NULL`. For `.mcool` files, which resolution
#'   to open. When `NULL`, the coarsest available resolution is used for the
#'   probe; you must set `resolution` explicitly in `$fetch()`.
#' @param max_fetch_bp Integer. Maximum genomic span per `$fetch()` call.
#'   Default `5e6` (5 Mb).
#' @return A `SeqH5` S3 object with `$fetch()` and `$fetch_binned()` methods.
#' @export
open_h5 <- function(path, resolution = NULL, max_fetch_bp = 5e6L) {

  if (!file.exists(path))
    stop("File not found: ", path, call. = FALSE)
  if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("open_h5() requires the 'rhdf5' package (Bioconductor).",
         call. = FALSE)

  # ---- Probe file structure ----
  h5_ls <- tryCatch(
    rhdf5::h5ls(path),
    error = function(e) stop("Cannot read H5 file: ", conditionMessage(e),
                              call. = FALSE)
  )

  # Detect .cool vs .mcool by checking for a 'resolutions' group
  is_mcool    <- any(h5_ls$name == "resolutions" & h5_ls$otype == "H5I_GROUP")
  resolutions <- NULL

  if (is_mcool) {
    res_group   <- h5_ls[h5_ls$group == "/resolutions", ]
    resolutions <- sort(as.integer(res_group$name))
    if (is.null(resolution)) resolution <- max(resolutions)   # coarsest default
    cool_root <- paste0("/resolutions/", as.integer(resolution))
  } else {
    cool_root <- "/"
  }

  # Read chromosome names and lengths from bins table
  chroms <- tryCatch({
    chrom_names   <- rhdf5::h5read(path, paste0(cool_root, "/chroms/name"))
    chrom_lengths <- rhdf5::h5read(path, paste0(cool_root, "/chroms/length"))
    data.frame(name   = as.character(chrom_names),
               length = as.integer(chrom_lengths),
               stringsAsFactors = FALSE)
  }, error = function(e)
    stop("Cannot read chromosome table from H5 file: ", conditionMessage(e),
         call. = FALSE)
  )

  # ---- Helper: fetch raw contacts from cool matrix ----
  .fetch_cool <- function(path, cool_root, chroms,
                          region1, region2, resolution) {

    res <- as.integer(resolution)

    chr1 <- as.character(GenomicRanges::seqnames(region1))
    s1   <- BiocGenerics::start(region1)
    e1   <- BiocGenerics::end(region1)
    chr2 <- as.character(GenomicRanges::seqnames(region2))
    s2   <- BiocGenerics::start(region2)
    e2   <- BiocGenerics::end(region2)

    # Chromosome index (0-based in cooler)
    ci1 <- match(chr1, chroms$name) - 1L
    ci2 <- match(chr2, chroms$name) - 1L
    if (is.na(ci1)) stop("Chromosome '", chr1, "' not in H5 file.",
                         call. = FALSE)
    if (is.na(ci2)) stop("Chromosome '", chr2, "' not in H5 file.",
                         call. = FALSE)

    # Bin indices for each chromosome (cooler bins are 0-based half-open)
    chrom_offsets <- rhdf5::h5read(path,
                                   paste0(cool_root, "/indexes/chrom_offset"))

    bin_offset1 <- chrom_offsets[ci1 + 1L]
    bin_offset2 <- chrom_offsets[ci2 + 1L]

    # Bin numbers within each chromosome
    b1_start <- as.integer((s1 - 1L) %/% res)
    b1_end   <- as.integer((e1 - 1L) %/% res)
    b2_start <- as.integer((s2 - 1L) %/% res)
    b2_end   <- as.integer((e2 - 1L) %/% res)

    # Global bin indices
    g1_start <- bin_offset1 + b1_start
    g1_end   <- bin_offset1 + b1_end
    g2_start <- bin_offset2 + b2_start
    g2_end   <- bin_offset2 + b2_end

    # Read pixels (bin1_id, bin2_id, count) for the region.
    # Cooler stores pixels in bin1_id order; use index to slice.
    bin1_offset <- rhdf5::h5read(path,
                                 paste0(cool_root, "/indexes/bin1_offset"))

    # Row range in pixel table for the bin1 range
    px_start <- bin1_offset[g1_start + 1L] + 1L   # 1-based
    px_end   <- bin1_offset[g1_end   + 2L]         # inclusive

    if (px_end < px_start) return(data.frame())

    px_idx   <- px_start:px_end
    bin2_ids <- rhdf5::h5read(path, paste0(cool_root, "/pixels/bin2_id"),
                              index = list(px_idx))
    counts   <- rhdf5::h5read(path, paste0(cool_root, "/pixels/count"),
                              index = list(px_idx))
    bin1_ids <- rhdf5::h5read(path, paste0(cool_root, "/pixels/bin1_id"),
                              index = list(px_idx))

    # Filter to bin2 range
    keep <- bin2_ids >= g2_start & bin2_ids <= g2_end
    if (!any(keep)) return(data.frame())

    bin1_ids <- bin1_ids[keep]
    bin2_ids <- bin2_ids[keep]
    counts   <- counts[keep]

    # Convert global bin indices back to genomic coords
    local1 <- bin1_ids - bin_offset1
    local2 <- bin2_ids - bin_offset2

    data.frame(
      seqnames1 = chr1,
      start1    = as.integer(local1 * res + 1L),
      end1      = as.integer(local1 * res + res),
      seqnames2 = chr2,
      start2    = as.integer(local2 * res + 1L),
      end2      = as.integer(local2 * res + res),
      score     = as.numeric(counts),
      stringsAsFactors = FALSE
    )
  }

  # ---- Build connection object ----
  obj <- list(
    path         = path,
    is_mcool     = is_mcool,
    cool_root    = cool_root,
    resolution   = resolution,
    resolutions  = resolutions,
    chromosomes  = chroms,
    max_fetch_bp = as.integer(max_fetch_bp)
  )

  obj$fetch <- function(region1, region2 = NULL,
                        resolution = obj$resolution,
                        ...) {

    if (is.null(resolution))
      stop("A `resolution` must be specified for $fetch().", call. = FALSE)
    if (!inherits(region1, "GRanges") || length(region1) != 1L)
      stop("`region1` must be a single-range GRanges.", call. = FALSE)
    if (is.null(region2)) region2 <- region1

    total_bp <- max(BiocGenerics::width(region1),
                    BiocGenerics::width(region2))
    if (total_bp > obj$max_fetch_bp)
      stop(sprintf(
        "Requested region spans %s bp but max_fetch_bp = %s. ",
        format(total_bp, big.mark = ","),
        format(obj$max_fetch_bp, big.mark = ",")),
        "Narrow your windows or increase max_fetch_bp with care.",
        call. = FALSE)

    # Resolve correct cool_root for this resolution in mcool
    cr <- if (obj$is_mcool)
      paste0("/resolutions/", as.integer(resolution))
    else
      obj$cool_root

    .fetch_cool(obj$path, cr, obj$chromosomes, region1, region2, resolution)
  }

  obj$fetch_binned <- function(region1, region2 = NULL,
                               bin_size   = NULL,
                               fun        = "mean",
                               resolution = obj$resolution,
                               ...) {

    if (is.null(resolution))
      stop("A `resolution` must be specified.", call. = FALSE)

    contacts <- obj$fetch(region1 = region1, region2 = region2,
                          resolution = resolution, ...)

    if (is.null(bin_size) || nrow(contacts) == 0L) return(contacts)

    bin_size <- as.integer(bin_size)
    if (bin_size <= as.integer(resolution))
      stop(sprintf(
        "`bin_size` (%d) must be larger than the native resolution (%d).",
        bin_size, as.integer(resolution)), call. = FALSE)

    agg_fun <- .resolve_bin_fun(fun)
    .rebin_contacts(contacts, bin_size, agg_fun)
  }

  class(obj) <- c("SeqH5", "SeqFileConn")
  obj
}

#' @export
print.SeqH5 <- function(x, ...) {
  res_str <- if (!is.null(x$resolutions))
    paste(utils::head(x$resolutions, 5), collapse = ", ")
  else
    as.character(x$resolution)
  cat(sprintf(
    "<%s>  %s\n  chromosomes: %d\n  resolution(s) (bp): %s\n  max_fetch_bp: %s\n  methods: $fetch(region, resolution)  $fetch_binned(..., bin_size, fun)\n",
    if (x$is_mcool) "SeqH5 [mcool]" else "SeqH5 [cool]",
    x$path,
    nrow(x$chromosomes),
    res_str,
    format(x$max_fetch_bp, big.mark = ",")))
  invisible(x)
}

#' Test whether `x` is a SeqPlotR file connection
#'
#' Connection objects from [open_bigwig()], [open_bam()], [open_hic()],
#' and [open_h5()] all carry the `SeqFileConn` class.
#'
#' @param x An object.
#' @return `TRUE` if `x` inherits from `SeqFileConn`, otherwise `FALSE`.
#' @export
is_seq_file_conn <- function(x) inherits(x, "SeqFileConn")
