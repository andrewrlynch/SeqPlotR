# ── Package datasets ─────────────────────────────────────────────────────────

#' UCSC hg38 cytoband table
#'
#' Cytogenetic band annotations for the hg38 human reference assembly,
#' pulled from the UCSC cytoBand track. Used by [seq_ideogram()] and
#' [load_cytobands()].
#'
#' @format A data frame with 1,549 rows and 5 columns:
#' \describe{
#'   \item{chrom}{Chromosome name (UCSC-style, with `"chr"` prefix).}
#'   \item{chromStart}{0-based start coordinate.}
#'   \item{chromEnd}{End coordinate (half-open).}
#'   \item{name}{Cytogenetic band name (e.g. `"p36.33"`).}
#'   \item{gieStain}{Giemsa stain intensity code (`"gneg"`, `"gpos25"`,
#'     `"gpos50"`, `"gpos75"`, `"gpos100"`, `"acen"`, `"stalk"`, `"gvar"`).}
#' }
#' @source UCSC Genome Browser — <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz>
"cytoband_hg38"
