# Load a cytoband annotation table

Loads UCSC-style cytoband annotations as a `GRanges`. With no arguments,
returns the built-in
[`cytoband_hg38`](http://andrewlynch.io/SeqPlotR/reference/cytoband_hg38.md)
dataset. Given a `path`, reads a TSV cytoband file with columns
`chrom chromStart chromEnd name gieStain` (the standard UCSC cytoBand
layout).

## Usage

``` r
load_cytobands(path = NULL, as_granges = TRUE)
```

## Arguments

- path:

  Optional path to a UCSC-style cytoband TSV. If `NULL` (default), the
  bundled `cytoband_hg38` table is used.

- as_granges:

  Logical. If `TRUE` (default), return a `GRanges` with `name` and
  `gieStain` as mcols — suitable for direct use with
  [`seq_ideogram()`](http://andrewlynch.io/SeqPlotR/reference/seq_ideogram.md).
  If `FALSE`, return the raw data frame.

## Value

A `GRanges` (default) or `data.frame` of cytobands.

## Examples

``` r
cb <- load_cytobands()
head(cb)
#> GRanges object with 6 ranges and 2 metadata columns:
#>       seqnames            ranges strand |        name    gieStain
#>          <Rle>         <IRanges>  <Rle> | <character> <character>
#>   [1]     chr1         1-2300000      * |      p36.33        gneg
#>   [2]     chr1   2300001-5300000      * |      p36.32      gpos25
#>   [3]     chr1   5300001-7100000      * |      p36.31        gneg
#>   [4]     chr1   7100001-9100000      * |      p36.23      gpos25
#>   [5]     chr1  9100001-12500000      * |      p36.22        gneg
#>   [6]     chr1 12500001-15900000      * |      p36.21      gpos50
#>   -------
#>   seqinfo: 711 sequences from an unspecified genome; no seqlengths
```
