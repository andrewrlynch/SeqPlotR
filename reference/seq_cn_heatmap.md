# Multi-sample copy-number heatmap

Builds a
[`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md)
with a single track arranging per-sample CN calls as a heatmap: samples
are placed on the y-axis (top to bottom according to `sample_order`),
genomic windows along the x-axis, and each bin is a tile coloured by CN
state.

## Usage

``` r
seq_cn_heatmap(
  data,
  windows,
  sample_col = NULL,
  cn_col = NULL,
  state_colors = NULL,
  sample_order = NULL,
  bins = NULL,
  track_height = 3,
  track_id = NULL,
  ...
)
```

## Arguments

- data:

  A `GRanges` (long format) or numeric matrix.

- windows:

  A `GRanges` of genomic view windows.

- sample_col:

  Name of the mcols column giving sample identity. Auto-detected from
  `c("sample", "sample_id", "Sample", "id")` when `NULL`.

- cn_col:

  Name of the mcols column giving integer CN state. Auto- detected from
  the same candidates as
  [`seq_copynumber()`](http://andrewlynch.io/SeqPlotR/reference/seq_copynumber.md).

- state_colors:

  Named character vector keyed by CN state string. Defaults to the
  [`seq_copynumber()`](http://andrewlynch.io/SeqPlotR/reference/seq_copynumber.md)
  palette.

- sample_order:

  Character vector of sample names in display order (top to bottom).
  When `NULL`, samples are sorted alphabetically.

- bins:

  Optional `GRanges` giving the genomic position of each matrix column;
  required when `data` is a matrix.

- track_height:

  Relative track height.

- track_id:

  Character `track_id` for the generated track.

- ...:

  Additional arguments forwarded to
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

## Value

A `SeqPlot` object.

## Details

Data may be supplied either as a long-format `GRanges` with one row per
(sample, bin) carrying `sample` and `cn` mcols columns, or as a numeric
matrix (rows = samples, cols = bins) plus a `bins` `GRanges` giving the
genomic position of each column (passed via `...`). Matrix input
requires the matrix argument to have sample names as `rownames`.

## Examples

``` r
library(GenomicRanges)
gr <- GRanges("chr1", IRanges(seq(1, 1e6, by = 2e4), width = 2e4),
              sample = rep(c("S1","S2","S3"), length.out = 50),
              cn     = sample(0:4, 50, replace = TRUE))
win <- GRanges("chr1", IRanges(1, 1e6))
seq_cn_heatmap(gr, windows = win)
```
