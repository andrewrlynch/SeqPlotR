# Hi-C contact matrix

Builds a
[`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md)
rendering a Hi-C contact matrix in one of four styles:

- `"full"`:

  Symmetric square heatmap with genomic position on both x- and y-axes.

- `"diagonal"`:

  Same coordinate system as `"full"` (kept as a separate keyword so call
  sites can switch styles without changing shape).

- `"triangle"`:

  Rotated 45 degrees; upper triangle only, y-axis is interaction
  distance in base pairs.

- `"rectangle"`:

  Same rotation as `"triangle"`, but y-axis is capped at `max_dist`,
  yielding a rectangle.

## Usage

``` r
seq_hic(
  data,
  windows,
  style = "triangle",
  max_dist = NULL,
  palette = "blues",
  na_color = "#FFFFD9",
  y_windows = NULL,
  combine_windows = FALSE,
  combine_y_windows = FALSE,
  flip_x = FALSE,
  flip_y = FALSE,
  track_height = 1,
  track_id = NULL,
  legend = NULL,
  show_legend = TRUE,
  ...
)
```

## Arguments

- data:

  Sparse `GRanges` (mcols: `i_start`, `i_end`, `j_start`, `j_end`,
  `score`) or a numeric matrix / data.frame whose row/column names
  encode bin positions. Each row must describe a well-formed contact
  where `i_end - i_start` and `j_end - j_start` both equal the Hi-C bin
  width (they define the tile's footprint on the position and distance
  axes). For cross-chromosomal contacts, include an optional `j_chrom`
  mcols column giving the j-bin's chromosome (the GRanges's `seqnames`
  is taken to be the i-bin's chromosome). When absent, both bins are
  assumed to live on the same chromosome.

- windows:

  `GRanges` defining the genomic region(s) to display on the x-axis.
  Multiple ranges produce side-by-side panels (one per range), useful
  for comparing several regions.

- style:

  One of `"full"`, `"diagonal"`, `"triangle"`, `"rectangle"`. Default
  `"triangle"`.

- max_dist:

  For `style = "rectangle"` only: cap the distance axis at this value
  (bp). Required for `"rectangle"`.

- palette:

  Colour scale palette for the tile fill. Passed to
  [`seq_scale_color_continuous()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_color_continuous.md).

- na_color:

  Colour for zero/NA contacts.

- y_windows:

  Optional `GRanges` for the genomic y-axis range in `"full"` /
  `"diagonal"` styles. Defaults to `windows` (square matrix). Pass a
  different `GRanges` to display an asymmetric region pair. Ignored for
  rotated styles.

- combine_windows:

  Logical; when `TRUE`, multi-region `windows` are concatenated into a
  single virtual track so cross-window data (e.g. inter-chromosomal
  contacts) renders continuously in one panel. Default `FALSE`.

- combine_y_windows:

  Symmetric to `combine_windows` for multi-region `y_windows` in the
  `full` / `diagonal` styles.

- flip_x, flip_y:

  Logical. Mirror the x or y axis. For the `triangle` style
  `flip_y = TRUE` produces a downward-pointing triangle; for `diagonal`
  it switches to the lower diagonal; for `full` it flips the matrix
  vertically (y) or horizontally (x). Tick labels follow the same
  orientation. Default `FALSE`.

- track_height:

  Relative track height.

- track_id:

  `track_id` for the generated track. Defaults to
  `paste0("hic_", style)`.

- legend:

  A `LegendKey` or `SeqLegendSpec` forwarded to the tile element. `NULL`
  (default) produces no legend entry.

- show_legend:

  Logical. When `FALSE`, the tile element contributes no legend. Default
  `TRUE`.

- ...:

  Additional arguments forwarded to
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

## Value

A `SeqPlot` with a single Hi-C track.

## Details

To show multiple styles side-by-side, call `seq_hic()` multiple times
and combine via
[`seq_resolve()`](http://andrewlynch.io/SeqPlotR/reference/seq_resolve.md)
or `%|%`.

## Examples

``` r
library(GenomicRanges)
set.seed(1)
n  <- 80
st <- sort(sample(seq(1, 1e6, by = 1e4), n))
gr <- GRanges("chr1", IRanges(st, width = 1e4),
              i_start = st, i_end = st + 1e4,
              j_start = st + sample(0:5e5, n, replace = TRUE),
              j_end   = st + sample(0:5e5, n, replace = TRUE) + 1e4,
              score   = rexp(n, rate = 0.5))
win <- GRanges("chr1", IRanges(1, 1e6))
seq_hic(gr, windows = win, style = "triangle")
```
