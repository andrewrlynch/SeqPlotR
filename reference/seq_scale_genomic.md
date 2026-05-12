# Genomic position scale

Creates a position scale based on genomic coordinates from a `GRanges`
object. Used for axes that represent physical genome positions (bp, kb,
Mb).

## Usage

``` r
seq_scale_genomic(
  windows = NULL,
  scale_factor = NULL,
  breaks = NULL,
  minor_breaks = NULL,
  expand = c(0, 0),
  cap = c("capped", "full", "exact", "ticks"),
  labels = NULL,
  oob = c("exclude", "perimeter"),
  pretty = NULL
)
```

## Arguments

- windows:

  A `GRanges` object defining the genomic windows.

- scale_factor:

  Numeric vector of per-window scale factors controlling the unit label
  (1e-6 = Mb, 1e-3 = kb, 1 = bp). If `NULL`, reads from
  `mcols(windows)$scale` or defaults to `1e-6`.

- breaks:

  Optional numeric vector of explicit break positions. When supplied,
  these are used instead of
  [`pretty()`](https://rdrr.io/r/base/pretty.html)-generated breaks.

- minor_breaks:

  Optional scalar (number of sub-divisions between major breaks) or
  numeric vector of explicit minor-break positions.

- expand:

  Length-2 `c(mul, add)` specifying multiplicative and additive padding
  around the data range. Default `c(0, 0)` for genomic (windows already
  set the visible range).

- cap:

  One of `"capped"` (axis line spans the break range), `"full"` (spans
  the expanded plot range), `"exact"` (spans the unexpanded data range),
  or `"ticks"` (no axis line, only ticks).

- labels:

  Optional character vector of tick labels (same length as `breaks`). If
  `NULL`, breaks are formatted as decimal numbers.

- oob:

  One of `"exclude"` (default — out-of-bounds data rows are dropped
  before the npc transform; a
  [`message()`](https://rdrr.io/r/base/message.html) reports the count)
  or `"perimeter"` (clamped to the axis limits and counted).

- pretty:

  Optional named list of arguments forwarded to base
  [`pretty()`](https://rdrr.io/r/base/pretty.html) (e.g.
  `list(min.n = 2, high.u.bias = 0)`). Names that are not base
  [`pretty()`](https://rdrr.io/r/base/pretty.html) arguments (such as
  `bounds` or `f.min`) will error from R.

## Value

A `SeqScaleGenomic` / `SeqPositionScale` object.
