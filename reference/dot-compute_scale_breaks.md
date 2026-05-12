# Compute breaks, labels, minor breaks, and axis_range for a scale

Given a scale (continuous / genomic / discrete) and the raw data range,
returns a list with the full axis-drawing metadata: `breaks`, `labels`,
`minor_breaks`, `axis_range` (where the axis line should span, governed
by `cap`), `plot_range` (the expanded visible range), and `data_range`.

## Usage

``` r
.compute_scale_breaks(scale, data_range, plot_range = NULL)
```

## Arguments

- scale:

  A `SeqPositionScale`.

- data_range:

  Length-2 numeric vector of the unexpanded data range. For genomic
  scales, this is the window range.

- plot_range:

  Optional pre-expanded plot range. When `NULL` (default), computed from
  `data_range` via
  [`.expand_limits()`](http://andrewlynch.io/SeqPlotR/reference/dot-expand_limits.md)
  using `scale$expand`.

## Value

A list with `breaks`, `labels`, `minor_breaks`, `axis_range`,
`plot_range`, `data_range`.

## Details

Algorithm ported from THEfunc's `make_axis_meta()`; uses base
[`pretty()`](https://rdrr.io/r/base/pretty.html) instead of
`scales::pretty_breaks()` so no new dependency is introduced.
