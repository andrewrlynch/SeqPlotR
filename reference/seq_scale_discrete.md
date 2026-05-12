# Discrete position scale

Creates a categorical position scale for axes with discrete levels
(e.g., cell types, sample names).

## Usage

``` r
seq_scale_discrete(
  levels = NULL,
  labels = NULL,
  expand = c(0, 0.5),
  cap = c("capped", "full", "exact", "ticks"),
  oob = c("exclude", "perimeter"),
  pretty = NULL
)
```

## Arguments

- levels:

  Character vector of category levels, in display order. If `NULL`,
  levels are auto-detected from element data.

- labels:

  Optional display labels (same length as `levels`). If `NULL`, level
  names are used as labels.

- expand:

  Length-2 `c(mul, add)` expansion. Default `c(0, 0.5)` (half a category
  of padding on each side).

- cap:

  Axis-line cap mode. See
  [`seq_scale_continuous()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_continuous.md).

## Value

A `SeqScaleDiscrete_Pos` / `SeqPositionScale` object.
