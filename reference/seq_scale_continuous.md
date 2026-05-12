# Continuous position scale

Creates a numeric position scale for axes displaying continuous data.

## Usage

``` r
seq_scale_continuous(
  limits = NULL,
  n_breaks = 5,
  breaks = NULL,
  minor_breaks = NULL,
  expand = c(0.025, 0),
  cap = c("capped", "full", "exact", "ticks"),
  labels = NULL,
  oob = c("exclude", "perimeter"),
  pretty = NULL
)
```

## Arguments

- limits:

  Optional numeric vector of length 2 clamping the axis range. If
  `NULL`, the range is auto-computed from element data.

- n_breaks:

  Target number of pretty breaks (default 5). Ignored when `breaks` is
  supplied.

- breaks:

  Optional numeric vector of explicit break positions.

- minor_breaks:

  Optional scalar (number of sub-divisions between major breaks) or
  numeric vector of explicit minor-break positions.

- expand:

  Length-2 `c(mul, add)` specifying padding around the data range. `mul`
  is multiplied by `diff(data_range)`; `add` is added in data units.
  Default `c(0.025, 0)` — a 2.5% breath.

- cap:

  One of `"capped"` (default — axis line spans the break range, so tick
  labels don't look stranded), `"full"` (spans the expanded plot range),
  `"exact"` (spans the unexpanded data range), or `"ticks"` (suppress
  the axis line entirely).

- labels:

  Optional character vector of tick labels (same length as `breaks`).
  When `NULL`, breaks are formatted automatically.

## Value

A `SeqScaleContinuous_Pos` / `SeqPositionScale` object.
