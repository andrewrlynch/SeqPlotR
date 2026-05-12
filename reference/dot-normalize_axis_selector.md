# Normalise `map(axis.x = ...)` / `map(axis.y = ...)` to a scalar integer

Axis selectors are per-element — they route the element to the primary
(1) or secondary (2) axis.
[`.resolve_mapping()`](http://andrewlynch.io/SeqPlotR/reference/dot-resolve_mapping.md)
broadcasts scalars to the data length; collapse back to a single integer
and validate.

## Usage

``` r
.normalize_axis_selector(v, which = "x")
```

## Arguments

- v:

  Raw resolved value (may be `NULL`, vector, or scalar).

- which:

  Either `"x"` or `"y"` — used for the error message.

## Value

Integer scalar in `{1L, 2L}`. Defaults to `1L` when `v` is NULL.
