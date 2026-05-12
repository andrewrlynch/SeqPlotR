# Apply an out-of-bounds policy to (x, y) coordinates

Filters or clamps `(x, y)` values that fall outside the (expanded) plot
range. The mode is set per-scale via `seq_scale_*(oob = ...)`.

## Usage

``` r
.apply_oob(
  x,
  y = NULL,
  plot_range_x,
  plot_range_y = NULL,
  mode = c("exclude", "perimeter"),
  label = ""
)
```

## Arguments

- x, y:

  Numeric vectors of equal length. `y` may be `NULL` for elements with
  only an x dimension.

- plot_range_x, plot_range_y:

  Length-2 numeric vectors of the expanded plot range. `plot_range_y`
  may be `NULL`.

- mode:

  One of `"exclude"` (drop OOB rows; default) or `"perimeter"` (clamp to
  the limit).

- label:

  Optional character used in the emitted message.

## Value

A list with components:

- `x`,`y`:

  the filtered or clamped coordinates

- `keep`:

  logical mask of which input rows survived (always `TRUE` for
  `perimeter`)

- `n_oob`:

  integer count of rows that were out of bounds
