# Convert genomic x and data y to canvas npc coordinates for a single panel

Internal helper used by link prep methods. Clamps to the panel's `inner`
rectangle so any value outside the data ranges lands on the boundary.

## Usage

``` r
.to_canvas(x_gen, y_data, panel_meta)
```

## Arguments

- x_gen:

  Numeric vector of genomic x positions.

- y_data:

  Numeric vector of data-scale y values (same length as `x_gen`).

- panel_meta:

  Panel metadata list (with `xscale`, `yscale`, and `inner` containing
  `x0`, `x1`, `y0`, `y1`).

## Value

List with `x` and `y` numeric vectors of canvas npc coordinates.
