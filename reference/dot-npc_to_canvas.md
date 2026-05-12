# Map unit-interval panel coordinates to canvas npc coordinates

Uses `panel_meta$inner` as the destination rectangle.

## Usage

``` r
.npc_to_canvas(u, v, panel_meta)
```

## Arguments

- u:

  Numeric vector of x in `[0, 1]` panel coordinates.

- v:

  Numeric vector of y in `[0, 1]` panel coordinates.

- panel_meta:

  A panel metadata list with element `inner` containing `x0`, `x1`,
  `y0`, `y1`.

## Value

A list with `x` and `y` numeric vectors of canvas npc coordinates.
