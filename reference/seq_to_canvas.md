# Map data coordinates into a panel's canvas npc coordinates

Maps genomic `x` and data `y` values into one panel of a track, the same
way the built-in elements do: `x` through the panel's x plot range and
`y` through its y plot range (falling back to the scale limits), into
the panel's `inner` rectangle. Use it in `prep()` of custom elements and
links, including links that place values on another track via
`layout_all_tracks[[track_id]]`.

## Usage

``` r
seq_to_canvas(x, y, panel, apply_transform = TRUE, clamp = TRUE)
```

## Arguments

- x:

  Numeric vector of genomic (or x-axis) positions.

- y:

  Numeric vector of y values in the track's axis units.

- panel:

  One panel-metadata list, e.g. `layout_track[[w]]` in an element's
  `prep()` or `layout_all_tracks[[track_id]][[w]]` in a link's.

- apply_transform:

  Logical. Apply the panel's `y_transform` (set by transformed scales
  such as
  [`seq_scale_count()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_count.md))
  to `y` first. Use `FALSE` for values that were already transformed by
  `resolve()`.

- clamp:

  Logical. Clamp the result to the panel's inner rectangle.

## Value

A list with numeric vectors `x` and `y` in canvas npc units. `NA` inputs
stay `NA`.

## Examples

``` r
panel <- list(xscale = c(0, 100), yscale = c(0, 1),
              inner = list(x0 = 0.1, x1 = 0.9, y0 = 0.2, y1 = 0.8))
seq_to_canvas(c(25, 50), c(0, 0.5), panel)
#> $x
#> [1] 0.3 0.5
#> 
#> $y
#> [1] 0.2 0.5
#> 
```
