# Two-channel colour key specification

Describes a smooth colour square over `[0, 1] x [0, 1]` for
[`seq_tile_mix()`](http://andrewlynch.io/SeqPlotR/reference/seq_tile_mix.md),
with channel `a` along x and channel `b` along y. Titles, tick positions
and tick labels are supplied by the caller.

## Usage

``` r
seq_legend_2d(
  color_a,
  color_b,
  x_title = NULL,
  y_title = NULL,
  x_breaks = c(0, 1),
  x_labels = NULL,
  y_breaks = c(0, 1),
  y_labels = NULL,
  note = NULL,
  size = 0.85,
  resolution = 128,
  position = c("canvas_margin", "track_margin", "inside"),
  side = "right",
  x = 0.08,
  y = 0.5
)
```

## Arguments

- color_a, color_b:

  Colours of the two channels.

- x_title, y_title:

  Axis titles, or `NULL`.

- x_breaks, y_breaks:

  Tick positions in `[0, 1]`.

- x_labels, y_labels:

  Tick labels. Default to the breaks.

- note:

  Optional text drawn below the key.

- size:

  Side length of the key in inches.

- resolution:

  Number of raster cells per side.

- position:

  One of `"canvas_margin"`, `"track_margin"`, `"inside"`.

- side:

  Margin side for margin positions.

- x:

  Horizontal anchor as a fraction of the target width.

- y:

  Vertical anchor of the key top as a fraction of the target height. In
  a canvas margin the key is placed below the other legend keys of the
  same side when they exist, and `y` is used otherwise.

## Value

A `SeqLegend2DSpec` object.

## Examples

``` r
seq_legend_2d("#D62728", "#2B59C3", x_title = "A", y_title = "B")
#> $color_a
#> [1] "#D62728"
#> 
#> $color_b
#> [1] "#2B59C3"
#> 
#> $x_title
#> [1] "A"
#> 
#> $y_title
#> [1] "B"
#> 
#> $x_breaks
#> [1] 0 1
#> 
#> $x_labels
#> [1] "0" "1"
#> 
#> $y_breaks
#> [1] 0 1
#> 
#> $y_labels
#> [1] "0" "1"
#> 
#> $note
#> NULL
#> 
#> $size
#> [1] 0.85
#> 
#> $resolution
#> [1] 128
#> 
#> $position
#> [1] "canvas_margin"
#> 
#> $side
#> [1] "right"
#> 
#> $x
#> [1] 0.08
#> 
#> $y
#> [1] 0.5
#> 
#> attr(,"class")
#> [1] "SeqLegend2DSpec"
```
