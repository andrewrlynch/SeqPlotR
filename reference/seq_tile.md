# Draw tiles (rectangles or rotated diamonds)

Composite element drawing one filled shape per genomic range. In the
default unrotated mode, each tile is a rectangle spanning `start..end`
on the x-axis at the row indicated by `y` (default 1). With
`aes(rotate = TRUE)` and a `data2` GRanges giving per-tile y-axis
coordinates, tiles are rendered as diamonds via a linear coordinate
transform in genomic space.

## Usage

``` r
seq_tile(data = NULL, mapping = NULL, aesthetics = aes(), data2 = NULL, ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Recognised fields: `x`, `y`, `fill`.

- aesthetics:

  Optional `SeqAes`. `rotate` toggles diamond mode; `fill` sets a
  constant fill color.

- data2:

  Optional `GRanges` providing y-axis genomic ranges for rotated mode.
  Must match `data` in length.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqTileR6` instance.

## Examples

``` r
seq_tile(map(x = start, fill = color))
#> <SeqTile>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: NULL
#>     data2: NULL
#>     draw: function () 
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), data2 = NULL, 
#>     legend: NULL
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
```
