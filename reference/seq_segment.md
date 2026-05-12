# Draw line segments

Primitive element drawing a straight segment from `(x, y)` to
`(x_end, y_end)` for each observation. When `map(color = col)` is used
with a non-color column, colors are auto-scaled (discrete → flexoki
palette; numeric → viridis gradient) and an auto-legend is generated.

## Usage

``` r
seq_segment(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `x_end`, `y`, `y_end`.
  Defaults: `x = start`, `x_end = end`, `y0 = 0`, `y1 = resolved y`.
  `color` is auto-scaled when mapped to a non-color column.

- aesthetics:

  Optional `SeqAes`. Supports `color`, `linewidth`, `linetype`, `alpha`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqSegmentR6` instance.

## Examples

``` r
seq_segment(map(x = start, x_end = end, y = score, y_end = score))
#> <SeqSegment>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: NULL
#>     draw: function () 
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), ...) 
#>     legend: NULL
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
```
