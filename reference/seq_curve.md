# Draw bezier curves

Primitive element drawing a cubic bezier curve from `(x, y)` to
`(x_end, y_end)` with a y-offset controlled by `curvature`.

## Usage

``` r
seq_curve(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y`, `x_end`, `y_end`.

- aesthetics:

  Optional `SeqAes`. Supports `curvature` (default 0.3), `color`,
  `linewidth`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqCurveR6` instance.

## Examples

``` r
seq_curve(map(x = start, y = score, x_end = end, y_end = score))
#> <SeqCurve>
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
