# Draw points

Primitive element drawing one glyph per observation at `(x, y)`.

## Usage

``` r
seq_point(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`. Falls back to the parent track's data.

- mapping:

  Optional `SeqMap` from
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md). Required
  fields: `x`, `y`. `x` defaults to `start(data)`; `y` defaults to 0.5
  when unmapped. `color`, `fill`, and `shape` are auto-scaled when
  mapped to non-color columns: discrete columns get the flexoki palette;
  numeric columns get a viridis gradient legend.

- aesthetics:

  Optional `SeqAes` from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md). Supports
  `color`, `fill`, `size`, `shape`, `alpha`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqPointR6` instance (S3 class `"SeqPoint" / "SeqElement"`).

## Examples

``` r
seq_point(map(x = start, y = score))
#> <SeqPoint>
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
seq_point(map(x = start, y = score, color = type))  # auto-legend
#> <SeqPoint>
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
