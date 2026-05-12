# Draw text labels

Primitive element drawing a text label at each observation's `(x, y)`.

## Usage

``` r
seq_text(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y`, `label`.

- aesthetics:

  Optional `SeqAes`. Supports `size`, `color`, `angle`, `hjust`,
  `vjust`, and a constant `label` when not supplied via
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqTextR6` instance.

## Examples

``` r
seq_text(map(x = start, y = score, label = name))
#> <SeqText>
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
