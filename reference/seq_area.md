# Draw filled area under a curve

Primitive element drawing a filled polygon formed by the data line
`(x, y)` and a return path along `baseline`.

## Usage

``` r
seq_area(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y`.

- aesthetics:

  Optional `SeqAes`. Supports `baseline` (default 0), `fill`, `color`,
  `alpha`, `linewidth`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqAreaR6` instance.

## Examples

``` r
seq_area(map(x = mid, y = score))
#> <SeqArea>
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
