# Draw a filled ribbon between two y series

Composite element filling the band between `y_min` and `y_max` at each
`x` value. Useful for confidence intervals and uncertainty bands.

## Usage

``` r
seq_ribbon(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y_min`, `y_max`.

- aesthetics:

  Optional `SeqAes`. Supports `fill`, `color`, `alpha`, `linewidth`.
  Defaults to `aes(fill = "grey60", alpha = 0.8)`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqRibbonR6` instance.

## Examples

``` r
seq_ribbon(map(x = start, y_min = lo, y_max = hi))
#> <SeqRibbon>
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
