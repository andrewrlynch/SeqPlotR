# Draw a line

Primitive element drawing a polyline through `(x, y)` for each interval.
`x` defaults to the interval midpoint; `y` defaults to 0.5.

## Usage

``` r
seq_line(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Recognized fields: `x`, `y`.

- aesthetics:

  Optional `SeqAes`. Recognized fields: `color`, `linewidth`,
  `linetype`, `alpha`, and `type` (`"step"` or `"s"` for step lines).

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqLineR6` instance.

## Examples

``` r
seq_line(map(x = mid, y = score))
#> <SeqLine>
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
