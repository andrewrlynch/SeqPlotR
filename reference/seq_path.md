# Draw connected paths

Primitive element drawing a polyline connecting observations in order,
optionally partitioned into separate paths by a grouping column.

## Usage

``` r
seq_path(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y`. Optional: `group`
  (column whose values split observations into separate paths).

- aesthetics:

  Optional `SeqAes`. Supports `color`, `linewidth`, `linetype`, `alpha`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqPathR6` instance.

## Examples

``` r
seq_path(map(x = mid, y = score))
#> <SeqPath>
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
