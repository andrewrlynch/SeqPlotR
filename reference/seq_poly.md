# Draw filled polygons

Primitive element drawing one or more filled polygons. Vertices come
from `(x, y)`; an optional `group` mapping partitions vertices across
multiple polygons.

## Usage

``` r
seq_poly(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required fields: `x`, `y`. Optional: `group`.

- aesthetics:

  Optional `SeqAes`. Supports `fill`, `color`, `alpha`, `linewidth`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqPolyR6` instance.

## Examples

``` r
seq_poly(map(x = start, y = score))
#> <SeqPoly>
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
