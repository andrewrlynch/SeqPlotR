# Draw lollipops

Composite element drawing a vertical stem from `baseline` (default 0) to
`y`, with a point placed at `y`.

## Usage

``` r
seq_lollipop(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Recognised fields: `x`, `y`.

- aesthetics:

  Optional `SeqAes`. Supports `baseline`, `color`, `linewidth`, `size`,
  `shape`, `alpha`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqLollipopR6` instance.

## Examples

``` r
seq_lollipop(map(x = start, y = score))
#> <SeqLollipop>
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
