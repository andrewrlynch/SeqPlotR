# Draw bars

Composite element drawing one filled rectangle per genomic range,
spanning `start..end` on the x-axis with height `y`. Supply a `group`
mapping to stack bars vertically at identical x positions.

## Usage

``` r
seq_bar(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`. Falls back to the parent track's data.

- mapping:

  Optional `SeqMap` from
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md). Recognised
  fields: `x`, `y` (bar height, default 1), `group` (stacking), `fill`
  (per-bar color).

- aesthetics:

  Optional `SeqAes` from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md). Supports
  `fill`, `color`, `linewidth`, `alpha`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqBarR6` instance.

## Examples

``` r
seq_bar(map(x = start, y = score))
#> <SeqBar>
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
