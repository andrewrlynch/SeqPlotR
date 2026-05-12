# Draw a kernel density estimate as a filled area

Composite element that computes
[`stats::density()`](https://rdrr.io/r/stats/density.html) of the
resolved `y` values and renders the distribution as a filled area. The
density evaluation axis (the distribution of y values) is drawn
horizontally across the panel using the track's `yscale`; densities are
normalised to the panel height.

## Usage

``` r
seq_density(data = NULL, mapping = NULL, aesthetics = aes(), ...)
```

## Arguments

- data:

  Optional `GRanges`.

- mapping:

  Optional `SeqMap`. Required field: `y`.

- aesthetics:

  Optional `SeqAes`. Supports `fill`, `color`, `alpha`, `linewidth`, and
  `bw` (kernel bandwidth, default `"nrd0"`).

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqDensityR6` instance.

## Examples

``` r
seq_density(map(y = score))
#> <SeqDensity>
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
