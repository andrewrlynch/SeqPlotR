# Draw two-channel colour tiles

Draws one rectangle per row over `[x0, x1] x [y_min, y_max]`. Two
intensities in `[0, 1]`, `a` and `b`, set the colour: the hue is
interpolated between `color_a` and `color_b` by `b / (a + b)` and the
depth from white is `max(a, b)`. Rows with both intensities zero or
missing are not drawn. Tiles carry a hairline border in their own colour
so adjacent tiles show no seams. Draw bound lines with
[`seq_segment()`](http://andrewlynch.io/SeqPlotR/reference/seq_segment.md).

## Usage

``` r
seq_tile_mix(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  legend = NULL,
  show_legend = TRUE,
  ...
)
```

## Arguments

- data:

  A `GRanges` with one row per tile.

- mapping:

  A `SeqMap`. Fields: `x0` (default `start`), `x1` (default `end`),
  `y_min`, `y_max`, `a`, `b`.

- aesthetics:

  A `SeqAes`. Recognised: `color_a` (default `"#D62728"`), `color_b`
  (default `"#2B59C3"`), `seam` (border line width, default `0.25`).

- legend:

  `NULL` for the default 2D key from
  [`seq_legend_2d()`](http://andrewlynch.io/SeqPlotR/reference/seq_legend_2d.md),
  `FALSE` for no key, or a `SeqLegend2DSpec` used as given.

- show_legend:

  Logical. `FALSE` suppresses the key.

- ...:

  Unused.

## Value

A `SeqTileMixR6` instance.

## Examples

``` r
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 101), c(100, 200)),
  lo = c(0, 0.2), hi = c(0.2, 0.4), g = c(1, 0.5), l = c(0, 0.5))
seq_tile_mix(gr, map(y_min = lo, y_max = hi, a = g, b = l))
#> <SeqTileMix>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: GRanges
#>     draw: function () 
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), legend = NULL, 
#>     legend: NULL
#>     legend_arg: NULL
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
#>     y_fields: y_min y_max
#>     y_transform: NULL
```
