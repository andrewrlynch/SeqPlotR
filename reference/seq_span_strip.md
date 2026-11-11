# Draw labelled genomic spans

Draws a coloured bar over each genomic span with a text label at the
span centre. Spans are clipped to each window. With `bar_height = 0`
only the labels are drawn.

## Usage

``` r
seq_span_strip(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  legend = NULL,
  legend_extra = NULL,
  show_legend = TRUE,
  legend_position = "canvas_margin",
  legend_side = "right",
  ...
)
```

## Arguments

- data:

  A `GRanges` with one row per span.

- mapping:

  A `SeqMap`. Fields: `x0` (default `start`), `x1` (default `end`),
  `label`, `fill` (concrete colours; `NA` gives `grey70`).

- aesthetics:

  A `SeqAes`. Recognised: `bar_height` (fraction of the panel height,
  default 0.35), `bar_y` (bar centre as a panel fraction, default 0.5),
  `label_cex` (0.6), `label_position` (`"above"` or `"inside"`), `alpha`
  (1).

- legend:

  `NULL` for no key, `FALSE` for none, or a `SeqLegendSpec` or list of
  `LegendKey` used as given.

- legend_extra:

  List of `LegendKey` appended to the legend.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the generated legend.

- ...:

  Unused.

## Value

A `SeqSpanStripR6` instance.

## Examples

``` r
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 500), c(400, 900)),
  name = c("a", "b"), col = c("tomato", "steelblue"))
seq_span_strip(gr, map(label = name, fill = col))
#> <SeqSpanStrip>
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
#>     legend_extra: NULL
#>     legend_position: canvas_margin
#>     legend_side: right
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
#>     y_fields: y y_end y_min y_max y0 y1 height
#>     y_transform: NULL
```
