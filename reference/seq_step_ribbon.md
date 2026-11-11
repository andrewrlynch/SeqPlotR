# Draw step ribbons with interval bars and direction triangles

Draws each row as a thin horizontal band over its genomic span at height
`y`, with an optional vertical interval line from `bar_lo` to `bar_hi`
at the span midpoint, a fixed-size triangle at the midpoint pointing up
or down, and an optional label right of the span. Band thickness scales
with the mapped `thickness`; triangle size is in mm, independent of span
width and axis range. A row overlapping several windows is drawn in
each, clipped to the window.

## Usage

``` r
seq_step_ribbon(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  legend = NULL,
  legend_labels = NULL,
  legend_extra = NULL,
  show_legend = TRUE,
  legend_position = "canvas_margin",
  legend_side = "right",
  ...
)
```

## Arguments

- data:

  Optional `GRanges`. Falls back to the parent track's data.

- mapping:

  Optional `SeqMap`. Fields: `x0` (default `start`), `x1` (default
  `end`), `y` (required), `thickness` (default 1), `color` (concrete
  colours, default `"grey30"`), `alpha` (default 1), `label`, `bar_lo`,
  `bar_hi`, `pointing` (`"up"`, `"down"` or `NA`). `y`, `bar_lo` and
  `bar_hi` are in y axis units.

- aesthetics:

  Optional `SeqAes`. `half_height` (axis units per thickness unit,
  default 0.006), `tri_half_mm` (1.1), `tri_height_mm` (1.9),
  `line_width` (1), `label_cex` (0.55), `label_mm` (1).

- legend:

  `NULL` for automatic keys (one line key per thickness value present,
  labelled `"thickness <k>"`), `FALSE` for none, or a `SeqLegendSpec` or
  list of `LegendKey` used as given.

- legend_labels:

  Named character vector relabelling the automatic keys, named by
  thickness value (e.g. `c("1" = "one copy")`).

- legend_extra:

  List of `LegendKey` appended to the legend keys.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the legend. Default `"canvas_margin"`, `"right"`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqStepRibbonR6` instance.

## Examples

``` r
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1e5, 5e5), width = 2e5),
                             t = c(0.3, 0.6), k = c(1L, 3L))
seq_step_ribbon(gr, map(y = t, thickness = k, pointing = "up"))
#> <SeqStepRibbon>
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
#>     legend_labels: NULL
#>     legend_position: canvas_margin
#>     legend_side: right
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
#>     y_fields: y bar_lo bar_hi
#>     y_transform: NULL
```
