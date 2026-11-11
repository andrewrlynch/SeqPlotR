# Draw paired ribbons for two values per span

Draws two thin bands per genomic span: one from the first value `y_up`
upward and one from the second value `y_down` downward, so equal values
touch rather than overlap. All second-value bands are drawn before all
first-value bands. Band thickness is a fraction of the top of the
track's y range, so it is the same in axis units for every row.

## Usage

``` r
seq_paired_ribbon(
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
  `end`), `y_up`, `y_down` (y axis units).

- aesthetics:

  Optional `SeqAes`. `fill_up` (default `"#D62728"`), `fill_down`
  (`"#1F77B4"`), `thickness` (fraction of the y range top, default
  0.012), `alpha`.

- legend:

  `NULL` for two automatic rect keys (`"first value"`,
  `"second value"`), `FALSE` for none, or a `SeqLegendSpec` or list of
  `LegendKey` used as given.

- legend_labels:

  Character relabelling the automatic keys, named `up` and `down` or
  given in that order.

- legend_extra:

  List of `LegendKey` appended to the legend keys.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the legend. Default `"canvas_margin"`, `"right"`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqPairedRibbonR6` instance.

## Examples

``` r
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 5e5), width = 5e5),
                             a = c(2, 3), b = c(1, 1))
seq_paired_ribbon(gr, map(y_up = a, y_down = b),
                  legend_labels = c(up = "A", down = "B"))
#> <SeqPairedRibbon>
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
#>     legend_labels: A B
#>     legend_position: canvas_margin
#>     legend_side: right
#>     mapping: SeqMap
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
#>     y_fields: y_up y_down
#>     y_transform: NULL
```
