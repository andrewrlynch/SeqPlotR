# Draw a node-link panel

Draws nodes, directed or plain edges and labelled group boxes at given
positions. No layout is computed: `x` and `y` are used as supplied, in
the track's x and y scales.

## Usage

``` r
seq_node_link(
  nodes = NULL,
  edges = NULL,
  groups = NULL,
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

- nodes, edges, groups:

  Data frames (or data tables) as described above.

- aesthetics:

  A `SeqAes`. Recognised: `node_mm` (3), `bar_linewidth` (3),
  `bar_alpha` (0.35), `arrow_mm` (2), `label_cex` (0.6),
  `edge_label_cex` (0.55), `box_r_mm` (1), `flag_color` (`"#D62728"`).

- legend:

  `NULL` for one key per node shape present, `FALSE` for none, or a
  `SeqLegendSpec` or list of `LegendKey` used as given.

- legend_labels:

  Named character vector relabelling the shape keys.

- legend_extra:

  List of `LegendKey` appended to the keys.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the generated legend.

- ...:

  Unused.

## Value

A `SeqNodeLinkR6` instance. Add it to a track with a `windows` range and
fixed scales.

## Details

Columns of `nodes`: `id`, `x`, `y` (required), `bar_lo`, `bar_hi` (a
translucent vertical bar), `shape` (`circle`, `triangle_up`,
`triangle_down`, `diamond` or `square`; default `circle`), `fill`
(default `grey60`), `alpha` (default 1), `label`.

Columns of `edges`: `from`, `to` (required node ids), `arrow` (default
`TRUE`), `color` (default `grey10`), `linetype` (default 1), `linewidth`
(default 1.4), `label`, `flag` (default `FALSE`; draws a marker at the
edge midpoint).

Columns of `groups`: `x`, `y`, `label` (required), `fill`, `color`.

## Examples

``` r
nodes <- data.frame(id = c("A", "B"), x = c(1, 2), y = c(0.2, 0.6))
edges <- data.frame(from = "A", to = "B", label = "fixed")
seq_node_link(nodes = nodes, edges = edges)
#> <SeqNodeLink>
#>   Inherits from: <SeqElement>
#>   Public:
#>     .infer_scale_y: function () 
#>     aesthetics: SeqAes
#>     auto_legend: NULL
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     coordCanvas: NULL
#>     data: data.frame
#>     draw: function () 
#>     edges: data.frame
#>     groups: NULL
#>     initialize: function (nodes = NULL, edges = NULL, groups = NULL, aesthetics = aes(), 
#>     legend: NULL
#>     legend_arg: NULL
#>     legend_extra: NULL
#>     legend_labels: NULL
#>     legend_position: canvas_margin
#>     legend_side: right
#>     mapping: SeqMap
#>     nodes: data.frame
#>     prep: function (layout_track, track_windows) 
#>     resolve: function (track_data = NULL, track_mapping = NULL) 
#>     resolved: NULL
#>     show_legend: TRUE
#>     y_fields: y y_end y_min y_max y0 y1 height
#>     y_transform: NULL
```
