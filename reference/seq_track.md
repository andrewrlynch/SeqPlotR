# Create a new seq_track

A track cell is partitioned into five nested zones:

1.  **Track outer margin** — outermost band on the track cell; reserved
    for axis titles derived from the track `mapping`.

2.  **Track inner margin** — separates titles from the window row.

3.  **Window outer margin** — per-window spacer inside the track plot
    region.

4.  **Window inner margin** — per-window band holding axis lines, ticks,
    and tick labels.

5.  **Plot area** — where elements render.

## Usage

``` r
seq_track(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  windows = NULL,
  track_id = NULL,
  direction = "right",
  track_width = 1,
  track_height = 1,
  scale_x = NULL,
  scale_y = NULL,
  scale_x2 = NULL,
  scale_y2 = NULL,
  y_windows = NULL,
  y_windows2 = NULL,
  track_outer_margin = 0.02,
  track_inner_margin = 0.02,
  window_outer_margin = 0,
  window_inner_margin = 0.02,
  window_margin = NULL,
  combine_windows = FALSE,
  combine_y_windows = FALSE,
  flip_x = FALSE,
  flip_y = FALSE,
  elements = list(),
  show_legend = TRUE,
  window_scale = NULL,
  ...
)
```

## Arguments

- data:

  A GRanges object providing the default data for elements in this
  track.

- mapping:

  A SeqMap object from
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) providing
  default aesthetic mappings.

- aesthetics:

  A SeqAes object from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md) providing
  constant aesthetics.

- windows:

  A GRanges object defining the genomic windows for this track. Elements
  cannot define their own windows — this is the only place windows are
  set.

- track_id:

  Character. Unique identifier for this track, used for patchwork layout
  matching and cross-track link references.

- direction:

  One of `"right"` (append to current row) or `"under"` (start new row).
  Ignored for the first track (always top-left) and when
  [`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md)
  is given a layout string.

- track_width:

  Relative width unit within its row. Default 1.

- track_height:

  Relative height unit within its column. Default 1.

- scale_x:

  Optional position scale for the primary x axis — a
  [`seq_scale_genomic()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_genomic.md)
  (default behaviour when `NULL`), `seq_scale_continuous(limits = ...)`
  for scalar data x, or `seq_scale_discrete(levels = ...)` for
  categorical x.

- scale_y:

  Optional position scale for the primary y axis. Pass
  `seq_scale_genomic(y_windows)` to flip the track so genomic position
  runs along y; continuous / discrete scales are supported too.

- scale_x2, scale_y2:

  Optional position scales for the secondary x and y axes. Elements
  routed to the secondary axis via `map(axis.x = 2)` or
  `map(axis.y = 2)` are positioned against these scales. When `NULL`,
  secondary scales are auto-inferred from the contributing elements at
  `layoutGrid()` time.

- y_windows:

  Optional `GRanges` used as the genomic y-axis range when flipping the
  track. Setting this (or a `SeqScaleGenomic` `scale_y`) auto-enables
  `uses_genomic_y`.

- y_windows2:

  Optional `GRanges` for a secondary genomic y axis.

- track_outer_margin:

  Per-track outer margin in npc units. Scalar, length-4
  `c(bottom, left, top, right)` (base-R `par(mar = ...)` order), or a
  named list with any of `top`, `right`, `bottom`, `left`. Default
  `0.02`. Axis titles draw here.

- track_inner_margin:

  Per-track inner margin. Same form as `track_outer_margin`. Default
  `0.02`.

- window_outer_margin:

  Per-window outer margin. Same form. Default `0`. Optional spacer
  around each window.

- window_inner_margin:

  Per-window inner margin. Same form. Default `0.02`. Axis ticks and
  tick labels draw here.

- window_margin:

  Deprecated. Use `aes("window.gap.width" = <value>)` on the plot or
  track aesthetics to control the inter-window gap. Passing a non-`NULL`
  value emits a deprecation warning and is otherwise ignored.

- combine_windows:

  Logical. When `TRUE`, multi-region `windows` are concatenated into a
  single virtual panel — useful for drawing data that spans regions
  (e.g. inter-chromosomal Hi-C contacts) in one continuous track. Each
  original window's axis labels and title are still drawn, separated by
  a thin boundary marker. Default `FALSE` (multi-region windows render
  as separate panels).

- combine_y_windows:

  Logical. Symmetric to `combine_windows` for the genomic y-axis (only
  relevant for tracks with multiple `y_windows`). Default `FALSE`.

- flip_x, flip_y:

  Logical. When `TRUE`, mirror the x or y axis so low data values render
  at the high end of the panel and vice versa. Tick labels follow the
  same orientation.

- elements:

  Optional list of `SeqElement` / `SeqLink` objects to pre-populate this
  track. Elements can also be added via `addElement()` or the `%+%`
  operator.

- show_legend:

  Logical. When `FALSE`, this track contributes no legend keys
  regardless of element `legend` fields. Default `TRUE`.

- window_scale:

  Numeric vector or `NULL`. Per-window x-axis scale factors (e.g. `1e-6`
  for Mb, `1e-3` for kb, `1` for bp). When `NULL` (default), the scale
  is inferred from the narrowest window. Length 1 applies the value to
  all windows; length equal to the number of windows uses values
  positionally; any other length triggers a warning and recycles with
  [`rep_len()`](https://rdrr.io/r/base/rep.html).

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqTrackR6` instance (S3 class `"SeqTrack"`).

## Details

Axis lines sit at the boundary between the plot area and the window
inner margin. Tick labels and axis titles pick up their text from the
track's `mapping`: the x-axis title is the expression assigned to `x`,
and the y-axis title is the expression assigned to `y`.

## Examples

``` r
seq_track(track_id = "A")
#> <SeqTrack>
#>   Public:
#>     addElement: function (elem) 
#>     aesthetics: SeqAes
#>     clone: function (deep = FALSE) 
#>     collect_legend_keys: function () 
#>     combine_windows: FALSE
#>     combine_y_windows: FALSE
#>     data: NULL
#>     direction: right
#>     elements: list
#>     flip_x: FALSE
#>     flip_y: FALSE
#>     has_axis_x2: FALSE
#>     has_axis_y2: FALSE
#>     initialize: function (data = NULL, mapping = NULL, aesthetics = aes(), windows = NULL, 
#>     mapping: NULL
#>     resolved_theme: NULL
#>     scale_x: NULL
#>     scale_x2: NULL
#>     scale_y: NULL
#>     scale_y2: NULL
#>     show_legend: TRUE
#>     track_height: 1
#>     track_id: A
#>     track_inner_margin: list
#>     track_outer_margin: list
#>     track_width: 1
#>     uses_genomic_y: FALSE
#>     window_inner_margin: list
#>     window_margin: NULL
#>     window_outer_margin: list
#>     window_scale: NULL
#>     windows: NULL
#>     y_windows: NULL
#>     y_windows2: NULL
```
