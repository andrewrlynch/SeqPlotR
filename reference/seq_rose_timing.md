# Draw a ROSE timing plot table set

Thin adapter that turns a `plot_timing` table set written by ROSE into a
[`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md)
built from the generic timing elements
([`seq_cross_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_cross_arch.md),
[`seq_step_ribbon()`](http://andrewlynch.io/SeqPlotR/reference/seq_step_ribbon.md),
[`seq_paired_ribbon()`](http://andrewlynch.io/SeqPlotR/reference/seq_paired_ribbon.md),
[`seq_interval_glyph()`](http://andrewlynch.io/SeqPlotR/reference/seq_interval_glyph.md),
[`seq_tile_mix()`](http://andrewlynch.io/SeqPlotR/reference/seq_tile_mix.md),
[`seq_node_link()`](http://andrewlynch.io/SeqPlotR/reference/seq_node_link.md),
[`seq_span_strip()`](http://andrewlynch.io/SeqPlotR/reference/seq_span_strip.md)).
The tables carry every position, interval, opacity, label and colour
category, so the adapter only selects rows, renames columns, looks up
palette colours and sets axes. Proxy rows are order bounds, not times.
Glyph rows without a position, such as a whole-genome doubling, are
drawn across every window at their given `y`.

## Usage

``` r
seq_rose_timing(
  x,
  layout = c("arches", "order", "heatmap", "poset"),
  poset_view = c("poset", "hasse", "groups"),
  palette = NULL,
  show_hidden = FALSE,
  plain = FALSE,
  cytobands = NULL,
  legend_width = 0.2,
  ...
)
```

## Arguments

- x:

  A directory holding the `plot_timing.*` files, or a list of data
  frames with a `schema_version` element or attribute. Read with
  [`seq_read_table_set()`](http://andrewlynch.io/SeqPlotR/reference/seq_read_table_set.md).

- layout:

  One of `"arches"` (SV arches over copy number with step ribbons),
  `"order"` (arches at their poset layer), `"heatmap"` (arches over a
  two-channel timing heatmap) or `"poset"` (node-link view).

- poset_view:

  For `layout = "poset"`: `"poset"` (all edges), `"hasse"` (reduction
  edges only) or `"groups"` (group boxes only).

- palette:

  Named list merged over the default palettes. Entries are `sv`, `step`,
  `cn`, `sbs`, `heatmap`, `glyph` and `node` (named colour vectors),
  `annotation` (named colours by `colour_key`) and `edge` (named lists
  of `color`, `linetype` and `linewidth`).

- show_hidden:

  Logical. Add a strip for glyph rows marked `hidden`.

- plain:

  Logical. Draw a genomic x axis on the last timing track instead of an
  ideogram.

- cytobands:

  Cytoband `GRanges` for the ideogram. Defaults to
  [`load_cytobands()`](http://andrewlynch.io/SeqPlotR/reference/load_cytobands.md).

- legend_width:

  Right margin of the page, as a fraction, left for the legends.

- ...:

  Passed to
  [`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md).

## Value

A `SeqPlot` that can be extended with `%+%`.

## Examples

``` r
d <- system.file("extdata", "rose", "timing", "sim_bfb_amplicon", package = "SeqPlotR")
if (nzchar(d)) {
  p <- seq_rose_timing(d, layout = "arches", plain = TRUE)
}
```
