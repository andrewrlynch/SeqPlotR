# Plot a derivative-chromosome reconstruction

Draws a `plot_derivatives` table set (schema version 1), as written by
the ROSE package, with
[`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md),
[`seq_connector()`](http://andrewlynch.io/SeqPlotR/reference/seq_connector.md),
[`seq_synteny()`](http://andrewlynch.io/SeqPlotR/reference/seq_synteny.md),
[`seq_track_label()`](http://andrewlynch.io/SeqPlotR/reference/seq_track_label.md)
and
[`seq_caption()`](http://andrewlynch.io/SeqPlotR/reference/seq_caption.md).
The adapter only selects, renames and maps table columns to element
arguments. One track is drawn per row: a reference row, derivative rows,
deleted copies and candidate strips. All rows share one base-pair scale.
Clusters with `cn_panels` also get copy number tracks above the rows.

## Usage

``` r
seq_rose_derivatives(
  x,
  fill = c("cytoband", "solid"),
  palette = NULL,
  block_colors = NULL,
  synteny = FALSE,
  label_inside = FALSE,
  flip_style = c("rotate", "overline", "minus"),
  show_loose = TRUE,
  cytobands = NULL,
  annotations = NULL,
  ...
)
```

## Arguments

- x:

  A directory holding `plot_derivatives.*.tsv` files, or a list of data
  frames with a `schema_version`. See
  [`seq_read_table_set()`](http://andrewlynch.io/SeqPlotR/reference/seq_read_table_set.md).

- fill:

  `"cytoband"` (default) or `"solid"`.

- palette:

  `NULL` (default hues), `"soft"`, or a character vector of colours,
  optionally named by `fill_key`.

- block_colors:

  Optional named colours overriding block fills by label, then by key.

- synteny:

  Logical. Draw ribbons between the reference and derivative rows.

- label_inside:

  Logical. Place block labels inside blocks when they fit.

- flip_style:

  How flipped blocks are marked: `"rotate"`, `"overline"` or `"minus"`.

- show_loose:

  Logical. Draw the loose derivative rows and the rows that continue
  them.

- cytobands:

  Cytoband table passed to
  [`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md).
  `NULL` uses the bundled hg38 bands.

- annotations:

  Optional data frame for
  [`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md)
  annotations, with an optional `track` column naming the row label to
  draw on. Rows without `track` go on the reference and derivative rows.

- ...:

  Arguments merged over the
  [`seq_plot()`](http://andrewlynch.io/SeqPlotR/reference/seq_plot.md)
  defaults.

## Value

A `SeqPlot`, extendable with `%+%`.
