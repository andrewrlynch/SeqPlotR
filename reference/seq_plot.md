# Create a new seq_plot

Create a new seq_plot

## Usage

``` r
seq_plot(
  layout = NULL,
  aesthetics = aes(),
  show_legend = TRUE,
  legend = NULL,
  windows = NULL,
  plot_margin = 0.02,
  highlight_margins = FALSE,
  ...
)
```

## Arguments

- layout:

  Either NULL (positional layout, default) or a multiline character
  string defining a patchwork layout. When a layout string is given,
  track positions are determined entirely by `track_id` matching —
  `direction` on
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
  is ignored.

- aesthetics:

  A SeqAes object from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md) for
  plot-wide aesthetics.

- show_legend:

  Logical. Global legend switch — when `FALSE` no legend is drawn
  regardless of element or track settings. Default `TRUE`.

- legend:

  Convenience alias: `legend = FALSE` is sugar for
  `show_legend = FALSE`. All other values are ignored.

- windows:

  Optional plot-wide `GRanges`. Tracks that do not set their own
  `windows` (and genomic scales nested in them) inherit from this value.

- plot_margin:

  Numeric. Outermost canvas margin in NPC units (default `0.02`).
  Per-side overrides via `aesthetics = aes(margins = list(top = …, …))`
  still win.

- highlight_margins:

  Logical. When `TRUE`, overlay translucent coloured rects on every
  margin band — a development aid. Default `FALSE`.

- ...:

  Additional arguments reserved for future use.

## Value

A `SeqPlotR6` instance (S3 class `"SeqPlot"`).

## Examples

``` r
seq_plot()

seq_plot(layout = "AB\nCC")
```
