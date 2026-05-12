# Create a color gradient (color-bar) legend specification

Constructs a `GradientLegendSpec` describing a continuous color-scale
legend (color bar). When `breaks` is supplied the bar is rendered as
discrete `LegendKey` entries instead — useful for heatmaps that need
only a handful of labelled stops.

## Usage

``` r
seq_gradient_legend(
  palette = "viridis",
  limits = c(0, 1),
  title = NULL,
  position = "inside",
  x = 0.5,
  y = 0.5,
  hjust = 0,
  orientation = NULL,
  side = NULL,
  breaks = NULL
)
```

## Arguments

- palette:

  One of `"viridis"`, `"plasma"`, `"magma"`, `"blues"`, `"reds"`.
  Default `"viridis"`.

- limits:

  Numeric vector of length 2. The data range the gradient spans. Default
  `c(0, 1)`.

- title:

  Optional character string. Legend group title.

- position:

  One of `"inside"`, `"track_margin"`, `"canvas_margin"`.

- x, y:

  Numeric in \[0, 1\]. Anchor of the bar within the target area. Default
  `0.5`.

- hjust:

  Horizontal justification of the bar relative to the anchor. Default
  `0`.

- orientation:

  `"horizontal"` or `"vertical"`. Bar direction. Inferred from `side`
  when `NULL`.

- side:

  For margin positions: one of `"top"`, `"bottom"`, `"left"`, `"right"`.
  Defaults to `"top"` for margin positions.

- breaks:

  `NULL` (continuous color bar), a single positive integer (that many
  evenly-spaced discrete keys), or a numeric vector of explicit break
  values.

## Value

A `GradientLegendSpec` S3 object.

## Examples

``` r
# Continuous color bar
seq_gradient_legend(palette = "viridis", limits = c(0, 100), title = "Score")
#> <GradientLegendSpec> palette=viridis  limits=[0, 100]  position=inside  side=NULL  breaks=NULL (continuous)  title="Score"

# Five discrete keys
seq_gradient_legend(palette = "plasma", limits = c(-2, 2), breaks = 5)
#> <GradientLegendSpec> palette=plasma  limits=[-2, 2]  position=inside  side=NULL  breaks=n=5 discrete
```
