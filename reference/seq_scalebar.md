# Reference scalebar element

Draws a horizontal scalebar of a specified genomic length with tick caps
and a label. Useful as a visual length reference when the x-axis is
hidden.

## Usage

``` r
seq_scalebar(
  length_bp,
  label = NULL,
  hjust = 0.95,
  vjust = 0.05,
  ticks = TRUE,
  tick_height = 0.15,
  bar_lwd = 1.2,
  bar_col = "#1C1B1A",
  label_cex = 0.65,
  label_col = "#1C1B1A",
  label_pad = 0.005
)
```

## Arguments

- length_bp:

  Numeric. Length of the scalebar in base pairs. Required.

- label:

  Character or `NULL`. Bar label. When `NULL`, auto- formatted from
  `length_bp` (e.g. `"50 kb"`).

- hjust:

  Numeric in 0–1. Horizontal position of the bar's right edge within the
  panel. Default `0.95`.

- vjust:

  Numeric in 0–1. Vertical position of the bar within the panel. Default
  `0.05` (near the bottom).

- ticks:

  Logical. Draw tick caps at each end. Default `TRUE`.

- tick_height:

  Numeric. Tick height as a fraction of panel height. Default `0.15`.

- bar_lwd:

  Numeric. Line width of the bar and ticks. Default `1.2`.

- bar_col:

  Character. Bar / tick color. Default `"#1C1B1A"`.

- label_cex:

  Numeric. Label character expansion. Default `0.65`.

- label_col:

  Character. Label color. Default `"#1C1B1A"`.

- label_pad:

  Numeric. NPC gap between bar top and label bottom. Default `0.005`.

## Value

A `SeqScalebarR6` instance.
