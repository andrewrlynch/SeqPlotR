# Create a legend placement specification

Constructs a `SeqLegendSpec` object that describes where and how a
legend group should be rendered. The spec captures placement intent only
— no rendering logic is executed here. It is consumed by the
legend-rendering layer in a later batch.

## Usage

``` r
seq_legend(
  keys,
  title = NULL,
  position = "inside",
  x = 0.5,
  y = 0.5,
  hjust = 0,
  orientation = NULL,
  nrow = NULL,
  ncol = NULL,
  side = NULL
)
```

## Arguments

- keys:

  A single `LegendKey`, or a named or unnamed list of `LegendKey`
  objects. Required.

- title:

  Optional character string. The legend group title drawn above or
  beside the key block.

- position:

  Character. One of `"inside"`, `"track_margin"`, or `"canvas_margin"`.
  Controls which area the legend occupies. Default `"inside"`.

- x:

  Numeric in \[0, 1\]. Horizontal position of the legend anchor within
  the target area. Default `0.5`.

- y:

  Numeric in \[0, 1\]. Vertical position of the legend anchor within the
  target area. Default `0.5`.

- hjust:

  Numeric in \[0, 1\]. Horizontal justification of the legend content
  block relative to the anchor. `0` = left-aligned, `1` = right-aligned.
  Default `0`.

- orientation:

  Character. One of `"horizontal"` or `"vertical"`. Controls whether
  keys are laid out as a row or a column. When `NULL` (default),
  orientation is inferred from `side`: `"vertical"` for
  `side %in% c("left", "right")`, `"horizontal"` otherwise.

- nrow:

  Integer or `NULL`. Number of rows in the key grid. When both `nrow`
  and `ncol` are `NULL`, defaults to 1 row (all keys in one row for
  horizontal orientation) or 1 column (all keys in one column for
  vertical orientation).

- ncol:

  Integer or `NULL`. Number of columns in the key grid.

- side:

  Character or `NULL`. For
  `position %in% c("track_margin", "canvas_margin")`, which margin to
  target. One of `"top"`, `"bottom"`, `"left"`, `"right"`. When `NULL`,
  defaults to `"top"` for margin positions and is left `NULL` for
  `"inside"`.

## Value

A `SeqLegendSpec` S3 object (a named list with class `"SeqLegendSpec"`).

## Examples

``` r
k1 <- LegendKey(label = "H3K27ac", color = "firebrick")
k2 <- LegendKey(label = "H3K4me3", color = "steelblue")

# Inside the track, top-left
spec <- seq_legend(list(k1, k2), title = "Marks", x = 0.02, y = 0.95)

# In the outer track margin, bottom
spec <- seq_legend(list(k1, k2), position = "track_margin", side = "bottom")
```
