# Filled trapezoid connecting two tracks

Draws a filled quadrilateral whose bottom edge spans a region in track
`t0` and whose top edge spans a homologous region in track `t1`. Use to
highlight syntenic blocks, orthologous segments, or otherwise paired
intervals across two stacked tracks.

## Usage

``` r
seq_synteny(
  data = NULL,
  mapping = NULL,
  t0 = NULL,
  t1 = NULL,
  aesthetics = aes(),
  height = NULL,
  smooth = FALSE,
  layer = c("above", "below"),
  ...
)
```

## Arguments

- data:

  Optional `GRanges` or `data.frame`.

- mapping:

  Optional [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).

- t0, t1:

  Track identifiers. Locked to the parent track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `fill` / `color` (fill color), `alpha` (default `0.4`), `linewidth`
  (border width).

- height:

  Optional block height as a fraction of the panel inner height. `NULL`
  (default) anchors at the panel edges.

- smooth:

  Logical. Smoothstep sides instead of straight ones.

- layer:

  `"above"` (default) draws after track elements; `"below"` draws before
  them, so ribbons sit beneath blocks.

- ...:

  Reserved.

## Value

A `SeqSyntenyR6` instance.

## Details

`data` may be a `GRanges` or a `data.frame`. Required
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) fields are
`x0` and `x1`; `x0_end` and `x1_end` default to one base past `x0`/`x1`
when absent. Numeric edges are used exactly (no rounding to whole
bases). `y0` and `y1` default to the bottom edge of t0's inner panel and
the top edge of t1's inner panel respectively — no scale conversion is
applied in that case. With `height` set, the defaults become the facing
edges of blocks of that height centred in each panel (the same geometry
as
[`seq_block()`](http://andrewlynch.io/SeqPlotR/reference/seq_block.md)
and
[`seq_connector()`](http://andrewlynch.io/SeqPlotR/reference/seq_connector.md)).

The optional logical
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) field
`twisted` pairs `x0` with `x1_end` and `x0_end` with `x1`, so the two
sides cross (an inverted copy). With `smooth = TRUE` each side follows a
smoothstep curve (16 points per side), which turns twisted ribbons into
S shapes.
