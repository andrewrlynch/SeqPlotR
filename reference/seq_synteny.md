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

- ...:

  Reserved.

## Value

A `SeqSyntenyR6` instance.

## Details

`data` may be a `GRanges` or a `data.frame`. Required
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) fields are
`x0` and `x1`; `x0_end` and `x1_end` default to one base past `x0`/`x1`
when absent. `y0` and `y1` default to the bottom edge of t0's inner
panel and the top edge of t1's inner panel respectively — no scale
conversion is applied in that case.
