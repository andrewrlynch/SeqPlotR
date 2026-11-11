# Draw interval glyphs: point estimates, one-sided bounds and intervals

Draws one glyph per row at genomic position `x` and value `y`, with a
vertical bar from `bar_lo` to `bar_hi`. The glyph type is read per row
from the `glyph` mapping:

## Usage

``` r
seq_interval_glyph(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  legend = NULL,
  legend_labels = NULL,
  legend_all = FALSE,
  legend_extra = NULL,
  show_legend = TRUE,
  legend_position = "canvas_margin",
  legend_side = "right",
  ...
)
```

## Arguments

- data:

  `GRanges`. Falls back to the parent track's data.

- mapping:

  `SeqMap` from
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md). Fields:
  `x` (default `start`), `y`, `bar_lo`, `bar_hi`, `glyph` (default
  `"point"`), `color`, `alpha`, `label` (printed right of the glyph when
  not `NA`) and `caution` (logical; draws a small caution triangle).

- aesthetics:

  `SeqAes` from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md): `size`
  (point size in mm, default 2), `bar_linewidth` (4), `bar_alpha` (1;
  multiplies the row alpha for bars), `cap_mm` (2.5), `glyph_points`
  (`TRUE`), `label_cex` (0.6), `label_mm` (1.5), `caution_mm` (1.6),
  `color`, `alpha`.

- legend:

  `NULL` (keys for the glyph types present, from
  [`seq_glyph_legend_keys()`](http://andrewlynch.io/SeqPlotR/reference/seq_glyph_legend_keys.md)),
  `FALSE` (no legend), a `SeqLegendSpec`, or a list of `LegendKey`
  objects used as given.

- legend_labels:

  Named character relabelling the auto keys by glyph type, e.g.
  `c(point = "estimate")`.

- legend_all:

  Logical. Show keys for all five drawn glyph types.

- legend_extra:

  Optional `LegendKey` or list of them appended to the legend.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the built legend. Default `"canvas_margin"`, `"right"`.

- ...:

  Unused.

## Value

A `SeqIntervalGlyph` element.

## Details

- `"point"`: filled point at `y` and a solid bar from `bar_lo` to
  `bar_hi`.

- `"latest"`: cap at `y` and an open bar from `y` down to `bar_lo`.

- `"earliest"`: cap at `y` and an open bar from `y` up to `bar_hi`.

- `"interval"`: bracket from `bar_lo` to `bar_hi`, plus a point at `y`
  when `y` is not `NA`.

- `"untimed"`: nothing is drawn.

- `"hidden"`: the row is dropped.

A row with `NA` `y` gets no point and no cap. All values are drawn as
given. `coordCanvas[[w]]` holds one data.frame per window with columns
`row, glyph, x, y, lo, hi, color, alpha, label, caution` (canvas npc).

## Examples

``` r
library(GenomicRanges)
gr <- GRanges("chr1", IRanges(c(2e5, 5e5, 8e5), width = 1),
              y = c(0.4, 0.6, NA), lo = c(0.3, 0.2, 0.5),
              hi = c(0.5, 0.6, 0.9),
              g = c("point", "latest", "interval"))
el <- seq_interval_glyph(gr, map(y = y, bar_lo = lo, bar_hi = hi,
                                 glyph = g))
```
