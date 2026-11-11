# Arches whose stems rise from one track to tops on another

Draws, per row, vertical stems at one or two genomic anchors starting at
`y0` on the track the element is added to, joined by a Bezier arch whose
tops sit at `top` on the `target` track. Both tracks must use the same
windows: x is shared and only y maps across tracks, through the target
track's panels (and its y transform, e.g.
[`seq_scale_count()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_count.md)).
Rows without a second anchor draw a single stem with a pennant. Optional
interval glyphs (see
[`seq_interval_glyph()`](http://andrewlynch.io/SeqPlotR/reference/seq_interval_glyph.md))
are drawn on both stems at the top, and an optional marker point on each
stem at `marker` in source track units. An anchor outside every window
is an error naming the row.

## Usage

``` r
seq_cross_arch(
  data = NULL,
  mapping = NULL,
  target,
  aesthetics = aes(),
  legend = NULL,
  legend_extra = NULL,
  show_legend = TRUE,
  legend_position = "canvas_margin",
  legend_side = "right",
  ...
)
```

## Arguments

- data:

  `GRanges` or `data.frame` with one row per arch.

- mapping:

  `SeqMap` from
  [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md). Fields:
  `x0` (default `start`), `chrom0` (default `seqnames` for GRanges),
  `x1` and `chrom1` (second anchor; `NA` or unmapped gives a single
  stem; `chrom1` defaults to `chrom0`), `top` (arch top in target track
  units; `NA` peaks at the target axis top), `y0` (stem base, default
  0), `color`, `alpha`, `glyph`, `bar_lo`, `bar_hi` (glyph type and bar
  in target track units), `marker` (point on each stem), `caution`
  (logical; draws a caution mark) and `id` (row label used in error
  messages).

- target:

  `track_id` of the track holding the arch tops.

- aesthetics:

  `SeqAes` from
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `linewidth` (1.2), `curve` (0.04 npc), `color`, `alpha`,
  `bar_linewidth` (4), `bar_alpha` (1), `glyph_points` (`TRUE`), `size`
  (glyph point size in mm, 2), `cap_mm` (2.5), `marker_size` (0.8),
  `marker_shape` (16), `single_mm` (pennant width, 1.6), `caution_mm`
  (1.6).

- legend:

  `NULL` (no keys of its own), `FALSE`, a `SeqLegendSpec` or a list of
  `LegendKey` objects used as given.

- legend_extra:

  Optional `LegendKey` or list of them appended.

- show_legend:

  Logical. `FALSE` suppresses the legend.

- legend_position, legend_side:

  Placement of the built legend.

- ...:

  Passed to the generator (e.g. `t0` for a plot-level link).

## Value

A `SeqCrossArch` link element; add it to the source track.

## Details

`coordCanvas` holds one list per row with canvas npc values
`row, x0, y0, x1, y1, top0, top1, single, color, lwd, glyph, lo0, hi0, lo1, hi1, ytop, marker0, marker1, caution`.

## Examples

``` r
library(GenomicRanges)
win <- GRanges("chr1", IRanges(1, 1e6))
sv <- data.frame(chr = "chr1", pos1 = 2e5, pos2 = 6e5, t = 0.4)
p <- seq_plot() %+%
  seq_track(track_id = "top", windows = win,
            scale_y = seq_scale_continuous(limits = c(0, 1))) %+%
  seq_track(track_id = "base", windows = win,
            scale_y = seq_scale_continuous(limits = c(0, 4))) %+%
  seq_cross_arch(sv, map(chrom0 = chr, x0 = pos1, x1 = pos2, top = t),
                 target = "top")
```
