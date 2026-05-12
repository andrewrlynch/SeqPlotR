# Chromosome ideogram from cytogenetic bands

Draws a chromosome ideogram: each cytoband becomes a filled rectangle
shaded by its Giemsa stain (`gpos25` darkens through `gpos100`; `gneg`
is white; `stalk` and `gvar` carry their conventional colors). Paired
`acen` bands within a window render as two inward-pointing red
triangles, marking the centromere.

## Usage

``` r
seq_ideogram(
  data = NULL,
  mapping = NULL,
  aesthetics = aes(),
  scope = "window",
  style = "block",
  highlight_range = NULL,
  ...
)
```

## Arguments

- data:

  Optional `GRanges` of cytobands. Falls back to the parent track's
  data. Must carry a `gieStain` mcol unless a `stain` mapping is
  supplied.

- mapping:

  Optional [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).
  Recognised: `stain`.

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `color` (band border color, default `"black"`), `linewidth` (band
  border width, default `0.1`), `outline` (nested
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md) controlling
  the chromosome's outer perimeter outline: sub-keys `col`, `lwd`,
  `visible`), `highlight` (nested
  [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md) controlling
  the `scope = "full"` highlight rectangle: sub-keys `fill`, `col`,
  `lwd`, `alpha`), `telomere.radius` (numeric; corner radius as a
  fraction of band height for `style = "rounded"`; `1.0` = full
  half-circle cap).

- scope:

  Character. One of `"window"` (default — only bands overlapping the
  track windows are drawn) or `"full"` (the whole chromosome is drawn
  rescaled to fill the panel; the current window region is overlaid as a
  translucent highlight rectangle).

- style:

  Character. One of `"block"` (default — rectangular bands) or
  `"rounded"` (rounded telomere caps on the leftmost and rightmost
  bands).

- highlight_range:

  Optional `GRanges`. Only honoured when `scope = "full"`. When set, the
  highlight rectangle marks this range instead of the parent track's
  `windows`. Lets you set the track's `windows` to span the full
  chromosome (so the x-axis reads chromosome coordinates) while still
  highlighting a sub-range.

- ...:

  Reserved.

## Value

A `SeqIdeogramR6` instance.

## Details

The simplest call supplies a `GRanges` of cytobands — use
[`load_cytobands()`](http://andrewlynch.io/SeqPlotR/reference/load_cytobands.md)
to load the bundled hg38 table:

    cb <- load_cytobands()
    seq_plot() %|%
      seq_track(track_id = "Ideo",
                windows = default_genome_windows()) %+%
      seq_ideogram(data = cb)

Map a non-standard stain column with `map(stain = my_col)`.
