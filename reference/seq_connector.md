# Bezier connector between interval edges

Draws a cubic Bezier from the left or right edge of one interval to the
left or right edge of another, with handles pointing away from each
edge. Intervals may sit on the same track or on different tracks. Use it
to join consecutive pieces of an assembly, a rearranged chromosome or a
read path.

## Usage

``` r
seq_connector(
  data = NULL,
  mapping = NULL,
  t0 = NULL,
  t1 = NULL,
  aesthetics = aes(),
  styles = NULL,
  height = 0.5,
  handle = 0.04,
  bulge_height = 0.45,
  layer = c("above", "below"),
  legend = NULL,
  show_legend = TRUE,
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
  `color`, `linetype`, `linewidth` for unstyled connectors.

- styles:

  Optional named list of `list(col, lty, lwd)`.

- height:

  Block height as a fraction of the panel inner height (default `0.5`).

- handle:

  Handle length as a fraction of the track inner width (default `0.04`;
  same-track handles use 30 percent of it).

- bulge_height:

  Same-track bulge as a fraction of the panel inner height (default
  `0.45`).

- layer:

  `"above"` (default) or `"below"` the track elements.

- legend:

  Optional `LegendKey`, list of keys or
  [`seq_legend()`](http://andrewlynch.io/SeqPlotR/reference/seq_legend.md).

- show_legend:

  Logical. When `style` is mapped, an automatic legend lists the keys
  present.

- ...:

  Reserved.

## Value

A `SeqConnectorR6` instance.

## Details

Same-track connectors leave from the top (or bottom) edge of blocks of
`height` centred in the panel and bulge up (or down); connectors between
tracks leave from the block centres with horizontal handles.

[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) fields: `x0`,
`x1` (required: the left edges of interval 0 and interval 1, in track
coordinates), `x0_end`, `x1_end` (right edges; default the left edges),
`chrom0`, `chrom1` (window seqnames; required for `data.frame` data),
`side0`, `side1` (`"left"` or `"right"`; default `"right"` and
`"left"`), `style` (a key into `styles`) and `bulge` (same-track only;
`TRUE` or `NA` up, `FALSE` down). Numeric edges are used exactly.

Default styles: `junction` solid blue, `inferred` red dashed, `unused`
grey. Entries of `styles` override these by name. A missing or unknown
key is drawn with `aes(color, linetype, linewidth)` (defaults
`"#3B82D6"`, `1`, `1.3`).

## Examples

``` r
w <- GenomicRanges::GRanges("chrA", IRanges::IRanges(0, 100))
df <- data.frame(chrom0 = "chrA", x0 = 10, x0_end = 30,
                 chrom1 = "chrA", x1 = 60, x1_end = 90)
p <- seq_plot() %+%
  seq_track(track_id = "a", windows = w) %+%
  seq_connector(df, map(x0 = x0, x0_end = x0_end, chrom0 = chrom0,
                        x1 = x1, x1_end = x1_end, chrom1 = chrom1),
                t0 = "a", t1 = "a")
```
