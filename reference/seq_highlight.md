# Cross-track region highlight band

Draws a continuous filled highlight band that passes through every track
from `t0` down to `t1` (inclusive), with per-track widths determined by
each track's own genomic scale. Adjacent tracks are bridged by
trapezoids in the inter-track gap, so the band fans / compresses
smoothly across tracks with different windows. Useful for ChIP-seq style
stacks where a single locus or region of interest should be highlighted
across many panels at once, or for an overview track stacked above a
zoomed detail.

## Usage

``` r
seq_highlight(
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

  Optional `GRanges` or `data.frame` of one or more highlight regions.

- mapping:

  Optional [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).

- t0, t1:

  Track identifiers (string `track_id` or integer index) bracketing the
  inclusive run of tracks. `t0` is required; `t1` defaults to `t0`
  (single-track highlight).

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `fill` (default `"grey50"`), `alpha` (default `0.25`), `color`
  (border; default `NA`), `linewidth` (default `0.5`).

- ...:

  Reserved.

## Value

A `SeqHighlightR6` instance.

## Details

Always used as a plot-level link — `t0` must be specified explicitly.
Set `t1 = NULL` (or omit it) to highlight a single track only.

Map vocabulary:

- `x0`, `x0_end`:

  required — genomic edges of the highlight region. The same genomic
  region is projected onto every track in `[t0..t1]` and naturally
  compresses or expands per-track based on each track's window range.

- `chrom0`:

  optional for GRanges (default `seqnames(data)`); required for
  `data.frame`.
