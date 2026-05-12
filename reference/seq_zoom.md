# Zoom / highlight polygon between two tracks

Projects a genomic region onto two tracks and connects the two
projections with a filled quadrilateral. Typically used to link an
overview track to a zoomed detail track (e.g. ideogram → locus). Always
used as a plot-level link — `t0` and `t1` must be specified explicitly.

## Usage

``` r
seq_zoom(
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

  Track identifiers (required).

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `fill` (default `"grey50"`), `alpha` (default `0.15`), `color`
  (border; default `NA`), `linewidth` (default `0.5`), `stemOffset`
  (default `0.01`).

- ...:

  Reserved.

## Value

A `SeqZoomR6` instance.

## Details

Map vocabulary:

- `x0`, `x0_end`:

  required — genomic edges of the region in `t0`.

- `x1`, `x1_end`:

  optional — genomic edges of the region in `t1`; default to `x0`,
  `x0_end` (same region projected onto both tracks).

- `chrom0`, `chrom1`:

  optional for GRanges (default `seqnames(data)`); required for
  `data.frame`.

Attachment edges auto-select so that the upper track attaches from its
bottom and the lower track from its top.
