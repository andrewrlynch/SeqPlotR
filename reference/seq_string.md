# Cubic-Bezier string link between two genomic loci

Draws a smooth Bezier curve connecting two anchors. When
`aes(type = "auto")` (the default), the C vs. S shape is inferred from
the resolved `strand0`/`strand1` fields — opposing strands (`+/-`,
`-/+`) produce an S-curve, matching strands produce a C-curve.

## Usage

``` r
seq_string(
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
  `color` (default `"purple"`), `linewidth` (default `1.5`), `alpha`
  (default `0.8`), `type` (`"auto"`, `"c"`, or `"s"`; default `"auto"`),
  `bulge`, `orientation`.

- ...:

  Reserved.

## Value

A `SeqStringR6` instance.

## Details

`data` may be a `GRanges` (anchor-1 columns in `mcols`) or a
`data.frame` (BEDPE-like). Required
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) fields are
`x0` and `x1`; `chrom0` and `chrom1` are required for `data.frame` input
and optional for `GRanges` (default `seqnames(data)`). `y0` and `y1`
default to `0`.
