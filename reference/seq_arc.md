# Bezier arch link between two genomic loci on the same track

Single Bezier arch with no stems. Anchors live on the same track: when
added inside a
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
via `%+%`, both `t0` and `t1` are locked to the parent track.

## Usage

``` r
seq_arc(
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

  Optional `GRanges` or `data.frame`. Falls back to the parent track's
  data.

- mapping:

  Optional [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md).

- t0, t1:

  Track identifiers. Locked to the parent track when added inside a
  [`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md).

- aesthetics:

  Optional [`aes()`](http://andrewlynch.io/SeqPlotR/reference/aes.md):
  `color`, `linewidth`, `orientation` (`"+"` / `"-"` / `"*"`), `curve`
  (`"length"`, `"equal"`, or numeric).

- ...:

  Reserved.

## Value

A `SeqArcR6` instance.

## Details

`data` may be a `GRanges` (with anchor-1 columns in `mcols`) or a
`data.frame` (BEDPE-like).
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) field names
are explicit:

- `x0`, `x1`:

  required — genomic position of each anchor.

- `chrom0`, `chrom1`:

  optional for GRanges (default `seqnames(data)`); required for
  `data.frame`.

- `y0`, `y1`:

  optional baselines (default `0`).

- `height`:

  optional arch peak height in data-scale units (default `1`).
