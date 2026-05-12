# Strand-classified SV reconstruction arches

Inherits everything from
[`seq_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_arch.md) and
additionally classifies each link by strand orientation and chromosome
pair, drawing each class on a fixed tier with a class-specific color.
Use to summarise structural variant calls.

## Usage

``` r
seq_recon(
  data = NULL,
  mapping = NULL,
  t0 = NULL,
  t1 = NULL,
  aesthetics = aes(),
  drawClasses = c("Inversion", "Dup/Del", "Translocation"),
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

- drawClasses:

  Character vector of tier labels to render (default:
  `c("Inversion", "Dup/Del", "Translocation")`).

- ...:

  Reserved.

## Value

A `SeqReconR6` instance.

## Details

Default colors come from
[`flexoki_palette()`](http://andrewlynch.io/SeqPlotR/reference/flexoki_palette.md)
and can be overridden via
`aes(h2hColor=, t2tColor=, dupColor=, delColor=, transColor=)`.

Both `strand0` and `strand1` must be specified in
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) — `seq_recon`
errors at `prep()` time if either is absent.
