# Bezier arch with vertical stems

Same anchor / `t0`-`t1` rules as
[`seq_arc()`](http://andrewlynch.io/SeqPlotR/reference/seq_arc.md), but
draws stems from the baseline up to the arch endpoints and renders stubs
(with a partner chromosome label) for half-visible links. Stub rendering
is controlled by `aes(plotStubs = TRUE)` (default).

## Usage

``` r
seq_arch(
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

A `SeqArchR6` instance.
