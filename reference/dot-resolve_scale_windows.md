# Resolve a genomic scale's windows from its enclosing track

When
[`seq_scale_genomic()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_genomic.md)
is used inside `seq_track(scale_x = ...)` without explicit `windows`,
the scale should fall back to the track's windows. This helper mutates
the scale in place if needed and recomputes `scale_factor` from the
resolved windows.

## Usage

``` r
.resolve_scale_windows(scale, track_windows)
```

## Arguments

- scale:

  A `SeqPositionScale` or `NULL`.

- track_windows:

  A `GRanges` of windows from the enclosing track, or `NULL`.

## Value

The (possibly updated) scale.
