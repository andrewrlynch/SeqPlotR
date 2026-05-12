# Genomic yscale for a track

Returns the yscale that should be used for a track's y-axis. When the
track uses a genomic y-axis (`uses_genomic_y == TRUE`), the range spans
the full extent of `track$y_windows`. Otherwise returns `c(0, 1)` as a
placeholder that callers may override using element data.

## Usage

``` r
.genomic_yscale(track)
```

## Arguments

- track:

  A `SeqTrackR6` instance.

## Value

Numeric length-2 vector: `c(y_min, y_max)`.
