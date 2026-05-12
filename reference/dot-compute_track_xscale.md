# Compute the xscale for a single track window

Dispatches on `track$scale_x`: continuous uses `limits` (or a `c(0, 1)`
placeholder when `limits` is NULL), discrete uses
`c(0.5, n_levels + 0.5)`, and everything else (including an explicit
[`seq_scale_genomic()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_genomic.md))
falls back to the genomic window range.

## Usage

``` r
.compute_track_xscale(track, window_gr, w)
```

## Arguments

- track:

  A `SeqTrackR6` instance.

- window_gr:

  The track's `windows` `GRanges`.

- w:

  Integer index of the window.

## Value

Numeric length-2 vector: `c(x_min, x_max)`.
