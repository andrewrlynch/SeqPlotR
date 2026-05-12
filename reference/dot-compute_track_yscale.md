# Compute the y-axis scale for a track

Used by the layout builders. Dispatches on (1) genomic y-windows, (2)
explicit `scale_y` type, and (3) a `c(0, 1)` placeholder fallback.

## Usage

``` r
.compute_track_yscale(track)
```

## Arguments

- track:

  A `SeqTrackR6` instance.

## Value

Numeric length-2 vector `c(y_min, y_max)`.
