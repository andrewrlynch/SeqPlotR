# Infer a continuous x-range for a track from its elements' resolved data

Used when the track has a non-genomic x mapping but no explicit
`scale_x`. Tries each element in turn: resolves against the track's
data + mapping and returns the range of the first numeric `x` vector.
Returns `NULL` when nothing usable is found.

## Usage

``` r
.infer_x_range(track)
```

## Arguments

- track:

  A `SeqTrackR6` instance.

## Value

Numeric length-2 vector or `NULL`.
