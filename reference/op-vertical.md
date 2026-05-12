# Start a new row (vertical stacking)

Convenience alias. Equivalent to adding a
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
with `direction = "under"`. All logic is in `%+%`.

## Usage

``` r
e1 %__% e2
```

## Arguments

- e1:

  A `SeqPlot` object.

- e2:

  A `SeqTrack` (direction forced to `"under"`) or any valid `%+%` RHS.

## Value

`e1`, invisibly modified in place.
