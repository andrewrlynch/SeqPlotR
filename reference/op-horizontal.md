# Append a track to the current row (horizontal stacking)

Convenience alias. Equivalent to adding a
[`seq_track()`](http://andrewlynch.io/SeqPlotR/reference/seq_track.md)
with `direction = "right"`. All logic is in `%+%`.

## Usage

``` r
e1 %|% e2
```

## Arguments

- e1:

  A `SeqPlot` object.

- e2:

  A `SeqTrack` (direction forced to `"right"`) or any valid `%+%` RHS.

## Value

`e1`, invisibly modified in place.
