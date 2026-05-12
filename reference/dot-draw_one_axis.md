# Draw a single axis (line, ticks, minor ticks, labels, title)

Draw a single axis (line, ticks, minor ticks, labels, title)

## Usage

``` r
.draw_one_axis(side, track, panels, slot)
```

## Arguments

- side:

  One of `"x1"`, `"x2"`, `"y1"`, `"y2"`.

- track:

  The `SeqTrackR6`.

- panels:

  The track's panel list.

- slot:

  The slot from
  [`.compute_axis_bands()`](http://andrewlynch.io/SeqPlotR/reference/dot-compute_axis_bands.md)
  for this axis.

## Value

Invisibly `NULL`.
