# Draw the track background / border and per-window panel chrome

Reads from `track$resolved_theme$chrome`.

## Usage

``` r
.draw_track_chrome(track, panels, tb)
```

## Arguments

- track:

  The `SeqTrackR6`.

- panels:

  Panel list for this track.

- tb:

  Track bounds list (`x0`, `x1`, `y0`, `y1`).

## Value

Invisibly `NULL`.
