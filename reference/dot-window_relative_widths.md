# Build per-window relative widths and resolve scale factors

Build per-window relative widths and resolve scale factors

## Usage

``` r
.window_relative_widths(win, track = NULL)
```

## Arguments

- win:

  A `GRanges` object of windows.

- track:

  Optional `SeqTrackR6`. When provided, `window_scale` is read from the
  track and the track-level inferred unit is used as the default. When
  `NULL`, falls back to `1e-6` (original behaviour).

## Value

A list with `rel` and `scale` numeric vectors of length `length(win)`.
