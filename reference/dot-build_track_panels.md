# Build the panel metadata list for one track in one (x, y) bounding box

Splits the track's x-band into per-window panels and packages each one
with `full`, `inner`, `xscale`, `yscale`, and other downstream-required
fields. Used by both layout builders.

## Usage

``` r
.build_track_panels(track, x0, x1, y0, y1, window_gap, track_key)
```

## Arguments

- track:

  A `SeqTrackR6` instance.

- x0, x1:

  Numeric npc x bounds for the track.

- y0, y1:

  Numeric npc y bounds for the track.

- window_gap:

  Numeric npc gap between adjacent window panels.

- track_key:

  The value to store in each panel's `track` field — integer index in
  positional mode, `track_id` string in patchwork mode.

## Value

A list of per-window panel metadata.
