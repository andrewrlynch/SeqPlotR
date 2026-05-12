# Build layout bounds for the patchwork-string mode

Looks up each track's `track_id` in the parsed layout regions; tracks
not present in the layout string are silently skipped.

## Usage

``` r
.build_patchwork_layout(tracks, layout_str, aesthetics)
```

## Arguments

- tracks:

  Flat list of `SeqTrackR6` objects.

- layout_str:

  The raw layout string.

- aesthetics:

  Merged plot-level aesthetics.

## Value

A list with `panelBounds` and `trackBounds`, each a *named* list keyed
by `track_id`.
