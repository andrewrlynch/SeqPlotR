# Layout the two axis slots available on each side of a track

For each side (top / bottom / left / right) there are up to two axes (x1
& x2 or y1 & y2) that could live there. The slot nearer the panel is
slot 1; the outer slot is slot 2. When both axes with the same dimension
share a side, x1 takes slot 1 and x2 takes slot 2; otherwise the single
axis takes slot 1 (the full band).

## Usage

``` r
.compute_axis_bands(resolved_theme, first_panel)
```

## Arguments

- resolved_theme:

  The track's resolved theme.

- first_panel:

  The track's first panel (source of band bounds).

## Value

A named list with one entry per side (`bottom`, `top`, `left`, `right`).
Each entry is a list of up to two slots; each slot has `anchor` (the npc
coordinate of the panel-adjacent edge), `outer` (the outer edge of the
slot), and `band_outer` (the outer edge of the whole track_inner margin
band).
