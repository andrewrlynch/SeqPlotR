# Remap a per-window panel list so xscale/yscale reflect the element's primary/secondary axis selection

When an element sets `map(axis.x = 2)` or `map(axis.y = 2)`, its
coordinate transforms must use `xscale2` / `yscale2` instead of the
primary scales. Every primitive's `prep()` reads `panel$xscale` /
`panel$yscale` directly, so this helper builds a per-element shallow
copy of the panel list with those fields swapped.

## Usage

``` r
.panels_for_element(layout_track, resolved)
```

## Arguments

- layout_track:

  List of per-window panel metadata.

- resolved:

  Element `resolved` list (contains `axis_x`/`axis_y`).

## Value

A layout_track with `xscale` / `yscale` rewritten if needed.
