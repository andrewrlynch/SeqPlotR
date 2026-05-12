# Expand a multi-y-window panel into one pseudo-panel per y sub-panel

Mirror of
[.expand_panels_for_combined_axis](http://andrewlynch.io/SeqPlotR/reference/dot-expand_panels_for_combined_axis.md)
but for the y-axis. Two modes:

- `y_sub_panels` (vertically stacked sub-panels): each entry gets its
  own npc band with the corresponding window's yscale.

- `virtual_map_y` (combined virtual y axis): the single y band is split
  into per-original-window sub-ranges along the virtual axis, each
  carrying the original genomic yscale for label rendering.

## Usage

``` r
.expand_panels_for_multi_y_axis(panels)
```

## Arguments

- panels:

  The track's panel list.

## Value

A list of pseudo-panels (one per y sub-panel); or the input panels
unchanged when no multi-y data is present.
