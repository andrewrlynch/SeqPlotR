# Expand a combined virtual panel into per-original-window sub-panels

When a track was laid out with `combine_windows = TRUE`, the layout
returns one panel whose xscale is virtual `c(1, virtual_total)`. For
axis drawing we want each original window to render its own ticks,
labels, and title at the corresponding virtual sub-range — so produce a
list of pseudo-panels covering each `vmap_x` entry, with `xscale` set to
the original genomic range and `inner.x0/x1` narrowed to the virtual
sub-range's npc extent.

## Usage

``` r
.expand_panels_for_combined_axis(panels)
```

## Arguments

- panels:

  The track's panel list as returned by
  [`.build_track_panels()`](http://andrewlynch.io/SeqPlotR/reference/dot-build_track_panels.md).

## Value

Either the input `panels` unchanged (when no virtual map is present) or
a list of virtual sub-panels suitable for x-axis drawing.
