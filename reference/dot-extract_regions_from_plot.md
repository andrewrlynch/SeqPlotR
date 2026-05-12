# Extract layout regions from a seq_plot object

Handles both patchwork (layout_str) and positional (rows) modes.

## Usage

``` r
.extract_regions_from_plot(plot_obj, margin = 0.05)
```

## Arguments

- plot_obj:

  A `SeqPlot` R6 object.

- margin:

  Numeric fractional margin.

## Value

Named list of npc bounding boxes per track ID.
