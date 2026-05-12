# Convert parsed layout to npc bounding boxes

Convert parsed layout to npc bounding boxes

## Usage

``` r
.regions_from_parse(parsed, margin = 0.05)
```

## Arguments

- parsed:

  Output of
  [`.parse_layout_string()`](http://andrewlynch.io/SeqPlotR/reference/dot-parse_layout_string.md).

- margin:

  Numeric fractional margin on each side.

## Value

Named list of `list(x0, x1, y0, y1)` per track ID.
