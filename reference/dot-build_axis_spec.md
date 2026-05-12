# Precompute a per-side axis spec from a flat theme

For one of `c("x1","x2","y1","y2")`, walk the hierarchy once per leaf
key and package the resolved values into a nested list the axis-draw
helpers consume directly. This avoids repeated chain walks at draw time.

## Usage

``` r
.build_axis_spec(flat_theme, side)
```

## Arguments

- flat_theme:

  Flattened theme (already merged plot+track).

- side:

  One of `"x1"`, `"x2"`, `"y1"`, `"y2"`.

## Value

A list with elements `position`, `line`, `ticks`, `text`, `title`,
`scale`, plus convenience scalars `axis_dim` and `axis_index`.
