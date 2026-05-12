# Merge theme-shortcut scale fields into an (optional) scale object

The theme can express scale settings via keys like
`axis.x1.scale.limits` / `.breaks` / `.expand` / `.cap`. This helper
builds or augments the canonical scale object:

## Usage

``` r
.merge_scale_with_theme(scale, axis_spec)
```

## Arguments

- scale:

  An existing scale or `NULL`.

- axis_spec:

  A per-side axis spec from
  [`.build_axis_spec()`](http://andrewlynch.io/SeqPlotR/reference/dot-build_axis_spec.md).

## Value

An updated scale, or `NULL` when neither source supplies one.

## Details

- If `scale` is `NULL` and the theme has at least one relevant entry, a
  [`seq_scale_continuous()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_continuous.md)
  is constructed from the theme values.

- If `scale` is non-`NULL`, any `NULL` fields on `scale` are filled from
  the theme. User-supplied fields always win.
