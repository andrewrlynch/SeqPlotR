# Resolve a theme key from a flat theme map

Walks
[`.axis_inheritance_chain()`](http://andrewlynch.io/SeqPlotR/reference/dot-axis_inheritance_chain.md)
and returns the value stored at the most specific key present. Falls
back to `default` if no entry exists.

## Usage

``` r
.resolve_theme(flat_theme, key_path, default = NULL)
```

## Arguments

- flat_theme:

  Named list of dotted keys (from
  [`.flatten_theme()`](http://andrewlynch.io/SeqPlotR/reference/dot-flatten_theme.md),
  usually merged with `default_theme`).

- key_path:

  Dotted key string (e.g. `"axis.x1.line.col"`).

- default:

  Value returned when no entry is found.

## Value

The resolved value, or `default`.
