# Build a `SeqMap` from a named list of column names

Helper used by the wrapper functions to construct mappings
programmatically —
[`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) captures
unevaluated expressions, which is awkward when the column name is
determined at runtime. This builds the same object shape directly from
string column names (or already-unquoted symbols).

## Usage

``` r
.map_from_names(...)
```

## Arguments

- ...:

  Named arguments. Strings are converted to bare R symbols; language /
  call objects are passed through unchanged.

## Value

A `SeqMap` object.
