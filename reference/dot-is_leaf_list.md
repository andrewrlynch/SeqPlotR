# Decide whether a named list should be treated as a theme leaf

Some leaf values are themselves named lists (e.g. the `margins` key on
the default theme stores `list(top, right, bottom, left)`). Heuristic:
if every element of the list is atomic and length 1 and the name set is
a subset of known leaf-list names, treat it as a leaf.

## Usage

``` r
.is_leaf_list(v)
```

## Arguments

- v:

  A named list.

## Value

`TRUE` if `v` should be emitted as a single leaf value.
