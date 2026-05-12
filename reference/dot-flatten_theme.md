# Flatten a (possibly nested) SeqAes to a dotted-key map

Recursively walks a `SeqAes` or named list, producing a single-level
named list where each name is a dot-separated path to a leaf value.
Leaves are anything that is not itself a `SeqAes` or a named `list`:
scalars, numeric vectors, functions, palettes, `SeqScale*` objects,
`NA`, etc. Unnamed elements are dropped.

## Usage

``` r
.flatten_theme(x, prefix = "")
```

## Arguments

- x:

  A `SeqAes`, named list, or `NULL`.

- prefix:

  Internal — used during recursion to accumulate the path.

## Value

A named list of leaf values, keyed by dotted paths.
