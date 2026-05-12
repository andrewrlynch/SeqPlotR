# Translate a SeqAes object to grid::gpar()

Internal helper. Maps SeqAes field names onto grid graphical parameter
names. `NULL` values are silently dropped by
[`grid::gpar()`](https://rdrr.io/r/grid/gpar.html) — this is
intentional, so absent fields fall through to grid defaults.

## Usage

``` r
.aes_to_gpar(a)
```

## Arguments

- a:

  A `SeqAes` object.

## Value

A `gpar` object.
