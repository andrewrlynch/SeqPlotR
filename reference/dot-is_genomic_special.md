# Is a `map()` expression a genomic special (bare symbol)?

Returns `TRUE` when `expr` is a bare symbol matching one of `start`,
`end`, `mid`, or `width` — the coordinate specials injected into the
mapping eval environment. Used to auto-detect genomic axes.

## Usage

``` r
.is_genomic_special(expr)
```

## Arguments

- expr:

  An R language object from a `SeqMap`.

## Value

Logical scalar.
