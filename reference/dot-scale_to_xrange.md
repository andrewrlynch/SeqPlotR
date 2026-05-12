# Compute an x-range from a scale object (helper shared by x1 and x2)

Compute an x-range from a scale object (helper shared by x1 and x2)

## Usage

``` r
.scale_to_xrange(sx, window_gr, w)
```

## Arguments

- sx:

  A `SeqPositionScale` or `NULL`.

- window_gr:

  A `GRanges` — used as the genomic fallback.

- w:

  Integer index of the window.

## Value

Numeric length-2 vector, or `NULL` when `sx` is NULL and no genomic
fallback is appropriate (caller decides).
