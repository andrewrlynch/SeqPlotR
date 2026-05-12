# Normalise a margin specification

Accepts a scalar (applied to all four sides), a length-4 numeric in
base-R order `c(bottom, left, top, right)`, or an already-normalised
list with fields `top`, `right`, `bottom`, `left`. Returns a list in the
`list(top, right, bottom, left)` canonical form.

## Usage

``` r
.normalize_margin(m)
```

## Arguments

- m:

  A scalar, length-4 numeric, or list.

## Value

A list with fields `top`, `right`, `bottom`, `left`.
