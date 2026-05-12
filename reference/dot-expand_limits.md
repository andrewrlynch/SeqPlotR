# Expand a data range by a `c(mul, add)` specification

Returns `c(lo - mul*span - add, hi + mul*span + add)` where
`span = diff(r)`. Used to introduce a small "breath" around the data so
points are not crammed against the panel edge.

## Usage

``` r
.expand_limits(r, expand = c(0, 0))
```

## Arguments

- r:

  Length-2 numeric vector (the data range).

- expand:

  Length-1 or length-2 numeric. Scalars become `c(mul, 0)`.

## Value

Length-2 numeric vector.
