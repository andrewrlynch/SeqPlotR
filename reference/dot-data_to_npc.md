# Normalise a data x value to the unit interval within a scale range

Normalise a data x value to the unit interval within a scale range

## Usage

``` r
.data_to_npc(x, xscale)
```

## Arguments

- x:

  Numeric vector of data x values.

- xscale:

  Numeric vector of length 2: `c(min, max)`.

## Value

Numeric vector of the same length as `x`, clamped to `[0, 1]`.
