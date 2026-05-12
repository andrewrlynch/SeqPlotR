# Normalise a data y value to the unit interval within a scale range

Normalise a data y value to the unit interval within a scale range

## Usage

``` r
.data_to_npc_y(y, yscale)
```

## Arguments

- y:

  Numeric vector of data y values.

- yscale:

  Numeric vector of length 2: `c(min, max)`.

## Value

Numeric vector of the same length as `y`, clamped to `[0, 1]`.
