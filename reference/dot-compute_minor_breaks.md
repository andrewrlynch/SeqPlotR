# Compute minor breaks between majors

If `user_val` is a numeric vector of length \> 1, it is returned as-is
(filtered to `plot_range`). If it is a single integer `N`, `N`
equally-spaced interior subdivisions are placed between adjacent major
breaks. If `NULL`, returns `NULL` (no minor breaks).

## Usage

``` r
.compute_minor_breaks(user_val, major_breaks, plot_range)
```

## Arguments

- user_val:

  User-supplied minor_breaks value, or `NULL`.

- major_breaks:

  Numeric vector of major-break positions.

- plot_range:

  Length-2 numeric vector of the visible plot range.

## Value

Numeric vector of minor-break positions, or `NULL`.
