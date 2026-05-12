# Format a numeric tick label

Uses decimal-comma-grouped formatting, no scientific notation, and trims
whitespace. A separate helper so the axis draw code reads clean.

## Usage

``` r
.axis_fmt(x)
```

## Arguments

- x:

  Numeric vector of tick values.

## Value

Character vector of formatted labels.
