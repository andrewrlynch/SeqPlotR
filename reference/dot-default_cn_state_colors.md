# Default CN state colour palette

Diverging palette for integer copy-number states: blue at CN=0, grey at
the diploid value (2), and increasingly warm colours above. Values
beyond the range fall back to `"grey80"`.

## Usage

``` r
.default_cn_state_colors()
```

## Value

A named character vector keyed by CN state (as a string).
