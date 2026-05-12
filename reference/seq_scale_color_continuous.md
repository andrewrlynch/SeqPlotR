# Continuous color scale for SeqAes mappings

Continuous color scale for SeqAes mappings

## Usage

``` r
seq_scale_color_continuous(
  palette = "viridis",
  limits = NULL,
  na_value = "grey80",
  midpoint = NULL
)

seq_scale_fill_continuous(...)
```

## Arguments

- palette:

  One of `"viridis"`, `"plasma"`, `"magma"`, `"blues"`, `"reds"`.

- limits:

  Optional numeric vector of length 2 clamping the scale range.

- na_value:

  Color for NA values (default `"grey80"`).

- midpoint:

  Optional midpoint for diverging scales (not yet implemented).

- ...:

  Arguments forwarded from `seq_scale_fill_continuous()` to
  `seq_scale_color_continuous()`.

## Value

A `SeqScaleContinuous` / `SeqScale` object.
