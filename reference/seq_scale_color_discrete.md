# Discrete color scale for SeqAes mappings

Discrete color scale for SeqAes mappings

## Usage

``` r
seq_scale_color_discrete(values = NULL, palette = NULL, na_value = "grey80")

seq_scale_fill_discrete(...)
```

## Arguments

- values:

  Optional named character vector mapping level names to colors.

- palette:

  Optional function `function(n)` returning `n` hex colors. Falls back
  to
  [`flexoki_palette()`](http://andrewlynch.io/SeqPlotR/reference/flexoki_palette.md)
  if `NULL`.

- na_value:

  Color for NA / unmatched values (default `"grey80"`).

- ...:

  Arguments forwarded from `seq_scale_fill_discrete()` to
  `seq_scale_color_discrete()`.

## Value

A `SeqScaleDiscrete` / `SeqScale` object.
