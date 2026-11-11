# Legend keys for interval glyph types

Builds one `LegendKey` per glyph type, each drawing a miniature of the
glyph through the key's `draw` function. Use it to label
[`seq_interval_glyph()`](http://andrewlynch.io/SeqPlotR/reference/seq_interval_glyph.md)
or
[`seq_cross_arch()`](http://andrewlynch.io/SeqPlotR/reference/seq_cross_arch.md)
glyphs in a legend.

## Usage

``` r
seq_glyph_legend_keys(
  glyphs = character(),
  labels = NULL,
  all = FALSE,
  color = "grey30"
)
```

## Arguments

- glyphs:

  Character vector of glyph types (`"point"`, `"latest"`, `"earliest"`,
  `"interval"`, `"untimed"`); others, including `"hidden"`, are ignored.
  Keys follow that fixed order.

- labels:

  Optional named character overriding the default labels
  (`point = "point estimate"`, `latest = "upper bound"`,
  `earliest = "lower bound"`, `interval = "interval"`,
  `untimed = "no estimate"`).

- all:

  Logical. Return keys for all five types regardless of `glyphs`.

- color:

  Glyph colour in the keys.

## Value

A named list of `LegendKey` objects (names are glyph types).

## Examples

``` r
keys <- seq_glyph_legend_keys(c("point", "latest"))
names(keys)
#> [1] "point"  "latest"
```
