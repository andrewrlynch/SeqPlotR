# Legend Key

Constructs a single legend key entry for a SeqPlotR element. A
`LegendKey` records the visual properties (colour, shape, fill, etc.),
the display label, and an optional group title for one row in a rendered
legend.

## Usage

``` r
LegendKey(
  label = NULL,
  title = NULL,
  color = "#1C1B1A",
  shape = "-",
  size = 1,
  alpha = 1,
  fill = NULL,
  lty = 1,
  ...
)
```

## Arguments

- label:

  Character or `NULL`. Text label shown beside the key glyph.

- title:

  Character or `NULL`. Legend group title this key belongs to.

- color:

  Character. Stroke / point colour. Default `"#1C1B1A"`.

- shape:

  Character. Glyph shape code (e.g. `"-"`, `"circle"`). Default `"-"`.

- size:

  Numeric. Relative glyph size. Default `1`.

- alpha:

  Numeric in `[0, 1]`. Opacity. Default `1`.

- fill:

  Character or `NULL`. Fill colour for filled glyphs. Default `NULL`.

- lty:

  Integer or character. Line type. Default `1`.

- ...:

  Additional fields stored verbatim in the `extra` sub-list.

## Value

A `LegendKey` S3 object (a named list with class `"LegendKey"`).
