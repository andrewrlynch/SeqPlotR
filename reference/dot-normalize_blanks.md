# Translate bare NA values on structural axis keys to `visible = FALSE`

Mirrors ggplot's `element_blank()`: when a user writes
`aes(axis.x.line = NA)` the intent is to hide that piece. After
flattening, scan the theme for the known structural sub-keys
(`axis.<side>.line`, `.title`, `.ticks`, `.labels`, `.text`,
`.gridline`) whose value is bare `NA`, and rewrite as
`<key>.visible = FALSE` (dropping the raw NA leaf so it does not confuse
downstream resolvers). Also recognises a bare string on
`axis.<side>.title` as a shorthand for `axis.<side>.title.label`.

## Usage

``` r
.normalize_blanks(flat)
```

## Arguments

- flat:

  A flat theme map (from
  [`.flatten_theme()`](http://andrewlynch.io/SeqPlotR/reference/dot-flatten_theme.md)).

## Value

The transformed flat theme.
