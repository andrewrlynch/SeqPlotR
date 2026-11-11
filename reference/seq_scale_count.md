# Count position scale with an optional pseudo-log transform

Builds a continuous y scale for copy-number style counts with integer
tick labels. In log mode values up to 1 stay linear and larger values
map to `1 + log2(v)`, so ticks fall at 0, 1, 2, 4, 8, ... When a track
uses this scale, elements added to it transform their mapped y fields
(see the `y_fields` field of
[`SeqElementR6`](http://andrewlynch.io/SeqPlotR/reference/SeqElementR6.md))
before mapping to the canvas. Constant aesthetics such as
`aes(height = )` are not transformed.

## Usage

``` r
seq_scale_count(
  max,
  log = max > 8,
  bottom = 0.018,
  top_pad = 0.5,
  expand = c(0.025, 0),
  oob = c("exclude", "perimeter")
)
```

## Arguments

- max:

  Largest count to show. Values below 1 are raised to 1.

- log:

  Logical. Use the pseudo-log transform. Default `max > 8`.

- bottom:

  Fraction of the top value reserved below 0, so that drawings offset
  beneath 0 stay visible. Default `0.018`.

- top_pad:

  Headroom above the top tick in linear mode. Default `0.5`.

- expand:

  Length-2 `c(mul, add)` padding, as in
  [`seq_scale_continuous()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_continuous.md).

- oob:

  Out-of-bounds handling, as in
  [`seq_scale_continuous()`](http://andrewlynch.io/SeqPlotR/reference/seq_scale_continuous.md).

## Value

A `SeqScaleCount` / `SeqScaleContinuous_Pos` / `SeqPositionScale` object
with extra fields `transform` and `inverse` (functions, or `NULL` in
linear mode) and `log`.

## Examples

``` r
s <- seq_scale_count(12, log = TRUE)
s$labels
#> [1] "0"  "1"  "2"  "4"  "8"  "16"
s$transform(12)
#> [1] 4.584963
seq_scale_count(8)$breaks
#> [1] 0 1 2 3 4 5 6 7 8
```
