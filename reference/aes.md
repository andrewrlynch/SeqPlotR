# Capture evaluated constant aesthetics

Captures static aesthetic values (colors, line widths, sizes) that do
not vary per observation. For per-observation, data-driven aesthetics,
use [`map()`](http://andrewlynch.io/SeqPlotR/reference/map.md) instead.

## Usage

``` r
aes(...)
```

## Arguments

- ...:

  Named values (e.g. `color = "blue"`, `linewidth = 1.5`).

## Value

A `SeqAes` object: a named list of evaluated values.

## Examples

``` r
a <- aes(color = "blue", linewidth = 1.5)
```
